"""On-disk SQLite cache for API responses with WAL mode, thread safety, and read-only URI support."""

from __future__ import annotations

import json
import random
import sqlite3
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Optional, TypeVar

from .schema import utc_now

DEFAULT_CACHE_PATH = Path("data/cache.sqlite")
DEFAULT_ONLINE_TTL_SECONDS = 30.0 * 86400.0  # 30 days
DEFAULT_BUSY_TIMEOUT_MS = 30_000

_T = TypeVar("_T")


def open_readonly_sqlite(path: Path, *, busy_timeout_ms: int = DEFAULT_BUSY_TIMEOUT_MS) -> sqlite3.Connection:
    """Open an existing SQLite database via a thread-safe read-only URI (`?mode=ro`).

    Falls back to `?mode=ro&immutable=1` when mounted on a strictly read-only directory
    where SQLite cannot create WAL shared-memory (`-shm`) sidecar files.
    """
    resolved = path.resolve()
    base_uri = resolved.as_uri()
    timeout_sec = max(1.0, float(busy_timeout_ms) / 1000.0)
    try:
        conn = sqlite3.connect(
            f"{base_uri}?mode=ro",
            uri=True,
            timeout=timeout_sec,
            check_same_thread=False,
        )
        conn.execute(f"PRAGMA busy_timeout = {int(busy_timeout_ms)}")
        conn.execute("SELECT 1 FROM sqlite_master LIMIT 1").fetchall()
        return conn
    except sqlite3.OperationalError:
        conn = sqlite3.connect(
            f"{base_uri}?mode=ro&immutable=1",
            uri=True,
            timeout=timeout_sec,
            check_same_thread=False,
        )
        conn.execute(f"PRAGMA busy_timeout = {int(busy_timeout_ms)}")
        return conn


def _is_lock_error(exc: sqlite3.OperationalError) -> bool:
    msg = str(exc).lower()
    return "locked" in msg or "busy" in msg


def _parse_utc_timestamp(ts: str) -> Optional[datetime]:
    text = (ts or "").strip()
    if not text:
        return None
    try:
        if text.endswith("Z"):
            text = text[:-1] + "+00:00"
        dt = datetime.fromisoformat(text)
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        return dt.astimezone(timezone.utc)
    except ValueError:
        return None


class ResponseCache:
    """Thread-safe persistent SQLite key-value cache storing JSON payloads and fetch timestamps."""

    def __init__(
        self,
        path: Path | str = DEFAULT_CACHE_PATH,
        *,
        max_age_seconds: Optional[float] = None,
        refresh: bool = False,
        readonly: bool = False,
        busy_timeout_ms: int = DEFAULT_BUSY_TIMEOUT_MS,
        max_lock_retries: int = 6,
    ) -> None:
        self.path = Path(path)
        self.max_age_seconds = max_age_seconds
        self.refresh = refresh
        self.busy_timeout_ms = max(100, int(busy_timeout_ms))
        self.max_lock_retries = max(0, int(max_lock_retries))
        self._lock = threading.RLock()
        self._session_written: set[str] = set()
        self.hits = 0
        self.misses = 0
        self._readonly = bool(readonly)
        self._has_table = False
        self._conn: Optional[sqlite3.Connection] = None

        if self._readonly:
            if self.path.exists():
                self._conn = open_readonly_sqlite(self.path, busy_timeout_ms=self.busy_timeout_ms)
                self._has_table = self._check_api_cache_table(self._conn)
        else:
            try:
                self.path.parent.mkdir(parents=True, exist_ok=True)
                timeout_sec = max(1.0, float(self.busy_timeout_ms) / 1000.0)
                conn = sqlite3.connect(
                    str(self.path),
                    timeout=timeout_sec,
                    check_same_thread=False,
                )
                self._conn = conn
                self._run_with_lock_retry(lambda: self._configure_rw_connection(conn))
                self._has_table = True
            except (sqlite3.OperationalError, OSError):
                if self._conn is not None:
                    try:
                        self._conn.close()
                    except Exception:
                        pass
                    self._conn = None
                if self.path.exists():
                    self._readonly = True
                    self._conn = open_readonly_sqlite(self.path, busy_timeout_ms=self.busy_timeout_ms)
                    self._has_table = self._check_api_cache_table(self._conn)
                else:
                    raise

    @staticmethod
    def _check_api_cache_table(conn: sqlite3.Connection) -> bool:
        row = conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'api_cache'"
        ).fetchone()
        return row is not None

    def _configure_rw_connection(self, conn: sqlite3.Connection) -> None:
        conn.execute(f"PRAGMA busy_timeout = {self.busy_timeout_ms}")
        conn.execute("PRAGMA journal_mode = WAL")
        conn.execute("PRAGMA synchronous = NORMAL")
        if not self._check_api_cache_table(conn):
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS api_cache (
                    cache_key TEXT PRIMARY KEY,
                    source_id TEXT NOT NULL,
                    payload_json TEXT NOT NULL,
                    fetched_utc TEXT NOT NULL
                )
                """
            )
            conn.commit()

    def _run_with_lock_retry(self, fn: Callable[[], _T]) -> _T:
        for attempt in range(self.max_lock_retries + 1):
            try:
                with self._lock:
                    return fn()
            except sqlite3.OperationalError as exc:
                if not _is_lock_error(exc) or attempt >= self.max_lock_retries:
                    raise
                delay = random.uniform(0.005, 0.020) * (2**attempt)
                time.sleep(min(delay, 0.50))
        raise RuntimeError("Unreachable SQLite retry exhaustion")

    def get(
        self,
        cache_key: str,
        *,
        max_age_seconds: Optional[float] = None,
    ) -> Optional[tuple[Any, str]]:
        """Return `(decoded_payload, fetched_utc)` or `None` if missing, bypassed, or expired."""
        with self._lock:
            if self.refresh and cache_key not in self._session_written:
                self.misses += 1
                return None
            if self._conn is None or not self._has_table:
                self.misses += 1
                return None

        def _query() -> Optional[tuple[str, str]]:
            assert self._conn is not None
            return self._conn.execute(
                "SELECT payload_json, fetched_utc FROM api_cache WHERE cache_key = ?",
                (cache_key,),
            ).fetchone()

        row = self._run_with_lock_retry(_query)
        with self._lock:
            if row is None:
                self.misses += 1
                return None

            payload_raw, fetched_utc = row
            effective_ttl = max_age_seconds if max_age_seconds is not None else self.max_age_seconds
            if effective_ttl is not None:
                if effective_ttl <= 0:
                    self.misses += 1
                    return None
                fetched_dt = _parse_utc_timestamp(fetched_utc)
                if fetched_dt is None:
                    self.misses += 1
                    return None
                age_sec = (datetime.now(timezone.utc) - fetched_dt).total_seconds()
                if age_sec > effective_ttl:
                    self.misses += 1
                    return None

            self.hits += 1
            return json.loads(payload_raw), fetched_utc

    def set(
        self,
        cache_key: str,
        source_id: str,
        payload: Any,
        fetched_utc: Optional[str] = None,
    ) -> str:
        ts = fetched_utc or utc_now()
        encoded = json.dumps(payload)
        with self._lock:
            if self._readonly or self._conn is None:
                return ts

        def _write() -> None:
            assert self._conn is not None
            self._conn.execute(
                """
                INSERT INTO api_cache (cache_key, source_id, payload_json, fetched_utc)
                VALUES (?, ?, ?, ?)
                ON CONFLICT(cache_key) DO UPDATE SET
                    source_id = excluded.source_id,
                    payload_json = excluded.payload_json,
                    fetched_utc = excluded.fetched_utc
                """,
                (cache_key, source_id, encoded, ts),
            )
            self._conn.commit()
            self._session_written.add(cache_key)

        self._run_with_lock_retry(_write)
        return ts

    def close(self) -> None:
        with self._lock:
            if self._conn is not None:
                self._conn.close()
                self._conn = None

    def __enter__(self) -> "ResponseCache":
        return self

    def __exit__(self, *args) -> None:
        self.close()
