"""On-disk SQLite cache for API responses with fetch timestamps and optional TTL."""

from __future__ import annotations

import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from .schema import utc_now

DEFAULT_CACHE_PATH = Path("data/cache.sqlite")
DEFAULT_ONLINE_TTL_SECONDS = 30.0 * 86400.0  # 30 days


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
    """Persistent SQLite key-value cache storing JSON payloads and fetch timestamps."""

    def __init__(
        self,
        path: Path | str = DEFAULT_CACHE_PATH,
        *,
        max_age_seconds: Optional[float] = None,
        refresh: bool = False,
    ) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(str(self.path))
        self._conn.execute(
            """
            CREATE TABLE IF NOT EXISTS api_cache (
                cache_key TEXT PRIMARY KEY,
                source_id TEXT NOT NULL,
                payload_json TEXT NOT NULL,
                fetched_utc TEXT NOT NULL
            )
            """
        )
        self._conn.commit()
        self.max_age_seconds = max_age_seconds
        self.refresh = refresh
        self._session_written: set[str] = set()
        self.hits = 0
        self.misses = 0

    def get(
        self,
        cache_key: str,
        *,
        max_age_seconds: Optional[float] = None,
    ) -> Optional[tuple[Any, str]]:
        """Return `(decoded_payload, fetched_utc)` or `None` if missing, bypassed, or expired."""
        if self.refresh and cache_key not in self._session_written:
            self.misses += 1
            return None

        row = self._conn.execute(
            "SELECT payload_json, fetched_utc FROM api_cache WHERE cache_key = ?",
            (cache_key,),
        ).fetchone()
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
        self._conn.execute(
            """
            INSERT INTO api_cache (cache_key, source_id, payload_json, fetched_utc)
            VALUES (?, ?, ?, ?)
            ON CONFLICT(cache_key) DO UPDATE SET
                source_id = excluded.source_id,
                payload_json = excluded.payload_json,
                fetched_utc = excluded.fetched_utc
            """,
            (cache_key, source_id, json.dumps(payload), ts),
        )
        self._conn.commit()
        self._session_written.add(cache_key)
        return ts

    def close(self) -> None:
        self._conn.close()

    def __enter__(self) -> "ResponseCache":
        return self

    def __exit__(self, *args) -> None:
        self.close()
