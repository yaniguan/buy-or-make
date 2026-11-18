"""On-disk SQLite cache for API responses with fetch timestamps."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any, Optional

from .schema import utc_now

DEFAULT_CACHE_PATH = Path("data/cache.sqlite")


class ResponseCache:
    """Persistent SQLite key-value cache storing JSON payloads and fetch timestamps."""

    def __init__(self, path: Path | str = DEFAULT_CACHE_PATH) -> None:
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
        self.hits = 0
        self.misses = 0

    def get(self, cache_key: str) -> Optional[tuple[Any, str]]:
        """Return (decoded_payload, fetched_utc) or None."""
        row = self._conn.execute(
            "SELECT payload_json, fetched_utc FROM api_cache WHERE cache_key = ?",
            (cache_key,),
        ).fetchone()
        if row is None:
            self.misses += 1
            return None
        self.hits += 1
        return json.loads(row[0]), row[1]

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
        return ts

    def close(self) -> None:
        self._conn.close()

    def __enter__(self) -> "ResponseCache":
        return self

    def __exit__(self, *args) -> None:
        self.close()
