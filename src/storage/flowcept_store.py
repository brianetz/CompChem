"""
Flowcept provenance helpers for the CompChem pipeline.

Centralizes how the rest of the codebase talks to Flowcept so individual
modules don't need to duplicate import-guards or know which capture mode
(offline JSONL buffer vs. online Mongo/Redis) is configured on this machine.

One-time setup per machine:
    pip install flowcept
    flowcept --init-settings
    flowcept --config-profile full-offline -y   # local dev: no Mongo/Redis needed

Docs: https://github.com/ORNL/flowcept
"""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, Optional

try:
    from flowcept import Flowcept as _Flowcept
    from flowcept import flowcept_task as _flowcept_task
    _FLOWCEPT_AVAILABLE = True
except ImportError:
    _Flowcept = None
    _FLOWCEPT_AVAILABLE = False

    def _flowcept_task(func=None, **_kwargs):
        """No-op fallback so the pipeline still runs if flowcept isn't installed yet."""
        if func is None:
            return lambda f: f
        return func


# Every pipeline module should import flowcept_task from here, not from
# `flowcept` directly - keeps the "is flowcept installed" fallback in one place.
flowcept_task = _flowcept_task


@contextmanager
def provenance_session() -> Iterator[None]:
    """Wrap a block of pipeline code so every @flowcept_task call inside it is
    captured under one Flowcept run. No-ops cleanly if flowcept isn't
    installed, so it's always safe to wrap code in this.
    """
    if not _FLOWCEPT_AVAILABLE:
        yield
        return
    with _Flowcept():
        yield


def read_captured_tasks(buffer_path: str = "flowcept_buffer.jsonl") -> Optional[Any]:
    """Read back captured tasks from the local offline JSONL buffer (the
    full-offline config profile). Returns a pandas DataFrame, or None if
    flowcept isn't installed or no buffer file exists yet.
    """
    if not _FLOWCEPT_AVAILABLE:
        return None
    if not Path(buffer_path).exists():
        return None
    return _Flowcept.read_buffer_file(buffer_path, return_df=True)


def query_online_tasks(**filters: Any) -> Optional[Any]:
    """Query captured tasks from the online DB (MongoDB), for when the
    full-online / online-database config profile is used instead of
    full-offline (e.g. later on HPC, or for shared/team use).
    """
    if not _FLOWCEPT_AVAILABLE:
        return None
    return _Flowcept.db.query(collection="tasks", filter=filters or None)
