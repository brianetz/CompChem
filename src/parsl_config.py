"""
Parsl execution configs for the CompChem pipeline.

get_local_config() is the Phase 1 config: single-machine, no scheduler,
multiple concurrent workers. Later HPC configs (Frontier/Defiant, using
SlurmProvider instead of LocalProvider) get added here as their own
functions - see the Phase 8 plan in the architecture doc - without any
pipeline code that calls parsl.load(...) needing to change.
"""
from __future__ import annotations

import os

from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.providers import LocalProvider


def get_local_config(max_workers: int | None = None) -> Config:
    """Local multi-worker config for running on this machine (no scheduler).

    max_workers defaults to (cpu_count - 1), so independent jobs (e.g.
    sibling BDE fragment calculations) actually run concurrently instead of
    one at a time, while leaving a core free for the OS/foreground work.
    """
    if max_workers is None:
        max_workers = max(1, (os.cpu_count() or 2) - 1)

    return Config(
        executors=[
            HighThroughputExecutor(
                label="htex_local",
                max_workers_per_node=max_workers,
                provider=LocalProvider(
                    init_blocks=1,
                    max_blocks=1,
                ),
            )
        ],
    )
