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
from parsl.addresses import address_by_interface
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SrunLauncher
from parsl.providers import LocalProvider, SlurmProvider


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


# Frontier compute nodes have 64 cores, 56 of them usable by jobs (the rest
# are reserved for system noise by default).
FRONTIER_USABLE_CORES = 56


def get_frontier_config(
    account: str,
    partition: str,
    nodes: int,
    walltime: str,
    workers_per_node: int,
    worker_init: str,
    interface: str = "hsn0",
) -> Config:
    """Frontier config: ONE Slurm block of `nodes` nodes, each running
    `workers_per_node` concurrent Parsl workers - so many tasks share a node
    instead of each task holding (and billing) a whole node to itself.

    The interchange runs wherever this is loaded (normally a login node);
    workers on the compute nodes connect back to it over `interface` (the
    Slingshot high-speed network, hsn0, reachable from compute nodes).
    """
    return Config(
        executors=[
            HighThroughputExecutor(
                label="htex_frontier",
                address=address_by_interface(interface),
                max_workers_per_node=workers_per_node,
                cores_per_worker=FRONTIER_USABLE_CORES / workers_per_node,
                provider=SlurmProvider(
                    account=account,
                    partition=partition,
                    nodes_per_block=nodes,
                    init_blocks=1,
                    # Not 0: the driver holds tasks back until the shared
                    # Ollama server is up, and with min_blocks=0 Parsl
                    # cancels the still-idle block after 2 min.
                    min_blocks=1,
                    max_blocks=1,
                    walltime=walltime,
                    worker_init=worker_init,
                    # One worker manager per node, given every usable core -
                    # srun's default would pin it to a single core.
                    launcher=SrunLauncher(overrides=f"--cpus-per-task={FRONTIER_USABLE_CORES}"),
                ),
            )
        ],
    )
