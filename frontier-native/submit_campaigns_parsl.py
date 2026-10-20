#!/usr/bin/env python
"""
Run agentic-campaign tasks in parallel on Frontier through Parsl, from a
Frontier login node: one planner per question (--question, repeatable, or
--questions-file), all running at once. --n-tasks K runs K copies of each
question instead - only useful to measure how repeatable the planner is.

Parsl counterpart to submit_campaigns.sh. Instead of one `sbatch` job (and so
one whole node) per task, Parsl requests ONE Slurm block of --nodes nodes and
packs --workers-per-node tasks onto each, running each task as a bash_app.
With the ollama backend, ONE shared `ollama serve` runs INSIDE that same
block, started by the block's worker_init on its first node before any
worker starts - so the server and the tasks always share one allocation and
one lifetime (no separate server job that could start, or time out, at a
different moment than the workers). It uses the node's GPU while the tasks
use its CPUs.

Each run gets its own directory, $ROOT_DIR/run_<timestamp>/, holding:
  task_<i>/agentic_campaign_report.md   (+ campaign.log, report_structures/,
                                          flowcept_buffer.jsonl)
  task_<i>.log                          each task's stdout/stderr
  ollama_serve.log                      the shared server's log
  parsl/                                Parsl's own logs and block scripts
                                        (block stdout/stderr in submit_scripts/)

This process must stay alive for the whole run (the Parsl interchange runs
inside it) - start it under tmux/screen or nohup on the login node.

Must be run with the campaign venv's python, e.g.:
    /path/to/venv/bin/python submit_campaigns_parsl.py --question "..." ...

Results database (--db) and staged mode (--staged):
  --db PATH   every task shares one results database (storage/results_db.py):
              a species already computed is reused, and tasks needing the same
              species at the same time compute it once.
  --staged    (needs --db) instead of N full pipelines (N = questions x
              --n-tasks), run three stages in the same block:
                1. plan   - N planning tasks, each writes task_<i>/spec.json
                2. compute - one task per UNIQUE species across all specs that
                            is not already in the database
                3. assemble - N tasks run each spec; every species is a cache
                            hit, so this is BDE math + critic + report only
  run_summary.json in the run directory records per-stage timings and counts.
"""
from __future__ import annotations

import argparse
import json
import os
import shlex
import sys
import time
from datetime import datetime
from pathlib import Path

USABLE_CORES = 56  # keep in sync with parsl_config.FRONTIER_USABLE_CORES


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--question", action="append", default=[],
                   help="Campaign question (commas are fine); repeat for a batch - one planner each")
    p.add_argument("--questions-file", default=None,
                   help="File with one question per line (blank lines and # comments ignored)")
    p.add_argument("--n-tasks", type=int, default=1,
                   help="Copies of EACH question (default: 1). >1 only to measure planner "
                        "repeatability: every copy plans the same question independently")
    p.add_argument("--root-dir", required=True, help="Shared output directory (created if missing)")
    p.add_argument("--repo-src", required=True, help="Absolute path to CompChem/src on Frontier")
    p.add_argument("--venv-path", required=True, help="Absolute path to the Frontier venv")
    p.add_argument("--partition", required=True, help="Slurm partition")
    p.add_argument("--account", default="STF053")
    p.add_argument("--qos", default=None,
                   help="Slurm QoS, e.g. develop (with --partition service: <=2h, 1 job)")
    p.add_argument("--nodes", type=int, default=1, help="Nodes in the Parsl block (default: 1)")
    p.add_argument("--workers-per-node", type=int, default=8,
                   help="Concurrent tasks per node; each gets 56/N cores (default: 8)")
    p.add_argument("--time-limit-minutes", type=int, default=60,
                   help="Walltime of the Parsl block - ALL tasks must finish within it (default: 60)")
    p.add_argument("--backend", choices=["ollama", "claude_cli", "anthropic"], default="ollama")
    p.add_argument("--ollama-model", default="nemotron-3-nano:latest")
    p.add_argument("--ollama-bin", default=os.path.expanduser("~/software/ollama/bin/ollama"))
    p.add_argument("--ollama-port", type=int, default=11434)
    p.add_argument("--ollama-models-dir", default=None, help="Shared pre-pulled Ollama model cache")
    p.add_argument("--ollama-num-parallel", type=int, default=None,
                   help="OLLAMA_NUM_PARALLEL (default: number of tasks). Each slot gets its own KV "
                        "cache, so keep it small for a model that nearly fills the GPUs")
    p.add_argument("--ollama-context-length", type=int, default=None,
                   help="OLLAMA_CONTEXT_LENGTH per slot (default: Ollama's VRAM-based choice, "
                        "up to 256k). Planner prompts are ~3k tokens; 8192 is enough")
    p.add_argument("--ollama-load-timeout", default=None,
                   help="OLLAMA_LOAD_TIMEOUT, e.g. 30m for a ~430 GB model read from Lustre "
                        "(Ollama default: 5m)")
    p.add_argument("--ollama-no-mmap", action="store_true",
                   help="Planners ask Ollama to load the model without mmap: a sequential "
                        "read (~2 GB/s from Lustre) instead of page faults (~60 MB/s)")
    p.add_argument("--plan-critic-skill", default=None,
                   help="Path to a plan-critic skill .md (e.g. <repo-src>/bde-plan-critic.md): "
                        "the planner's model reviews each validated BDE plan with it and a FAIL "
                        "earns the planner's one retry. Needs a --repo-src that has the critic")
    p.add_argument("--plan-chat", action="store_true",
                   help="Plan each question as a group chat of Academy agents (Planner, "
                        "Validator, plan critic with --plan-critic-skill) under a supervising "
                        "round-robin moderator (agents/plan_discussion.py)")
    p.add_argument("--plan-chat-rounds", type=int, default=3,
                   help="Maximum discussion rounds with --plan-chat (default: 3)")
    p.add_argument("--plan-board", action="store_true",
                   help="Build each plan module by module on a shared board reviewed by the "
                        "Validator and the plan critic (agents/plan_board.py)")
    p.add_argument("--plan-board-attempts", type=int, default=3,
                   help="Attempts per module with --plan-board (default: 3)")
    p.add_argument("--claude-model", default=None, help="Only used with --backend claude_cli")
    p.add_argument("--db", default=None,
                   help="Shared results database (SQLite on Lustre), e.g. "
                        "/lustre/orion/stf053/scratch/$USER/compchem-db/results.sqlite")
    p.add_argument("--staged", action="store_true",
                   help="Plan all -> compute unique missing species -> assemble (needs --db)")
    p.add_argument("--no-rag", action="store_true",
                   help="Don't give the planner prior results retrieved from --db")
    p.add_argument("--interface", default="hsn0",
                   help="Network interface compute nodes use to reach this host (default: hsn0)")
    args = p.parse_args()
    if args.staged and not args.db:
        p.error("--staged needs --db")
    questions = list(args.question)
    if args.questions_file:
        for line in Path(args.questions_file).read_text().splitlines():
            if line.strip() and not line.lstrip().startswith("#"):
                questions.append(line.strip())
    if not questions:
        p.error("give at least one --question or a --questions-file")
    if args.n_tasks < 1:
        p.error("--n-tasks must be >= 1")
    args.questions = questions
    # task index (1-based) -> question; copies of one question are adjacent
    args.task_questions = {i + 1: q for i, q in enumerate(q for q in questions for _ in range(args.n_tasks))}
    args.total_tasks = len(args.task_questions)
    return args


def ollama_worker_init(args: argparse.Namespace, run_dir: Path) -> list[str]:
    """worker_init lines that start the shared Ollama server on the block's
    first node and point every worker at it. They run in the block's batch
    script before workers are launched, so no task can start before the
    server is up and the model is confirmed present; on failure the block
    exits instead of running tasks that would all fail with a 404."""
    ollama = shlex.quote(args.ollama_bin)
    model = shlex.quote(args.ollama_model)
    port = args.ollama_port
    local = f"OLLAMA_HOST=127.0.0.1:{port}"
    lines = []
    if args.ollama_models_dir:
        lines.append(f"export OLLAMA_MODELS={shlex.quote(args.ollama_models_dir)}")
    server_env = f"OLLAMA_HOST=0.0.0.0:{port} OLLAMA_NUM_PARALLEL={args.ollama_num_parallel or args.total_tasks}"
    if args.ollama_context_length:
        server_env += f" OLLAMA_CONTEXT_LENGTH={args.ollama_context_length}"
    if args.ollama_load_timeout:
        server_env += f" OLLAMA_LOAD_TIMEOUT={shlex.quote(args.ollama_load_timeout)}"
    lines += [
        # Append: when Parsl replaces a timed-out block, the new server must
        # not overwrite the first block's log.
        f"{server_env} "
        f"{ollama} serve >> {shlex.quote(str(run_dir / 'ollama_serve.log'))} 2>&1 &",
        f"for _ in $(seq 1 120); do {local} {ollama} list >/dev/null 2>&1 && break; sleep 1; done",
        f"{local} {ollama} show {model} >/dev/null 2>&1 || "
        f"{{ echo \"Ollama server/model {args.ollama_model} not available on $(hostname)\" >&2; exit 1; }}",
        # Workers (on this node or others in the block) inherit this.
        f"export OLLAMA_HOST=$(hostname):{port}",
        *([f"export OLLAMA_CHAT_OPTIONS={shlex.quote(json.dumps({'use_mmap': False}))}"]
          if args.ollama_no_mmap else []),
        'echo "Shared Ollama server ready at $OLLAMA_HOST"',
    ]
    return lines


def backend_args(args: argparse.Namespace) -> list[str]:
    if args.backend == "ollama":
        return ["--ollama", "--model", args.ollama_model]
    if args.backend == "claude_cli":
        return ["--claude-cli"] + (["--model", args.claude_model] if args.claude_model else [])
    return []


def main() -> None:
    args = parse_args()
    sys.path.insert(0, args.repo_src)
    # Parsl launches its interchange as `interchange.py` from PATH, which only
    # has the venv's bin/ on it if the venv was activated - don't rely on that.
    os.environ["PATH"] = f"{args.venv_path}/bin:{os.environ.get('PATH', '')}"

    import parsl
    from parsl import bash_app
    from parsl_config import get_frontier_config

    run_dir = Path(args.root_dir) / f"run_{datetime.now():%Y%m%d_%H%M%S}"
    run_dir.mkdir(parents=True)
    print(f"Run directory: {run_dir}")

    threads = max(1, USABLE_CORES // args.workers_per_node)
    worker_init = "\n".join([
        f"source {shlex.quote(args.venv_path)}/bin/activate",
        f"export PYTHONPATH={shlex.quote(args.repo_src)}:${{PYTHONPATH:-}}",
        # PySCF/BLAS threads per task, so co-located tasks don't oversubscribe the node.
        f"export OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} MKL_NUM_THREADS={threads}",
        *([f"export COMPCHEM_DB={shlex.quote(args.db)}"] if args.db else []),
        f"export COMPCHEM_RUN_ID={run_dir.name}",
        *([f"export COMPCHEM_PLAN_CRITIC_SKILL={shlex.quote(args.plan_critic_skill)}"]
          if args.plan_critic_skill else []),
        *(ollama_worker_init(args, run_dir) if args.backend == "ollama" else []),
    ])
    config = get_frontier_config(
        account=args.account,
        partition=args.partition,
        nodes=args.nodes,
        walltime=f"{args.time_limit_minutes // 60:02d}:{args.time_limit_minutes % 60:02d}:00",
        workers_per_node=args.workers_per_node,
        worker_init=worker_init,
        interface=args.interface,
        qos=args.qos,
    )
    config.run_dir = str(run_dir / "parsl")

    summary = {"run_dir": str(run_dir), "mode": "staged" if args.staged else ("db" if args.db else "baseline"),
               "n_tasks": args.total_tasks, "copies_per_question": args.n_tasks,
               "nodes": args.nodes, "workers_per_node": args.workers_per_node, "backend": args.backend,
               "question": args.questions[0] if len(args.questions) == 1 else args.questions,
               "task_questions": args.task_questions, "stages": {}}
    db = None
    if args.db:
        os.environ["COMPCHEM_DB"] = args.db
        from storage.results_db import ResultsDB
        db = ResultsDB(args.db)
        summary["db"] = args.db
        summary["db_before"] = db.stats()

    results, dfk = {}, None
    t_start = time.time()
    try:
        dfk = parsl.load(config)
        print(f"Parsl block requested: {args.nodes} node(s) x {args.workers_per_node} workers")

        @bash_app
        def run_cmd(task_dir: str, cmd: str, stdout=None, stderr=None):
            return f"mkdir -p {task_dir} && cd {task_dir} && {cmd}"

        def wait_all(futures: dict, label: str) -> dict:
            out = {}
            for i, fut in futures.items():
                try:
                    fut.result()
                    out[i] = "COMPLETED"
                except Exception as e:  # BashExitFailure, worker loss, walltime, ...
                    out[i] = f"FAILED ({type(e).__name__}: {e})"
                print(f"[{label} {i}] {out[i]}")
            return out

        def campaign(i: int, extra: list[str], log: str, subdir: str = "") -> object:
            # Stage 1 runs in task_<i>/plan/ so its provenance buffer isn't
            # overwritten by the assemble stage running in task_<i>/.
            task_dir = run_dir / f"task_{i}" / subdir if subdir else run_dir / f"task_{i}"
            cmd = shlex.join(["python", "-m", "agents.run_agentic_campaign", "--question", args.task_questions[i],
                              *backend_args(args), *(["--no-rag"] if args.no_rag else []),
                              *(["--plan-chat", "--plan-chat-rounds", str(args.plan_chat_rounds)]
                                if args.plan_chat else []),
                              *(["--plan-board", "--plan-board-attempts", str(args.plan_board_attempts)]
                                if args.plan_board else []),
                              "--auto-approve", *extra])
            return run_cmd(shlex.quote(str(task_dir)), cmd, stdout=log, stderr=log)

        tasks = range(1, args.total_tasks + 1)
        if not args.staged:
            t0 = time.time()
            futures = {i: campaign(i, [], str(run_dir / f"task_{i}.log")) for i in tasks}
            print(f"Submitted {args.total_tasks} tasks; waiting (Ctrl-C cancels the whole run)...")
            results = wait_all(futures, "task")
            summary["stages"]["full"] = {"tasks": args.total_tasks, "wall_s": round(time.time() - t0, 1)}
        else:
            from agents.compute_species import collect_requests, split_cached

            # Stage 1: plan
            t0 = time.time()
            futures = {i: campaign(i, ["--plan-only", "../spec.json"], str(run_dir / f"plan_{i}.log"), "plan")
                       for i in tasks}
            print(f"[stage 1] {args.total_tasks} planning tasks submitted "
                  f"({len(args.questions)} question(s) x {args.n_tasks})")
            plan_status = wait_all(futures, "plan")
            specs = {}
            for i in tasks:
                path = run_dir / f"task_{i}" / "spec.json"
                if path.exists():
                    specs[i] = json.loads(path.read_text())
            summary["stages"]["plan"] = {"tasks": args.total_tasks, "specs_ok": len(specs),
                                         "wall_s": round(time.time() - t0, 1)}

            # Stage 2: compute unique species not already in the database
            t0 = time.time()
            spec_ids = list(specs)
            requests, total, skipped = collect_requests([specs[i] for i in spec_ids])
            hits, missing = split_cached(requests, db)
            req_file = run_dir / "species_requests.json"
            req_file.write_text(json.dumps(missing, indent=2))
            (run_dir / "species").mkdir(exist_ok=True)
            print(f"[stage 2] species needed: {total} across specs -> {len(requests)} unique -> "
                  f"{len(hits)} already in DB, {len(missing)} to compute")
            futures = {
                j: run_cmd(shlex.quote(str(run_dir / "species" / f"species_{j}")),
                           shlex.join(["python", "-m", "agents.compute_species",
                                       "--requests", str(req_file), "--index", str(j)]),
                           stdout=str(run_dir / "species" / f"species_{j}.log"),
                           stderr=str(run_dir / "species" / f"species_{j}.log"))
                for j in range(len(missing))
            }
            species_status = wait_all(futures, "species")
            # A species task killed mid-compute (worker loss, walltime) leaves its
            # claim 'running'; fail it now so assemble tasks don't wait on it.
            released = db.release_run_claims(run_dir.name, "claim abandoned during compute stage")
            summary["stages"]["compute"] = {
                "species_requested": total, "species_unique": len(requests), "db_hits": len(hits),
                "tasks": len(missing),
                "tasks_ok": sum(v == "COMPLETED" for v in species_status.values()),
                "claims_released": released,
                "specs_not_deduplicated": [spec_ids[k] for k in skipped],
                "wall_s": round(time.time() - t0, 1),
                "species": [{"smiles": r["smiles"], "charge": r["charge"], "multiplicity": r["multiplicity"],
                             "method": f"{r['xc']}/{r['basis']}", "needed_by": r["needed_by"],
                             "status": species_status.get(j, "?")}
                            for j, r in enumerate(missing)]
                           + [{"smiles": r["smiles"], "charge": r["charge"], "multiplicity": r["multiplicity"],
                               "method": f"{r['xc']}/{r['basis']}", "needed_by": r["needed_by"],
                               "status": "DB_HIT"} for r in hits],
            }

            # Stage 3: assemble every planned campaign from the database
            t0 = time.time()
            futures = {i: campaign(i, ["--from-spec", "spec.json"], str(run_dir / f"task_{i}.log"))
                       for i in specs}
            print(f"[stage 3] {len(futures)} assemble tasks submitted")
            results = wait_all(futures, "task")
            for i in tasks:
                results.setdefault(i, f"NO SPEC (planning {plan_status.get(i)})")
            summary["stages"]["assemble"] = {"tasks": len(futures), "wall_s": round(time.time() - t0, 1)}
    finally:
        if dfk is not None:
            dfk.cleanup()
        summary["total_wall_s"] = round(time.time() - t_start, 1)
        summary["task_status"] = results
        if db is not None:
            summary["db_after"] = db.stats()
        (run_dir / "run_summary.json").write_text(json.dumps(summary, indent=2, default=str))

    print("\nReport paths:")
    for i in range(1, args.total_tasks + 1):
        print(f"  {run_dir}/task_{i}/agentic_campaign_report.md  [{results.get(i, 'UNKNOWN')}]")
    print(f"Summary: {run_dir}/run_summary.json")


if __name__ == "__main__":
    main()
