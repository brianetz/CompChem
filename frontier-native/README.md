# Native-sbatch deployment: parallel agentic campaigns on Frontier

No S3M, no REST API, no token. Run entirely from a Frontier **login node**
you're SSH'd into, using plain `sbatch`/`sacct`. Use this while sorting out
S3M/Frontier API issues (see `deploy/s3m-frontier/README.md` for the open
items there); switch back to the S3M version later if/when that path is
validated - the two are independent, kept in separate folders on purpose.

Copy this whole `deploy/frontier-native/` folder to Frontier (e.g. as part
of cloning the repo) and run `submit_campaigns.sh` from a login node.

## Why native instead of S3M

- `submit_campaigns.sh` runs ON Frontier, so it can create `--root-dir`
  itself - no more "create this yourself first" step.
- `sbatch` always writes a real `slurm-<job_id>.out` log file. The Defiant/
  S3M debugging session found jobs that failed in ~1s via the REST path
  with **no log file ever created**, which made the root cause very hard to
  pin down. That failure mode is specific to S3M's REST submission path and
  doesn't apply to native `sbatch`.
- No `environment` field is required the way Slurm's REST API demands it -
  `sbatch` naturally inherits your login shell's environment; the wrapper
  explicitly passes `--export=ALL` so the exported `QUESTION`/`REPO_SRC`/etc.
  reach the job.
- No token to manage, rotate, or accidentally leak.

## Prerequisites (one-time setup, via your SSH session on Frontier)

1. Clone/copy the repo onto Frontier's scratch filesystem:
   ```bash
   git clone <repo-url> /path/to/scratch/CompChem
   ```
2. Build a Python venv there (same package list validated on Defiant):
   ```bash
   python3 -m venv /path/to/scratch/venvs/CompChem
   source /path/to/scratch/venvs/CompChem/bin/activate
   pip install --upgrade pip
   pip install rdkit ase pyscf cclib pandas numpy scipy
   pip install "flowcept>=1.0.2" academy-py ollama anthropic "tblite[ase]"
   ```
   (Quote `pkg>=x.y.z` version pins - unquoted `>` is parsed as shell output
   redirection by zsh/bash.)
3. Confirm Ollama is installed and runnable:
   ```bash
   ls -l $HOME/software/ollama/bin/ollama   # the script's default path
   $HOME/software/ollama/bin/ollama --version
   ```
   If missing/stale, reinstall via the no-sudo tarball method (same as
   Defiant's - see `deploy/s3m/README.md` prerequisite #4), then pass
   `--ollama-bin`/`--ollama-models-dir` explicitly if you install it
   somewhere other than the default path.
4. Confirm the partition name and that your account has access to it:
   ```bash
   sinfo
   sacctmgr show assoc where user=$USER format=account,user,partition,maxjobs,maxsubmit
   ```
   Frontier's compute nodes are GPU nodes (4x MI250X each) - this pipeline
   is CPU-only and requests no GPUs, so each task occupies (and burns
   allocation hours on) a full GPU node. If Frontier has a CPU-only
   partition, prefer it; otherwise be aware of the cost before running many
   parallel tasks.

You do **not** need to pre-create `--root-dir` - `submit_campaigns.sh` does
this itself (`mkdir -p`), since it runs on the same filesystem.

## Usage

```bash
./submit_campaigns.sh \
  --question "What's the O-H bond dissociation energy for water in kcal/mol, using a cheap xtb prescreen and a DFT refinement?" \
  --n-tasks 8 \
  --partition <confirmed partition> \
  --root-dir /path/to/scratch/campaign_runs \
  --repo-src /path/to/scratch/CompChem/src \
  --venv-path /path/to/scratch/venvs/CompChem \
  --ollama-model nemotron-3-nano:latest
```

Run `./submit_campaigns.sh --help` for the full flag list. `--account`
defaults to `STF053`; everything else required has no default on purpose
(no path has been confirmed on Frontier this session - see
`deploy/s3m-frontier/README.md`'s open items, which apply here too for the
venv/Ollama paths, just not for the S3M-specific ones).

This submits `--n-tasks` separate `sbatch` jobs (one per task, each
`--nodes=1 --ntasks=1`), prints each job ID and its expected report path,
then polls `sacct` every `--poll-interval-s` (default 30s) until all reach
a terminal state. Pass `--no-poll` to submit and check back later with
`sacct -j <id> --format=JobID,State,ExitCode,Reason -p` directly.

To use `--claude-cli` instead of local Ollama (needs a login/compute node
with internet egress and a logged-in `claude` CLI):
```bash
  --backend claude_cli
```

## Files

- `job_template.sh` - the batch script run inside each `sbatch` job.
  Forked from `deploy/s3m-frontier/job_template.sh` with S3M-specific
  framing removed; the actual backend-selection/run logic is identical.
- `submit_campaigns.sh` - the native-Slurm counterpart to
  `deploy/s3m-frontier/submit_campaigns.py`: creates `--root-dir`, submits
  N `sbatch` jobs with `--export=ALL` (after exporting the needed env vars
  itself), polls via `sacct`, and prints each task's report path.
- `ollama_server.sh` - with the `ollama` backend, ONE shared `ollama serve`
  job that all tasks connect to via `OLLAMA_HOST` (with
  `OLLAMA_NUM_PARALLEL=N_TASKS`), instead of one server per task. It
  hands its address to the tasks through `$ROOT_DIR/ollama_host`, fails
  fast if the model isn't pulled, and exits once every task has touched
  `$ROOT_DIR/done_<jobid>`. Each request carries its own message history,
  so tasks don't share context.
- `submit_campaigns_parsl.py` - Parsl alternative to `submit_campaigns.sh`:
  one Slurm block of `--nodes` nodes running `--workers-per-node` tasks
  each (default 1 node x 8), instead of one node per task. Uses
  `get_frontier_config()` in `CompChem/src/parsl_config.py`. The shared
  `ollama serve` runs inside that same block (started by its
  `worker_init`, GPU for Ollama, CPUs for the tasks), so server and tasks
  always share one allocation; `ollama_server.sh` is only used by
  `submit_campaigns.sh`. Run with the venv's python, under tmux/nohup (the
  Parsl interchange lives in this process). Outputs go to
  `$ROOT_DIR/run_<timestamp>/task_<i>/`.

  One planner per question: pass `--question` once per question, or
  `--questions-file` (one per line), and every question is planned once,
  all in parallel. `--n-tasks K` (default 1) runs K copies of each question
  instead - only for measuring how repeatable the planner is.

  ```bash
  /path/to/venv/bin/python submit_campaigns_parsl.py \
    --questions-file questions.txt --staged --db /path/to/results.sqlite \
    --partition service --qos develop --time-limit-minutes 120 \
    --root-dir /path/to/campaign_runs --repo-src /path/to/CompChem/src \
    --venv-path /path/to/venv --ollama-model nemotron-3-nano:latest
  ```
