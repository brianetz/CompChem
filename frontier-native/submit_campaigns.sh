#!/bin/bash
# Submit N parallel agentic-campaign tasks on Frontier via native `sbatch`,
# run directly from a Frontier login node (no S3M, no token, no HTTP).
#
# This is the native-Slurm counterpart to deploy/s3m-frontier/submit_campaigns.py.
# Key differences, all upsides of running natively instead of through S3M:
#   - This script runs ON Frontier, so it CAN create --root-dir itself
#     (the S3M version's local machine has no filesystem access at all and
#     requires you to pre-create it separately).
#   - `sbatch` always writes a real stdout/stderr log file
#     ($ROOT_DIR/slurm-<job_id>.out) - the S3M/Defiant debugging session
#     found cases where submitted jobs failed in ~1s with NO log file ever
#     created, which made root-causing much harder; that failure mode is
#     specific to the REST submission path and does not apply here.
#   - No `environment` field is required (unlike Slurm's REST API) - sbatch
#     naturally inherits your login shell's environment, plus we pass
#     `--export=ALL` explicitly so the exported QUESTION/REPO_SRC/etc. below
#     reach the job.
#
# QUESTION is passed via an exported shell variable, not Slurm's
# `--export=KEY=VALUE,...` comma-list syntax, because a real campaign
# question will very likely contain commas (e.g. "a cheap xtb prescreen,
# then DFT refinement"), which would corrupt that comma-separated form.
# `--export=ALL` instead forwards the whole environment, commas and all.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
JOB_TEMPLATE="$SCRIPT_DIR/job_template.sh"

# ---- defaults ----
ACCOUNT="STF053"
TIME_LIMIT_MINUTES=60
JOB_NAME_PREFIX="agentic-campaign-frontier"
LLM_BACKEND="ollama"
OLLAMA_MODEL="nemotron-3-nano:latest"
OLLAMA_BIN=""
OLLAMA_MODELS=""
CLAUDE_MODEL=""
POLL=1
POLL_INTERVAL_S=30

QUESTION=""
N_TASKS=""
ROOT_DIR=""
REPO_SRC=""
VENV_PATH=""
PARTITION=""

usage() {
  cat <<'EOF'
Usage: submit_campaigns.sh --question "..." --n-tasks N --root-dir DIR \
         --repo-src DIR --venv-path DIR --partition NAME [options]

Required:
  --question TEXT          Campaign question (quote it; commas are fine)
  --n-tasks N               Number of parallel tasks to submit
  --root-dir DIR            Shared output directory (created if missing -
                             this script runs ON Frontier, so it can do this)
  --repo-src DIR            Absolute path to CompChem/src on Frontier
  --venv-path DIR           Absolute path to the Frontier venv to activate
  --partition NAME          Slurm partition (confirm with `sinfo` first)

Options:
  --account NAME            Slurm account (default: STF053)
  --time-limit-minutes N    Per-task wall time limit (default: 60)
  --job-name-prefix NAME    Job name prefix (default: agentic-campaign-frontier)
  --backend NAME            ollama | claude_cli | anthropic (default: ollama)
  --ollama-model NAME       (default: nemotron-3-nano:latest)
  --ollama-bin PATH         Override job_template.sh's default $OLLAMA_BIN
  --ollama-models-dir DIR   Shared pre-pulled Ollama model cache
  --claude-model NAME       Only used with --backend claude_cli
  --poll-interval-s N       Seconds between sacct polls (default: 30)
  --no-poll                 Submit and exit without polling
  -h, --help                Show this help
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --question) QUESTION="$2"; shift 2 ;;
    --n-tasks) N_TASKS="$2"; shift 2 ;;
    --root-dir) ROOT_DIR="$2"; shift 2 ;;
    --repo-src) REPO_SRC="$2"; shift 2 ;;
    --venv-path) VENV_PATH="$2"; shift 2 ;;
    --partition) PARTITION="$2"; shift 2 ;;
    --account) ACCOUNT="$2"; shift 2 ;;
    --time-limit-minutes) TIME_LIMIT_MINUTES="$2"; shift 2 ;;
    --job-name-prefix) JOB_NAME_PREFIX="$2"; shift 2 ;;
    --backend) LLM_BACKEND="$2"; shift 2 ;;
    --ollama-model) OLLAMA_MODEL="$2"; shift 2 ;;
    --ollama-bin) OLLAMA_BIN="$2"; shift 2 ;;
    --ollama-models-dir) OLLAMA_MODELS="$2"; shift 2 ;;
    --claude-model) CLAUDE_MODEL="$2"; shift 2 ;;
    --poll-interval-s) POLL_INTERVAL_S="$2"; shift 2 ;;
    --no-poll) POLL=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

for name in QUESTION N_TASKS ROOT_DIR REPO_SRC VENV_PATH PARTITION; do
  if [ -z "${!name}" ]; then
    echo "Missing required argument: --$(echo "$name" | tr '[:upper:]_' '[:lower:]-')" >&2
    usage >&2
    exit 1
  fi
done

case "$LLM_BACKEND" in
  ollama|claude_cli|anthropic) ;;
  *) echo "Unknown --backend: $LLM_BACKEND (must be ollama|claude_cli|anthropic)" >&2; exit 1 ;;
esac

mkdir -p "$ROOT_DIR"

export LLM_BACKEND REPO_SRC ROOT_DIR QUESTION VENV_PATH
[ -n "$OLLAMA_MODEL" ] && export OLLAMA_MODEL
[ -n "$OLLAMA_BIN" ] && export OLLAMA_BIN
[ -n "$OLLAMA_MODELS" ] && export OLLAMA_MODELS
[ -n "$CLAUDE_MODEL" ] && export CLAUDE_MODEL

# ollama backend: ONE shared server job (ollama_server.sh) serves every task,
# instead of each task running its own `ollama serve`. Tasks are held until it
# has started (--dependency=after) and then wait for its ready file.
SERVER_JOB_ID=""
DEPENDENCY_ARGS=()
if [ "$LLM_BACKEND" = "ollama" ]; then
  # Clear handshake files left over from a previous run in the same ROOT_DIR.
  rm -f "$ROOT_DIR/ollama_host" "$ROOT_DIR/ollama_failed" "$ROOT_DIR"/done_*
  export N_TASKS
  SERVER_JOB_ID=$(sbatch --parsable \
    --account="$ACCOUNT" \
    --partition="$PARTITION" \
    --time="$((TIME_LIMIT_MINUTES + 20))" \
    --nodes=1 \
    --ntasks=1 \
    --job-name="${JOB_NAME_PREFIX}-ollama" \
    --output="$ROOT_DIR/slurm-%j.out" \
    --export=ALL \
    "$SCRIPT_DIR/ollama_server.sh")
  echo "Submitted shared Ollama server -> job $SERVER_JOB_ID"
  DEPENDENCY_ARGS=(--dependency="after:$SERVER_JOB_ID")
fi

job_ids=()
for i in $(seq 1 "$N_TASKS"); do
  job_id=$(sbatch --parsable \
    "${DEPENDENCY_ARGS[@]}" \
    --account="$ACCOUNT" \
    --partition="$PARTITION" \
    --time="$TIME_LIMIT_MINUTES" \
    --nodes=1 \
    --ntasks=1 \
    --job-name="${JOB_NAME_PREFIX}-${i}" \
    --output="$ROOT_DIR/slurm-%j.out" \
    --export=ALL \
    "$JOB_TEMPLATE")
  echo "Submitted task $i -> job $job_id"
  job_ids+=("$job_id")
done

echo
echo "Submitted ${#job_ids[@]} jobs: ${job_ids[*]}"
for jid in "${job_ids[@]}"; do
  echo "  job $jid -> $ROOT_DIR/task_${jid}/agentic_campaign_report.md"
done

if [ "$POLL" -eq 0 ]; then
  [ -n "$SERVER_JOB_ID" ] && echo "Shared Ollama server job $SERVER_JOB_ID stops after all tasks finish (or: scancel $SERVER_JOB_ID)"
  exit 0
fi

echo
echo "Polling (Ctrl-C to stop polling without cancelling jobs)..."
declare -A states
for jid in "${job_ids[@]}"; do states["$jid"]="PENDING"; done

is_terminal() {
  case "$1" in
    COMPLETED|FAILED|CANCELLED|TIMEOUT|NODE_FAIL|OUT_OF_MEMORY|DEADLINE|BOOT_FAIL) return 0 ;;
    *) return 1 ;;
  esac
}

while true; do
  all_done=1
  for jid in "${job_ids[@]}"; do
    if is_terminal "${states[$jid]}"; then
      continue
    fi
    all_done=0
    # sacct retains history after a job leaves squeue, unlike squeue itself -
    # use it so this still works once a fast-failing job has already finished.
    state=$(sacct -j "$jid" --format=State --noheader --parsable2 2>/dev/null | head -n1 | tr -d ' ')
    [ -n "$state" ] && states["$jid"]="$state"
    echo "[poll] job $jid: ${states[$jid]}"
  done
  if [ "$all_done" -eq 1 ]; then
    break
  fi
  sleep "$POLL_INTERVAL_S"
done

if [ -n "$SERVER_JOB_ID" ]; then
  # The server exits on its own once every task has touched its done_ file;
  # this just covers tasks that died without running their EXIT trap.
  scancel "$SERVER_JOB_ID" 2>/dev/null || true
fi

echo
echo "Final states:"
for jid in "${job_ids[@]}"; do
  echo "  job $jid: ${states[$jid]}"
done

echo
echo "Report paths (check directly - you're already on this filesystem):"
for jid in "${job_ids[@]}"; do
  echo "  $ROOT_DIR/task_${jid}/agentic_campaign_report.md"
done
