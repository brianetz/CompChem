#!/bin/bash
# Native-sbatch batch script for one agentic-campaign task on FRONTIER.
#
# Forked from deploy/s3m-frontier/job_template.sh, with S3M-specific framing
# removed - this is submitted with plain `sbatch` from a Frontier login
# node (see submit_campaigns.sh), not via the S3M REST API. Account,
# partition, and time limit are passed as `sbatch` CLI flags by the wrapper
# script, not baked in here as #SBATCH directives, so they stay easy to
# change per run without editing this file.
#
# $ROOT_DIR is ONE shared directory, created directly by submit_campaigns.sh
# (which, unlike the S3M version, runs ON Frontier and so CAN touch this
# filesystem itself - no more "create this yourself first" step). Each job
# then makes and cd's into its own subdirectory named by its own
# $SLURM_JOB_ID (always unique, known only once Slurm assigns it), which is
# what keeps each parallel task's agentic_campaign_report.md,
# report_structures/, and flowcept_buffer.jsonl from colliding with any
# other task's.

set -euo pipefail

: "${LLM_BACKEND:=ollama}"            # ollama (default for HPC) | claude_cli | anthropic
# $HOME/software/ollama/bin/ollama was the ORIGINAL confirmed Ollama install
# path on Frontier from before the Defiant work in this repo's history -
# kept as the default, but confirm it still exists (`ls -l` it) before
# relying on it, or pass OLLAMA_BIN explicitly via submit_campaigns.sh.
: "${OLLAMA_BIN:=$HOME/software/ollama/bin/ollama}"
: "${OLLAMA_MODEL:=nemotron-3-nano:latest}"
: "${CLAUDE_MODEL:=}"                 # only used when LLM_BACKEND=claude_cli
: "${VENV_PATH:?VENV_PATH must be set - build a venv on Frontier's scratch filesystem first, see README}"
: "${REPO_SRC:?REPO_SRC must point at the CompChem src/ directory}"
: "${ROOT_DIR:?ROOT_DIR must point at an EXISTING shared directory on Frontier's filesystem}"
: "${OLLAMA_MODELS:=}"                # shared pre-pulled model cache dir, ollama backend only
: "${QUESTION:?QUESTION must be set to the campaign question}"

source "$VENV_PATH/bin/activate"
export PYTHONPATH="$REPO_SRC:${PYTHONPATH:-}"
TASK_DIR="$ROOT_DIR/task_${SLURM_JOB_ID}"
mkdir -p "$TASK_DIR"
cd "$TASK_DIR"

case "$LLM_BACKEND" in
  ollama)
    # Talk to the ONE shared server started by ollama_server.sh rather than
    # running a private `ollama serve` here. Tell the server we're done on
    # any exit, so it can shut down once every task has finished.
    trap 'touch "$ROOT_DIR/done_${SLURM_JOB_ID}"' EXIT
    : "${OLLAMA_WAIT_S:=900}"
    waited=0
    while [ ! -f "$ROOT_DIR/ollama_host" ]; do
      if [ -f "$ROOT_DIR/ollama_failed" ]; then
        echo "Shared Ollama server failed: $(cat "$ROOT_DIR/ollama_failed")" >&2
        exit 1
      fi
      if [ "$waited" -ge "$OLLAMA_WAIT_S" ]; then
        echo "Shared Ollama server not ready after ${OLLAMA_WAIT_S}s" >&2
        exit 1
      fi
      sleep 5
      waited=$((waited + 5))
    done
    # The ollama Python client (llm/planner_llm.py) reads OLLAMA_HOST.
    export OLLAMA_HOST="$(cat "$ROOT_DIR/ollama_host")"
    echo "Using shared Ollama server at $OLLAMA_HOST"
    BACKEND_ARGS=(--ollama --model "$OLLAMA_MODEL")
    ;;
  claude_cli)
    BACKEND_ARGS=(--claude-cli)
    if [ -n "$CLAUDE_MODEL" ]; then
      BACKEND_ARGS+=(--model "$CLAUDE_MODEL")
    fi
    ;;
  anthropic)
    BACKEND_ARGS=()
    ;;
  *)
    echo "Unknown LLM_BACKEND: $LLM_BACKEND" >&2
    exit 1
    ;;
esac

python -m agents.run_agentic_campaign \
  --question "$QUESTION" \
  "${BACKEND_ARGS[@]}" \
  --auto-approve
