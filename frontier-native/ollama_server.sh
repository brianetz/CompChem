#!/bin/bash
# Native-sbatch batch script for the ONE shared Ollama server that all
# campaign tasks (job_template.sh) talk to over the network, instead of each
# task starting its own private `ollama serve`. Submitted by
# submit_campaigns.sh before the tasks.
#
# Ollama is stateless between requests - every `ollama.chat` call carries its
# client's full message history - so tasks sharing this server never see each
# other's context. OLLAMA_NUM_PARALLEL lets their requests run concurrently
# instead of queueing.
#
# Handshake with the tasks, all through files in $ROOT_DIR:
#   ollama_host   written once the server is up AND the model is confirmed
#                 present; contains host:port for tasks to export as OLLAMA_HOST
#   ollama_failed written instead if the server/model check fails, so tasks
#                 fail fast instead of waiting out their timeout
#   done_<jobid>  touched by each task on exit; the server shuts down once all
#                 $N_TASKS are present (or at its own time limit)

set -euo pipefail

: "${OLLAMA_BIN:=$HOME/software/ollama/bin/ollama}"
: "${OLLAMA_MODEL:=nemotron-3-nano:latest}"
: "${OLLAMA_MODELS:=}"
: "${OLLAMA_PORT:=11434}"
: "${ROOT_DIR:?ROOT_DIR must point at the shared campaign directory}"
: "${N_TASKS:?N_TASKS must be set to the number of client tasks}"

SERVER_DIR="$ROOT_DIR/ollama_server_${SLURM_JOB_ID}"
mkdir -p "$SERVER_DIR"
cd "$SERVER_DIR"

fail() {
  echo "$1" >&2
  echo "$1" >"$ROOT_DIR/ollama_failed"
  exit 1
}

if [ -n "$OLLAMA_MODELS" ]; then
  export OLLAMA_MODELS
fi
export OLLAMA_NUM_PARALLEL="${OLLAMA_NUM_PARALLEL:-$N_TASKS}"
export OLLAMA_HOST="0.0.0.0:$OLLAMA_PORT"

"$OLLAMA_BIN" serve >ollama_serve.log 2>&1 &
OLLAMA_PID=$!
trap 'kill "$OLLAMA_PID" 2>/dev/null || true' EXIT

ready=0
for _ in $(seq 1 60); do
  if "$OLLAMA_BIN" list >/dev/null 2>&1; then
    ready=1
    break
  fi
  sleep 1
done
[ "$ready" -eq 1 ] || fail "ollama serve did not come up on $(hostname); see $SERVER_DIR/ollama_serve.log"

# Compute nodes have no internet, so the model must already be pulled - catch
# that here once, rather than as a 404 inside every task.
"$OLLAMA_BIN" show "$OLLAMA_MODEL" >/dev/null 2>&1 \
  || fail "model '$OLLAMA_MODEL' not found on $(hostname) - run \`ollama pull $OLLAMA_MODEL\` on a login node first"

echo "$(hostname):$OLLAMA_PORT" >"$ROOT_DIR/ollama_host"
echo "Ollama server ready at $(hostname):$OLLAMA_PORT (OLLAMA_NUM_PARALLEL=$OLLAMA_NUM_PARALLEL)"

while kill -0 "$OLLAMA_PID" 2>/dev/null; do
  n_done=$(find "$ROOT_DIR" -maxdepth 1 -name 'done_*' | wc -l)
  if [ "$n_done" -ge "$N_TASKS" ]; then
    echo "All $N_TASKS tasks finished; shutting down."
    exit 0
  fi
  sleep 15
done
fail "ollama serve exited unexpectedly; see $SERVER_DIR/ollama_serve.log"
