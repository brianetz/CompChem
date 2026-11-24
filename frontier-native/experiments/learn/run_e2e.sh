#!/bin/bash
# Fair-comparison END-TO-END run: plan -> DFT compute -> assemble on the 17-task subset (tasks_e2e.json),
# first-try plans only (no oracle hint), each arm from a fresh copy of the same frozen DB.
# Everything except MODEL / RAG / CRITIC is fixed here so both arms run the same protocol.
# usage: run_e2e.sh <label> <model> <rag on|off> <critic on|off> - - <partition> <qos> [extra...]
H=${COMPCHEM_HUB:-/lustre/orion/stf053/scratch/akhilesh/compchem}; L=$H/learn; SRC=$H/code/src
LABEL=$1 MODEL=$2 RAG=$3 CRITIC=$4 OFFSET=$5 LIMIT=$6 PART=$7 QOS=$8; shift 8
FROZEN=$H/db/rag_eval_v1.sqlite            # sha256 c339a259e96365e5..., read-only
DB=$L/dbs/$LABEL.sqlite; cp $FROZEN $DB; chmod u+w $DB   # each run gets its own fresh copy
SW=(); [ "$RAG" = off ] && SW+=(--no-rag); [ "$CRITIC" = off ] && SW+=(--no-plan-critic)
echo "=== $LABEL start $(date -Is) $MODEL rag=$RAG critic=$CRITIC e2e17 on $PART/$QOS ${SW[*]} $*"
cd $H/driver
$H/venv/bin/python submit_campaigns_parsl.py --tasks $L/tasks_e2e.json --split eval \
  --workers-per-node 8 --nodes 1 --time-limit-minutes 120 --partition $PART --qos $QOS --account STF053 \
  --root-dir $L/runs --repo-src $SRC --venv-path $H/venv \
  --backend ollama --ollama-model $MODEL --ollama-models-dir $H/models \
  --ollama-servers 8 --ollama-num-parallel 1 --ollama-context-length 8192 --ollama-load-timeout 30m --ollama-no-mmap \
  --staged --db $DB --plan-board --plan-board-deadline-min 30 --plan-critic-skill $SRC/bde-plan-critic.md \
  --trace "${SW[@]}" "$@" > $L/runs/driver_$LABEL.log 2>&1 < /dev/null
echo "=== $LABEL exit $? $(date -Is) $(grep -a 'Run directory' $L/runs/driver_$LABEL.log)"
