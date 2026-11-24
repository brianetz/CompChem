#!/bin/bash
# Planning-only campaigns over task_gen tasks, graded by the oracle.
# usage: run_learn.sh <label> <model> <servers> <split> <offset> <limit> <partition> <qos> [extra driver args...]
H=${COMPCHEM_HUB:-/lustre/orion/stf053/scratch/akhilesh/compchem}; L=$H/learn; SRC=$H/code/src
LABEL=$1 MODEL=$2 SERVERS=$3 SPLIT=$4 OFFSET=$5 LIMIT=$6 PART=$7 QOS=$8; shift 8
echo "=== $LABEL start $(date -Is) $MODEL x$SERVERS $SPLIT[$OFFSET:+$LIMIT] on $PART/$QOS $*"
cd $H/driver
$H/venv/bin/python submit_campaigns_parsl.py --tasks $L/tasks.json --split $SPLIT --task-offset $OFFSET --task-limit $LIMIT \
  --workers-per-node 8 --nodes 1 --time-limit-minutes 120 --partition $PART --qos $QOS --account STF053 \
  --root-dir $L/runs --repo-src $SRC --venv-path $H/venv \
  --backend ollama --ollama-model $MODEL --ollama-models-dir $H/models \
  --ollama-servers $SERVERS --ollama-num-parallel 1 --ollama-context-length 8192 --ollama-load-timeout 30m --ollama-no-mmap \
  --staged --db $L/learn.sqlite --no-rag --plan-board --plan-board-deadline-min 30 --plan-critic-skill $SRC/bde-plan-critic.md \
  --plan-stage-only --trace "$@" > $L/runs/driver_$LABEL.log 2>&1 < /dev/null
echo "=== $LABEL exit $? $(date -Is) $(grep -a 'Run directory' $L/runs/driver_$LABEL.log)"
