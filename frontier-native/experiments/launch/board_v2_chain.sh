#!/bin/bash
# Plan board v2 (commit 3a64568 planner code), retrieval off: the stage 7/8 runs with any Ollama model
# and any number of Ollama servers per node, each run from its own reset copy of the results DB.
# usage: board_v2_chain.sh <model> <servers> <db-tag> <partition> <qos> <run>...
#   run = 7T 7M 8a 8b M  (toluene x8, methanol x8, set 1, set 1 again, set 2)
#   DB copies are $H/db/pre_<N7_toluene|N7_methanol|N8a|N8b|NM>_<db-tag>.sqlite
H=${COMPCHEM_HUB:-/lustre/orion/stf053/scratch/akhilesh/compchem}; SRC=$H/code/src
TOLUENE="What's the benzylic C-H bond dissociation energy of toluene in kcal/mol, using a cheap xtb prescreen and a DFT refinement?"
METHANOL="What's the O-H bond dissociation energy of methanol (CH3OH) in kcal/mol, using a cheap xtb prescreen and a DFT refinement?"
MODEL=$1 SERVERS=$2 TAG=$3 PART=$4 QOS=$5; shift 5
run() {  # run <label> <db> <question args...>
  local label=$1 db=$2; shift 2
  echo "=== $label start $(date -Is) on $PART/$QOS ($MODEL, $SERVERS servers)"
  cd $H/driver
  $H/venv/bin/python submit_campaigns_parsl.py "$@" \
    --workers-per-node 8 --nodes 1 --time-limit-minutes 120 --partition $PART --qos $QOS --account STF053 \
    --root-dir $H/runs/campaigns --repo-src $SRC --venv-path $H/venv \
    --backend ollama --ollama-model $MODEL --ollama-models-dir $H/models \
    --ollama-servers $SERVERS --ollama-num-parallel $(( 8 / SERVERS )) --ollama-context-length 8192 \
    --ollama-load-timeout 30m --ollama-no-mmap \
    --staged --no-rag --plan-board --plan-board-deadline-min 45 --plan-critic-skill $SRC/bde-plan-critic.md \
    --db $H/db/$db > $H/runs/campaigns/driver_$label.log 2>&1 < /dev/null
  echo "=== $label exit $? $(date -Is) $(grep -a 'Run directory' $H/runs/campaigns/driver_$label.log)"
}
for r in "$@"; do case $r in
  7T) run 7T_$TAG pre_N7_toluene_$TAG.sqlite  --question "$TOLUENE"  --n-tasks 8 ;;
  7M) run 7M_$TAG pre_N7_methanol_$TAG.sqlite --question "$METHANOL" --n-tasks 8 ;;
  8a) run 8a_$TAG pre_N8a_$TAG.sqlite --questions-file $H/inputs/new_elements_20261021.txt ;;
  8b) run 8b_$TAG pre_N8b_$TAG.sqlite --questions-file $H/inputs/new_elements_20261021.txt ;;
  M)  run M_$TAG  pre_NM_$TAG.sqlite  --questions-file $H/inputs/new_molecules_20261028.txt ;;
esac; done
echo "=== ALL DONE $(date -Is)"
