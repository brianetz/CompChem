#!/bin/bash
# nemotron-3-super:120b (Q4_K_M, 120B total / 12B active) on plan board v2 (commit 3a64568), retrieval off:
# repeats of the nemotron-3-nano runs N7-T, N7-M, N8-a, N8-b and NM (2026-10-28), same DB starting state.
# usage: board_v2_super_chain.sh <partition> <qos> <run>...   run = S7T S7M S8a S8b SM
H=${COMPCHEM_HUB:-/lustre/orion/stf053/scratch/akhilesh/compchem}; SRC=$H/code/src
TOLUENE="What's the benzylic C-H bond dissociation energy of toluene in kcal/mol, using a cheap xtb prescreen and a DFT refinement?"
METHANOL="What's the O-H bond dissociation energy of methanol (CH3OH) in kcal/mol, using a cheap xtb prescreen and a DFT refinement?"
PART=$1 QOS=$2; shift 2
run() {  # run <label> <db> <question args...>
  local label=$1 db=$2; shift 2
  echo "=== $label start $(date -Is) on $PART/$QOS"
  cd $H/driver
  $H/venv/bin/python submit_campaigns_parsl.py "$@" \
    --workers-per-node 8 --nodes 1 --time-limit-minutes 120 --partition $PART --qos $QOS --account STF053 \
    --root-dir $H/runs/campaigns --repo-src $SRC --venv-path $H/venv \
    --backend ollama --ollama-model nemotron-3-super:120b --ollama-models-dir $H/models \
    --ollama-servers 4 --ollama-num-parallel 2 --ollama-context-length 8192 --ollama-load-timeout 30m --ollama-no-mmap \
    --staged --no-rag --plan-board --plan-board-deadline-min 45 --plan-critic-skill $SRC/bde-plan-critic.md \
    --db $H/db/$db > $H/runs/campaigns/driver_$label.log 2>&1 < /dev/null
  echo "=== $label exit $? $(date -Is) $(grep -a 'Run directory' $H/runs/campaigns/driver_$label.log)"
}
for r in "$@"; do case $r in
  S7T) run S7T_toluene_super  pre_N7_toluene_super.sqlite  --question "$TOLUENE"  --n-tasks 8 ;;
  S7M) run S7M_methanol_super pre_N7_methanol_super.sqlite --question "$METHANOL" --n-tasks 8 ;;
  S8a) run S8a_super pre_N8a_super.sqlite --questions-file $H/inputs/new_elements_20261021.txt ;;
  S8b) run S8b_super pre_N8b_super.sqlite --questions-file $H/inputs/new_elements_20261021.txt ;;
  SM)  run SM_super  pre_NM_super.sqlite  --questions-file $H/inputs/new_molecules_20261028.txt ;;
esac; done
echo "=== ALL DONE $(date -Is)"
