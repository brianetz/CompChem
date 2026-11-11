#!/bin/bash
# nemotron-3-nano on plan board v2 (commit 3a64568), retrieval off: the stage 7 and 8 runs
# (B1v2, B2v2, NE1/NE2 with Llama 3.1 405B) repeated with the lighter model.
H=${COMPCHEM_HUB:-/lustre/orion/stf053/scratch/akhilesh/compchem}
SRC=$H/code/src
TOLUENE="What's the benzylic C-H bond dissociation energy of toluene in kcal/mol, using a cheap xtb prescreen and a DFT refinement?"
METHANOL="What's the O-H bond dissociation energy of methanol (CH3OH) in kcal/mol, using a cheap xtb prescreen and a DFT refinement?"

run() {  # run <label> <db> <question args...>
  local label=$1 db=$2; shift 2
  echo "=== $label start $(date -Is)"
  cd $H/driver
  $H/venv/bin/python submit_campaigns_parsl.py "$@" \
    --workers-per-node 8 --nodes 1 --time-limit-minutes 120 --partition service --qos develop --account STF053 \
    --root-dir $H/runs/campaigns --repo-src $SRC --venv-path $H/venv \
    --backend ollama --ollama-model nemotron-3-nano:latest --ollama-models-dir $H/models \
    --staged --no-rag --plan-board --plan-board-deadline-min 45 --plan-critic-skill $SRC/bde-plan-critic.md \
    --db $H/db/$db > $H/runs/campaigns/driver_$label.log 2>&1 < /dev/null
  echo "=== $label exit $? $(date -Is) $(grep -a 'Run directory' $H/runs/campaigns/driver_$label.log)"
}

run N7_toluene_nemo  pre_N7_toluene_nemo.sqlite  --question "$TOLUENE"  --n-tasks 8
run N7_methanol_nemo pre_N7_methanol_nemo.sqlite --question "$METHANOL" --n-tasks 8
run N8a_nemo pre_N8a_nemo.sqlite --questions-file $H/inputs/new_elements_20261021.txt
run N8b_nemo pre_N8b_nemo.sqlite --questions-file $H/inputs/new_elements_20261021.txt
echo "=== ALL DONE $(date -Is)"
