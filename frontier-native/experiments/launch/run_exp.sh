#!/bin/bash
# usage: QUESTION="..." NTASKS=56 run_exp.sh <label> <db-path|none> [extra driver flags...]
H=${COMPCHEM_HUB:-/lustre/orion/stf053/scratch/akhilesh/compchem}
LABEL=$1; DB=$2; shift 2
DEFAULT_Q="What's the O-H bond dissociation energy for water in kcal/mol, using a cheap xtb prescreen and a DFT refinement?"
Q=${QUESTION:-$DEFAULT_Q}
cd $H/driver
DBARG=(); [ "$DB" != none ] && DBARG=(--db "$DB")
exec $H/venv/bin/python submit_campaigns_parsl.py \
  --question "$Q" \
  --n-tasks "${NTASKS:-56}" --workers-per-node 8 --nodes 1 --time-limit-minutes 60 \
  --partition service --qos develop --account STF053 \
  --root-dir $H/runs/campaigns \
  --repo-src $H/code/src \
  --venv-path $H/venv \
  --backend ollama --ollama-model nemotron-3-nano:latest \
  --ollama-models-dir $H/models \
  "${DBARG[@]}" "$@" > "$H/runs/campaigns/driver_${LABEL}.log" 2>&1
