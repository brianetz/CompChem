#!/bin/bash
# Run several run_learn.sh runs one after another on one queue.
# usage: learn_chain.sh <partition> <qos> <model> <servers> "<label> <split> <offset> <limit> [extra...]" ...
L=/lustre/orion/stf053/scratch/akhilesh/compchem/learn
PART=$1 QOS=$2 MODEL=$3 SERVERS=$4; shift 4
for spec in "$@"; do
  set -- $spec
  label=$1 split=$2 off=$3 lim=$4; shift 4
  bash $L/run_learn.sh $label $MODEL $SERVERS $split $off $lim $PART $QOS "$@"
done
echo "=== CHAIN DONE $(date -Is)"
