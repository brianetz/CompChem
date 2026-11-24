#!/bin/bash
# One repeat of one arm = the 65 eval tasks in three parts, run back to back on one queue.
# usage: arm_chain.sh <arm A|B> <repeat> <model> <rag> <critic> <partition> <qos> [extra...]
L=/lustre/orion/stf053/scratch/akhilesh/compchem/learn
ARM=$1 REP=$2 MODEL=$3 RAG=$4 CRITIC=$5 PART=$6 QOS=$7; shift 7
for p in "1 0 22" "2 22 22" "3 44 21"; do set -- $p "$@"; n=$1 off=$2 lim=$3; shift 3
  bash $L/run_arm.sh arm${ARM}_r${REP}_p$n $MODEL $RAG $CRITIC $off $lim $PART $QOS "$@"; done
echo "=== CHAIN arm$ARM r$REP DONE $(date -Is)"
