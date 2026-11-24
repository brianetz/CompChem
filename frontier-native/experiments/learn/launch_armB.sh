#!/bin/bash
# Wait for the registration job, check the v2 model exists, then start Arm B (3 repeats + E2E), same queues as Arm A.
L=/lustre/orion/stf053/scratch/akhilesh/compchem/learn; REG=$1
while :; do s=$(sacct -nXj $REG -o State | head -1 | tr -d ' '); case $s in COMPLETED) break;; FAILED|CANCELLED*|TIMEOUT|OUT_OF_ME*) echo "register $REG $s"; exit 1;; esac; sleep 120; done
M=/lustre/orion/stf053/scratch/akhilesh/ollama-models/manifests/registry.ollama.ai/library/nemotron-3-nano-learn-v2/latest
[ -f $M ] || { echo "model manifest missing: $M"; exit 1; }
echo "=== launching Arm B $(date -Is)"; cd $L
setsid nohup bash arm_chain.sh B 1 nemotron-3-nano-learn-v2:latest on on service develop > chains/armB_r1.log 2>&1 < /dev/null &
sleep 3; setsid nohup bash arm_chain.sh B 2 nemotron-3-nano-learn-v2:latest on on batch normal > chains/armB_r2.log 2>&1 < /dev/null &
sleep 3; setsid nohup bash arm_chain.sh B 3 nemotron-3-nano-learn-v2:latest on on batch normal > chains/armB_r3.log 2>&1 < /dev/null &
sleep 3; setsid nohup bash run_e2e.sh armB_e2e nemotron-3-nano-learn-v2:latest on on - - batch normal > chains/armB_e2e.log 2>&1 < /dev/null &
wait
