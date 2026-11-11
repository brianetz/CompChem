#!/bin/bash
# Check that 4 Ollama servers on one node each get their own 2 GCDs (ROCR_VISIBLE_DEVICES) and answer.
#SBATCH -A STF053 -p batch -N 1 -t 00:25:00 -J super4test -o test_multi_ollama_%j.out
O=${OLLAMA_BIN:-~/software/ollama/bin/ollama}; D=${COMPCHEM_HUB:-/lustre/orion/stf053/scratch/akhilesh/compchem}/jobs/test_multi_ollama_$SLURM_JOB_ID; mkdir -p $D
export OLLAMA_MODELS=/lustre/orion/stf053/scratch/akhilesh/ollama-models
echo "slurm ROCR=$ROCR_VISIBLE_DEVICES HIP=$HIP_VISIBLE_DEVICES"
for s in 0 1 2 3; do g=$((2*s)),$((2*s+1))
  env -u HIP_VISIBLE_DEVICES ROCR_VISIBLE_DEVICES=$g OLLAMA_HOST=0.0.0.0:$((11434+s)) OLLAMA_NUM_PARALLEL=2 OLLAMA_CONTEXT_LENGTH=8192 OLLAMA_LOAD_TIMEOUT=30m $O serve > $D/serve_$s.log 2>&1 &
  SERVERS+=($!)
done
sleep 20
for s in 0 1 2 3; do
  ( t=$(date +%s); curl -s -m 1200 http://127.0.0.1:$((11434+s))/api/chat -d '{"model":"nemotron-3-super:120b","stream":false,"think":false,"options":{"use_mmap":false,"num_predict":20},"messages":[{"role":"user","content":"Say OK."}]}' | head -c 300; echo " <- server $s in $(( $(date +%s)-t )) s" ) &
  REQS+=($!)
done; wait "${REQS[@]}"
for s in 0 1 2 3; do echo "server $s: $(grep -a 'inference compute' $D/serve_$s.log | grep -c ROCm) ROCm GPUs; $(grep -a -o 'total_vram="[^"]*"' $D/serve_$s.log | head -1)"; done
rocm-smi --showmemuse 2>/dev/null | grep -i "GPU\[" | head -16
kill "${SERVERS[@]}"
