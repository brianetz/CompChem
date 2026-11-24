#!/bin/bash
#SBATCH -A STF053 -p batch -q normal -N 1 -t 00:45:00 -J n3reg2 -o /lustre/orion/stf053/scratch/akhilesh/compchem/learn/train/register_v2_%j.out
# Convert adapter_v2 to GGUF and register nemotron-3-nano-learn-v2 (base nemotron-3-nano:latest + ADAPTER), then smoke-test it.
set -e
source /lustre/orion/stf053/proj-shared/akhilesh/venvs/gpus_notebook/bin/activate
cd /lustre/orion/stf053/scratch/akhilesh/compchem/learn/tools/llama.cpp
python convert_lora_to_gguf.py --base /lustre/orion/stf053/proj-shared/akhilesh/hf_cache/NVIDIA-Nemotron-3-Nano-30B-A3B-BF16 \
  --outfile /lustre/orion/stf053/scratch/akhilesh/compchem/learn/train/adapter_v2/adapter.gguf --outtype f16 /lustre/orion/stf053/scratch/akhilesh/compchem/learn/train/adapter_v2 2>&1 | tail -5
printf 'FROM nemotron-3-nano:latest\nADAPTER /lustre/orion/stf053/scratch/akhilesh/compchem/learn/train/adapter_v2/adapter.gguf\n' > /lustre/orion/stf053/scratch/akhilesh/compchem/learn/train/adapter_v2/Modelfile
export OLLAMA_MODELS=/lustre/orion/stf053/scratch/akhilesh/ollama-models OLLAMA_HOST=127.0.0.1:11558 ROCR_VISIBLE_DEVICES=0
unset HIP_VISIBLE_DEVICES CUDA_VISIBLE_DEVICES
rsync -a --exclude 'cuda*' --exclude 'vulkan*' ~/software/ollama/ /tmp/ollama-$SLURM_JOB_ID/
O=/tmp/ollama-$SLURM_JOB_ID/bin/ollama
$O serve > /lustre/orion/stf053/scratch/akhilesh/compchem/learn/train/adapter_v2/create_serve.log 2>&1 & sleep 15
$O create nemotron-3-nano-learn-v2 -f /lustre/orion/stf053/scratch/akhilesh/compchem/learn/train/adapter_v2/Modelfile
$O run nemotron-3-nano-learn-v2 "Reply with the single word OK." --think=false 2>&1 | tail -3
$O list | grep nemotron
