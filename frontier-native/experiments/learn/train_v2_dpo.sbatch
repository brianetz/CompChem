#!/bin/bash
#SBATCH -A STF053 -p batch -q normal -N 1 -t 02:00:00 -J n3dpo2 -o /lustre/orion/stf053/scratch/akhilesh/compchem/learn/train/train_v2_dpo_%j.out
# Round 2, stage 2: DPO starting from the SFT adapter (reference = base, LoRA disabled), same settings as v1.
set -e
source /lustre/orion/stf053/proj-shared/akhilesh/venvs/gpus_notebook/bin/activate
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 PYTHONPATH=/ccs/proj/stf053/akhilesh/CompChem/src
export PYTORCH_HIP_ALLOC_CONF=expandable_segments:True
export MIOPEN_USER_DB_PATH=/tmp/miopen-$SLURM_JOB_ID MIOPEN_CUSTOM_CACHE_DIR=/tmp/miopen-$SLURM_JOB_ID; mkdir -p /tmp/miopen-$SLURM_JOB_ID
python -m learn.train_lora --base /lustre/orion/stf053/proj-shared/akhilesh/hf_cache/NVIDIA-Nemotron-3-Nano-30B-A3B-BF16 \
  --dpo /lustre/orion/stf053/scratch/akhilesh/compchem/learn/data_v2/dpo.jsonl --init-adapter /lustre/orion/stf053/scratch/akhilesh/compchem/learn/train/adapter_v2_sft --out /lustre/orion/stf053/scratch/akhilesh/compchem/learn/train/adapter_v2 \
  --max-len 7168 --max-gib 20 --dpo-epochs 1 --accum 8 --dpo-lr 2e-5 --beta 0.1 --sft-alpha 0.2
