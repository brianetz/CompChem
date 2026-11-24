#!/bin/bash
#SBATCH -A STF053 -p batch -q normal -N 1 -t 02:00:00 -J n3sft2 -o /lustre/orion/stf053/scratch/akhilesh/compchem/learn/train/train_v2_sft_%j.out
# Round 2, stage 1: episodes v2 (v1 training runs + the col_e rerun) -> data_v2 -> SFT-only adapter.
set -e
source /lustre/orion/stf053/proj-shared/akhilesh/venvs/gpus_notebook/bin/activate
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 PYTHONPATH=/ccs/proj/stf053/akhilesh/CompChem/src
export PYTORCH_HIP_ALLOC_CONF=expandable_segments:True
export MIOPEN_USER_DB_PATH=/tmp/miopen-$SLURM_JOB_ID MIOPEN_CUSTOM_CACHE_DIR=/tmp/miopen-$SLURM_JOB_ID; mkdir -p /tmp/miopen-$SLURM_JOB_ID
CPY=/lustre/orion/stf053/scratch/akhilesh/compchem/venv/bin/python
$CPY -m learn.episodes /lustre/orion/stf053/scratch/akhilesh/compchem/learn/runs/run_20261122_114551 /lustre/orion/stf053/scratch/akhilesh/compchem/learn/runs/run_20261122_120236 /lustre/orion/stf053/scratch/akhilesh/compchem/learn/runs/run_20261122_133036 /lustre/orion/stf053/scratch/akhilesh/compchem/learn/runs/run_20261122_152308 \
  /lustre/orion/stf053/scratch/akhilesh/compchem/learn/runs/run_20261122_120241 /lustre/orion/stf053/scratch/akhilesh/compchem/learn/runs/run_20261122_222828 /lustre/orion/stf053/scratch/akhilesh/compchem/learn/runs/run_20261122_233424 /lustre/orion/stf053/scratch/akhilesh/compchem/learn/runs/run_20261123_005148 /lustre/orion/stf053/scratch/akhilesh/compchem/learn/runs/run_20261123_014931 \
  --tasks /lustre/orion/stf053/scratch/akhilesh/compchem/learn/tasks.json --out /lustre/orion/stf053/scratch/akhilesh/compchem/learn/episodes_train_v2.jsonl
$CPY -m learn.build_dataset /lustre/orion/stf053/scratch/akhilesh/compchem/learn/episodes_train_v2.jsonl --tasks /lustre/orion/stf053/scratch/akhilesh/compchem/learn/tasks.json --out-dir /lustre/orion/stf053/scratch/akhilesh/compchem/learn/data_v2
cat /lustre/orion/stf053/scratch/akhilesh/compchem/learn/data_v2/dataset_stats.json
python -m learn.train_lora --base /lustre/orion/stf053/proj-shared/akhilesh/hf_cache/NVIDIA-Nemotron-3-Nano-30B-A3B-BF16 \
  --sft /lustre/orion/stf053/scratch/akhilesh/compchem/learn/data_v2/sft.jsonl --out /lustre/orion/stf053/scratch/akhilesh/compchem/learn/train/adapter_v2_sft \
  --max-len 7168 --max-gib 20 --sft-epochs 1 --accum 8 --lr 1e-4
