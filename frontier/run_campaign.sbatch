#!/bin/bash
# CompChem agentic campaign on one Frontier compute node.
#
# Submit (from anywhere):
#   sbatch frontier/run_campaign.sbatch "What's the O-H BDE of water in kcal/mol?"
#   sbatch frontier/run_campaign.sbatch --smiles CCO        # rule-based, no LLM
# Extra campaign flags can follow, e.g. --model claude-haiku-4-5-20251001
#
#SBATCH -A stf053
#SBATCH -J compchem
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -q debug
#SBATCH -o /lustre/orion/stf053/scratch/%u/compchem-runs/%x-%j.out

set -euo pipefail

REPO=/ccs/proj/stf053/akhilesh/CompChem
VENV=/ccs/proj/stf053/akhilesh/venvs/compchem
RUN_DIR=/lustre/orion/stf053/scratch/$USER/compchem-runs/$SLURM_JOB_ID

if [[ $# -eq 0 ]]; then
    echo "usage: sbatch $0 \"<question>\" | --smiles <SMILES> [campaign flags]" >&2
    exit 2
fi

# Python env + code. PYTHONPATH replaces `cd src`, so outputs land in RUN_DIR.
source "$VENV/bin/activate"
export PYTHONPATH=$REPO/src
export PATH=$HOME/.local/bin:$PATH          # `claude` for --claude-cli

# Compute nodes only reach the internet through the OLCF proxy.
export http_proxy=http://proxy.ccs.ornl.gov:3128
export https_proxy=http://proxy.ccs.ornl.gov:3128
export no_proxy=localhost,127.0.0.1,.ccs.ornl.gov,.olcf.ornl.gov

# xtb and PySCF are OpenMP-threaded; give them the whole node.
export OMP_NUM_THREADS=56
export OMP_STACKSIZE=1G

mkdir -p "$RUN_DIR" && cd "$RUN_DIR"
echo "job $SLURM_JOB_ID on $(hostname), branch $(git -C "$REPO" rev-parse --abbrev-ref HEAD)@$(git -C "$REPO" rev-parse --short HEAD)"
echo "run dir: $RUN_DIR"

if [[ $1 == --smiles ]]; then
    shift; SMILES=$1; shift
    CAMPAIGN=("$SMILES" "$@")
else
    QUESTION=$1; shift
    CAMPAIGN=(--question "$QUESTION" --claude-cli "$@")
fi

# Batch jobs have no terminal, so the human approval gates must be skipped.
srun -n1 -c56 --cpu-bind=threads \
    python -m agents.run_agentic_campaign "${CAMPAIGN[@]}" --auto-approve

echo "done: $RUN_DIR/agentic_campaign_report.md"
