# Experiment scripts

Launchers and analysis for the BDE campaign experiments run on Frontier with
`../submit_campaigns_parsl.py`. They assume the hub layout on Lustre
(`$COMPCHEM_HUB`, default `/lustre/orion/stf053/scratch/akhilesh/compchem`):

```
$COMPCHEM_HUB/code    -> this repository
$COMPCHEM_HUB/venv    -> Python 3.11 venv with requirements installed
$COMPCHEM_HUB/driver  -> code/frontier-native
$COMPCHEM_HUB/models  -> shared Ollama model cache (OLLAMA_MODELS)
$COMPCHEM_HUB/db      -> results databases (one copy per experiment)
$COMPCHEM_HUB/inputs  -> question files (one question per line; copies in inputs/)
$COMPCHEM_HUB/runs/campaigns/run_<date>_<time>/   one folder per run
```

## launch/

| Script | What it runs |
|---|---|
| `run_exp.sh` | One run of K copies of a question (default water O–H, nemotron-3-nano). `QUESTION="..." NTASKS=8 run_exp.sh <label> <db or none> [driver flags]` |
| `board_v2_nano_chain.sh` | nemotron-3-nano on plan board v2, retrieval off: toluene × 8, methanol × 8, molecule set 1 twice |
| `board_v2_super_chain.sh` | The same runs with nemotron-3-super:120b on 4 Ollama servers per node. `board_v2_super_chain.sh <partition> <qos> S7T S7M S8a S8b SM` |
| `test_multi_ollama.sbatch` | 25-min check that 4 Ollama servers on one node each get their own 2 GCDs and answer |

Each run starts from its own copy of the results database (`--db`), so runs
being compared start from the same cached species.

Several Ollama servers per node: Slurm sets `ROCR_VISIBLE_DEVICES=0-7`, so a
server is restricted to its GCDs by replacing that list, not by adding
`HIP_VISIBLE_DEVICES` (which leaves all but the first server on the CPU).

## analysis/

Run with the repo's `src/` on `PYTHONPATH` and the hub venv's Python.

| Script | Output |
|---|---|
| `analyze_run.py <run_dir>` | Verdicts, BDEs, species computed vs reused, block time of one run |
| `classify_board.py <run_dir>...` | Per campaign: correct / wrong answer passed / failure reason; writes `<run_dir>/.board_scores.json` |
| `board_metrics.py <label>=<run_dir>...` | Correct, board rounds, BDE and method per campaign, planning and run time (after `classify_board.py`) |
| `classify_repeatability.py <run_dir>...` | Methanol O–H repeatability runs: same, wrong and failed answers |
| `stage_times.py <out.json>` | Time per run split into queue, startup, planning, plan critic, execution, result generation; from `run_summary.json`, Slurm, Ollama logs and Flowcept records |
| `make_time_section.py`, `make_error_section.py` | Charts and tables of time per stage and of failure reasons per run, from the JSON above |

"Correct" means the results critic gave PASS and the parent and fragments are
the requested bond's homolysis (`EXPECTED` in `classify_board.py`).
