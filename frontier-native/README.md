# Running BDE campaigns on Frontier

`submit_campaigns_parsl.py` runs a batch of agentic bond-dissociation-energy
(BDE) campaigns inside **one Slurm job** on Frontier, from a login node:

- the job's first node starts one or more **Ollama** servers on its GPUs, which
  serve the LLM planner, plan critic and supervisor;
- **Parsl** runs 8 campaign workers per node on the CPUs (7 cores each);
- with `--staged` every campaign is planned first, each unique species is then
  computed once (GFN2-xTB geometry + PySCF DFT energy), and every campaign is
  assembled from the stored energies;
- a shared **results database** (SQLite on Lustre) reuses species already
  computed by earlier runs, and can feed prior results back to the planner;
- **Flowcept** records every planner/critic call and every computation
  (provenance) next to the outputs.

## 1. One-time setup

Keep the code where it is backed up and the outputs on Lustre:

| What | Where (example) | Why |
|---|---|---|
| Repository | `/ccs/proj/<project>/<user>/CompChem` | backed up, not purged |
| Runs, databases, question files | `/lustre/orion/<project>/scratch/<user>/compchem/` | fast, shared with compute nodes; **purged after ~90 days unused** |
| Ollama models | `/lustre/orion/<project>/scratch/<user>/ollama-models` | 24–431 GB per model |

**Python environment** (Python 3.11; the system Python 3.6 has no RDKit wheels):

```bash
/opt/cray/pe/python/3.11.7/bin/python3.11 -m venv /ccs/proj/<project>/<user>/venvs/compchem
source /ccs/proj/<project>/<user>/venvs/compchem/bin/activate
pip install --upgrade pip
pip install rdkit ase pyscf cclib pandas numpy scipy parsl "flowcept>=1.0.2" academy-py ollama anthropic "tblite[ase]"
```

**Ollama**: install the release tarball natively under `~/software/ollama`
(no container: the Apptainer build crashed under ROCm). Then pull models from a
**login node** (compute nodes have no direct internet) into the shared model
directory, using a throw-away server on a spare port:

```bash
export OLLAMA_MODELS=/lustre/orion/<project>/scratch/<user>/ollama-models
OLLAMA_HOST=127.0.0.1:11555 ~/software/ollama/bin/ollama serve &
OLLAMA_HOST=127.0.0.1:11555 ~/software/ollama/bin/ollama pull nemotron-3-super:120b
OLLAMA_HOST=127.0.0.1:11555 ~/software/ollama/bin/ollama pull nemotron-3-nano:latest
kill %1
```

The results database needs no setup: pass a path with `--db` and it is
created on first use. Copy it before an experiment you want to compare
against a fixed starting state.

## 2. Run a batch of campaigns

Put one question per line in a file (blank lines and `#` comments are
ignored); see `experiments/inputs/` for two 8-molecule sets. Then, from a login
node, inside `tmux` or with `nohup` (the Parsl interchange lives in the driver
process, so it must stay alive until the run ends):

```bash
H=/lustre/orion/<project>/scratch/<user>/compchem
SRC=/ccs/proj/<project>/<user>/CompChem/src
VENV=/ccs/proj/<project>/<user>/venvs/compchem
cd /ccs/proj/<project>/<user>/CompChem/frontier-native

$VENV/bin/python submit_campaigns_parsl.py \
  --questions-file $H/inputs/new_elements_20261021.txt \
  --nodes 1 --workers-per-node 8 --time-limit-minutes 120 \
  --partition service --qos develop --account <PROJECT> \
  --root-dir $H/runs/campaigns --repo-src $SRC --venv-path $VENV \
  --backend ollama --ollama-model nemotron-3-super:120b --ollama-models-dir $H/models \
  --ollama-servers 4 --ollama-num-parallel 2 --ollama-context-length 8192 \
  --ollama-load-timeout 30m --ollama-no-mmap \
  --staged --db $H/db/results.sqlite \
  --plan-board --plan-board-deadline-min 45 --plan-critic-skill $SRC/bde-plan-critic.md
```

The driver prints the run folder (`Run directory: .../run_<date>_<time>`),
waits for the job, and writes `run_summary.json` at the end. To ask the same
question several times (to measure repeatability), use `--question "..."` with
`--n-tasks 8` instead of `--questions-file`.

### Choosing the model and serving

Ollama serves the nemotron models one request at a time per server, so with
one server every campaign queues behind the others. `--ollama-servers N`
starts N servers on the first node, each on 8/N GPU dies and its own port;
campaign *i* uses server *(i−1) mod N*. Measured on plan board v2, retrieval
off, 8 campaigns per run:

| Model (`--ollama-model`) | Size | Serving | Correct | Run of 8 campaigns |
|---|---|---|---|---|
| `nemotron-3-super:120b` | 87 GB, 2 dies per copy | `--ollama-servers 4` | 40 / 40 | 11–15 min |
| `nemotron-3-nano:latest` | 24 GB, 1 die | `--ollama-servers 8` | 39 / 40 | 12–17 min typical, up to 48 |
| `nemotron-3-nano:latest` | 24 GB, 1 die | `--ollama-servers 4` | 40 / 40 | 16–22 min typical, up to 38 |
| `nemotron-3-nano:latest` | 24 GB, 1 die | 1 server | 38 / 40 | 46–61 min |
| `llama3.1:405b-instruct-q8_0` | 431 GB, all 8 dies | 1 server, `--ollama-num-parallel 4` | 32 / 32 | 25–30 min |

`nemotron-3-super` on 4 servers is the recommended default: the most accurate
and the most consistent run time. With 405B use `--ollama-no-mmap` and
`--ollama-load-timeout 30m` (loading 431 GB takes about 4 minutes).

### Main options

| Option | Effect |
|---|---|
| `--questions-file F` / `--question Q` | One campaign per question; `--n-tasks K` repeats each question K times |
| `--staged --db PATH` | Plan all, compute each unique missing species once, assemble from the database (recommended) |
| `--no-rag` | Don't give the planner prior results from the database. Retrieval on keeps the DFT method consistent between campaigns; off gives an independent test of the planner |
| `--plan-board` | Plan module by module (target, fragments, method, rationale), each checked by code and reviewed by the plan critic |
| `--plan-board-deadline-min M` | Planning budget per campaign (45 is what the experiments used) |
| `--plan-critic-skill PATH` | The plan critic's instructions (`src/bde-plan-critic.md`) |
| `--ollama-servers N` | Ollama servers per node (1, 2, 4 or 8) |
| `--partition/--qos` | `service`/`develop` (≤ 2 h; used for most runs) or `extended`/`normal` for a second batch at the same time |

`--help` lists everything.

## 3. What a run writes

```
run_<date>_<time>/
  run_summary.json            stage timings, species computed vs reused, per-campaign status
  plan_<i>.log                planning output of campaign i, including the board history
                              (each round's code check and critic verdict)
  task_<i>.log                assemble output of campaign i
  task_<i>/spec.json          the accepted plan (missing if planning failed)
  task_<i>/agentic_campaign_report.md   energies, BDE and QC verdict
  task_<i>/plan/flowcept_buffer.jsonl   provenance of planning: every planner, critic and
                                        supervisor call with its model, input and output
  task_<i>/flowcept_buffer.jsonl        provenance of assembly: database lookups, executor,
                                        results critic, reporter
  species/species_<j>/flowcept_buffer.jsonl   provenance of one species' xtb and DFT run
  ollama_serve.log | ollama_serve_<s>.log     model server logs (GPUs found, per-call timings)
  parsl/                      Parsl, interchange, worker and Slurm block logs
```

Flowcept writes a campaign's buffer when that campaign finishes, so a campaign
killed mid-planning has no provenance file; look in `parsl/` (interchange log,
block `.stderr`) and the Ollama logs instead.

To read a campaign's agent decisions:

```bash
cd run_<...>/task_1/plan && PYTHONPATH=$SRC $VENV/bin/python $SRC/verify_agent_provenance.py
```

## 4. Checking the results

The scripts in `experiments/analysis/` score and time runs (use the venv's
Python with `PYTHONPATH=$SRC`):

```bash
cd experiments/analysis
$VENV/bin/python classify_board.py  $H/runs/campaigns/run_<...>          # correct / wrong / failure reason per campaign
$VENV/bin/python board_metrics.py   label=$H/runs/campaigns/run_<...>    # rounds, BDE and method per campaign
$VENV/bin/python run_times.py       label=$H/runs/campaigns/run_<...>    # queue, planning, critic, compute; model work vs slot wait
```

"Correct" means the results critic passed the campaign and its parent and
fragments are the requested bond's homolysis (`EXPECTED` in `classify_board.py`;
add new molecules there).

## 5. Repeating the experiments

`experiments/launch/board_v2_chain.sh` runs the standard comparison set
(toluene benzylic C–H × 8, methanol O–H × 8, molecule set 1 twice, set 2)
with any model and server count, each run from its own copy of the results
database:

```bash
experiments/launch/board_v2_chain.sh nemotron-3-super:120b 4 <db-tag> service develop 7T 7M 8a 8b M
```

It expects `$COMPCHEM_HUB/db/pre_<N7_toluene|N7_methanol|N8a|N8b|NM>_<db-tag>.sqlite`;
see `experiments/README.md` for the other launchers and the expected folder layout.

## 6. Known problems

- **GPU selection with several servers.** Slurm sets `ROCR_VISIBLE_DEVICES=0-7`;
  the driver restricts each server by replacing that list. Adding
  `HIP_VISIBLE_DEVICES` on top leaves all but the first server on the CPU
  (`total_vram="0 B"` in `ollama_serve_<s>.log`), and their campaigns hang.
  Test a new setup with `experiments/launch/test_multi_ollama.sbatch`.
- **No timeout on LLM calls.** The planning budget is checked between calls, so
  a call to a stuck server holds its campaign until the job ends.
- **The staged driver waits for the slowest plan.** Compute starts only after
  every campaign has planned; one campaign that needs many board rounds sets
  the run time, and if the job ends first, finished plans are lost too.
- **Method drift with retrieval off.** Planners choose different basis sets for
  the same question (e.g. HF: 127.6 kcal/mol at 6-31G*, 132.7 at def2-SVP);
  turn retrieval on, or compare BDEs at the same method.
- **Scratch purge.** Runs and databases on Lustre scratch are deleted after
  ~90 days without access; copy anything you need to keep to `/ccs/proj` or
  project storage.
- **Driver log names.** The experiment launchers name `driver_<label>.log` by
  label, so rerunning a label overwrites the earlier driver log.
- **Planning dominates.** 85–99% of a run's wall time is planning and the plan
  critic; xtb and DFT take seconds for these small molecules.

## Older path: one `sbatch` job per campaign

`submit_campaigns.sh` (with `job_template.sh` and `ollama_server.sh`) submits
one single-node job per campaign plus one shared Ollama job, and polls them
with `sacct`. It predates the Parsl driver, uses a whole GPU node per campaign
and has no results database, staged mode or plan board; prefer
`submit_campaigns_parsl.py`. `./submit_campaigns.sh --help` lists its options.
