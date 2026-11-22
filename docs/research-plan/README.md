# Research plan: learning from the provenance of mistakes to improve the planner

**Question.** Can the provenance that Flowcept records for every planning
decision, including the wrong ones, be turned into a planner that makes fewer
mistakes, needs fewer plan-board rounds and fails less often, using only local
models on Frontier?

Pitch (two slides, for a non-specialist audience): [`pitch/bde-provenance-pitch.pdf`](pitch/bde-provenance-pitch.pdf),
editable with speaker notes: [`pitch/bde-provenance-pitch.pptx`](pitch/bde-provenance-pitch.pptx).

## What the logs already hold

Counted over all plan-board runs on Frontier (2026-10-20 to 2026-11-11):

| | |
|---|---|
| LLM calls recorded by Flowcept | 3,001: 1,673 planner, 1,296 plan critic, 32 supervisor (1,728 nemotron-3-nano, 706 Llama 3.1 405B, 567 nemotron-3-super) |
| Plan-board campaigns | 212 |
| Module attempts accepted / rejected | 887 / 625, of which 140 rejected as unusable replies (not a tool call) |
| Rejected → accepted pairs (same module and campaign) | 188 |
| Campaigns correct since plan board v2 | 189 of 192; the three misses ran out of planning time |

## The loop

1. **Run.** The local planner builds a plan on the board (code checks + LLM critic); xtb and DFT compute it.
2. **Record.** Flowcept logs every planner decision. Veritas (code.ornl.gov/ae4s/veritas) logs how the run was launched (command, git commit, environment, output hashes) and checks the result against a rule written in advance. One run ID joins the two.
3. **Blame and label.** Fixed rules settle clear cases (a server with no GPU is infrastructure; a parent that is not the named molecule is the planner). A person decides the rest. An LLM may draft a diagnosis, never the planner's own model.
4. **Learn** from verified planner mistakes only: a mistake bank in the prompt, new deterministic checks on the board, retraining the model.
5. **Gate.** A new planner version is used only if it passes at least as many verification items as the current one.

## Aims

- **A0 · Complete the provenance.** Record each code-check verdict and the driver's stage transitions, timeouts and lost tasks as Flowcept tasks, linked to the decision they judge, in one queryable store. Today the verdicts live only in `plan_<i>.log`, and run-level failures leave no trace in Flowcept.
- **A1 · Learn from mistakes without training.** (a) Mistake bank: retrieve past rejected attempts and their fixes for the same module and chemistry into the planner's prompt. (b) Learned checks: cluster failure signatures; keep a proposed check only if it flags every known-bad and passes every known-good plan in the history (e.g. the parent's formula must match the molecule named; a fragment's label must match its SMILES). (c) Root-cause reopen: when a module fails the same way against the same accepted upstream module, reopen the upstream one.
- **A2 · Retrain the local model.** LoRA fine-tuning on the logs: supervised examples from accepted modules, preference training (DPO) on the rejected → accepted pairs, distillation from nemotron-3-super into nemotron-3-nano.
- **A3 · Predict and escalate.** Predict budget timeouts from a campaign's first 3–5 rounds; hand predicted failures to a larger model or give them more time. Score the LLM critic against ground truth.

## Why Flowcept and Veritas are both needed

| | Flowcept (provenance) | Veritas (verification) |
|---|---|---|
| Watches | the planner agents inside a run | the operator agent outside it (Claude Code's commands) |
| Records | decisions and computations as they happen | command, commit, environment, output hashes |
| Judges | nothing: stores what the checks and critic said | each result, against a pass/fail rule written in advance, by rerun or a person |
| Answers | which decision led here? | is this result true, and whose fault is a failure? |

Uses for the planner: correct blame (only planner-caused failures are learned
from), verified training and retrieval data, Pact specs as the planner's
contract and as board checks, trust-weighted retrieval, human judgement on hard
cases, and a promotion gate. Needed first: a shared run ID, a
`bde_campaign.pact` spec, CompChem registered in Veritas, and a trust-status
column in the results database.

## Retraining

Export the provenance → split by molecule → LoRA fine-tune with PyTorch on
ROCm (one node, 8 × 64 GB) → merge → convert to GGUF → `ollama create` → rerun
`frontier-native/experiments/launch/board_v2_chain.sh`. nemotron-3-nano needs
its Hugging Face bf16 weights (≈ 60 GB) and its Mamba layers may lack fast ROCm
kernels; the fallback is Qwen2.5-7B-Instruct. The current data fixes reply
format and habits but is thin for chemistry judgement; nemotron-3-super can
generate more (8 campaigns in 11–15 min). Test molecules never appear in
training, retrieval or learned checks.

## Evaluation and timeline

- **Benchmark:** 30–50 molecules with reference BDEs, including multi-site molecules; held-out molecules and one held-out bond class.
- **Metrics:** correct rate, first-try acceptance per module, board rounds, budget timeouts, unusable replies, run time, error vs reference.
- **Baselines:** nano and super on board v2 (2026-11-11), retrieval off and on.
- **Weeks 1–2** A0 and critic scoring · **3–4** A1 · **5–7** A2 · **8** A3 and write-up.
- **Risks:** ROCm support for nemotron's Mamba layers; too little training data; gains hidden by an easy benchmark; learned checks over-fitting the history; Lustre scratch purge after 90 days.

## Worked examples (logs under `/lustre/orion/stf053/scratch/akhilesh/compchem/runs/campaigns/`)

| Case | Campaign | What happened | Provenance |
|---|---|---|---|
| A · calculations ran | B1 campaign 4, toluene benzylic C–H, Llama 3.1 405B | plan accepted in 6 rounds; 95.10 kcal/mol, QC PASS | `run_20261020_182457/task_4/plan/`, `species/species_2/`, `species/species_4/`, `task_4/` |
| B · wrong answer passed | B1 campaign 6 | fragment `Cc1[c]cccc1` labelled "benzyl" (a ring C–H); every check passed; 117.45 kcal/mol | `run_20261020_182457/task_6/` |
| C · no calculation (planner) | N8-b-8s campaign 4, hydrogen peroxide, nemotron-3-nano | parent written as `[OH]O[OH]` (H₂O₃), accepted after the critic errored; correct `[OH]` + `[OH]` rejected 12 times; budget ran out after 38 rounds | `run_20261111_123714/task_4/plan/flowcept_buffer.jsonl`, `plan_4.log` |
| D · no calculation (infrastructure) | first S8b attempt, campaign 1, ammonia, nemotron-3-super | plan accepted in 6 rounds; 3 of 4 Ollama servers had no GPU; compute never started | `run_20261102_175722/task_1/plan/`; failure only in `ollama_serve_1..3.log`, `parsl/000/htex_frontier/interchange.log`, Slurm job 5837388 |
