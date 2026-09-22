---
name: compchem-workflow
description: Use when working in the CompChem repo - running or demoing the agentic plan/approve/execute/critique/approve/report campaign pipeline, explaining how it works, adding a new analysis skill, or troubleshooting a run.
---

# CompChem agentic workflow

This repo is an agentic computational chemistry tool: a plan → approve → execute → critique → approve → report loop that turns a molecule (or a plain-language question) into a computed result, with a human approving at two points and every step captured as provenance. Full narrative detail lives in `README.md` at the repo root - read it once if you need more than this skill gives you, but this file should be enough to run and explain the workflow without re-deriving anything.

## Two generations of code - know which one you're in

- **`src/agents/`, `src/skills/`, `src/engines/`, `src/llm/`, `src/storage/`** - the active, current workflow. Everything below is about this.
- **`src/conformer.py`, `fragmenter.py`, `nwchem.py`, `bde_workflow.py`, `run_all.py`, `parser.py`** - the original hand-written BDE pipeline this project started from. Still functional, still Parsl-orchestrated, but not where new work lands. Don't confuse a question about "the workflow" with this older code unless the user specifically asks about the original BDE scripts.

## How the pipeline works

Four real Academy actor agents, orchestrated by `agents/run_agentic_campaign.py`:

1. **PlannerAgent** - a bare SMILES gets a rule-based spec; a natural-language `--question` gets a real Claude API call that reasons about what's needed (e.g. recognizing that a bond-dissociation-energy question needs a parent *and* fragments, with correct charge/multiplicity - not just one molecule).
2. **Human approval gate 1** - pauses before anything computes. Approve or abort.
3. **ExecutorAgent** - dispatches to a pluggable skill under `src/skills/` (currently `single_point`, `bde`): GFN2-xTB (xtb) prescreen optimization, then a DFT refinement (PySCF by default; NWChem is wired in as an ASE-calculator alternate but not yet planner-selectable or runtime-validated - don't imply it works end-to-end yet).
4. **CriticAgent** - PASS / WARN / FAIL on convergence, missing energies, spin-state sanity.
5. **Human approval gate 2** - approve (finalize) / reprioritize (feedback goes back through the planner for a revised spec, a real Claude call, bounded by `--max-replan-rounds`) / new_experiment (fresh SMILES or question) / abort.
6. **ReporterAgent** - writes `agentic_campaign_report.md` with every species' energy and a rendered structure image.

Every agent decision and every LLM call, including each human gate choice, is captured as its own Flowcept-backed provenance record (`agent_tool` / `ai_model_invocation`), not just the underlying calculation tasks.

## Running it

Always `cd src` first - the modules import as `agents.*`/`skills.*` relative to that directory, so `python -m agents.run_agentic_campaign` from the repo root fails with `ModuleNotFoundError: No module named 'agents'`. This is the single most common way a first run goes wrong; check the working directory before troubleshooting anything else.

```bash
conda activate compchem-agent   # env from environment.yml at repo root
cd src
```

Quick rule-based smoke test, gates skipped:
```bash
python -m agents.run_agentic_campaign --auto-approve
```

Natural-language question, LLM-planned, gates live (needs `ANTHROPIC_API_KEY` set in the shell):
```bash
python -m agents.run_agentic_campaign --question "What's the O-H bond dissociation energy for water in kcal/mol, using a cheap xtb prescreen and a DFT refinement?" --model claude-haiku-4-5-20251001
```

At the gates: plan gate is `y`/`n`; post-analysis gate is `a` (approve) / `r` (reprioritize, then feedback text) / `n` (new_experiment, then a new question or SMILES) / `x` (abort).

Flags: `--question "..."`, `--model <name>` (default `claude-sonnet-5`, or `$ANTHROPIC_MODEL`), `--auto-approve` (skip both gates), `--max-replan-rounds <N>` (default 2).

Output: terminal shows the spec/energies/verdict live; `agentic_campaign_report.md` + `report_structures/*.png` (both gitignored, regenerated per run) are the human-facing artifacts; `flowcept_buffer.jsonl` (gitignored) holds full provenance, inspectable with `python verify_agent_provenance.py`.

## Adding a new analysis type

New module under `src/skills/` exposing `NAME`, `DESCRIPTION`, `run(spec)`, `review(results)`, `describe(spec, results, verdict)`, registered in `skills/__init__.py`. `PlannerAgent`'s LLM tool schema and `ExecutorAgent`/`CriticAgent`/`ReporterAgent` all pick it up automatically from the registry - don't hand-edit those three classes for a new analysis type, that defeats the point of the registry.

## Current status / known gaps (say this plainly if asked "does X work")

- Local path (RDKit → xtb → PySCF, both analysis types, both gates) is complete and working.
- NWChem is wired as an ASE-calculator refinement engine but is **not yet selectable by either planner** and **not yet runtime-validated** against a real binary through the agent pipeline - it compiles/imports cleanly, that's all that's confirmed so far.
- HPC deployment (OLCF Defiant) is in progress outside this repo's own code: NWChem 7.2.3 built successfully via Spack; a Slurm functional test of that binary is still being worked through. A Globus Compute endpoint exists and is confirmed running, paused pending that validation. None of this is reachable or runnable from a normal dev machine - don't suggest running against Defiant as part of a quick demo.

## Troubleshooting

- `ModuleNotFoundError: No module named 'agents'` → not `cd`'d into `src/` (see above).
- A run seems to hang with no output → an approval gate is waiting on `input()` at the terminal (gates are on by default); either respond to the prompt or re-run with `--auto-approve`.
- `--model <alias>` not found → try the fully-dated model name (e.g. `claude-haiku-4-5-20251001` instead of `claude-haiku-4-5`).
- Missing structure images in the report → `pymol` isn't on `PATH` in this environment (the renderer shells out to it); check `which pymol`, image rendering fails gracefully (report still writes, just without that image) with a one-line diagnostic printed to the terminal.
