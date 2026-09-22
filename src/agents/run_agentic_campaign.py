"""
Phase 4/5/6/10 driver + smoke test: wires PlannerAgent -> ExecutorAgent ->
CriticAgent -> ReporterAgent together through Academy's Manager, running
entirely locally (LocalExchangeFactory + ThreadPoolExecutor - no Redis, no
HPC). Phase 4 gave every stage a real Academy `Agent`; Phase 5 added an
LLM-backed planning path (PlannerAgent.plan_campaign_from_prompt, via the
Anthropic Claude API) alongside the original rule-based one, selected here
by which CLI argument you pass. Phase 6 made Executor/Critic/Reporter
skill-agnostic (src/skills/ - architecture plan Sec.16): this driver no
longer branches on campaign/analysis type at all - whichever analysis_type
the planner picked (single_point, bde, or any future registered skill) is
dispatched generically via run_analysis/review_analysis/write_analysis_report.

Phase 10 adds two human-approval gates (agents/human_gate.py) around that
same pipeline: one right after planning (before ANY compute runs), one
after the critic's review (before the report is finalized). The
post-analysis gate can send the spec back through the planner for revision
("reprioritize"), start over with a fresh question/SMILES
("new_experiment"), finalize as-is ("approve"), or stop ("abort") - this is
the first real branching this driver has ever had; previously the critic's
verdict was computed but never acted on.

Usage:
    python -m agents.run_agentic_campaign                     # default: water ("O"), rule-based planner
    python -m agents.run_agentic_campaign CCO                  # explicit SMILES, rule-based planner
    python -m agents.run_agentic_campaign --question "..."     # natural-language question, LLM-backed planner (Phase 5)
        [--model claude-opus-5]                                # optional: override ANTHROPIC_MODEL for this run
        [--auto-approve]                                       # skip both human-approval gates (Phase 10) - scripted/test runs
        [--max-replan-rounds 2]                                # cap on 'reprioritize' round-trips through the planner (default: 2)

Without --auto-approve, this campaign now PAUSES TWICE for a terminal
prompt: once after planning, once after the critic's review. Existing
scripts/tests that ran this non-interactively need --auto-approve added to
keep running unattended, exactly as before this feature existed.

The LLM-backed path requires `pip install anthropic` and an ANTHROPIC_API_KEY
in the environment, and makes a real (billed) network call to the Claude API.
"""
from __future__ import annotations

import argparse
import asyncio
from concurrent.futures import ThreadPoolExecutor

from academy.exchange import LocalExchangeFactory
from academy.manager import Manager

from agents.planner import PlannerAgent
from agents.executor import ExecutorAgent
from agents.critic import CriticAgent
from agents.reporter import ReporterAgent
from agents.human_gate import request_plan_approval, request_analysis_approval
from storage.flowcept_store import provenance_session


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "smiles",
        nargs="?",
        default=None,
        help='SMILES string for the rule-based planner (default: "O", water). '
        "Ignored if --question is given.",
    )
    parser.add_argument(
        "--question",
        default=None,
        help="Natural-language research question for the LLM-backed planner (Phase 5).",
    )
    parser.add_argument(
        "--model",
        default=None,
        help="Override the ANTHROPIC_MODEL default for --question runs.",
    )
    parser.add_argument(
        "--auto-approve",
        action="store_true",
        help="Skip both human-approval gates (Phase 10, agents/human_gate.py) "
        "- for scripted/unattended/test runs. Without this flag, the "
        "campaign pauses twice and prompts on the terminal.",
    )
    parser.add_argument(
        "--max-replan-rounds",
        type=int,
        default=2,
        help="Cap on how many times the post-analysis 'reprioritize' choice "
        "can send the spec back through the planner before the campaign "
        "just finalizes with whatever it currently has (default: 2).",
    )
    return parser.parse_args()


def _abort_report(
    spec: dict, decision, *, stage: str, results: dict | None = None, verdict: dict | None = None
) -> str:
    """Phase 10: a short markdown report for a campaign stopped at a human-
    approval gate, so main()'s normal report-writing path doesn't need any
    special-casing for the abort outcome."""
    lines = ["# Campaign Aborted", "", f"Aborted by human decision at the **{stage}** approval gate."]
    if decision.feedback:
        lines += ["", f"Reason given: {decision.feedback}"]
    lines += ["", "## Spec at time of abort", "", "```", str(spec), "```"]
    if results is not None:
        lines += ["", "## Results at time of abort", "", "```", str(results), "```"]
    if verdict is not None:
        lines += ["", "## Verdict at time of abort", "", "```", str(verdict), "```"]
    return "\n".join(lines)


async def run_campaign(args: argparse.Namespace) -> str:
    async with await Manager.from_exchange_factory(
        factory=LocalExchangeFactory(),
        executors=ThreadPoolExecutor(),
    ) as manager:
        planner = await manager.launch(PlannerAgent())
        executor = await manager.launch(ExecutorAgent())
        critic = await manager.launch(CriticAgent())
        reporter = await manager.launch(ReporterAgent())

        if args.question:
            spec = await planner.plan_campaign_from_prompt(
                question=args.question, model=args.model
            )
        else:
            spec = await planner.plan_campaign(smiles=args.smiles or "O")
        print("Planner spec:", spec)

        results: dict | None = None
        verdict: dict | None = None

        # Phase 10: this loop is the first real branching this driver has
        # ever had - the critic's verdict used to be computed and printed
        # but never acted on. Bounded by --max-replan-rounds so a
        # 'reprioritize' choice can never loop forever (same "cap recovery
        # attempts" philosophy as the LLM self-correction retry in
        # llm/planner_llm.py).
        for round_num in range(args.max_replan_rounds + 1):
            plan_decision = await request_plan_approval(spec, auto_approve=args.auto_approve)
            if plan_decision.action == "abort":
                return _abort_report(spec, plan_decision, stage="plan")

            # Phase 6: skill-agnostic dispatch - whichever analysis_type the
            # planner chose (single_point, bde, or any future registered skill
            # in src/skills/) is looked up and run generically. No per-type
            # branching lives here anymore; see src/skills/__init__.py.
            results = await executor.run_analysis(spec)
            print("Executor results:", results)

            verdict = await critic.review_analysis(results)
            print("Critic verdict:", verdict)

            analysis_decision = await request_analysis_approval(
                spec, results, verdict, auto_approve=args.auto_approve
            )

            if analysis_decision.action == "approve":
                break
            if analysis_decision.action == "abort":
                return _abort_report(
                    spec, analysis_decision, stage="analysis", results=results, verdict=verdict
                )
            if analysis_decision.action == "new_experiment":
                if analysis_decision.new_question:
                    spec = await planner.plan_campaign_from_prompt(
                        question=analysis_decision.new_question, model=args.model
                    )
                else:
                    spec = await planner.plan_campaign(smiles=analysis_decision.new_smiles or "O")
                print("Planner spec (new experiment):", spec)
                continue
            if analysis_decision.action == "reprioritize":
                if round_num >= args.max_replan_rounds:
                    print(
                        f"Reached --max-replan-rounds={args.max_replan_rounds}; "
                        "finalizing report with current results instead of replanning again."
                    )
                    break
                spec = await planner.replan_campaign(
                    previous_spec=spec, previous_results=results, verdict=verdict,
                    feedback=analysis_decision.feedback or "", model=args.model,
                )
                print("Planner spec (revised):", spec)
                continue

        report = await reporter.write_analysis_report(spec, results, verdict)

    # No explicit per-agent .shutdown() calls here: the Manager's own
    # __aexit__ (triggered by the `async with` block above ending) already
    # tears down every launched agent. Calling .shutdown() on each handle
    # too, before that, raced with the manager's own teardown and produced
    # harmless-but-noisy "Response received without corresponding request"
    # warnings (one per agent) as each shutdown ack arrived after its
    # mailbox was already gone.
    return report


def main() -> None:
    args = parse_args()
    with provenance_session():
        report = asyncio.run(run_campaign(args))

    out_path = "agentic_campaign_report.md"
    with open(out_path, "w") as f:
        f.write(report)

    print("\n" + report)
    print(f"\nReport written to {out_path}")


if __name__ == "__main__":
    main()
