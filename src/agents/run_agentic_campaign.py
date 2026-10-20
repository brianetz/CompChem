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
        [--ollama]                                             # use a local Ollama model instead of the Anthropic API
        [--claude-cli]                                         # use the logged-in Claude Code CLI instead of the Anthropic API
        [--auto-approve]                                       # skip both human-approval gates (Phase 10) - scripted/test runs
        [--max-replan-rounds 2]                                # cap on 'reprioritize' round-trips through the planner (default: 2)
        [--plan-only spec.json]                                # stop after the plan gate; write the spec (staged driver, stage 1)
        [--from-spec spec.json]                                # skip planning; execute this spec (staged driver, stage 3)
        [--no-rag]                                             # don't add retrieved prior results to the planner prompt

With $COMPCHEM_DB set, every species goes through the shared results
database (storage/results_db.py): already-computed species are reused, and
concurrent campaigns needing the same species compute it once.

Without --auto-approve, this campaign now PAUSES TWICE for a terminal
prompt: once after planning, once after the critic's review. Existing
scripts/tests that ran this non-interactively need --auto-approve added to
keep running unattended, exactly as before this feature existed.

The LLM-backed path (--question, or 'reprioritize' at the post-analysis
gate) needs one of three backends (--ollama and --claude-cli are mutually
exclusive; neither given means the Anthropic default):
  - Anthropic (default): `pip install anthropic` and an ANTHROPIC_API_KEY
    in the environment. Makes a real (billed) network call to the Claude API.
  - Ollama (--ollama flag): `pip install ollama` and a local `ollama serve`
    running, with the chosen model already pulled (`ollama pull llama3.1`
    or similar). No API key, no network call, no billing - see
    llm/planner_llm.py's module docstring for the tool-calling caveats of
    this path (it's less reliable than Claude's forced tool_choice).
  - Claude Code CLI (--claude-cli flag): Claude Code itself installed and
    logged in (`claude` on PATH). Shells out to `claude -p`, authenticating
    with the CLI's own login instead of ANTHROPIC_API_KEY - for a user with
    a Claude Code subscription but no separate API key. Still counts
    against that login's usage/plan, and is slower than a direct API call
    since each attempt spawns a real (tool-less) Claude Code session - see
    llm/planner_llm.py's module docstring for how the tool schema is
    conveyed (plain-text instructions, not native tool-calling) and its
    own failure modes.
"""
from __future__ import annotations

import argparse
import asyncio
import json
import os
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
        help="Override the default model for --question runs: the ANTHROPIC_MODEL "
        "default (used for both the Anthropic and --claude-cli backends), or "
        "with --ollama, the OLLAMA_MODEL default.",
    )
    backend_group = parser.add_mutually_exclusive_group()
    backend_group.add_argument(
        "--ollama",
        action="store_true",
        help="Use a local Ollama model instead of the Anthropic Claude API for "
        "LLM-backed planning (--question, and any 'reprioritize' replanning). "
        "No ANTHROPIC_API_KEY needed. Requires `pip install ollama` and a "
        "running `ollama serve` with the model already pulled (default model: "
        "llama3.1, override with --model or $OLLAMA_MODEL).",
    )
    backend_group.add_argument(
        "--claude-cli",
        action="store_true",
        help="Use the logged-in Claude Code CLI (`claude -p`) instead of the "
        "Anthropic Claude API for LLM-backed planning. No ANTHROPIC_API_KEY "
        "needed - authenticates via Claude Code's own login instead. Requires "
        "Claude Code installed and logged in (`claude` on PATH). Still counts "
        "against that login's usage/plan, and is slower than a direct API "
        "call (see module docstring).",
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
    parser.add_argument(
        "--plan-only",
        metavar="SPEC_JSON",
        default=None,
        help="Plan (and pass the plan gate), write the approved spec to this "
        "JSON file, and stop - no compute. Stage 1 of the staged driver.",
    )
    parser.add_argument(
        "--from-spec",
        metavar="SPEC_JSON",
        default=None,
        help="Skip planning and execute the spec in this JSON file (written by "
        "--plan-only). Stage 3 of the staged driver.",
    )
    parser.add_argument(
        "--no-rag",
        action="store_true",
        help="Don't add prior results retrieved from $COMPCHEM_DB to the "
        "planner's prompt (storage/rag.py).",
    )
    parser.add_argument(
        "--plan-chat",
        action="store_true",
        help="Plan --question as a group chat of Academy agents (Planner, "
        "Validator, and the plan critic if $COMPCHEM_PLAN_CRITIC_SKILL is set) "
        "under a round-robin moderator that supervises the discussion "
        "(agents/plan_discussion.py), instead of the planner's single retry.",
    )
    parser.add_argument(
        "--plan-chat-rounds",
        type=int,
        default=3,
        help="Maximum discussion rounds with --plan-chat (default: 3).",
    )
    parser.add_argument(
        "--plan-board",
        action="store_true",
        help="Build the --question plan module by module on a shared board: the "
        "Planner publishes each module, the Validator and the plan critic "
        "($COMPCHEM_PLAN_CRITIC_SKILL, default src/bde-plan-critic.md) review "
        "it, and a supervising moderator accepts modules and the final plan "
        "(agents/plan_board.py).",
    )
    parser.add_argument(
        "--plan-board-attempts",
        type=int,
        default=3,
        help="Attempts per module with --plan-board (default: 3).",
    )
    return parser.parse_args()


async def _plan_on_board(manager, args: argparse.Namespace, backend: str) -> dict:
    """Launch the plan-board agents for one question and return the plan
    they agreed on (raises ValueError if none passed)."""
    from agents.plan_board import (
        BoardCriticAgent, BoardPlannerAgent, BoardValidatorAgent, PlanBoard, PlanBoardModerator,
    )
    from storage.rag import build_retrieval_context

    board = await manager.launch(PlanBoard())
    planner = await manager.launch(BoardPlannerAgent, kwargs={"board": board, "backend": backend, "model": args.model})
    validator = await manager.launch(BoardValidatorAgent, kwargs={"board": board})
    critic = await manager.launch(BoardCriticAgent, kwargs={"board": board, "backend": backend, "model": args.model})
    moderator = await manager.launch(PlanBoardModerator, kwargs={
        "board": board, "planner": planner, "validator": validator, "critic": critic,
        "backend": backend, "model": args.model, "max_attempts": args.plan_board_attempts})
    context, retrieval = build_retrieval_context(args.question)
    result = await moderator.plan(args.question, context)
    print("Plan board:", result["outcome"], "after", len(result["rounds"]), "reviewed round(s);",
          "supervisor:", result["supervisor"] or "no action")
    if result["spec"] is None:
        raise ValueError(f"plan board ended without a valid plan ({result['outcome']}): {result['rounds']}")
    spec = result["spec"]
    spec["retrieval"] = retrieval
    spec["plan_board"] = {k: result[k] for k in ("outcome", "rounds", "supervisor", "final_reviews", "board")}
    return spec


async def _plan_by_discussion(manager, args: argparse.Namespace, backend: str) -> dict:
    """Launch the plan-discussion agents for one question and return the
    spec the discussion agreed on (raises ValueError if none passed)."""
    from agents.plan_discussion import (
        CriticChatAgent, PlanDiscussionManager, PlannerChatAgent, ValidatorChatAgent,
    )
    from llm.planner_llm import _resolve_model
    from llm.plan_critic_llm import skill_path
    from storage.rag import build_retrieval_context

    model = _resolve_model(backend, args.model)
    participants = [
        await manager.launch(PlannerChatAgent(backend, args.model)),
        await manager.launch(ValidatorChatAgent(model)),
    ]
    if skill_path():
        participants.append(await manager.launch(CriticChatAgent(backend, args.model)))
    moderator = await manager.launch(
        PlanDiscussionManager,
        kwargs={"participants": participants, "backend": backend, "model": args.model,
                "max_rounds": args.plan_chat_rounds},
    )
    context, retrieval = build_retrieval_context(args.question)
    result = await moderator.plan(args.question, context)
    print("Plan discussion:", result["outcome"], "after", len(result["rounds"]), "round(s);",
          "supervisor:", result["supervisor"] or "no action")
    if result["spec"] is None:
        raise ValueError(f"plan discussion ended without a valid plan ({result['outcome']}): "
                         f"{result['rounds']}")
    spec = result["spec"]
    spec["retrieval"] = retrieval
    spec["plan_discussion"] = {k: result[k] for k in ("outcome", "rounds", "supervisor", "history")}
    return spec


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


async def run_campaign(args: argparse.Namespace) -> str | None:
    if args.ollama:
        backend = "ollama"
    elif args.claude_cli:
        backend = "claude_cli"
    else:
        backend = "anthropic"

    async with await Manager.from_exchange_factory(
        factory=LocalExchangeFactory(),
        executors=ThreadPoolExecutor(),
    ) as manager:
        planner = await manager.launch(PlannerAgent())
        executor = await manager.launch(ExecutorAgent())
        critic = await manager.launch(CriticAgent())
        reporter = await manager.launch(ReporterAgent())

        if args.from_spec:
            with open(args.from_spec) as f:
                spec = json.load(f)
        elif args.question and args.plan_board:
            spec = await _plan_on_board(manager, args, backend)
        elif args.question and args.plan_chat:
            spec = await _plan_by_discussion(manager, args, backend)
        elif args.question:
            spec = await planner.plan_campaign_from_prompt(
                question=args.question, model=args.model, backend=backend
            )
        else:
            spec = await planner.plan_campaign(smiles=args.smiles or "O")
        print("Planner spec:", spec)

        if args.plan_only:
            plan_decision = await request_plan_approval(spec, auto_approve=args.auto_approve)
            if plan_decision.action == "abort":
                return _abort_report(spec, plan_decision, stage="plan")
            with open(args.plan_only, "w") as f:
                json.dump(spec, f, indent=2, default=str)
            print(f"Spec written to {args.plan_only}")
            return None

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
                        question=analysis_decision.new_question, model=args.model, backend=backend
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
                    feedback=analysis_decision.feedback or "", model=args.model, backend=backend,
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
    if args.no_rag:
        os.environ["COMPCHEM_RAG"] = "0"
    with provenance_session():
        report = asyncio.run(run_campaign(args))
    if report is None:  # --plan-only: the spec file is the output
        return

    out_path = "agentic_campaign_report.md"
    with open(out_path, "w") as f:
        f.write(report)

    print("\n" + report)
    print(f"\nReport written to {out_path}")


if __name__ == "__main__":
    main()
