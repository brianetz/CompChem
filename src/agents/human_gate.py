"""
Human-in-the-loop approval gates for run_agentic_campaign.py (Phase 10).

Two gates sit in the plan -> execute -> critique -> report pipeline:

1. request_plan_approval - right after the planner produces a spec, before
   ExecutorAgent.run_analysis dispatches ANY compute. This deliberately
   combines what could be two separate checkpoints ("is this plan right"
   and "ok to actually run it") into one, since today those are the exact
   same moment in the code - ExecutorAgent.run_analysis (agents/executor.py)
   is a single synchronous call with no distinct "dispatch" sub-step to
   gate separately. Revisit this once Phase 8's Globus Compute swap makes
   "submission" a materially different, higher-stakes action (leaving the
   laptop and landing on Defiant) - at that point a second, later gate
   inside ExecutorAgent, right before the `agent_run_sync(skill.run, spec)`
   call the architecture doc already earmarks for that swap, would be the
   natural place to split them.

2. request_analysis_approval - after the critic's review, before the
   reporter finalizes anything. Offers four outcomes: approve (finalize as
   normal), reprioritize (PlannerAgent.replan_campaign revises the spec
   using the previous spec/results/verdict plus free-text human feedback,
   then the loop re-runs execute -> critique - bounded by
   --max-replan-rounds so this can never loop forever, mirroring the
   existing "cap recovery attempts" pattern already used for the LLM
   self-correction retry in llm/planner_llm.py), new_experiment (start over
   from a fresh question/SMILES), or abort (write a short abort report and
   stop).

Both gates are plain blocking `input()` prompts, not some new async
signaling mechanism - there is nothing else to pause/resume on anywhere in
this codebase (no asyncio.Event, no callback hook), and
run_agentic_campaign.py is a synchronous CLI tool one person runs
interactively, so a blocking prompt is the right level of engineering here.
`--auto-approve` bypasses both gates entirely (always returns "approve",
instantly, no prompt) so scripted/unattended/test runs keep working exactly
as they did before this feature existed - this is NOT the default; gates
are on unless you pass that flag.

Every decision - human-entered or auto-approved - is logged through the
same PROV-AGENT "AgentTool" mechanism (storage/agent_provenance.py) every
other agent decision in this pipeline already uses, under agent_id
"HumanApprovalGate", so an approval/rejection/revision sits in the same
provenance graph as everything else, not off to the side.
"""
from __future__ import annotations

from dataclasses import dataclass

from storage.agent_provenance import agent_tool_task

_AGENT_ID = "HumanApprovalGate"


@dataclass
class ApprovalDecision:
    action: str  # "approve" | "abort" | "reprioritize" | "new_experiment"
    feedback: str | None = None
    new_question: str | None = None
    new_smiles: str | None = None


def _prompt_choice(prompt: str, choices: dict[str, str]) -> str:
    """Blocking prompt, re-asks on an unrecognized answer. `choices` maps
    the single-letter answer to its full action name."""
    options = "/".join(f"{letter}={name}" for letter, name in choices.items())
    while True:
        answer = input(f"{prompt} [{options}]: ").strip().lower()
        if answer in choices:
            return choices[answer]
        print(f"Unrecognized answer {answer!r} - please enter one of: {', '.join(choices)}")


async def request_plan_approval(spec: dict, *, auto_approve: bool = False) -> ApprovalDecision:
    """Gate 1+2 combined (see module docstring): approve the plan / abort."""
    with agent_tool_task(_AGENT_ID, "approve_plan", used={"spec": spec}) as t:
        if auto_approve:
            decision = ApprovalDecision(action="approve")
            t.end(generated={"decision": decision.action, "auto_approved": True})
            return decision

        print("\n=== Human approval requested: campaign plan ===")
        print(spec)
        action = _prompt_choice(
            "Approve this plan and run it?", {"y": "approve", "n": "abort"},
        )
        reason = None
        if action == "abort":
            reason = input("Reason (optional, for the record): ").strip() or None
        decision = ApprovalDecision(action=action, feedback=reason)
        t.end(generated={"decision": decision.action, "reason": reason, "auto_approved": False})
        return decision


async def request_analysis_approval(
    spec: dict, results: dict, verdict: dict, *, auto_approve: bool = False
) -> ApprovalDecision:
    """Gate 3: approve the finished analysis, send it back for
    reprioritization, start a brand-new experiment, or abort."""
    with agent_tool_task(
        _AGENT_ID, "approve_analysis",
        used={"analysis_type": spec.get("analysis_type"), "qc_flag": verdict.get("qc_flag")},
    ) as t:
        if auto_approve:
            decision = ApprovalDecision(action="approve")
            t.end(generated={"decision": decision.action, "auto_approved": True})
            return decision

        print("\n=== Human approval requested: analysis complete ===")
        print("Results:", results)
        print("Verdict:", verdict)
        action = _prompt_choice(
            "What next?",
            {"a": "approve", "r": "reprioritize", "n": "new_experiment", "x": "abort"},
        )
        decision = ApprovalDecision(action=action)
        if action == "reprioritize":
            decision.feedback = input("What should change for the revised plan? ").strip()
        elif action == "new_experiment":
            question = input(
                "New research question (leave blank to give a SMILES instead): "
            ).strip()
            if question:
                decision.new_question = question
            else:
                decision.new_smiles = input("New SMILES: ").strip() or "O"
        elif action == "abort":
            decision.feedback = input("Reason (optional, for the record): ").strip() or None

        t.end(generated={
            "decision": decision.action,
            "feedback": decision.feedback,
            "new_question": decision.new_question,
            "new_smiles": decision.new_smiles,
            "auto_approved": False,
        })
        return decision
