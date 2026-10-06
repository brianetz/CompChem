"""
PlannerAgent - Phase 4 gave this a simple rule-based campaign-spec builder
(`plan_campaign`, explicit molecule/method parameters in, structured spec
out). Phase 5 adds a second, LLM-backed action (`plan_campaign_from_prompt`)
that turns a natural-language research question into the same spec shape
via the Anthropic Claude API (see `llm.planner_llm`) - Executor/Critic/
Reporter don't need to know or care which action produced the spec they're
handed. Phase 9 (architecture plan Sec.18) wraps both actions in a
PROV-AGENT 'AgentTool' provenance record (storage/agent_provenance.py),
attributing each spec to this agent's own identity in the provenance graph.
Phase 10 adds a third action, `replan_campaign`, invoked from the
post-analysis human-approval gate's 'reprioritize' choice
(agents/human_gate.py) - it revises a previous spec in light of its own
results, the critic's verdict, and free-text human feedback.
"""
from __future__ import annotations

from academy.agent import Agent, action

from llm.planner_llm import (
    generate_campaign_spec_llm,
    generate_revised_campaign_spec_llm,
)
from storage.agent_provenance import agent_tool_task

_AGENT_ID = "PlannerAgent"


class PlannerAgent(Agent):
    """Turns a molecule specification - explicit or natural-language - into
    a structured campaign spec."""

    @action
    async def plan_campaign(
        self,
        smiles: str,
        charge: int = 0,
        multiplicity: int = 1,
        basis: str = "6-31g*",
        xc: str = "b3lyp",
    ) -> dict:
        """Rule-based planning (Phase 4): build a spec directly from
        explicit parameters. No LLM call, no network dependency. Always the
        'single_point' analysis skill (src/skills/single_point.py, Phase 6)
        - this action has no way to express a multi-species analysis."""
        with agent_tool_task(
            _AGENT_ID, "plan_campaign",
            used={"smiles": smiles, "charge": charge, "multiplicity": multiplicity, "basis": basis, "xc": xc},
        ) as t:
            spec = {
                "analysis_type": "single_point",
                "smiles": smiles,
                "charge": charge,
                "multiplicity": multiplicity,
                "method_tiers": ["xtb-GFN2", f"pyscf-{xc}/{basis}"],
                "basis": basis,
                "xc": xc,
            }
            t.end(generated=spec)
        return spec

    @action
    async def plan_campaign_from_prompt(
        self, question: str, model: str | None = None, backend: str = "anthropic"
    ) -> dict:
        """LLM-backed planning (Phase 5): turns a natural-language research
        question into the same campaign spec shape `plan_campaign` returns.
        `backend` selects which LLM answers the tool call - 'anthropic'
        (default, the hosted Claude API), 'ollama' (a local model via a
        running `ollama serve`, CLI: --ollama), or 'claude_cli' (the
        logged-in Claude Code CLI via `claude -p`, CLI: --claude-cli, for
        no-API-key use). Either way this is a blocking network/subprocess
        call, so it's dispatched through
        agent_run_sync like Executor's engine calls. The actual
        prompt/response cycle(s) - including any Phase-6-added bounded
        self-correction retry - are captured separately, one level down,
        inside generate_campaign_spec_llm itself (llm/planner_llm.py) as
        PROV-AGENT 'AIModelInvocation' records; this outer task instead
        records the PlannerAgent's own tool invocation as a whole."""
        with agent_tool_task(
            _AGENT_ID, "plan_campaign_from_prompt",
            used={"question": question, "model": model, "backend": backend},
        ) as t:
            spec = await self.agent_run_sync(
                generate_campaign_spec_llm, question, model=model, backend=backend,
            )
            t.end(generated=spec)
        return spec

    @action
    async def replan_campaign(
        self,
        previous_spec: dict,
        previous_results: dict,
        verdict: dict,
        feedback: str,
        model: str | None = None,
        backend: str = "anthropic",
    ) -> dict:
        """Phase 10: revise a previous campaign spec in light of its own
        results, the critic's verdict, and free-text human feedback given at
        the post-analysis approval gate (agents/human_gate.py,
        'reprioritize' choice). Always goes through the LLM path
        (generate_revised_campaign_spec_llm) regardless of whether the
        original spec came from plan_campaign or plan_campaign_from_prompt -
        interpreting free-text feedback needs language understanding either
        way. `backend` selects 'anthropic' (default), 'ollama', or
        'claude_cli', same as plan_campaign_from_prompt above. Same
        agent_run_sync dispatch and
        PROV-AGENT 'AgentTool' wrapping as the other two actions above."""
        with agent_tool_task(
            _AGENT_ID, "replan_campaign",
            used={
                "previous_analysis_type": previous_spec.get("analysis_type"),
                "verdict": verdict,
                "feedback": feedback,
                "model": model,
                "backend": backend,
            },
        ) as t:
            spec = await self.agent_run_sync(
                generate_revised_campaign_spec_llm,
                previous_spec, previous_results, verdict, feedback, model=model, backend=backend,
            )
            t.end(generated=spec)
        return spec
