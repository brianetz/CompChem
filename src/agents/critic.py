"""
CriticAgent - applies whichever analysis skill's own review() to
ExecutorAgent's results, by looking the skill up in the skill registry
(src/skills/ - Phase 6, architecture plan Sec.16). Deliberately
skill-agnostic, same rationale as ExecutorAgent.

Phase 9 (architecture plan Sec.18) wraps review_analysis in a PROV-AGENT
'AgentTool' provenance record (storage/agent_provenance.py), attributing
the PASS/WARN/FAIL verdict to this agent's own identity in the provenance
graph.
"""
from __future__ import annotations

from academy.agent import Agent, action

from skills import get_skill
from storage.agent_provenance import agent_tool_task

_AGENT_ID = "CriticAgent"


class CriticAgent(Agent):
    """Applies physical-plausibility checks to ExecutorAgent's campaign
    results before they're allowed into a report."""

    @action
    async def review_analysis(self, results: dict) -> dict:
        skill = get_skill(results["analysis_type"])
        with agent_tool_task(
            _AGENT_ID, "review_analysis", used={"analysis_type": results.get("analysis_type")},
        ) as t:
            verdict = skill.review(results)
            t.end(generated=verdict)
        return verdict
