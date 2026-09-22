"""
ExecutorAgent - runs whichever analysis skill a campaign spec names
(spec['analysis_type']), by looking it up in the skill registry
(src/skills/ - Phase 6, architecture plan Sec.16) and running it. This agent
is deliberately skill-agnostic: adding a new analysis type (IR frequencies,
NMR shifts, UV-Vis/TDDFT, HOMO-LUMO gap, redox potentials, Marcus theory,
...) means adding a new skill module under src/skills/, not touching this
class.

Phase 9 (architecture plan Sec.18) wraps run_analysis in a PROV-AGENT
'AgentTool' provenance record (storage/agent_provenance.py), attributing
the skill dispatch to this agent's own identity in the provenance graph -
skills themselves stay agent-agnostic, same rationale as everywhere else in
this file.
"""
from __future__ import annotations

from academy.agent import Agent, action

from skills import get_skill
from storage.agent_provenance import agent_tool_task

_AGENT_ID = "ExecutorAgent"


class ExecutorAgent(Agent):
    """Executes a campaign spec end-to-end and returns JSON-safe results."""

    @action
    async def run_analysis(self, spec: dict) -> dict:
        skill = get_skill(spec["analysis_type"])
        with agent_tool_task(
            _AGENT_ID, "run_analysis", used={"analysis_type": spec.get("analysis_type")},
        ) as t:
            results = await self.agent_run_sync(skill.run, spec)
            t.end(generated={"analysis_type": results.get("analysis_type"), "status": results.get("status")})
        return results
