"""
ReporterAgent - turns a campaign spec + Executor results + Critic verdict
into a short, human-readable markdown report. The header (title, timestamp,
verdict, research question, planner rationale) is generic; the body is
each analysis skill's own describe() (src/skills/ - Phase 6, architecture
plan Sec.16). Deliberately skill-agnostic, same rationale as ExecutorAgent/
CriticAgent.

Phase 7 (architecture plan Sec.17) added an image_dir parameter, passed
straight through to the skill's describe(): skills may render a PyMOL PNG
of each species' structure under that directory and embed it in the report
body. ReporterAgent itself stays agnostic to rendering, same as it's
agnostic to every other analysis-specific detail.

Phase 9 (architecture plan Sec.18) wraps write_analysis_report in a
PROV-AGENT 'AgentTool' provenance record (storage/agent_provenance.py),
attributing the finished report text itself to this agent's own identity in
the provenance graph - directly serving this project's original vision
(architecture plan Sec.1) of a report "traceable back to the exact code,
inputs, and model calls that produced it."
"""
from __future__ import annotations

from datetime import datetime, timezone

from academy.agent import Agent, action

from skills import get_skill
from storage.agent_provenance import agent_tool_task

_AGENT_ID = "ReporterAgent"


class ReporterAgent(Agent):
    @action
    async def write_analysis_report(
        self,
        spec: dict,
        results: dict,
        verdict: dict,
        image_dir: str = "report_structures",
    ) -> str:
        skill = get_skill(results["analysis_type"])

        with agent_tool_task(
            _AGENT_ID, "write_analysis_report",
            used={"analysis_type": results.get("analysis_type"), "qc_flag": verdict.get("qc_flag")},
        ) as t:
            title_smiles = spec.get("smiles") or (spec.get("parent") or {}).get("smiles", "?")

            verdict_line = f"QC verdict: **{verdict['qc_flag']}**"
            if verdict.get("qc_notes"):
                verdict_line += f" ({verdict['qc_notes']})"

            lines = [
                f"# {skill.NAME} campaign report: {title_smiles}",
                "",
                f"Generated: {datetime.now(timezone.utc).isoformat()}",
                "",
                verdict_line,
                "",
            ]

            if spec.get("source_question"):
                lines += ["## Research question", spec["source_question"], ""]
            if spec.get("rationale"):
                lines += [
                    f"## Planner rationale ({spec.get('planner_model', 'LLM')})",
                    spec["rationale"],
                    "",
                ]

            lines += skill.describe(spec, results, verdict, image_dir)
            lines.append("")
            report = "\n".join(lines)

            t.end(generated={"report": report})
        return report
