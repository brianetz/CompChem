"""
LLM half of the BDE plan critic: the judgment-based review in
src/bde-plan-critic.md, run by the same model/backend as the planner after
the deterministic validator (analysis/bde_plan_validator.py) has passed.

Enabled by COMPCHEM_PLAN_CRITIC_SKILL=<path to the skill .md>; off when
unset. The skill text is the critic's system prompt as written; the
validator report is handed over so the model judges only what code cannot
(rationale completeness, bond choice against the question, spin-state
justification) and never re-derives - or overrules - a code result.
"""
from __future__ import annotations

import json
import os
from typing import Any

from llm.plan_checks import _decode
from storage.agent_provenance import agent_llm_invocation_task

_VERDICTS = ("PASS", "WARN", "FAIL")

_TOOL = {
    "name": "emit_plan_review",
    "description": "Return the plan review: verdict, failed checks with exact fixes, and warnings.",
    "input_schema": {
        "type": "object",
        "properties": {
            "verdict": {"type": "string", "enum": list(_VERDICTS)},
            "failed_checks": {
                "type": "array", "items": {"type": "string"},
                "description": "One string per failed check: quoted evidence, then the exact fix.",
            },
            "warnings": {
                "type": "array", "items": {"type": "string"},
                "description": "One string per warning, with a one-line suggestion.",
            },
        },
        "required": ["verdict", "failed_checks", "warnings"],
    },
}

_CALL_INSTRUCTIONS = """

---
You are the plan critic. The deterministic validator has ALREADY run; its
report is in the user message - treat it as authoritative and do not
recompute electron counts, parity, atom balance or the cleaved bond. Apply
the rest of the procedure and the rationale rubric above, then call the
emit_plan_review tool exactly once. Use FAIL only for problems that require
the plan to change (wrong bond for the question, wrong or unjustified spin
state, method inconsistent with the spec, a rationale that states wrong
numbers or implies a 298 K enthalpy). Rubric items a rationale merely omits
are WARN. Keep each entry short enough to send back to the planner."""


def skill_path() -> str | None:
    return os.environ.get("COMPCHEM_PLAN_CRITIC_SKILL") or None


def _load_skill(path: str) -> str:
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def _user_message(spec: dict, question: str, validator: dict) -> str:
    shown = {k: spec.get(k) for k in ("analysis_type", "parent", "fragments", "basis", "xc", "rationale")}
    report = {
        "verdict": validator["verdict"],
        "warnings": validator["warns"],
        "cleaved_bond": validator.get("cleaved_bond"),
        "symmetry_unique_bonds": [
            {"bond": c["bond"], "fragments": c["fragments"]} for c in validator.get("candidates", [])
        ],
    }
    return (
        f"User question:\n{question}\n\n"
        f"Planner spec:\n{json.dumps(shown, indent=1)}\n\n"
        f"Deterministic validator report:\n{json.dumps(report, indent=1)}"
    )


def _as_list(value: Any) -> list[str]:
    value = _decode(value)
    if value is None:
        return []
    if isinstance(value, str):
        return [value] if value.strip() else []
    if isinstance(value, list):
        return [v if isinstance(v, str) else json.dumps(v) for v in value]
    return [json.dumps(value)]


def review_plan(spec: dict, question: str, validator: dict, *, backend: str, model: str) -> dict:
    """One critic call. Returns {'verdict', 'failed_checks', 'warnings'};
    verdict 'ERROR' (with 'error') if the model's reply can't be used - an
    unusable review never blocks a plan the validator passed."""
    from llm.planner_llm import _BACKENDS

    attempt_fn, _ = _BACKENDS[backend]
    system_prompt = _load_skill(skill_path()) + _CALL_INSTRUCTIONS
    messages = [{"role": "user", "content": _user_message(spec, question, validator)}]
    with agent_llm_invocation_task(
        "PlanCriticAgent", "review_plan",
        used={"question": question, "model": model, "backend": backend,
              "skill": skill_path(), "validator_verdict": validator["verdict"]},
    ) as t:
        try:
            raw, _, _ = attempt_fn(model, system_prompt, [_TOOL], messages)
        except Exception as exc:  # no tool call, HTTP error, ...
            review = {"verdict": "ERROR", "failed_checks": [], "warnings": [],
                      "error": f"{type(exc).__name__}: {exc}"}
            t.end(generated=review)
            return review
        verdict = str(raw.get("verdict", "")).strip().upper()
        review = {
            "verdict": verdict if verdict in _VERDICTS else "ERROR",
            "failed_checks": _as_list(raw.get("failed_checks")),
            "warnings": _as_list(raw.get("warnings")),
        }
        if review["verdict"] == "ERROR":
            review["error"] = f"unusable verdict {raw.get('verdict')!r}"
        t.end(generated={"tool_input": raw, "review": review})
    return review
