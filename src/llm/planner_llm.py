"""
LLM-backed campaign planning (Phase 5) via the Anthropic Claude API.

Turns a natural-language research question into a structured campaign spec.
For the 'single_point' analysis this matches the flat shape
`agents.planner.PlannerAgent.plan_campaign` already produces (Phase 4), so
ExecutorAgent/CriticAgent/ReporterAgent handle it unchanged; for any other
registered analysis skill (src/skills/ - Phase 6, architecture plan Sec.16)
it returns the nested parent/fragments shape those skills expect.

Which analysis types Claude is even allowed to choose is pulled dynamically
from the skill registry (`skills.list_skills()`) rather than hand-listed
here - registering a new skill module in src/skills/ automatically makes it
choosable by the planner, with no prompt edit required in this file. What
IS still hand-maintained here is the input SHAPE (flat 'single_point'-like
vs nested 'parent + fragments'-like) - a future skill needing a genuinely
different shape (neither of those two) will need a new branch in
generate_campaign_spec_llm below, even though ExecutorAgent/CriticAgent/
ReporterAgent themselves need no changes regardless (see architecture plan
Sec.16 for this known boundary).

Backend decision (2026-09-15): gpt-oss:20b via local Ollama (the original
Phase 5 plan) couldn't be gotten running locally; switched to a hosted API
instead. Chose the Anthropic Claude API over OpenAI's, using Claude's
native tool-calling (forced via `tool_choice`) to get a reliable structured
spec back rather than parsing free text. Trade-off worth remembering for
Phase 8: this requires outbound network access and an API key, which
Frontier/Defiant compute nodes likely won't have - the HPC phase may need
to revisit this (e.g. call the API from a login/gateway node, or fall back
to a local model there).

Requires:
    pip install anthropic
    export ANTHROPIC_API_KEY=...   (never hardcode a key in code or files)

Model choice is configurable via the ANTHROPIC_MODEL env var. Default here
is claude-sonnet-5 (a strong balance of reasoning quality and latency/cost
for this small, structured task) - see
https://docs.claude.com/en/docs/about-claude/models for the current model
lineup if this default ever needs updating.

Phase 10 (2026-09-17) adds a sibling function,
generate_revised_campaign_spec_llm, used by PlannerAgent.replan_campaign
when a human chooses 'reprioritize' at the post-analysis approval gate
(agents/human_gate.py). It shares this module's system prompt, tool
schema, bounded self-correction retry, and output shape - the only real
difference is the initial user message, which embeds the previous
spec/results/verdict plus the human's free-text feedback instead of a
bare question string.

Bounded self-correction (found 2026-09-15, real run): a forced tool call
does not GUARANTEE the model actually fills in every field the system
prompt asks for - a real run against a cheaper/faster model alias returned
analysis_type='bde' with an empty 'fragments' list, despite the prompt
saying that's required, and despite JSON-schema hints in the tool
definition (see _build_tool_schema). Rather than crash on the first
mistake, generate_campaign_spec_llm now gives the model exactly one
corrective turn: it sends the bad tool call back as an is_error tool_result
naming exactly what was missing, and lets the model retry once, still
forced through tool_choice. Only a second consecutive failure raises. This
mirrors the "cap recovery attempts" pitfall (architecture plan Sec.7) -
bounded at 2 total attempts, not unbounded retrying.
"""
from __future__ import annotations

import json
import os
from typing import Any

import anthropic

from storage.flowcept_store import flowcept_task
from storage.agent_provenance import agent_llm_invocation_task
from skills import list_skills

_DEFAULT_MODEL = "claude-sonnet-5"
_MAX_ATTEMPTS = 2


def _build_system_prompt() -> str:
    skill_lines = "\n".join(f"- '{name}': {desc}" for name, desc in list_skills())
    return f"""You are a computational chemistry campaign planner. Given a
plain-language research question, decide which registered analysis type it
needs:

{skill_lines}

CRITICAL for any analysis that needs 'fragments' (e.g. a bond dissociation
energy): a fragment is what you get from breaking exactly ONE bond of the
parent, with every OTHER bond of the parent left intact. Never fully
atomize the molecule (breaking every bond) when only one bond/transformation
was asked about.

Worked example - the O-H bond in water (H-O-H): breaking one O-H bond
leaves the other O-H bond intact, giving a hydroxyl radical (SMILES
'[OH]', charge 0, multiplicity 2) and a hydrogen atom (SMILES '[H]',
charge 0, multiplicity 2). The WRONG answer (a common mistake) is a bare
oxygen atom plus two separate hydrogen atoms ('[O]' + '[H]' + '[H]') -
that is full atomization (breaking BOTH O-H bonds), not the single-bond
homolysis that was asked about, and also is not what "the parent minus one
bond" actually looks like structurally. Apply the same logic to any other
molecule or transformation: identify exactly what changes, then write the
SMILES for exactly what's left on each side.

For every molecule/species you specify, determine:
- its SMILES string
- its net charge and spin multiplicity (neutral closed-shell molecules are
  charge 0, multiplicity 1; a radical fragment from homolytic bond cleavage
  is typically multiplicity 2; adjust charge/multiplicity for ions/other
  radical types as chemically appropriate). The number of valence electrons
  implied by your SMILES/charge must be able to support your multiplicity
  (multiplicity 2 needs an ODD valence-electron count) - double check this
  arithmetic before answering, since an inconsistent combination will make
  the downstream calculation fail outright.

Also choose a reasonable DFT basis set and exchange-correlation functional
for a refinement step run with PySCF, after a cheap GFN2-xTB prescreen for
each species.

Always call the emit_campaign_spec tool with your answer - never answer in
plain text. If the chosen analysis type needs more than the parent molecule,
you MUST include 'fragments' as a non-empty array (at least one entry) in
that SAME tool call - a bond dissociation energy (or any other
fragments-based analysis) with an empty or missing 'fragments' list cannot
be computed at all and will be rejected. If anything is ambiguous, make the
most chemically reasonable assumption and say so in `rationale`."""


def _build_tool_schema() -> dict:
    analysis_types = [name for name, _ in list_skills()]
    return {
        "name": "emit_campaign_spec",
        "description": (
            "Emit a structured computational chemistry campaign spec. "
            "'analysis_type' must be one of the registered analysis types "
            "(see the system prompt for what each one means). For an "
            "analysis type that needs more than the parent molecule (e.g. "
            "a bond dissociation energy), also provide 'fragments': the "
            "additional species that analysis needs. Never emit "
            "analysis_type='bde' (or any other fragments-based type) with "
            "an empty 'fragments' list."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "analysis_type": {
                    "type": "string",
                    "enum": analysis_types,
                    "description": "Which registered analysis type this question needs.",
                },
                "parent": {
                    "type": "object",
                    "description": "The main/reference molecule for this analysis.",
                    "properties": {
                        "smiles": {"type": "string"},
                        "charge": {"type": "integer"},
                        "multiplicity": {"type": "integer"},
                    },
                    "required": ["smiles", "charge", "multiplicity"],
                },
                "fragments": {
                    "type": "array",
                    "description": (
                        "Additional related species this analysis needs beyond "
                        "the parent (e.g. dissociation fragments for a bond "
                        "dissociation energy, oxidized/reduced forms for a redox "
                        "potential, product/reactant states for Marcus theory). "
                        "Leave empty ONLY for an analysis_type that truly needs "
                        "nothing beyond the parent (e.g. 'single_point') - for "
                        "'bde' and similar multi-species analyses this must be "
                        "a non-empty array."
                    ),
                    "items": {
                        "type": "object",
                        "properties": {
                            "label": {
                                "type": "string",
                                "description": "Short human-readable label, e.g. 'OH radical'",
                            },
                            "smiles": {"type": "string"},
                            "charge": {"type": "integer"},
                            "multiplicity": {"type": "integer"},
                        },
                        "required": ["smiles", "charge", "multiplicity"],
                    },
                },
                "basis": {
                    "type": "string",
                    "description": "PySCF-compatible basis set name, e.g. '6-31g*'",
                },
                "xc": {
                    "type": "string",
                    "description": "PySCF-compatible exchange-correlation functional, e.g. 'b3lyp'",
                },
                "rationale": {
                    "type": "string",
                    "description": "Brief explanation of the choices above",
                },
            },
            "required": ["analysis_type", "parent", "basis", "xc", "rationale"],
            "if": {"properties": {"analysis_type": {"const": "bde"}}},
            "then": {"required": ["fragments"]},
        },
    }


def _require_species_dict(obj: Any, what: str) -> dict:
    """Validate that `obj` is shaped like {'smiles': ..., 'charge': ...,
    'multiplicity': ...}. Bug found 2026-09-16 (Sec.18): the tool schema
    marks 'parent' (and each fragment) as an object, but Claude's forced
    tool_choice does not hard-enforce nested shapes - a real run returned
    analysis_type='bde' with 'parent' as a bare SMILES string instead of an
    object, crashing deep in dict-key access with an opaque TypeError.
    Raising ValueError here instead routes the mistake through the same
    bounded self-correction retry in generate_campaign_spec_llm already
    used for the missing-fragments case, rather than crashing the run."""
    if not isinstance(obj, dict) or "smiles" not in obj:
        raise ValueError(
            f"{what} must be an object with a 'smiles' field (and "
            f"optionally 'charge'/'multiplicity'), got {obj!r}"
        )
    return obj


def _spec_from_raw(raw: dict[str, Any], question: str, resolved_model: str) -> dict:
    """Build the final campaign-spec dict from one validated tool-call
    payload. Raises ValueError if a fragments-requiring analysis_type came
    back with no fragments, or if 'parent'/a fragment entry isn't shaped as
    the expected object - the caller (generate_campaign_spec_llm) decides
    whether that's worth a corrective retry."""
    analysis_type = raw.get("analysis_type", "single_point")
    parent = _require_species_dict(raw.get("parent"), "'parent'")

    if analysis_type == "single_point":
        return {
            "analysis_type": "single_point",
            "smiles": parent["smiles"],
            "charge": int(parent.get("charge", 0)),
            "multiplicity": int(parent.get("multiplicity", 1)),
            "method_tiers": ["xtb-GFN2", f"pyscf-{raw['xc']}/{raw['basis']}"],
            "basis": raw["basis"],
            "xc": raw["xc"],
            "rationale": raw.get("rationale"),
            "source_question": question,
            "planner_model": resolved_model,
        }

    # Every other currently-registered analysis type (just 'bde' for now)
    # uses the nested parent/fragments shape. See module docstring: a
    # future skill needing a genuinely different shape needs its own branch
    # here.
    fragments_raw = raw.get("fragments") or []
    if not fragments_raw:
        raise ValueError(
            f"analysis_type={analysis_type!r} but no fragments were provided"
        )
    fragments = [
        _require_species_dict(frag, f"fragments[{i}]")
        for i, frag in enumerate(fragments_raw)
    ]
    return {
        "analysis_type": analysis_type,
        "parent": {
            "smiles": parent["smiles"],
            "charge": int(parent.get("charge", 0)),
            "multiplicity": int(parent.get("multiplicity", 1)),
        },
        "fragments": [
            {
                "label": frag.get("label", frag["smiles"]),
                "smiles": frag["smiles"],
                "charge": int(frag.get("charge", 0)),
                "multiplicity": int(frag.get("multiplicity", 1)),
            }
            for frag in fragments
        ],
        "basis": raw["basis"],
        "xc": raw["xc"],
        "rationale": raw.get("rationale"),
        "source_question": question,
        "planner_model": resolved_model,
    }


@flowcept_task
def generate_campaign_spec_llm(question: str, model: str | None = None) -> dict:
    """Ask Claude to turn a natural-language question into a campaign spec.

    Returns a JSON-safe dict. For analysis_type='single_point' this matches
    the flat shape PlannerAgent.plan_campaign already produces (plus extra
    provenance fields: 'rationale', 'source_question', 'planner_model'), so
    ExecutorAgent.run_analysis/CriticAgent.review_analysis/
    ReporterAgent.write_analysis_report handle it unchanged. For any other
    analysis type it instead returns the nested parent/fragments shape
    those skills expect (see module docstring for the known boundary here).

    Gives the model one bounded corrective retry (see module docstring) if
    it chooses a fragments-requiring analysis_type but returns no fragments;
    raises ValueError only if that still fails on the second attempt."""
    resolved_model = model or os.environ.get("ANTHROPIC_MODEL", _DEFAULT_MODEL)
    client = anthropic.Anthropic()  # reads ANTHROPIC_API_KEY from the environment

    system_prompt = _build_system_prompt()
    tools = [_build_tool_schema()]
    messages: list[dict[str, Any]] = [{"role": "user", "content": question}]

    last_error: Exception | None = None
    for attempt in range(1, _MAX_ATTEMPTS + 1):
        # One PROV-AGENT 'AIModelInvocation' record per actual API call -
        # architecture plan Sec.18. A bounded self-correction retry (below)
        # shows up as a second, separate record here, making the retry
        # itself part of the auditable trail rather than an invisible
        # implementation detail.
        with agent_llm_invocation_task(
            "PlannerAgent", "generate_campaign_spec_llm",
            used={"question": question, "model": resolved_model, "attempt": attempt},
        ) as t:
            response = client.messages.create(
                model=resolved_model,
                max_tokens=1024,
                system=system_prompt,
                tools=tools,
                tool_choice={"type": "tool", "name": "emit_campaign_spec"},
                messages=messages,
            )

            tool_use = next(block for block in response.content if block.type == "tool_use")
            raw: dict[str, Any] = dict(tool_use.input)
            t.end(generated={"tool_input": raw})

        try:
            return _spec_from_raw(raw, question, resolved_model)
        except ValueError as exc:
            last_error = exc
            if attempt >= _MAX_ATTEMPTS:
                break
            # Bounded self-correction: hand the mistake back as a tool
            # error and let the model try again once, still forced through
            # tool_choice (architecture plan Sec.7/Sec.17 - cap recovery
            # attempts rather than retry unboundedly).
            messages = messages + [
                {"role": "assistant", "content": response.content},
                {
                    "role": "user",
                    "content": [
                        {
                            "type": "tool_result",
                            "tool_use_id": tool_use.id,
                            "content": (
                                f"Invalid emit_campaign_spec call: {exc}. Call "
                                "emit_campaign_spec again, this time including a "
                                "non-empty 'fragments' array."
                            ),
                            "is_error": True,
                        }
                    ],
                },
            ]

    raise ValueError(
        f"Claude failed to produce a valid campaign spec after {_MAX_ATTEMPTS} "
        f"attempts, for question: {question!r} ({last_error})"
    )



def _build_revision_user_message(
    previous_spec: dict, previous_results: dict, verdict: dict, feedback: str
) -> str:
    """Build the initial user message for a replan request: the same
    information a human reviewing this campaign would have in front of
    them - what was planned, what happened, what the critic said, and what
    the human wants changed."""
    return (
        "You previously planned this computational chemistry campaign:\n"
        f"{json.dumps(previous_spec, indent=2, default=str)}\n\n"
        "It was executed, producing these results:\n"
        f"{json.dumps(previous_results, indent=2, default=str)}\n\n"
        "The critic's review of those results:\n"
        f"{json.dumps(verdict, indent=2, default=str)}\n\n"
        "A human reviewed this and gave the following feedback on what "
        "should change for a revised attempt:\n"
        f"{feedback}\n\n"
        "Call emit_campaign_spec again with a REVISED campaign spec that "
        "addresses this feedback. You may change analysis_type, method "
        "tiers, basis/xc, charge/multiplicity, or the parent/fragments "
        "themselves if the feedback calls for it - do not just repeat the "
        "previous spec unchanged unless the feedback genuinely calls for "
        "no changes."
    )


@flowcept_task
def generate_revised_campaign_spec_llm(
    previous_spec: dict,
    previous_results: dict,
    verdict: dict,
    feedback: str,
    model: str | None = None,
) -> dict:
    """Ask Claude to revise a campaign spec in light of its own prior
    results, the critic's verdict, and free-text human feedback (Phase 10 -
    the 'reprioritize' choice at the post-analysis human-approval gate,
    agents/human_gate.py). Shares generate_campaign_spec_llm's tool schema,
    bounded self-correction retry (see module docstring), and output shape;
    the only difference is what goes into the initial user message (see
    _build_revision_user_message above) and that two extra provenance
    fields, 'revised_from' and 'revision_feedback', are stamped onto the
    returned spec."""
    resolved_model = model or os.environ.get("ANTHROPIC_MODEL", _DEFAULT_MODEL)
    client = anthropic.Anthropic()

    system_prompt = _build_system_prompt()
    tools = [_build_tool_schema()]
    user_message = _build_revision_user_message(previous_spec, previous_results, verdict, feedback)
    messages: list[dict[str, Any]] = [{"role": "user", "content": user_message}]
    question_for_spec = (
        f"{previous_spec.get('source_question', previous_spec.get('smiles', ''))} "
        f"| REVISION FEEDBACK: {feedback}"
    )

    last_error: Exception | None = None
    for attempt in range(1, _MAX_ATTEMPTS + 1):
        with agent_llm_invocation_task(
            "PlannerAgent", "generate_revised_campaign_spec_llm",
            used={"feedback": feedback, "model": resolved_model, "attempt": attempt},
        ) as t:
            response = client.messages.create(
                model=resolved_model,
                max_tokens=1024,
                system=system_prompt,
                tools=tools,
                tool_choice={"type": "tool", "name": "emit_campaign_spec"},
                messages=messages,
            )
            tool_use = next(block for block in response.content if block.type == "tool_use")
            raw: dict[str, Any] = dict(tool_use.input)
            t.end(generated={"tool_input": raw})

        try:
            spec = _spec_from_raw(raw, question_for_spec, resolved_model)
            spec["revised_from"] = previous_spec.get("analysis_type")
            spec["revision_feedback"] = feedback
            return spec
        except ValueError as exc:
            last_error = exc
            if attempt >= _MAX_ATTEMPTS:
                break
            messages = messages + [
                {"role": "assistant", "content": response.content},
                {
                    "role": "user",
                    "content": [
                        {
                            "type": "tool_result",
                            "tool_use_id": tool_use.id,
                            "content": (
                                f"Invalid emit_campaign_spec call: {exc}. Call "
                                "emit_campaign_spec again, this time including a "
                                "non-empty 'fragments' array."
                            ),
                            "is_error": True,
                        }
                    ],
                },
            ]

    raise ValueError(
        f"Claude failed to produce a valid revised campaign spec after "
        f"{_MAX_ATTEMPTS} attempts ({last_error})"
    )
