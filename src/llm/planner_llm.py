"""
LLM-backed campaign planning (Phase 5) via the Anthropic Claude API, with
two later-added alternate backends: a local-Ollama backend and a
Claude-Code-CLI backend (see "Backend decision" below for both).

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
_spec_from_raw below, even though ExecutorAgent/CriticAgent/ReporterAgent
themselves need no changes regardless (see architecture plan Sec.16 for this
known boundary).

Backend decision (2026-09-15): gpt-oss:20b via local Ollama (the original
Phase 5 plan) couldn't be gotten running locally; switched to a hosted API
instead. Chose the Anthropic Claude API over OpenAI's, using Claude's
native tool-calling (forced via `tool_choice`) to get a reliable structured
spec back rather than parsing free text. Trade-off worth remembering for
Phase 8: this requires outbound network access and an API key, which
Frontier/Defiant compute nodes likely won't have - the HPC phase may need
to revisit this (e.g. call the API from a login/gateway node, or fall back
to a local model there).

Ollama backend (later addition): re-added as an opt-in `backend="ollama"`
path (CLI: `--ollama`, agents/run_agentic_campaign.py), for offline/no-API-
key use with whatever model the caller already has pulled locally. Unlike
Anthropic, Ollama's tool-calling API has no forced `tool_choice` - the model
decides whether to call a tool at all, so a model that ignores the
instruction to always call `emit_campaign_spec` surfaces as a ValueError
("did not call the tool") rather than silently falling back to free text.
Pick a model documented as tool-call-capable (e.g. llama3.1, qwen2.5) and
expect it to need the same bounded self-correction retries, or more, that
Claude occasionally needs (see below) - smaller local models are more prone
to malformed tool calls, not less.

Claude-Code-CLI backend (later addition): a third opt-in `backend=
"claude_cli"` path (CLI: `--claude-cli`), for a user with Claude Code
installed and logged in but no separate ANTHROPIC_API_KEY - it shells out
to `claude -p` (the CLI's own headless/print mode), which authenticates
with the CLI's existing login (subscription/OAuth) rather than reading
ANTHROPIC_API_KEY from the environment. `claude -p` has no native
tool-calling of its own (no tool_choice, no tools= param), so the tool
schema is instead spelled out as plain-text instructions
(_json_output_instructions) and the reply is parsed as raw JSON
(_extract_json_object) - a real test call came back wrapped in a markdown
code fence despite being told not to, so that gets stripped rather than
failing. `--allowedTools ""` (an empty whitelist) keeps the spawned CLI
session from touching the filesystem or running anything - it's used
purely as a text-completion engine here, not as an agent. Conversation
state for the bounded retry lives server-side via `--resume <session_id>`
rather than replayed messages, since that's how `claude -p` sessions work.
Worth knowing before reaching for this backend: it still counts against
whatever Claude Code plan/usage the CLI is logged into (not free just
because no API key is set), and each call spawns a real (if tool-less)
Claude Code session, so it's slower than a direct API call.

Requires (Anthropic backend, the default):
    pip install anthropic
    export ANTHROPIC_API_KEY=...   (never hardcode a key in code or files)

Requires (Ollama backend, `backend="ollama"` / CLI `--ollama`):
    pip install ollama
    A running local `ollama serve`, with the chosen model already pulled
    (`ollama pull llama3.1` or similar).

Requires (Claude-Code-CLI backend, `backend="claude_cli"` / CLI
`--claude-cli`):
    Claude Code itself installed and logged in (`claude` on PATH,
    `claude /login` already run) - no ANTHROPIC_API_KEY needed.

Model choice is configurable via the ANTHROPIC_MODEL env var (Anthropic or
Claude-Code-CLI backend - both select a Claude model by the same names) or
the OLLAMA_MODEL env var (Ollama backend), or the `--model` CLI flag either
way. Default here is claude-sonnet-5 for Anthropic/Claude-Code-CLI (a
strong balance of reasoning quality and latency/cost for this small,
structured task - see https://docs.claude.com/en/docs/about-claude/models
for the current model lineup if this default ever needs updating) and
llama3.1 for Ollama (a widely-available tool-calling-capable default; swap
for whatever's already pulled locally).

Phase 10 (2026-09-17) adds a sibling function,
generate_revised_campaign_spec_llm, used by PlannerAgent.replan_campaign
when a human chooses 'reprioritize' at the post-analysis approval gate
(agents/human_gate.py). It shares this module's system prompt, tool
schema, bounded self-correction retry, and output shape - the only real
difference is the initial user message, which embeds the previous
spec/results/verdict plus the human's free-text feedback instead of a
bare question string. Both functions delegate their actual backend call
and retry loop to the shared `_generate_spec_via_llm` below.

Bounded self-correction (found 2026-09-15, real run): a forced tool call
does not GUARANTEE the model actually fills in every field the system
prompt asks for - a real run against a cheaper/faster model alias returned
analysis_type='bde' with an empty 'fragments' list, despite the prompt
saying that's required, and despite JSON-schema hints in the tool
definition (see _build_tool_schema). Rather than crash on the first
mistake, _generate_spec_via_llm now gives the model exactly one corrective
turn: it sends the bad tool call back as a tool-error message and lets the
model retry once, still asked to call the same tool. Only a second
consecutive failure raises. This mirrors the "cap recovery attempts"
pitfall (architecture plan Sec.7) - bounded at 2 total attempts, not
unbounded retrying.
"""
from __future__ import annotations

import json
import os
import re
from typing import Any

import anthropic

from storage.flowcept_store import flowcept_task
from storage.agent_provenance import agent_llm_invocation_task
from skills import list_skills

_DEFAULT_MODEL = "claude-sonnet-5"
_DEFAULT_OLLAMA_MODEL = "llama3.1"
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
    bounded self-correction retry in _generate_spec_via_llm already used
    for the missing-fragments case, rather than crashing the run."""
    if not isinstance(obj, dict) or "smiles" not in obj:
        raise ValueError(
            f"{what} must be an object with a 'smiles' field (and "
            f"optionally 'charge'/'multiplicity'), got {obj!r}"
        )
    return obj


def _spec_from_raw(raw: dict[str, Any], question: str, resolved_model: str) -> dict:
    """Build the final campaign-spec dict from one validated tool-call
    payload. Raises ValueError if a fragments-requiring analysis_type came
    back with no fragments, if 'parent'/a fragment entry isn't shaped as
    the expected object, or if 'basis'/'xc' are missing - the caller
    (_generate_spec_via_llm) decides whether that's worth a corrective
    retry.

    Bug found 2026-09-24 (Ollama backend): the tool schema marks 'basis'
    and 'xc' required, but unlike Anthropic's forced tool_choice, Ollama's
    tool-calling has no hard schema enforcement - a real run against a
    local model returned a valid-looking BDE tool call with both fields
    simply absent. The original code read raw['basis']/raw['xc'] directly,
    which raised a bare KeyError that bypassed the bounded self-correction
    retry below (it only catches ValueError) and crashed the whole call.
    Checking for them explicitly here, before either return branch, routes
    that failure through the same retry as every other malformed-tool-call
    case instead."""
    analysis_type = raw.get("analysis_type", "single_point")
    parent = _require_species_dict(raw.get("parent"), "'parent'")
    missing = [field for field in ("basis", "xc") if not raw.get(field)]
    if missing:
        raise ValueError(f"missing required field(s): {', '.join(missing)}")

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


def _resolve_model(backend: str, model: str | None) -> str:
    if backend == "ollama":
        return model or os.environ.get("OLLAMA_MODEL", _DEFAULT_OLLAMA_MODEL)
    return model or os.environ.get("ANTHROPIC_MODEL", _DEFAULT_MODEL)


def _anthropic_attempt(
    model: str, system_prompt: str, tools: list[dict], messages: list[dict[str, Any]]
) -> tuple[dict, Any, Any]:
    """One Claude tool-call attempt. Returns (raw tool input, the
    assistant-turn entry to replay into history on a retry, an opaque
    'retry ref' - here the tool_use_id the retry's tool_result must echo)."""
    client = anthropic.Anthropic()  # reads ANTHROPIC_API_KEY from the environment
    response = client.messages.create(
        model=model,
        max_tokens=1024,
        system=system_prompt,
        tools=tools,
        tool_choice={"type": "tool", "name": "emit_campaign_spec"},
        messages=messages,
    )
    tool_use = next(block for block in response.content if block.type == "tool_use")
    raw: dict[str, Any] = dict(tool_use.input)
    return raw, {"role": "assistant", "content": response.content}, tool_use.id


def _anthropic_retry_messages(
    messages: list[dict[str, Any]], assistant_entry: Any, retry_ref: Any, error_text: str
) -> list[dict[str, Any]]:
    return messages + [
        assistant_entry,
        {
            "role": "user",
            "content": [
                {
                    "type": "tool_result",
                    "tool_use_id": retry_ref,
                    "content": error_text,
                    "is_error": True,
                }
            ],
        },
    ]


def _ollama_tool_schema(anthropic_tool: dict) -> dict:
    """Ollama's `tools=` param wants OpenAI-shaped function tools
    ({'type': 'function', 'function': {'name', 'description',
    'parameters'}}), not Anthropic's {'name', 'description',
    'input_schema'} shape. The JSON-schema body itself (input_schema /
    parameters) is reused as-is - only the wrapper differs."""
    return {
        "type": "function",
        "function": {
            "name": anthropic_tool["name"],
            "description": anthropic_tool["description"],
            "parameters": anthropic_tool["input_schema"],
        },
    }


def _ollama_attempt(
    model: str, system_prompt: str, tools: list[dict], messages: list[dict[str, Any]]
) -> tuple[dict, Any, Any]:
    """One local-Ollama tool-call attempt. Unlike Anthropic, Ollama has no
    forced tool_choice - the model decides whether to call a tool at all,
    so a model that ignores the system prompt's "always call the tool"
    instruction surfaces here as a ValueError (routed through the same
    bounded retry as a malformed tool call), not a silent free-text
    fallback."""
    import ollama  # local import: optional dependency, only needed for this backend

    full_messages = [{"role": "system", "content": system_prompt}] + messages
    response = ollama.chat(
        model=model,
        messages=full_messages,
        tools=[_ollama_tool_schema(t) for t in tools],
    )
    msg = response.message
    if not msg.tool_calls:
        raise ValueError(
            "local Ollama model did not call the emit_campaign_spec tool "
            f"(returned plain text instead: {msg.content!r}); this model may "
            "not support reliable tool calling - try a model documented as "
            "tool-capable (e.g. llama3.1, qwen2.5)"
        )
    call = msg.tool_calls[0]
    raw: dict[str, Any] = dict(call.function.arguments)
    return raw, msg, call.function.name


def _ollama_retry_messages(
    messages: list[dict[str, Any]], assistant_entry: Any, retry_ref: Any, error_text: str
) -> list[dict[str, Any]]:
    return messages + [
        assistant_entry,
        {"role": "tool", "tool_name": retry_ref, "content": error_text},
    ]


def _extract_json_object(text: str) -> dict:
    """Pull the campaign-spec JSON object out of a `claude -p` text reply.

    Despite being told to respond with ONLY the JSON object, real runs have
    done two different things `json.loads` alone can't handle: wrapped the
    object in a markdown code fence, and/or prefixed it with an explanatory
    sentence before the fence (e.g. "No `emit_campaign_spec` tool exists in
    this environment, so per the schema instructions I'm responding with
    the JSON object directly.\\n\\n```json\\n{...}"). A leading-fence-only
    strip (the original approach here) misses that second case entirely -
    this instead searches the whole text for a fenced block (fence can
    start anywhere, not just at position 0), then falls back to the first
    '{' in whatever's left; `json.JSONDecoder.raw_decode` is used instead
    of a naive `json.loads` so trailing prose AFTER the object (also seen
    in practice) doesn't fail the parse either."""
    fence_match = re.search(r"```(?:json)?\s*\n(.*?)```", text, re.DOTALL)
    candidate = fence_match.group(1) if fence_match else text

    start = candidate.find("{")
    if start == -1:
        raise ValueError(
            f"claude -p did not return any JSON object; raw response: {text[:500]!r}"
        )
    try:
        parsed, _ = json.JSONDecoder().raw_decode(candidate[start:])
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"claude -p did not return valid JSON ({exc}); raw response: {text[:500]!r}"
        )
    if not isinstance(parsed, dict):
        raise ValueError(f"claude -p returned valid JSON but not an object: {parsed!r}")
    return parsed


def _json_output_instructions(anthropic_tool: dict) -> str:
    """`claude -p` (Claude Code's headless CLI mode) has no native forced
    tool-calling - unlike the Anthropic Messages API's tool_choice or
    Ollama's tools= param, there's no way to make the CLI guarantee a
    structured call. Instead the schema is spelled out as plain
    instructions appended to the system prompt, and the reply is parsed as
    raw JSON (_extract_json_object) - the bounded self-correction retry
    that already handles the other two backends' failure modes also covers
    "didn't return valid JSON" and "wrapped it in a markdown fence anyway"
    here."""
    return (
        "\n\nRespond with ONLY a single JSON object - no markdown code "
        "fences, no prose before or after it - matching exactly this JSON "
        f"Schema:\n\n{json.dumps(anthropic_tool['input_schema'], indent=2)}\n\n"
        f"The object's top-level fields must match the schema above for "
        f"the '{anthropic_tool['name']}' call."
    )


def _claude_cli_attempt(
    model: str, system_prompt: str, tools: list[dict], messages: list[dict[str, Any]]
) -> tuple[dict, Any, Any]:
    """One `claude -p` (Claude Code's own CLI, headless/print mode)
    attempt. Uses the CLI's own logged-in credentials (subscription/OAuth)
    instead of a raw ANTHROPIC_API_KEY - for a user who has Claude Code
    installed and logged in but no separate API key. `--allowedTools ""`
    is an empty tool whitelist, so the spawned session can't touch the
    filesystem, run commands, or do anything besides answer in text - it's
    used purely as a text-completion engine here, not as an agent.

    Conversation state lives server-side (`--resume <session_id>`), unlike
    the other two backends which replay the full message history
    themselves - so `messages` here only ever needs its last entry (the
    next prompt to send, and the session_id to resume, if any)."""
    import shutil
    import subprocess

    claude_bin = shutil.which("claude")
    if claude_bin is None:
        raise RuntimeError(
            "claude CLI not found on PATH - the --claude-cli backend needs "
            "Claude Code itself installed and logged in "
            "(see https://claude.com/claude-code)"
        )

    last = messages[-1]
    prompt = last["content"]
    session_id = last.get("session_id")

    cmd = [
        claude_bin, "-p", prompt,
        "--output-format", "json",
        "--model", model,
        "--allowedTools", "",
    ]
    if session_id:
        cmd += ["--resume", session_id]
    else:
        # System prompt is only settable on the first turn of a session -
        # --resume replays the recorded one automatically (system-prompt
        # snapshotting, on by default), so it's omitted on retries.
        cmd += ["--system-prompt", system_prompt + _json_output_instructions(tools[0])]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        raise RuntimeError(f"claude -p exited {result.returncode}: {result.stderr.strip()[:1000]}")

    payload = json.loads(result.stdout)
    if payload.get("is_error"):
        raise ValueError(f"claude -p returned an error: {payload.get('result')!r}")

    raw = _extract_json_object(payload["result"])
    return raw, None, payload["session_id"]


def _claude_cli_retry_messages(
    messages: list[dict[str, Any]], assistant_entry: Any, retry_ref: Any, error_text: str
) -> list[dict[str, Any]]:
    # retry_ref is the previous attempt's session_id; --resume in
    # _claude_cli_attempt continues that same server-side conversation, so
    # (unlike the other two backends) no prior turns need replaying here.
    return [{"content": error_text, "session_id": retry_ref}]


_BACKENDS = {
    "anthropic": (_anthropic_attempt, _anthropic_retry_messages),
    "ollama": (_ollama_attempt, _ollama_retry_messages),
    "claude_cli": (_claude_cli_attempt, _claude_cli_retry_messages),
}


def _generate_spec_via_llm(
    *,
    question_for_spec: str,
    initial_messages: list[dict[str, Any]],
    model: str | None,
    backend: str,
    provenance_action: str,
) -> dict:
    """Shared attempt/retry loop behind both generate_campaign_spec_llm and
    generate_revised_campaign_spec_llm: resolve the model, dispatch to the
    chosen backend's one-attempt call, validate the result via
    _spec_from_raw, and give the model one bounded corrective retry (see
    module docstring) before raising."""
    if backend not in _BACKENDS:
        raise ValueError(f"Unknown LLM backend {backend!r}; expected one of {sorted(_BACKENDS)}")
    attempt_fn, retry_fn = _BACKENDS[backend]

    resolved_model = _resolve_model(backend, model)
    system_prompt = _build_system_prompt()
    tools = [_build_tool_schema()]
    messages = list(initial_messages)

    last_error: Exception | None = None
    for attempt in range(1, _MAX_ATTEMPTS + 1):
        # One PROV-AGENT 'AIModelInvocation' record per actual LLM call -
        # architecture plan Sec.18. A bounded self-correction retry (below)
        # shows up as a second, separate record here, making the retry
        # itself part of the auditable trail rather than an invisible
        # implementation detail.
        with agent_llm_invocation_task(
            "PlannerAgent", provenance_action,
            used={
                "question": question_for_spec,
                "model": resolved_model,
                "backend": backend,
                "attempt": attempt,
            },
        ) as t:
            raw, assistant_entry, retry_ref = attempt_fn(resolved_model, system_prompt, tools, messages)
            t.end(generated={"tool_input": raw})

        try:
            return _spec_from_raw(raw, question_for_spec, resolved_model)
        except ValueError as exc:
            last_error = exc
            if attempt >= _MAX_ATTEMPTS:
                break
            # Bounded self-correction: hand the mistake back as a tool
            # error and let the model try again once (architecture plan
            # Sec.7/Sec.17 - cap recovery attempts rather than retry
            # unboundedly). The message names the ACTUAL problem (missing
            # fragments, missing basis/xc, malformed parent/fragment, ...)
            # rather than a fixed "add fragments" hint that would mislead
            # the retry for any other failure (found 2026-09-24: a local
            # Ollama model tripped the missing-basis/xc case, not fragments).
            error_text = (
                f"Invalid emit_campaign_spec call: {exc}. Call "
                "emit_campaign_spec again, fixing exactly that problem and "
                "keeping everything else about your answer the same."
            )
            messages = retry_fn(messages, assistant_entry, retry_ref, error_text)

    raise ValueError(
        f"{backend} backend failed to produce a valid campaign spec after "
        f"{_MAX_ATTEMPTS} attempts, for question: {question_for_spec!r} ({last_error})"
    )


@flowcept_task
def generate_campaign_spec_llm(
    question: str, model: str | None = None, backend: str = "anthropic"
) -> dict:
    """Ask an LLM (Claude via the Anthropic API by default, a local Ollama
    model if backend='ollama', or Claude via the logged-in Claude Code CLI
    if backend='claude_cli') to turn a natural-language question into a campaign
    spec.

    Returns a JSON-safe dict. For analysis_type='single_point' this matches
    the flat shape PlannerAgent.plan_campaign already produces (plus extra
    provenance fields: 'rationale', 'source_question', 'planner_model'), so
    ExecutorAgent.run_analysis/CriticAgent.review_analysis/
    ReporterAgent.write_analysis_report handle it unchanged. For any other
    analysis type it instead returns the nested parent/fragments shape
    those skills expect (see module docstring for the known boundary here).

    Gives the model one bounded corrective retry (see module docstring) if
    it chooses a fragments-requiring analysis_type but returns no fragments,
    or (Ollama only) doesn't call the tool at all; raises ValueError only if
    that still fails on the second attempt."""
    return _generate_spec_via_llm(
        question_for_spec=question,
        initial_messages=[{"role": "user", "content": question}],
        model=model,
        backend=backend,
        provenance_action="generate_campaign_spec_llm",
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
    backend: str = "anthropic",
) -> dict:
    """Ask an LLM (Claude via the Anthropic API by default, a local Ollama
    model if backend='ollama', or Claude via the logged-in Claude Code CLI
    if backend='claude_cli') to revise a campaign spec in light of its own prior
    results, the critic's verdict, and free-text human feedback (Phase 10 -
    the 'reprioritize' choice at the post-analysis human-approval gate,
    agents/human_gate.py). Shares generate_campaign_spec_llm's tool schema,
    backend dispatch, and bounded self-correction retry (see module
    docstring) via _generate_spec_via_llm; the only difference is what goes
    into the initial user message (see _build_revision_user_message above)
    and that two extra provenance fields, 'revised_from' and
    'revision_feedback', are stamped onto the returned spec."""
    question_for_spec = (
        f"{previous_spec.get('source_question', previous_spec.get('smiles', ''))} "
        f"| REVISION FEEDBACK: {feedback}"
    )
    user_message = _build_revision_user_message(previous_spec, previous_results, verdict, feedback)

    spec = _generate_spec_via_llm(
        question_for_spec=question_for_spec,
        initial_messages=[{"role": "user", "content": user_message}],
        model=model,
        backend=backend,
        provenance_action="generate_revised_campaign_spec_llm",
    )
    spec["revised_from"] = previous_spec.get("analysis_type")
    spec["revision_feedback"] = feedback
    return spec
