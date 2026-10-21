"""
Plan board: a modular BDE plan built on a shared board (blackboard
pattern) by Academy agents, under a supervising moderator.

The plan is split into modules, built in dependency order:
  target     analysis type, parent molecule, and the bond being asked about
  fragments  the species that breaking that bond gives
  method     DFT functional and basis set
  rationale  the written justification

  PlanBoard        Academy agent holding the shared state: every published
                   module version, every review, and which version of each
                   module is accepted.
  Planner (LLM)    publishes one module at a time to the board, seeing the
                   question, the accepted modules and the reviews of its
                   earlier attempts at this module.
  Validator        deterministic checks of the latest version of a module;
                   posts PASS or FAIL with exact fixes to the board.
  PlanCritic (LLM) the bde-plan-critic skill as its instructions; reviews
                   each module the Validator passed (ACCEPT or REVISE), and
                   the assembled plan at the end, where it may reopen
                   modules.
  Moderator        runs the rounds: a module is accepted only when the
                   Validator passes it and the critic accepts it; a module
                   is retried up to max_attempts times; the plan passes
                   only after every module is accepted and the final whole-
                   plan review accepts it (or reopened modules are rebuilt,
                   once). A supervise control loop pauses the board when the
                   Planner repeats a Validator-rejected mistake, appends one
                   corrective instruction to the Planner's prompt, and
                   terminates if the mistake survives that edit.

A plan therefore passes only after at least five reviewed rounds (four
modules and the final review). If the critic never accepts a module the
Validator passed, the last validated version is used and marked
critic_unresolved in the spec, so a plan is never lost to wording alone.
"""
from __future__ import annotations

import asyncio
import json
import time
import logging
import os
from typing import Any

from academy.agent import Agent, action, loop
from academy.handle import Handle

from agents.plan_discussion import _EDIT_PROMPT, _INSTRUCTION_TOOL
from analysis.bde_plan_validator import critique_plan, enumerate_bde_candidates, problems_as_text
from analysis.species_checks import electron_count, parse_smiles, spin_allowed
from llm import planner_llm as P
from llm.plan_checks import _decode, _fix_species
from llm.plan_critic_llm import _as_list
from skills import list_skills
from storage.agent_provenance import agent_llm_invocation_task
from storage.rag import bond_types_in

logger = logging.getLogger(__name__)

MODULES = ("target", "fragments", "method", "rationale")
DEFAULT_SKILL = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bde-plan-critic.md")

_SPECIES = {
    "type": "object",
    "properties": {"label": {"type": "string"},
                   "smiles": {"type": "string", "description": "A SMILES string, e.g. 'CO' (methanol) or "
                              "'C[CH2]' (ethyl radical) - never a name or a formula such as 'CH3OH'"},
                   "charge": {"type": "integer"}, "multiplicity": {"type": "integer"}},
    "required": ["smiles", "charge", "multiplicity"],
}


def module_tool(module: str) -> dict:
    """The Planner's tool for publishing one module."""
    props: dict[str, Any]
    if module == "target":
        props = {"analysis_type": {"type": "string", "enum": [n for n, _ in list_skills()]},
                 "parent": _SPECIES,
                 "bond": {"type": "string", "description": "The bond the question asks about, e.g. 'O-H' or 'benzylic C-H'"}}
    elif module == "fragments":
        props = {"fragments": {"type": "array", "items": _SPECIES,
                               "description": "The species that breaking exactly that one bond gives"}}
    elif module == "method":
        props = {"xc": {"type": "string", "description": "PySCF functional, e.g. 'b3lyp'"},
                 "basis": {"type": "string", "description": "PySCF basis set, e.g. '6-31g*'"}}
    else:
        props = {"rationale": {"type": "string", "description": "Why this plan answers the question"}}
    return {"name": f"publish_{module}", "description": f"Publish the {module} module of the plan to the board.",
            "input_schema": {"type": "object", "properties": props, "required": list(props)}}


_BOARD_INSTRUCTIONS = """

---
BOARD MODE. The plan is built one module at a time on a shared board that
a Validator and a plan critic review. Ignore the instruction above to call
emit_campaign_spec: call ONLY the publish_<module> tool you are given, for
the module you are asked for, once. Keep accepted modules as they are, and
when an earlier attempt at this module was rejected, fix exactly what the
reviews say."""

_MODULE_REVIEW_TOOL = {
    "name": "emit_module_review",
    "description": "Review one module of the plan on the board.",
    "input_schema": {"type": "object", "properties": {
        "verdict": {"type": "string", "enum": ["ACCEPT", "REVISE"]},
        "issues": {"type": "array", "items": {"type": "string"},
                   "description": "What must change, each with the exact fix; empty for ACCEPT"}},
        "required": ["verdict", "issues"]},
}
_FINAL_REVIEW_TOOL = {
    "name": "emit_final_review",
    "description": "Review the assembled plan.",
    "input_schema": {"type": "object", "properties": {
        "verdict": {"type": "string", "enum": ["ACCEPT", "REVISE"]},
        "reopen": {"type": "array", "items": {"type": "string", "enum": list(MODULES)},
                   "description": "Modules that must be rebuilt; empty for ACCEPT"},
        "issues": {"type": "array", "items": {"type": "string"}}},
        "required": ["verdict", "reopen", "issues"]},
}
_MODULE_FOCUS = {
    "target": "Is this the molecule and the bond the question asks about (e.g. benzylic vs aromatic C-H), and is the parent's charge/spin right?",
    "fragments": "Are these exactly the radicals from breaking THAT bond? The Validator lists the parent's symmetry-unique splits and which one these fragments are; check it is the site the question names.",
    "method": "Is the functional/basis appropriate for open-shell main-group thermochemistry and usable for every species? Prefer a method already cached unless the question asks otherwise.",
    "rationale": "Apply the rationale rubric. REVISE only for wrong statements or rubric items that are missing.",
}


# Found 2026-10-20 (B2): the supervisor's example repeated the Planner's
# uppercase "SMILES" key and a formula ("CH3OH") as SMILES.
_EDIT_FIELD_NAMES = (" Any example you give must use the publish tool's exact field names - lowercase smiles, "
                     "charge, multiplicity - and valid SMILES (methanol is 'CO', not 'CH3OH').")


def _critic_instructions(module: str | None) -> str:
    if module is None:
        return ("\n\n---\nYou are the plan critic for a plan assembled on a shared board. Every module has "
                "already been accepted by the Validator and by you. Review the plan as a whole against the "
                "question. ACCEPT unless something must change; then list the modules to reopen. "
                "Call emit_final_review once.")
    return (f"\n\n---\nYou are the plan critic reviewing ONE module of a plan on a shared board: "
            f"'{module}'. The Validator's code checks have passed it; never re-argue them. "
            f"{_MODULE_FOCUS[module]} ACCEPT unless this module must change. Call emit_module_review once.")


def _skill_text() -> str:
    path = os.environ.get("COMPCHEM_PLAN_CRITIC_SKILL") or DEFAULT_SKILL
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def _species_problems(sp: Any, name: str) -> list[str]:
    sp = _fix_species(sp)
    if not isinstance(sp, dict) or not isinstance(sp.get("smiles"), str):
        return [f"{name} must be a JSON object with smiles, charge and multiplicity"]
    mol = parse_smiles(sp["smiles"])
    if mol is None:
        return [f"{name}: SMILES {sp['smiles']!r} is not valid"]
    n_e = electron_count(mol, int(sp.get("charge", 0)))
    if not spin_allowed(n_e, int(sp.get("multiplicity", 1))):
        return [f"{name}: {n_e} electrons cannot have multiplicity {sp.get('multiplicity')}"]
    return []


def _split_label(c: dict) -> str:
    """'[CH2]c1ccccc1 + [H] (benzylic C-H)' - a split with its site."""
    return " + ".join(c["fragments"]) + f" ({(c.get('site') + ' ') if c.get('site') else ''}{c['bond']})"


# A well-formed example of each module, for molecules other than the one
# asked about, so the Planner copies the format, not the chemistry
# (2026-10-21: 7 of 39 rejected module versions in B1/B2 were format
# errors and 20 were rationales the critic sent back).
_MODULE_EXAMPLES = {
    "target": {"analysis_type": "bde", "parent": {"label": "water", "smiles": "O", "charge": 0, "multiplicity": 1},
               "bond": "O-H"},
    "fragments": {"fragments": [{"label": "ethyl radical", "smiles": "C[CH2]", "charge": 0, "multiplicity": 2},
                                {"label": "hydrogen atom", "smiles": "[H]", "charge": 0, "multiplicity": 2}]},
    "method": {"xc": "b3lyp", "basis": "6-31g*"},
}
_FALLBACK_RATIONALE = (
    "Breaking one O-H bond of water (H-O-H) leaves the other O-H bond intact, giving hydroxyl radical ([OH], "
    "neutral doublet: 9 total electrons (7 valence), odd, consistent with multiplicity 2) plus a hydrogen atom "
    "([H], neutral doublet: 1 electron, odd, consistent with multiplicity 2), not full atomization into O + 2H. "
    "Parent water is a closed-shell singlet (charge 0, multiplicity 1). Workflow: GFN2-xTB geometry optimization "
    "of parent and both fragments as a cheap prescreen, then a B3LYP/def2-TZVP DFT single point at those "
    "geometries (unrestricted for the radicals). BDE = E(OH) + E(H) - E(H2O), converted to kcal/mol. This is an "
    "electronic-energy BDE without ZPE or thermal corrections.")


def _reference_rationale() -> str:
    """The skill's own 'Good (water O-H)' rationale, the critic's yardstick."""
    import re

    m = re.search(r"Good \(water O-H\):\s*>\s*(.+?)\n\s*\n", _skill_text(), re.S)
    return m.group(1).strip() if m else _FALLBACK_RATIONALE


def module_example(module: str) -> str:
    """The format guidance shown to the Planner for one module."""
    if module == "rationale":
        return ("Reference rationale, for water: \"" + _reference_rationale() + "\"\nWrite yours with the same "
                "structure and the same items - the bond broken and what stays intact, each species' charge, total "
                "electron count and multiplicity with the parity argument, the parent's spin state, single-bond "
                "homolysis versus full atomization, each workflow stage and why, the functional and basis of THIS "
                "plan's accepted method module, the BDE formula and units, and that it is an electronic BDE "
                "without ZPE - but for this question's molecule, bond and fragments.")
    note = " Prefer the method the retrieved results above recommend, if any." if module == "method" else ""
    return ("Example of a well-formed '" + module + "' module, for a DIFFERENT molecule - copy its format "
            "(field names, JSON objects, SMILES), not its chemistry: " + json.dumps(_MODULE_EXAMPLES[module]) + note)


def validate_module(module: str, content: dict, accepted: dict, question: str) -> dict:
    """Deterministic review of one module, given the accepted ones.
    Returns {'verdict', 'problems', 'notes', 'signature'}."""
    problems: list[str] = []
    notes: list[str] = []
    if "error" in content:
        problems = [f"your reply could not be used ({content['error']}); call publish_{module} once"]
    elif module == "target":
        if content.get("analysis_type") != "bde" and bond_types_in(question):
            problems.append("the question asks for a bond dissociation energy: analysis_type must be 'bde'")
        problems += _species_problems(content.get("parent"), "parent")
        if not problems:
            parent = _fix_species(content["parent"])
            declared = bond_types_in(str(content.get("bond", ""))) or bond_types_in(question)
            asked = bond_types_in(question)
            if asked and declared and not set(declared) & set(asked):
                problems.append(f"the bond {content.get('bond')!r} is not the {'/'.join(asked)} bond the question asks about")
            cands = [c for c in enumerate_bde_candidates(parent["smiles"]) if not declared or c["bond"] in declared]
            if declared and not cands:
                problems.append(f"the parent {parent['smiles']!r} has no single {'/'.join(declared)} bond that can be broken")
            elif cands:
                notes.append(f"symmetry-unique {'/'.join(declared) or 'single'} bonds of {parent['smiles']}: "
                             + "; ".join(_split_label(c) for c in cands[:6]))
    elif module == "fragments":
        tgt = accepted.get("target") or {}
        spec = {"analysis_type": tgt.get("analysis_type", "bde"), "parent": _fix_species(tgt.get("parent")),
                "fragments": [_fix_species(f) for f in (_decode(content.get("fragments")) or [])], "rationale": ""}
        if not spec["fragments"] or not all(isinstance(f, dict) and "smiles" in f for f in spec["fragments"]):
            problems.append("fragments must be a JSON array of objects with smiles, charge and multiplicity")
        else:
            q = f"{question} ({tgt.get('bond', '')})"
            report = critique_plan(spec, q)
            problems = problems_as_text(report["fails"])
            if report.get("cleaved_bond"):
                notes.append("these fragments break " + report["cleaved_bond"]["bond"] + ": "
                             + _split_label(report["cleaved_bond"]))
            wanted = bond_types_in(q)
            opts = [_split_label(c) for c in report.get("candidates", []) if not wanted or c["bond"] in wanted]
            if opts:
                notes.append(f"the parent's symmetry-unique {'/'.join(wanted) or 'single'}-bond splits: " + "; ".join(opts[:6]))
    elif module == "method":
        try:
            from pyscf import dft, gto
            gto.basis.load(str(content.get("basis")), "C")
            dft.libxc.parse_xc(str(content.get("xc")))
        except Exception as exc:
            problems.append(f"PySCF cannot use xc={content.get('xc')!r} / basis={content.get('basis')!r} "
                            f"({type(exc).__name__})")
    elif module == "rationale":
        if len(str(content.get("rationale", "")).strip()) < 40:
            problems.append("the rationale is missing or too short")
    # What was wrong, for the supervisor's repeat detection.
    signature = [module] + sorted(p[:80] for p in problems) if problems else None
    return {"verdict": "FAIL" if problems else "PASS", "problems": problems, "notes": notes,
            "signature": signature}


class PlanBoard(Agent):
    """The shared board: module versions, reviews and accepted modules."""

    def __init__(self) -> None:
        super().__init__()
        self.reset_state("", "")

    def reset_state(self, question: str, context: str) -> None:
        self.question, self.context = question, context
        self.entries: list[dict] = []
        self.accepted: dict[str, int] = {}
        self.flags: dict[str, str] = {}

    @action
    async def reset(self, question: str, context: str) -> None:
        self.reset_state(question, context)

    @action
    async def post(self, module: str, content: dict, author: str) -> int:
        entry = {"id": len(self.entries), "module": module, "author": author, "content": content,
                 "reviews": [], "status": "proposed"}
        self.entries.append(entry)
        return entry["id"]

    @action
    async def review(self, entry_id: int, author: str, verdict: str, notes: list[str]) -> None:
        self.entries[entry_id]["reviews"].append({"author": author, "verdict": verdict, "notes": notes})

    @action
    async def accept(self, module: str, entry_id: int, flag: str | None = None) -> None:
        self.accepted[module] = entry_id
        self.entries[entry_id]["status"] = "accepted"
        if flag:
            self.flags[module] = flag

    @action
    async def reopen(self, modules: list[str]) -> None:
        for m in modules:
            if m in self.accepted:
                self.entries[self.accepted.pop(m)]["status"] = "reopened"

    @action
    async def snapshot(self) -> dict:
        return {"question": self.question, "context": self.context, "entries": self.entries,
                "accepted": {m: self.entries[i]["content"] for m, i in self.accepted.items()},
                "flags": dict(self.flags)}


def render_board(snap: dict, module: str | None) -> str:
    """The board as the Planner and the critic read it."""
    lines = [f"Question: {snap['question']}"]
    if snap.get("context"):
        lines += ["", snap["context"]]
    lines += ["", "Accepted modules:"]
    lines += [f"- {m}: {json.dumps(c)}" for m, c in snap["accepted"].items()] or ["- (none yet)"]
    if module:
        tries = [e for e in snap["entries"] if e["module"] == module and e["status"] != "accepted"]
        if tries:
            lines += ["", f"Earlier attempts at '{module}' and their reviews:"]
            for e in tries[-3:]:
                lines.append(f"- {json.dumps(e['content'])}")
                lines += [f"    {r['author']}: {r['verdict']} {' | '.join(r['notes'])}" for r in e["reviews"]]
    return "\n".join(lines)


class BoardPlannerAgent(Agent):
    """LLM planner that publishes one module at a time to the board."""

    def __init__(self, board: Handle[PlanBoard], backend: str, model: str | None = None) -> None:
        super().__init__()
        self.board, self.backend = board, backend
        self.model = P._resolve_model(backend, model)
        self.prompt = P._build_system_prompt() + _BOARD_INSTRUCTIONS

    @action
    async def get_prompt(self) -> str:
        return self.prompt

    @action
    async def update_prompt(self, new_prompt: str) -> None:
        self.prompt = new_prompt

    @action
    async def publish(self, module: str) -> int:
        snap = await self.board.snapshot()
        message = (render_board(snap, module) + "\n\n" + module_example(module)
                   + f"\n\nPublish the '{module}' module now.")
        raw = await self.agent_run_sync(self._call, module, message)
        content = {k: _decode(v) for k, v in raw.items()} if isinstance(raw, dict) else {"error": repr(raw)}
        return await self.board.post(module, content, "Planner")

    def _call(self, module: str, message: str) -> Any:
        attempt_fn, _ = P._BACKENDS[self.backend]
        with agent_llm_invocation_task("PlannerAgent", f"publish_{module}",
                                       used={"model": self.model, "backend": self.backend, "module": module}) as t:
            try:
                raw, _, _ = attempt_fn(self.model, self.prompt, [module_tool(module)],
                                       [{"role": "user", "content": message}])
            except Exception as exc:
                raw = {"error": f"{type(exc).__name__}: {exc}"[:600]}
            t.end(generated={"tool_input": raw})
        return raw


class BoardValidatorAgent(Agent):
    """Deterministic reviewer: posts PASS/FAIL for the latest version of a module."""

    def __init__(self, board: Handle[PlanBoard]) -> None:
        super().__init__()
        self.board = board

    @action
    async def evaluate(self, entry_id: int) -> dict:
        snap = await self.board.snapshot()
        entry = snap["entries"][entry_id]
        result = validate_module(entry["module"], entry["content"], snap["accepted"], snap["question"])
        await self.board.review(entry_id, "Validator", result["verdict"], result["problems"] + result["notes"])
        return result


class BoardCriticAgent(Agent):
    """LLM critic with the bde-plan-critic skill: reviews validated modules
    and the assembled plan."""

    def __init__(self, board: Handle[PlanBoard], backend: str, model: str | None = None) -> None:
        super().__init__()
        self.board, self.backend = board, backend
        self.model = P._resolve_model(backend, model)

    @action
    async def evaluate(self, entry_id: int) -> str:
        snap = await self.board.snapshot()
        entry = snap["entries"][entry_id]
        message = (render_board(snap, None) + f"\n\nModule under review: {entry['module']}\n"
                   + json.dumps(entry["content"]) + "\n\nValidator: "
                   + " | ".join(r["verdict"] + " " + " | ".join(r["notes"]) for r in entry["reviews"]))
        raw = await self.agent_run_sync(self._call, entry["module"], message, _MODULE_REVIEW_TOOL)
        verdict = str(raw.get("verdict", "")).upper()
        verdict = verdict if verdict in ("ACCEPT", "REVISE") else "ERROR"
        await self.board.review(entry_id, "PlanCritic", verdict, _as_list(raw.get("issues")))
        return verdict

    @action
    async def final_review(self, plan: dict) -> dict:
        snap = await self.board.snapshot()
        message = render_board(snap, None) + "\n\nAssembled plan:\n" + json.dumps(plan)
        raw = await self.agent_run_sync(self._call, None, message, _FINAL_REVIEW_TOOL)
        verdict = str(raw.get("verdict", "")).upper()
        reopen = [m for m in _as_list(raw.get("reopen")) if m in MODULES]
        return {"verdict": verdict if verdict in ("ACCEPT", "REVISE") else "ERROR",
                "reopen": reopen, "issues": _as_list(raw.get("issues"))}

    def _call(self, module: str | None, message: str, tool: dict) -> dict:
        attempt_fn, _ = P._BACKENDS[self.backend]
        with agent_llm_invocation_task("PlanCriticAgent", f"review_{module or 'plan'}",
                                       used={"model": self.model, "backend": self.backend, "module": module}) as t:
            try:
                raw, _, _ = attempt_fn(self.model, _skill_text() + _critic_instructions(module), [tool],
                                       [{"role": "user", "content": message}])
            except Exception as exc:
                raw = {"verdict": "ERROR", "error": f"{type(exc).__name__}: {exc}"[:300]}
            t.end(generated={"tool_input": raw})
        return raw if isinstance(raw, dict) else {"verdict": "ERROR"}


class PlanBoardModerator(Agent):
    """Runs the board: module rounds, acceptance, final review, supervision."""

    def __init__(self, board: Handle[PlanBoard], planner: Handle[BoardPlannerAgent],
                 validator: Handle[BoardValidatorAgent], critic: Handle[BoardCriticAgent], *,
                 backend: str, model: str | None = None, max_attempts: int | None = None,
                 deadline_s: float | None = 3 * 60, critic_revise_limit: int = 3,
                 reopen_budget: int = 1) -> None:
        super().__init__()
        self.board, self.planner, self.validator, self.critic = board, planner, validator, critic
        self.backend = backend
        self.model = P._resolve_model(backend, model)
        # max_attempts None/0: keep going until the module is accepted or the
        # planning time budget (deadline_s) runs out.
        self.max_attempts, self.reopen_budget = max_attempts or None, reopen_budget
        self.deadline_s, self.critic_revise_limit = deadline_s, critic_revise_limit
        self.started = time.monotonic()
        self.timed_out = False
        self.rounds: list[dict] = []
        self.supervisor_log: list[dict] = []
        self.running = asyncio.Event()
        self.running.set()
        self.terminate = asyncio.Event()
        self.round_done = asyncio.Event()
        self.round_checked = asyncio.Event()

    def _out_of_time(self) -> bool:
        return self.deadline_s is not None and time.monotonic() - self.started > self.deadline_s

    async def _build_module(self, module: str) -> bool:
        last_validated, revisions, attempt = None, 0, 0
        self.terminate.clear()
        while self.max_attempts is None or attempt < self.max_attempts:
            if self._out_of_time():
                self.timed_out = True
                break
            attempt += 1
            await self.running.wait()
            entry_id = await self.planner.publish(module)
            val = await self.validator.evaluate(entry_id)
            crit = None
            if val["verdict"] == "PASS":
                last_validated = entry_id
                crit = await self.critic.evaluate(entry_id)
            self.rounds.append({"round": len(self.rounds) + 1, "module": module, "attempt": attempt,
                                "validator": val["verdict"], "critic": crit, "signature": val["signature"]})
            if val["verdict"] == "PASS" and crit in ("ACCEPT", "ERROR"):
                await self.board.accept(module, entry_id)
                return True
            if val["verdict"] == "PASS":
                # The critic's veto is bounded: after critic_revise_limit
                # rejections of validated versions, keep the last one.
                revisions += 1
                if revisions >= self.critic_revise_limit:
                    break
            self.round_checked.clear()
            self.round_done.set()
            await self.round_checked.wait()
            if self.terminate.is_set():
                break
        if last_validated is not None:
            await self.board.accept(module, last_validated, flag="critic_unresolved")
            return True
        return False

    @action
    async def plan(self, question: str, context: str = "") -> dict:
        """Build the plan for one question. Returns {'spec' (None on failure),
        'outcome', 'rounds', 'supervisor', 'final_reviews', 'board'}."""
        await self.board.reset(question, context)
        self.started, self.timed_out = time.monotonic(), False
        todo, final_reviews, outcome, spec = list(MODULES), [], "failed", None
        while True:
            for module in todo:
                if not await self._build_module(module):
                    snap = await self.board.snapshot()
                    why = " (planning time budget used up)" if self.timed_out else ""
                    return {"spec": None, "outcome": f"no valid {module}{why}", "rounds": self.rounds,
                            "supervisor": self.supervisor_log, "final_reviews": final_reviews, "board": snap}
            snap = await self.board.snapshot()
            acc = snap["accepted"]
            raw = {**acc["target"], "parent": _fix_species(acc["target"].get("parent")),
                   "fragments": [_fix_species(f) for f in _decode(acc["fragments"]["fragments"])],
                   **acc["method"], "rationale": acc["rationale"].get("rationale")}
            raw.pop("bond", None)
            spec = P._spec_from_raw(raw, question, self.model)
            full = critique_plan(spec, f"{question} ({acc['target'].get('bond', '')})")
            if full["fails"]:
                return {"spec": None, "outcome": "assembled plan failed validation: "
                        + "; ".join(problems_as_text(full["fails"])), "rounds": self.rounds,
                        "supervisor": self.supervisor_log, "final_reviews": final_reviews, "board": snap}
            review = await self.critic.final_review(raw)
            final_reviews.append(review)
            self.rounds.append({"round": len(self.rounds) + 1, "module": "final", "attempt": len(final_reviews),
                                "validator": "PASS", "critic": review["verdict"], "signature": None})
            if review["verdict"] == "REVISE" and review["reopen"] and self.reopen_budget > 0:
                self.reopen_budget -= 1
                todo = [m for m in MODULES if m in review["reopen"]]
                await self.board.reopen(todo)
                continue
            outcome = "passed" if review["verdict"] != "REVISE" else "passed_final_review_unresolved"
            break
        spec["bond_requested"] = acc["target"].get("bond")
        spec["critic_unresolved"] = sorted(snap["flags"])
        return {"spec": spec, "outcome": outcome, "rounds": self.rounds, "supervisor": self.supervisor_log,
                "final_reviews": final_reviews, "board": await self.board.snapshot()}

    @loop
    async def supervise(self, shutdown: asyncio.Event) -> None:
        """After every rejected round: a repeated Validator failure on the
        same module means the Planner is going in circles - pause and
        correct its prompt. A repeat after that correction escalates to the
        Validator's exact fixes (unbounded mode) or terminates the module
        (with a max_attempts cap)."""
        while not shutdown.is_set():
            try:
                await asyncio.wait_for(self.round_done.wait(), timeout=1.0)
            except asyncio.TimeoutError:
                continue
            self.round_done.clear()
            try:
                await self._check_progress()
            finally:
                self.round_checked.set()

    async def _check_progress(self) -> None:
        sigs = [r["signature"] for r in self.rounds if r["module"] == self.rounds[-1]["module"]]
        if len(sigs) < 2 or not sigs[-1] or sigs[-1] != sigs[-2]:
            return
        edits = [e for e in self.supervisor_log if e["signature"] == sigs[-1] and e["action"].startswith("PAUSE")]
        if edits and self.max_attempts is not None:
            self.supervisor_log.append({"round": self.rounds[-1]["round"], "action": "TERMINATE",
                                        "signature": sigs[-1], "reason": "same mistake after a correction"})
            self.terminate.set()
            return
        if len(edits) >= 2:
            return  # already escalated for this mistake; keep going until the time budget
        self.running.clear()  # PAUSE
        try:
            snap = await self.board.snapshot()
            if edits:  # second time: no more LLM wording, the Validator's own fixes
                last = [e for e in snap["entries"] if e["module"] == self.rounds[-1]["module"]][-1]
                instruction = ("You repeated this mistake after a correction. Fix exactly what the Validator "
                               "says: " + " | ".join(n for r in last["reviews"] for n in r["notes"])[:1500]
                               + " Use the publish tool's exact field names: smiles, charge, multiplicity.")
                diagnosis, action = "repeated after the first correction", "PAUSE+ESCALATE"
            else:
                instruction, diagnosis = await self.agent_run_sync(self._write_instruction, snap)
                action = "PAUSE+EDIT"
            prompt = await self.planner.get_prompt()
            await self.planner.update_prompt(prompt + "\n\nSupervisor instruction (added after a repeated "
                                             "mistake): " + instruction)
            self.supervisor_log.append({"round": self.rounds[-1]["round"], "action": action,
                                        "signature": sigs[-1], "diagnosis": diagnosis, "instruction": instruction})
        finally:
            self.running.set()

    def _write_instruction(self, snap: dict) -> tuple[str, str]:
        module = self.rounds[-1]["module"]
        board = render_board(snap, module)
        attempt_fn, _ = P._BACKENDS[self.backend]
        with agent_llm_invocation_task("PlanBoardModerator", "supervisor_edit",
                                       used={"model": self.model, "backend": self.backend, "module": module}) as t:
            try:
                raw, _, _ = attempt_fn(self.model, _EDIT_PROMPT + _EDIT_FIELD_NAMES, [_INSTRUCTION_TOOL],
                                       [{"role": "user", "content": "Validator reviews are on the board.\n\n"
                                                                    + board[-12000:]}])
                instruction = str(raw.get("instruction") or "").strip()
                diagnosis = str(raw.get("diagnosis") or "").strip()
            except Exception as exc:
                instruction, diagnosis = "", f"supervisor LLM failed: {type(exc).__name__}: {exc}"[:300]
            if not instruction:
                last = [e for e in snap["entries"] if e["module"] == module][-1]
                instruction = ("Fix exactly this, which you have now got wrong twice: "
                               + " | ".join(n for r in last["reviews"] for n in r["notes"])[:1500])
            t.end(generated={"instruction": instruction, "diagnosis": diagnosis})
        return instruction, diagnosis
