"""
Plan discussion: BDE planning as an Academy group chat.

Follows the pattern of Academy's examples/08-discussion (GroupChatAgent +
RoundRobinGroupChatManager): every participant is an Academy agent that
holds handles to its peers and broadcasts what it says to them, so they
share one conversation; a moderator agent runs the rounds, decides when to
stop, and supervises the conversation from a background control loop.

Participants, in turn order:
  Planner     LLM. Proposes a campaign spec through the emit_campaign_spec
              tool, seeing everything said so far.
  Validator   Deterministic (analysis/bde_plan_validator.py). Checks the
              latest proposal and says PASS or FAIL with exact fixes; on a
              FAIL it also lists the symmetry-unique splits of the bond the
              question asks about.
  PlanCritic  LLM, with the bde-plan-critic skill as its instructions
              (only when COMPCHEM_PLAN_CRITIC_SKILL is set). Speaks only
              after a Validator PASS: code results are never re-argued.

Moderator (PlanDiscussionManager):
  - stops as soon as the Validator says PASS and the critic (if any) does
    not say FAIL; otherwise runs another round, up to max_rounds;
  - supervise (@loop): after every round, if the Validator rejected the
    same thing two rounds running - the planner repeating its mistake -
    the moderator PAUSES the chat, asks the LLM for one corrective
    instruction for the Planner and appends it to the Planner's system
    prompt (update_prompt), then resumes. If the same mistake comes back
    after that edit, it TERMINATES the discussion.

Differences from the example, on purpose: the stopping decision is
deterministic (the Validator's verdict) instead of an LLM call, and the
moderator appends an instruction to the Planner's prompt instead of
letting an LLM rewrite the whole prompt.
"""
from __future__ import annotations

import asyncio
import json
import logging
from typing import Any

from academy.agent import Agent, action, loop
from academy.handle import Handle

from analysis.bde_plan_validator import critique_plan, problems_as_text
from llm import planner_llm as P
from llm.plan_checks import check_spec, normalize_raw
from llm.plan_critic_llm import review_plan, skill_path
from storage.agent_provenance import agent_llm_invocation_task
from storage.rag import bond_types_in

logger = logging.getLogger(__name__)

_INSTRUCTION_TOOL = {
    "name": "emit_instruction",
    "description": "Return one corrective instruction for the Planner.",
    "input_schema": {
        "type": "object",
        "properties": {
            "diagnosis": {"type": "string", "description": "How the Planner keeps failing, in one sentence."},
            "instruction": {"type": "string", "description": "One instruction to add to the Planner's system prompt."},
        },
        "required": ["diagnosis", "instruction"],
    },
}

_EDIT_PROMPT = """You are supervising a discussion between a Planner (an LLM
that writes computational chemistry campaign specs), a deterministic
Validator and a plan critic. The Planner has made the SAME mistake in two
consecutive rounds even though the Validator explained it. Read the
discussion, identify exactly what the Planner keeps getting wrong, and
write ONE short, concrete instruction that will be added to the Planner's
system prompt so that its next proposal fixes it. Use the Validator's own
fixes and the bond splits it listed. Call emit_instruction once."""


def _merge_turns(messages: list[dict[str, str]]) -> list[dict[str, str]]:
    """Merge consecutive same-role turns (several peers speaking in a row),
    which some chat APIs require."""
    out: list[dict[str, str]] = []
    for m in messages:
        if out and out[-1]["role"] == m["role"]:
            out[-1] = {"role": m["role"], "content": out[-1]["content"] + "\n\n" + m["content"]}
        else:
            out.append(dict(m))
    return out


def validate_proposal(raw: Any, question: str, model: str) -> dict:
    """The Validator's whole judgment of one Planner proposal, as data:
    {'verdict', 'spec', 'problems' (instructions), 'report' (bde only),
    'signature' (what was wrong, for repeat detection)}."""
    if not isinstance(raw, dict) or "error" in raw:
        err = raw.get("error") if isinstance(raw, dict) else repr(raw)
        return {"verdict": "FAIL", "spec": None, "report": None,
                "problems": [f"your reply could not be used ({err}); call emit_campaign_spec with JSON objects"],
                "signature": ("unusable reply",)}
    fixed, repairs = normalize_raw(raw)
    try:
        spec = P._spec_from_raw(fixed, question, model)
    except ValueError as exc:
        return {"verdict": "FAIL", "spec": None, "report": None, "problems": [str(exc)],
                "signature": ("shape", str(exc)[:80])}
    if repairs:
        spec["planner_repairs"] = repairs
    if spec.get("analysis_type") != "bde":
        problems = check_spec(spec, question)
        return {"verdict": "FAIL" if problems else "PASS", "spec": spec, "report": None,
                "problems": problems, "signature": tuple(problems) or None}

    report = critique_plan(spec, question)
    problems = problems_as_text(report["fails"])
    if problems:
        # Show the right answer's shape: the symmetry-unique splits of the
        # bond the question names (the wrong-bond and non-homolysis fixes
        # already list them; parity and SMILES failures did not).
        wanted = bond_types_in(question)
        options = [" + ".join(c["fragments"]) for c in report.get("candidates", [])
                   if not wanted or c["bond"] in wanted]
        if options and not any("homolyses" in p for p in problems):
            problems.append(f"The parent's symmetry-unique {'/'.join(wanted) or 'single'}-bond "
                            f"splits are: {'; '.join(options[:6])}")
    frags = tuple(sorted(str(f.get("smiles")) for f in spec.get("fragments") or []))
    signature = (tuple(sorted(f["check"] for f in report["fails"])), frags) if report["fails"] else None
    return {"verdict": report["verdict"] if report["verdict"] == "FAIL" else "PASS",
            "spec": spec, "report": report, "problems": problems, "signature": signature}


class PlanChatAgent(Agent):
    """One participant in the plan discussion; the example's GroupChatAgent
    primitives (peers, receive, respond, prompt get/update, clear)."""

    def __init__(self, role: str, prompt: str) -> None:
        super().__init__()  # agent_run_sync needs Agent's own state
        self.role = role
        self.prompt = prompt
        self.peers: list[Handle[Any]] = []
        self.messages: list[dict[str, str]] = []

    @action
    async def add_peer(self, peer: Handle[Any]) -> None:
        """Peers receive every message this agent sends."""
        self.peers.append(peer)

    @action
    async def get_role(self) -> str:
        return self.role

    @action
    async def get_prompt(self) -> str:
        return self.prompt

    @action
    async def update_prompt(self, new_prompt: str) -> None:
        self.prompt = new_prompt

    @action
    async def receive(self, role: str, content: str) -> None:
        self.messages.append({"role": "user", "content": f"{role}: {content}"})
        self._on_receive(role, content)

    def _on_receive(self, role: str, content: str) -> None:
        pass

    @action
    async def clear(self) -> None:
        self.messages = []

    async def _say(self, content: str) -> str:
        self.messages.append({"role": "assistant", "content": content})
        for peer in self.peers:
            await peer.receive(self.role, content)
        return content


class PlannerChatAgent(PlanChatAgent):
    """LLM planner: proposes a spec (as the emit_campaign_spec tool input,
    JSON) given the whole discussion."""

    def __init__(self, backend: str, model: str | None = None) -> None:
        super().__init__("Planner", P._build_system_prompt())
        self.backend = backend
        self.model = P._resolve_model(backend, model)
        self.turn = 0

    @action
    async def respond(self) -> str:
        self.turn += 1
        raw = await self.agent_run_sync(self._propose)
        return await self._say(json.dumps(raw))

    def _propose(self) -> Any:
        attempt_fn, _ = P._BACKENDS[self.backend]
        with agent_llm_invocation_task(
            "PlannerAgent", "plan_discussion_turn",
            used={"model": self.model, "backend": self.backend, "attempt": self.turn},
        ) as t:
            try:
                raw, _, _ = attempt_fn(self.model, self.prompt, [P._build_tool_schema()],
                                       _merge_turns(self.messages))
            except Exception as exc:  # unusable reply, HTTP 500, ...
                raw = {"error": f"{type(exc).__name__}: {exc}"[:600]}
            t.end(generated={"tool_input": raw})
        return raw


class ValidatorChatAgent(PlanChatAgent):
    """Deterministic validator: judges the latest Planner proposal."""

    def __init__(self, model_label: str = "") -> None:
        super().__init__("Validator", "deterministic plan validator")
        self.question = ""
        self.proposal: Any = None
        self.model_label = model_label
        self.last: dict = {}

    def _on_receive(self, role: str, content: str) -> None:
        if role == "User" and not self.question:
            self.question = content.split("\n\n", 1)[0]
        elif role == "Planner":
            try:
                self.proposal = json.loads(content)
            except ValueError:
                self.proposal = {"error": "reply was not JSON"}

    @action
    async def respond(self) -> str:
        self.last = validate_proposal(self.proposal, self.question, self.model_label)
        if self.last["verdict"] == "PASS":
            warns = [f"{w['evidence']} ({w['fix']})" for w in (self.last["report"] or {}).get("warns", [])]
            msg = "PASS" + (". Warnings: " + "; ".join(warns) if warns else "")
        else:
            msg = "FAIL: " + " | ".join(self.last["problems"])
        return await self._say(msg)

    @action
    async def latest(self) -> dict:
        """The last judgment, as data (spec, verdict, signature, report)."""
        out = dict(self.last)
        out["signature"] = list(out["signature"]) if out.get("signature") else None
        return out


class CriticChatAgent(PlanChatAgent):
    """LLM plan critic with the bde-plan-critic skill as its instructions.
    Reviews only plans the Validator passed."""

    def __init__(self, backend: str, model: str | None = None) -> None:
        super().__init__("PlanCritic", f"skill: {skill_path()}")
        self.backend = backend
        self.model = P._resolve_model(backend, model)
        self.question = ""
        self.proposal: Any = None

    def _on_receive(self, role: str, content: str) -> None:
        if role == "User" and not self.question:
            self.question = content.split("\n\n", 1)[0]
        elif role == "Planner":
            try:
                self.proposal = json.loads(content)
            except ValueError:
                self.proposal = None

    @action
    async def respond(self) -> str:
        judged = validate_proposal(self.proposal, self.question, self.model)
        if judged["verdict"] != "PASS" or judged["report"] is None:
            return await self._say("SKIP: the Validator has not passed this plan")
        review = await self.agent_run_sync(
            review_plan, judged["spec"], self.question, judged["report"],
            backend=self.backend, model=self.model,
        )
        parts = [review["verdict"]]
        if review.get("failed_checks"):
            parts.append("Failed: " + "; ".join(review["failed_checks"]))
        if review.get("warnings"):
            parts.append("Warnings: " + "; ".join(review["warnings"]))
        return await self._say(". ".join(parts))


class PlanDiscussionManager(Agent):
    """Moderator: wires the participants, runs round-robin rounds, stops on
    a Validator PASS (and no critic FAIL), and supervises the discussion
    from a control loop."""

    def __init__(
        self,
        participants: list[Handle[Any]],
        *,
        backend: str,
        model: str | None = None,
        max_rounds: int = 3,
    ) -> None:
        super().__init__()
        self.participants = participants
        self.backend = backend
        self.model = P._resolve_model(backend, model)
        self.max_rounds = max_rounds
        self.history: list[dict[str, str]] = []
        self.rounds: list[dict] = []
        self.supervisor_log: list[dict] = []
        self.running = asyncio.Event()
        self.running.set()
        self.terminate = asyncio.Event()
        self.round_done = asyncio.Event()
        self.round_checked = asyncio.Event()

    async def agent_on_startup(self) -> None:
        for participant in self.participants:
            for peer in self.participants:
                if peer != participant:
                    await participant.add_peer(peer)
            await participant.add_peer(Handle(self.agent_id))

    @action
    async def receive(self, role: str, content: str) -> None:
        logger.info("%s says: %s", role, content[:300])
        self.history.append({"role": role, "content": content})

    @action
    async def plan(self, question: str, context: str = "") -> dict:
        """Run the discussion for one question. Returns {'spec' (None if no
        plan passed), 'outcome', 'rounds', 'supervisor', 'history'}."""
        first = f"{question}\n\n{context}" if context else question
        self.history.append({"role": "User", "content": first})
        for participant in self.participants:
            await participant.receive("User", first)

        outcome, spec, last_passed = "max_rounds", None, None
        for rnd in range(1, self.max_rounds + 1):
            judged: dict = {}
            critic_verdict = None
            for participant in self.participants:
                await self.running.wait()
                if self.terminate.is_set():
                    break
                role = await participant.get_role()
                said = await participant.respond()
                if role == "Validator":
                    judged = await participant.latest()
                elif role == "PlanCritic" and not said.startswith("SKIP"):
                    critic_verdict = said.split(".", 1)[0].split(":", 1)[0].strip().upper()
            if self.terminate.is_set():
                outcome = "terminated"
                break
            self.rounds.append({"round": rnd, "validator": judged.get("verdict"),
                                "signature": judged.get("signature"), "critic": critic_verdict,
                                "fragments": [f.get("smiles") for f in (judged.get("spec") or {}).get("fragments") or []]})
            if judged.get("verdict") == "PASS":
                last_passed = judged["spec"]
                if critic_verdict != "FAIL":
                    spec, outcome = last_passed, "passed"
                    break
            # Hand the round to the supervisor and wait for its check, so the
            # next round starts after any prompt edit (deterministic order).
            self.round_checked.clear()
            self.round_done.set()
            await self.round_checked.wait()
            if self.terminate.is_set():
                outcome = "terminated"
                break
        if outcome != "passed" and last_passed is not None:
            # Never lose a plan the Validator passed: the critic only judges
            # wording and bond choice, and a later round may have broken it.
            spec, outcome = last_passed, outcome + "_kept_last_validated_plan"
        return {"spec": spec, "outcome": outcome, "rounds": self.rounds,
                "supervisor": self.supervisor_log, "history": self.history}

    @loop
    async def supervise(self, shutdown: asyncio.Event) -> None:
        """Background check after every round: a repeated Validator failure
        means the Planner is going in circles - pause and correct it, or
        terminate if a correction already failed."""
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
        sigs = [r["signature"] for r in self.rounds]
        if len(sigs) < 2 or not sigs[-1] or sigs[-1] != sigs[-2]:
            return
        if any(e["signature"] == sigs[-1] for e in self.supervisor_log):
            self.supervisor_log.append({"round": len(sigs), "action": "TERMINATE", "signature": sigs[-1],
                                        "reason": "same mistake after a prompt correction"})
            self.terminate.set()
            return
        self.running.clear()  # PAUSE
        try:
            planner = None
            for participant in self.participants:
                if await participant.get_role() == "Planner":
                    planner = participant
            if planner is None:
                return
            instruction, diagnosis = await self.agent_run_sync(self._write_instruction)
            prompt = await planner.get_prompt()
            await planner.update_prompt(prompt + "\n\nSupervisor instruction (added after a repeated "
                                        "mistake): " + instruction)
            self.supervisor_log.append({"round": len(sigs), "action": "PAUSE+EDIT", "signature": sigs[-1],
                                        "diagnosis": diagnosis, "instruction": instruction})
        finally:
            self.running.set()

    def _write_instruction(self) -> tuple[str, str]:
        transcript = "\n".join(f"{m['role']}: {m['content']}" for m in self.history)[-12000:]
        attempt_fn, _ = P._BACKENDS[self.backend]
        with agent_llm_invocation_task(
            "PlanDiscussionManager", "supervisor_edit",
            used={"model": self.model, "backend": self.backend},
        ) as t:
            try:
                raw, _, _ = attempt_fn(self.model, _EDIT_PROMPT, [_INSTRUCTION_TOOL],
                                       [{"role": "user", "content": transcript}])
                instruction = str(raw.get("instruction") or "").strip()
                diagnosis = str(raw.get("diagnosis") or "").strip()
            except Exception as exc:
                instruction, diagnosis = "", f"supervisor LLM failed: {type(exc).__name__}: {exc}"[:300]
            if not instruction:
                # Fall back to the Validator's own last fixes.
                last_validator = next((m["content"] for m in reversed(self.history)
                                       if m["role"] == "Validator"), "")
                instruction = "Fix exactly this, which you have now got wrong twice: " + last_validator[:1500]
            t.end(generated={"instruction": instruction, "diagnosis": diagnosis})
        return instruction, diagnosis
