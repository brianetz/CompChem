"""
Smoke test for the plan board (agents/plan_board.py): the real Academy
agents in a local Manager, with a scripted stand-in for the LLM.

Run:
    python src/smoke_test_plan_board.py
"""
from __future__ import annotations

import asyncio
import time
from concurrent.futures import ThreadPoolExecutor

from academy.exchange import LocalExchangeFactory
from academy.manager import Manager

import llm.planner_llm as P
from agents.plan_board import (
    BoardCriticAgent, BoardPlannerAgent, BoardValidatorAgent, PlanBoard, PlanBoardModerator,
)
from storage.flowcept_store import provenance_session

TQ = "What's the benzylic C-H bond dissociation energy of toluene in kcal/mol?"
TARGET = {"analysis_type": "bde", "parent": {"smiles": "Cc1ccccc1", "charge": 0, "multiplicity": 1},
          "bond": "benzylic C-H"}
H = {"smiles": "[H]", "charge": 0, "multiplicity": 2}
METHOD = {"xc": "b3lyp", "basis": "6-31g*"}
RATIONALE = {"rationale": "Breaking the benzylic C-H bond of toluene gives the benzyl radical and an H atom; "
                          "electronic BDE = E(benzyl) + E(H) - E(toluene) in kcal/mol."}


def frags(s):
    return {"fragments": [{"smiles": s, "charge": 0, "multiplicity": 2, "label": "benzyl"}, H]}


S: dict = {}
PROMPTS: list[str] = []


def stub(model, system_prompt, tools, messages):
    name, msg = tools[0]["name"], messages[0]["content"]
    if name.startswith("publish_"):
        PROMPTS.append(system_prompt)
        time.sleep(S.get("delay", 0))
        return S[name[len("publish_"):]].pop(0), None, None
    if name == "emit_module_review":
        assert "BDE plan critic" in system_prompt
        if "Module under review: fragments" in msg and "Cc1ccc[c]c1" in msg.split("Module under review")[1]:
            return {"verdict": "REVISE", "issues": ["this is an aromatic C-H, not benzylic"]}, None, None
        if "Module under review: rationale" in msg and S.get("rationale_revise"):
            return {"verdict": "REVISE", "issues": ["missing electron counts"]}, None, None
        return {"verdict": "ACCEPT", "issues": []}, None, None
    if name == "emit_final_review":
        return S["final"].pop(0), None, None
    if name == "emit_instruction":
        return {"diagnosis": "brackets a ring carbon", "instruction": "Write benzyl as [CH2]c1ccccc1."}, None, None
    raise AssertionError(name)


P._BACKENDS["stub"] = (stub, lambda *a: a[0])
OK = {"verdict": "ACCEPT", "reopen": [], "issues": []}


async def run(script, **moderator):
    S.clear(); S.update({k: list(v) if isinstance(v, list) else v for k, v in script.items()})
    S.setdefault("final", [OK])
    PROMPTS.clear()
    async with await Manager.from_exchange_factory(factory=LocalExchangeFactory(), executors=ThreadPoolExecutor()) as m:
        board = await m.launch(PlanBoard())
        kw = {"board": board, "backend": "stub", "model": "m"}
        planner = await m.launch(BoardPlannerAgent, kwargs=kw)
        validator = await m.launch(BoardValidatorAgent, kwargs={"board": board})
        critic = await m.launch(BoardCriticAgent, kwargs=kw)
        mod = await m.launch(PlanBoardModerator, kwargs={"board": board, "planner": planner, "validator": validator,
                                                         "critic": critic, "backend": "stub", "model": "m",
                                                         **moderator})
        return await mod.plan(TQ, "")


def show(tag, r):
    print(tag, r["outcome"], [(x["module"], x["attempt"], x["validator"], x["critic"]) for x in r["rounds"]],
          [s["action"] for s in r["supervisor"]])


async def main():
    base = {"target": [TARGET], "method": [METHOD], "rationale": [RATIONALE]}

    # A: fragments wrong twice -> supervisor edits the Planner prompt -> fixed.
    r = await run({**base, "fragments": [frags("C[c]1ccccc1"), frags("C[c]1ccccc1"), frags("[CH2]c1ccccc1")]})
    show("A", r)
    assert r["outcome"] == "passed" and r["spec"]["fragments"][0]["smiles"] == "[CH2]c1ccccc1"
    assert [s["action"] for s in r["supervisor"]] == ["PAUSE+EDIT"]
    assert any("Write benzyl as [CH2]c1ccccc1." in p for p in PROMPTS)
    assert len(r["rounds"]) == 7  # target, 3x fragments, method, rationale, final review
    notes = [n for e in r["board"]["entries"] if e["module"] == "fragments" for rv in e["reviews"] for n in rv["notes"]]
    assert any("[CH2]c1ccccc1 + [H]" in n for n in notes), notes  # Validator shows the benzylic split

    # B: Validator passes a ring C-H; the critic sends it back; benzylic accepted.
    r = await run({**base, "fragments": [frags("Cc1ccc[c]c1"), frags("[CH2]c1ccccc1")]})
    show("B", r)
    fr = [x for x in r["rounds"] if x["module"] == "fragments"]
    assert [(x["validator"], x["critic"]) for x in fr] == [("PASS", "REVISE"), ("PASS", "ACCEPT")]
    assert r["spec"]["fragments"][0]["smiles"] == "[CH2]c1ccccc1"

    # C: the final review reopens the method module once.
    r = await run({**base, "fragments": [frags("[CH2]c1ccccc1")], "method": [METHOD, METHOD],
                   "final": [{"verdict": "REVISE", "reopen": ["method"], "issues": ["x"]}, OK]})
    show("C", r)
    assert r["outcome"] == "passed" and [x["module"] for x in r["rounds"]].count("method") == 2

    # D: the critic never accepts the rationale -> kept, flagged.
    r = await run({**base, "fragments": [frags("[CH2]c1ccccc1")], "rationale": [RATIONALE] * 3,
                   "rationale_revise": True})
    show("D", r)
    assert r["spec"] is not None and r["spec"]["critic_unresolved"] == ["rationale"]

    # E: with a cap of 3 attempts, no valid target -> no plan.
    bad = {"analysis_type": "bde", "parent": {"smiles": "[CH3O]", "charge": 0, "multiplicity": 1}, "bond": "C-H"}
    r = await run({"target": [bad] * 3}, max_attempts=3)
    show("E", r)
    assert r["spec"] is None and r["outcome"] == "no valid target"

    # F: no cap (default): wrong 4 times, then right -> edit, then escalate, then a plan.
    r = await run({**base, "target": [bad] * 4 + [TARGET], "fragments": [frags("[CH2]c1ccccc1")]})
    show("F", r)
    assert r["outcome"] == "passed" and [s["action"] for s in r["supervisor"]] == ["PAUSE+EDIT", "PAUSE+ESCALATE"]
    assert any("Use the publish tool's exact field names" in p for p in PROMPTS)

    # G: never valid -> stops when the time budget runs out, not after N attempts.
    r = await run({"target": [bad] * 1000, "delay": 0.05}, deadline_s=0.5)
    show("G", (r | {"rounds": r["rounds"][-2:]}))
    assert r["spec"] is None and r["outcome"] == "no valid target (planning time budget used up)"
    assert len(r["rounds"]) > 3

    print("\nOK: plan board behaves as expected.")


if __name__ == "__main__":
    with provenance_session():
        asyncio.run(main())
