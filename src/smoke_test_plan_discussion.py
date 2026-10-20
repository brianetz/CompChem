"""
Smoke test for the plan discussion (agents/plan_discussion.py): the real
Academy agents and moderator in a local Manager, with a scripted stand-in
for the LLM, so no model or GPU is needed.

Run:
    python src/smoke_test_plan_discussion.py
"""
from __future__ import annotations

import asyncio
import os
from concurrent.futures import ThreadPoolExecutor

from academy.exchange import LocalExchangeFactory
from academy.manager import Manager

import llm.planner_llm as P
from agents.plan_discussion import (
    CriticChatAgent, PlanDiscussionManager, PlannerChatAgent, ValidatorChatAgent,
)
from storage.flowcept_store import provenance_session

TQ = "What's the benzylic C-H bond dissociation energy of toluene in kcal/mol?"
H = {"smiles": "[H]", "charge": 0, "multiplicity": 2, "label": "H atom"}


def plan(frag_smiles):
    return {"analysis_type": "bde", "parent": {"smiles": "Cc1ccccc1", "charge": 0, "multiplicity": 1},
            "fragments": [{"smiles": frag_smiles, "charge": 0, "multiplicity": 2, "label": "benzyl"}, H],
            "basis": "6-31g*", "xc": "b3lyp", "rationale": "benzylic C-H homolysis"}


SCRIPT: dict = {}
SEEN_PLANNER_PROMPTS: list[str] = []


def stub_attempt(model, system_prompt, tools, messages):
    name = tools[0]["name"]
    if name == "emit_campaign_spec":
        SEEN_PLANNER_PROMPTS.append(system_prompt)
        return SCRIPT["planner"].pop(0), None, None
    if name == "emit_plan_review":
        return SCRIPT["critic"].pop(0), None, None
    if name == "emit_instruction":
        assert "Validator" in messages[0]["content"]
        return {"diagnosis": "brackets the ring carbon", "instruction": "Write benzyl as [CH2]c1ccccc1."}, None, None
    raise AssertionError(name)


P._BACKENDS["stub"] = (stub_attempt, lambda *a: a[0])


async def discuss(planner_replies, critic_replies=(), critic=False, rounds=3):
    SCRIPT["planner"] = list(planner_replies)
    SCRIPT["critic"] = list(critic_replies)
    SEEN_PLANNER_PROMPTS.clear()
    if critic:
        os.environ["COMPCHEM_PLAN_CRITIC_SKILL"] = os.path.join(os.path.dirname(__file__), "bde-plan-critic.md")
    else:
        os.environ.pop("COMPCHEM_PLAN_CRITIC_SKILL", None)
    async with await Manager.from_exchange_factory(
        factory=LocalExchangeFactory(), executors=ThreadPoolExecutor(),
    ) as manager:
        parts = [await manager.launch(PlannerChatAgent("stub", "m")),
                 await manager.launch(ValidatorChatAgent("m"))]
        if critic:
            parts.append(await manager.launch(CriticChatAgent("stub", "m")))
        mod = await manager.launch(PlanDiscussionManager, kwargs={
            "participants": parts, "backend": "stub", "model": "m", "max_rounds": rounds})
        result = await mod.plan(TQ, "")
    return result


def summary(result):
    return result["outcome"], [(r["round"], r["validator"], r["critic"]) for r in result["rounds"]], \
        [s["action"] for s in result["supervisor"]]


async def main():
    # A: same mistake twice -> supervisor edits the Planner's prompt -> fixed.
    r = await discuss([plan("C[c]1ccccc1"), plan("C[c]1ccccc1"), plan("[CH2]c1ccccc1")])
    print("A", summary(r))
    assert r["outcome"] == "passed" and len(r["rounds"]) == 3
    assert [s["action"] for s in r["supervisor"]] == ["PAUSE+EDIT"]
    assert "Write benzyl as [CH2]c1ccccc1." in SEEN_PLANNER_PROMPTS[-1]
    assert "Write benzyl" not in SEEN_PLANNER_PROMPTS[0]
    val = [m["content"] for m in r["history"] if m["role"] == "Validator"]
    assert "[CH2]c1ccccc1 + [H]" in val[0], val[0]  # the Validator shows the benzylic split

    # B: correct first time, critic WARN -> one round.
    r = await discuss([plan("[CH2]c1ccccc1")], [{"verdict": "WARN", "failed_checks": [], "warnings": ["no ZPE"]}],
                      critic=True)
    print("B", summary(r))
    assert r["outcome"] == "passed" and len(r["rounds"]) == 1 and r["rounds"][0]["critic"] == "WARN"

    # C: the mistake survives the correction -> TERMINATE, no plan.
    r = await discuss([plan("C[c]1ccccc1")] * 3, rounds=4)
    print("C", summary(r))
    assert r["outcome"] == "terminated" and r["spec"] is None
    assert [s["action"] for s in r["supervisor"]] == ["PAUSE+EDIT", "TERMINATE"]

    # D: critic FAIL on a valid plan, then the planner breaks it -> keep the
    # last plan the Validator passed.
    r = await discuss([plan("[CH2]c1ccccc1"), plan("CCO"), plan("CCO")],
                      [{"verdict": "FAIL", "failed_checks": ["rationale"], "warnings": []}], critic=True)
    print("D", summary(r))
    assert r["spec"] is not None and r["spec"]["fragments"][0]["smiles"] == "[CH2]c1ccccc1"
    assert r["outcome"].endswith("_kept_last_validated_plan")

    print("\nOK: plan discussion behaves as expected.")


if __name__ == "__main__":
    with provenance_session():
        asyncio.run(main())
