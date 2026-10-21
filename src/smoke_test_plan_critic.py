"""
Smoke test for the bde-plan-critic skill wiring: the deterministic
validator (analysis/bde_plan_validator.py) on known good/bad plans, and the
planner's retry loop with the LLM critic enabled, using a stub backend in
place of a real model.

Run:
    python src/smoke_test_plan_critic.py
"""
from __future__ import annotations

import os

from analysis.bde_plan_validator import critique_plan, enumerate_bde_candidates
import llm.planner_llm as planner
from storage.flowcept_store import provenance_session

MQ = "What's the O-H bond dissociation energy of methanol (CH3OH) in kcal/mol?"
TQ = "What's the benzylic C-H bond dissociation energy of toluene in kcal/mol?"
H = {"smiles": "[H]", "charge": 0, "multiplicity": 2}


def spec(parent, frags, basis="6-31g*", rationale="r"):
    return {"analysis_type": "bde", "parent": {"smiles": parent, "charge": 0, "multiplicity": 1},
            "fragments": frags, "basis": basis, "xc": "b3lyp", "rationale": rationale}


def rad(smiles, mult=2, charge=0):
    return {"smiles": smiles, "charge": charge, "multiplicity": mult}


def main():
    # --- validator ------------------------------------------------------
    cases = [
        ("methanol C[O] + [H]", spec("CO", [rad("C[O]"), H]), MQ, "PASS", None),
        ("methanol C[O:1] + [H]", spec("CO", [rad("C[O:1]"), H]), MQ, "PASS", None),
        ("methanol [CH3O]", spec("CO", [rad("[CH3O]"), H]), MQ, "FAIL", "smiles_parse"),
        ("methanol CO doublet", spec("CO", [rad("CO"), H]), MQ, "FAIL", "parity"),
        ("methanol C-H break", spec("CO", [rad("[CH2]O"), H]), MQ, "FAIL", "requested_bond"),
        ("methanol OH + H", spec("CO", [rad("[OH]"), H]), MQ, "FAIL", "atom_balance"),
        ("toluene benzylic", spec("Cc1ccccc1", [rad("[CH2]c1ccccc1"), H]), TQ, "PASS", None),
        ("toluene ring C-H", spec("Cc1ccccc1", [rad("Cc1ccc[c]c1"), H]), TQ, "FAIL", "site"),
        ("toluene, no site word", spec("Cc1ccccc1", [rad("Cc1ccc[c]c1"), H]), "C-H BDE of toluene?", "WARN", None),
        ("singlet O atom", spec("O", [rad("[O]", mult=1), rad("[H][H]", mult=1)]), None, "FAIL", "atom_ground_state"),
        ("charge mismatch", spec("CO", [rad("C[O-]", mult=1), rad("[H+]", mult=1, charge=1)]), MQ, "FAIL", "formal_charge"),
    ]
    for name, sp, q, want, check in cases:
        rep = critique_plan(sp, q)
        got_checks = [f["check"] for f in rep["fails"]]
        print(f"{name:24} {rep['verdict']:4} fails={got_checks} warns={[w['check'] for w in rep['warns']]}")
        assert rep["verdict"] == want, (name, rep)
        assert check is None or check in got_checks, (name, got_checks)

    bonds = sorted(c["bond"] for c in enumerate_bde_candidates("Cc1ccccc1"))
    print("toluene symmetry-unique bonds:", bonds)
    assert bonds.count("C-H") == 4, bonds  # benzylic + ortho/meta/para ring C-H

    # A C-H break of methanol for an O-H question names the right answer in the fix.
    fix = critique_plan(spec("CO", [rad("[CH3]"), rad("[OH]")]), MQ)["fails"][0]["fix"]
    assert "C[O] + [H]" in fix or "[H] + C[O]" in fix, fix

    # --- planner loop with the critic on (stub backend) -----------------
    skill = os.path.join(os.path.dirname(__file__), "bde-plan-critic.md")
    os.environ["COMPCHEM_PLAN_CRITIC_SKILL"] = skill
    plan_replies = iter([
        {**spec("CO", [rad("[CH3O]"), H])},            # attempt 1: validator FAIL
        {**spec("CO", [rad("C[O]"), H])},              # attempt 2: valid
    ])
    critic_replies = iter([{"verdict": "WARN", "failed_checks": [], "warnings": ["no ZPE statement"]}])

    def stub_attempt(model, system_prompt, tools, messages):
        if tools[0]["name"] == "emit_plan_review":
            assert "BDE plan critic" in system_prompt
            return next(critic_replies), None, None
        return next(plan_replies), {"role": "assistant", "content": "x"}, "ref"

    def stub_retry(messages, entry, ref, error_text):
        return messages + [{"role": "user", "content": error_text}]

    planner._BACKENDS["stub"] = (stub_attempt, stub_retry)
    out = planner._generate_spec_via_llm(
        question_for_spec=MQ, initial_messages=[{"role": "user", "content": MQ}],
        model="m", backend="stub", provenance_action="test",
    )
    stages = [(r["attempt"], r["stage"], r["verdict"]) for r in out["plan_review_log"]]
    print("review log:", stages)
    assert stages == [(1, "validator", "FAIL"), (2, "validator", "PASS"), (2, "critic", "WARN")], stages

    # Critic FAIL on attempt 1 earns the retry; FAIL on the last attempt keeps the plan.
    plan_replies = iter([spec("CO", [rad("C[O]"), H]), spec("CO", [rad("C[O]"), H])])
    critic_replies = iter([
        {"verdict": "FAIL", "failed_checks": ["rationale says '7 total electrons' for OH"], "warnings": []},
        {"verdict": "FAIL", "failed_checks": ["still wrong"], "warnings": []},
    ])
    out = planner._generate_spec_via_llm(
        question_for_spec=MQ, initial_messages=[{"role": "user", "content": MQ}],
        model="m", backend="stub", provenance_action="test",
    )
    stages = [(r["attempt"], r["stage"], r["verdict"]) for r in out["plan_review_log"]]
    print("review log:", stages)
    assert stages[-1] == (2, "critic", "FAIL") and out["fragments"][0]["smiles"] == "C[O]", stages

    # Stringified arrays from the critic are decoded.
    critic_replies = iter([{"verdict": "pass", "failed_checks": "[]", "warnings": '["a", "b"]'}])
    plan_replies = iter([spec("CO", [rad("C[O]"), H])])
    out = planner._generate_spec_via_llm(
        question_for_spec=MQ, initial_messages=[{"role": "user", "content": MQ}],
        model="m", backend="stub", provenance_action="test",
    )
    crit = out["plan_review_log"][-1]
    assert crit["verdict"] == "PASS" and crit["warnings"] == ["a", "b"], crit

    print("\nOK: plan validator and LLM plan critic behave as expected.")


if __name__ == "__main__":
    with provenance_session():
        main()
