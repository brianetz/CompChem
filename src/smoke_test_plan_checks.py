"""
Smoke test for the plan checks (llm/plan_checks.py) and the BDE results
critic's chemistry checks (skills/bde.review), replayed on the tool calls
and results seen in the 2026-10-19 methanol runs (nemotron-3-nano and
Llama 3.1 405B). No LLM or compute needed.

Run:
    python src/smoke_test_plan_checks.py
"""
from __future__ import annotations

from llm.planner_llm import _spec_from_raw
from llm.plan_checks import check_spec, normalize_raw
from skills.bde import review

Q = ("What's the O-H bond dissociation energy of methanol (CH3OH) in kcal/mol, "
     "using a cheap xtb prescreen and a DFT refinement?")
BASE = {"analysis_type": "bde", "basis": "6-31g*", "xc": "b3lyp", "rationale": "r"}
H = {"smiles": "[H]", "charge": 0, "multiplicity": 2, "label": "H atom"}
PARENT = {"smiles": "CO", "charge": 0, "multiplicity": 1}


def plan(raw):
    """(spec or None, repairs, problems) - what _generate_spec_via_llm sees."""
    fixed, repairs = normalize_raw(raw)
    try:
        spec = _spec_from_raw(fixed, Q, "test")
    except ValueError as exc:
        return None, repairs, [f"shape: {exc}"]
    return spec, repairs, check_spec(spec, Q)


def frag(smiles, mult=2, label="methoxy radical"):
    return {"smiles": smiles, "charge": 0, "multiplicity": mult, "label": label}


def main():
    # --- shape repairs (Llama 3.1 405B) ---------------------------------
    cases = {
        "JSON string": {**BASE, "parent": '{"smiles": "CO", "charge": 0, "multiplicity": 1}',
                        "fragments": '[{"smiles": "C[O]", "charge": 0, "multiplicity": 2}, '
                                     '{"smiles": "[H]", "charge": 0, "multiplicity": 2}]'},
        "Python repr": {**BASE, "parent": "{'smiles': 'CO', 'charge': 0, 'multiplicity': 1}",
                        "fragments": "[{'smiles': '[H]', 'charge': 0, 'multiplicity': 2}, "
                                     "{'smiles': 'C[O]', 'charge': 0, 'multiplicity': 2}]"},
        "flattened parent": {**BASE, "parent_smiles": "CO", "parent_charge": 0,
                             "parent_multiplicity": 1, "fragments": [frag("C[O]"), H]},
    }
    for name, raw in cases.items():
        spec, repairs, problems = plan(raw)
        print(f"{name:18} repairs={repairs} problems={problems}")
        assert spec is not None and repairs and not problems, name

    truncated = {**BASE, "parent": PARENT, "fragments": "[{'smiles': '[H]"}
    spec, _, problems = plan(truncated)
    assert spec is None and "JSON array of objects (not a string)" in problems[0], problems

    # Bare 'H' is mapped to '[H]' instead of failing at RDKit.
    spec, _, problems = plan({**BASE, "parent": PARENT,
                              "fragments": [frag("C[O]"), {**H, "smiles": "H"}]})
    assert not problems and spec["fragments"][1]["smiles"] == "[H]", problems

    # --- chemistry problems (both models) -------------------------------
    bad = {
        "[CH3O] (invalid)": ([frag("[CH3O]"), H], "is not valid"),
        "CO as doublet": ([frag("CO"), H], "cannot have multiplicity 2"),
        "OH + H, no CH3": ([frag("[OH]", label="OH"), H], "add up to"),
        "C-H break": ([frag("[CH2]O", label="hydroxymethyl"), H], "C-H bond"),
        "C-O break": ([frag("[CH3]", label="methyl"), frag("[OH]", label="OH")], "O-C bond"),
        "formula SMILES": ([frag("CH3O"), H], "is not valid"),
    }
    for name, (frags, expect) in bad.items():
        spec, _, problems = plan({**BASE, "parent": PARENT, "fragments": frags})
        print(f"{name:18} -> {problems[0][:110] if problems else 'NO PROBLEM'}")
        assert any(expect in p for p in problems), (name, problems)

    for good in ("C[O]", "[O]C", "C[O:1]"):
        _, _, problems = plan({**BASE, "parent": PARENT, "fragments": [frag(good), H]})
        assert not problems, (good, problems)
    _, _, problems = plan({**BASE, "parent": {"smiles": "O", "charge": 0, "multiplicity": 1},
                           "fragments": [frag("[OH]", label="OH"), H]})
    assert not problems, problems  # water O-H still fine, for an O-H question

    # --- results critic -------------------------------------------------
    ok = {"status": "ok", "xtb": {"converged": True},
          "refinement": {"converged": True, "energy_hartree": -1.0}}

    def results(bde, balance=None, bond="O-H"):
        return {"parent": ok, "fragments": [{**ok, "label": "a"}, {**ok, "label": "b"}],
                "bde_kcal_mol": bde, "atom_balance": balance, "bond_type": bond,
                "requested_bonds": ["O-H"]}

    verdicts = {
        "good": (review(results(103.99)), "PASS"),
        "24 763 kcal/mol": (review(results(24763.28, "fragments add up to HO2 but ...", None)), "FAIL"),
        "C-H for O-H question": (review(results(101.75, bond="C-H")), "FAIL"),
    }
    for name, (v, want) in verdicts.items():
        print(f"critic {name:22} {v['qc_flag']:4} {v['qc_notes']}")
        assert v["qc_flag"] == want, (name, v)

    print("\nOK: plan checks and results critic behave as expected.")


if __name__ == "__main__":
    main()
