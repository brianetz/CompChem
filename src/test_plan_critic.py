"""Tests for agents/plan_critic.py (no Academy, no LLM, no network: stub planner).

Run from src/:  python -m unittest test_plan_critic
"""
import asyncio
import unittest

from agents.plan_critic import augment_verdict, critique_spec
from analysis import bde_plan_validator as v
from test_bde_plan_validator import water_spec


class StubPlanner:
    def __init__(self, replacement):
        self.replacement = replacement
        self.calls = []

    async def replan_campaign(self, **kwargs):
        self.calls.append(kwargs)
        return self.replacement


def run(coro):
    return asyncio.run(coro)


class PlanCritic(unittest.TestCase):
    def test_non_bde_spec_is_passed_through_untouched(self):
        planner = StubPlanner(None)
        spec = {"analysis_type": "single_point", "smiles": "O"}
        out, report = run(critique_spec(spec, planner))
        self.assertIs(out, spec)
        self.assertIsNone(report)
        self.assertEqual(planner.calls, [])

    def test_failing_plan_triggers_exactly_one_corrective_replan(self):
        bad = water_spec(parent="O")            # malformed parent -> FAIL (no RDKit needed)
        planner = StubPlanner(water_spec())
        out, report = run(critique_spec(bad, planner, model="m"))
        self.assertEqual(len(planner.calls), 1)
        call = planner.calls[0]
        self.assertEqual(call["model"], "m")
        self.assertIs(call["previous_spec"], bad)
        self.assertIn("parent_shape", call["feedback"])
        self.assertEqual(out, planner.replacement)
        if v.HAVE_RDKIT:
            self.assertNotEqual(report.verdict, v.FAIL, report.feedback())

    def test_retry_is_bounded_when_planner_keeps_failing(self):
        bad = water_spec(parent="O")
        planner = StubPlanner(bad)              # planner never fixes it
        _, report = run(critique_spec(bad, planner, max_retries=1))
        self.assertEqual(len(planner.calls), 1)
        self.assertEqual(report.verdict, v.FAIL)

    def test_max_retries_zero_never_calls_planner(self):
        bad = water_spec(parent="O")
        planner = StubPlanner(water_spec())
        run(critique_spec(bad, planner, max_retries=0))
        self.assertEqual(planner.calls, [])

    def test_passing_plan_makes_no_planner_call(self):
        planner = StubPlanner(None)
        _, report = run(critique_spec(water_spec(), planner))
        self.assertEqual(planner.calls, [])
        if v.HAVE_RDKIT:
            self.assertNotEqual(report.verdict, v.FAIL)


class PostRun(unittest.TestCase):
    PASS_VERDICT = {"qc_flag": "PASS", "qc_notes": None}

    def test_non_bde_results_unchanged(self):
        r = {"analysis_type": "single_point"}
        self.assertIs(augment_verdict(r, self.PASS_VERDICT), self.PASS_VERDICT)

    def test_negative_bde_escalates_to_fail_and_keeps_existing_notes(self):
        r = {"analysis_type": "bde", "bde_kcal_mol": -5.0, "parent": {}, "fragments": []}
        out = augment_verdict(r, {"qc_flag": "WARN", "qc_notes": "x_xtb_not_converged"})
        self.assertEqual(out["qc_flag"], "FAIL")
        self.assertIn("x_xtb_not_converged", out["qc_notes"])
        self.assertIn("bde_value", out["qc_notes"])

    def test_never_relaxes_an_existing_fail(self):
        r = {"analysis_type": "bde", "bde_kcal_mol": 100.0, "parent": {}, "fragments": []}
        out = augment_verdict(r, {"qc_flag": "FAIL", "qc_notes": "a"})
        self.assertEqual(out["qc_flag"], "FAIL")

    def test_plausible_bde_adds_no_notes(self):
        r = {"analysis_type": "bde", "bde_kcal_mol": 118.0, "parent": {}, "fragments": []}
        out = augment_verdict(r, self.PASS_VERDICT)
        self.assertEqual((out["qc_flag"], out["qc_notes"]), ("PASS", None))

    def test_failed_species_are_skipped_not_crashed(self):
        r = {"analysis_type": "bde", "bde_kcal_mol": None,
             "parent": {"status": "failed", "smiles": "O"}, "fragments": []}
        out = augment_verdict(r, {"qc_flag": "FAIL", "qc_notes": "parent_calculation_failed:e"})
        self.assertEqual(out["qc_flag"], "FAIL")

    def _oh(self, s2, atom_spin, mult=2):
        return {"analysis_type": "bde", "bde_kcal_mol": 118.0, "parent": {}, "fragments": [
            {"label": "OH", "smiles": "[OH]", "multiplicity": mult,
             "refinement": {"s2": s2, "atom_spin": atom_spin}}]}

    def test_clean_radical_spin_passes(self):
        out = augment_verdict(self._oh(0.7512, [0.93, 0.07]), self.PASS_VERDICT)
        self.assertEqual(out["qc_flag"], "PASS")

    def test_spin_contamination_fails_the_verdict(self):
        out = augment_verdict(self._oh(1.30, None), self.PASS_VERDICT)
        self.assertEqual(out["qc_flag"], "FAIL")
        self.assertIn("spin_contamination", out["qc_notes"])

    @unittest.skipUnless(v.HAVE_RDKIT, "RDKit not installed")
    def test_radical_on_wrong_atom_fails(self):
        out = augment_verdict(self._oh(0.7512, [0.05, 0.95]), self.PASS_VERDICT)
        self.assertEqual(out["qc_flag"], "FAIL")
        self.assertIn("radical_location", out["qc_notes"])

    def test_engine_without_spin_data_is_skipped(self):
        r = self._oh(None, None)
        self.assertEqual(augment_verdict(r, self.PASS_VERDICT)["qc_flag"], "PASS")

    @unittest.skipUnless(v.HAVE_RDKIT, "RDKit not installed")
    def test_connectivity_flags_dissociated_geometry(self):
        stretched = {"symbols": ["O", "H", "H"],
                     "positions": [[0, 0, 0.117], [0, 0.757, -0.469], [0, -3.5, -0.469]]}
        good = {"symbols": ["O", "H", "H"],
                "positions": [[0, 0, 0.117], [0, 0.757, -0.469], [0, -0.757, -0.469]]}
        mk = lambda xtb: {"analysis_type": "bde", "bde_kcal_mol": 118.0,
                          "parent": {"smiles": "O", "xtb": xtb}, "fragments": []}
        self.assertEqual(augment_verdict(mk(good), self.PASS_VERDICT)["qc_flag"], "PASS")
        self.assertEqual(augment_verdict(mk(stretched), self.PASS_VERDICT)["qc_flag"], "WARN")


if __name__ == "__main__":
    unittest.main()
