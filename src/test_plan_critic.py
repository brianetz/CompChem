"""Tests for agents/plan_critic.py (no Academy, no LLM, no network: stub planner).

Run from src/:  python -m unittest test_plan_critic
"""
import asyncio
import unittest

from agents.plan_critic import critique_spec
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


if __name__ == "__main__":
    unittest.main()
