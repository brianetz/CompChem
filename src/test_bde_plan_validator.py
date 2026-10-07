"""Tests for bde_plan_validator (plain unittest; pytest collects these too).

Run from src/:  python -m pytest tests/test_bde_plan_validator.py
RDKit-dependent tests are skipped automatically when RDKit is not importable.
"""
import unittest

try:  # repo layout
    from analysis import bde_plan_validator as v
except ImportError:  # flat layout (file next to the module)
    import bde_plan_validator as v

RATIONALE_1 = (
    "Breaking one O-H bond of water (H-O-H) leaves the other O-H bond intact, giving hydroxyl "
    "radical ([OH], neutral doublet: 9 valence electrons, odd, consistent with multiplicity 2) plus a "
    "hydrogen atom ([H], neutral doublet: 1 valence electron, odd, consistent with multiplicity 2) "
    "- not full atomization into O + 2H. Parent water is closed-shell singlet (charge 0, multiplicity 1). "
    "Workflow: cheap GFN2-xTB geometry optimization/prescreen for parent and both fragments, followed by "
    "a DFT single-point refinement with B3LYP/def2-TZVP (a standard, reliable combination for main-group "
    "thermochemistry and radical species) to compute BDE = E(OH) + E(H) - E(H2O), converted to kcal/mol."
)
RATIONALE_1_FIXED = RATIONALE_1.replace("9 valence electrons", "9 total electrons (7 valence)")
RATIONALE_2 = (
    "Breaking one O-H bond in water leaves the other O-H bond intact, giving hydroxyl radical ([OH], "
    "doublet) plus a hydrogen atom ([H], doublet) rather than full atomization into O+H+H. A GFN2-xTB "
    "prescreen gives cheap geometries, followed by B3LYP/def2-TZVP single-point refinement, a standard "
    "accurate combination for BDEs of small main-group molecules."
)

WATER_COUNTS = dict(total_counts=[10, 9, 1], valence_counts=[8, 7, 1])


def water_spec(**overrides):
    spec = {
        "analysis_type": "bde",
        "parent": {"smiles": "O", "charge": 0, "multiplicity": 1},
        "fragments": [
            {"smiles": "[OH]", "charge": 0, "multiplicity": 2},
            {"smiles": "[H]", "charge": 0, "multiplicity": 2},
        ],
        "basis": "def2-TZVP",
        "xc": "B3LYP",
        "rationale": RATIONALE_1_FIXED,
    }
    spec.update(overrides)
    return spec


class PureLogic(unittest.TestCase):
    def test_electron_counts(self):
        self.assertEqual(v.electron_counts(["O", "H", "H"], 0)["total"], 10)
        self.assertEqual(v.electron_counts(["O", "H", "H"], 0)["valence"], 8)
        oh = v.electron_counts(["O", "H"], 0)
        self.assertEqual((oh["total"], oh["valence"]), (9, 7))
        self.assertEqual(v.electron_counts(["O", "H"], -1)["total"], 10)
        self.assertEqual(v.electron_counts(["Xx"], 0)["unsupported"], ["Xx"])

    def test_parity(self):
        self.assertTrue(v.parity_ok(9, 2))
        self.assertFalse(v.parity_ok(8, 2))   # the original "[O] written as O" failure shape
        self.assertTrue(v.parity_ok(10, 1))
        self.assertTrue(v.parity_ok(16, 3))   # O2 triplet
        self.assertFalse(v.parity_ok(10, 2))
        self.assertFalse(v.parity_ok(1, 4))   # more unpaired electrons than electrons

    def test_ground_state_table(self):
        self.assertEqual(v.ground_state_multiplicity("O"), 3)
        self.assertEqual(v.ground_state_multiplicity("H"), 2)
        self.assertIsNone(v.ground_state_multiplicity("Xx"))

    def test_rationale_1_mislabeled_valence_count_is_caught(self):
        issues = v.check_electron_claims(RATIONALE_1, **WATER_COUNTS)
        self.assertEqual(len(issues), 1)
        self.assertIn("9", issues[0])
        self.assertIn("TOTAL", issues[0])

    def test_rationale_1_fixed_is_clean(self):
        self.assertEqual(v.check_electron_claims(RATIONALE_1_FIXED, **WATER_COUNTS), [])

    def test_unqualified_wrong_number_is_caught(self):
        issues = v.check_electron_claims("OH has 11 electrons", **WATER_COUNTS)
        self.assertEqual(len(issues), 1)

    def test_spin_contamination_bands(self):
        status = lambda s2: v.check_spin("OH", 2, s2)[0].status
        self.assertEqual(status(0.7512), v.PASS)
        self.assertEqual(status(0.84), v.WARN)
        self.assertEqual(status(0.95), v.FAIL)
        self.assertEqual(v.check_spin("H2O", 1, None), [])

    def test_spin_location(self):
        by_name = lambda checks: {c.name: c.status for c in checks}
        good = v.check_spin("OH", 2, 0.751, atom_spin=[0.93, 0.07], expected_atoms=[0])
        self.assertEqual(by_name(good)["OH_radical_location"], v.PASS)
        wrong_site = v.check_spin("X", 2, 0.751, atom_spin=[0.05, 0.95], expected_atoms=[0])
        self.assertEqual(by_name(wrong_site)["X_radical_location"], v.FAIL)
        delocalized = v.check_spin("Y", 2, 0.751, atom_spin=[0.45, 0.55], expected_atoms=[0])
        self.assertEqual(by_name(delocalized)["Y_radical_location"], v.WARN)
        bad_sum = v.check_spin("Z", 2, 0.751, atom_spin=[0.5, 0.2], expected_atoms=[0])
        self.assertEqual(by_name(bad_sum)["Z_spin_sum"], v.WARN)

    def test_bde_value_band(self):
        self.assertEqual(v.check_bde_value(116.2).status, v.PASS)
        self.assertEqual(v.check_bde_value(-3.0).status, v.FAIL)
        self.assertEqual(v.check_bde_value(250.0).status, v.WARN)
        self.assertEqual(v.check_bde_value(None).status, v.FAIL)

    def test_malformed_parent_is_reported_not_raised(self):
        # the Phase 9 failure mode: parent given as a bare string
        report = v.critique_plan(water_spec(parent="O"))
        self.assertEqual(report.verdict, v.FAIL)
        self.assertIn("parent_shape", [c.name for c in report.checks])

    def test_report_qc_shape(self):
        report = v.critique_plan(water_spec(parent="O"))
        qc = report.to_qc()
        self.assertEqual(qc["qc_flag"], v.FAIL)
        self.assertTrue(qc["qc_notes"])


@unittest.skipUnless(v.HAVE_RDKIT, "RDKit not installed")
class WithRDKit(unittest.TestCase):
    def test_water_homolysis(self):
        # atom 0 is O, atoms 1 and 2 are H in AddHs ordering
        keys = sorted(v.homolysis_fragment_keys("O", 0, 1))
        self.assertEqual(keys, sorted([v.canonical_key("[OH]"), v.canonical_key("[H]")]))

    def test_ring_bond_gives_single_diradical(self):
        # a C-C bond in cyclohexane (atoms 0 and 1)
        self.assertEqual(len(v.homolysis_fragment_keys("C1CCCCC1", 0, 1)), 1)

    def test_infer_bond_for_water(self):
        hits = v.infer_cleaved_bonds("O", ["[OH]", "[H]"])
        self.assertEqual(len(hits), 2)                    # both O-H bonds
        self.assertEqual(len({h["class"] for h in hits}), 1)  # one symmetry class

    def test_impossible_fragments_are_rejected(self):
        # closed-shell ethanol + H atom cannot come from any single-bond homolysis
        self.assertEqual(v.infer_cleaved_bonds("CCO", ["CCO", "[H]"]), [])

    def test_ethanol_oh_cleavage_is_found_and_distinct_from_ch(self):
        oh = v.infer_cleaved_bonds("CCO", ["CC[O]", "[H]"])
        self.assertEqual(len(oh), 1)
        self.assertIn("O2", oh[0]["label"])
        ch = v.infer_cleaved_bonds("CCO", ["[CH2]CO", "[H]"])
        self.assertTrue(ch)
        self.assertNotEqual(oh[0]["class"], ch[0]["class"])

    def test_declared_bond_that_contradicts_fragments_fails(self):
        spec = {
            "analysis_type": "bde",
            "parent": {"smiles": "CCO", "charge": 0, "multiplicity": 1},
            "fragments": [
                {"smiles": "CC[O]", "charge": 0, "multiplicity": 2},
                {"smiles": "[H]", "charge": 0, "multiplicity": 2},
            ],
            "basis": "def2-TZVP", "xc": "B3LYP", "rationale": RATIONALE_1_FIXED,
            "bond_cleaved": {"atom_indices": [0, 3]},  # a C-H bond, not the O-H the fragments imply
        }
        report = v.critique_plan(spec)
        self.assertIn("bond_cleaved", [c.name for c in report.checks if c.status == v.FAIL])

    def test_original_O_for_OH_bug_is_a_fail(self):
        spec = water_spec(fragments=[
            {"smiles": "O", "charge": 0, "multiplicity": 2},
            {"smiles": "[H]", "charge": 0, "multiplicity": 2},
        ])
        report = v.critique_plan(spec)
        self.assertEqual(report.verdict, v.FAIL)
        names = [c.name for c in report.checks if c.status == v.FAIL]
        self.assertIn("fragment1_electron_parity", names)

    def test_good_plan_has_no_failures(self):
        report = v.critique_plan(water_spec())
        fails = [c for c in report.checks if c.status == v.FAIL]
        self.assertEqual(fails, [], report.feedback())
        bond = [c for c in report.checks if c.name == "bond_cleavage"][0]
        self.assertEqual(bond.status, v.PASS)

    def test_rationale_1_as_written_fails_on_electron_number(self):
        report = v.critique_plan(water_spec(rationale=RATIONALE_1))
        self.assertIn("rationale_electron_numbers", [c.name for c in report.checks if c.status == v.FAIL])

    def test_rationale_2_triggers_rubric_warnings(self):
        report = v.critique_plan(water_spec(rationale=RATIONALE_2))
        warned = {c.name for c in report.checks if c.status == v.WARN}
        self.assertTrue({"rationale_parent_state", "rationale_formula", "rationale_units"} <= warned)

    def test_atom_conservation_failure(self):
        spec = water_spec(fragments=[
            {"smiles": "[OH]", "charge": 0, "multiplicity": 2},
            {"smiles": "[H][H]", "charge": 0, "multiplicity": 1},
        ])
        report = v.critique_plan(spec)
        self.assertIn("atom_conservation", [c.name for c in report.checks if c.status == v.FAIL])

    def test_candidates_cover_symmetry_unique_ch_bonds(self):
        # propane: two unique C-H environments (primary, secondary)
        cands = v.enumerate_bde_candidates("CCC", elements=("C", "H"))
        self.assertEqual(len(cands), 2)

    def test_connectivity_check(self):
        # water geometry, planned as water: PASS; planned as OH + H together: mismatch
        symbols = ["O", "H", "H"]
        pos = [[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]]
        self.assertEqual(v.check_connectivity("water", symbols, pos, "O").status, v.PASS)
        stretched = [[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -3.5, -0.469]]
        self.assertNotEqual(v.check_connectivity("water", symbols, stretched, "O").status, v.PASS)

    def test_expected_radical_atoms(self):
        self.assertEqual(v.expected_radical_atoms("[OH]"), [0])


if __name__ == "__main__":
    unittest.main()
