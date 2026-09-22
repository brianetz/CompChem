"""
Phase 2 smoke test: prove build_bde_table/save_bde_table actually gate on
qc_flag (PASS only) instead of silently including WARN/FAIL records, and
that excluded records are reported (not just dropped).

Run:
    python src/smoke_test_critic.py
"""
from __future__ import annotations

from storage.records import CalculationRecord
from storage.flowcept_store import provenance_session
from analysis.qc_checks import apply_qc_checks
from analysis.bde import build_bde_table, save_bde_table


def make_record(calc_id, molecule_id, role, energy=None, fragment_id=None,
                 broken_bond=None, success=True, n_atoms=5, imaginary_freq_count=0):
    return CalculationRecord(
        calc_id=calc_id, molecule_id=molecule_id, fragment_id=fragment_id, role=role,
        smiles=None, parent_smiles=None, broken_bond=broken_bond, formula="CH4",
        charge=0, multiplicity=1 if role == "parent" else 2, n_atoms=n_atoms,
        software="synthetic", job_type="opt_freq", method="B3LYP", basis="6-31G*",
        input_file=None, output_file=f"{calc_id}.out", xyz_file=None,
        job_id=None, hostname=None, timestamp=None,
        success=success, optimization_converged=True, imaginary_freq_count=imaginary_freq_count,
        scf_energy_hartree=energy,
    )


def main():
    with provenance_session():
        records = [
            # Clean PASS bond: parent + 2 good fragments -> should end up in the BDE table
            apply_qc_checks(make_record("parent", "mol1", "parent", -40.5)),
            apply_qc_checks(make_record("good_1", "mol1", "fragment", -20.1, "good_1", "C-H_0")),
            apply_qc_checks(make_record("good_2", "mol1", "fragment", -20.2, "good_2", "C-H_0")),

            # WARN bond: one fragment has an imaginary frequency (not a true minimum)
            apply_qc_checks(make_record("warn_1", "mol1", "fragment", -20.0, "warn_1", "C-C_1", imaginary_freq_count=1)),
            apply_qc_checks(make_record("warn_2", "mol1", "fragment", -20.3, "warn_2", "C-C_1")),

            # FAIL bond: one fragment's job failed outright
            apply_qc_checks(make_record("fail_1", "mol1", "fragment", None, "fail_1", "O-H_2", success=False)),
            apply_qc_checks(make_record("fail_2", "mol1", "fragment", -20.4, "fail_2", "O-H_2")),
        ]

        print("qc_flags:", {r.calc_id: r.qc_flag for r in records})

        bde_df, excluded_df = build_bde_table(records)

        print("\n--- BDE table (PASS only) ---")
        print(bde_df[["bond_label", "bde_electronic_kcal_mol"]] if len(bde_df) else "(empty)")

        print("\n--- excluded records ---")
        print(excluded_df[["calc_id", "qc_flag", "qc_notes"]])

        included_bonds = set(bde_df["bond_label"]) if len(bde_df) else set()
        assert "C-H_0" in included_bonds, "expected the clean bond to be included"
        assert "C-C_1" not in included_bonds, "WARN bond (imaginary freq) should be excluded"
        assert "O-H_2" not in included_bonds, "FAIL bond should be excluded"
        print("\nOK: qc_flag gating behaves as expected.")

        save_bde_table(records, "critic_smoketest.csv")


if __name__ == "__main__":
    main()
