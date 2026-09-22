"""
Phase 0 smoke test: prove Flowcept actually captures task provenance and that
it can be queried back out, using the new parse_and_analyze building blocks
(no real NWChem output files needed - uses synthetic CalculationRecords).

One-time setup:
    pip install flowcept cclib ase pandas pyarrow
    flowcept --init-settings
    flowcept --config-profile full-offline -y

Run:
    python src/smoke_test_flowcept.py

Note: @flowcept_task-decorated functions should only return JSON-safe values
(dataclasses, dicts with str keys, lists, primitives) or None - Flowcept's
orjson-based serializer can't handle a raw pandas DataFrame's default integer
row-index as a captured return value. That's why build_bde_table (returns a
DataFrame) is left undecorated, while save_bde_table (writes CSV, returns
None) carries the decorator instead.
"""
from __future__ import annotations

import pandas as pd

from storage.flowcept_store import provenance_session, read_captured_tasks
from storage.records import CalculationRecord
from analysis.qc_checks import apply_qc_checks
from analysis.bde import save_bde_table


def make_fake_record(calc_id, molecule_id, role, energy, fragment_id=None, broken_bond=None):
    return CalculationRecord(
        calc_id=calc_id, molecule_id=molecule_id, fragment_id=fragment_id, role=role,
        smiles=None, parent_smiles=None, broken_bond=broken_bond, formula="CH4",
        charge=0, multiplicity=1 if role == "parent" else 2, n_atoms=5,
        software="synthetic", job_type="opt_freq", method="B3LYP", basis="6-31G*",
        input_file=None, output_file=f"{calc_id}.out", xyz_file=None,
        job_id=None, hostname=None, timestamp=None,
        success=True, optimization_converged=True, imaginary_freq_count=0,
        scf_energy_hartree=energy,
    )


def main():
    buffer_path = "flowcept_buffer.jsonl"
    bde_outfile = "bde_smoketest.csv"

    with provenance_session():
        parent = apply_qc_checks(make_fake_record("parent", "mol1", "parent", -40.5))
        f1 = apply_qc_checks(make_fake_record("f1", "mol1", "fragment", -20.1, "f1", "C-H_0"))
        f2 = apply_qc_checks(make_fake_record("f2", "mol1", "fragment", -20.2, "f2", "C-H_0"))
        save_bde_table([parent, f1, f2], bde_outfile)

    print("QC flags:", parent.qc_flag, f1.qc_flag, f2.qc_flag)
    print(pd.read_csv(bde_outfile))

    print("\n--- reading back captured Flowcept tasks ---")
    df = read_captured_tasks(buffer_path)
    if df is None or len(df) == 0:
        print(f"No captured tasks found at ./{buffer_path} - check that "
              f"'flowcept --config-profile full-offline -y' was run, and that "
              f"this script is being run from the same directory as the buffer file.")
    else:
        print(f"Captured {len(df)} task records.")
        cols = [c for c in ("activity_id", "status", "started_at", "ended_at") if c in df.columns]
        print(df[cols] if cols else df.head())


if __name__ == "__main__":
    main()
