"""
Phase 3 smoke test: prove PySCF and xtb (via tblite) both work locally as
real alternate QC engines - no NWChem needed. Demonstrates the intended
tiered workflow: RDKit provides a starting geometry, xtb (cheap) optimizes
it, PySCF (higher accuracy) computes a single-point energy at that geometry.

Run:
    python src/smoke_test_engines.py
"""
from __future__ import annotations

from conformer import get_lowest_energy_conformer
from engines.rdkit_to_ase import rdkit_mol_to_ase_atoms
from engines.xtb_engine import optimize_with_xtb
from engines.pyscf_engine import singlepoint_energy_pyscf
from storage.flowcept_store import provenance_session


def main():
    with provenance_session():
        conf_result = get_lowest_energy_conformer(smiles="O", num_confs=5)  # water
        atoms = rdkit_mol_to_ase_atoms(conf_result["lowest_energy_conformer"])
        print(f"Starting geometry: {len(atoms)} atoms, {atoms.get_chemical_formula()}")

        xtb_result = optimize_with_xtb(atoms, charge=0, multiplicity=1)
        print(f"xtb (GFN2-xTB) optimized energy: {xtb_result['energy_hartree']:.6f} hartree "
              f"(converged={xtb_result['converged']})")

        pyscf_result = singlepoint_energy_pyscf(
            symbols=xtb_result["symbols"],
            positions=xtb_result["positions"],
            charge=0, multiplicity=1, basis="6-31g*", xc="b3lyp",
        )
        print(f"PySCF (B3LYP/6-31G*) single-point energy at xtb geometry: "
              f"{pyscf_result['energy_hartree']:.6f} hartree (converged={pyscf_result['converged']})")

    assert xtb_result["converged"], "xtb geometry optimization should converge for water"
    assert xtb_result["energy_hartree"] < 0, "xtb energy should be negative"
    assert pyscf_result["converged"], "PySCF SCF should converge for water"
    assert pyscf_result["energy_hartree"] < 0, "PySCF energy should be negative"
    print("\nOK: both engines produced real, converged energies locally - no NWChem involved.")


if __name__ == "__main__":
    main()
