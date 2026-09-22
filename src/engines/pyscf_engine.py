"""
PySCF DFT single-point energy - a higher-accuracy tier than xtb, meant to be
run at a geometry already optimized more cheaply (e.g. by
engines.xtb_engine.optimize_with_xtb). No geometry optimization or frequency
analysis here yet (see architecture plan Phase 3 scope note) - just a real,
local, converged DFT energy at a fixed geometry.
"""
from __future__ import annotations

from typing import List, Tuple

from pyscf import gto, dft

from storage.flowcept_store import flowcept_task


def _spin_from_multiplicity(multiplicity: int) -> int:
    """PySCF's `spin` is 2S (number of unpaired electrons), not the chemist's
    multiplicity (2S+1) - e.g. a doublet radical (multiplicity=2) is spin=1."""
    return multiplicity - 1


@flowcept_task
def singlepoint_energy_pyscf(symbols: List[str], positions: List[Tuple[float, float, float]],
                              charge: int = 0, multiplicity: int = 1,
                              basis: str = "6-31g*", xc: str = "b3lyp") -> dict:
    """Single-point DFT energy for a fixed geometry, given as plain
    symbols/positions lists (JSON-safe inputs, and directly compatible with
    what xtb_engine.optimize_with_xtb returns)."""
    atom_spec = list(zip(symbols, [tuple(p) for p in positions]))
    mol = gto.M(atom=atom_spec, basis=basis, charge=charge,
                spin=_spin_from_multiplicity(multiplicity), unit="Angstrom")

    mf = dft.UKS(mol) if mol.spin != 0 else dft.RKS(mol)
    mf.xc = xc
    energy_hartree = mf.kernel()

    return {
        "engine": f"pyscf-{xc}/{basis}",
        "energy_hartree": float(energy_hartree),
        "converged": bool(mf.converged),
    }
