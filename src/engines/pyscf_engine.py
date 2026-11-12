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

    result = {
        "engine": f"pyscf-{xc}/{basis}",
        "energy_hartree": float(energy_hartree),
        "converged": bool(mf.converged),
    }
    # Additive diagnostics for the opt-in post-run spin checks
    # (agents/plan_critic.py): <S^2> and per-atom Mulliken spin populations
    # (atom order = input order). Only for open-shell (UKS) runs; never allowed
    # to affect the energy result.
    if mol.spin != 0:
        try:
            result.update(_spin_diagnostics(mol, mf))
        except Exception:
            pass
    return result


def _spin_diagnostics(mol, mf) -> dict:
    s2, _ = mf.spin_square()
    dm_a, dm_b = mf.make_rdm1()
    spin_per_ao = ((dm_a - dm_b) * mf.get_ovlp()).sum(axis=0)  # diag of (Da-Db)S
    atom_spin = [float(spin_per_ao[p0:p1].sum()) for _, _, p0, p1 in mol.aoslice_by_atom()]
    # 'atom_spin' (plain list, input atom order) is what check_spin consumes;
    # 'atom_spin_labeled' is the same data keyed "<element><index>" for humans.
    labeled = {f"{mol.atom_symbol(i)}{i}": round(x, 4) for i, x in enumerate(atom_spin)}
    return {"s2": float(s2), "atom_spin": atom_spin, "atom_spin_labeled": labeled}
