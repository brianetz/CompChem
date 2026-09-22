"""
GFN2-xTB via tblite's ASE calculator - the cheap/fast prescreening tier
named in the architecture plan (Phase 3). This runs real semiempirical
tight-binding quantum chemistry locally, not a mock.
"""
from __future__ import annotations

from ase.atoms import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite

from utils.units import ev_to_hartree
from storage.flowcept_store import flowcept_task


@flowcept_task
def optimize_with_xtb(atoms: Atoms, charge: int = 0, multiplicity: int = 1,
                       fmax: float = 0.05, steps: int = 200) -> dict:
    """Optimize a geometry with GFN2-xTB. Returns a JSON-safe dict (energy,
    convergence flag, final coordinates) rather than the ase.Atoms object
    itself, since Flowcept can't serialize arbitrary objects (see Phase 0
    lesson in the architecture doc) - callers that need the optimized
    geometry for a follow-up calculation should rebuild an Atoms object from
    'symbols'/'positions'.
    """
    working = atoms.copy()
    working.calc = TBLite(method="GFN2-xTB", charge=charge, multiplicity=multiplicity)

    opt = BFGS(working, logfile=None)
    converged = opt.run(fmax=fmax, steps=steps)

    energy_hartree = float(ev_to_hartree(working.get_potential_energy()))

    return {
        "engine": "xtb-GFN2",
        "energy_hartree": energy_hartree,
        "converged": bool(converged),
        "symbols": working.get_chemical_symbols(),
        "positions": working.get_positions().tolist(),
    }
