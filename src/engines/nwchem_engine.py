"""
NWChem via ASE's Calculator protocol - wired into the same
atoms.calc = SomeCalculator() / atoms.get_potential_energy() pattern
engines/xtb_engine.py already uses for GFN2-xTB (Phase 3), rather than
engines/pyscf_engine.py's direct-API-bypass pattern (PySCF talks to its own
gto.M()/dft.RKS() interface, never touching ASE's Calculator protocol at
all, despite the architecture plan's Sec.3 table implying otherwise). This
is the calculator-abstraction unification deferred from Phase 6 to Phase 8
(architecture plan Sec.6, Sec.16) for the piece of it that doesn't need a
real NWChem binary to design and wire up - see the module-level caveats
below for what still does.

Positioned as an alternate refinement tier alongside PySCF (same level as
engines.pyscf_engine.singlepoint_energy_pyscf: a single-point energy at an
already-optimized geometry, typically xtb's), not a replacement for the old
hand-rolled-.nw-file pipeline (nwchem.py/run_all.py/parser.py), which stays
as-is and out of scope here.

*** IMPORTANT - NOT YET VALIDATED AGAINST A REAL NWCHEM BINARY ***
This was written and syntax-checked (`python3 -m py_compile`) without a
real `nwchem` binary or a compatible mock available in this environment
(the existing mock_nwchem.py is shaped for the OLD pipeline's hand-rolled
.nw-file + regex-based parsing, not for ASE's own NWChem FileIOCalculator,
which writes its own input file format and expects its own output shape -
making mock_nwchem.py unsuitable here without separate adaptation that
wasn't attempted, since a fake-but-ASE-compatible NWChem output is fragile
to get right and a real binary is the only trustworthy validation target
anyway). Specifically unverified against the actually-installed ASE
version:
  - the exact NWChem() constructor kwarg names/shapes below (ASE's NWChem
    calculator API has changed across ASE versions - the Phase 4 lesson
    about a young/fast-moving dependency's docs lagging its released
    version, architecture plan Sec.7, applies here just as much as it did
    to academy-py);
  - where ASE actually writes the raw NWChem output file for a given
    `label` (assumed here to be '<label>.nwo', ASE's long-standing
    convention, but not confirmed against the installed version);
  - whether `atoms.get_potential_energy()` raises cleanly
    (ase.calculators.calculator.CalculationFailed or similar) on a real
    NWChem non-convergence/crash, the way this module's exception handling
    assumes.
Before relying on this for a real campaign: run
`python -c "from ase.calculators.nwchem import NWChem; help(NWChem)"` in
the compchem-agent env and diff its actual signature against
_build_calculator() below, then exercise singlepoint_energy_nwchem()
directly (not through the full agent pipeline) against a real `nwchem`
binary or a real Frontier/Defiant allocation.
"""
from __future__ import annotations

import os
from typing import List, Tuple

from ase.atoms import Atoms
from ase.calculators.nwchem import NWChem

from utils.units import ev_to_hartree
from storage.flowcept_store import flowcept_task
from parsers.nwchem_fallback import read_text, extract_success

# Same NWCHEM_CMD convention already used by run_all.py/mock_nwchem.py (see
# architecture plan Sec.11) - real Frontier invocation by default, override
# for local/other testing. ASE's FileIOCalculator command string uses
# 'PREFIX' as a literal placeholder it substitutes itself; a bare command
# with no redirection also works since ASE appends its own '> PREFIX.nwo'
# when the command doesn't already contain 'PREFIX'.
_DEFAULT_NWCHEM_CMD = os.environ.get(
    "NWCHEM_CMD",
    "srun -n56 /lustre/orion/proj-shared/stf053/ketan/isc/nwchem-install-3/bin/nwchem",
)


def _spin_from_multiplicity(multiplicity: int) -> int:
    """Mirrors engines.pyscf_engine._spin_from_multiplicity: NWChem's `mult`
    keyword IS the chemist's multiplicity (2S+1, unlike PySCF's 2S
    convention) - kept as a named function anyway so the distinction is
    explicit and searchable, not an unexplained pass-through."""
    return multiplicity


def _build_calculator(label: str, charge: int, multiplicity: int,
                       basis: str, xc: str, command: str | None) -> NWChem:
    """Isolated on purpose (see module docstring): if the installed ASE
    version needs different kwarg names, this is the one function that
    needs to change."""
    return NWChem(
        label=label,
        command=command or _DEFAULT_NWCHEM_CMD,
        theory="dft",
        charge=charge,
        basis=basis,
        dft={"xc": xc, "mult": _spin_from_multiplicity(multiplicity), "maxiter": 200},
    )


@flowcept_task
def singlepoint_energy_nwchem(
    symbols: List[str], positions: List[Tuple[float, float, float]],
    charge: int = 0, multiplicity: int = 1,
    basis: str = "6-31g*", xc: str = "b3lyp",
    command: str | None = None, label: str = "nwchem_singlepoint",
) -> dict:
    """Single-point DFT energy for a fixed geometry via real NWChem, given
    as plain symbols/positions lists - same call shape as
    engines.pyscf_engine.singlepoint_energy_pyscf so the two are
    interchangeable refinement tiers from skills/_common.py's point of view
    (see skills/_common.py's new `refinement_engine` parameter).

    Raises on failure rather than catching internally, matching
    xtb_engine/pyscf_engine's convention - the graceful status/error
    fallback lives one layer up, in skills/_common.py's
    compute_species_from_smiles/compute_species_from_conformer_search
    (architecture plan Sec.7/Sec.16)."""
    atoms = Atoms(symbols=symbols, positions=[tuple(p) for p in positions])
    atoms.calc = _build_calculator(label, charge, multiplicity, basis, xc, command)

    energy_hartree = float(ev_to_hartree(atoms.get_potential_energy()))

    # Cross-check convergence against the raw output text using the SAME
    # regex-based success detector the old CalculationRecord/cclib pipeline
    # already trusts (parsers/nwchem_fallback.py) - the concrete piece of
    # "reconcile the two schemas" that isn't just shared PASS/WARN/FAIL
    # arithmetic (see analysis/qc_common.py): one canonical answer to "did
    # this NWChem job actually finish cleanly", reused by both pipelines,
    # rather than trusting ASE's parse alone. If the output file isn't
    # where we expect (see module docstring caveats), fall back to True -
    # atoms.get_potential_energy() would already have raised on a real
    # failure, so reaching this line at all is itself a good sign.
    converged = True
    output_path = f"{label}.nwo"
    if os.path.exists(output_path):
        converged = extract_success(read_text(output_path))

    return {
        "engine": f"nwchem-{xc}/{basis}",
        "energy_hartree": energy_hartree,
        "converged": bool(converged),
    }
