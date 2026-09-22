"""
Shared helpers for analysis skills (src/skills/*.py) - small, reusable
building blocks composed by more than one skill. Kept out of any single
skill module and out of the agent classes, which should stay skill-agnostic
(Phase 6 - architecture plan Sec.16).
"""
from __future__ import annotations

import os
import re

from conformer import get_lowest_energy_conformer
from engines.rdkit_to_ase import rdkit_mol_to_ase_atoms, smiles_to_ase_atoms
from engines.xtb_engine import optimize_with_xtb
from engines.pyscf_engine import singlepoint_energy_pyscf
from engines.nwchem_engine import singlepoint_energy_nwchem

# Refinement tier dispatch (Phase 8, partial - architecture plan Sec.19):
# pyscf and nwchem are interchangeable single-point refinement engines run
# at whatever geometry xtb already optimized. 'pyscf' stays the default so
# every existing caller (Phases 5-9) is unaffected; NWChem is opt-in via
# refinement_engine='nwchem' until it's been validated against a real
# binary (see engines/nwchem_engine.py's module docstring) and, later,
# made choosable by the LLM/rule-based planners rather than only by a
# direct caller of these two functions.
_REFINEMENT_ENGINES = {
    "pyscf": singlepoint_energy_pyscf,
    "nwchem": singlepoint_energy_nwchem,
}


def _run_refinement(engine_name: str, symbols, positions, charge, multiplicity, basis, xc) -> dict:
    try:
        engine_fn = _REFINEMENT_ENGINES[engine_name]
    except KeyError:
        raise ValueError(
            f"Unknown refinement_engine {engine_name!r} - choose one of "
            f"{sorted(_REFINEMENT_ENGINES)}"
        ) from None
    return engine_fn(symbols, positions, charge=charge, multiplicity=multiplicity, basis=basis, xc=xc)


def compute_species_from_smiles(
    smiles: str, charge: int, multiplicity: int, basis: str, xc: str,
    refinement_engine: str = "pyscf",
) -> dict:
    """xtb -> refinement engine (PySCF by default, NWChem opt-in - see
    _REFINEMENT_ENGINES above) for one species given directly by SMILES (a
    single RDKit-embedded starting geometry, no multi-conformer
    force-field search) - the right tool for a small/rigid/radical species
    (a BDE fragment, a redox state, ...) where xtb will re-optimize the
    geometry regardless of how good the starting guess is.

    Never raises: an engine failure (including a chemically invalid
    charge/multiplicity combination an LLM planner proposed) comes back as
    an explicit status/error field instead of crashing the caller - see
    architecture plan Sec.7/Sec.15."""
    try:
        atoms = smiles_to_ase_atoms(smiles)
        xtb_result = optimize_with_xtb(atoms, charge=charge, multiplicity=multiplicity)
        refinement_result = _run_refinement(
            refinement_engine, xtb_result["symbols"], xtb_result["positions"],
            charge, multiplicity, basis, xc,
        )
        return {
            "smiles": smiles, "charge": charge, "multiplicity": multiplicity,
            "status": "ok", "error": None, "refinement_engine": refinement_engine,
            "xtb": xtb_result, "refinement": refinement_result,
        }
    except Exception as exc:
        return {
            "smiles": smiles, "charge": charge, "multiplicity": multiplicity,
            "status": "failed", "error": f"{type(exc).__name__}: {exc}",
            "refinement_engine": refinement_engine,
            "xtb": None, "refinement": None,
        }


def compute_species_from_conformer_search(
    smiles: str, charge: int, multiplicity: int, basis: str, xc: str,
    refinement_engine: str = "pyscf",
) -> dict:
    """xtb -> refinement engine (PySCF by default, NWChem opt-in) for one
    species via conformer.get_lowest_energy_conformer's full multi-conformer
    MMFF/UFF search - the right tool for a flexible 'parent' molecule, where
    a good starting geometry actually matters more than for a small rigid
    fragment (see compute_species_from_smiles above).

    Same never-raises contract as compute_species_from_smiles."""
    try:
        conf_result = get_lowest_energy_conformer(smiles)
        mol = conf_result["lowest_energy_conformer"]
        atoms = rdkit_mol_to_ase_atoms(mol)
        xtb_result = optimize_with_xtb(atoms, charge=charge, multiplicity=multiplicity)
        refinement_result = _run_refinement(
            refinement_engine, xtb_result["symbols"], xtb_result["positions"],
            charge, multiplicity, basis, xc,
        )
        return {
            "smiles": smiles, "charge": charge, "multiplicity": multiplicity,
            "status": "ok", "error": None, "refinement_engine": refinement_engine,
            "conformer_energy_ff": conf_result["energy"],
            "xtb": xtb_result, "refinement": refinement_result,
        }
    except Exception as exc:
        return {
            "smiles": smiles, "charge": charge, "multiplicity": multiplicity,
            "status": "failed", "error": f"{type(exc).__name__}: {exc}",
            "refinement_engine": refinement_engine,
            "conformer_energy_ff": None, "xtb": None, "refinement": None,
        }


def safe_filename(text: str) -> str:
    """Turn an arbitrary label/SMILES (which may contain '[', ']', '/', etc,
    all unsafe in a filename) into a filesystem-safe filename stem for a
    rendered structure image (Phase 7 - architecture plan Sec.17)."""
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_")
    return stem or "species"


def render_species_image(label: str, xtb_result: dict | None, image_dir: str) -> str | None:
    """Render one species' xtb-optimized geometry to a PNG (Phase 7 -
    engines/render_engine.py, PyMOL) and return a ready-to-embed markdown
    image line, or None if rendering wasn't possible.

    A missing geometry (species failed before xtb ran), a missing PyMOL
    install, or any rendering failure are all reasons to just omit the
    image - never a reason to fail the whole report, matching the
    graceful-degradation approach used throughout this repo (architecture
    plan Sec.7)."""
    if not xtb_result or not xtb_result.get("symbols") or not xtb_result.get("positions"):
        return None

    from engines.render_engine import render_structure_png

    out_path = os.path.join(image_dir, f"{safe_filename(label)}.png")
    result = render_structure_png(
        xtb_result["symbols"], xtb_result["positions"], out_path, label=label,
    )
    if result["status"] != "ok":
        print(f"[report] structure image skipped for {label!r}: {result.get('error')}")
        return None
    return f"![{label}]({result['path']})"
