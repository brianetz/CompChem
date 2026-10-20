"""
bde analysis skill: bond dissociation energy for a parent molecule plus one
or more dissociation fragments (BDE = sum(fragment energies) - parent
energy). Phase 5's original multi-species capability, now the second entry
in the analysis-skill registry (Phase 6 - architecture plan Sec.16). Phase 7
(architecture plan Sec.17) added rendered structure images for the parent
and every fragment to describe()'s report body.
"""
from __future__ import annotations

from skills._common import compute_species_from_smiles, render_species_image
from analysis.qc_common import classify, energy_present
from utils.units import delta_hartree_to_kcal
from storage.results_db import get_db, detect_broken_bond, canonical_smiles, method_string
from storage.rag import bond_types_in
from analysis.species_checks import atom_balance

NAME = "bde"
# Electronic homolytic BDEs of real bonds: weak peroxide/halogen bonds sit
# near 30-50 kcal/mol, the strongest (C#O, N#N) near 226-257. Outside this,
# the species do not describe one bond breaking.
_BDE_RANGE_KCAL = (0.0, 270.0)
DESCRIPTION = (
    "A bond dissociation energy question: needs a parent molecule plus "
    "one or more dissociation fragments - the species the SPECIFIC bond "
    "breaks into (NOT full atomization of the molecule)."
)


def species_requests(spec: dict) -> list[dict]:
    """Every species this spec needs computed, as keyword arguments for
    skills._common.compute_species_from_smiles - lets a driver deduplicate
    species across many specs before submitting any compute
    (agents/compute_species.py)."""
    basis, xc = spec["basis"], spec["xc"]
    engine = spec.get("refinement_engine", "pyscf")
    species = [spec["parent"], *(spec.get("fragments") or [])]
    return [
        {"path": "smiles", "smiles": s["smiles"], "charge": s.get("charge", 0),
         "multiplicity": s.get("multiplicity", 1), "basis": basis, "xc": xc,
         "refinement_engine": engine}
        for s in species
    ]


def _record_bde(spec: dict, results: dict) -> None:
    """Index a finished BDE by bond type in the results database (no-op
    when caching is off or any species lacks a cache key)."""
    db = get_db()
    parent = results["parent"]
    frags = results["fragments"]
    keys = [(s.get("cache") or {}).get("key") for s in [parent, *frags]]
    if db is None or not all(keys):
        return
    frag_smiles = [canonical_smiles(f["smiles"]) for f in frags]
    parent_smiles = canonical_smiles(parent["smiles"])
    results["bond_type"] = detect_broken_bond(parent_smiles, frag_smiles)
    results["bde_record_id"] = db.record_bde(
        parent_key=keys[0], fragment_keys=keys[1:], parent_smiles=parent_smiles,
        fragment_smiles=frag_smiles, bond_type=results["bond_type"],
        method=method_string(results["refinement_engine"], spec["xc"], spec["basis"], "smiles-embed"),
        bde_hartree=results["bde_hartree"], bde_kcal_mol=results["bde_kcal_mol"],
        qc_flag=review(results)["qc_flag"], source_question=spec.get("source_question"),
    )


def run(spec: dict) -> dict:
    parent_spec = spec["parent"]
    fragment_specs = spec.get("fragments") or []
    if not fragment_specs:
        raise ValueError("bde analysis requires at least one entry in spec['fragments']")

    basis, xc = spec["basis"], spec["xc"]
    refinement_engine = spec.get("refinement_engine", "pyscf")

    parent_result = compute_species_from_smiles(
        parent_spec["smiles"], parent_spec.get("charge", 0),
        parent_spec.get("multiplicity", 1), basis, xc,
        refinement_engine=refinement_engine,
    )

    fragment_results = []
    for frag in fragment_specs:
        frag_result = compute_species_from_smiles(
            frag["smiles"], frag.get("charge", 0), frag.get("multiplicity", 1), basis, xc,
            refinement_engine=refinement_engine,
        )
        frag_result["label"] = frag.get("label", frag["smiles"])
        fragment_results.append(frag_result)

    parent_refinement = parent_result.get("refinement")
    fragment_refinements = [f.get("refinement") for f in fragment_results]

    if parent_refinement is None or any(f is None for f in fragment_refinements):
        # At least one species failed (see compute_species_from_smiles) -
        # no BDE to report; review() below will flag this FAIL from each
        # species' 'status'/'error' fields.
        bde_hartree = None
        bde_kcal_mol = None
    else:
        parent_energy = parent_refinement["energy_hartree"]
        fragments_energy_sum = sum(f["energy_hartree"] for f in fragment_refinements)
        bde_hartree = fragments_energy_sum - parent_energy
        bde_kcal_mol = delta_hartree_to_kcal(bde_hartree)

    frag_smiles = [f["smiles"] for f in fragment_specs]
    results = {
        "analysis_type": NAME,
        "refinement_engine": refinement_engine,
        "parent": parent_result,
        "fragments": fragment_results,
        "bde_hartree": bde_hartree,
        "bde_kcal_mol": bde_kcal_mol,
        # For review(): what the fragments add up to, which bond they come
        # from, and which bond the question asked about.
        "atom_balance": atom_balance(parent_spec["smiles"], frag_smiles),
        "bond_type": detect_broken_bond(parent_spec["smiles"], frag_smiles),
        "requested_bonds": bond_types_in(spec.get("source_question") or ""),
    }
    try:
        _record_bde(spec, results)
    except Exception as exc:  # indexing is a side effect; never fail the BDE over it
        print(f"[results-db] BDE not indexed: {type(exc).__name__}: {exc}")
    return results


def review(results: dict) -> dict:
    """PASS/WARN/FAIL arithmetic now shared with analysis/qc_checks.py's
    CalculationRecord-based checks via analysis.qc_common.classify (Phase 8,
    partial - architecture plan Sec.19); the per-species notes stay specific
    to this pipeline. '_refinement_not_converged'/'_missing_energy' replace
    the old hard-coded '_pyscf_not_converged' naming so this is correct
    whichever refinement_engine ran (engine identity is still visible in
    the report body via describe(), just not baked into the note names)."""
    notes = []

    def _check(species: dict, label: str) -> None:
        if species.get("status") == "failed":
            notes.append(f"{label}_calculation_failed:{species.get('error')}")
            return
        # species.get("xtb", {}) would return None (not {}) when the key is
        # present but explicitly None - `or {}` normalizes both cases.
        xtb = species.get("xtb") or {}
        refinement = species.get("refinement") or {}
        if not xtb.get("converged"):
            notes.append(f"{label}_xtb_not_converged")
        if not refinement.get("converged"):
            notes.append(f"{label}_refinement_not_converged")
        if not energy_present(refinement.get("energy_hartree")):
            notes.append(f"{label}_missing_energy")

    _check(results.get("parent", {}), "parent")
    for frag in results.get("fragments", []):
        _check(frag, frag.get("label", frag.get("smiles", "fragment")))

    # Found 2026-10-19: a methanol "BDE" from [OH] + [H] (no CH3) came out
    # at 24 763 kcal/mol and passed, and a C-H break passed for an O-H
    # question. Convergence alone does not make the answer right.
    if results.get("atom_balance"):
        notes.append(f"atom_balance_mismatch:{results['atom_balance']}")
    bde = results.get("bde_kcal_mol")
    if bde is not None and not _BDE_RANGE_KCAL[0] <= bde <= _BDE_RANGE_KCAL[1]:
        notes.append(f"bde_implausible:{bde:.1f} kcal/mol outside {_BDE_RANGE_KCAL}")
    requested, broken = results.get("requested_bonds") or [], results.get("bond_type")
    if requested and broken and broken not in requested:
        notes.append(f"wrong_bond:{broken} broken, question asks {'/'.join(requested)}")
    elif requested and not broken and not results.get("atom_balance"):
        notes.append("bond_undetermined")

    hard_fail_markers = ("_refinement_not_converged", "_missing_energy", "_calculation_failed",
                         "atom_balance_mismatch", "bde_implausible", "wrong_bond")
    return classify(notes, hard_fail_markers)


def describe(
    spec: dict, results: dict, verdict: dict, image_dir: str = "report_structures"
) -> list[str]:
    parent = results.get("parent", {})
    fragments = results.get("fragments", [])
    engine_name = results.get("refinement_engine", "pyscf")

    def _species_line(prefix: str, species: dict) -> str:
        if species.get("status") == "failed":
            return f"- {prefix}: {species.get('smiles')} -> FAILED: {species.get('error')}"
        refinement = species.get("refinement") or {}
        outcome = (species.get("cache") or {}).get("outcome")
        cache_note = f" [cache: {outcome}]" if outcome and outcome != "disabled" else ""
        return (
            f"- {prefix}: {species.get('smiles')} "
            f"(charge={species.get('charge')}, multiplicity={species.get('multiplicity')}) -> "
            f"{refinement.get('energy_hartree')} hartree (converged={refinement.get('converged')})"
            f"{cache_note}"
        )

    lines = [
        "## Method tiers",
        "- xtb (GFN2-xTB) geometry optimization",
        f"- {engine_name} {spec['xc']}/{spec['basis']} single-point at each optimized geometry",
        "",
        "## Species",
        _species_line("Parent", parent),
    ]
    for frag in fragments:
        lines.append(_species_line(f"Fragment ({frag.get('label')})", frag))

    image_lines = []
    parent_image = render_species_image(
        f"parent_{parent.get('smiles', '?')}", parent.get("xtb"), image_dir
    )
    if parent_image:
        image_lines.append(parent_image)
    for frag in fragments:
        frag_image = render_species_image(
            frag.get("label", frag.get("smiles", "fragment")), frag.get("xtb"), image_dir
        )
        if frag_image:
            image_lines.append(frag_image)
    if image_lines:
        lines += ["", "## Structures"] + image_lines

    bde_kcal_mol = results.get("bde_kcal_mol")
    bde_line = (
        f"- {bde_kcal_mol:.2f} kcal/mol ({results.get('bde_hartree')} hartree)"
        if bde_kcal_mol is not None
        else "- BDE not computed"
    )
    lines += ["", "## Bond dissociation energy", bde_line]
    if results.get("bond_type"):
        lines.append(f"- bond: {results['bond_type']}")
    return lines
