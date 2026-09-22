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

NAME = "bde"
DESCRIPTION = (
    "A bond dissociation energy question: needs a parent molecule plus "
    "one or more dissociation fragments - the species the SPECIFIC bond "
    "breaks into (NOT full atomization of the molecule)."
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

    return {
        "analysis_type": NAME,
        "refinement_engine": refinement_engine,
        "parent": parent_result,
        "fragments": fragment_results,
        "bde_hartree": bde_hartree,
        "bde_kcal_mol": bde_kcal_mol,
    }


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

    hard_fail_markers = ("_refinement_not_converged", "_missing_energy", "_calculation_failed")
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
        return (
            f"- {prefix}: {species.get('smiles')} "
            f"(charge={species.get('charge')}, multiplicity={species.get('multiplicity')}) -> "
            f"{refinement.get('energy_hartree')} hartree (converged={refinement.get('converged')})"
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
    return lines
