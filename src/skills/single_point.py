"""
single_point analysis skill: energy/geometry for ONE molecule (xtb GFN2-xTB
optimization -> PySCF DFT single-point). The original Phase 3/4 capability,
now the first entry in the analysis-skill registry (Phase 6 - architecture
plan Sec.16) instead of a bespoke ExecutorAgent method. Phase 7 (architecture
plan Sec.17) added an optional rendered structure image to describe()'s
report body.
"""
from __future__ import annotations

from skills._common import compute_species_from_conformer_search, render_species_image
from analysis.qc_common import classify, energy_present

NAME = "single_point"
DESCRIPTION = (
    "A plain energy/geometry question about ONE molecule - not a bond "
    "dissociation, redox, or any question needing more than one species."
)


def run(spec: dict) -> dict:
    result = compute_species_from_conformer_search(
        spec["smiles"],
        spec.get("charge", 0),
        spec.get("multiplicity", 1),
        spec["basis"],
        spec["xc"],
        refinement_engine=spec.get("refinement_engine", "pyscf"),
    )
    result["analysis_type"] = NAME
    return result


def review(results: dict) -> dict:
    """PASS/WARN/FAIL arithmetic now shared with analysis/qc_checks.py's
    CalculationRecord-based checks via analysis.qc_common.classify (Phase 8,
    partial - architecture plan Sec.19); only the notes this function builds
    (from a live engine results dict, not a parsed record) stay specific to
    this pipeline. 'refinement' replaces the old hard-coded 'pyscf' key so
    this works unchanged whichever refinement_engine ran (§19)."""
    engine_name = results.get("refinement_engine", "pyscf")
    notes = []

    if results.get("status") == "failed":
        notes.append(f"calculation_failed:{results.get('error')}")
    else:
        xtb = results.get("xtb") or {}
        refinement = results.get("refinement") or {}
        if not xtb.get("converged"):
            notes.append("xtb_not_converged")
        if not refinement.get("converged"):
            notes.append(f"{engine_name}_scf_not_converged")
        if not energy_present(refinement.get("energy_hartree")):
            notes.append("missing_refinement_energy")

    hard_fail_markers = (f"{engine_name}_scf_not_converged", "missing_refinement_energy", "calculation_failed")
    return classify(notes, hard_fail_markers)


def describe(
    spec: dict, results: dict, verdict: dict, image_dir: str = "report_structures"
) -> list[str]:
    engine_name = results.get("refinement_engine", "pyscf")
    lines = [
        "## Method tiers",
        "- xtb (GFN2-xTB) geometry optimization",
        f"- {engine_name} {spec['xc']}/{spec['basis']} single-point at the optimized geometry",
        "",
        "## Results",
    ]
    if results.get("status") == "failed":
        lines.append(f"- FAILED: {results.get('error')}")
    else:
        xtb = results.get("xtb") or {}
        refinement = results.get("refinement") or {}
        lines += [
            f"- xtb energy: {xtb.get('energy_hartree')} hartree (converged={xtb.get('converged')})",
            f"- {engine_name} energy: {refinement.get('energy_hartree')} hartree (converged={refinement.get('converged')})",
        ]
        image_line = render_species_image(results.get("smiles", "species"), xtb, image_dir)
        if image_line:
            lines += ["", "## Structure", image_line]
    return lines
