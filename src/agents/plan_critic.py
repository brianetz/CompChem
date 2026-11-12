"""
Opt-in pre-run critique of an LLM-generated campaign spec (BDE only).

Glue between the driver (agents/run_agentic_campaign.py, ``--critique-plan``)
and the deterministic validator (analysis/bde_plan_validator.py). Nothing in
the default pipeline imports this module: it only runs when the flag is given.

``critique_spec`` runs ``critique_plan`` on the spec and prints the findings
(so they appear right above the human plan-approval gate). If any check FAILs
it makes at most ``max_retries`` (default 1) corrective round-trip(s) through
``PlannerAgent.replan_campaign``, feeding the validator's findings back as the
feedback text. It never loops, and never decides on its own whether to run
anything - the human gate (or the driver's --auto-approve policy) still does.

Specs for analysis types the validator does not cover (e.g. single_point)
are returned untouched with ``report=None``.
"""
from __future__ import annotations

from analysis import bde_plan_validator as v
from storage.agent_provenance import agent_tool_task

_AGENT_ID = "PlanCritic"


def _print_report(report: v.PlanReport, label: str) -> None:
    print(f"\n--- Plan critique ({label}): {report.verdict} ---")
    for c in report.checks:
        if c.status != v.PASS:
            print(f"  [{c.status}] {c.name}: {c.message}")
    if report.verdict == v.PASS:
        print("  all checks passed")


async def critique_spec(
    spec: dict, planner, *, model: str | None = None, max_retries: int = 1
) -> tuple[dict, v.PlanReport | None]:
    """Return ``(spec, report)``; ``spec`` may be a corrected replacement."""
    if spec.get("analysis_type") != "bde":
        return spec, None

    report = v.critique_plan(spec)
    _print_report(report, "initial")

    attempts = 0
    while report.verdict == v.FAIL and attempts < max_retries:
        attempts += 1
        with agent_tool_task(
            _AGENT_ID, "corrective_replan",
            used={"attempt": attempts, "qc": report.to_qc()},
        ) as t:
            spec = await planner.replan_campaign(
                previous_spec=spec,
                previous_results={"note": "Not run. The plan was rejected by the deterministic pre-run validator."},
                verdict=report.to_qc(),
                feedback=(
                    "Fix every FAIL below (and address WARNs where sensible). Keep the same "
                    "scientific question and the same parent molecule.\n" + report.feedback()
                ),
                model=model,
            )
            print("Planner spec (corrected after plan critique):", spec)
            report = v.critique_plan(spec)
            t.end(generated={"qc": report.to_qc()})
        _print_report(report, f"after correction {attempts}")

    return spec, report


_RANK = {v.PASS: 0, v.WARN: 1, v.FAIL: 2}


def augment_verdict(results: dict, verdict: dict) -> dict:
    """Post-run: add BDE sanity-band and optimized-geometry connectivity checks
    to the critic's verdict. Returns a NEW dict; the skill's own verdict is only
    ever made stricter, never relaxed. Non-BDE results are returned unchanged.

    Spin checks (``check_spin``) use the <S^2> / per-atom spin populations the
    PySCF engine returns for open-shell species; engines that do not report
    them (e.g. NWChem today) are simply skipped.
    """
    if results.get("analysis_type") != "bde":
        return verdict

    checks = []
    if results.get("bde_kcal_mol") is not None:
        checks.append(v.check_bde_value(results["bde_kcal_mol"]))

    species = [("parent", results.get("parent") or {})]
    species += [(f.get("label") or f.get("smiles", "fragment"), f) for f in results.get("fragments", [])]

    for label, sp in species:
        refinement = sp.get("refinement") or {}
        if sp.get("status") == "failed" or refinement.get("s2") is None:
            continue
        try:
            expected = v.expected_radical_atoms(sp["smiles"]) if v.HAVE_RDKIT else None
            checks += v.check_spin(
                str(label), int(sp.get("multiplicity", 1)), refinement["s2"],
                atom_spin=refinement.get("atom_spin"), expected_atoms=expected or None,
            )
        except Exception as exc:
            checks.append(v.Check(f"{label}_spin", v.WARN, f"check could not run: {exc}"))

    if v.HAVE_RDKIT:
        for label, sp in species:
            xtb = sp.get("xtb") or {}
            if sp.get("status") == "failed" or not xtb.get("symbols") or not xtb.get("positions"):
                continue
            try:
                checks.append(v.check_connectivity(str(label), xtb["symbols"], xtb["positions"], sp["smiles"]))
            except Exception as exc:  # a diagnostic must never crash the campaign
                checks.append(v.Check(f"{label}_connectivity", v.WARN, f"check could not run: {exc}"))

    notes = [f"plan_critic:{c.name}:{c.message}".replace(";", ",") for c in checks if c.status != v.PASS]
    worst = max((c.status for c in checks), key=_RANK.get, default=v.PASS)
    flag = max(verdict.get("qc_flag", v.PASS), worst, key=_RANK.get)
    old = verdict.get("qc_notes")
    merged = ";".join(([old] if old else []) + notes) or None
    print("Post-run critique:", {"qc_flag": flag, "added_notes": notes})
    return {**verdict, "qc_flag": flag, "qc_notes": merged}
