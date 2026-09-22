"""
Analysis skill registry (Phase 6 - architecture plan Sec.16).

An "analysis skill" is a self-contained module in this package exposing:
    NAME: str                                       - the analysis_type key
    DESCRIPTION: str                                 - one-line description,
        used to build PlannerAgent's LLM tool schema dynamically, so a new
        skill is automatically choosable by the planner once registered
        below - no hand-editing llm/planner_llm.py's prompt each time.
    run(spec: dict) -> dict                          - blocking; does the
        actual calculation(s), embeds 'analysis_type' in its return dict.
    review(results: dict) -> dict                    - {'qc_flag', 'qc_notes'}
    describe(spec, results, verdict,
              image_dir='report_structures') -> list[str]
        - the analysis-specific body of the markdown report (the header -
        title/timestamp/verdict/question/rationale - is generic and lives
        in ReporterAgent). Phase 7 (architecture plan Sec.17) added the
        optional image_dir parameter: a skill may call
        skills._common.render_species_image to embed a PyMOL-rendered PNG
        of a species' optimized structure under that directory - purely a
        report-layer nicety, never required, and never a reason to fail
        the report if PyMOL isn't installed or rendering fails.

ExecutorAgent, CriticAgent, and ReporterAgent each have exactly ONE generic
action that looks a skill up here by analysis_type and calls it. Adding a
new analysis (IR frequencies, NMR shifts, UV-Vis/TDDFT, HOMO-LUMO gap, redox
potentials, Marcus theory, ...) means writing one new module with this same
shape and adding it below - not touching any agent class.
"""
from __future__ import annotations

from skills import single_point, bde

_REGISTRY = {
    single_point.NAME: single_point,
    bde.NAME: bde,
}


def get_skill(analysis_type: str):
    try:
        return _REGISTRY[analysis_type]
    except KeyError:
        raise ValueError(
            f"Unknown analysis_type {analysis_type!r} - registered skills: "
            f"{sorted(_REGISTRY)}"
        ) from None


def list_skills() -> list[tuple[str, str]]:
    """(name, description) for every registered skill - used to build
    PlannerAgent's LLM tool schema (see llm/planner_llm.py) without
    hand-editing it every time a skill is added."""
    return [(mod.NAME, mod.DESCRIPTION) for mod in _REGISTRY.values()]
