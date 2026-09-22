"""
Shared PASS/WARN/FAIL classification logic (Phase 8, partial - architecture
plan Sec.19). Before this module existed, the same "given a list of QC
problem notes, decide PASS vs WARN vs FAIL" block was written out three
times independently: once in analysis/qc_checks.py (the old
CalculationRecord/cclib-based pipeline), and once each in
skills/single_point.py and skills/bde.py (the Phase 6 skill-registry
pipeline). All three used the same shape (a notes list, a set of substrings
that force FAIL, everything else WARN) but had drifted into three separate
copies that could disagree given the same underlying facts. This module is
the one place that decision now lives; every caller only supplies its own
domain-specific notes and its own definition of which notes are hard
failures for its situation.

Deliberately NOT unifying the *notes themselves* (qc_checks.py's
"job_not_successful"/"missing_geometry" vocabulary vs. the skills'
"xtb_not_converged"/"<engine>_scf_not_converged" vocabulary) - those come
from genuinely different underlying data (a parsed CalculationRecord vs. a
live engine-call results dict) and forcing one shared vocabulary would
obscure which pipeline actually produced which note. What's shared is only
the classification arithmetic itself.
"""
from __future__ import annotations

from typing import Iterable, Optional


def classify(notes: Iterable[str], hard_fail_markers: Iterable[str]) -> dict:
    """Turn a list of problem notes into a {'qc_flag', 'qc_notes'} verdict.

    No notes -> PASS. Any note containing one of `hard_fail_markers` as a
    substring -> FAIL. Any other non-empty notes -> WARN. `qc_notes` is the
    semicolon-joined notes list, or None if there were none - matching the
    contract every caller (old and new pipeline alike) already expects.
    """
    notes = list(notes)
    hard_fail_markers = tuple(hard_fail_markers)

    if not notes:
        qc_flag = "PASS"
    elif any(marker in note for note in notes for marker in hard_fail_markers):
        qc_flag = "FAIL"
    else:
        qc_flag = "WARN"

    return {"qc_flag": qc_flag, "qc_notes": ";".join(notes) if notes else None}


def energy_present(energy_hartree: Optional[float]) -> bool:
    """Tiny named predicate so call sites read as intent
    ('not energy_present(x)') rather than a bare None-check repeated
    everywhere - purely a readability aid, no behavior of its own."""
    return energy_hartree is not None
