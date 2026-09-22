from storage.records import CalculationRecord
from storage.flowcept_store import flowcept_task
from analysis.qc_common import classify, energy_present

# Notes that force FAIL rather than WARN for this pipeline specifically - a
# job that didn't run, or ran but produced no energy, is unusable; a
# stray imaginary frequency or a missing geometry is worth flagging but the
# record may still be salvageable (see analysis.qc_common's module
# docstring for why this stays pipeline-specific rather than shared).
_HARD_FAIL_MARKERS = ("job_not_successful", "missing_energy")


@flowcept_task
def apply_qc_checks(record: CalculationRecord) -> CalculationRecord:
    """Unchanged verdict logic from before Phase 8's partial unification
    (architecture plan Sec.19) - only the PASS/WARN/FAIL arithmetic itself
    was factored out into analysis.qc_common.classify, shared with the
    skill-registry pipeline's per-skill review() functions. The notes this
    function builds (job_not_successful/missing_energy/N_imaginary_freqs/
    missing_geometry) are unchanged."""
    notes = []

    if not record.success:
        notes.append("job_not_successful")

    if not energy_present(record.scf_energy_hartree):
        notes.append("missing_energy")

    if record.imaginary_freq_count is not None:
        if record.role != "ts" and record.imaginary_freq_count > 0:
            notes.append(f"{record.imaginary_freq_count}_imaginary_freqs")

    if record.n_atoms is None or record.n_atoms == 0:
        notes.append("missing_geometry")

    verdict = classify(notes, _HARD_FAIL_MARKERS)
    record.qc_flag = verdict["qc_flag"]
    record.qc_notes = verdict["qc_notes"]

    return record