"""
Phase 2 smoke test: prove run_individual_bde (old/Parsl pipeline) no longer
returns misleading -1 sentinels on failure - it should return None values
plus a clear status/error field instead.

Run:
    python src/smoke_test_bde_hardening.py
"""
import os
import tempfile

from bde_workflow import run_individual_bde
from storage.flowcept_store import provenance_session

FAKE_OUT = """Total DFT energy  =  {energy:.6f}
Zero-Point correction to Energy  =    25.0000 kcal/mol
Thermal correction to Enthalpy  =    27.0000 kcal/mol
Total Entropy  =    70.0000 cal/mol-K
Task  times  cpu:        0.1s     wall:        0.2s
"""


def main():
    with tempfile.TemporaryDirectory() as outdir, provenance_session():
        # Good case: both fragment .out files present and parseable
        with open(os.path.join(outdir, "good_1.out"), "w") as f:
            f.write(FAKE_OUT.format(energy=-20.1))
        with open(os.path.join(outdir, "good_2.out"), "w") as f:
            f.write(FAKE_OUT.format(energy=-20.2))

        good_result = run_individual_bde(
            e0=-40.5, frags={"label": "good"}, h0=0.03, outdir=outdir, s0=0.1, z0=0.03,
        )
        print("good case:", good_result)
        assert good_result["status"] == "ok"
        assert good_result["bd_energy"] is not None

        # Bad case: only ONE of the two fragment .out files exists (the other
        # job "failed" and never produced output) - this used to silently
        # return {bd_energy: -1, ...} with no indication anything went wrong.
        with open(os.path.join(outdir, "bad_1.out"), "w") as f:
            f.write(FAKE_OUT.format(energy=-20.1))
        # bad_2.out intentionally not written

        bad_result = run_individual_bde(
            e0=-40.5, frags={"label": "bad"}, h0=0.03, outdir=outdir, s0=0.1, z0=0.03,
        )
        print("bad case:", bad_result)
        assert bad_result["status"] == "failed"
        assert bad_result["bd_energy"] is None
        assert bad_result["error"]

    print("\nOK: failures now come back as None + status/error, not a misleading -1.")


if __name__ == "__main__":
    main()
