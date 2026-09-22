#!/usr/bin/env python3
"""
Stand-in for the real NWChem binary, used ONLY for local Parsl-orchestration
testing (Phase 1) when the real binary/Slurm allocation isn't available (e.g.
on a laptop). Reads a .nw input file, fabricates a plausible-but-fake energy
(deterministic per-input, so reruns are reproducible), and prints
NWChem-output-shaped text to stdout so the existing parsers (parser.py's
line-based parse_nwchem_output, and nwchem_fallback.py's regexes) work
end-to-end without modification.

This does NOT do any real chemistry - it exists purely to prove the Parsl
orchestration layer (concurrent submission, futures, failure handling)
actually works, decoupled from whether a real QC engine is available.

Usage:   python mock_nwchem.py <input.nw>
To use:  export NWCHEM_CMD="python $(pwd)/mock_nwchem.py"   (see run_all.py)
"""
import hashlib
import random
import sys
import time


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: mock_nwchem.py <input.nw>", file=sys.stderr)
        return 1

    nw_path = sys.argv[1]
    with open(nw_path) as f:
        nw_text = f.read()

    # Deterministic-but-varied fake numbers, seeded off the input file's
    # content so the same molecule always gets the same fake energy.
    seed = int(hashlib.sha256(nw_text.encode()).hexdigest(), 16) % (2**32)
    rng = random.Random(seed)

    energy_hartree = -40.0 - rng.random() * 20
    zpe_kcal = 20 + rng.random() * 20
    enthalpy_kcal = zpe_kcal + rng.random() * 5
    entropy_cal = 60 + rng.random() * 30

    time.sleep(rng.uniform(0.5, 2.0))  # simulate real jobs taking varying time

    print(f"Job input: {nw_path}")
    print()
    print(f"Total DFT energy  =  {energy_hartree:.6f}")
    print(f"Zero-Point correction to Energy  =    {zpe_kcal:.4f} kcal/mol")
    print(f"Thermal correction to Enthalpy  =    {enthalpy_kcal:.4f} kcal/mol")
    print(f"Total Entropy  =    {entropy_cal:.4f} cal/mol-K")
    print()
    print("Task  times  cpu:        0.1s     wall:        0.2s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
