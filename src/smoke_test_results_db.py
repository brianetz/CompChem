"""
Results-database smoke test (storage/results_db.py): canonical keys, bond
typing, cross-process claim deduplication, and a real water BDE computed
twice - the second time entirely from cache.

Run from src/ (or with PYTHONPATH=src), on a filesystem that supports POSIX
locks (Orion does; it's mounted with flock):
    python smoke_test_results_db.py [db_dir]
"""
from __future__ import annotations

import multiprocessing as mp
import os
import sys
import tempfile
import time

from storage.flowcept_store import provenance_session
from storage.results_db import (
    ResultsDB, canonical_smiles, detect_broken_bond, get_or_compute_species, method_string, species_key,
)


def check_keys():
    m = method_string("pyscf", "B3LYP", "6-31G*", "smiles-embed")
    assert species_key("O", 0, 1, m) == species_key("[H]O[H]", 0, 1, m), "same species, different spelling"
    assert species_key("O", 0, 1, m) != species_key("O", 0, 3, m), "multiplicity must change the key"
    assert species_key("O", 0, 1, m) != species_key("O", 0, 1, method_string("pyscf", "b3lyp", "def2-svp", "smiles-embed"))
    assert species_key("O", 0, 1, m) == species_key("O", 0, 1, method_string("pyscf", "b3lyp", "6-31G(d)", "smiles-embed"))
    assert method_string("pyscf", "b3lyp", "def2SVP", "x") == method_string("pyscf", "B3LYP", "def2-svp", "x")
    assert canonical_smiles("not a smiles((").startswith("raw:")
    print("OK keys")


def check_bonds():
    cases = [
        ("O", ["[OH]", "[H]"], "O-H"),
        ("CO", ["[CH2]O", "[H]"], "C-H"),
        ("CO", ["C[O]", "[H]"], "O-H"),
        ("CO", ["[CH3]", "[OH]"], "O-C"),
        ("CC", ["[CH3]", "[CH3]"], "C-C"),
        ("CCc1ccccc1", ["C[CH]c1ccccc1", "[H]"], "C-H"),
        ("O", ["[OH-]", "[H+]"], "O-H"),          # heterolysis: connectivity fallback
        ("O", ["[O]", "[H]"], None),              # not a single-bond break
        ("CCO", ["CCO", "[H]"], None),            # malformed planner output: no X-H cut fits
        ("CO", ["CO", "[H]"], None),
        ("C=O", ["[CH2]", "[O]"], "O=C"),         # double bond, connectivity fallback
    ]
    for parent, frags, want in cases:
        got = detect_broken_bond(parent, frags)
        assert got == want, f"{parent} -> {frags}: expected {want}, got {got}"
    print("OK bond typing")


def check_claims(tmp):
    db = ResultsDB(os.path.join(tmp, "claims.sqlite"), claim_ttl_s=2)
    state, tok_a = db.claim("k1", "C", 0, 1, "m")
    assert state == "claimed"
    assert db.claim("k1", "C", 0, 1, "m")[0] == "running", "live claim must be respected"
    time.sleep(2.2)  # A stops heartbeating -> stale
    state, tok_b = db.claim("k1", "C", 0, 1, "m")
    assert state == "claimed" and tok_b != tok_a, "stale claim must be taken over"
    assert not db.store_failure("k1", tok_a, "late failure from A", 1.0), "superseded claimant must not write"
    assert not db.heartbeat("k1", tok_a)
    assert db.store_result("k1", tok_b, {"refinement": {"energy_hartree": -1.0, "converged": True}}, 1.0)
    assert db.get_species("k1")["status"] == "done"
    os.environ["COMPCHEM_RUN_ID"] = "run_test"
    db.claim("k2", "CC", 0, 1, "m")
    assert db.release_run_claims("run_test", "test") == 1 and db.get_species("k2")["status"] == "failed"
    del os.environ["COMPCHEM_RUN_ID"]
    print("OK claims: takeover, token guard, heartbeat, release")


def check_long_compute_heartbeat(tmp):
    """A compute longer than the TTL stays claimed (heartbeats) - waiters don't steal it."""
    db = ResultsDB(os.path.join(tmp, "hb.sqlite"), claim_ttl_s=1.2)
    calls = []

    def slow():
        calls.append(1)
        time.sleep(3)
        return {"status": "ok", "refinement": {"energy_hartree": -2.0, "converged": True}}

    import threading
    out = {}
    t = threading.Thread(target=lambda: out.setdefault("a", get_or_compute_species(
        smiles="C", charge=0, multiplicity=1, engine="pyscf", xc="b3lyp", basis="hb",
        protocol="smiles-embed", compute=slow, db=db)))
    t.start(); time.sleep(0.5)
    b = get_or_compute_species(smiles="C", charge=0, multiplicity=1, engine="pyscf", xc="b3lyp",
                               basis="hb", protocol="smiles-embed", compute=slow, db=db)
    t.join()
    assert len(calls) == 1 and b["cache"]["outcome"] == "waited", (calls, b["cache"])
    print("OK heartbeat: 3 s compute with 1.2 s TTL computed once")


def check_rag_and_fallback(tmp):
    from storage.rag import bond_types_in
    assert bond_types_in("the C=O bond and the O–H bond and C#N") == ["O=C", "O-H", "N#C"], bond_types_in("the C=O bond and the O–H bond and C#N")
    blocker = os.path.join(tmp, "a_file")
    open(blocker, "w").close()
    os.environ["COMPCHEM_DB"] = os.path.join(blocker, "results.sqlite")  # parent is a file: unusable
    from skills._common import _with_cache
    r = _with_cache("O", 0, 1, "sto-3g", "b3lyp", "pyscf", "smiles-embed",
                    lambda: {"status": "ok", "refinement": {"energy_hartree": -1.0, "converged": True}})
    assert r["cache"]["outcome"] in ("bypassed", "computed"), r["cache"]
    del os.environ["COMPCHEM_DB"]
    print(f"OK rag bond regex; DB failure falls back ({r['cache']['outcome']})")


def _racer(db_path, counter_path, out_q):
    db = ResultsDB(db_path)

    def compute():
        with open(counter_path, "a") as f:
            f.write(f"{os.getpid()}\n")
        time.sleep(3)
        return {"smiles": "[OH]", "charge": 0, "multiplicity": 2, "status": "ok", "error": None,
                "refinement_engine": "pyscf", "xtb": {"converged": True},
                "refinement": {"energy_hartree": -75.7, "converged": True}}

    r = get_or_compute_species(smiles="[OH]", charge=0, multiplicity=2, engine="pyscf", xc="b3lyp",
                               basis="race", protocol="smiles-embed", compute=compute, db=db)
    out_q.put(r["cache"]["outcome"])


def check_race(tmp):
    db_path, counter = os.path.join(tmp, "race.sqlite"), os.path.join(tmp, "computed.txt")
    ResultsDB(db_path)
    q = mp.Queue()
    procs = [mp.Process(target=_racer, args=(db_path, counter, q)) for _ in range(8)]
    [p.start() for p in procs]
    [p.join() for p in procs]
    outcomes = sorted(q.get() for _ in procs)
    n_computed = len(open(counter).read().split())
    assert n_computed == 1, f"expected exactly one compute, got {n_computed}"
    assert outcomes.count("computed") == 1, outcomes
    print(f"OK race: 8 processes -> 1 compute; outcomes {outcomes}")


def check_water_bde(tmp):
    os.environ["COMPCHEM_DB"] = os.path.join(tmp, "results.sqlite")
    from skills import bde

    spec = {"analysis_type": "bde", "parent": {"smiles": "O", "charge": 0, "multiplicity": 1},
            "fragments": [{"label": "OH", "smiles": "[OH]", "charge": 0, "multiplicity": 2},
                          {"label": "H", "smiles": "[H]", "charge": 0, "multiplicity": 2}],
            "basis": "6-31g*", "xc": "b3lyp", "source_question": "smoke test: O-H BDE of water"}
    t0 = time.time(); first = bde.run(spec); t1 = time.time(); second = bde.run(spec); t2 = time.time()
    outcomes1 = [s["cache"]["outcome"] for s in [first["parent"], *first["fragments"]]]
    outcomes2 = [s["cache"]["outcome"] for s in [second["parent"], *second["fragments"]]]
    assert outcomes1 == ["computed"] * 3, outcomes1
    assert outcomes2 == ["hit"] * 3, outcomes2
    assert abs(first["bde_kcal_mol"] - second["bde_kcal_mol"]) < 1e-9
    assert first["bond_type"] == "O-H", first.get("bond_type")
    from storage.results_db import get_db
    rows = get_db().query_bdes(bond_types=["O-H"])
    assert len(rows) == 1 and rows[0]["qc_flag"] == "PASS", rows
    print(f"OK water BDE {first['bde_kcal_mol']:.2f} kcal/mol: cold {t1 - t0:.1f}s, warm {t2 - t1:.2f}s; "
          f"stats {get_db().stats()}")


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else None
    with tempfile.TemporaryDirectory(dir=base) as tmp:
        check_keys()
        check_bonds()
        check_race(tmp)
        check_claims(tmp)
        check_long_compute_heartbeat(tmp)
        with provenance_session():
            check_rag_and_fallback(tmp)
            check_water_bde(tmp)
    print("\nOK: results database behaves as expected.")


if __name__ == "__main__":
    main()
