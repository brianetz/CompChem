"""
Species-level compute for the staged driver (target workflow Aim 2): instead
of every campaign task running its whole pipeline, the driver plans all
campaigns first, collects the species they need, drops duplicates and
anything already in the results database, and runs ONE task per remaining
species through this module.

    python -m agents.compute_species --requests requests.json --index 3

Each task goes through the same cached path the pipeline uses
(skills._common.compute_species_*), so its result lands in $COMPCHEM_DB and
the campaigns that need it later find it as a cache hit.
"""
from __future__ import annotations

import argparse
import json
import sys

from skills import get_skill
from storage.results_db import ResultsDB, method_string, species_key

_PROTOCOL = {"smiles": "smiles-embed", "conformer": "conformer-search"}


def request_key(req: dict) -> str:
    method = method_string(req["refinement_engine"], req["xc"], req["basis"], _PROTOCOL[req["path"]])
    return species_key(req["smiles"], req["charge"], req["multiplicity"], method)[0]


def collect_requests(specs: list[dict]) -> tuple[list[dict], int, list[int]]:
    """Unique species requests across specs (keyed exactly as the cache keys
    them), each annotated with 'key' and 'needed_by' (how many campaigns
    asked for it). Also returns the total before deduplication, and the
    indices of specs whose species couldn't be listed (malformed spec, or a
    skill without species_requests) - those campaigns still run in the
    assemble stage, computing their own species the ordinary way."""
    unique: dict[str, dict] = {}
    total, skipped = 0, []
    for idx, spec in enumerate(specs):
        try:
            skill = get_skill(spec["analysis_type"])
            reqs = [{**r, "key": request_key(r)} for r in skill.species_requests(spec)]
        except Exception as exc:
            print(f"[compute_species] spec {idx} not deduplicated: {type(exc).__name__}: {exc}")
            skipped.append(idx)
            continue
        for req in reqs:
            total += 1
            if req["key"] in unique:
                unique[req["key"]]["needed_by"] += 1
            else:
                unique[req["key"]] = {**req, "needed_by": 1}
    return list(unique.values()), total, skipped


def split_cached(requests: list[dict], db: ResultsDB) -> tuple[list[dict], list[dict]]:
    """(already done in the database, still to compute)."""
    hits, missing = [], []
    for req in requests:
        (hits if db.lookup(req["key"]) else missing).append(req)
    return hits, missing


def compute(req: dict) -> dict:
    from skills._common import compute_species_from_smiles, compute_species_from_conformer_search

    fn = compute_species_from_smiles if req["path"] == "smiles" else compute_species_from_conformer_search
    return fn(req["smiles"], req["charge"], req["multiplicity"], req["basis"], req["xc"],
              refinement_engine=req["refinement_engine"])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--requests", required=True, help="JSON list written by the driver")
    p.add_argument("--index", type=int, required=True, help="Which request in the list to compute")
    args = p.parse_args()

    from storage.flowcept_store import provenance_session

    with open(args.requests) as f:
        req = json.load(f)[args.index]
    with provenance_session():
        result = compute(req)
    summary = {"key": req["key"], "smiles": req["smiles"], "status": result.get("status"),
               "error": result.get("error"), "cache": result.get("cache"),
               "energy_hartree": (result.get("refinement") or {}).get("energy_hartree")}
    print("SPECIES_RESULT " + json.dumps(summary))
    sys.exit(0 if result.get("status") == "ok" else 1)


if __name__ == "__main__":
    main()
