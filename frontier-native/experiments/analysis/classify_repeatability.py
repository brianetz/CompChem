"""Classify every campaign of a methanol O-H repeatability run.

usage: classify_repeatability.py <run_dir> [<run_dir> ...]
"""
import glob
import json
import os
import re
import sys
from collections import Counter

from storage.results_db import canonical_smiles, detect_broken_bond, normalize_basis

# Independent references (login node, no database), B3LYP, kcal/mol
REF = {
    "6-31g*":    {"O-H": 103.99, "C-H": 101.75, "O-C": 96.39},
    "def2-svp":  {"O-H": 105.87, "C-H": 99.20,  "O-C": 97.76},
    "def2-tzvp": {"O-H": 107.49, "C-H": 100.18, "O-C": 93.50},
}
CORRECT_FRAGS = sorted([(canonical_smiles("C[O]"), 0, 2), (canonical_smiles("[H]"), 0, 2)])


def classify(task_dir):
    spec_path = os.path.join(task_dir, "spec.json")
    if not os.path.exists(spec_path):
        return {"outcome": "planning failed"}
    spec = json.load(open(spec_path))
    rep_path = os.path.join(task_dir, "agentic_campaign_report.md")
    report = open(rep_path).read() if os.path.exists(rep_path) else ""
    verdict = (re.search(r"QC verdict: \*\*(\w+)\*\*", report) or [None, "NO REPORT"])[1]
    basis = normalize_basis(spec.get("basis", "?"))
    xc = (spec.get("xc") or "?").lower()
    method = f"{xc}/{basis}"
    info = {"method": method, "verdict": verdict}
    if spec.get("analysis_type") != "bde" or not spec.get("fragments"):
        info.update(outcome="wrong result", why=f"analysis_type={spec.get('analysis_type')}")
        return info
    parent = spec["parent"]
    frags = spec["fragments"]
    try:
        bond = detect_broken_bond(parent["smiles"], [f["smiles"] for f in frags])
    except Exception:
        bond = None
    info["bond"] = bond
    if verdict != "PASS":
        m = re.search(r"\((.*)\)", report.split("QC verdict:")[1].splitlines()[0]) if "QC verdict:" in report else None
        notes = m.group(1) if m else verdict
        why = ("SMILES not parseable" if "parse SMILES" in notes else
               "spin/electron mismatch" if "spin" in notes.lower() or "electrons" in notes.lower() else notes[:80])
        info.update(outcome="failed", why=why)
        return info
    bde = float(re.search(r"- ([-\d.]+) kcal/mol", report).group(1))
    info["bde"] = bde
    got = sorted((canonical_smiles(f["smiles"]), int(f.get("charge", 0)), int(f.get("multiplicity", 1))) for f in frags)
    right_species = canonical_smiles(parent["smiles"]) == canonical_smiles("CO") and got == CORRECT_FRAGS
    ref = REF.get(basis, {}).get("O-H") if xc == "b3lyp" else None
    if not right_species:
        info.update(outcome="wrong result", why=f"broke {bond or 'unknown bond'}: {' + '.join(s for s, _, _ in got)}")
    elif ref is not None and abs(bde - ref) > 0.01:
        info.update(outcome="wrong result", why=f"{bde:.2f} != reference {ref:.2f}")
    else:
        info.update(outcome="correct", ref=ref)
    return info


def main():
    for run_dir in sys.argv[1:]:
        tasks = sorted([p for p in glob.glob(f"{run_dir}/task_*") if os.path.isdir(p)], key=lambda p: int(p.rsplit("_", 1)[1]))
        rows = [classify(t) for t in tasks]
        correct_methods = Counter(r["method"] for r in rows if r["outcome"] == "correct")
        modal = correct_methods.most_common(1)[0][0] if correct_methods else None
        for r in rows:
            if r["outcome"] == "correct":
                r["outcome"] = "same result" if r["method"] == modal else "correct, other method"
        summary = json.load(open(f"{run_dir}/run_summary.json"))
        comp = summary["stages"].get("compute", {})
        out = {
            "run": os.path.basename(run_dir), "n": len(rows), "modal_method": modal,
            "outcomes": dict(Counter(r["outcome"] for r in rows)),
            "pass": sum(r.get("verdict") == "PASS" for r in rows),
            "values": {f"{m} {b:.2f} [{o}]": n for (m, b, o), n in
                       Counter((r["method"], r["bde"], r["outcome"]) for r in rows if "bde" in r).most_common()},
            "wrong_why": dict(Counter(r["why"] for r in rows if r["outcome"] == "wrong result")),
            "failed_why": dict(Counter(r["why"] for r in rows if r["outcome"] == "failed")),
            "methods_all": dict(Counter(r.get("method") for r in rows if "method" in r)),
            "species_unique": comp.get("species_unique"), "species_tasks": comp.get("tasks"),
            "species_ok": comp.get("tasks_ok"), "db_hits": comp.get("db_hits"),
            "stages": {k: v.get("wall_s") for k, v in summary["stages"].items()},
            "total_wall_s": summary.get("total_wall_s"),
        }
        print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
