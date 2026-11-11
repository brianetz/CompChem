"""Score board-v2 runs: per campaign, correct / wrong answer / failure reason.

usage: classify_board.py <run_dir> [...]
"""
import glob
import json
import os
import re
import sys
from collections import Counter

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def can(s):
    m = Chem.MolFromSmiles(s or "")
    return Chem.MolToSmiles(m) if m else f"?{s}"


# keyword in the question -> (parent, fragments) of the correct homolysis
EXPECTED = {
    "toluene": ("Cc1ccccc1", ["[CH2]c1ccccc1", "[H]"]),
    "methanol": ("CO", ["C[O]", "[H]"]),
    "ammonia": ("N", ["[NH2]", "[H]"]),
    "methanethiol": ("CS", ["C[S]", "[H]"]),
    "chloromethane": ("CCl", ["[CH3]", "[Cl]"]),
    "hydrogen peroxide": ("OO", ["[OH]", "[OH]"]),
    "ethane": ("CC", ["[CH3]", "[CH3]"]),
    "phenol": ("Oc1ccccc1", ["[O]c1ccccc1", "[H]"]),
    "propene": ("C=CC", ["[CH2]C=C", "[H]"]),
    "hydrogen fluoride": ("F", ["[F]", "[H]"]),
    "methane": ("C", ["[CH3]", "[H]"]),
    "methylamine": ("CN", ["C[NH]", "[H]"]),
    "formic acid": ("O=CO", ["O=C[O]", "[H]"]),
    "bromomethane": ("CBr", ["[CH3]", "[Br]"]),
    "propane": ("CCC", ["[CH3]", "[CH2]C"]),
    "dimethyl disulfide": ("CSSC", ["C[S]", "C[S]"]),
    "acetaldehyde": ("CC=O", ["C[C]=O", "[H]"]),
    "hydrogen chloride": ("Cl", ["[Cl]", "[H]"]),
}


def expected_for(question):
    q = question.lower()
    for key in sorted(EXPECTED, key=len, reverse=True):
        if key in q:
            p, f = EXPECTED[key]
            return key, can(p), sorted(can(x) for x in f)
    return None, None, None


def classify(task_dir, question):
    key, ep, ef = expected_for(question)
    out = {"molecule": key}
    spec_p = os.path.join(task_dir, "spec.json")
    plan_log = os.path.join(os.path.dirname(task_dir), f"plan_{task_dir.rsplit('_', 1)[1]}.log")
    if not os.path.exists(spec_p):
        txt = open(plan_log, errors="ignore").read() if os.path.exists(plan_log) else ""
        if "planning time budget used up" in txt:
            why = "planning budget ran out"
        elif "Plan board" not in txt and "Spec written" not in txt:
            why = "worker lost"  # killed (block walltime) before the planner printed anything
        elif "ResponseError" in txt or "status code: 500" in txt:
            why = "LLM output format (Ollama 500)"
        elif "ManagerLost" in txt or "walltime" in txt.lower():
            why = "worker lost"
        else:
            why = "no plan produced"
        return {**out, "outcome": "no plan", "why": why}
    spec = json.load(open(spec_p))
    rep_p = os.path.join(task_dir, "agentic_campaign_report.md")
    report = open(rep_p).read() if os.path.exists(rep_p) else ""
    verdict = (re.search(r"QC verdict: \*\*(\w+)\*\*", report) or [None, "NO REPORT"])[1]
    bde = re.search(r"- ([-\d.]+) kcal/mol", report)
    out.update(verdict=verdict, bde=float(bde.group(1)) if bde else None,
               method=f"{spec.get('xc')}/{spec.get('basis')}".lower())
    parent = can((spec.get("parent") or {}).get("smiles"))
    frags = sorted(can(f.get("smiles")) for f in spec.get("fragments") or [])
    right = parent == ep and frags == ef
    out["fragments"] = " + ".join(frags)
    if verdict == "PASS":
        out["outcome"] = "correct" if right else "wrong answer passed"
    else:
        notes = report.split("QC verdict:", 1)[1].splitlines()[0] if "QC verdict:" in report else ""
        out["outcome"] = "failed"
        out["why"] = ("invalid SMILES" if "parse SMILES" in notes else
                      "impossible spin" if "multiplicity" in notes.lower() or "spin" in notes.lower() else notes[:90])
    return out


def main():
    for run in sys.argv[1:]:
        s = json.load(open(f"{run}/run_summary.json"))
        tq = s.get("task_questions") or {}
        rows = []
        for t in sorted(glob.glob(f"{run}/task_*/"), key=lambda p: int(p.rstrip("/").rsplit("_", 1)[1])):
            i = t.rstrip("/").rsplit("_", 1)[1]
            q = tq.get(i) or tq.get(int(i)) if isinstance(tq, dict) else None
            q = q or s.get("question") or ""
            rows.append({"task": i, **classify(t.rstrip("/"), q)})
        c = Counter(r["outcome"] for r in rows)
        print(f"== {os.path.basename(run)}  {dict(c)}  plan wall {s['stages'].get('plan', {}).get('wall_s')} s, total {s.get('total_wall_s')} s")
        for r in rows:
            if r["outcome"] != "correct":
                print("   ", r)
        print("    methods:", dict(Counter(r.get("method") for r in rows if r.get("method"))))
        json.dump(rows, open(f"{run}/.board_scores.json", "w"), indent=1)


if __name__ == "__main__":
    main()
