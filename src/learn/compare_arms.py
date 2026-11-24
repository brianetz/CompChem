"""Fair two-arm comparison: summarise each repeat of each arm, then mean / min / max over repeats.

A repeat is the three eval parts (65 tasks) of one arm; each arm is given as NAME=r1a,r1b,r1c;r2a,...
End-to-end runs (plan -> DFT -> assemble) are summarised per task with --e2e NAME=run_dir.

usage: python -m learn.compare_arms --tasks tasks.json --out compare.json \
           A="runA1p1,runA1p2,runA1p3;runA2p1,..." B="..." [--e2e A=run --e2e B=run]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import statistics as st

from learn.eval_report import summarise

SCALARS = ["first_try_correct", "after_retry_correct", "wrong_plan_caught_by_oracle", "no_plan",
           "format_failure_rate", "content_mistake_rate", "attempts_per_campaign", "plan_wall_min"]


def repeat_stats(runs: list[str], tasks: dict) -> dict:
    s = summarise(runs, tasks)
    s["no_plan"] = s["failure_classes"].get("no_plan", 0)
    wall = [json.load(open(os.path.join(r, "run_summary.json")))["stages"]["plan"]["wall_s"] for r in runs]
    s["plan_wall_min"] = round(sum(wall) / 60, 1)   # summed over the three parts
    s["runs"] = [os.path.basename(r) for r in runs]
    return s


def agg(reps: list[dict]) -> dict:
    out = {}
    for k in SCALARS:
        v = [r[k] for r in reps]
        out[k] = {"mean": round(st.mean(v), 3), "min": min(v), "max": max(v), "values": v}
    mods = sorted({m for r in reps for m in r["first_attempt_accept_by_module"]})
    out["first_attempt_accept_by_module"] = {
        m: [sum(r["first_attempt_accept_by_module"].get(m, [0, 0])[0] for r in reps),
            sum(r["first_attempt_accept_by_module"].get(m, [0, 0])[1] for r in reps)] for m in mods}
    classes = sorted({c for r in reps for c in r["by_bond_class"]})
    out["by_bond_class"] = {c: [sum(r["by_bond_class"].get(c, [0, 0])[0] for r in reps),
                                sum(r["by_bond_class"].get(c, [0, 0])[1] for r in reps)] for c in classes}
    fc = {}
    for r in reps:
        for c, n in r["failure_classes"].items():
            fc[c] = fc.get(c, 0) + n
    out["failure_classes_total"] = fc
    # tasks solved first try in every repeat / in none
    sets = [set(r["tasks_correct"]) for r in reps]
    out["always_correct"] = len(set.intersection(*sets)) if sets else 0
    out["solved_by_any"] = len(set.union(*sets)) if sets else 0
    out["repeats"] = reps
    return out


_BDE = re.compile(r"(-?\d+(?:\.\d+)?)\s*kcal/mol")


def e2e(run: str, tasks: dict) -> dict:
    o = json.load(open(os.path.join(run, "oracle.json")))
    summ = json.load(open(os.path.join(run, "run_summary.json")))
    comp = summ["stages"].get("compute", {})
    rows = []
    for k, v in sorted(o.items(), key=lambda x: int(x[0])):
        d = os.path.join(run, f"task_{k}")
        rep = os.path.join(d, "agentic_campaign_report.md")
        txt = open(rep).read() if os.path.exists(rep) else ""
        m = re.search(r"QC verdict: \*\*(\w+)\*\*", txt)
        verdict = m.group(1) if m else None
        sec = txt.split("## Bond dissociation energy", 1)[1] if "## Bond dissociation energy" in txt else ""
        bde = _BDE.findall(sec)
        converged = txt.count("converged=True"), txt.count("converged=False")
        rows.append({"task": v["task"], "plan_correct": v["first"]["correct"],
                     "error_class": v["first"]["error_class"],
                     "scf_converged": converged[0], "scf_not_converged": converged[1],
                     "report": bool(txt), "status": summ.get("task_status", {}).get(k),
                     "results_verdict": verdict, "bde_kcal": float(bde[0]) if bde else None})
    return {"run": os.path.basename(run), "tasks": len(rows),
            "plan_correct": sum(r["plan_correct"] for r in rows),
            "reports": sum(r["report"] for r in rows),
            "completed_correct": sum(r["plan_correct"] and r["report"] and r["bde_kcal"] is not None for r in rows),
            "species_tasks": comp.get("tasks"), "species_ok": comp.get("tasks_ok"), "db_hits": comp.get("db_hits"),
            "wall_min": {s: round(x.get("wall_s", 0) / 60, 1) for s, x in summ["stages"].items() if isinstance(x, dict)},
            "rows": rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--e2e", action="append", default=[])
    ap.add_argument("arms", nargs="+")
    a = ap.parse_args()
    tasks = {t["id"]: t for t in json.load(open(a.tasks))["tasks"]}
    out = {"arms": {}, "e2e": {}}
    for spec in a.arms:
        name, reps = spec.split("=", 1)
        out["arms"][name] = agg([repeat_stats(r.split(","), tasks) for r in reps.split(";")])
        s = out["arms"][name]
        print(name, {k: s[k]["values"] for k in SCALARS})
    for spec in a.e2e:
        name, run = spec.split("=", 1)
        out["e2e"][name] = e = e2e(run, tasks)
        print(name, "e2e plan", e["plan_correct"], "/", e["tasks"], "reports", e["reports"], "completed+correct", e["completed_correct"])
    json.dump(out, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
