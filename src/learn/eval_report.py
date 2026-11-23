"""Summarise evaluation runs (oracle.json + episodes) per condition.

usage: python -m learn.eval_report --tasks tasks.json --out report.json  NAME=run_dir,run_dir ... [NAME=...]
"""
from __future__ import annotations

import argparse
import json
import os
from collections import Counter, defaultdict

from learn.episodes import run_episodes


def summarise(runs: list[str], tasks: dict) -> dict:
    rows, eps = [], []
    for r in runs:
        o = json.load(open(os.path.join(r, "oracle.json")))
        rows += [dict(v, run=os.path.basename(r), index=int(k)) for k, v in o.items()]
        eps += [e for e in run_episodes(r, tasks) if e["phase"] == "first"]
    n = len(rows)
    first = sum(r["first"]["correct"] for r in rows)
    final = sum(r["first"]["correct"] or bool(r["retry"] and r["retry"]["correct"]) for r in rows)
    why = Counter(r["first"]["error_class"] for r in rows if not r["first"]["correct"])
    att = [e for e in eps if e["blame"] not in ("infrastructure", "checker_false_reject")]
    by_mod = defaultdict(lambda: [0, 0])
    for e in eps:
        if e.get("attempt") in (1, "1"):
            by_mod[e["module"]][1] += 1
            by_mod[e["module"]][0] += e["blame"] == "accepted"
    rounds = Counter((e["run"], e["task_index"]) for e in eps)
    per_class = defaultdict(lambda: [0, 0])
    for r in rows:
        t = tasks.get(r["task"], {})
        per_class[t.get("bond_class", "?")][1] += 1
        per_class[t.get("bond_class", "?")][0] += r["first"]["correct"]
    return {"campaigns": n, "first_try_correct": first, "after_retry_correct": final,
            "wrong_plan_caught_by_oracle": sum(1 for r in rows if r["first"]["error_class"] not in (None, "no_plan")),
            "failure_classes": dict(why),
            "attempts": len(att), "format_failure_rate": round(sum(e["blame"] == "planner_format" for e in att) / max(1, len(att)), 3),
            "content_mistake_rate": round(sum(e["blame"] == "planner_content" for e in att) / max(1, len(att)), 3),
            "first_attempt_accept_by_module": {m: [a, b] for m, (a, b) in sorted(by_mod.items())},
            "attempts_per_campaign": round(sum(rounds.values()) / max(1, len(rounds)), 1),
            "by_bond_class": {k: v for k, v in sorted(per_class.items(), key=lambda x: -x[1][1])},
            "tasks_correct": sorted(r["task"] for r in rows if r["first"]["correct"])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("conds", nargs="+")
    a = ap.parse_args()
    tasks = {t["id"]: t for t in json.load(open(a.tasks))["tasks"]}
    out = {}
    for c in a.conds:
        name, runs = c.split("=", 1)
        out[name] = summarise(runs.split(","), tasks)
        s = out[name]
        print(f"{name}: first {s['first_try_correct']}/{s['campaigns']}, retry {s['after_retry_correct']}, "
              f"format {s['format_failure_rate']}, content {s['content_mistake_rate']}, attempts/campaign "
              f"{s['attempts_per_campaign']}, fails {s['failure_classes']}, first-attempt accept {s['first_attempt_accept_by_module']}")
    json.dump(out, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
