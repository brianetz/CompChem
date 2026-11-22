"""Export plan-board provenance (Flowcept records written with COMPCHEM_TRACE=1)
into attempt-level episodes for lessons and training.

One episode = one Planner attempt at one module, with the prompt it saw, its
reply, how the board judged it, the oracle's grade of the whole campaign, and a
blame label. Only attempts with blame 'planner_content', 'planner_format' or
'oracle_wrong' are mistakes to learn from; 'infrastructure' and 'blocked'
attempts are excluded.

usage: python -m learn.episodes <run_dir> [<run_dir> ...] --out episodes.jsonl
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re

INFRA = re.compile(r"ResponseError|status code: 5\d\d|timed? ?out|Connection|EOF|Broken pipe|ManagerLost", re.I)
FORMAT = re.compile(r"did not call the|could not be used|JSON", re.I)
ORACLE_MODULE = {"wrong_parent": "target", "fragments_do_not_add_up": "fragments", "not_a_homolysis": "fragments",
                 "wrong_bond": "fragments", "spin_or_charge": "fragments", "bad_method": "method"}


def _records(path: str) -> list[dict]:
    out = []
    if not os.path.exists(path):
        return out
    for line in open(path, errors="ignore"):
        try:
            d = json.loads(line)
        except Exception:
            continue
        if d.get("type") == "task":
            out.append(d)
    return sorted(out, key=lambda d: d.get("started_at") or 0)


def blame(reply: dict, verdict: dict | None) -> str:
    err = str(reply.get("error", "")) if isinstance(reply, dict) else ""
    if err:
        return "infrastructure" if INFRA.search(err) and not FORMAT.search(err) else "planner_format"
    if verdict is None:
        return "unknown"
    if verdict.get("validator") != "PASS":
        return "planner_content"
    if verdict.get("critic") == "REVISE":
        return "critic_revise"
    return "accepted"


def campaign_episodes(plan_dir: str, meta: dict) -> list[dict]:
    recs = _records(os.path.join(plan_dir, "flowcept_buffer.jsonl"))
    calls = [r for r in recs if r.get("agent_id") == "PlannerAgent" and str(r.get("activity_id", "")).startswith("publish_")]
    verdicts = [r for r in recs if r.get("activity_id") == "round_result"]
    by_mod: dict[str, list[dict]] = {}
    for v in verdicts:
        by_mod.setdefault((v.get("used") or {}).get("module"), []).append(v)
    seen: dict[str, int] = {}
    eps = []
    for c in calls:
        used, gen = c.get("used") or {}, c.get("generated") or {}
        mod = used.get("module")
        k = seen.get(mod, 0)
        seen[mod] = k + 1
        v = by_mod.get(mod, [])[k] if k < len(by_mod.get(mod, [])) else None
        vg = (v or {}).get("generated") or {}
        reply = gen.get("tool_input") if isinstance(gen.get("tool_input"), dict) else {"error": str(gen)}
        eps.append({**meta, "thinking": gen.get("thinking"), "module": mod, "attempt": used.get("attempt", k + 1), "round": ((v or {}).get("used") or {}).get("round"),
                    "system_prompt": used.get("system_prompt"), "user_message": used.get("user_message"),
                    "board": used.get("board"), "lessons_in_prompt": bool(used.get("lessons")),
                    "reply": reply, "validator": vg.get("validator"), "problems": vg.get("problems", []),
                    "critic": vg.get("critic"), "critic_notes": vg.get("critic_notes", []),
                    "blame": blame(reply, vg or None), "flowcept_task_id": c.get("task_id"),
                    "model": used.get("model")})
    outcome = next((r.get("generated") or {} for r in recs if r.get("activity_id") == "plan_outcome"), {})
    grade = meta.get("oracle") or {}
    # Root cause: the oracle names the module behind a wrong accepted plan.
    if grade and not grade.get("correct") and grade.get("error_class") in ORACLE_MODULE:
        bad = ORACLE_MODULE[grade["error_class"]]
        accepted = [e for e in eps if e["module"] == bad and e["blame"] in ("accepted", "critic_revise")]
        if accepted:
            accepted[-1]["blame"] = "oracle_wrong"
        if grade["error_class"] == "wrong_parent":
            # later fragment rejections that only disagree with the wrong parent are not the Planner's fault
            for e in eps:
                if e["module"] == "fragments" and e["blame"] == "planner_content" and \
                        any("but the parent is" in str(p) for p in e["problems"]):
                    e["blame"] = "blocked"
    # A rejected fragments attempt that is the task's correct answer was the board check's
    # mistake (e.g. site-word bugs fixed 2026-11-22), not the Planner's: never a negative.
    exp = meta.get("expected")
    if exp:
        from learn.task_gen import canon
        want = sorted(f["smiles"] for f in exp["fragments"])
        for e in eps:
            if e["module"] == "fragments" and e["blame"] == "planner_content":
                frs = e["reply"].get("fragments") if isinstance(e["reply"], dict) else None
                if isinstance(frs, str):
                    try:
                        frs = json.loads(frs)
                    except Exception:
                        frs = None
                got = sorted(canon(str(f.get("smiles", ""))) or "?" for f in frs) if isinstance(frs, list) and \
                    all(isinstance(f, dict) for f in frs) else None
                if got == want:
                    e["blame"] = "checker_false_reject"
    for e in eps:
        e["plan_outcome"] = outcome.get("outcome")
    return eps


def run_episodes(run_dir: str, tasks: dict | None = None) -> list[dict]:
    oracle = json.load(open(os.path.join(run_dir, "oracle.json"))) if os.path.exists(os.path.join(run_dir, "oracle.json")) else {}
    out = []
    for tdir in sorted(glob.glob(os.path.join(run_dir, "task_*")), key=lambda p: int(p.rsplit("_", 1)[1])):
        i = tdir.rsplit("_", 1)[1]
        o = oracle.get(i) or oracle.get(int(i)) or {}
        base = {"run": os.path.basename(run_dir), "task_index": int(i), "task": o.get("task"), "question": o.get("question"),
                "expected": (tasks or {}).get(o.get("task"), {}).get("expected")}
        out += campaign_episodes(os.path.join(tdir, "plan"), {**base, "phase": "first", "oracle": o.get("first")})
        if os.path.isdir(os.path.join(tdir, "plan_retry")):
            out += campaign_episodes(os.path.join(tdir, "plan_retry"), {**base, "phase": "retry", "oracle": o.get("retry")})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tasks", default=None, help="tasks.json, to recognise correct answers the board check rejected")
    a = ap.parse_args()
    tasks = {t["id"]: t for t in json.load(open(a.tasks))["tasks"]} if a.tasks else None
    n = 0
    with open(a.out, "w") as fh:
        for r in a.runs:
            for e in run_episodes(r, tasks):
                fh.write(json.dumps(e) + "\n")
                n += 1
    print(f"{n} episodes -> {a.out}")


if __name__ == "__main__":
    main()
