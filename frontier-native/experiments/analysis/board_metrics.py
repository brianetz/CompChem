"""Per-run board metrics: correct, rounds per campaign, BDE per campaign, plan/total wall.
usage: board_metrics.py <label>=<run_dir> ...   (run classify_board.py first)"""
import json, re, sys, glob, os, statistics as st
out = {}
for arg in sys.argv[1:]:
    label, run = arg.split("=", 1)
    s = json.load(open(f"{run}/run_summary.json"))
    rows = json.load(open(f"{run}/.board_scores.json"))
    for r in rows:
        txt = open(f"{run}/plan_{r['task']}.log", errors="ignore").read()
        m = re.findall(r"after (\d+) reviewed round", txt)
        r["rounds"] = int(m[-1]) if m else None
    rd = [r["rounds"] for r in rows if r["rounds"]]
    out[label] = {"run": os.path.basename(run), "n": len(rows),
                  "correct": sum(r["outcome"] == "correct" for r in rows),
                  "rounds_mean": round(st.mean(rd), 1) if rd else None,
                  "plan_min": round(s["stages"]["plan"]["wall_s"] / 60, 1),
                  "total_min": round(s["total_wall_s"] / 60, 1),
                  "rows": [{k: r.get(k) for k in ("task", "molecule", "outcome", "why", "bde", "method", "rounds")} for r in rows]}
    o = out[label]
    print(f"{label:6} {o['run']}  {o['correct']}/{o['n']}  rounds {o['rounds_mean']}  plan {o['plan_min']} min  total {o['total_min']} min")
    for r in o["rows"]:
        print(f"     {r['task']:>2} {str(r['molecule']):20} {r['outcome']:8} {r['bde']}  {r['method']}  r{r['rounds']}  {r.get('why') or ''}")
json.dump(out, open(os.environ.get("OUT", "/dev/null"), "w"), indent=1)
