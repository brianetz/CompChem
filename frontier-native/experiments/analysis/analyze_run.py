"""Summarize one campaign run directory: verdicts, BDEs, species computes, block time."""
import glob
import json
import os
import re
import subprocess
import sys
from collections import Counter


def flowcept_counts(run_dir):
    c = Counter()
    for buf in glob.glob(f"{run_dir}/task_*/flowcept_buffer.jsonl") + glob.glob(f"{run_dir}/species/flowcept_buffer.jsonl"):
        for line in open(buf):
            try:
                d = json.loads(line)
            except Exception:
                continue
            a = d.get("activity_id")
            if a in ("singlepoint_energy_pyscf", "optimize_with_xtb", "generate_campaign_spec_llm"):
                if d.get("subtype") == "ai_model_invocation":
                    c["llm_calls"] += 1
                elif a != "generate_campaign_spec_llm":
                    c[a] += 1
            if a == "species_cache_event":
                c["cache_" + (d.get("generated") or {}).get("outcome", "?")] += 1
    return c


def block_time(run_dir):
    scripts = glob.glob(f"{run_dir}/parsl/*/submit_scripts/parsl.htex_frontier.block-*")
    names = {os.path.basename(s).split(".stdout")[0].split(".stderr")[0] for s in scripts if not s.endswith(".sh")}
    total = 0
    for n in names:
        out = subprocess.run(["sacct", "-n", "-X", "--name", n, "-o", "ElapsedRaw,NNodes", "-S", "2026-10-01"],
                             capture_output=True, text=True).stdout.split()
        if len(out) >= 2:
            total += int(out[0]) * int(out[1])
    return total


def main(run_dir):
    reports = sorted(glob.glob(f"{run_dir}/task_*/agentic_campaign_report.md"))
    verdicts, bdes = Counter(), Counter()
    for r in reports:
        t = open(r).read()
        m = re.search(r"QC verdict: \*\*(\w+)\*\*", t)
        verdicts[m.group(1) if m else "ABORT/OTHER"] += 1
        b = re.search(r"- ([\d.]+) kcal/mol", t)
        meth = re.search(r"pyscf (\S+) single-point", t)
        if b and m and m.group(1) == "PASS":
            bdes[(meth.group(1) if meth else "?", round(float(b.group(1)), 2))] += 1
    n_tasks = len(glob.glob(f"{run_dir}/task_*.log")) or len(glob.glob(f"{run_dir}/task_*"))
    summary_path = f"{run_dir}/run_summary.json"
    rs = json.load(open(summary_path)) if os.path.exists(summary_path) else {}
    out = {
        "run": os.path.basename(run_dir), "mode": rs.get("mode", "baseline"), "tasks": n_tasks,
        "reports": len(reports), "verdicts": dict(verdicts),
        "bde_pass_by_method": {f"{k[0]} {k[1]}": v for k, v in sorted(bdes.items(), key=lambda kv: -kv[1])},
        "counts": dict(flowcept_counts(run_dir)), "block_node_s": block_time(run_dir),
        "stages": rs.get("stages"), "total_wall_s": rs.get("total_wall_s"),
        "db_before": rs.get("db_before"), "db_after": rs.get("db_after"),
    }
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main(sys.argv[1].rstrip("/"))
