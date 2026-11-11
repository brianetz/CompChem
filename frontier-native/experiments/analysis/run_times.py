"""Time breakdown of staged board runs, for the 'Where the time goes' charts (section 12).

Wall clock per run: Slurm queue (run folder created -> block start), Ollama start (block start
-> first server listening), planning and plan critic (planning wave, split by the two agents'
share of LLM-call time), execution (compute wave), result generation (assemble wave).

Per campaign: every Flowcept LLM call is matched to the Ollama request that served it (same end
time); the request's prompt-eval + eval time is the model working, the rest of the call is
waiting for a free model slot. Works with one or several Ollama servers (ollama_serve*.log).

usage: run_times.py <label>=<run_dir> ...  [> json on stdout with --json]
"""
import glob
import json
import os
import re
import subprocess
import sys
from datetime import datetime

TIMING = re.compile(r"(prompt eval time|\beval time) =\s+([\d.]+) ms")
GIN = re.compile(r"\[GIN\] (\d{4}/\d\d/\d\d - \d\d:\d\d:\d\d) \| (\d+) \|\s+(\S+) \|.*\"/api/chat\"")
LISTEN = re.compile(r'time=(\S+) level=INFO source=routes.go:\d+ msg="Listening on')


def ollama_requests(run):
    """[(end_ts, working_s)] for every /api/chat request in every server log."""
    reqs, listen = [], []
    for log in glob.glob(f"{run}/ollama_serve*.log"):
        pending = 0.0
        for line in open(log, errors="ignore"):
            m = TIMING.search(line)
            if m:
                pending += float(m.group(2)) / 1000
                continue
            m = GIN.search(line)
            if m:
                end = datetime.strptime(m.group(1), "%Y/%m/%d - %H:%M:%S").timestamp()
                reqs.append((end, pending))
                pending = 0.0
                continue
            m = LISTEN.search(line)
            if m:
                listen.append(datetime.fromisoformat(m.group(1)).timestamp())
    return sorted(reqs), (min(listen) if listen else None)


def llm_calls(run):
    calls = []
    for b in glob.glob(f"{run}/task_*/plan/flowcept_buffer.jsonl") + glob.glob(f"{run}/task_*/flowcept_buffer.jsonl"):
        for line in open(b, errors="ignore"):
            try:
                d = json.loads(line)
            except Exception:
                continue
            if d.get("subtype") == "ai_model_invocation" and d.get("started_at") and d.get("ended_at"):
                calls.append({"agent": d.get("agent_id"), "start": d["started_at"], "end": d["ended_at"]})
    return calls


def block_start(run):
    for s in glob.glob(f"{run}/parsl/*/submit_scripts/parsl.htex_frontier.block-*"):
        if s.endswith((".stdout", ".stderr", ".sh")):
            continue
        out = subprocess.run(["sacct", "-n", "-X", "--name", os.path.basename(s), "-o", "Start", "-S", "2026-10-01"],
                             capture_output=True, text=True).stdout.split()
        if out:
            return datetime.fromisoformat(out[0]).timestamp()
    return None


def analyse(run):
    s = json.load(open(f"{run}/run_summary.json"))
    n = s["n_tasks"]
    run_start = datetime.strptime(os.path.basename(run)[4:19], "%Y%m%d_%H%M%S").timestamp()
    bstart = block_start(run)
    reqs, listen = ollama_requests(run)
    calls = llm_calls(run)
    # match calls to requests: nearest end time (GIN is to the second), each request used once
    used = set()
    work = {"plan": 0.0, "critic": 0.0}
    dur = {"plan": 0.0, "critic": 0.0}
    matched = 0
    for c in sorted(calls, key=lambda c: c["end"]):
        kind = "critic" if c["agent"] == "PlanCriticAgent" else "plan"
        dur[kind] += c["end"] - c["start"]
        best = min((i for i in range(len(reqs)) if i not in used), key=lambda i: abs(reqs[i][0] - c["end"]), default=None)
        if best is not None and abs(reqs[best][0] - c["end"]) <= 2.5:
            used.add(best)
            work[kind] += reqs[best][1]
            matched += 1
    llm = dur["plan"] + dur["critic"]
    queue = max(0.0, bstart - run_start) if bstart else 0.0
    start = max(0.0, listen - bstart) if (listen and bstart) else 0.0
    st = s["stages"]
    decision = max(0.0, st["plan"]["wall_s"] - queue - start)
    cshare = dur["critic"] / llm if llm else 0.0
    wall = {"queue": queue, "startup": start, "plan": decision * (1 - cshare), "critic": decision * cshare,
            "exec": st["compute"]["wall_s"], "result": st["assemble"]["wall_s"]}
    wall = {k: round(v / 60, 2) for k, v in wall.items()}
    exec_per = 0.0  # execution seconds per campaign: species compute summed over workers / campaigns
    for b in glob.glob(f"{run}/species/species_*/flowcept_buffer.jsonl"):
        for line in open(b, errors="ignore"):
            try:
                d = json.loads(line)
            except Exception:
                continue
            if d.get("activity_id") in {"optimize_with_xtb", "singlepoint_energy_pyscf"} and d.get("ended_at"):
                exec_per += d["ended_at"] - d["started_at"]
    per = {"plan_work": work["plan"] / n / 60, "critic_work": work["critic"] / n / 60,
           "slot_wait": max(0.0, llm - work["plan"] - work["critic"]) / n / 60, "exec_s": exec_per / n}
    per = {k: round(v, 2) for k, v in per.items()}
    return {"run": os.path.basename(run), "n": n, "wall_min": wall, "used_min": round(sum(wall.values()), 1),
            "per_campaign": per, "per_total_min": round(llm / n / 60, 1), "calls": len(calls), "matched": matched}


if __name__ == "__main__":
    out = {}
    for arg in [a for a in sys.argv[1:] if a != "--json"]:
        label, run = arg.split("=", 1)
        out[label] = r = analyse(run)
        if "--json" not in sys.argv:
            print(f"{label:9} {r['run']} used {r['used_min']} min  wall {r['wall_min']}  per {r['per_campaign']} "
                  f"(calls {r['per_total_min']} min/campaign)  matched {r['matched']}/{r['calls']}")
    if "--json" in sys.argv:
        json.dump(out, sys.stdout, indent=1)
