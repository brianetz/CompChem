"""Per-run time breakdown for every CompChem campaign run.

For each run: setup (model, planner mode, critic, RAG, campaigns, time limit),
wall-clock zones (startup, plan, plan critic, execution, result generation,
unused time limit) and per-campaign averages by stage, from run_summary.json,
the Parsl block script, task command lines and Flowcept records.

usage: stage_times.py <out.json>
"""
import glob
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from datetime import datetime

ROOT = "/lustre/orion/stf053/scratch/akhilesh/compchem/runs/campaigns"

LABELS = {
    "run_20261014_162020": "Baseline", "run_20261018_224316": "A", "run_20261018_230942": "B",
    "run_20261018_233407": "D", "run_20261018_235544": "C", "run_20261019_002148": "V",
    "run_20261019_081257": "R1", "run_20261019_092019": "R2", "run_20261019_121942": "L405",
    "run_20261020_100129": "L405v2", "run_20261020_102220": "F1", "run_20261020_114059": "F2",
    "run_20261020_122720": "F3", "run_20261020_112921": "G3", "run_20261020_120916": "G1",
    "run_20261020_124021": "G2", "run_20261020_182457": "B1", "run_20261020_190936": "B2",
    "run_20261021_115950": "B1v2", "run_20261021_122925": "B2v2", "run_20261021_132444": "NE1",
    "run_20261021_141911": "NE2",
}
EXEC = {"optimize_with_xtb", "singlepoint_energy_pyscf", "singlepoint_energy_nwchem", "get_lowest_energy_conformer"}
RESULT = {"review_analysis", "write_analysis_report", "render_structure_png"}


def task_command(run):
    for f in sorted(glob.glob(f"{run}/plan_*.log")) + sorted(glob.glob(f"{run}/task_*.log")):
        txt = open(f, errors="ignore").read(20000)
        m = re.search(r"executable follows <--\n(.*?)\n--> end executable", txt, re.S)
        if m:
            return m.group(1)
    return ""


def block_info(run):
    out = {"limit_min": None, "partition": None}
    for s in glob.glob(f"{run}/parsl/*/submit_scripts/parsl.htex_frontier.block-*"):
        if s.endswith((".stdout", ".stderr", ".sh")):
            continue
        txt = open(s, errors="ignore").read()
        m = re.search(r"#SBATCH --time=(\d+)", txt)
        out["limit_min"] = int(m.group(1)) if m else None
        m = re.search(r"#SBATCH --partition=(\S+)", txt)
        out["partition"] = m.group(1) if m else None
        name = os.path.basename(s)
        acct = subprocess.run(["sacct", "-n", "-X", "--name", name, "-o", "ElapsedRaw,Start", "-S", "2026-10-01"],
                              capture_output=True, text=True).stdout.split()
        if acct:
            out["block_elapsed_s"] = int(acct[0])
            out["block_start"] = acct[1]
        break
    return out


def flowcept_records(run):
    recs = []
    for b in glob.glob(f"{run}/**/flowcept_buffer.jsonl", recursive=True):
        where = "species" if "/species/" in b else ("plan" if "/plan/" in b else "task")
        for line in open(b, errors="ignore"):
            try:
                d = json.loads(line)
            except Exception:
                continue
            if d.get("started_at") and d.get("ended_at"):
                d["_where"] = where
                d["_task"] = re.search(r"/(task_\d+|species_\d+)/", b).group(1) if re.search(r"/(task_\d+|species_\d+)/", b) else "?"
                recs.append(d)
    return recs


def zone_of(d):
    a, sub, agent = d.get("activity_id"), d.get("subtype"), d.get("agent_id")
    if sub == "ai_model_invocation":
        return "critic" if agent == "PlanCriticAgent" else "plan"
    if a in EXEC:
        return "exec"
    if a in RESULT:
        return "result"
    return None


def analyse(run_dir):
    name = os.path.basename(run_dir)
    summ = json.load(open(f"{run_dir}/run_summary.json")) if os.path.exists(f"{run_dir}/run_summary.json") else {}
    cmd = task_command(run_dir)
    model = (re.search(r"--model (\S+)", cmd) or [None, "?"])[1]
    model = "Llama 3.1 405B" if "405b" in model else ("nemotron-3-nano" if "nemotron" in model else model)
    planner = "board" if "--plan-board" in cmd else "chat" if "--plan-chat" in cmd else "single"
    critic = "--plan-critic-skill" in cmd or "--plan-critic" in cmd
    rag = "--no-rag" not in cmd and bool(summ.get("db"))
    blk = block_info(run_dir)
    recs = flowcept_records(run_dir)
    n_campaigns = summ.get("n_tasks") or len(glob.glob(f"{run_dir}/task_*/"))

    # per-campaign summed seconds by zone (the top-level LLM invocation records
    # only, so nested wrappers aren't double counted)
    per = defaultdict(float)
    for d in recs:
        z = zone_of(d)
        if z:
            per[z] += d["ended_at"] - d["started_at"]
    per_campaign = {z: round(v / max(1, n_campaigns), 1) for z, v in per.items()}

    run_start = datetime.strptime(name[4:], "%Y%m%d_%H%M%S").timestamp()
    first = min((d["started_at"] for d in recs), default=None)
    startup = (first - run_start) if first else 0.0
    st = summ.get("stages", {})
    total = summ.get("total_wall_s") or blk.get("block_elapsed_s") or 0.0
    llm = per.get("plan", 0) + per.get("critic", 0)
    critic_share = per.get("critic", 0) / llm if llm else 0.0
    if "plan" in st:  # staged: stage walls are measured directly
        decision = max(0.0, st["plan"]["wall_s"] - startup)
        execution = st["compute"]["wall_s"]
        result = st["assemble"]["wall_s"]
        split = "measured"
    else:  # one-pass runs: stages interleave across workers; split by share of summed time
        rest = max(0.0, total - startup)
        tot = sum(per.get(z, 0) for z in ("plan", "critic", "exec", "result")) or 1
        decision = rest * (per.get("plan", 0) + per.get("critic", 0)) / tot
        execution = rest * per.get("exec", 0) / tot
        result = rest - decision - execution
        split = "estimated"
    wall = {
        "startup": round(startup, 1),
        "plan": round(decision * (1 - critic_share), 1),
        "critic": round(decision * critic_share, 1),
        "exec": round(execution, 1),
        "result": round(result, 1),
    }
    limit_s = (blk.get("limit_min") or 0) * 60
    used = sum(wall.values())
    wall["unused"] = round(max(0.0, limit_s - used), 1) if limit_s else None
    reports = len(glob.glob(f"{run_dir}/task_*/agentic_campaign_report.md"))
    return {
        "run": name, "label": LABELS.get(name, name), "model": model, "planner": planner, "critic": critic,
        "rag": rag, "mode": summ.get("mode", "baseline"), "campaigns": n_campaigns, "reports": reports,
        "limit_min": blk.get("limit_min"), "partition": blk.get("partition"),
        "total_s": round(total, 1), "wall": wall, "split": split, "per_campaign": per_campaign,
        "question": (summ.get("question") or "")[:60],
    }


def main():
    out = []
    for run in sorted(glob.glob(f"{ROOT}/run_*")):
        if os.path.basename(run) not in LABELS:
            continue
        out.append(analyse(run))
    json.dump(out, open(sys.argv[1], "w"), indent=1)
    for r in out:
        print(f"{r['label']:8s} {r['model']:16s} {r['planner']:6s} critic={r['critic']!s:5s} rag={r['rag']!s:5s} "
              f"n={r['campaigns']:>3} limit={r['limit_min']} total={r['total_s']:>7} {r['split']:9s} wall={r['wall']} per={r['per_campaign']}")


if __name__ == "__main__":
    main()
