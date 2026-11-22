"""Build training data from episodes (learn/episodes.py).

sft.jsonl  module attempts the board accepted in campaigns the oracle graded
           correct: the model imitates what worked.
dpo.jsonl  (chosen, rejected) pairs for the same module and campaign: chosen = an
           accepted attempt from an oracle-correct campaign (or its hinted retry),
           rejected = a planner-caused mistake (rejected by the code checks, or
           accepted but shown wrong by the oracle). Infrastructure failures and
           attempts blocked by an earlier wrong module are never used.

Every prompt is rebuilt as a first try: the board as it stood (question and
accepted modules) without earlier attempts, reviews, oracle hints or lessons,
and the base system prompt without supervisor edits, so the model learns to be
right the first time.

usage: python -m learn.build_dataset episodes.jsonl [...] --tasks tasks.json --out-dir DIR
"""
from __future__ import annotations

import argparse
import json
import os
import re
from collections import defaultdict

from agents.plan_board import module_example, module_tool, render_board
from llm import planner_llm as P

MISTAKE = {"planner_content", "oracle_wrong", "planner_format"}
_PLAIN = re.compile(r"returned plain text instead: (['\"])(.*)\1\); this model", re.S)


def rejected_reply(ep: dict) -> dict | None:
    """The mistaken reply as the model wrote it: a tool call with wrong content, or,
    for a format failure, the plain text it wrote instead of calling the tool."""
    r = ep.get("reply")
    if not isinstance(r, dict):
        return None
    if "error" not in r:
        return {"name": module_tool(ep["module"])["name"], "arguments": r, "thinking": ep.get("thinking") or ""}
    m = _PLAIN.search(r["error"])
    if not m:
        return None
    text = m.group(2).encode().decode("unicode_escape", errors="ignore") if "\\n" in m.group(2) else m.group(2)
    return {"text": text[:6000], "thinking": ep.get("thinking") or ""}


def first_try_prompt(ep: dict) -> tuple[str, str]:
    board = ep.get("board") or {}
    snap = {"question": board.get("question", ep.get("question", "")), "context": "",
            "accepted": board.get("accepted", {}), "entries": []}
    user = render_board(snap, ep["module"]) + "\n\n" + module_example(ep["module"]) \
        + f"\n\nPublish the '{ep['module']}' module now."
    system = (ep.get("system_prompt") or "").split("\n\nSupervisor instruction")[0]
    return system, user


def tool_spec(module: str) -> dict:
    return P._ollama_tool_schema(module_tool(module))


def example(ep: dict) -> dict:
    system, user = first_try_prompt(ep)
    return {"module": ep["module"], "task": ep.get("task"), "run": ep.get("run"),
            "source_task_id": ep.get("flowcept_task_id"),
            "messages": [{"role": "system", "content": system}, {"role": "user", "content": user}],
            "tools": [tool_spec(ep["module"])],
            "reply": {"name": module_tool(ep["module"])["name"], "arguments": ep["reply"],
                      "thinking": ep.get("thinking") or ""}}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("episodes", nargs="+")
    ap.add_argument("--tasks", required=True, help="tasks.json: eval-split tasks are never used for training")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    split = {t["id"]: t["split"] for t in json.load(open(a.tasks))["tasks"]}
    eps = [json.loads(l) for f in a.episodes for l in open(f) if l.strip()]
    eps = [e for e in eps if e.get("system_prompt") and e.get("module") in ("target", "fragments", "method", "rationale")]
    leak = [e for e in eps if split.get(e.get("task")) != "train"]
    eps = [e for e in eps if split.get(e.get("task")) == "train"]
    by_cm = defaultdict(list)
    for e in eps:
        by_cm[(e["run"], e["task_index"], e["module"])].append(e)

    def good(e):
        return e["blame"] == "accepted" and bool((e.get("oracle") or {}).get("correct"))

    sft, dpo = [], []
    for key, group in by_cm.items():
        chosen = [e for e in group if good(e)]
        bad = [e for e in group if e["blame"] in MISTAKE and rejected_reply(e)]
        for c in chosen[:1]:
            sft.append(example(c))
            for r in bad[:4]:   # at most 4 mistakes per (campaign, module): no single campaign dominates
                ex = example(c)
                dpo.append({**{k: ex[k] for k in ("module", "task", "run", "messages", "tools")},
                            "chosen": ex["reply"], "rejected": rejected_reply(r),
                            "rejected_blame": r["blame"], "rejected_problems": r.get("problems", [])[:3],
                            "source_task_ids": [c.get("flowcept_task_id"), r.get("flowcept_task_id")]})
    os.makedirs(a.out_dir, exist_ok=True)
    for name, rows in (("sft", sft), ("dpo", dpo)):
        with open(os.path.join(a.out_dir, f"{name}.jsonl"), "w") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
    stats = {"episodes": len(eps), "excluded_not_train": len(leak), "sft": len(sft), "dpo": len(dpo),
             "dpo_by_module": {m: sum(1 for d in dpo if d["module"] == m) for m in ("target", "fragments", "method", "rationale")},
             "blame": {b: sum(1 for e in eps if e["blame"] == b) for b in sorted({e["blame"] for e in eps})}}
    json.dump(stats, open(os.path.join(a.out_dir, "dataset_stats.json"), "w"), indent=1)
    print(json.dumps(stats))


if __name__ == "__main__":
    main()
