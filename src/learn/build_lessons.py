"""Turn the Planner's fixed mistakes into short, general lessons (mistake memory).

A mistake qualifies when the same campaign later fixed it: a module attempt the
code checks rejected followed by an accepted one, or a plan the oracle showed
wrong followed by a correct hinted retry. A different, stronger model than the
one that erred writes one general sentence from (wrong reply, what the check
said, corrected reply). Only training-split tasks are used.

Runs where Ollama servers are up: OLLAMA_HOSTS=host:port,host:port,...
usage: python -m learn.build_lessons episodes.jsonl [...] --tasks tasks.json --model nemotron-3-super:120b --out lessons.jsonl
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import re
from collections import defaultdict

from learn.task_gen import _bond_class

PROMPT = """A planning AI sets up bond dissociation energy (BDE) calculations, one module at a time
(target, fragments, method, rationale). It made a mistake that a check caught, then fixed it.

Question: {question}
Module: {module}
Wrong attempt: {wrong}
What the check said: {feedback}
Corrected attempt: {right}

Write ONE lesson, at most 30 words, that would stop this mistake on SIMILAR questions about other
molecules. State the general rule (about the bond, the site, SMILES, radicals, spin, or the method), not
this molecule's answer. Reply with the lesson sentence only."""


def _short(x, n=600):
    s = json.dumps(x) if not isinstance(x, str) else x
    return s if len(s) <= n else s[:n] + "..."


def candidates(eps: list[dict], split: dict) -> list[dict]:
    by_cm = defaultdict(list)
    for e in eps:
        if split.get(e.get("task")) == "train":
            by_cm[(e["run"], e["task_index"], e["phase"], e["module"])].append(e)
    out = []
    for (run, ti, phase, mod), group in by_cm.items():
        group.sort(key=lambda e: e.get("attempt") or 0)
        fixed = next((e for e in group if e["blame"] == "accepted"), None)
        for e in group:
            if fixed and e["blame"] == "planner_content" and (e.get("attempt") or 0) < (fixed.get("attempt") or 0):
                out.append({"module": mod, "question": e.get("question"), "wrong": e.get("reply"),
                            "feedback": "; ".join(map(str, e.get("problems") or []))[:600], "right": fixed.get("reply"),
                            "kind": "code_check", "sources": [e.get("flowcept_task_id"), fixed.get("flowcept_task_id")]})
    # oracle-wrong first plans fixed by the hinted retry
    first = {(e["run"], e["task_index"], e["module"]): e for e in eps if e["blame"] == "oracle_wrong" and e["phase"] == "first"}
    for (run, ti, mod), bad in first.items():
        if split.get(bad.get("task")) != "train":
            continue
        retry_ok = [e for e in eps if (e["run"], e["task_index"], e["module"]) == (run, ti, mod) and e["phase"] == "retry"
                    and e["blame"] == "accepted" and (e.get("oracle") or {}).get("correct")]
        if retry_ok:
            out.append({"module": mod, "question": bad.get("question"), "wrong": bad.get("reply"),
                        "feedback": (retry_ok[-1].get("board") or {}).get("context") or (bad.get("oracle") or {}).get("error_class"),
                        "right": retry_ok[-1].get("reply"), "kind": "oracle",
                        "sources": [bad.get("flowcept_task_id"), retry_ok[-1].get("flowcept_task_id")]})
    return out


def write_lesson(c: dict, host: str, model: str) -> str | None:
    import ollama
    client = ollama.Client(host=host)
    r = client.chat(model=model, messages=[{"role": "user", "content": PROMPT.format(
        question=c["question"], module=c["module"], wrong=_short(c["wrong"]), feedback=_short(c["feedback"], 700),
        right=_short(c["right"]))}], think=False, options={"temperature": 0.2})
    text = (r.message.content or "").strip().splitlines()
    text = next((l.strip(" -*\"") for l in text if l.strip()), "")
    return text if 8 <= len(text) <= 400 else None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("episodes", nargs="+")
    ap.add_argument("--tasks", required=True)
    ap.add_argument("--model", default="nemotron-3-super:120b")
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-module", type=int, default=80)
    a = ap.parse_args()
    split = {t["id"]: t["split"] for t in json.load(open(a.tasks))["tasks"]}
    eps = [json.loads(l) for f in a.episodes for l in open(f) if l.strip()]
    cands = candidates(eps, split)
    hosts = [h for h in os.environ.get("OLLAMA_HOSTS", os.environ.get("OLLAMA_HOST", "127.0.0.1:11434")).split(",") if h]
    print(f"{len(cands)} fixed mistakes; writing lessons with {a.model} on {len(hosts)} server(s)")
    lessons, seen, per_mod = [], set(), defaultdict(int)
    with cf.ThreadPoolExecutor(max_workers=len(hosts)) as pool:
        futs = {pool.submit(write_lesson, c, hosts[i % len(hosts)], a.model): c for i, c in enumerate(cands)}
        for f in cf.as_completed(futs):
            c = futs[f]
            try:
                text = f.result()
            except Exception as exc:
                print("lesson failed:", exc)
                continue
            key = re.sub(r"\W+", " ", (text or "").lower()).strip()
            if not text or key in seen or per_mod[c["module"]] >= a.per_module:
                continue
            seen.add(key)
            per_mod[c["module"]] += 1
            lessons.append({"module": c["module"], "lesson": text, "question": c["question"], "kind": c["kind"],
                            "bond_class": _bond_class(re.sub(r".*the (.*?) bond dissociation.*", r"\1", c["question"] or "")),
                            "sources": c["sources"]})
    with open(a.out, "w") as fh:
        for les in lessons:
            fh.write(json.dumps(les) + "\n")
    print(f"{len(lessons)} lessons -> {a.out}: {dict(per_mod)}")


if __name__ == "__main__":
    main()
