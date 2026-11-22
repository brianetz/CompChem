"""Lessons from past mistakes, retrieved for the Planner's prompt.

A lesson is one short, general sentence written from a mistake and the attempt
that fixed it (learn/build_lessons.py). Retrieval filters by plan-board module
and ranks by word overlap with the question, with a bonus for the same bond
class, so a lesson about benzylic C-H reaches the next benzylic question.
"""
from __future__ import annotations

import json
import math
import re
from collections import Counter

_WORD = re.compile(r"[A-Za-z][A-Za-z0-9\-]+")
_STOP = set("the of a an in and to for with what s is kcal mol using cheap xtb prescreen dft refinement "
            "bond dissociation energy".split())


def _bag(text: str) -> Counter:
    return Counter(w.lower() for w in _WORD.findall(text or "") if w.lower() not in _STOP)


def _bond_class(text: str) -> str | None:
    m = re.search(r"\b([A-Z][a-z]?-[A-Z][a-z]?)\b", text or "")
    return m.group(1) if m else None


class LessonMemory:
    def __init__(self, lessons: list[dict]):
        self.lessons = lessons

    @classmethod
    def load(cls, path: str) -> "LessonMemory":
        return cls([json.loads(l) for l in open(path) if l.strip()])

    def retrieve(self, module: str, question: str, k: int = 3) -> list[dict]:
        q = _bag(question)
        qn = math.sqrt(sum(v * v for v in q.values())) or 1.0
        qc = _bond_class(question)
        scored = []
        for les in self.lessons:
            if les.get("module") != module:
                continue
            b = _bag(les.get("question", "") + " " + les.get("lesson", ""))
            bn = math.sqrt(sum(v * v for v in b.values())) or 1.0
            s = sum(q[w] * b[w] for w in q) / (qn * bn)
            if qc and les.get("bond_class") == qc:
                s += 0.25
            scored.append((s, les))
        scored.sort(key=lambda x: -x[0])
        out, seen = [], set()
        for s, les in scored:
            if les["lesson"] in seen or s <= 0:
                continue
            seen.add(les["lesson"])
            out.append(les)
            if len(out) == k:
                break
        return out
