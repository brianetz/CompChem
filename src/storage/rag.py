"""
Retrieval-augmented planning: before the planner LLM answers a question,
pull relevant prior results out of the results database (storage/results_db.py)
and hand them to it as context - which BDEs were already computed for the
bond types / molecules the question is about, with which method, and which
species are already cached.

Two effects: the planner sees real prior numbers ("results inform the next
plan"), and it is steered toward a method that is already cached, which
both turns more species into cache hits and stops the run-to-run basis-set
drift seen on Frontier (6-31G* / def2-SVP / def2-TZVP for one question).

Retrieval is deliberately simple and transparent: exact bond-type matches
parsed from the question ('O-H', 'C–H', ...) plus token overlap with the
questions earlier BDEs were computed for. Disabled with COMPCHEM_RAG=0
(CLI: --no-rag) or when no database is configured.
"""
from __future__ import annotations

import json
import os
import re
from collections import Counter

from storage.results_db import ResultsDB, get_db

_BOND_RE = re.compile(r"\b([A-Z][a-z]?)\s*([-–—‐=#≡])\s*([A-Z][a-z]?)\b")
_BOND_ORDER = {"=": "=", "#": "#", "≡": "#"}  # every dash variant means a single bond
_TOKEN_RE = re.compile(r"[a-z0-9]+")
_STOP = {"the", "a", "an", "of", "for", "in", "and", "with", "what", "whats", "is", "s", "using",
         "energy", "bond", "dissociation", "kcal", "mol", "cheap", "dft", "xtb", "refinement",
         "prescreen", "compute", "calculate", "please", "to", "on", "by", "at"}


def _normalize_bond(a: str, b: str, sym: str = "-") -> str | None:
    from rdkit import Chem

    table = Chem.GetPeriodicTable()
    try:
        za, zb = table.GetAtomicNumber(a), table.GetAtomicNumber(b)
    except Exception:
        return None
    if za == 0 or zb == 0:
        return None
    return f"{a}{sym}{b}" if za >= zb else f"{b}{sym}{a}"


def bond_types_in(question: str) -> list[str]:
    out = []
    for a, sym, b in _BOND_RE.findall(question):
        bt = _normalize_bond(a, b, _BOND_ORDER.get(sym, "-"))
        if bt and bt not in out:
            out.append(bt)
    return out


def _tokens(text: str | None) -> set[str]:
    return {t for t in _TOKEN_RE.findall((text or "").lower()) if t not in _STOP}


def _method_label(method: str) -> str:
    # 'pyscf:b3lyp/6-31g*|xtb-gfn2|smiles-embed' -> 'b3lyp/6-31g* (pyscf)'
    head = method.split("|", 1)[0]
    engine, _, xc_basis = head.partition(":")
    return f"{xc_basis} ({engine})"


def retrieve(question: str, db: ResultsDB, k: int = 8) -> list[dict]:
    """Top-k prior PASS BDE rows for this question, best first."""
    bonds = bond_types_in(question)
    q_tokens = _tokens(question)
    candidates = db.query_bdes(qc_flag="PASS", limit=500)
    scored = []
    for row in candidates:
        overlap = q_tokens & _tokens(row["source_question"])
        jaccard = len(overlap) / max(1, len(q_tokens | _tokens(row["source_question"])))
        score = (1.0 if row["bond_type"] and row["bond_type"] in bonds else 0.0) + jaccard
        if score >= 0.2:
            scored.append((score, row["created_at"], row))
    scored.sort(key=lambda t: (t[0], t[1]), reverse=True)
    return [row | {"score": round(score, 3)} for score, _, row in scored[:k]]


def build_retrieval_context(question: str, db: ResultsDB | None = None) -> tuple[str, dict]:
    """(context text for the planner prompt, provenance metadata). Empty text
    when retrieval is off or nothing relevant is stored."""
    meta = {"enabled": False, "bond_types": [], "record_ids": [], "recommended_method": None}
    if os.environ.get("COMPCHEM_RAG", "1") == "0":
        return "", meta
    db = db if db is not None else get_db()
    if db is None:
        return "", meta
    meta["enabled"] = True
    meta["bond_types"] = bond_types_in(question)
    rows = retrieve(question, db)
    if not rows:
        return "", meta

    # Recommend the method whose cached species have been reused most - the
    # de facto standard for this kind of question - not just the most recent.
    methods = Counter()
    for r in rows:
        parent = db.get_species(r["parent_key"]) or {}
        methods[r["method"]] += 1 + (parent.get("hits") or 0)
    top_method, _ = methods.most_common(1)[0]
    meta["record_ids"] = [r["id"] for r in rows]
    meta["recommended_method"] = _method_label(top_method)

    lines = [
        "Context retrieved from the shared results database (prior calculations "
        "relevant to this question, most relevant first):",
    ]
    for r in rows:
        frags = " + ".join(json.loads(r["fragment_smiles"]))
        lines.append(
            f"- {r['bond_type'] or 'bond ?'} in {r['parent_smiles']} -> {frags}: "
            f"{r['bde_kcal_mol']:.2f} kcal/mol at {_method_label(r['method'])}, critic {r['qc_flag']}"
            f" (asked as: {r['source_question']!r})"
        )
    lines += [
        "",
        f"Every species above is already cached at its listed method, so reusing that method costs "
        f"no new compute and keeps results comparable. Unless the question explicitly asks for a "
        f"different functional or basis set, use {meta['recommended_method'].split(' ')[0]} "
        f"(xc/basis) for this campaign. Still specify every species explicitly.",
    ]
    return "\n".join(lines), meta
