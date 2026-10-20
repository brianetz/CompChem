"""
Plan checks: repair and vet a planner's tool call before any compute runs
(the "plan critic" box of the target workflow).

Two parts:

- normalize_raw(): undo shape quirks that carry a correct answer in the
  wrong container. Found 2026-10-19 with Llama 3.1 405B on Ollama, which
  has no schema enforcement: 6 of 8 plans sent 'parent'/'fragments' as a
  string holding the object as JSON or as a Python repr, or flattened the
  parent into 'parent_smiles'/'parent_charge'/'parent_multiplicity'.
- check_spec(): chemistry problems that would make a species fail at RDKit
  or xtb, or a BDE answer the wrong question. Each problem is phrased as an
  instruction, because _generate_spec_via_llm hands it back to the model
  for its one corrective retry.
"""
from __future__ import annotations

import ast
import json
from typing import Any

from analysis.species_checks import (
    atom_balance, electron_count, parse_smiles, spin_allowed,
)

# Spellings that are unambiguous but not valid SMILES.
_SMILES_ALIASES = {"H": "[H]"}


def _decode(value: Any) -> Any:
    """A JSON or Python-literal string -> the object it holds; anything else
    unchanged. A truncated string stays a string (and is rejected later)."""
    if not isinstance(value, str) or value.strip()[:1] not in ("{", "["):
        return value
    for loader in (json.loads, ast.literal_eval):
        try:
            return loader(value)
        except (ValueError, SyntaxError, TypeError, MemoryError, RecursionError):
            continue
    return value


def _fix_species(obj: Any) -> Any:
    obj = _decode(obj)
    if isinstance(obj, dict) and isinstance(obj.get("smiles"), str):
        smiles = obj["smiles"].strip()
        obj = {**obj, "smiles": _SMILES_ALIASES.get(smiles, smiles)}
    return obj


def normalize_raw(raw: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    """Return (repaired copy of the tool input, list of repairs made)."""
    raw = dict(raw)
    repairs = []

    if raw.get("parent") is None and "parent_smiles" in raw:
        raw["parent"] = {
            field: raw.pop(f"parent_{field}")
            for field in ("smiles", "charge", "multiplicity", "label")
            if f"parent_{field}" in raw
        }
        repairs.append("regrouped parent_* fields into 'parent'")

    if isinstance(raw.get("parent"), str) and isinstance(_decode(raw["parent"]), dict):
        repairs.append("decoded 'parent' from a string")
    raw["parent"] = _fix_species(raw.get("parent"))

    frags = raw.get("fragments")
    if isinstance(frags, str) and isinstance(_decode(frags), list):
        repairs.append("decoded 'fragments' from a string")
    frags = _decode(frags)
    if isinstance(frags, dict):  # a single fragment not wrapped in a list
        frags = [frags]
    if isinstance(frags, list):
        if any(isinstance(f, str) and isinstance(_decode(f), dict) for f in frags):
            repairs.append("decoded fragment entries from strings")
        raw["fragments"] = [_fix_species(f) for f in frags]
    return raw, repairs


def _species_problems(sp: dict, name: str, parent_canon: str | None) -> list[str]:
    from rdkit import Chem

    smiles = sp["smiles"]
    mol = parse_smiles(smiles)
    if mol is None:
        return [
            f"{name}: SMILES {smiles!r} is not valid. A bracket atom holds ONE "
            "element plus its own hydrogens ('[CH3]', '[OH]', '[NH2]'); a radical "
            "centre inside a larger molecule is written by bracketing only that "
            "atom, e.g. ethoxy radical 'CC[O]', hydroxymethyl radical '[CH2]O'. "
            "Write formulas as SMILES ('O' for water, not 'H2O')"
        ]
    charge, mult = int(sp.get("charge", 0)), int(sp.get("multiplicity", 1))
    problems = []
    n_e = electron_count(mol, charge)
    if not spin_allowed(n_e, mult):
        problems.append(
            f"{name}: {smiles!r} with charge {charge} has {n_e} electrons "
            f"({'odd' if n_e % 2 else 'even'}), which cannot have multiplicity {mult}"
            + (". A radical fragment has one hydrogen fewer than the parent at the "
               "broken bond - its SMILES must show that"
               if mult == 2 and n_e % 2 == 0 else "")
        )
    if parent_canon is not None and Chem.MolToSmiles(Chem.RemoveHs(mol)) == parent_canon:
        problems.append(f"{name}: {smiles!r} is the parent molecule itself, not a fragment")
    return problems


def check_spec(spec: dict, question: str) -> list[str]:
    """Chemistry problems in a built campaign spec (empty list = OK)."""
    from rdkit import Chem
    from storage.rag import bond_types_in
    from storage.results_db import detect_broken_bond

    if spec.get("analysis_type") == "single_point":
        return _species_problems(spec, "the molecule", None)

    parent = spec["parent"]
    fragments = spec.get("fragments") or []
    problems = _species_problems(parent, "parent", None)
    parent_mol = parse_smiles(parent["smiles"])
    parent_canon = Chem.MolToSmiles(Chem.RemoveHs(parent_mol)) if parent_mol is not None else None
    for i, frag in enumerate(fragments):
        problems += _species_problems(frag, f"fragment {frag.get('label') or i}", parent_canon)
    if problems:
        return problems

    frag_smiles = [f["smiles"] for f in fragments]
    mismatch = atom_balance(parent["smiles"], frag_smiles)
    if mismatch:
        problems.append(
            f"{mismatch}. The fragments must contain exactly the parent's atoms, "
            "split at the one bond that breaks"
        )
    if spec.get("analysis_type") == "bde" and not mismatch:
        wanted = bond_types_in(question)
        broken = detect_broken_bond(parent["smiles"], frag_smiles)
        if wanted and broken and broken not in wanted:
            problems.append(
                f"these fragments come from breaking a {broken} bond, but the "
                f"question asks about {' or '.join(wanted)}"
            )
    return problems
