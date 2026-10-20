"""
Chemistry checks on planned species, shared by the plan checks
(llm/plan_checks.py, before any compute) and the BDE results critic
(skills/bde.review, after compute).

Found 2026-10-19 (methanol runs): every failed species died before DFT -
at RDKit (unparseable SMILES such as '[CH3O]' or bare 'H') or at xtb (a
doublet with an even electron count, e.g. methanol 'CO' passed off as the
methoxy radical) - and the results critic passed a "BDE" whose fragments
did not add up to the parent (24 763 kcal/mol). These checks catch all of
that from the SMILES alone.
"""
from __future__ import annotations

from collections import Counter


def parse_smiles(smiles: str):
    """RDKit Mol with explicit hydrogens, or None if RDKit can't read it.
    Same relaxed fallback as engines.rdkit_to_ase.smiles_to_ase_atoms, so a
    SMILES accepted here is one the compute step can build."""
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
        except Exception:
            return None
    return Chem.AddHs(mol)


def element_counts(mol) -> Counter:
    return Counter(atom.GetSymbol() for atom in mol.GetAtoms())


def electron_count(mol, charge: int) -> int:
    return sum(atom.GetAtomicNum() for atom in mol.GetAtoms()) - int(charge)


def spin_allowed(n_electrons: int, multiplicity: int) -> bool:
    """An odd electron count needs an even multiplicity, and vice versa."""
    return multiplicity >= 1 and (n_electrons % 2) != (multiplicity % 2)


def formula(counts: Counter) -> str:
    """Hill-order formula ('CH4O'), for readable messages."""
    def part(el):
        n = counts[el]
        return el + (str(n) if n > 1 else "")

    if "C" in counts:
        order = ["C", "H"] + sorted(el for el in counts if el not in ("C", "H"))
    else:
        order = sorted(counts)
    return "".join(part(el) for el in order if counts[el])


def atom_balance(parent_smiles: str, fragment_smiles: list[str]) -> str | None:
    """None if the fragments contain exactly the parent's atoms (including
    hydrogens), else a description of the mismatch. Also None if any SMILES
    is unparseable - that is reported separately."""
    parent = parse_smiles(parent_smiles)
    frags = [parse_smiles(s) for s in fragment_smiles]
    if parent is None or any(f is None for f in frags):
        return None
    want = element_counts(parent)
    got = sum((element_counts(f) for f in frags), Counter())
    if got == want:
        return None
    missing, extra = want - got, got - want
    bits = []
    if missing:
        bits.append(f"missing {formula(missing)}")
    if extra:
        bits.append(f"extra {formula(extra)}")
    return (f"fragments add up to {formula(got)} but the parent is {formula(want)} "
            f"({', '.join(bits)})")
