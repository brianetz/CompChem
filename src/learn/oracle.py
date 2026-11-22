"""Grade a BDE plan against the computed answer of its task (see task_gen).

Independent of the plan board's own checks: it knows the right answer, so it
also catches plans the board accepted but that break the wrong bond (case B).
No DFT is run. Each error class has a hint that says what is wrong without
giving the answer, for one retry.
"""
from __future__ import annotations

from collections import Counter

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

from learn.task_gen import canon, homolysis

RDLogger.DisableLog("rdApp.*")

HINTS = {
    "no_plan": None,  # nothing to correct: the board produced no plan
    "wrong_parent": ("An independent check rejected the plan: the parent molecule is not the molecule named "
                     "in the question. Check the parent's SMILES against the molecular formula {formula}."),
    "fragments_do_not_add_up": ("An independent check rejected the plan: the two fragments do not add up to the "
                                "parent molecule. Breaking one bond must keep every atom of the parent."),
    "not_a_homolysis": ("An independent check rejected the plan: the fragments contain the parent's atoms but are "
                        "not what breaking a single bond gives. Each fragment keeps one electron of the broken bond."),
    "wrong_bond": ("An independent check rejected the plan: the fragments come from breaking a different bond "
                   "from the one the question asks for ({site}). Find that bond in the parent and break it."),
    "spin_or_charge": ("An independent check rejected the plan: a species has the wrong charge or spin "
                       "multiplicity. Neutral radicals from a homolysis are doublets."),
    "bad_method": ("An independent check rejected the plan: PySCF cannot run the chosen method ({method}). "
                   "Use a standard functional and basis set."),
}


def _species(sp) -> tuple[str | None, int | None, int | None]:
    if not isinstance(sp, dict):
        return None, None, None
    return canon(str(sp.get("smiles", ""))), sp.get("charge"), sp.get("multiplicity")


def _formula(smiles: str | None) -> str | None:
    m = Chem.MolFromSmiles(smiles or "")
    return rdMolDescriptors.CalcMolFormula(Chem.AddHs(m)) if m else None


def valid_splits(parent_smiles: str) -> set[tuple[str, ...]]:
    mol_h = Chem.AddHs(Chem.MolFromSmiles(parent_smiles))
    out = set()
    for b in mol_h.GetBonds():
        if b.GetBondTypeAsDouble() == 1 and not b.IsInRing():
            frags = homolysis(mol_h, b.GetIdx())
            if frags:
                out.add(tuple(sorted(canon(f) for f in frags)))
    return out


def method_runnable(xc: str | None, basis: str | None, smiles: str | None) -> bool:
    try:
        from pyscf import gto
        from pyscf.dft import libxc
        libxc.parse_xc(str(xc))
        elements = {a.GetSymbol() for a in Chem.AddHs(Chem.MolFromSmiles(smiles)).GetAtoms()}
        for el in elements:
            if not gto.basis.load(str(basis), el):
                return False
        return True
    except Exception:
        return False


def grade(spec: dict | None, task: dict) -> dict:
    exp = task["expected"]
    fmt = {"formula": task.get("formula"), "site": task.get("site")}
    if not spec:
        return {"correct": False, "error_class": "no_plan", "hint": None}
    p_smi, p_q, p_m = _species(spec.get("parent"))
    frags = [_species(f) for f in spec.get("fragments") or []]
    exp_frags = sorted((f["smiles"], f["charge"], f["multiplicity"]) for f in exp["fragments"])
    details = {"parent": p_smi, "fragments": sorted(f[0] or "?" for f in frags),
               "method": f"{spec.get('xc')}/{spec.get('basis')}"}
    frag_key = tuple(sorted(f[0] or "?" for f in frags))
    if p_smi != exp["parent"]["smiles"]:
        cls = "wrong_parent"
    elif len(frags) != 2 or "?" in frag_key or not _adds_up(list(frag_key), p_smi):
        cls = "fragments_do_not_add_up"
    elif frag_key != tuple(f[0] for f in exp_frags):
        cls = "wrong_bond" if frag_key in valid_splits(p_smi) else "not_a_homolysis"
    elif (p_q, p_m) != (0, 1) or sorted(frags) != exp_frags:
        cls = "spin_or_charge"
    elif not method_runnable(spec.get("xc"), spec.get("basis"), p_smi):
        cls = "bad_method"
    else:
        return {"correct": True, "error_class": None, "hint": None, "details": details}
    hint = HINTS[cls].format(method=details["method"], **fmt) if HINTS[cls] else None
    return {"correct": False, "error_class": cls, "hint": hint, "details": details}


def _adds_up(frag_smiles: list[str], parent_smiles: str) -> bool:
    """Do the fragments' atoms sum to the parent's (H included)?"""
    total = Counter()
    for s in frag_smiles:
        m = Chem.MolFromSmiles(s or "")
        if m is None:
            return False
        total.update(a.GetSymbol() for a in Chem.AddHs(m).GetAtoms())
    pm = Chem.MolFromSmiles(parent_smiles or "")
    return pm is not None and total == Counter(a.GetSymbol() for a in Chem.AddHs(pm).GetAtoms())
