"""BDE questions with computed answers.

For every molecule in MOLECULES, break each non-ring single bond, keep one bond
per set of equivalent bonds (same fragments), describe its site in words, and
keep it only if the description picks out that bond alone. The expected answer
(parent, fragments, charges, multiplicities) comes from RDKit, so a plan can be
graded without any DFT.

usage: python -m learn.task_gen --out tasks.json [--eval-frac 0.25] [--holdout-class C-Cl] [--seed 7]
"""
from __future__ import annotations

import argparse
import json
import random
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

# (name, SMILES). Small closed-shell molecules of H, C, N, O, F, S, Cl, Br (no iodine:
# 6-31G* has no I basis). The 16 molecules of the earlier board v2 test sets are included.
MOLECULES = [
    ("hydrogen", "[H][H]"), ("water", "O"), ("hydrogen fluoride", "F"), ("hydrogen chloride", "Cl"),
    ("hydrogen bromide", "Br"), ("hydrogen sulfide", "S"), ("ammonia", "N"), ("hydrogen peroxide", "OO"),
    ("hydrazine", "NN"), ("hydroxylamine", "NO"),
    ("methane", "C"), ("ethane", "CC"), ("propane", "CCC"), ("n-butane", "CCCC"), ("isobutane", "CC(C)C"),
    ("neopentane", "CC(C)(C)C"), ("cyclopropane", "C1CC1"), ("cyclopentane", "C1CCCC1"), ("cyclohexane", "C1CCCCC1"),
    ("ethylene", "C=C"), ("propene", "C=CC"), ("1-butene", "C=CCC"), ("isobutylene", "C=C(C)C"),
    ("1,3-butadiene", "C=CC=C"), ("cyclohexene", "C1=CCCCC1"), ("acetylene", "C#C"), ("propyne", "C#CC"),
    ("benzene", "c1ccccc1"), ("toluene", "Cc1ccccc1"), ("ethylbenzene", "CCc1ccccc1"), ("p-xylene", "Cc1ccc(C)cc1"),
    ("cumene", "CC(C)c1ccccc1"), ("styrene", "C=Cc1ccccc1"), ("phenol", "Oc1ccccc1"), ("aniline", "Nc1ccccc1"),
    ("anisole", "COc1ccccc1"), ("benzaldehyde", "O=Cc1ccccc1"), ("benzyl alcohol", "OCc1ccccc1"),
    ("thiophenol", "Sc1ccccc1"), ("chlorobenzene", "Clc1ccccc1"), ("fluorobenzene", "Fc1ccccc1"),
    ("bromobenzene", "Brc1ccccc1"), ("pyridine", "c1ccncc1"), ("benzyl chloride", "ClCc1ccccc1"),
    ("methanol", "CO"), ("ethanol", "CCO"), ("1-propanol", "CCCO"), ("2-propanol", "CC(C)O"),
    ("tert-butanol", "CC(C)(C)O"), ("allyl alcohol", "C=CCO"), ("ethylene glycol", "OCCO"),
    ("dimethyl ether", "COC"), ("diethyl ether", "CCOCC"), ("tetrahydrofuran", "C1CCOC1"),
    ("formaldehyde", "C=O"), ("acetaldehyde", "CC=O"), ("propanal", "CCC=O"), ("acetone", "CC(C)=O"),
    ("butanone", "CCC(C)=O"), ("formic acid", "OC=O"), ("acetic acid", "CC(=O)O"), ("methyl formate", "COC=O"),
    ("methyl acetate", "COC(C)=O"), ("acetyl chloride", "CC(=O)Cl"), ("acetamide", "CC(N)=O"), ("formamide", "NC=O"),
    ("methylamine", "CN"), ("dimethylamine", "CNC"), ("trimethylamine", "CN(C)C"), ("ethylamine", "CCN"),
    ("acetonitrile", "CC#N"), ("hydrogen cyanide", "C#N"),
    ("methanethiol", "CS"), ("ethanethiol", "CCS"), ("dimethyl sulfide", "CSC"), ("dimethyl disulfide", "CSSC"),
    ("fluoromethane", "CF"), ("chloromethane", "CCl"), ("bromomethane", "CBr"), ("dichloromethane", "ClCCl"),
    ("chloroform", "ClC(Cl)Cl"), ("trifluoromethane", "FC(F)F"), ("chloroethane", "CCCl"), ("bromoethane", "CCBr"),
    ("vinyl chloride", "C=CCl"), ("methyl hydroperoxide", "COO"), ("dimethyl peroxide", "COOC"),
]

QUESTION = ("What's the {site} bond dissociation energy of {name} ({formula}) in kcal/mol, "
            "using a cheap xtb prescreen and a DFT refinement?")
QUESTION_BOND = ("What's the bond dissociation energy of the {site} in {name} ({formula}), in kcal/mol, "
                 "using a cheap xtb prescreen and a DFT refinement?")


def canon(smiles_or_mol) -> str | None:
    """Canonical SMILES used for every comparison (explicit H atoms folded in)."""
    mol = Chem.MolFromSmiles(smiles_or_mol) if isinstance(smiles_or_mol, str) else smiles_or_mol
    if mol is None:
        return None
    return Chem.MolToSmiles(Chem.RemoveHs(mol, sanitize=False) if mol.GetNumAtoms() > 1 else mol)


def homolysis(mol_h: Chem.Mol, bond_idx: int) -> list[Chem.Mol] | None:
    """Break one single bond homolytically: each end keeps one electron."""
    rw = Chem.RWMol(mol_h)
    b = rw.GetBondWithIdx(bond_idx)
    a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
    rw.RemoveBond(a1, a2)
    for a in (a1, a2):
        at = rw.GetAtomWithIdx(a)
        at.SetNumRadicalElectrons(at.GetNumRadicalElectrons() + 1)
        at.SetNoImplicit(True)
    m = rw.GetMol()
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return None
    frags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True)
    return list(frags) if len(frags) == 2 else None


def _multiplicity(mol: Chem.Mol) -> int:
    return sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()) + 1


def _group(atom: Chem.Atom) -> str:
    """Words for the group an atom belongs to (for heavy-atom bond descriptions)."""
    if atom.GetIsAromatic():
        return "ring"
    sym = atom.GetSymbol()
    if sym == "C":
        nh = atom.GetTotalNumHs(includeNeighbors=True)
        if any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(atom).GetSymbol() == "O" for b in atom.GetBonds()):
            return "carbonyl"
        return {3: "methyl", 2: "methylene", 1: "methine"}.get(nh, "carbon")
    return {"O": "oxygen", "N": "nitrogen", "S": "sulfur"}.get(sym, sym)


def _ch_descriptors(c: Chem.Atom, mol: Chem.Mol) -> list[str]:
    """Increasingly specific names for a C-H site."""
    nbrs = [n for n in c.GetNeighbors() if n.GetSymbol() != "H"]
    out = []
    if c.GetIsAromatic():
        ring = mol.GetRingInfo()
        subst = [a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic() and a.GetSymbol() == "C"
                 and any(n.GetSymbol() != "H" and not n.GetIsAromatic() for n in a.GetNeighbors())]
        out.append("aromatic ring C-H")
        if len(subst) == 1:
            d = len(Chem.GetShortestPath(mol, subst[0], c.GetIdx())) - 1
            pos = {1: "ortho", 2: "meta", 3: "para"}.get(d)
            if pos:
                out.append(f"{pos} ring C-H")
        hetero = [a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic() and a.GetSymbol() != "C"]
        if len(hetero) == 1 and not subst:
            d = len(Chem.GetShortestPath(mol, hetero[0], c.GetIdx())) - 1
            pos = {1: "2-position", 2: "3-position", 3: "4-position"}.get(d)
            if pos:
                out.append(f"ring C-H at the {pos}")
        return out
    hyb = c.GetHybridization()
    if any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(c).GetSymbol() == "O" for b in c.GetBonds()):
        return ["aldehydic C-H" if any(n.GetSymbol() == "C" for n in nbrs) else "formyl C-H"]
    if hyb == Chem.HybridizationType.SP:
        return ["acetylenic C-H"]
    if hyb == Chem.HybridizationType.SP2:
        return ["vinylic C-H"]
    if any(n.GetIsAromatic() for n in nbrs):
        return ["benzylic C-H"]
    if any(n.GetHybridization() == Chem.HybridizationType.SP2 and n.GetSymbol() == "C"
           and any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(n).GetSymbol() == "C" for b in n.GetBonds())
           for n in nbrs):
        return ["allylic C-H"]
    kind = {3: "methyl", 2: "methylene", 1: "methine"}.get(c.GetTotalNumHs(includeNeighbors=True), "")
    if not kind:
        return ["C-H"]
    out = [f"{kind} C-H"]
    for n in nbrs:
        g = _group(n)
        if g in ("oxygen", "nitrogen", "sulfur"):
            out.append(f"{kind} C-H on the {g}")
        elif g == "carbonyl":
            out.append(f"{kind} C-H next to the carbonyl")
    return out


def _descriptors(mol_h: Chem.Mol, bond: Chem.Bond) -> list[str]:
    a, b = bond.GetBeginAtom(), bond.GetEndAtom()
    if a.GetSymbol() == "H" and b.GetSymbol() == "H":
        return ["H-H"]
    if b.GetSymbol() == "H":
        a, b = b, a
    if a.GetSymbol() == "H":
        x = b
        if x.GetSymbol() == "C":
            return _ch_descriptors(x, mol_h)
        if x.GetDegree() == 1:  # diatomic hydride
            return [f"H-{x.GetSymbol()}"]
        return [f"{x.GetSymbol()}-H"]
    s1, s2 = sorted([a.GetSymbol(), b.GetSymbol()], key=lambda s: (s != "C", s))
    base = f"{s1}-{s2}"
    g1, g2 = _group(a), _group(b)
    if g1 == g2:
        return [base, f"central {base}"]
    g1, g2 = sorted([g1, g2])
    return [base, f"{base} bond between the {g1} group and the {g2}".replace("the ring group", "the ring")
            .replace("the oxygen group", "the oxygen").replace("the nitrogen group", "the nitrogen")
            .replace("the sulfur group", "the sulfur")]


def _bond_class(site: str) -> str:
    """'methyl C-H on the oxygen' -> 'C-H'; 'C-O bond between ...' -> 'C-O'."""
    if "C-H" in site:
        return "C-H"
    return next((tok for tok in site.replace("central ", "").split() if "-" in tok), site)


def tasks_for(name: str, smiles: str) -> list[dict]:
    mol = Chem.MolFromSmiles(smiles)
    mol_h = Chem.AddHs(mol)
    parent = canon(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    classes: dict[tuple, dict] = {}
    for bond in mol_h.GetBonds():
        if bond.GetBondTypeAsDouble() != 1 or bond.IsInRing():
            continue
        frags = homolysis(mol_h, bond.GetIdx())
        if not frags:
            continue
        key = tuple(sorted(canon(f) for f in frags))
        c = classes.setdefault(key, {"fragments": [{"smiles": canon(f), "charge": 0, "multiplicity": _multiplicity(f)}
                                                   for f in frags], "descs": []})
        c["descs"].append(_descriptors(mol_h, bond))
    # a description is usable when no other bond class shares it
    usage = defaultdict(set)
    for key, c in classes.items():
        for descs in c["descs"]:
            for d in descs:
                usage[d].add(key)
    out = []
    for key, c in classes.items():
        site = next((d for d in c["descs"][0] if usage[d] == {key}), None)
        if site is None:
            continue
        out.append({"id": f"{name.replace(' ', '_')}:{site.replace(' ', '_')}", "molecule": name, "smiles": smiles,
                    "formula": formula, "site": site, "bond_class": _bond_class(site),
                    "question": (QUESTION_BOND if " bond between " in site else QUESTION).format(
                        site=site, name=name, formula=formula),
                    "expected": {"parent": {"smiles": parent, "charge": 0, "multiplicity": 1},
                                 "fragments": sorted(c["fragments"], key=lambda f: f["smiles"])}})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--eval-frac", type=float, default=0.25)
    ap.add_argument("--holdout-class", default="C-Cl")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    tasks = [t for name, smi in MOLECULES for t in tasks_for(name, smi)]
    rng = random.Random(a.seed)
    mols = sorted({t["molecule"] for t in tasks})
    rng.shuffle(mols)
    eval_mols = set(mols[: int(len(mols) * a.eval_frac)])
    for t in tasks:
        t["split"] = "eval" if (t["molecule"] in eval_mols or t["bond_class"] == a.holdout_class) else "train"
        t["holdout_class"] = t["bond_class"] == a.holdout_class
    json.dump({"seed": a.seed, "holdout_class": a.holdout_class, "tasks": tasks}, open(a.out, "w"), indent=1)
    n = defaultdict(int)
    for t in tasks:
        n[t["split"]] += 1
    print(f"{len(tasks)} tasks from {len(mols)} molecules: {dict(n)}; held-out class {a.holdout_class}: "
          f"{sum(t['holdout_class'] for t in tasks)}")


if __name__ == "__main__":
    main()
