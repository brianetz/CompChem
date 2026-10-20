"""
Deterministic half of the BDE plan critic (src/bde-plan-critic.md): the
checks the skill says code must decide before any judgment is applied.

critique_plan(spec, question) checks every species (SMILES parses,
electron parity vs multiplicity, radical centres vs multiplicity, formal
charge vs stated charge, ground-state multiplicity of isolated atoms),
conservation of atoms and charge from parent to fragments, and that the
fragments are a single-bond homolysis of the parent - of the bond the
question asks about, when it names one. A FAIL means the plan must be
revised; WARNs are for the LLM critic / a human to judge.

Not implemented yet (post-compute checks in the skill): <S^2> vs s(s+1),
spin-population location, and connectivity of the optimized geometry.
The BDE value band lives in skills/bde.review.
"""
from __future__ import annotations

from collections import Counter

from analysis.species_checks import electron_count, element_counts, formula, parse_smiles, spin_allowed

# Ground-state spin multiplicity of isolated neutral atoms.
_ATOM_GROUND_MULT = {
    "H": 2, "He": 1, "Li": 2, "Be": 1, "B": 2, "C": 3, "N": 4, "O": 3, "F": 2, "Ne": 1,
    "Na": 2, "Mg": 1, "Al": 2, "Si": 3, "P": 4, "S": 3, "Cl": 2, "Ar": 1, "Br": 2, "I": 2,
}
_BOND_SYMBOL = {1.0: "-", 1.5: ":", 2.0: "=", 3.0: "#"}


def _canonical(mol) -> str:
    """Canonical SMILES without explicit H or atom-map numbers, so 'C[O:1]',
    '[O]C' and 'C[O]' compare equal."""
    from rdkit import Chem

    mol = Chem.Mol(mol)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(Chem.RemoveHs(mol, sanitize=False))


def _item(check: str, evidence: str, fix: str, species: str | None = None) -> dict:
    return {"check": check, "species": species, "evidence": evidence, "fix": fix}


def enumerate_bde_candidates(parent_smiles: str) -> list[dict]:
    """Every symmetry-unique single bond of the parent (explicit H), with
    the two radical fragments its homolysis gives. Ring bonds are skipped:
    cutting one gives a diradical, not two fragments."""
    from rdkit import Chem

    parent = parse_smiles(parent_smiles)
    if parent is None:
        return []
    ranks = list(Chem.CanonicalRankAtoms(parent, breakTies=False))
    seen, out = set(), []
    for bond in parent.GetBonds():
        if bond.IsInRing():
            continue
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        sym_class = tuple(sorted((ranks[i], ranks[j])))
        if sym_class in seen:
            continue
        seen.add(sym_class)
        rw = Chem.RWMol(parent)
        for atom in rw.GetAtoms():
            atom.SetNoImplicit(True)
        rw.RemoveBond(i, j)
        for idx in (i, j):
            atom = rw.GetAtomWithIdx(idx)
            atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons() + 1)
        pieces = Chem.GetMolFrags(rw.GetMol(), asMols=True, sanitizeFrags=False)
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        if a1.GetAtomicNum() < a2.GetAtomicNum():
            a1, a2 = a2, a1
        out.append({
            "atom_indices": [i, j],
            "bond": f"{a1.GetSymbol()}{_BOND_SYMBOL.get(bond.GetBondTypeAsDouble(), '~')}{a2.GetSymbol()}",
            "fragments": sorted(_canonical(p) for p in pieces),
        })
    return out


def _check_species(sp: dict, name: str, rationale: str) -> tuple[list[dict], list[dict], object]:
    fails, warns = [], []
    smiles = sp.get("smiles")
    mol = parse_smiles(smiles)
    if mol is None:
        fails.append(_item(
            "smiles_parse", f"SMILES {smiles!r} is not valid",
            "a bracket atom holds ONE element plus its own hydrogens ('[CH3]', '[OH]'); "
            "bracket only the radical centre of a larger radical ('CC[O]', '[CH2]O'); "
            "write SMILES, not formulas ('O', not 'H2O')", name))
        return fails, warns, None

    charge, mult = int(sp.get("charge", 0)), int(sp.get("multiplicity", 1))
    n_e = electron_count(mol, charge)
    if not spin_allowed(n_e, mult):
        fails.append(_item(
            "parity", f"{smiles!r}, charge {charge}: {n_e} total electrons "
            f"({'odd' if n_e % 2 else 'even'}) cannot have multiplicity {mult}",
            "a radical fragment has one hydrogen fewer than the parent at the broken bond; "
            "its SMILES must show that radical centre", name))

    formal = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    if formal != charge:
        fails.append(_item(
            "formal_charge", f"{smiles!r} carries formal charge {formal:+d} but the stated charge is {charge}",
            "make the SMILES charges and the stated charge agree", name))

    if mol.GetNumAtoms() == 1:  # explicit H included, so this is a lone atom
        el = mol.GetAtomWithIdx(0).GetSymbol()
        ground = _ATOM_GROUND_MULT.get(el)
        if charge == 0 and ground and mult != ground:
            fails.append(_item(
                "atom_ground_state", f"isolated {el} atom given multiplicity {mult}",
                f"the ground state of a neutral {el} atom has multiplicity {ground}", name))
    else:
        radicals = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
        if radicals == 1 and mult != 2:
            fails.append(_item(
                "radical_centres", f"{smiles!r} has one radical centre but multiplicity {mult}",
                "one unpaired electron means a doublet (multiplicity 2)", name))
        elif radicals >= 2:
            allowed = range(radicals % 2 + 1, radicals + 2, 2)
            if mult not in allowed:
                fails.append(_item(
                    "radical_centres", f"{smiles!r} has {radicals} radical electrons but multiplicity {mult}",
                    f"use one of multiplicities {list(allowed)}", name))
            elif not any(w in rationale.lower() for w in ("singlet", "triplet")):
                warns.append(_item(
                    "radical_centres", f"{smiles!r} has {radicals} radical electrons",
                    "justify the singlet/triplet choice and say which is the ground state", name))
        elif radicals == 0 and mult > 1:
            warns.append(_item(
                "radical_centres", f"{smiles!r} shows no radical centre but has multiplicity {mult}",
                "check that this open-shell state is intended (e.g. triplet O2)", name))
    return fails, warns, mol


def critique_plan(spec: dict, question: str | None = None) -> dict:
    """{'verdict': PASS|WARN|FAIL, 'fails': [...], 'warns': [...],
    'cleaved_bond': ..., 'candidates': [...]} for one BDE spec."""
    from storage.rag import bond_types_in

    rationale = spec.get("rationale") or ""
    parent = spec.get("parent") or {}
    fragments = spec.get("fragments") or []
    fails, warns = [], []

    f, w, parent_mol = _check_species(parent, "parent", rationale)
    fails += f
    warns += w
    if parent_mol is not None and int(parent.get("multiplicity", 1)) != 1:
        warns.append(_item("parent_state", f"parent multiplicity {parent.get('multiplicity')}",
                           "a stable closed-shell parent is a singlet; justify otherwise", "parent"))
    frag_mols = []
    for i, frag in enumerate(fragments):
        f, w, mol = _check_species(frag, f"fragment {frag.get('label') or i}", rationale)
        fails += f
        warns += w
        frag_mols.append(mol)

    report = {"cleaved_bond": None, "candidates": []}
    if parent_mol is not None and all(m is not None for m in frag_mols):
        want, got = element_counts(parent_mol), sum((element_counts(m) for m in frag_mols), Counter())
        if want != got:
            fails.append(_item(
                "atom_balance", f"fragments add up to {formula(got)} but the parent is {formula(want)}",
                "the fragments must contain exactly the parent's atoms, split at the one bond that breaks"))
        charge_sum = sum(int(fr.get("charge", 0)) for fr in fragments)
        if charge_sum != int(parent.get("charge", 0)):
            fails.append(_item(
                "charge_balance", f"fragment charges sum to {charge_sum}, parent charge is {parent.get('charge', 0)}",
                "homolysis keeps every fragment's charge summing to the parent's"))

        candidates = enumerate_bde_candidates(parent["smiles"])
        report["candidates"] = candidates
        planned = sorted(_canonical(m) for m in frag_mols)
        matches = [c for c in candidates if c["fragments"] == planned]
        wanted = bond_types_in(question or "")
        if not fails and not matches:
            hint = [c for c in candidates if not wanted or c["bond"] in wanted]
            fails.append(_item(
                "homolysis", f"fragments {planned} are not the homolysis of any single bond of "
                f"{parent['smiles']!r} (wrong bond, radical on the wrong atom, or not one bond breaking)",
                "symmetry-unique homolyses" + (f" of {'/'.join(wanted)} bonds" if wanted else "")
                + ": " + "; ".join(" + ".join(c["fragments"]) for c in hint[:6])))
        elif matches:
            report["cleaved_bond"] = matches[0]
            if wanted and matches[0]["bond"] not in wanted:
                options = [" + ".join(c["fragments"]) for c in candidates if c["bond"] in wanted]
                fails.append(_item(
                    "requested_bond", f"these fragments break a {matches[0]['bond']} bond",
                    f"the question asks about {' or '.join(wanted)}"
                    + (f"; its symmetry-unique homolyses: {'; '.join(options[:6])}" if options else "")))
            same_type = [c for c in candidates if c["bond"] == matches[0]["bond"]]
            if len(same_type) > 1 and not spec.get("bond_cleaved"):
                warns.append(_item(
                    "bond_choice", f"the parent has {len(same_type)} symmetry-distinct "
                    f"{matches[0]['bond']} bonds",
                    "confirm the planned one is the bond the question means"))

    if any(int(s.get("charge", 0)) < 0 for s in [parent, *fragments]):
        basis = str(spec.get("basis", "")).lower()
        if "+" not in basis and "aug" not in basis and not basis.startswith("ma-"):
            warns.append(_item("method", f"an anion with basis {spec.get('basis')!r}",
                               "use diffuse functions (e.g. 6-31+g*, aug-cc-pvdz, ma-def2-svp)"))

    verdict = "FAIL" if fails else ("WARN" if warns else "PASS")
    return {"verdict": verdict, "fails": fails, "warns": warns, **report}


def problems_as_text(items: list[dict]) -> list[str]:
    """Validator items as one-line instructions for the planner's retry."""
    return [
        (f"{it['species']}: " if it.get("species") else "") + f"{it['evidence']}. Fix: {it['fix']}"
        for it in items
    ]
