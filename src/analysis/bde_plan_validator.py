"""Deterministic validation of LLM-generated BDE campaign specs.

Intended home in the repo: ``src/analysis/bde_plan_validator.py``.

Two layers:

* ``critique_plan(spec)``: runs BEFORE any compute. Re-derives the electronic
  structure bookkeeping and the bond cleavage independently of the LLM and
  compares them with what the planner wrote (species table, rationale text).
* ``check_spin`` / ``check_connectivity`` / ``check_bde_value``: run AFTER the
  calculations, to confirm the radical ended up where the plan said it would.

Everything that needs RDKit imports it lazily, so the pure-Python pieces
(electron counting, rationale number checks, spin and value checks) work
without it. Atom indices throughout refer to ``Chem.AddHs(Chem.MolFromSmiles(s))``
ordering (heavy atoms in SMILES order, then hydrogens). Verify that this matches
the ordering produced by ``engines/rdkit_to_ase.py`` before relying on the
spin-localization check.

Status: the pure-Python parts are covered by the unit tests shipped alongside
this file. The RDKit-dependent parts were written against the documented RDKit
API but have NOT been run; run the test file in the ``compchem-agent`` env first.
"""
from __future__ import annotations

import json
import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from typing import Any, Iterable, Optional, Sequence

try:  # RDKit is optional at import time
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    HAVE_RDKIT = True
except Exception:  # ImportError, or a broken install
    Chem = None  # type: ignore[assignment]
    rdDetermineBonds = None  # type: ignore[assignment]
    HAVE_RDKIT = False

PASS, WARN, FAIL = "PASS", "WARN", "FAIL"

# symbol -> (atomic number, main-group valence electrons, neutral-atom ground-state multiplicity)
_ELEMENTS: dict[str, tuple[int, int, int]] = {
    "H": (1, 1, 2), "He": (2, 2, 1),
    "Li": (3, 1, 2), "Be": (4, 2, 1), "B": (5, 3, 2), "C": (6, 4, 3),
    "N": (7, 5, 4), "O": (8, 6, 3), "F": (9, 7, 2), "Ne": (10, 8, 1),
    "Na": (11, 1, 2), "Mg": (12, 2, 1), "Al": (13, 3, 2), "Si": (14, 4, 3),
    "P": (15, 5, 4), "S": (16, 6, 3), "Cl": (17, 7, 2), "Ar": (18, 8, 1),
    "Se": (34, 6, 3), "Br": (35, 7, 2), "I": (53, 7, 2),
}


# --------------------------------------------------------------------------
# Result containers
# --------------------------------------------------------------------------
@dataclass
class Check:
    name: str
    status: str  # PASS | WARN | FAIL
    message: str


@dataclass
class PlanReport:
    checks: list[Check] = field(default_factory=list)

    def add(self, name: str, status: str, message: str) -> None:
        self.checks.append(Check(name, status, message))

    @property
    def verdict(self) -> str:
        statuses = {c.status for c in self.checks}
        if FAIL in statuses:
            return FAIL
        if WARN in statuses:
            return WARN
        return PASS

    def to_qc(self) -> dict[str, Any]:
        """Same shape as the repo's qc_common.classify() output."""
        notes = [f"{c.name}: {c.message}" for c in self.checks if c.status != PASS]
        return {"qc_flag": self.verdict, "qc_notes": notes}

    def feedback(self) -> str:
        """Text suitable for a bounded corrective turn back to the planner."""
        bad = [c for c in self.checks if c.status != PASS]
        if not bad:
            return "No problems found."
        return "\n".join(f"[{c.status}] {c.name}: {c.message}" for c in bad)


# --------------------------------------------------------------------------
# Pure-Python electron bookkeeping
# --------------------------------------------------------------------------
def electron_counts(symbols: Sequence[str], charge: int) -> dict[str, Any]:
    """Total and valence electron counts for a species given by element symbols."""
    unsupported = sorted({s for s in symbols if s not in _ELEMENTS})
    if unsupported:
        return {"total": None, "valence": None, "unsupported": unsupported}
    total = sum(_ELEMENTS[s][0] for s in symbols) - charge
    valence = sum(_ELEMENTS[s][1] for s in symbols) - charge
    return {"total": total, "valence": valence, "unsupported": []}


def parity_ok(total_electrons: int, multiplicity: int) -> bool:
    """Multiplicity is consistent with the electron count (odd N <=> even multiplicity)."""
    unpaired = multiplicity - 1
    return 0 <= unpaired <= total_electrons and (total_electrons - unpaired) % 2 == 0


def ground_state_multiplicity(symbol: str) -> Optional[int]:
    entry = _ELEMENTS.get(symbol)
    return entry[2] if entry else None


_CLAIM_RE = re.compile(r"\b(\d+)\s+(valence|total|core)?\s*electrons?\b", re.IGNORECASE)


def check_electron_claims(
    text: str, total_counts: Iterable[int], valence_counts: Iterable[int]
) -> list[str]:
    """Compare 'N (valence|total) electrons' claims in prose with computed counts.

    Counts are pooled across species, so this catches wrong or mislabeled numbers
    but cannot attribute a claim to a particular species.
    """
    totals, valences = set(total_counts), set(valence_counts)
    issues = []
    for m in _CLAIM_RE.finditer(text or ""):
        n, qualifier = int(m.group(1)), (m.group(2) or "").lower()
        if qualifier == "valence":
            if n not in valences:
                hint = " (that matches a TOTAL electron count; the label is wrong)" if n in totals else ""
                issues.append(f"'{m.group(0)}': no species has {n} valence electrons{hint}")
        elif n not in totals and n not in valences:
            issues.append(f"'{m.group(0)}': no species has {n} electrons")
    return issues


# --------------------------------------------------------------------------
# RDKit helpers
# --------------------------------------------------------------------------
def _require_rdkit() -> None:
    if not HAVE_RDKIT:
        raise ImportError("RDKit is required for this check (conda env 'compchem-agent').")


def _hmol(smiles: str):
    _require_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit cannot parse SMILES {smiles!r}")
    return Chem.AddHs(mol)


def _key_from_hmol(hmol) -> str:
    """Canonical key; both planned and derived fragments go through this same path."""
    key = Chem.MolToSmiles(Chem.RemoveHs(hmol))
    return key or Chem.MolToSmiles(hmol)


def canonical_key(smiles: str) -> str:
    return _key_from_hmol(_hmol(smiles))


def _homolyse(hmol, a: int, b: int):
    """Homolytically cleave single bond a-b; return fragment molecules (with Hs)."""
    bond = hmol.GetBondBetweenAtoms(a, b)
    if bond is None:
        raise ValueError(f"no bond between atoms {a} and {b}")
    if bond.GetBondType() != Chem.BondType.SINGLE:
        raise ValueError("only single bonds are supported for homolysis")
    rw = Chem.RWMol(hmol)
    rw.RemoveBond(a, b)
    for i in (a, b):
        atom = rw.GetAtomWithIdx(i)
        atom.SetNoImplicit(True)  # stop RDKit from refilling the lost valence with H
        atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons() + 1)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)


def homolysis_fragment_keys(parent_smiles: str, a: int, b: int) -> list[str]:
    """Canonical keys of the pieces produced by homolysing bond a-b.

    One key back (instead of two) means the bond is in a ring: the result is a
    single diradical, not two fragments.
    """
    return [_key_from_hmol(f) for f in _homolyse(_hmol(parent_smiles), a, b)]


def _bond_label(hmol, a: int, b: int) -> str:
    return f"{hmol.GetAtomWithIdx(a).GetSymbol()}{a}-{hmol.GetAtomWithIdx(b).GetSymbol()}{b}"


def infer_cleaved_bonds(parent_smiles: str, planned_fragments: Sequence[str]) -> list[dict]:
    """Every single bond whose homolysis reproduces the planned fragments exactly.

    Radical location is part of the comparison (canonical SMILES encodes which
    atom carries the radical), so a plan that puts the radical on the wrong atom
    returns no hits. Symmetry-equivalent bonds share a 'class' value.
    """
    hmol = _hmol(parent_smiles)
    target = Counter(canonical_key(s) for s in planned_fragments)
    ranks = list(Chem.CanonicalRankAtoms(hmol, breakTies=False))
    hits = []
    for bond in hmol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        try:
            keys = Counter(_key_from_hmol(f) for f in _homolyse(hmol, a, b))
        except Exception:
            continue
        if keys == target:
            hits.append({
                "atoms": (a, b),
                "class": tuple(sorted((ranks[a], ranks[b]))),
                "label": _bond_label(hmol, a, b),
            })
    return hits


def enumerate_bde_candidates(
    parent_smiles: str, elements: Optional[tuple[str, str]] = None
) -> list[dict]:
    """One representative per symmetry-unique single bond (optionally filtered by element pair).

    Use for 'weakest X-Y bond' questions: the plan should cover every entry.
    """
    hmol = _hmol(parent_smiles)
    ranks = list(Chem.CanonicalRankAtoms(hmol, breakTies=False))
    wanted = tuple(sorted(elements)) if elements else None
    groups: dict[tuple, dict] = {}
    for bond in hmol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        pair = tuple(sorted((hmol.GetAtomWithIdx(a).GetSymbol(), hmol.GetAtomWithIdx(b).GetSymbol())))
        if wanted and pair != wanted:
            continue
        cls = tuple(sorted((ranks[a], ranks[b])))
        if cls in groups:
            groups[cls]["n_equivalent"] += 1
            continue
        try:
            keys = [_key_from_hmol(f) for f in _homolyse(hmol, a, b)]
        except Exception:
            continue
        groups[cls] = {
            "atoms": (a, b),
            "label": _bond_label(hmol, a, b),
            "fragments": keys,
            "in_ring": len(keys) == 1,
            "n_equivalent": 1,
        }
    return sorted(groups.values(), key=lambda g: g["label"])


def expected_radical_atoms(smiles: str) -> list[int]:
    """Indices (AddHs ordering) of atoms carrying radical electrons in this species."""
    hmol = _hmol(smiles)
    return [a.GetIdx() for a in hmol.GetAtoms() if a.GetNumRadicalElectrons() > 0]


# --------------------------------------------------------------------------
# Pre-execution plan critique
# --------------------------------------------------------------------------
def _species(obj: Any) -> Optional[dict]:
    if not isinstance(obj, dict) or not isinstance(obj.get("smiles"), str):
        return None
    return {
        "smiles": obj["smiles"],
        "charge": int(obj.get("charge", 0) or 0),
        "multiplicity": int(obj.get("multiplicity", 1) or 1),
        "raw": obj,
    }


def _check_species(report: PlanReport, label: str, sp: dict) -> Optional[dict]:
    """Electron bookkeeping for one species; returns counts for later prose checks."""
    smiles, charge, mult = sp["smiles"], sp["charge"], sp["multiplicity"]
    try:
        hmol = _hmol(smiles)
    except ValueError as exc:
        report.add(f"{label}_smiles", FAIL, str(exc))
        return None
    symbols = [a.GetSymbol() for a in hmol.GetAtoms()]
    counts = electron_counts(symbols, charge)
    if counts["total"] is None:
        report.add(f"{label}_elements", WARN, f"unsupported elements {counts['unsupported']}; electron checks skipped")
        return None

    formal = Chem.GetFormalCharge(hmol)
    if formal != charge:
        report.add(f"{label}_charge", FAIL, f"{smiles} has formal charge {formal:+d} but the spec says {charge:+d}")

    total = counts["total"]
    if not parity_ok(total, mult):
        report.add(
            f"{label}_electron_parity", FAIL,
            f"{smiles}: {total} electrons is incompatible with multiplicity {mult}"
            " (odd electron count needs an even multiplicity, and vice versa)",
        )
    else:
        report.add(f"{label}_electron_parity", PASS, f"{smiles}: {total} electrons, multiplicity {mult}")

    n_rad = sum(a.GetNumRadicalElectrons() for a in hmol.GetAtoms())
    if n_rad == 1 and mult != 2:
        report.add(f"{label}_radical_spin", FAIL, f"{smiles} is a monoradical, so multiplicity must be 2 (got {mult})")
    elif n_rad == 0 and mult != 1:
        report.add(f"{label}_radical_spin", WARN,
                   f"{smiles} has no radical centres but multiplicity is {mult}; valid only for known high-spin ground states (e.g. O2, CH2)")
    elif n_rad == 2 and mult in (1, 3):
        report.add(f"{label}_radical_spin", WARN,
                   f"{smiles} is a diradical; the rationale must justify singlet vs triplet (got {mult})")
    elif n_rad >= 3:
        report.add(f"{label}_radical_spin", WARN, f"{smiles} has {n_rad} radical electrons; confirm the spin state")
    else:
        report.add(f"{label}_radical_spin", PASS, f"{n_rad} radical electron(s), multiplicity {mult}")

    if len(symbols) == 1 and charge == 0:
        expected = ground_state_multiplicity(symbols[0])
        if expected is not None and mult != expected:
            report.add(f"{label}_atom_ground_state", WARN,
                       f"isolated {symbols[0]} atom: ground-state multiplicity is {expected}, plan uses {mult}")
    return {"symbols": symbols, **counts}


def _norm(text: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (text or "").lower())


def _check_rationale(report: PlanReport, spec: dict, counts: list[dict], n_fragments: int) -> None:
    text = spec.get("rationale") or ""
    if not text.strip():
        report.add("rationale_present", FAIL, "no rationale text")
        return
    low = text.lower()

    if counts:  # without computed counts every claim would look wrong
        issues = check_electron_claims(
            text, [c["total"] for c in counts], [c["valence"] for c in counts]
        )
        for issue in issues:
            report.add("rationale_electron_numbers", FAIL, issue)

    for key, label in (("basis", "basis set"), ("xc", "functional")):
        value = spec.get(key)
        if value and _norm(str(value)) not in _norm(text):
            report.add(f"rationale_{key}", WARN, f"rationale never mentions the planned {label} ({value})")

    rubric = [
        ("rationale_names_bond", bool(re.search(r"\b[A-Z][a-z]?\s*[-–]\s*[A-Z][a-z]?\b", text)) or "bond" in low,
         "does not say which bond breaks"),
        ("rationale_spin_states", any(w in low for w in ("singlet", "doublet", "triplet", "multiplicity")),
         "does not state spin states / multiplicities"),
        ("rationale_parent_state", any(w in low for w in ("closed-shell", "closed shell", "singlet")),
         "does not state the parent's spin state"),
        ("rationale_formula", bool(re.search(r"bde\s*=|e\(", low)),
         "does not write out the BDE formula (sum of fragment energies minus parent)"),
        ("rationale_units", "kcal" in low, "does not state the units (kcal/mol)"),
        ("rationale_target_quantity", any(w in low for w in ("zpe", "zero-point", "thermal", "enthalpy", "electronic")),
         "does not say whether the BDE is electronic-only or includes ZPE/thermal corrections"),
        ("rationale_workflow", "xtb" in low and ("single-point" in low or "single point" in low),
         "does not describe the xtb-then-DFT-single-point stages"),
    ]
    if n_fragments == 2:
        rubric.append(("rationale_contrast", any(w in low for w in ("atomization", "atomisation")),
                       "does not contrast single-bond homolysis with full atomization"))
    for name, ok, msg in rubric:
        if not ok:
            report.add(name, WARN, msg)


def critique_plan(spec: dict) -> PlanReport:
    """Validate a BDE campaign spec before any compute runs."""
    report = PlanReport()

    if spec.get("analysis_type") != "bde":
        report.add("analysis_type", FAIL, f"expected 'bde', got {spec.get('analysis_type')!r}")

    parent = _species(spec.get("parent"))
    fragments_raw = spec.get("fragments")
    fragments = [_species(f) for f in fragments_raw] if isinstance(fragments_raw, list) else []
    if parent is None:
        report.add("parent_shape", FAIL, "parent must be an object with smiles/charge/multiplicity")
    if not fragments or any(f is None for f in fragments):
        report.add("fragments_shape", FAIL, "fragments must be a non-empty list of objects with smiles/charge/multiplicity")
    if parent is None or not fragments or any(f is None for f in fragments):
        return report

    for key in ("basis", "xc"):
        if not spec.get(key):
            report.add(f"spec_{key}", FAIL, f"spec has no {key}")
    for sp in [parent, *fragments]:
        for key in ("basis", "xc"):
            if sp["raw"].get(key) and sp["raw"][key] != spec.get(key):
                report.add("method_consistency", FAIL,
                           f"{sp['smiles']} overrides {key} ({sp['raw'][key]} vs {spec.get(key)}); every species must use the same level")

    if len(fragments) != 2:
        report.add("fragment_count", WARN,
                   f"{len(fragments)} fragments: not a simple single-bond homolysis; the rationale must say what is being computed")

    counts: list[dict] = []
    if HAVE_RDKIT:
        for label, sp in [("parent", parent)] + [(f"fragment{i + 1}", f) for i, f in enumerate(fragments)]:
            info = _check_species(report, label, sp)
            if info:
                counts.append(info)

        # atom and charge conservation
        try:
            parent_atoms = Counter(a.GetSymbol() for a in _hmol(parent["smiles"]).GetAtoms())
            frag_atoms: Counter = Counter()
            for f in fragments:
                frag_atoms.update(a.GetSymbol() for a in _hmol(f["smiles"]).GetAtoms())
            if parent_atoms != frag_atoms:
                report.add("atom_conservation", FAIL,
                           f"parent has {dict(parent_atoms)} but fragments sum to {dict(frag_atoms)}")
            else:
                report.add("atom_conservation", PASS, "fragments contain exactly the parent's atoms")
            if parent["charge"] != sum(f["charge"] for f in fragments):
                report.add("charge_conservation", FAIL, "fragment charges do not sum to the parent charge")
        except ValueError:
            pass  # already reported as a SMILES failure above

        # independent derivation of the cleavage
        if len(fragments) == 2:
            _check_cleavage(report, spec, parent, fragments)
    else:
        report.add("rdkit_unavailable", WARN, "RDKit not installed; structural and electron checks skipped")

    _check_rationale(report, spec, counts, len(fragments))
    return report


def _check_cleavage(report: PlanReport, spec: dict, parent: dict, fragments: list[dict]) -> None:
    planned = [f["smiles"] for f in fragments]
    declared = spec.get("bond_cleaved")
    if isinstance(declared, dict) and declared.get("atom_indices"):
        a, b = declared["atom_indices"]
        try:
            got = Counter(homolysis_fragment_keys(parent["smiles"], a, b))
        except Exception as exc:
            report.add("bond_cleaved", FAIL, f"declared bond {a}-{b} cannot be homolysed: {exc}")
            return
        want = Counter(canonical_key(s) for s in planned)
        if got == want:
            report.add("bond_cleaved", PASS, f"declared bond {a}-{b} reproduces the planned fragments")
        else:
            report.add("bond_cleaved", FAIL,
                       f"declared bond {a}-{b} gives {sorted(got)} but the plan lists {sorted(want)}")
        return

    hits = infer_cleaved_bonds(parent["smiles"], planned)
    classes = {h["class"] for h in hits}
    if not hits:
        ring_note = ""
        try:
            if any(len(c["fragments"]) == 1 for c in enumerate_bde_candidates(parent["smiles"])):
                ring_note = " (the parent has ring bonds; cleaving one gives a single diradical, not two fragments)"
        except Exception:
            pass
        report.add("bond_cleavage", FAIL,
                   "no single-bond homolysis of the parent reproduces the planned fragments, radical positions included"
                   f"{ring_note}")
    elif len(classes) > 1:
        report.add("bond_cleavage", WARN,
                   "planned fragments match more than one bond class: " + ", ".join(h["label"] for h in hits)
                   + "; have the planner declare the bond")
    else:
        report.add("bond_cleavage", PASS,
                   f"fragments (radical positions included) match cleaving {hits[0]['label']}"
                   f" ({len(hits)} symmetry-equivalent bond(s))")


# --------------------------------------------------------------------------
# Post-calculation checks
# --------------------------------------------------------------------------
def check_spin(
    label: str,
    multiplicity: int,
    s2: Optional[float],
    atom_spin: Optional[Sequence[float]] = None,
    expected_atoms: Optional[Sequence[int]] = None,
) -> list[Check]:
    """<S^2> contamination and radical-location checks.

    s2: <S^2> from the UKS run (PySCF: mf.spin_square()[0]); None for restricted runs.
    atom_spin: per-atom Mulliken spin populations (PySCF: mf.mulliken_spin_pop()[1]).
    expected_atoms: where the plan put the radical (see expected_radical_atoms).
    Thresholds are conventional, not universal; adjust to taste.
    """
    out: list[Check] = []
    spin = (multiplicity - 1) / 2
    ideal = spin * (spin + 1)
    if s2 is not None:
        dev = abs(s2 - ideal)
        warn_tol = max(0.10 * ideal, 0.05)
        fail_tol = max(0.25 * ideal, 0.15)
        if dev > fail_tol:
            out.append(Check(f"{label}_spin_contamination", FAIL, f"<S^2>={s2:.3f}, ideal {ideal:.3f}"))
        elif dev > warn_tol:
            out.append(Check(f"{label}_spin_contamination", WARN, f"<S^2>={s2:.3f}, ideal {ideal:.3f}"))
        else:
            out.append(Check(f"{label}_spin_contamination", PASS, f"<S^2>={s2:.3f}, ideal {ideal:.3f}"))

    if atom_spin is not None:
        net = sum(atom_spin)
        if abs(net - (multiplicity - 1)) > 0.05:
            out.append(Check(f"{label}_spin_sum", WARN,
                             f"spin populations sum to {net:.2f}, expected {multiplicity - 1}"))
        if expected_atoms:
            on_expected = sum(atom_spin[i] for i in expected_atoms)
            heaviest = max(range(len(atom_spin)), key=lambda i: abs(atom_spin[i]))
            if on_expected >= 0.6 and heaviest in expected_atoms:
                out.append(Check(f"{label}_radical_location", PASS,
                                 f"{on_expected:.2f} of the unpaired spin sits on the planned radical atom(s)"))
            elif on_expected >= 0.3:
                out.append(Check(f"{label}_radical_location", WARN,
                                 f"only {on_expected:.2f} of the spin is on the planned atom(s); delocalized radical or wrong site"))
            else:
                out.append(Check(f"{label}_radical_location", FAIL,
                                 f"only {on_expected:.2f} of the spin is on the planned atom(s); largest spin is on atom {heaviest}"))
    return out


def _graph_key(atomic_nums: Sequence[int], edges: Iterable[tuple[int, int]]) -> str:
    """Element-labelled connectivity fingerprint, independent of bond orders."""
    rw = Chem.RWMol()
    for z in atomic_nums:
        atom = Chem.Atom(int(z))
        atom.SetNoImplicit(True)
        rw.AddAtom(atom)
    for i, j in edges:
        if i != j and rw.GetBondBetweenAtoms(i, j) is None:
            rw.AddBond(i, j, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    return Chem.MolToSmiles(mol)


def check_connectivity(
    label: str,
    symbols: Sequence[str],
    positions: Sequence[Sequence[float]],
    planned_smiles: str,
    severity: str = WARN,
) -> Check:
    """Did the optimized geometry keep the connectivity of the planned species?

    Catches rearrangements (1,2-H shifts, ring closures) and further dissociation
    during the xtb optimization. Bonds are perceived from covalent radii, so a very
    stretched bond can trigger a false alarm; hence WARN by default.
    """
    _require_rdkit()
    block = "\n".join(
        [str(len(symbols)), label]
        + [f"{s} {x:.6f} {y:.6f} {z:.6f}" for s, (x, y, z) in zip(symbols, positions)]
    ) + "\n"
    geom = Chem.MolFromXYZBlock(block)
    if geom is None:
        return Check(f"{label}_connectivity", WARN, "could not read the geometry")
    rdDetermineBonds.DetermineConnectivity(geom)
    found = _graph_key(
        [a.GetAtomicNum() for a in geom.GetAtoms()],
        [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in geom.GetBonds()],
    )
    planned = _hmol(planned_smiles)
    expected = _graph_key(
        [a.GetAtomicNum() for a in planned.GetAtoms()],
        [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in planned.GetBonds()],
    )
    if found == expected:
        return Check(f"{label}_connectivity", PASS, "optimized geometry matches the planned connectivity")
    return Check(f"{label}_connectivity", severity,
                 f"perceived connectivity {found} differs from planned {expected}; possible rearrangement or dissociation")


def check_bde_value(bde_kcal: Optional[float]) -> Check:
    """Coarse sanity band for a single-bond BDE in kcal/mol."""
    if bde_kcal is None:
        return Check("bde_value", FAIL, "no BDE was computed")
    if bde_kcal < 0:
        return Check("bde_value", FAIL, f"negative BDE ({bde_kcal:.1f} kcal/mol): fragments lower than parent")
    if bde_kcal < 20 or bde_kcal > 150:
        return Check("bde_value", WARN, f"{bde_kcal:.1f} kcal/mol is outside the usual 20-150 range for a single bond")
    return Check("bde_value", PASS, f"{bde_kcal:.1f} kcal/mol")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python -m analysis.bde_plan_validator spec.json")
    with open(sys.argv[1]) as fh:
        rep = critique_plan(json.load(fh))
    for c in rep.checks:
        print(f"[{c.status}] {c.name}: {c.message}")
    print(f"\nVERDICT: {rep.verdict}")
