"""
Bridge between RDKit (used for conformer search/fragmentation) and ASE (used
to drive the alternate QC engines in this package) - the common currency
between the cheminformatics layer and the QC-engine layer.
"""
from __future__ import annotations

from ase import Atoms
from rdkit import Chem
from rdkit.Chem import AllChem


def rdkit_mol_to_ase_atoms(mol) -> Atoms:
    """Convert an RDKit Mol with an embedded conformer (e.g. the
    'lowest_energy_conformer' from conformer.get_lowest_energy_conformer)
    into an ase.Atoms object."""
    conf = mol.GetConformer()
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    positions = [tuple(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    return Atoms(symbols=symbols, positions=positions)


def smiles_to_ase_atoms(smiles: str) -> Atoms:
    """Build an ase.Atoms object directly from a SMILES string, with a
    single embedded conformer rather than conformer.get_lowest_energy_conformer's
    full multi-conformer MMFF/UFF force-field search.

    Added in Phase 5 (§14 architecture plan) for BDE campaign species -
    dissociation fragments (radicals, lone atoms like '[H]') that a
    100-conformer force-field search isn't designed for and that xtb will
    geometry-optimize from this starting guess anyway, so a single cheap
    embed is all that's needed. Handles the single-atom case explicitly
    (RDKit's distance-geometry embedder has nothing to do with zero bonds).
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        # Some radical SMILES (e.g. explicit open-valence atoms) don't
        # survive RDKit's default sanitization; retry with a relaxed pass.
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            raise ValueError(f"Could not parse SMILES: {smiles!r}")
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)

    mol = Chem.AddHs(mol)

    if mol.GetNumAtoms() == 1:
        conf = Chem.Conformer(1)
        conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
        mol.AddConformer(conf, assignId=True)
    else:
        AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE)

    return rdkit_mol_to_ase_atoms(mol)
