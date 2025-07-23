from rdkit import Chem
from flowcept import Flowcept, flowcept_task

@flowcept_task
def break_bonds(smiles):
    """
    Returns list of dictionaries: 
    {
        'label': bond_label,
        'fragment1': frag1_smi,
        'fragment2': frag2_smi
    }
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    frags = []

    for bond in mol.GetBonds():
        if bond.IsInRing():
            continue  # skip ring bonds

        if bond.GetBondTypeAsDouble() > 1.9999:
            continue # skip double/triple bonds
        try:
            frag = break_individual_bond(bond, mol)
        except:
            frag = {}
        frags.append(frag)

    return {"fragments": frags}

@flowcept_task
def break_individual_bond(bond, mol):
    atom1 = bond.GetBeginAtom()
    atom2 = bond.GetEndAtom()
    idx1 = atom1.GetIdx()
    idx2 = atom2.GetIdx()
    bond_idx = bond.GetIdx()
    label = f"{atom1.GetSymbol()}-{atom2.GetSymbol()}_{bond_idx}"
    rwmol = Chem.RWMol(mol)
    rwmol.RemoveBond(idx1, idx2)
    # Set radicals at both ends
    for idx in [idx1, idx2]:
        atom = rwmol.GetAtomWithIdx(idx)
        h_count = atom.GetNumExplicitHs()
        if h_count > 0:
            atom.SetNumExplicitHs(h_count - 1)
        atom.SetNoImplicit(True)
        atom.SetNumRadicalElectrons(1)
    Chem.SanitizeMol(rwmol)
    fragments = Chem.GetMolFrags(rwmol, asMols=True)
    if len(fragments) == 2:
        smi1 = Chem.MolToSmiles(fragments[0], isomericSmiles=True)
        smi2 = Chem.MolToSmiles(fragments[1], isomericSmiles=True)
        return {
            'label': label,
            'fragment1': smi1,
            'fragment2': smi2
        }
    else:
        return {
            'label': label,
            'fragment1': "not_enough_fragments",
            'fragment2': "not_enough_fragments"
        }
