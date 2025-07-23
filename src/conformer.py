import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from flowcept import Flowcept, flowcept_task

@flowcept_task
def get_lowest_energy_conformer(smiles, num_confs=20):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: '{smiles}'")

    mol = Chem.AddHs(mol)
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=AllChem.ETKDG()))
    energies = []

    use_mmff = AllChem.MMFFHasAllMoleculeParams(mol)

    for conf_id in ids:
        if use_mmff:
            ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf_id)
        else:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
        ff.Minimize()
        energies.append(ff.CalcEnergy())

    best_index = int(np.argmin(energies))
    best_id = ids[int(np.argmin(energies))]
    best_energy = energies[best_index]
    best_conf = mol.GetConformer(best_id)

    # Keep only best conformer
    new_mol = Chem.Mol(mol)
    new_mol.RemoveAllConformers()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        pos = best_conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, pos)
    new_mol.AddConformer(conf, assignId=True)

    return {
        "lowest_energy_conformer": new_mol,
        "conformer_id": best_id,
        "energy": best_energy,
        "smiles": Chem.MolToSmiles(Chem.RemoveHs(new_mol))
    }

