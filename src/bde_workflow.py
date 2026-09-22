import os
from fragmenter import break_bonds
from conformer import get_lowest_energy_conformer
from nwchem import mol_to_xyz, xyz_dict_to_str, write_nwchem_input
from parser import parse_nwchem_output
import pandas as pd
from storage.flowcept_store import flowcept_task

@flowcept_task
def write_species_files(smiles, name, outdir, charge, mult):
    mol_data = get_lowest_energy_conformer(smiles=smiles, num_confs=20)
    mol = mol_data["lowest_energy_conformer"]
    xyz_data = mol_to_xyz(mol=mol)
    xyz_str = xyz_dict_to_str(xyz_data)

    # Set the correct multiplicity for the molecule
    if name == "parent":
        charge, mult = 0, 1
    else:
        charge, mult = 0, 2

    nwchem_dict = write_nwchem_input(xyz_str=xyz_str, job_name=name, charge=charge, mult=mult)
    nwchem_input = nwchem_dict["input_text"]

    os.makedirs(outdir, exist_ok=True)

    xyz_path = os.path.join(outdir, f"{name}.xyz")
    nw_path = os.path.join(outdir, f"{name}.nw")
    
    with open(os.path.join(outdir, f"{name}.xyz"), 'w') as f:
        f.write(xyz_str)
    with open(os.path.join(outdir, f"{name}.nw"), 'w') as f:
        f.write(nwchem_input)

    return {
        "name": name,
        "xyz_path": xyz_path,
        "nw_path": nw_path,
        "charge": charge,
        "multiplicity": mult,
        "smiles": mol_data["smiles"],
        "energy": mol_data["energy"]
    }

@flowcept_task
def run_bde(smiles, outdir="bde_calc", generate_inputs=True, parse_outputs=True):
    print("Starting RUN_BDE exec")
    os.makedirs(outdir, exist_ok=True)
    all_species = {}
    bond_data = []

    if generate_inputs:
        all_species["parent"] = write_species_files(smiles=smiles, name="parent", outdir=outdir, charge=0, mult=1)
        fragments = break_bonds(smiles).get('fragments')
        for frags in fragments:
            label = frags["label"]
            smi1= frags["fragment1"]
            smi2 = frags["fragment2"]
            write_species_files(smiles=smi1, name=f"{label}_1", outdir=outdir, charge=0, mult=2)
            write_species_files(smiles=smi2, name=f"{label}_2", outdir=outdir, charge=0, mult=2)

    if parse_outputs:
        try:
            parent_data = parse_nwchem_output(filename=f"{outdir}/parent.out")
            e0 = parent_data["energy"]
            z0 = parent_data["zpe"]
            h0 = parent_data["enthalpy"]
            s0 = parent_data["entropy"]
        except Exception as exc:
            print(f"Parent output missing or failed: {exc}")
            return {"species": all_species, "bde_data": [], "error": f"parent output missing or failed: {exc}"}

        fragments = break_bonds(smiles=smiles).get("fragments")
        for frags in fragments:
            print("Going to run run_individual_bde")
            individual_bond_data = run_individual_bde(e0=e0, frags=frags, h0=h0, outdir=outdir, s0=s0, z0=z0)
            bond_data.append(individual_bond_data)

    # save results
    df = pd.DataFrame(bond_data)
    df.to_csv(os.path.join(outdir, "bde_results.csv"), index=False)
    print(df)

    return {
        "species": all_species,
        "bde_data": bond_data,
        "output_csv": os.path.join(outdir, "bde_results.csv")
    }

@flowcept_task
def run_individual_bde(e0, frags, h0, outdir, s0, z0):
    label = frags["label"]
    try:
        frag1_data = parse_nwchem_output(filename=f"{outdir}/{label}_1.out")
        frag2_data = parse_nwchem_output(filename=f"{outdir}/{label}_2.out")
        e1 = frag1_data["energy"]
        z1 = frag1_data["zpe"]
        h1 = frag1_data["enthalpy"]
        s1 = frag1_data["entropy"]

        e2 = frag2_data["energy"]
        z2 = frag2_data["zpe"]
        h2 = frag2_data["enthalpy"]
        s2 = frag2_data["entropy"]

        # Fail loudly (with a specific reason) rather than letting a missing
        # value silently blow up into the generic except-clause below, which
        # used to mask *every* kind of failure behind the same -1 sentinel.
        missing = [
            name for name, val in [
                ("frag1_zpe", z1), ("frag1_enthalpy", h1), ("frag1_entropy", s1),
                ("frag2_zpe", z2), ("frag2_enthalpy", h2), ("frag2_entropy", s2),
                ("parent_zpe", z0), ("parent_enthalpy", h0), ("parent_entropy", s0),
            ] if val is None
        ]
        if missing:
            raise ValueError(f"missing thermochemistry value(s): {', '.join(missing)}")

        bd_energy = (e1 + z1 + e2 + z2) - (e0 + z0)
        bd_enthalpy = (e1 + h1 + e2 + h2) - (e0 + h0)
        bd_free_energy = bd_enthalpy * 627.509 - (298.15 * (s2 + s1 - s0))  # G = H - T*S T=temperature
        return {
            "bond_id": label,
            "bd_energy": bd_energy * 627.509,
            "bd_enthalpy": bd_enthalpy * 627.509,
            "bd_free_energy": bd_free_energy,
            "status": "ok",
            "error": None,
        }
    except Exception as exc:
        # None (not -1) so a failed row is unmistakably missing data rather
        # than a plausible-looking-but-fake number, and status/error say why.
        print(f"BDE calc failed for bond {label}: {exc}")
        return {
            "bond_id": label,
            "bd_energy": None,
            "bd_enthalpy": None,
            "bd_free_energy": None,
            "status": "failed",
            "error": str(exc),
        }

