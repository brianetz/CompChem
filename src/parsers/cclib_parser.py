from __future__ import annotations

from pathlib import Path
from typing import Optional, Dict, Any

from cclib.io import ccread
from ase import Atoms
from ase.io import write as ase_write

from storage.records import CalculationRecord
from parsers.nwchem_fallback import (
    read_text,
    extract_method_basis,
    extract_charge_multiplicity,
    extract_success,
    extract_job_metadata,
    extract_final_energy_hartree,
)
from utils.paths import ensure_dir
from storage.flowcept_store import flowcept_task


def count_imaginary_frequencies(vibfreqs) -> Optional[int]:
    if vibfreqs is None:
        return None
    return sum(1 for f in vibfreqs if f < 0.0)


def infer_formula(atomnos: Optional[list[int]]) -> Optional[str]:
    """
    Very lightweight formula builder.
    You can replace this later with RDKit or ASE chemistry helpers if desired.
    """
    if not atomnos:
        return None

    periodic = {
        1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 15: "P", 16: "S", 17: "Cl",
        35: "Br", 53: "I"
    }

    counts = {}
    for z in atomnos:
        sym = periodic.get(z, f"Z{z}")
        counts[sym] = counts.get(sym, 0) + 1

    # crude Hill-style-ish ordering
    ordered = []
    if "C" in counts:
        ordered.append(("C", counts.pop("C")))
    if "H" in counts:
        ordered.append(("H", counts.pop("H")))
    ordered.extend(sorted(counts.items()))

    return "".join(sym if n == 1 else f"{sym}{n}" for sym, n in ordered)


def save_final_geometry(atomnos, atomcoords, xyz_path: str) -> Optional[str]:
    if atomnos is None or atomcoords is None:
        return None
    if len(atomcoords) == 0:
        return None

    atoms = Atoms(numbers=atomnos, positions=atomcoords[-1])
    ensure_dir(Path(xyz_path).parent)
    ase_write(xyz_path, atoms)
    return xyz_path


@flowcept_task
def parse_output_to_record(
    output_file: str,
    calc_id: str,
    molecule_id: str,
    role: str,
    fragment_id: Optional[str] = None,
    smiles: Optional[str] = None,
    parent_smiles: Optional[str] = None,
    broken_bond: Optional[str] = None,
    job_type: Optional[str] = "opt_freq",
    xyz_dir: str = "geometries",
    extra_metadata: Optional[Dict[str, Any]] = None,
) -> CalculationRecord:
    """
    Parse a computational chemistry output file into a standardized record.
    """

    output_path = Path(output_file)
    text = read_text(output_file)

    # --- cclib parse ---
    data = ccread(output_file)

    # --- fallback / program-specific parsing ---
    method, basis = extract_method_basis(text)
    charge_fb, multiplicity_fb = extract_charge_multiplicity(text)
    success = extract_success(text)
    job_meta = extract_job_metadata(text)

    atomnos = getattr(data, "atomnos", None)
    atomcoords = getattr(data, "atomcoords", None)
    vibfreqs = getattr(data, "vibfreqs", None)

    # cclib energies are often in eV for some attrs, but scfenergies are usually array-like
    # We'll store final SCF energy as-is here only if present and convert if needed later.
    scf_energy = extract_final_energy_hartree(text)
    #scf_energy = None
    #if hasattr(data, "scfenergies") and data.scfenergies is not None and len(data.scfenergies) > 0:
        # NOTE: cclib scfenergies are often in eV.
        # You may want to convert these if needed depending on parser behavior for NWChem.
    #    scf_energy = float(data.scfenergies[-1])

    charge = getattr(data, "charge", None)
    multiplicity = getattr(data, "mult", None)

    if charge is None:
        charge = charge_fb
    if multiplicity is None:
        multiplicity = multiplicity_fb

    xyz_path = str(Path(xyz_dir) / f"{calc_id}.xyz")
    xyz_file = save_final_geometry(atomnos, atomcoords, xyz_path)

    record = CalculationRecord(
        calc_id=calc_id,
        molecule_id=molecule_id,
        fragment_id=fragment_id,
        role=role,

        smiles=smiles,
        parent_smiles=parent_smiles,
        broken_bond=broken_bond,
        formula=infer_formula(atomnos),
        charge=charge,
        multiplicity=multiplicity,
        n_atoms=len(atomnos) if atomnos is not None else None,

        software="NWChem",
        job_type=job_type,
        method=method,
        basis=basis,
        input_file=None,
        output_file=str(output_path),
        xyz_file=xyz_file,
        job_id=job_meta.get("job_id"),
        hostname=job_meta.get("hostname"),
        timestamp=None,

        success=success,
        optimization_converged=None,  # can improve later
        imaginary_freq_count=count_imaginary_frequencies(vibfreqs),

        scf_energy_hartree=scf_energy,  # WARNING: verify unit for your NWChem/cclib combo
        zpve_hartree=None,
        enthalpy_corr_hartree=None,
        gibbs_corr_hartree=None,
        entropy_cal_molK=None,

        vibfreqs=vibfreqs.tolist() if hasattr(vibfreqs, "tolist") else vibfreqs,
        atomnos=atomnos.tolist() if hasattr(atomnos, "tolist") else atomnos,
        final_coords=atomcoords[-1].tolist() if atomcoords is not None else None,

        extra=extra_metadata or {},
    )

    return record
