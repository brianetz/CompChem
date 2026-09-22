from dataclasses import dataclass, asdict, field
from typing import Optional, List, Dict, Any


@dataclass
class CalculationRecord:
    # --- identity ---
    calc_id: str
    molecule_id: str
    fragment_id: Optional[str]
    role: str  # "parent", "fragment", "conformer", etc.

    # --- chemistry ---
    smiles: Optional[str]
    parent_smiles: Optional[str]
    broken_bond: Optional[str]
    formula: Optional[str]
    charge: Optional[int]
    multiplicity: Optional[int]
    n_atoms: Optional[int]

    # --- provenance ---
    software: str
    job_type: Optional[str]
    method: Optional[str]
    basis: Optional[str]
    input_file: Optional[str]
    output_file: str
    xyz_file: Optional[str]
    job_id: Optional[str]
    hostname: Optional[str]
    timestamp: Optional[str]

    # --- status / QC ---
    success: bool
    optimization_converged: Optional[bool]
    imaginary_freq_count: Optional[int]
    qc_flag: Optional[str] = None
    qc_notes: Optional[str] = None

    # --- energies ---
    scf_energy_hartree: Optional[float] = None
    zpve_hartree: Optional[float] = None
    enthalpy_corr_hartree: Optional[float] = None
    gibbs_corr_hartree: Optional[float] = None
    entropy_cal_molK: Optional[float] = None

    # --- optional parsed arrays / metadata ---
    vibfreqs: Optional[List[float]] = None
    atomnos: Optional[List[int]] = None
    final_coords: Optional[List[List[float]]] = None

    # --- extra metadata bucket ---
    extra: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)