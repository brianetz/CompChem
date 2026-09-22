import re
from pathlib import Path
from typing import Optional, Dict, Any


def read_text(outfile: str) -> str:
    return Path(outfile).read_text(errors="ignore")


def extract_method_basis(text: str) -> tuple[Optional[str], Optional[str]]:
    method = None
    basis = None

    # Example patterns - adjust based on your NWChem inputs
    m = re.search(r"xc\s+([A-Za-z0-9\-\_]+)", text, re.IGNORECASE)
    if m:
        method = m.group(1)

    b = re.search(r"basis\s+.*?\n(.*?)\nend", text, re.IGNORECASE | re.DOTALL)
    if b:
        basis_block = b.group(1).strip().splitlines()
        if basis_block:
            basis = basis_block[0].strip()

    return method, basis


def extract_charge_multiplicity(text: str) -> tuple[Optional[int], Optional[int]]:
    charge = None
    multiplicity = None

    c = re.search(r"charge\s+(-?\d+)", text, re.IGNORECASE)
    if c:
        charge = int(c.group(1))

    m = re.search(r"(?:mult|multiplicity)\s+(\d+)", text, re.IGNORECASE)
    if m:
        multiplicity = int(m.group(1))

    return charge, multiplicity


def extract_success(text: str) -> bool:
    success_markers = [
        "Total times  cpu:",
        "CITATION",
        "Task  times  cpu:",
    ]
    return any(marker in text for marker in success_markers)


def extract_job_metadata(text: str) -> Dict[str, Any]:
    metadata = {}

    host = re.search(r"hostname\s*=\s*(\S+)", text, re.IGNORECASE)
    if host:
        metadata["hostname"] = host.group(1)

    job = re.search(r"(?:jobid|job id)\s*[:=]\s*(\S+)", text, re.IGNORECASE)
    if job:
        metadata["job_id"] = job.group(1)

    return metadata

def extract_final_energy_hartree(text: str) -> Optional[float]:
    """
    Extract final DFT/SCF total energy from NWChem output in Hartree.
    Adjust regex if your output formatting differs.
    """
    patterns = [
        r"Total DFT energy\s*=\s*(-?\d+\.\d+)",
        r"Total SCF energy\s*=\s*(-?\d+\.\d+)",
    ]

    for pattern in patterns:
        matches = re.findall(pattern, text, re.IGNORECASE)
        if matches:
            return float(matches[-1])

    return None
