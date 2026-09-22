import re
from flowcept import Flowcept, flowcept_task

@flowcept_task
def parse_nwchem_output(filename):
    energy = zpe = enthalpy = entropy = None

    with open(filename) as f:
        for line in f:
            if "Total DFT energy" in line:
                energy = float(line.split('=')[1].strip())
            elif "Zero-Point correction to Energy" in line:
                # Extract only the kcal/mol value and convert to Hartree
                zpe_kcal = float(line.split('=')[1].split("kcal")[0].strip())
                zpe = zpe_kcal / 627.509
            elif "Thermal correction to Enthalpy" in line:
                enthalpy_kcal = float(line.split('=')[1].split("kcal")[0].strip())
                enthalpy = enthalpy_kcal / 627.509
            elif "Total Entropy" in line:
                entropy_cal = float(line.split('=')[1].split("cal")[0].strip())  # in cal/mol·K
                entropy = entropy_cal / 1000

    if energy is None:
        raise ValueError(f"Total energy not found in {filename}")
    if zpe is None:
        print(f"Warning: ZPE not found in {filename}")
        #zpe = 0.0
    if enthalpy is None:
        print(f"Warning: Enthalpy not found in {filename}")
        #enthalpy = 0.0
    if entropy is None:
        print(f"Warning: Entropy not found in {filename}")
        #entropy = 0.0

    return {
        "energy": energy,
        "zpe": zpe,
        "enthalpy": enthalpy,
        "entropy": entropy
    }
