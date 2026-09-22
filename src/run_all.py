import os

import parsl
from parsl import bash_app

from bde_workflow import run_bde
from parsl_config import get_local_config
from storage.flowcept_store import flowcept_task, provenance_session

# Real Frontier invocation by default; override locally for testing, e.g.:
#   export NWCHEM_CMD="python $(pwd)/mock_nwchem.py"
NWCHEM_CMD = os.environ.get(
    "NWCHEM_CMD",
    "srun -n56 /lustre/orion/proj-shared/stf053/ketan/isc/nwchem-install-3/bin/nwchem",
)


@bash_app
def run_nwchem_job(nw_path: str, nwchem_cmd: str, stdout=None, stderr=None):
    # Parsl bash_apps run in a separate worker process, so they can't see
    # module-level globals from the process that defined them - nwchem_cmd
    # has to come in as an explicit argument (not read from NWCHEM_CMD at
    # call time), or every job fails with "name 'NWCHEM_CMD' is not defined".
    return f"{nwchem_cmd} {nw_path}"


@flowcept_task
def run_nwchem_jobs(output_dir):
    """Submit every .nw file in output_dir to Parsl concurrently, then wait
    for all of them. This replaces the old subprocess-per-file loop (which
    was effectively serial - each subprocess.run() call blocked before the
    next job was even submitted) and the separate time.sleep(10) polling
    loop that used to wait for .out files to appear: Parsl's futures give us
    both real concurrency and a real "is it done" signal for free.
    """
    files = [f for f in os.listdir(output_dir) if f.endswith(".nw")]

    futures = {}
    for nw_file in files:
        base = os.path.splitext(nw_file)[0]
        nw_path = os.path.join(output_dir, nw_file)
        out_path = os.path.join(output_dir, base + ".out")
        err_path = os.path.join(output_dir, base + ".err")
        # Submitting here (not awaiting) is what makes every job in this
        # loop run concurrently rather than one at a time.
        futures[nw_file] = run_nwchem_job(nw_path, NWCHEM_CMD, stdout=out_path, stderr=err_path)

    job_results = []
    for nw_file, future in futures.items():
        base = os.path.splitext(nw_file)[0]
        try:
            future.result()  # blocks until THIS job is done; raises on failure
            status = "Completed"
        except Exception as exc:
            status = f"failed: {exc}"
        job_results.append({
            "input_file": nw_file,
            "output_file": base + ".out",
            "status": status,
        })
    return job_results


@flowcept_task
def main():
    with provenance_session():
        smiles = input("Enter SMILES string for parent molecule: ").strip()
        output_dir = "bde_calc"

        # Step 1: generate input files
        bde_input = run_bde(smiles=smiles, outdir=output_dir, generate_inputs=True, parse_outputs=False)

        # Step 2: run NWChem jobs concurrently via Parsl (no manual polling)
        job_results = run_nwchem_jobs(output_dir)

        # Step 3: parse outputs + compute BDEs
        print("Parsing outputs...")
        bde_results = run_bde(smiles=smiles, outdir=output_dir, generate_inputs=False, parse_outputs=True)
        print("Results saved to bde_results.csv")

        return {
            "input_species": bde_input["species"],
            "job_status": job_results,
            "bde_results": bde_results,
        }


if __name__ == "__main__":
    parsl.load(get_local_config())
    try:
        result = main()
        print(result)
    finally:
        parsl.dfk().cleanup()
