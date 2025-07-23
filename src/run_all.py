import os
import subprocess
import time
from bde_workflow import run_bde
from flowcept import Flowcept, flowcept_task

@flowcept_task
def run_nwchem_job(nw_file, output_dir):
    log_file = os.path.splitext(nw_file)[0] + ".out"
    nw_path = os.path.join(output_dir, nw_file)
    log_path = os.path.join(output_dir, log_file)
    print(f"Running: {nw_file}")
    try:
        with open(log_path, "w") as log:
            subprocess.run(
                ["srun", "-n56", "/lustre/orion/proj-shared/stf053/ketan/isc/nwchem-install-3/bin/nwchem", nw_path],
                stdout=log, stderr=subprocess.STDOUT, check=True)
        status = "Completed"
    except subprocess.CalledProcessError:
        status = "failed"
    job_result = {
        "input_file": nw_file,
        "output_file": log_file,
        "status": status
    }
    return job_result


@flowcept_task
def run_nwchem_jobs(output_dir):
    job_results = []
    files = [f for f in os.listdir(output_dir) if f.endswith(".nw")]
    for nw_file in files:
        job_result = run_nwchem_job(nw_file, output_dir)
        job_results.append(job_result)
    return job_results

@flowcept_task
def wait_for_jobs(output_dir):
    print("\nWaiting for all jobs to finish...")
    while True:
        running = [f for f in os.listdir(output_dir) if f.endswith(".nw") and not os.path.exists(f"{output_dir}/{f[:-3]}.out")]
        if not running:
            print("All jobs finished.\n")
            return {"status": "complete", "unfinished": []}
        print(f"{len(running)} jobs still running...")
        time.sleep(10)

@flowcept_task
def main():
    with Flowcept():
        smiles = input("Enter SMILES string for parent molecule: ").strip()
        output_dir = "bde_calc"

        # Step 1: generate input files
        bde_input = run_bde(smiles=smiles, outdir=output_dir, generate_inputs=True, parse_outputs=False)
    
        # Step 2: run NWChem jobs
        job_results = run_nwchem_jobs(output_dir)

        # Step 3: wait until all output files exist
        wait_result = wait_for_jobs(output_dir=output_dir)

        # Step 4: reparse results
        print("Parsing outputs...")
        bde_results = run_bde(smiles=smiles, outdir=output_dir, generate_inputs=False, parse_outputs=True)
        print("Results saved to bde_results.csv")

        return {
            "input_species": bde_input["species"],
            "job_status": job_results,
            "wait_status": wait_result,
            "bde_results": bde_results
        }

if __name__ == "__main__":
    main()
