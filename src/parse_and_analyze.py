from pathlib import Path

from parsers.cclib_parser import parse_output_to_record
from analysis.qc_checks import apply_qc_checks
from storage.csv_store import save_records_csv, save_records_parquet
from analysis.bde import save_bde_table
from utils.paths import ensure_dir
from storage.flowcept_store import provenance_session


def infer_role_from_filename(path: Path) -> str:
    name = path.stem.lower()
    if "frag" in name:
        return "fragment"
    return "parent"


def infer_fragment_id(path: Path) -> str | None:
    name = path.stem.lower()
    if "frag" in name:
        return path.stem
    return None


def infer_molecule_id(path: Path) -> str:
    # Replace this with your real naming convention
    return path.parent.name


def main():
    output_dir = Path("outputs")
    results_dir = ensure_dir("results")
    geometry_dir = ensure_dir("geometries")

    records = []

    with provenance_session():
        for outfile in output_dir.rglob("*.out"):
            role = infer_role_from_filename(outfile)
            fragment_id = infer_fragment_id(outfile)
            molecule_id = infer_molecule_id(outfile)
            calc_id = outfile.stem

            record = parse_output_to_record(
                output_file=str(outfile),
                calc_id=calc_id,
                molecule_id=molecule_id,
                role=role,
                fragment_id=fragment_id,
                xyz_dir=str(geometry_dir),
            )

            record = apply_qc_checks(record)
            records.append(record)

        save_records_csv(records, str(results_dir / "calculations.csv"))
        #save_records_parquet(records, str(results_dir / "calculations.parquet"))
        save_bde_table(records, str(results_dir / "bde_results.csv"))

    print(f"Parsed {len(records)} calculations.")
    print(f"Saved results to {results_dir}")


if __name__ == "__main__":
    main()
