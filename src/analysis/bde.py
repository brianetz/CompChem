from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple
import pandas as pd

from storage.records import CalculationRecord
from utils.units import delta_hartree_to_kcal
from storage.flowcept_store import flowcept_task


def records_to_df(records: List[CalculationRecord]) -> pd.DataFrame:
    return pd.DataFrame([r.to_dict() for r in records])


def compute_bde_kcal(parent_e: float, frag_a_e: float, frag_b_e: float) -> float:
    delta = frag_a_e + frag_b_e - parent_e
    return delta_hartree_to_kcal(delta)


def build_bde_table(records: List[CalculationRecord]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Assumes fragments belonging to a bond break share the same broken_bond label
    and molecule_id.

    Returns (bde_df, excluded_df): only records whose qc_flag is "PASS" (set
    by analysis.qc_checks.apply_qc_checks) are used to compute BDEs - this is
    the single source of truth for what's trustworthy, rather than this
    function re-deriving its own success/energy check that could drift out
    of sync with the critic. Everything else (WARN or FAIL - e.g. imaginary
    frequencies, missing energy, a failed job) is returned separately in
    excluded_df with its qc_notes, so nothing is silently dropped.
    """

    df = records_to_df(records)

    if "qc_flag" not in df.columns:
        raise ValueError(
            "records have no qc_flag - run analysis.qc_checks.apply_qc_checks "
            "on every record before calling build_bde_table."
        )

    excluded = df[df["qc_flag"] != "PASS"].copy()
    df = df[df["qc_flag"] == "PASS"].copy()

    parents = df[df["role"] == "parent"].copy()
    frags = df[df["role"] == "fragment"].copy()

    results = []

    for _, parent in parents.iterrows():
        mol_id = parent["molecule_id"]
        broken_bond = parent.get("broken_bond", None)

        # In many workflows, parent may not have broken_bond.
        # So instead, find all fragment pairs for that molecule.
        frag_subset = frags[frags["molecule_id"] == mol_id]

        if "broken_bond" in frag_subset.columns and frag_subset["broken_bond"].notna().any():
            grouped = frag_subset.groupby("broken_bond")
        else:
            grouped = [("UNKNOWN", frag_subset)]

        for bond_label, group in grouped:
            if len(group) != 2:
                continue

            frag_rows = group.sort_values("fragment_id")
            frag_a = frag_rows.iloc[0]
            frag_b = frag_rows.iloc[1]

            try:
                bde = compute_bde_kcal(
                    parent["scf_energy_hartree"],
                    frag_a["scf_energy_hartree"],
                    frag_b["scf_energy_hartree"],
                )
            except Exception:
                continue

            results.append({
                "molecule_id": mol_id,
                "bond_label": bond_label,
                "parent_calc_id": parent["calc_id"],
                "frag_a_calc_id": frag_a["calc_id"],
                "frag_b_calc_id": frag_b["calc_id"],
                "parent_energy_hartree": parent["scf_energy_hartree"],
                "frag_a_energy_hartree": frag_a["scf_energy_hartree"],
                "frag_b_energy_hartree": frag_b["scf_energy_hartree"],
                "bde_electronic_kcal_mol": bde,
            })

    return pd.DataFrame(results), excluded


@flowcept_task
def save_bde_table(records: List[CalculationRecord], outfile: str) -> None:
    df, excluded = build_bde_table(records)
    df.to_csv(outfile, index=False)

    if len(excluded) > 0:
        excluded_path = str(Path(outfile).with_name(Path(outfile).stem + "_excluded.csv"))
        cols = [c for c in ("calc_id", "molecule_id", "role", "broken_bond", "qc_flag", "qc_notes") if c in excluded.columns]
        excluded[cols].to_csv(excluded_path, index=False)
        print(f"WARNING: {len(excluded)} record(s) excluded from BDE table by QC flag - see {excluded_path}")
