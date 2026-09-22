from pathlib import Path
from typing import List

import pandas as pd

from storage.records import CalculationRecord
from utils.paths import ensure_dir
from storage.flowcept_store import flowcept_task


def records_to_dataframe(records: List[CalculationRecord]) -> pd.DataFrame:
    return pd.DataFrame([r.to_dict() for r in records])


@flowcept_task
def save_records_csv(records: List[CalculationRecord], outfile: str) -> None:
    ensure_dir(Path(outfile).parent)
    df = records_to_dataframe(records)
    df.to_csv(outfile, index=False)


def save_records_parquet(records: List[CalculationRecord], outfile: str) -> None:
    ensure_dir(Path(outfile).parent)
    df = records_to_dataframe(records)
    df.to_parquet(outfile, index=False)