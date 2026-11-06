"""
make_natural_benchmark.py
==========================
Build a "natural-only" benchmark by removing every D-Phenylalanine ('f')
tagged peptide from train.csv, val.csv, and test.csv.

D-Phe is the only non-proteinogenic residue in this dataset (encoded as
lowercase 'f'); L-Phe-tagged and label-free peptides are built entirely from
natural L-amino acids. This benchmark lets models be evaluated on natural
sequences alone, without any stereochemical outlier.

Run after split_data.py has produced train.csv / val.csv / test.csv:
    python make_natural_benchmark.py

Three output files are produced:
  curated_data/natural_train.csv
  curated_data/natural_val.csv
  curated_data/natural_test.csv
"""

import pandas as pd
from pathlib import Path


def make_natural_split(input_path: str | Path, output_path: str | Path) -> pd.DataFrame:
    """Drop D-Phe-tagged peptides (sequences containing lowercase 'f')."""
    input_path = Path(input_path)
    output_path = Path(output_path)

    df = pd.read_csv(input_path)
    natural = df[~df["Peptide"].str.contains("f")].reset_index(drop=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    natural.to_csv(output_path, index=False)
    print(f"{input_path.name}: {len(df)} -> {len(natural)} rows (dropped {len(df) - len(natural)} D-Phe peptides) -> {output_path}")
    return natural


if __name__ == "__main__":
    for split in ["train", "val", "test"]:
        make_natural_split(
            input_path=f"curated_data/{split}.csv",
            output_path=f"curated_data/natural_{split}.csv",
        )
