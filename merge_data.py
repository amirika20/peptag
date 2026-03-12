"""
merge_data.py
=============
Merge all per-library curated CSVs into a single dataset file, enriched with:
  - Sequence length (number of residues).
  - Isomeric SMILES with explicit alpha-carbon stereochemistry.

Cross-library deduplication: if the same peptide sequence appears in more
than one library file, all rows are collapsed into a single entry by:
  - averaging B (%B at elution) and M (observed mass),
  - taking the mode of Z (precursor charge state).

Run this script after data_processing.py has populated curated_data/:
    python merge_data.py

Two output files are produced:
  curated_data/merged.csv  – train + val pool (all libraries except 'test')
  curated_data/test.csv    – held-out benchmark (only the 'test' library)
"""

import pandas as pd
from pathlib import Path

from data_processing import peptide_to_smiles


def merge_libraries(
    curated_path: str | Path,
    output_path: str | Path,
    folders: list[str] | None = None,
) -> pd.DataFrame:
    """Merge curated CSV files into one output file.

    Parameters
    ----------
    curated_path : path to the curated_data directory
    output_path  : destination CSV file
    folders      : if given, only these subdirectory names are included;
                   if None, all subdirectories except 'test' are included
    """
    curated_path = Path(curated_path)
    output_path = Path(output_path)

    frames = []
    for lib_dir in sorted(curated_path.iterdir()):
        if not lib_dir.is_dir():
            continue
        # Respect the folder filter: skip 'test' when building the train+val pool,
        # or skip any directory not in the explicit allow-list when one is given.
        if folders is None:
            if lib_dir.name == "test":
                continue
        else:
            if lib_dir.name not in folders:
                continue

        for csv_file in sorted(lib_dir.glob("*.csv")):
            df = pd.read_csv(csv_file)
            # Drop raw RT (seconds) — only the calibrated %B column is kept
            df = df.drop(columns=["RT"])
            frames.append(df)
            print(f"  loaded {lib_dir.name}/{csv_file.name}  ({len(df)} rows)")

    if not frames:
        raise RuntimeError(f"No curated CSV files found under {curated_path}")

    merged = pd.concat(frames, ignore_index=True)
    print(f"\nTotal rows before deduplication: {len(merged)}")

    # ---- add sequence length ----------------------------------------------
    merged["Length"] = merged["Peptide"].str.len()

    # ---- generate isomeric SMILES for each unique sequence ----------------
    # SMILES are generated once per unique sequence to avoid redundant work.
    # Failures (e.g. unknown residue codes) are stored as None and reported.
    unique_seqs = merged["Peptide"].unique()
    print(f"Generating SMILES for {len(unique_seqs)} unique sequences...")

    smiles_map = {}
    errors = []
    for seq in unique_seqs:
        try:
            smiles_map[seq] = peptide_to_smiles(seq)
        except Exception as e:
            smiles_map[seq] = None
            errors.append((seq, str(e)))

    if errors:
        print(f"  WARNING: SMILES generation failed for {len(errors)} sequence(s):")
        for seq, err in errors:
            print(f"    {seq}: {err}")

    merged["SMILES"] = merged["Peptide"].map(smiles_map)

    # ---- cross-library deduplication -------------------------------------
    # If a peptide appears in multiple library files (e.g. the same sequence
    # synthesised in two different batches), collapse to one row by:
    #   B  → mean   (average %B across observations)
    #   M  → mean   (average observed mass)
    #   Z  → mode   (most common charge state; ties broken by first occurrence)
    def mode_first(s):
        m = s.mode()
        return m.iloc[0] if not m.empty else s.iloc[0]

    merged = (
        merged.groupby("Peptide", sort=False)
        .agg(
            B=("B", "mean"),
            M=("M", "mean"),
            Z=("Z", mode_first),
            Length=("Length", "first"),
            SMILES=("SMILES", "first"),
        )
        .reset_index()
    )

    merged = merged[["Peptide", "B", "M", "Z", "Length", "SMILES"]]
    merged = merged.sort_values("Peptide").reset_index(drop=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_path, index=False)
    print(f"\nSaved {len(merged)} unique peptides → {output_path}")
    return merged


if __name__ == "__main__":
    merge_libraries(
        curated_path="curated_data",
        output_path="curated_data/merged.csv",
    )
    merge_libraries(
        curated_path="curated_data",
        output_path="curated_data/test.csv",
        folders=["test"],
    )
