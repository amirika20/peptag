"""
split_data.py
=============
Split trainval.csv into stratified train and validation sets.

Stratification key: (terminal_type, phe_type, length)
  - terminal_type : 'K' if the peptide ends in K / Kf / KF (Lys C-terminus)
                    'R' if the peptide ends in R (Arg C-terminus)
  - phe_type      : 'D'  – D-Phe tag (sequence contains lowercase 'f')
                    'L'  – L-Phe tag (sequence starts with uppercase 'F')
                    'LF' – label-free (no phenylalanine tag)
  - length        : total sequence length

Stratifying on all three factors guarantees that every combination of
terminal type, tagging variant, and length is proportionally represented
in both train and val.  This is important because each stratum has a
systematically different RT distribution.

For the single stratum with only 3 sequences (R_LF_L7), all 3 are kept;
2 go to train and 1 to val — an acceptable imbalance given the tiny size.

Run after merge_data.py:
    python split_data.py
"""

import pandas as pd
from pathlib import Path
from sklearn.model_selection import train_test_split


def _terminal(seq: str) -> str:
    """Return 'R' for Arg-terminated peptides, 'K' for Lys-terminated."""
    return "R" if seq[-1] == "R" else "K"


def _phe_type(seq: str) -> str:
    """Return the phenylalanine tagging variant of a sequence.

    Returns
    -------
    'D'  : sequence contains lowercase 'f' (D-Phenylalanine tag)
    'L'  : sequence starts with 'F' (L-Phenylalanine tag)
    'LF' : no phenylalanine tag (label-free)
    """
    if "f" in seq:
        return "D"
    if seq[0] == "F":
        return "L"
    return "LF"


def split_data(
    input_path: str | Path,
    train_path: str | Path,
    val_path: str | Path,
    val_fraction: float = 0.1,
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Stratified train/val split of a merged peptide CSV.

    Parameters
    ----------
    input_path   : path to trainval.csv (or any file with a 'Peptide' column)
    train_path   : output path for the training set
    val_path     : output path for the validation set
    val_fraction : fraction of data to use for validation (default 0.1 = 10 %)
    seed         : random seed for reproducibility
    """
    df = pd.read_csv(input_path)

    # Build a compound stratum label, e.g. "K_D_L8" for an 8-residue
    # D-Phe K-term peptide.  Used only for stratification; dropped afterwards.
    df["_term"]    = df["Peptide"].apply(_terminal)
    df["_phe"]     = df["Peptide"].apply(_phe_type)
    df["_stratum"] = df["_term"] + "_" + df["_phe"] + "_L" + df["Length"].astype(str)

    train_idx, val_idx = train_test_split(
        df.index,
        test_size=val_fraction,
        random_state=seed,
        stratify=df["_stratum"],   # ensures proportional representation in each stratum
    )

    train_df = df.loc[train_idx].drop(columns=["_term", "_phe", "_stratum"])
    val_df   = df.loc[val_idx].drop(columns=["_term", "_phe", "_stratum"])

    train_df = train_df.sort_values("Peptide").reset_index(drop=True)
    val_df   = val_df.sort_values("Peptide").reset_index(drop=True)

    Path(train_path).parent.mkdir(parents=True, exist_ok=True)
    Path(val_path).parent.mkdir(parents=True, exist_ok=True)
    train_df.to_csv(train_path, index=False)
    val_df.to_csv(val_path, index=False)

    # ---- summary -----------------------------------------------------------
    print(f"Train : {len(train_df):>6}  ({100*(1-val_fraction):.0f}%)")
    print(f"Val   : {len(val_df):>6}  ({100*val_fraction:.0f}%)")
    print()

    # Per-stratum check: term × phe × length all present in both splits
    for split_name, split_df in [("Train", train_df), ("Val", val_df)]:
        terms  = set(split_df["Peptide"].apply(_terminal))
        phes   = set(split_df["Peptide"].apply(_phe_type))
        lengths = set(split_df["Length"])
        print(f"{split_name}:")
        print(f"  terminal types : {sorted(terms)}")
        print(f"  phe types      : {sorted(phes)}")
        print(f"  lengths        : {sorted(lengths)}")
        print()

    return train_df, val_df


if __name__ == "__main__":
    split_data(
        input_path="curated_data/trainval.csv",
        train_path="curated_data/train.csv",
        val_path="curated_data/val.csv",
        val_fraction=0.1,
        seed=42,
    )
