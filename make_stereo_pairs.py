"""
Build a CSV of D-Phe / L-Phe stereo pairs from test.csv.

For every sequence containing D-Phe ('f'), the L-Phe counterpart is obtained
by replacing all 'f' with 'F'.  Only pairs where both versions exist in the
dataset AND whose %B difference meets the minimum threshold are kept.

A pair is considered meaningful only when |B_f - B_F| >= min_delta_b.
Pairs below this threshold are discarded because the retention shift is
within experimental noise and does not represent a reliable stereochemical
effect on chromatographic behaviour.

Output columns:
    Sequence_f  – peptide sequence with D-Phe (f)
    Sequence_F  – peptide sequence with L-Phe (F)
    SMILES_f    – SMILES string for the D-Phe version
    SMILES_F    – SMILES string for the L-Phe version
    B_f         – %B at elution for the D-Phe version
    B_F         – %B at elution for the L-Phe version
    delta_B     – signed difference B_f − B_F (positive = D-Phe elutes later)
"""

import pandas as pd
from pathlib import Path


def make_stereo_pairs(
    input_path: str | Path,
    output_path: str | Path,
    min_delta_b: float = 1.0,
) -> pd.DataFrame:
    """Extract D-Phe / L-Phe stereo pairs with a minimum %B difference filter.

    Parameters
    ----------
    input_path  : CSV file to read (must have 'Peptide' and 'B' columns).
    output_path : Destination CSV file.
    min_delta_b : Minimum absolute difference in %B required to keep a pair.
                  Pairs with |B_f - B_F| < min_delta_b are discarded as
                  not meaningfully separated.  Default: 1.0 %B.

    Returns
    -------
    pd.DataFrame with columns: Sequence_f, Sequence_F, SMILES_f, SMILES_F, B_f, B_F, delta_B.
    """
    df = pd.read_csv(input_path)

    b_map = df.set_index("Peptide")["B"].to_dict()
    smiles_map = df.set_index("Peptide")["SMILES"].to_dict()

    # Start from the D-Phe rows and look up the L-Phe counterpart
    f_rows = df[df["Peptide"].str.contains("f", regex=False)].copy()
    f_rows["counterpart"] = f_rows["Peptide"].str.replace("f", "F", regex=False)

    # Keep only pairs where the L-Phe counterpart exists in the dataset
    mask = f_rows["counterpart"].isin(b_map)
    pairs = f_rows[mask][["Peptide", "SMILES", "B", "counterpart"]].copy()
    pairs["B_F"] = pairs["counterpart"].map(b_map)
    pairs["SMILES_F"] = pairs["counterpart"].map(smiles_map)

    pairs = pairs.rename(columns={
        "Peptide": "Sequence_f",
        "SMILES":  "SMILES_f",
        "B":       "B_f",
        "counterpart": "Sequence_F",
    })

    pairs = pairs[["Sequence_f", "Sequence_F", "SMILES_f", "SMILES_F", "B_f", "B_F"]]

    # Compute signed retention shift: positive means D-Phe elutes later than L-Phe
    pairs["delta_B"] = pairs["B_f"] - pairs["B_F"]

    # Apply minimum meaningful threshold: discard pairs where the retention
    # difference is too small to be distinguished from experimental noise
    n_before = len(pairs)
    pairs = pairs[pairs["delta_B"].abs() >= min_delta_b]
    n_dropped = n_before - len(pairs)

    pairs = pairs.sort_values("Sequence_f").reset_index(drop=True)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(output_path, index=False)

    print(f"Total pairs found       : {n_before}")
    print(f"Dropped (|ΔB| < {min_delta_b:.1f} %B): {n_dropped}")
    print(f"Kept (|ΔB| ≥ {min_delta_b:.1f} %B)   : {len(pairs)}")
    print(f"Saved → {output_path}")
    return pairs


if __name__ == "__main__":
    make_stereo_pairs(
        input_path="curated_data/test.csv",
        output_path="curated_data/stereo_pairs.csv",
        min_delta_b=1.0,
    )
