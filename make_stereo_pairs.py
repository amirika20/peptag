"""
Build a CSV of D-Phe / L-Phe stereo pairs from test.csv.

For every sequence containing D-Phe ('f'), the L-Phe counterpart is obtained
by replacing all 'f' with 'F'.  Only pairs where both versions exist in the
dataset are kept.

Output columns:
    Sequence_f  – peptide sequence with D-Phe (f)
    Sequence_F  – peptide sequence with L-Phe (F)
    B_f         – B percentage for the D-Phe version
    B_F         – B percentage for the L-Phe version
"""

import pandas as pd
from pathlib import Path


def make_stereo_pairs(input_path: str | Path, output_path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(input_path)

    b_map = df.set_index("Peptide")["B"].to_dict()

    # Start from the D-Phe rows and look up the L-Phe counterpart
    f_rows = df[df["Peptide"].str.contains("f", regex=False)].copy()
    f_rows["counterpart"] = f_rows["Peptide"].str.replace("f", "F", regex=False)

    # Keep only pairs where the L-Phe counterpart exists
    mask = f_rows["counterpart"].isin(b_map)
    pairs = f_rows[mask][["Peptide", "B", "counterpart"]].copy()
    pairs["B_F"] = pairs["counterpart"].map(b_map)

    pairs = pairs.rename(columns={
        "Peptide": "Sequence_f",
        "B":       "B_f",
        "counterpart": "Sequence_F",
    })

    pairs = pairs[["Sequence_f", "Sequence_F", "B_f", "B_F"]]
    pairs = pairs.sort_values("Sequence_f").reset_index(drop=True)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(output_path, index=False)
    print(f"Saved {len(pairs)} stereo pairs → {output_path}")
    return pairs


if __name__ == "__main__":
    make_stereo_pairs(
        input_path="curated_data/test.csv",
        output_path="curated_data/stereo_pairs.csv",
    )
