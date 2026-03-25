"""
Build a CSV of single-residue substitution pairs from test.csv.

Two peptides form a substitution pair when they have the same length and differ
at exactly one sequence position.  All standard (non-tag) residues as well as
tagged (f/F-prefixed) sequences are included, so e.g. fAPGIHR vs fAPTIHR is a
valid pair.

Only pairs where |B_1 - B_2| >= min_delta_b are kept.

Output columns:
    Sequence_1  – first peptide sequence (alphabetically earlier)
    Sequence_2  – second peptide sequence
    Position    – 1-indexed position where the two sequences differ
    Residue_1   – residue in Sequence_1 at that position
    Residue_2   – residue in Sequence_2 at that position
    SMILES_1    – SMILES string for Sequence_1
    SMILES_2    – SMILES string for Sequence_2
    B_1         – %B at elution for Sequence_1
    B_2         – %B at elution for Sequence_2
    delta_B     – signed difference B_1 − B_2 (positive = Sequence_1 elutes later)
"""

import pandas as pd
from itertools import combinations
from pathlib import Path


def _count_differences(s1: str, s2: str) -> tuple[int, int]:
    """Return (number of differing positions, index of the differing position).

    Returns (-1, -1) if the sequences differ at more than one position.
    The index is 0-based.
    """
    diff_count = 0
    diff_pos = -1
    for i, (a, b) in enumerate(zip(s1, s2)):
        if a != b:
            diff_count += 1
            if diff_count == 1:
                diff_pos = i
            else:
                return -1, -1
    return diff_count, diff_pos


def make_substitution_pairs(
    input_path: str | Path,
    output_path: str | Path,
    min_delta_b: float = 1.0,
) -> pd.DataFrame:
    """Extract single-residue substitution pairs with a minimum %B difference filter.

    Parameters
    ----------
    input_path  : CSV file to read (must have 'Peptide', 'B', and 'SMILES' columns).
    output_path : Destination CSV file.
    min_delta_b : Minimum absolute difference in %B required to keep a pair.
                  Default: 1.0 %B.

    Returns
    -------
    pd.DataFrame with columns: Sequence_1, Sequence_2, Position, Residue_1,
                                Residue_2, SMILES_1, SMILES_2, B_1, B_2, delta_B.
    """
    df = pd.read_csv(input_path)

    b_map = df.set_index("Peptide")["B"].to_dict()
    smiles_map = df.set_index("Peptide")["SMILES"].to_dict()

    # Group sequences by length for efficient comparison
    from collections import defaultdict
    by_length: dict[int, list[str]] = defaultdict(list)
    for seq in df["Peptide"]:
        by_length[len(seq)].append(seq)

    records = []
    for length, seqs in by_length.items():
        for s1, s2 in combinations(seqs, 2):
            n_diff, diff_pos = _count_differences(s1, s2)
            if n_diff != 1:
                continue
            # Ensure s1 is alphabetically first for consistency
            if s1 > s2:
                s1, s2 = s2, s1
                diff_pos = diff_pos  # position unchanged
            records.append({
                "Sequence_1": s1,
                "Sequence_2": s2,
                "Position":   diff_pos + 1,  # 1-indexed
                "Residue_1":  s1[diff_pos],
                "Residue_2":  s2[diff_pos],
                "SMILES_1":   smiles_map[s1],
                "SMILES_2":   smiles_map[s2],
                "B_1":        b_map[s1],
                "B_2":        b_map[s2],
            })

    pairs = pd.DataFrame(records)

    if pairs.empty:
        print("No pairs found.")
        return pairs

    pairs["delta_B"] = pairs["B_1"] - pairs["B_2"]

    n_before = len(pairs)
    pairs = pairs[pairs["delta_B"].abs() >= min_delta_b]
    n_dropped = n_before - len(pairs)

    pairs = pairs.sort_values(["Sequence_1", "Sequence_2"]).reset_index(drop=True)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(output_path, index=False)

    print(f"Total pairs found       : {n_before}")
    print(f"Dropped (|ΔB| < {min_delta_b:.1f} %B): {n_dropped}")
    print(f"Kept (|ΔB| ≥ {min_delta_b:.1f} %B)   : {len(pairs)}")
    print(f"Saved → {output_path}")
    return pairs


if __name__ == "__main__":
    make_substitution_pairs(
        input_path="curated_data/test.csv",
        output_path="curated_data/substitution_pairs.csv",
        min_delta_b=1.0,
    )
