"""
Build a CSV of tagged / untagged pairs from test.csv.

"Tagged" peptides have f (D-Phe) or F (L-Phe) prepended as a terminal tag.
If the core sequence ends with K, the same tag is also appended at the end.

    untagged core      : APGIHR          → tagged : FAPGIHR
    untagged core (K)  : CADIEHWNK       → tagged : FCADIEHWNKF  (or f...f)

The untagged core is the sequence with the terminal f/F residues stripped.
Only pairs where both the core and its tagged counterpart exist in the
dataset AND whose %B difference meets the minimum threshold are kept.

A pair is considered meaningful only when |B_untagged - B_tagged| >= min_delta_b.

Output columns:
    Sequence_untagged – peptide sequence without terminal f/F tags (the core)
    Sequence_tagged   – peptide sequence with terminal f/F tag(s) added
    Tag               – tag type used: 'f' (D-Phe) or 'F' (L-Phe)
    SMILES_untagged   – SMILES string for the untagged version
    SMILES_tagged     – SMILES string for the tagged version
    B_untagged        – %B at elution for the untagged version
    B_tagged          – %B at elution for the tagged version
    delta_B           – signed difference B_untagged − B_tagged
                        (positive = untagged elutes later)
"""

import pandas as pd
from pathlib import Path


def _get_core(seq: str) -> str:
    """Strip leading (and matching trailing) f/F tag to recover the core sequence."""
    if not seq or seq[0] not in "fF":
        return seq
    tag = seq[0]
    core = seq[1:]
    if core and core[-1] == tag:
        core = core[:-1]
    return core


def make_tag_pairs(
    input_path: str | Path,
    output_path: str | Path,
    min_delta_b: float = 1.0,
) -> pd.DataFrame:
    """Extract tagged / untagged pairs with a minimum %B difference filter.

    Parameters
    ----------
    input_path  : CSV file to read (must have 'Peptide', 'B', and 'SMILES' columns).
    output_path : Destination CSV file.
    min_delta_b : Minimum absolute difference in %B required to keep a pair.
                  Pairs with |B_untagged - B_tagged| < min_delta_b are discarded.
                  Default: 1.0 %B.

    Returns
    -------
    pd.DataFrame with columns: Sequence_untagged, Sequence_tagged, Tag,
                                SMILES_untagged, SMILES_tagged,
                                B_untagged, B_tagged, delta_B.
    """
    df = pd.read_csv(input_path)

    b_map = df.set_index("Peptide")["B"].to_dict()
    smiles_map = df.set_index("Peptide")["SMILES"].to_dict()

    # Tagged rows: sequence starts with f or F (terminal Phe tag)
    tagged_rows = df[df["Peptide"].str[0].isin(["f", "F"])].copy()

    records = []
    for _, row in tagged_rows.iterrows():
        seq_tagged = row["Peptide"]
        tag = seq_tagged[0]
        core = _get_core(seq_tagged)

        if core in b_map:
            records.append({
                "Sequence_untagged": core,
                "Sequence_tagged":   seq_tagged,
                "Tag":               tag,
                "SMILES_untagged":   smiles_map[core],
                "SMILES_tagged":     row["SMILES"],
                "B_untagged":        b_map[core],
                "B_tagged":          row["B"],
            })

    pairs = pd.DataFrame(records)

    if pairs.empty:
        print("No pairs found.")
        return pairs

    pairs["delta_B"] = pairs["B_untagged"] - pairs["B_tagged"]

    n_before = len(pairs)
    pairs = pairs[pairs["delta_B"].abs() >= min_delta_b]
    n_dropped = n_before - len(pairs)

    pairs = pairs.sort_values(["Sequence_untagged", "Tag"]).reset_index(drop=True)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(output_path, index=False)

    print(f"Total pairs found       : {n_before}")
    print(f"Dropped (|ΔB| < {min_delta_b:.1f} %B): {n_dropped}")
    print(f"Kept (|ΔB| ≥ {min_delta_b:.1f} %B)   : {len(pairs)}")
    print(f"  f-tagged pairs        : {(pairs['Tag'] == 'f').sum()}")
    print(f"  F-tagged pairs        : {(pairs['Tag'] == 'F').sum()}")
    print(f"Saved → {output_path}")
    return pairs


if __name__ == "__main__":
    make_tag_pairs(
        input_path="curated_data/test.csv",
        output_path="curated_data/tag_pairs.csv",
        min_delta_b=1.0,
    )
