"""
Smoke-test that the StereoPep dataset loads correctly from Hugging Face Hub.
"""

from datasets import load_dataset


REPO_ID = "amirka20/StereoPep"
EXPECTED_COLUMNS = {"Peptide", "B", "M", "Z", "Length", "SMILES"}
EXPECTED_STEREO_COLUMNS = {"Sequence_f", "Sequence_F", "SMILES_f", "SMILES_F", "B_f", "B_F", "delta_B"}
EXPECTED_TAG_COLUMNS = {"Sequence_untagged", "Sequence_tagged", "Tag", "SMILES_untagged", "SMILES_tagged", "B_untagged", "B_tagged", "delta_B"}
EXPECTED_SUBSTITUTION_COLUMNS = {"Sequence_1", "Sequence_2", "Position", "Residue_1", "Residue_2", "SMILES_1", "SMILES_2", "B_1", "B_2", "delta_B"}


def check(condition: bool, msg: str) -> None:
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {msg}")


print("=" * 60)
print(f"Testing dataset: {REPO_ID}")
print("=" * 60)

# --- config: StereoPep ---
print("\nConfig 'StereoPep':")
ds = load_dataset(REPO_ID, "StereoPep")

for split, expected_rows in [("train", 41_455), ("val", 4_607), ("test", 2_726), ("trainval", 46_062)]:
    split_ds = ds[split]
    check(set(split_ds.column_names) == EXPECTED_COLUMNS, f"{split}: columns match")
    check(len(split_ds) == expected_rows, f"{split}: {len(split_ds):,} rows (expected {expected_rows:,})")

# --- config: diastereomer_pairs ---
print("\nConfig 'diastereomer_pairs':")
stereo = load_dataset(REPO_ID, "diastereomer_pairs")["diastereomer_pairs"]
check(set(stereo.column_names) == EXPECTED_STEREO_COLUMNS, "diastereomer_pairs: columns match")
check(len(stereo) == 542, f"diastereomer_pairs: {len(stereo):,} rows (expected 542)")

# --- config: terminal_tag_pairs ---
print("\nConfig 'terminal_tag_pairs':")
tag = load_dataset(REPO_ID, "terminal_tag_pairs")["terminal_tag_pairs"]
check(set(tag.column_names) == EXPECTED_TAG_COLUMNS, "terminal_tag_pairs: columns match")
check(len(tag) == 1_549, f"terminal_tag_pairs: {len(tag):,} rows (expected 1,549)")

# --- config: point_mutant_pairs ---
print("\nConfig 'point_mutant_pairs':")
sub = load_dataset(REPO_ID, "point_mutant_pairs")["point_mutant_pairs"]
check(set(sub.column_names) == EXPECTED_SUBSTITUTION_COLUMNS, "point_mutant_pairs: columns match")
check(len(sub) == 6_942, f"point_mutant_pairs: {len(sub):,} rows (expected 6,942)")

print("\nDone.")
