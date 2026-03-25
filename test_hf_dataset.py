"""
Smoke-test that the peptag dataset loads correctly from Hugging Face Hub.
"""

from datasets import load_dataset


REPO_ID = "amirka20/peptag"
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

# --- config: peptag ---
print("\nConfig 'peptag':")
ds = load_dataset(REPO_ID, "peptag")

for split, expected_rows in [("train", 41_455), ("val", 4_607), ("test", 2_726), ("trainval", 46_062)]:
    split_ds = ds[split]
    check(set(split_ds.column_names) == EXPECTED_COLUMNS, f"{split}: columns match")
    check(len(split_ds) == expected_rows, f"{split}: {len(split_ds):,} rows (expected {expected_rows:,})")

# --- config: stereo_pairs ---
print("\nConfig 'stereo_pairs':")
stereo = load_dataset(REPO_ID, "stereo_pairs")["stereo_pairs"]
check(set(stereo.column_names) == EXPECTED_STEREO_COLUMNS, "stereo_pairs: columns match")
check(len(stereo) == 542, f"stereo_pairs: {len(stereo):,} rows (expected 542)")

# --- config: tag_pairs ---
print("\nConfig 'tag_pairs':")
tag = load_dataset(REPO_ID, "tag_pairs")["tag_pairs"]
check(set(tag.column_names) == EXPECTED_TAG_COLUMNS, "tag_pairs: columns match")
check(len(tag) == 1_549, f"tag_pairs: {len(tag):,} rows (expected 1,549)")

# --- config: substitution_pairs ---
print("\nConfig 'substitution_pairs':")
sub = load_dataset(REPO_ID, "substitution_pairs")["substitution_pairs"]
check(set(sub.column_names) == EXPECTED_SUBSTITUTION_COLUMNS, "substitution_pairs: columns match")
check(len(sub) == 6_942, f"substitution_pairs: {len(sub):,} rows (expected 6,942)")

print("\nDone.")
