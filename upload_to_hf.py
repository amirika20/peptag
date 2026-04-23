"""
Upload StereoPep datasets to Hugging Face Hub.

Loads trainval, stereo_pairs, stereo_pairs_trainval, tag_pairs, substitution_pairs, test, train, and val CSVs from curated_data/.

Because stereo_pairs, tag_pairs, and substitution_pairs have different columns from the other splits,
each group is pushed as a separate named config (subset) on the same repo:

  - config "StereoPep"           → train / val / test / trainval splits
  - config "stereo_pairs"        → stereo_pairs / stereo_pairs_trainval splits
  - config "tag_pairs"           → tag_pairs split
  - config "substitution_pairs"  → substitution_pairs split
"""

import pandas as pd
from datasets import Dataset, DatasetDict
from huggingface_hub import login

REPO_ID = "amirka20/StereoPep"


def load_csv(path: str) -> Dataset:
    df = pd.read_csv(path)
    print(f"  {path}: {len(df):,} rows, columns: {list(df.columns)}")
    return Dataset.from_pandas(df, preserve_index=False)


if __name__ == "__main__":
    login()

    data_dir = "curated_data"

    # --- config: StereoPep (shared schema) ---
    print("\nLoading main splits ...")
    main_ds = DatasetDict({
        "train":  load_csv(f"{data_dir}/train.csv"),
        "val":    load_csv(f"{data_dir}/val.csv"),
        "test":   load_csv(f"{data_dir}/test.csv"),
        "trainval": load_csv(f"{data_dir}/trainval.csv"),
    })
    print(f"\nPushing config 'StereoPep' to {REPO_ID} ...")
    main_ds.push_to_hub(REPO_ID, config_name="StereoPep", private=False)

    # --- config: stereo_pairs (different schema) ---
    print("\nLoading stereo pairs ...")
    stereo_ds = DatasetDict({
        "stereo_pairs":          load_csv(f"{data_dir}/stereo_pairs.csv"),
        "stereo_pairs_trainval": load_csv(f"{data_dir}/stereo_pairs_trainval.csv"),
    })
    print(f"\nPushing config 'stereo_pairs' to {REPO_ID} ...")
    stereo_ds.push_to_hub(REPO_ID, config_name="stereo_pairs", private=False)

    # --- config: tag_pairs (different schema) ---
    print("\nLoading tag pairs ...")
    tag_ds = DatasetDict({
        "tag_pairs": load_csv(f"{data_dir}/tag_pairs.csv"),
    })
    print(f"\nPushing config 'tag_pairs' to {REPO_ID} ...")
    tag_ds.push_to_hub(REPO_ID, config_name="tag_pairs", private=False)

    # --- config: substitution_pairs (different schema) ---
    print("\nLoading substitution pairs ...")
    substitution_ds = DatasetDict({
        "substitution_pairs": load_csv(f"{data_dir}/substitution_pairs.csv"),
    })
    print(f"\nPushing config 'substitution_pairs' to {REPO_ID} ...")
    substitution_ds.push_to_hub(REPO_ID, config_name="substitution_pairs", private=False)

    print("\nDone.")
