"""
Upload peptag datasets to Hugging Face Hub.

Loads trainval, stereo_pairs, test, train, and val CSVs from curated_data/.

Because stereo_pairs has different columns from the other splits, each group
is pushed as a separate named config (subset) on the same repo:

  - config "peptag"        → train / val / test / trainval splits
  - config "stereo_pairs"  → stereo_pairs split
"""

import pandas as pd
from datasets import Dataset, DatasetDict
from huggingface_hub import login

REPO_ID = "amirka20/peptag"


def load_csv(path: str) -> Dataset:
    df = pd.read_csv(path)
    print(f"  {path}: {len(df):,} rows, columns: {list(df.columns)}")
    return Dataset.from_pandas(df, preserve_index=False)


if __name__ == "__main__":
    login()

    data_dir = "curated_data"

    # --- config: peptag (shared schema) ---
    print("\nLoading main splits ...")
    main_ds = DatasetDict({
        "train":  load_csv(f"{data_dir}/train.csv"),
        "val":    load_csv(f"{data_dir}/val.csv"),
        "test":   load_csv(f"{data_dir}/test.csv"),
        "trainval": load_csv(f"{data_dir}/trainval.csv"),
    })
    print(f"\nPushing config 'peptag' to {REPO_ID} ...")
    main_ds.push_to_hub(REPO_ID, config_name="peptag", private=False)

    # --- config: stereo_pairs (different schema) ---
    print("\nLoading stereo pairs ...")
    stereo_ds = DatasetDict({
        "stereo_pairs": load_csv(f"{data_dir}/stereo_pairs.csv"),
    })
    print(f"\nPushing config 'stereo_pairs' to {REPO_ID} ...")
    stereo_ds.push_to_hub(REPO_ID, config_name="stereo_pairs", private=False)

    print("\nDone.")
