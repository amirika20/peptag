"""
Upload StereoPep datasets to Hugging Face Hub.

Loads trainval, diastereomer_pairs, diastereomer_pairs_trainval, terminal_tag_pairs, point_mutant_pairs, test, train, and val CSVs from curated_data/.

Because diastereomer_pairs, terminal_tag_pairs, and point_mutant_pairs have different columns from the other splits,
each group is pushed as a separate named config (subset) on the same repo:

  - config "StereoPep"            → train / val / test / trainval splits
  - config "diastereomer_pairs"   → diastereomer_pairs / diastereomer_pairs_trainval splits
  - config "terminal_tag_pairs"   → terminal_tag_pairs split
  - config "point_mutant_pairs"   → point_mutant_pairs split
"""

import os
import pandas as pd
from datasets import Dataset, DatasetDict
from huggingface_hub import login

REPO_ID = "stereopep-ano/stereopep"
HF_TOKEN = os.environ.get("HF_TOKEN")  # set via: export HF_TOKEN=hf_...


def load_csv(path: str) -> Dataset:
    df = pd.read_csv(path)
    print(f"  {path}: {len(df):,} rows, columns: {list(df.columns)}")
    return Dataset.from_pandas(df, preserve_index=False)


if __name__ == "__main__":
    login(token=HF_TOKEN)

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

    # --- config: diastereomer_pairs (different schema) ---
    print("\nLoading diastereomer pairs ...")
    stereo_ds = DatasetDict({
        "diastereomer_pairs":          load_csv(f"{data_dir}/diastereomer_pairs.csv"),
        "diastereomer_pairs_trainval": load_csv(f"{data_dir}/diastereomer_pairs_trainval.csv"),
    })
    print(f"\nPushing config 'diastereomer_pairs' to {REPO_ID} ...")
    stereo_ds.push_to_hub(REPO_ID, config_name="diastereomer_pairs", private=False)

    # --- config: terminal_tag_pairs (different schema) ---
    print("\nLoading terminal tag pairs ...")
    tag_ds = DatasetDict({
        "terminal_tag_pairs": load_csv(f"{data_dir}/terminal_tag_pairs.csv"),
    })
    print(f"\nPushing config 'terminal_tag_pairs' to {REPO_ID} ...")
    tag_ds.push_to_hub(REPO_ID, config_name="terminal_tag_pairs", private=False)

    # --- config: point_mutant_pairs (different schema) ---
    print("\nLoading point mutant pairs ...")
    substitution_ds = DatasetDict({
        "point_mutant_pairs": load_csv(f"{data_dir}/point_mutant_pairs.csv"),
    })
    print(f"\nPushing config 'point_mutant_pairs' to {REPO_ID} ...")
    substitution_ds.push_to_hub(REPO_ID, config_name="point_mutant_pairs", private=False)

    print("\nDone.")
