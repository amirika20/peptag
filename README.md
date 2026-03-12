---
license: mit
task_categories:
  - tabular-regression
language:
  - en
tags:
  - peptide
  - retention-time
  - HPLC
  - mass-spectrometry
  - cheminformatics
  - SMILES
  - stereochemistry
pretty_name: PepTag — Peptide Retention Time Dataset
size_categories:
  - 10K<n<100K
---

# PepTag — Peptide Reversed-Phase HPLC Retention Time Dataset

PepTag is a curated dataset of **~46 000 unique synthetic peptides** with experimentally measured reversed-phase HPLC (RP-HPLC) retention times, expressed as the percentage of organic modifier (acetonitrile, %B) at elution. The dataset is designed for training and benchmarking **retention-time prediction models** and includes full SMILES with chirality centres, enabling cheminformatics and graph-neural-network approaches.

---

## Motivation

Predicting peptide retention time is a core problem in proteomics and peptidomics. Accurate RT prediction aids targeted mass-spectrometry workflows, library searching, and the design of synthetic peptide libraries. Most public RT datasets cover only standard L-amino-acid sequences of a single length. PepTag extends the landscape by including:

- Peptides of **6 – 15 residues** and variable-length libraries.
- Both **K-terminus** (Lys C-terminus) and **R-terminus** (Arg C-terminus) scaffolds.
- Three **N/C-terminal tagging variants**: D-Phenylalanine, L-Phenylalanine, and label-free (no tag).
- A dedicated **stereo-pair subset** pairing each D-Phe peptide with its L-Phe counterpart for studying the stereochemical effect on retention.
- Programmatically generated **isomeric SMILES** with explicit chirality at every alpha-carbon.

---

## Dataset at a Glance

| Split | File | Rows |
|---|---|---|
| Full merged | `curated_data/merged.csv` | 46 062 |
| Train (90 %) | `curated_data/train.csv` | 41 455 |
| Validation (10 %) | `curated_data/val.csv` | 4 607 |
| Test (held-out library) | `curated_data/test.csv` | 2 726 |
| Stereo pairs | `curated_data/stereo_pairs.csv` | 783 |

The test set comes from an independent library (`rawdata/test/`) that was not used to construct the train/val split, making it a true held-out benchmark.

---

## Experimental Setup

### Peptide Libraries

Peptides were synthesised as one-bead-one-compound (OBOC) combinatorial libraries. Each library targets a fixed sequence length or a variable-length design:

| Library | Length | Notes |
|---|---|---|
| 6mer, 6mer2 | 6 | Two independent synthesis batches |
| 8mer, 8mer2 | 8 | Two independent synthesis batches |
| 9mer | 9 | |
| 10mer | 10 | |
| 12mer | 12 | |
| 13mer | 13 | |
| 15mer | 15 | |
| varLen, varLen2 | Mixed | Variable-length; two batches |
| test | Mixed | Independent held-out benchmark library |

### Terminal Scaffolds

Two C-terminal residues anchor the library, determining how beads are cleaved and identified:

- **K-term (Lys C-terminus):** The peptide C-terminus is a lysine residue with a free alpha-carboxyl. An optional phenylalanine *tag* can be conjugated as an amide to the epsilon-amine of this Lys.
- **R-term (Arg C-terminus):** The peptide C-terminus is an arginine residue with a free alpha-carboxyl.

### Phenylalanine Tags

An N-terminal and/or C-terminal phenylalanine is used as a mass-spectrometry ionisation and identification tag. Three variants exist:

| Variant | Code | Description |
|---|---|---|
| D-Phenylalanine | `f` (lowercase) | Non-proteinogenic D-configured Phe at N-term (and appended to K epsilon-amine for K-term libs) |
| L-Phenylalanine | `F` (uppercase) | Proteinogenic L-configured Phe, same positions as D-Phe |
| Label-Free | `LF` | No phenylalanine tag; sequence starts with a random residue |

The `2plex` files (where both D- and L-Phe tags are present simultaneously) are excluded from the curated dataset because co-eluting stereo-isomers cannot be cleanly deconvolved.

### LC-MS Acquisition

Peptides were separated by RP-HPLC using a C18 column and identified by tandem mass spectrometry (LC-MS/MS). Raw PSM (peptide-spectrum match) files were exported from a database search engine (FragPipe / MSFragger) as TSV tables containing one row per spectrum-peptide match.

### HPLC Gradients

Three acetonitrile gradient methods were used across libraries. Raw retention times (in seconds) are converted to **%B at elution** — the percentage of acetonitrile in the mobile phase at the moment of peptide elution — using a piecewise-linear interpolation that accounts for the column dead volume (3.257 min):

| Method | Used for | Gradient (time min → %B) |
|---|---|---|
| method_1 | 6mer, 9mer | 0→6, 10→11, 40→21, 55→31, 60→61 |
| method_2 | 12mer; 15mer LF | 0→6, 10→11, 40→21, 55→31, 70→91 |
| method_3 | all others | 0→6, 10→11, 40→21, 55→31, 70→91, 80→91 |

Converting to %B makes retention values comparable across runs recorded under different gradient programmes.

---

## Curation Pipeline

The raw PSM files are processed by four scripts in sequence:

```
rawdata/
  <library>/psm/*.tsv
        │
        ▼  data_processing.py  (clean_psm)
curated_data/
  <library>/*.csv
        │
        ▼  merge_data.py  (merge_libraries)
curated_data/merged.csv   ←── train+val pool
curated_data/test.csv     ←── held-out test pool
        │
        ├─▶  split_data.py  →  train.csv, val.csv
        │
        └─▶  make_stereo_pairs.py  →  stereo_pairs.csv
```

### Step 1 — `data_processing.py`: PSM Cleaning

For each raw TSV file, `clean_psm()` performs the following:

1. **Column selection:** retains `Protein`, `Peptide`, `Retention`, `Intensity`, `Observed Mass`, `Delta Mass`, `Charge`.
2. **Quality filters:**
   - Remove PSMs with zero intensity.
   - Remove PSMs with `|ΔMass| ≥ 0.1 Da` (mass accuracy filter).
3. **Sequence validation:** filters peptides by expected terminal residues based on the file name (e.g. K-term D-Phe files must start with `F` and end with `KF`).
4. **D-Phe encoding:** converts the terminal phenylalanine residues to lowercase `f` for D-Phe libraries, distinguishing them from L-Phe (`F`).
5. **RT → %B conversion:** applies the gradient-appropriate piecewise-linear function; rows outside the gradient window are dropped.
6. **Deduplication per protein:** for peptides observed in multiple spectra, keeps the highest-intensity PSM.
7. **Output:** `Peptide`, `RT` (raw seconds), `B` (%B), `M` (observed mass), `Z` (charge state).

### Step 2 — `merge_data.py`: Library Merging and SMILES Generation

`merge_libraries()` concatenates all per-library curated CSVs (excluding `test`) into a single file, then:

1. **SMILES generation:** calls `peptide_to_smiles()` from `data_processing.py` for every unique sequence. SMILES include:
   - Explicit `[C@@H]` (S config) for L-amino acids at the alpha-carbon.
   - Explicit `[C@H]` (R config) for D-Phenylalanine (`f`).
   - Special handling for Glycine (no chiral centre), Proline (pyrrolidine ring), Isoleucine (two chiral centres), and Threonine (two chiral centres).
   - K-terminus Phe tag encoded as an amide branch on the Lys side chain.
2. **Cross-library deduplication:** if the same sequence appears in multiple libraries, rows are collapsed by averaging `B` and `M`, taking the mode of `Z`.
3. Adds a `Length` column (sequence length).

`test.csv` is produced by a separate call to `merge_libraries()` restricted to the `test` sub-directory.

### Step 3 — `split_data.py`: Stratified Train/Val Split

`split_data()` splits `merged.csv` into `train.csv` (90 %) and `val.csv` (10 %) using a **stratified split**. The stratification key is the combination of:

- **Terminal type:** `K` or `R`.
- **Phe type:** `D` (contains `f`), `L` (starts with `F`), or `LF` (neither).
- **Length:** sequence length.

This guarantees that every combination of terminal type, tagging variant, and length is represented proportionally in both splits. The random seed is fixed at 42 for reproducibility.

### Step 4 — `make_stereo_pairs.py`: Stereo Pair Extraction

`make_stereo_pairs()` extracts pairs of sequences that differ only in the stereochemistry of phenylalanine:

- Start from every D-Phe sequence in `test.csv`.
- Generate the putative L-Phe counterpart by replacing all `f` → `F`.
- Keep only pairs where both members exist in the dataset.

The resulting `stereo_pairs.csv` enables direct measurement of the **retention shift induced by phenylalanine stereoinversion**.

---

## Data Format

### `merged.csv` / `train.csv` / `val.csv` / `test.csv`

| Column | Type | Description |
|---|---|---|
| `Peptide` | string | One-letter amino acid sequence. Lowercase `f` = D-Phenylalanine; all other uppercase letters = standard L-amino acids. |
| `B` | float | Elution point expressed as % acetonitrile (%B) in the mobile phase. This is the primary regression target. |
| `M` | float | Observed monoisotopic mass (Da) from the mass spectrometer. |
| `Z` | int | Precursor ion charge state. |
| `Length` | int | Sequence length (number of residues). |
| `SMILES` | string | Isomeric SMILES with explicit stereochemistry at every alpha-carbon. |

### `stereo_pairs.csv`

| Column | Type | Description |
|---|---|---|
| `Sequence_f` | string | Peptide sequence containing D-Phe (`f`). |
| `Sequence_F` | string | Same sequence with all `f` replaced by `F` (L-Phe counterpart). |
| `B_f` | float | Elution %B for the D-Phe version. |
| `B_F` | float | Elution %B for the L-Phe version. |

---

## Amino Acid Alphabet

| Code | Residue | Notes |
|---|---|---|
| A | Alanine | |
| C | Cysteine | |
| D | Aspartate | |
| E | Glutamate | |
| F | L-Phenylalanine | Used as N/C-terminal ionisation tag |
| f | D-Phenylalanine | Non-proteinogenic; encoded lowercase |
| G | Glycine | No alpha-carbon chirality |
| H | Histidine | |
| I | L-Isoleucine | Two chiral centres: (2S,3S) |
| K | Lysine | C-terminus anchor for K-term libraries |
| L | Leucine | |
| M | Methionine | |
| N | Asparagine | |
| P | L-Proline | Pyrrolidine ring; special backbone topology |
| Q | Glutamine | |
| R | Arginine | C-terminus anchor for R-term libraries |
| S | Serine | |
| T | L-Threonine | Two chiral centres: (2S,3R) |
| V | Valine | |
| W | Tryptophan | |
| Y | Tyrosine | |

---

## Usage

### Loading with pandas

```python
import pandas as pd

train = pd.read_csv("curated_data/train.csv")
val   = pd.read_csv("curated_data/val.csv")
test  = pd.read_csv("curated_data/test.csv")

# Example: filter to 8-residue K-term label-free peptides
subset = train[(train["Length"] == 8) & train["Peptide"].str.endswith("K")]
print(subset.head())
```

### Loading from Hugging Face

```python
from datasets import load_dataset

ds = load_dataset("amirika20/peptag")
train = ds["train"].to_pandas()
```

### Reproducing SMILES for a sequence

```python
from data_processing import peptide_to_smiles

# L-amino acid peptide
print(peptide_to_smiles("ACDEFGHIK"))

# D-Phe tagged K-term peptide (starts and ends with f/Kf)
print(peptide_to_smiles("fACDEFGHIKf"))
```

### Reproducing the full pipeline

```bash
# 1. Curate raw PSMs
python data_processing.py

# 2. Merge libraries and generate SMILES
python merge_data.py

# 3. Split into train / val
python split_data.py

# 4. Extract stereo pairs from test set
python make_stereo_pairs.py
```

---

## Repository Structure

```
peptag/
├── rawdata/                    # Raw PSM TSV files from LC-MS/MS searches
│   ├── 6mer/psm/               # 6-residue library, batch 1
│   ├── 6mer2/psm/              # 6-residue library, batch 2
│   ├── 8mer/psm/               # 8-residue library, batch 1
│   ├── 8mer2/psm/              # 8-residue library, batch 2
│   ├── 9mer/psm/
│   ├── 10mer/psm/
│   ├── 12mer/psm/
│   ├── 13mer/psm/
│   ├── 15mer/psm/
│   ├── varLen/psm/             # Variable-length library, batch 1
│   ├── varLen2/psm/            # Variable-length library, batch 2
│   └── test/psm/               # Held-out benchmark library
│
├── curated_data/               # Processed outputs
│   ├── <library>/*.csv         # Per-library curated CSVs (intermediate)
│   ├── merged.csv              # All train+val peptides (deduplicated, with SMILES)
│   ├── train.csv               # 90% stratified training split
│   ├── val.csv                 # 10% stratified validation split
│   ├── test.csv                # Held-out test set (independent library)
│   └── stereo_pairs.csv        # D-Phe / L-Phe stereo pairs from test set
│
├── data_processing.py          # PSM cleaning + SMILES generation utilities
├── merge_data.py               # Library merging and deduplication
├── split_data.py               # Stratified train/val split
├── make_stereo_pairs.py        # Stereo-pair extraction
├── LICENSE
└── README.md
```

---

## Citation

If you use PepTag in your research, please cite:

```bibtex
@dataset{kazeminia2026peptag,
  author    = {Kazeminia, Amirabbas},
  title     = {PepTag: Peptide Reversed-Phase HPLC Retention Time Dataset},
  year      = {2026},
  publisher = {Hugging Face},
  url       = {https://huggingface.co/datasets/amirika20/peptag}
}
```

---

## License

This dataset and codebase are released under the [MIT License](LICENSE).
