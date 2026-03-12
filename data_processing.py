"""
data_processing.py
==================
Core utilities for the PepTag curation pipeline.

Two responsibilities:
  1. **SMILES generation** – convert a one-letter peptide sequence (including the
     non-proteinogenic D-Phenylalanine encoded as lowercase 'f') into an isomeric
     SMILES string with explicit alpha-carbon stereochemistry.

     SMILES are generated at the **pH 3 protonation state** matching the RP-HPLC
     mobile phase conditions (solvent A: water + 0.1% formic acid, pH ≈ 3).
     At pH 3 the following groups are protonated:
       - N-terminus alpha-amine  → [NH3+]   (pKa ~8)
       - Lys epsilon-amine       → [NH3+]   (pKa ~10.5)
       - Arg guanidinium         → [NH2+]=  (pKa ~12.5)
       - His imidazolium         → [nH+]    (pKa ~6)
       - Asp/Glu carboxyl        → COOH     (pKa ~3.7/4.1; neutral form in SMILES)
       - C-terminus carboxyl     → COOH     (neutral form in SMILES)
     Groups that are already neutral and unchanged: Cys thiol (–SH), Tyr phenol
     (–OH), Ser/Thr hydroxyl (–OH), peptide bond amides (–CO–NH–).

  2. **PSM cleaning** – read a raw PSM (peptide-spectrum match) TSV exported from
     a database search engine (FragPipe / MSFragger), apply quality filters, convert
     raw retention time to %B (% acetonitrile at elution), and write a clean CSV.

Run this file directly to curate all libraries under rawdata/:
    python data_processing.py
"""

import pandas as pd
from pathlib import Path


# ---------------------------------------------------------------------------
# SMILES generation
# ---------------------------------------------------------------------------

# Side chain SMILES fragments attached at the alpha-carbon branch point.
# These are written as if the alpha-C reads its substituents in the order:
#   backbone-N  →  alpha-C  →  side-chain  →  backbone-C(=O)
# so that the @/@@ stereo tokens produce the correct CIP configuration.
# Glycine and Proline are NOT in this dict; they are handled separately.
#
# Protonation state: pH 3 (RP-HPLC conditions).
#   K  → epsilon-NH3+   ([NH3+], pKa ~10.5 → fully protonated at pH 3)
#   R  → guanidinium    ([NH2+]=, pKa ~12.5 → fully protonated at pH 3)
#   H  → imidazolium    ([nH+], pKa ~6 → fully protonated at pH 3)
#   D/E → carboxylic acid (COOH neutral, pKa ~3.7/4.1 → mostly protonated at pH 3)
#   C  → thiol (–SH, neutral, pKa ~8.3 → protonated at pH 3)
#   Y  → phenol (–OH, neutral, pKa ~10.1 → protonated at pH 3)
_SIDE_CHAINS = {
    'A': 'C',
    'C': 'CS',
    'D': 'CC(=O)O',
    'E': 'CCC(=O)O',
    'F': 'Cc1ccccc1',
    'f': 'Cc1ccccc1',           # D-Phe — same side chain, D config at alpha-C
    'H': 'Cc1c[nH]c[nH+]1',    # imidazolium: N-delta neutral, N-epsilon protonated
    'I': '[C@@H](CC)C',         # L-Ile: beta-C is (S)
    'K': 'CCCC[NH3+]',          # epsilon-ammonium (protonated at pH 3)
    'L': 'CC(C)C',
    'M': 'CCSC',
    'N': 'CC(=O)N',
    'Q': 'CCC(=O)N',
    'R': 'CCCNC(=[NH2+])N',     # guanidinium (protonated at pH 3)
    'S': 'CO',
    'T': '[C@@H](O)C',          # L-Thr: beta-C is (R)
    'V': 'C(C)C',
    'W': 'Cc1c[nH]c2ccccc12',
    'Y': 'Cc1ccc(O)cc1',
}

# Alpha-carbon SMILES token for L vs D amino acids.
# f (D-Phe) → [C@H] (R config); all other non-Gly/Pro → [C@@H] (S config).
def _alpha_c(aa):
    if aa == 'G':
        return 'C'   # no chirality — CH2
    if aa == 'f':
        return '[C@H]'   # D-Phe → (R)
    return '[C@@H]'      # L-amino acid → (S)


def peptide_to_smiles(sequence: str) -> str:
    """Build an isomeric SMILES at pH 3 protonation state for a peptide sequence.

    Conventions used throughout this codebase:
      - Uppercase letters  = L-amino acids
      - 'f'                = D-Phenylalanine
      - Sequence ending in 'Kf' or 'KF'  → K-term Phe: K is the C-terminus
        (free alpha-COOH); the trailing f/F is attached as an amide to the
        epsilon-amine of K's side chain.
      - Sequence ending in 'K' (no trailing f/F) → K-term LF: normal K at
        C-terminus with free epsilon-NH3+.
      - Sequence ending in 'R' → R-term: linear peptide, free COOH at R.

    Protonation (pH 3, matching RP-HPLC mobile phase):
      - N-terminus alpha-amine → [NH3+]
      - Lys epsilon-amine      → [NH3+]
      - Arg guanidinium        → [NH2+]=
      - His imidazolium        → [nH+]
      - C-terminus             → COOH (already neutral in SMILES notation)
      - Asp / Glu              → COOH (already neutral in SMILES notation)
    """
    # Determine K-term Phe topology
    k_side = None
    if len(sequence) >= 2 and sequence[-2] == 'K' and sequence[-1] in ('f', 'F'):
        tag = sequence[-1]
        if tag == 'f':
            # D-Phe tag (R config at alpha-C).
            # The tag is read C→N (from-atom = C=O), which inverts @/@@
            # compared to the standard N→C backbone direction.
            # Verified: C(=O)[C@@H](Cc1ccccc1)N → from CIP, 1→2→3 CW = R = D-Phe.
            # Terminal amine of the Phe tag is [NH3+] at pH 3.
            k_side = 'CCCCNC(=O)[C@@H](Cc1ccccc1)[NH3+]'
        else:
            # L-Phe tag (S config at alpha-C).
            # Same reversal: C(=O)[C@H](Cc1ccccc1)N → S = L-Phe.
            k_side = 'CCCCNC(=O)[C@H](Cc1ccccc1)[NH3+]'
        core = sequence[:-1]   # drop the trailing Phe tag; K is now last
    else:
        core = sequence

    return _build_backbone_smiles(core, k_side)


def _build_backbone_smiles(sequence: str, k_side: str | None) -> str:
    """Assemble backbone SMILES from N-terminus to C-terminus at pH 3.

    The N-terminal nitrogen is protonated:
      - Primary amine (all residues except Pro): [NH3+]
      - Secondary amine (Pro at N-terminus):     [NH2+]
    Internal backbone amide nitrogens are NOT protonated (peptide bond pKa << 0).
    """
    n = len(sequence)
    parts = []

    for i, aa in enumerate(sequence):
        is_first = (i == 0)
        is_last = (i == n - 1)
        cterm = 'O' if is_last else ''   # append to C(=O) → C(=O)O at C-terminus

        if aa == 'G':
            # Glycine: no side chain, alpha-C is CH2 (no chirality).
            # N-terminus → [NH3+]; internal amide → N.
            nterm = '[NH3+]' if is_first else 'N'
            parts.append(f'{nterm}CC(=O){cterm}')

        elif aa == 'P':
            # Proline: pyrrolidine ring fuses the backbone N into the side chain.
            # N1 opens/closes the ring; [C@@H]1 is the alpha-C (L-Pro = S).
            # At N-terminus, Pro N is a secondary amine → [NH2+] when protonated.
            # In mid-chain, Pro N is a tertiary amide (no free H, stays as N).
            nterm = '[NH2+]' if is_first else 'N'
            parts.append(f'{nterm}1CCC[C@@H]1C(=O){cterm}')

        else:
            ac = _alpha_c(aa)
            # N-terminus → [NH3+]; internal amide → N
            nterm = '[NH3+]' if is_first else 'N'

            # Choose side chain: K at C-terminus gets optional Phe tag
            if aa == 'K' and is_last and k_side is not None:
                sc = k_side
            else:
                sc = _SIDE_CHAINS[aa]

            parts.append(f'{nterm}{ac}({sc})C(=O){cterm}')

    return ''.join(parts)


# ---------------------------------------------------------------------------
# Gradient RT → %B conversion
# ---------------------------------------------------------------------------
# Raw retention times (seconds) are converted to %B (% acetonitrile in the
# mobile phase at the moment of elution) using a piecewise-linear
# interpolation over the programmed HPLC gradient.
#
# The dead volume of 3.257 min is subtracted from the corrected RT before
# mapping onto the gradient table, so that t=0 corresponds to the first
# peptide entering the column.  Rows where the corrected RT falls outside
# the defined gradient window return None and are subsequently dropped.
#
# Each method corresponds to a distinct gradient programme used across
# different library batches (see README for the mapping).

_DEAD_VOLUME_MIN = 3.25666666666666  # column dead volume in minutes


def _rt_to_b(RT_seconds, time_min, gradient_pct):
    """Convert a raw retention time (seconds) to %B using a piecewise gradient.

    Parameters
    ----------
    RT_seconds   : retention time in seconds as reported by the search engine
    time_min     : list of gradient time knots (minutes after injection)
    gradient_pct : list of %B values at the corresponding knots

    Returns None if RT falls outside the defined gradient window.
    """
    corrected_rt = RT_seconds / 60 - _DEAD_VOLUME_MIN
    for i in range(len(time_min) - 1):
        if time_min[i] < corrected_rt < time_min[i + 1]:
            slope = (gradient_pct[i + 1] - gradient_pct[i]) / (time_min[i + 1] - time_min[i])
            return gradient_pct[i] + slope * (corrected_rt - time_min[i])
    return None


def method_1(RT):
    """Gradient method for 6mer and 9mer libraries.

    Gradient: 6→11→21→31→61 %B over 0→10→40→55→60 min.
    """
    return _rt_to_b(RT, [0, 10, 40, 55, 60], [6, 11, 21, 31, 61])


def method_2(RT):
    """Gradient method for 12mer and 15mer-LF libraries.

    Gradient: 6→11→21→31→91 %B over 0→10→40→55→70 min.
    """
    return _rt_to_b(RT, [0, 10, 40, 55, 70], [6, 11, 21, 31, 91])


def method_3(RT):
    """Gradient method for all remaining libraries (default).

    Gradient: 6→11→21→31→91→91 %B over 0→10→40→55→70→80 min.
    The plateau at 91 %B from 70–80 min ensures all late-eluting peptides
    are flushed before re-equilibration.
    """
    return _rt_to_b(RT, [0, 10, 40, 55, 70, 80], [6, 11, 21, 31, 91, 91])


def should_skip(file_name):
    """Return True if the file should not be curated.

    Excluded file types:
    - '2plex' / '2-plex': both D-Phe and L-Phe tags were present in the same
      run.  Co-eluting stereo-isomers cannot be cleanly resolved by RT alone.
    - '13C': isotope-labelled variants used for quantification, not RT measurement.
    """
    return "2plex" in file_name or "2-plex" in file_name or "13C" in file_name


def get_method(lib_name, file_name):
    """Return the RT→%B conversion function appropriate for this library/file.

    Mapping:
      6mer, 9mer          → method_1  (shorter gradient)
      12mer               → method_2  (longer gradient)
      15mer + LF          → method_2
      everything else     → method_3  (longest gradient, with 91 %B plateau)
    """
    if lib_name in ("6mer", "9mer"):
        return method_1
    if lib_name == "12mer":
        return method_2
    if lib_name == "15mer" and "LF" in file_name:
        return method_2
    return method_3


def clean_psm(input_path, output_path, method_fn):
    """Read a raw PSM file, apply quality filters, and write a clean CSV.

    Parameters
    ----------
    input_path : str or Path
        Raw PSM file (TSV or CSV) exported from FragPipe / MSFragger.
    output_path : str or Path
        Destination CSV file.
    method_fn : callable
        One of method_1 / method_2 / method_3; converts RT (seconds) → %B.

    Returns
    -------
    int
        Number of rows written to output_path.

    Processing steps
    ----------------
    1. Load PSM file and retain only the columns needed for downstream use.
    2. Quality filters: positive intensity AND |ΔMass| < 0.1 Da.
    3. Sequence validation: require the expected terminal residues inferred
       from the file name (K-term vs R-term, D-Phe vs L-Phe vs LF).
    4. D-Phe encoding: replace the terminal phenylalanine characters with
       lowercase 'f' so D-Phe ('f') and L-Phe ('F') are distinguishable.
    5. RT → %B conversion; rows outside the gradient window are dropped.
    6. Per-protein deduplication: keep the highest-intensity PSM.
    7. Write sorted output CSV with columns: Peptide, RT, B, M, Z.
    """
    input_path = str(input_path)

    if input_path.endswith(".tsv"):
        df = pd.read_csv(input_path, sep="\t")
    elif input_path.endswith(".csv"):
        df = pd.read_csv(input_path, sep=",")
    else:
        raise ValueError(f"Unsupported file format: {input_path}")

    # Step 1: retain only the columns required for downstream processing
    df = df[["Protein", "Peptide", "Retention", "Intensity", "Observed Mass", "Delta Mass", "Charge"]]

    # Step 2: quality filters
    df = df[df["Intensity"] > 0]           # exclude zero-intensity spectra
    df = df[abs(df["Delta Mass"]) < 0.1]   # exclude mis-assigned PSMs (>0.1 Da mass error)

    # Infer labelling scheme from the file name
    is_k_term = "K-term" in input_path
    is_r_term = "R-term" in input_path
    is_d_phe = "D-Phe" in input_path
    is_l_phe = "L-Phe" in input_path
    is_lf = "LF" in input_path

    # Step 3: sequence validation – keep only correctly terminated peptides.
    # Phe-tagged libraries: N-term F (or f) + expected C-terminus
    if is_d_phe or is_l_phe:
        df = df[df["Peptide"].str[0] == "F"]   # N-terminal tag (D or L config — both stored as 'F' in PSMs)
        if is_k_term:
            df = df[df["Peptide"].str[-2:] == "KF"]   # C-terminal: ...K-Phe amide
        elif is_r_term:
            df = df[df["Peptide"].str[-1] == "R"]     # C-terminal: ...R
    elif is_lf:
        # Label-free: no Phe tag, but still anchor-terminated
        if is_k_term:
            df = df[df["Peptide"].str[-1] == "K"]
        elif is_r_term:
            df = df[df["Peptide"].str[-1] == "R"]

    # Step 4: D-Phe encoding.
    # The mass spectrometer cannot distinguish D-Phe from L-Phe; both appear
    # as 'F' in the PSM file.  We re-encode terminal 'F' residues as 'f'
    # (D-Phenylalanine) based on the file name label.
    if is_d_phe and is_k_term:
        # Pattern: F<core>KF  →  f<core>Kf
        df = df.copy()
        df["Peptide"] = df["Peptide"].apply(lambda x: "f" + x[1:-1] + "f")
    elif is_d_phe and is_r_term:
        # Pattern: F<core>R  →  f<core>R
        df = df.copy()
        df["Peptide"] = df["Peptide"].apply(lambda x: "f" + x[1:])

    # Step 5: RT → %B conversion; rows outside gradient window become NaN and are dropped
    df["B"] = df["Retention"].apply(method_fn)
    df = df.dropna(subset=["B"])

    # Step 6: per-protein deduplication (each "Protein" ID maps to one unique peptide)
    # Retain only the highest-intensity observation for each protein.
    df = df.loc[df.groupby("Protein")["Intensity"].idxmax()]

    # Sanity check: warn if two different proteins share the same peptide sequence
    # (should not happen, but flags potential library contamination or search errors)
    duplicated_peptides = df[df["Peptide"].duplicated(keep=False)]["Peptide"].unique()
    if len(duplicated_peptides) > 0:
        print(f"  WARNING: peptide sequence shared across multiple proteins in {input_path}:")
        for pep in duplicated_peptides:
            print(f"    {pep}")

    # Step 7: write clean output
    df = df.sort_values(by="Protein")
    df = df[["Peptide", "Retention", "B", "Observed Mass", "Charge"]]
    df = df.rename(columns={"Retention": "RT", "Observed Mass": "M", "Charge": "Z"})
    df.to_csv(output_path, index=False)
    return len(df)


if __name__ == "__main__":
    base_path = Path("/home/amirabbas-kazeminia/Projects/peptag")
    rawdata_path = base_path / "rawdata"
    curated_path = base_path / "curated_data"

    lib_dirs = sorted([d for d in rawdata_path.iterdir() if d.is_dir()])

    for lib_dir in lib_dirs:
        lib_name = lib_dir.name
        psm_dir = lib_dir / "psm"

        if not psm_dir.exists():
            print(f"Skipping {lib_name}: no psm/ subdirectory found.")
            continue

        output_dir = curated_path / lib_name
        output_dir.mkdir(parents=True, exist_ok=True)

        input_files = sorted([f for f in psm_dir.iterdir() if f.is_file()])
        print(f"\n[{lib_name}] {len(input_files)} file(s) found")

        for file in input_files:
            if should_skip(file.name):
                print(f"  SKIP  {file.name}")
                continue

            method_fn = get_method(lib_name, file.name)
            output_file = output_dir / (file.stem + ".csv")

            try:
                n = clean_psm(str(file), str(output_file), method_fn)
                print(f"  OK    {file.name}  ({method_fn.__name__}, {n} rows -> {output_file.name})")
            except Exception as e:
                print(f"  ERROR {file.name}: {e}")
