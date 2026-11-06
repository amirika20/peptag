"""
make_croissant_rai.py
=====================
Rebuild the two Croissant metadata files from the live Hugging Face dataset.

Hugging Face generates a Croissant record describing the dataset's files, configs
and field types. That generated record is authoritative for structure, but it has
no place for Responsible AI documentation, so the RAI fields are maintained here
and merged on top of a freshly fetched copy.

Two files are written:

  croissant_metadata.json       the generated record, unmodified
  croissant_rai_StereoPep.json  the same record plus rai:* and prov:* fields

Regenerate whenever the dataset's configs or splits change, so the metadata does
not drift from what is actually published.

Usage:
    python make_croissant_rai.py                 # fetch and write both files
    python make_croissant_rai.py --dry-run       # report, write nothing
    python make_croissant_rai.py --repo other/ds
"""

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_REPO = "amirka20/StereoPep"

PLAIN_OUT = REPO_ROOT / "croissant_metadata.json"
RAI_OUT = REPO_ROOT / "croissant_rai_StereoPep.json"

CURATION_REPO = "https://github.com/amirika20/StereoPep-Curation"

# JSON-LD prefixes the RAI and provenance terms need.
EXTRA_CONTEXT = {
    "rai": "http://mlcommons.org/croissant/RAI/",
    "prov": "http://www.w3.org/ns/prov#",
}

RAI_FIELDS = {
    "rai:dataLimitations": (
        "This dataset was designed and generated in a controlled chemistry laboratory "
        "setting using a specific set of synthetic peptide libraries, and therefore "
        "reflects the design choices and experimental constraints of those libraries "
        "rather than the full diversity of peptide chemical space. Stereochemical "
        "variation is introduced through D-/L-phenylalanine tagging at the N-terminus "
        "and, for lysine-terminated sequences, at the side-chain N-epsilon-amine; "
        "stereocenters at internal positions are not varied, so stereochemical coverage "
        "is incomplete. The data are reversed-phase liquid chromatography retention-time "
        "measurements collected under a single column and gradient configuration, "
        "expressed as percent acetonitrile at elution, and may not generalize to other "
        "chromatographic systems or assays. The dataset was constructed for developing "
        "and benchmarking machine learning models, particularly for peptide and molecular "
        "representation learning and for evaluating sensitivity to stereochemical "
        "perturbations. It is not recommended for downstream applications requiring broad "
        "biochemical generalization, such as drug discovery pipelines, without further "
        "validation, nor for tasks that depend on accurate three-dimensional structural "
        "modeling or physical simulation. It should not be used for clinical, medical, or "
        "safety-critical decision-making."
    ),
    "rai:dataBiases": (
        "The dataset exhibits selection bias due to its construction from a predefined set "
        "of synthetic peptide libraries, leading to overrepresentation of specific sequence "
        "motifs, amino acid compositions, and stereochemical configurations that arise from "
        "the design process described in the accompanying paper. Sequences terminate in "
        "arginine or lysine to mimic tryptic digestion products, so the composition "
        "distribution differs from peptides observed in shotgun proteomics. Other regions of "
        "peptide space, including diverse modifications and broader stereochemical variants, "
        "are underrepresented or absent. Label bias may be present due to experimental noise "
        "and systematic effects in chromatography measurements, such as gradient conditions "
        "and instrument calibration. The paired benchmarks emphasize relative retention-time "
        "differences within matched peptide pairs, which may bias models toward learning "
        "local comparative signals rather than globally generalizable patterns. In the "
        "held-out diastereomeric pair benchmark, 64.0 percent of pairs have the D-phenylalanine "
        "form eluting later, so pairwise accuracy carries a majority-class baseline of 0.640 "
        "and threshold-free metrics such as AUC should be preferred. Models trained on this "
        "dataset may perform well on similar synthetic distributions but fail to generalize "
        "to natural peptides or unseen chemical modifications."
    ),
    "rai:personalSensitiveInformation": (
        "This dataset does not contain any personal or sensitive information. It does not "
        "include data related to gender, socio-economic status, geography, language, age, "
        "culture, experience or seniority, health or medical information, or political or "
        "religious beliefs, as it consists entirely of chemical measurements on synthetic "
        "peptides produced in a laboratory setting."
    ),
    "rai:dataUseCases": (
        "This dataset is intended to represent peptide physicochemical behavior, specifically "
        "chromatographic retention time as a proxy for hydrophobicity, as well as "
        "stereochemical effects arising from differences such as D- versus L-amino acid "
        "substitutions. It is designed to evaluate protein and molecular models for "
        "representation learning and their ability to capture stereochemistry-dependent "
        "changes in measurable properties. Alongside retention-time regression, it provides "
        "three structured pair benchmarks covering matched diastereomers, single-residue "
        "substitutions, and terminal tag effects, and a canonical-amino-acid-only control "
        "subset. The dataset has been validated for sequence-based prediction tasks involving "
        "retention time and pairwise stereochemical comparisons, as demonstrated in the "
        "accompanying experiments. It has not been validated for three-dimensional structural "
        "modeling, biological function prediction such as binding affinity or toxicity, or "
        "applications requiring in vivo relevance, and its use in such contexts may not be "
        "appropriate."
    ),
    "rai:dataSocialImpact": (
        "The dataset has the potential to positively impact research by enabling more rigorous "
        "evaluation of machine learning models in chemistry, particularly in capturing subtle "
        "stereochemical effects, and may contribute to the development of more robust and "
        "generalizable models for molecular design. However, misuse of the dataset without "
        "accounting for its limitations could lead to incorrect scientific conclusions or poor "
        "model generalization, especially if applied outside its intended domain. There is also "
        "a risk that over-reliance on this dataset could reinforce narrow modeling assumptions "
        "tied to synthetic laboratory conditions. These risks are mitigated by clearly "
        "documenting the dataset's scope, biases, and limitations, and by encouraging users to "
        "validate models on complementary datasets and broader settings before applying them in "
        "real-world scenarios."
    ),
    # The peptides are chemically synthesized, but every label is a real physical
    # measurement, so no records are synthetically generated in the ML sense.
    "rai:hasSyntheticData": False,
    "prov:wasGeneratedBy": [
        {
            "@type": "prov:Activity",
            "prov:type": {"@id": "https://www.wikidata.org/wiki/Q4929239"},
            "prov:label": "fast flow peptide synthesis",
        }
    ],
    "prov:wasDerivedFrom": [
        {
            "@id": CURATION_REPO,
            "prov:label": "StereoPep curation pipeline",
            "sc:license": "MIT License",
        }
    ],
}

# The Hub API endpoint returns the complete record: the structural part
# (distribution + recordSet) plus the descriptive header (creator, url, license,
# keywords, identifier). The datasets-server "croissant-crumbs" endpoint returns
# only the structural part, so it is not used here.
CROISSANT_URL = "https://huggingface.co/api/datasets/{}/croissant"


def _get(url: str) -> dict:
    """GET a JSON body, returning the parsed payload even on an HTTP error.

    While the dataset viewer is reindexing, Hugging Face answers with HTTP 500 and
    a JSON body explaining that the response is not ready. That is a wait-and-retry
    condition, not a failure, so the body has to be read rather than raised past.
    """
    req = urllib.request.Request(url, headers={"User-Agent": "stereopep-croissant/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:  # noqa: S310
            return json.load(resp)
    except urllib.error.HTTPError as exc:
        try:
            return json.load(exc)
        except Exception:
            return {"error": f"HTTP {exc.code} {exc.reason}"}


def fetch(repo: str, attempts: int = 10, delay: int = 30) -> dict:
    url = CROISSANT_URL.format(repo)

    for attempt in range(1, attempts + 1):
        data = _get(url)

        if data.get("recordSet"):
            return data

        reason = data.get("error") or "response has no recordSet yet"
        if attempt == attempts:
            print(f"ERROR: gave up after {attempts} attempts: {reason}", file=sys.stderr)
            print("The dataset viewer is still indexing. Retry in a few minutes.", file=sys.stderr)
            sys.exit(1)

        print(f"  attempt {attempt}/{attempts}: {reason}; retrying in {delay}s ...")
        time.sleep(delay)

    raise AssertionError("unreachable")


def build_rai(base: dict) -> dict:
    doc = json.loads(json.dumps(base))  # deep copy

    ctx = doc.get("@context")
    if not isinstance(ctx, dict):
        print("ERROR: unexpected @context shape; cannot merge RAI terms.", file=sys.stderr)
        sys.exit(1)
    for prefix, iri in EXTRA_CONTEXT.items():
        ctx.setdefault(prefix, iri)

    doc.update(RAI_FIELDS)
    return doc


def summarize(doc: dict, label: str) -> None:
    configs = [r.get("@id") or r.get("name") for r in doc.get("recordSet", [])]
    configs = [c for c in configs if c]
    print(f"  {label}")
    print(f"    creator    : {doc.get('creator', {}).get('name')}")
    print(f"    url        : {doc.get('url')}")
    print(f"    recordSets : {len(doc.get('recordSet', []))}")
    for c in configs:
        print(f"       - {c}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("Usage:")[0].strip())
    parser.add_argument("--repo", default=DEFAULT_REPO)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    print(f"Fetching Croissant for {args.repo} ...")
    base = fetch(args.repo)
    rai = build_rai(base)

    print()
    summarize(base, f"{PLAIN_OUT.name} (generated record, unmodified)")
    print()
    summarize(rai, f"{RAI_OUT.name} (+ RAI and provenance)")
    added = [k for k in RAI_FIELDS if k in rai]
    print(f"    added      : {len(added)} fields -> {', '.join(added)}")

    if args.dry_run:
        print("\n--dry-run set: no files written.")
        return

    PLAIN_OUT.write_text(json.dumps(base, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    RAI_OUT.write_text(json.dumps(rai, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"\nWrote {PLAIN_OUT.name} ({PLAIN_OUT.stat().st_size:,} bytes)")
    print(f"Wrote {RAI_OUT.name} ({RAI_OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
