"""
cancer_cooccurrence_pipeline.py

Among HCC cases who are carriers of rare variants in each of the 6 cancer
genes (BRCA2, MSH6, PMS2, CHEK2, BRIP1, FANCA), tests whether they have
a co-occurring diagnosis of breast, ovarian, prostate, or pancreatic cancer
compared to non-carriers, using Fisher's exact test.

This script reads the carrier_status.tsv files already produced by
carrier_comorbidity_pipeline.py, run that script first.

Steps
-----
  1. Read carrier status files (one per gene x annotation category)
  2. Scan the condition file for cancer ICD codes in HCC cases
  3. Run Fisher's exact test for each gene x category x cancer type
  4. Write summary-level output CSV  

OUTPUT
------
  carrier_cancer_cooccurrence.csv
      One row per gene x annotation category x cancer type.
      Contains only counts, percentages, odds ratios, and p-values 

USAGE
-----
  python cancer_cooccurrence_pipeline.py \\
    --cases-file     /path/to/cases_only.txt \\
    --carrier-dir    /path/to/carrier_status_files/ \\
    --condition-file /path/to/condition_occurrence.txt \\
    --out-dir        /path/to/output_directory

  Run with --help to see all options.

NOTES
-----
  - --cases-file: two-column tab-separated file, no header; person_id is column 0.
    This is the cases_only.txt produced by carrier_comorbidity_pipeline.py.
  - --carrier-dir: directory containing {GENE}_{CATEGORY}_carrier_status.tsv files
    produced by carrier_comorbidity_pipeline.py. Run that script first.
  - --condition-file is assumed to be tab-separated with columns:
    [person_id, condition_source_value] containing ICD-9/ICD-10 codes.
    *** If your biobank uses a different format, adapt load_cancer_codes(). ***
"""

import argparse
import sys
import pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact


# ---------------------------------------------------------------------------
# ICD code definitions — edit here if needed
# ---------------------------------------------------------------------------

CANCER_CODES = {
    "Breast_Cancer": [
        "C50", "C50.0", "C50.1", "C50.2", "C50.3", "C50.4",
        "C50.5", "C50.6", "C50.8", "C50.9",
        "174", "174.0", "174.1", "174.2", "174.3", "174.4",
        "174.5", "174.6", "174.8", "174.9",
        "175", "175.0", "175.9",
    ],
    "Ovarian_Cancer": [
        "C56", "C56.1", "C56.2", "C56.9",
        "183", "183.0",
    ],
    "Prostate_Cancer": [
        "C61", "185",
    ],
    "Pancreatic_Cancer": [
        "C25", "C25.0", "C25.1", "C25.2", "C25.3", "C25.4",
        "C25.7", "C25.8", "C25.9",
        "157", "157.0", "157.1", "157.2", "157.3", "157.4",
        "157.8", "157.9",
    ],
}

GENES      = ["BRCA2", "MSH6", "PMS2", "CHEK2", "BRIP1", "FANCA"]
CATEGORIES = ["pLoF", "damaging_missense", "pLoF_damaging_missense"]


# ---------------------------------------------------------------------------
# ARGUMENT PARSING
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Cancer co-occurrence pipeline: tests whether HCC cases who carry "
            "rare variants in BRCA2/MSH6/PMS2/CHEK2/BRIP1/FANCA have enriched "
            "co-occurrence of breast, ovarian, prostate, or pancreatic cancer "
            "vs non-carriers. Outputs summary-level CSV only."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  python cancer_cooccurrence_pipeline.py \\\n"
            "    --cases-file     /path/to/cases_only.txt \\\n"
            "    --carrier-dir    /path/to/carrier_status_files/ \\\n"
            "    --condition-file /path/to/condition_occurrence.txt \\\n"
            "    --out-dir        /path/to/output_directory\n"
        ),
    )

    parser.add_argument(
        "--cases-file", required=True, type=Path,
        metavar="PATH",
        help=(
            "Two-column tab-separated file with no header; person_id is column 0. "
            "This is the cases_only.txt produced by carrier_comorbidity_pipeline.py."
        ),
    )
    parser.add_argument(
        "--carrier-dir", required=True, type=Path,
        metavar="PATH",
        help=(
            "Directory containing {GENE}_{CATEGORY}_carrier_status.tsv files "
            "produced by carrier_comorbidity_pipeline.py."
        ),
    )
    parser.add_argument(
        "--condition-file", required=True, type=Path,
        metavar="PATH",
        help=(
            "Tab-separated file with columns [person_id, condition_source_value] "
            "containing ICD-9/ICD-10 codes. "
            "Adapt load_cancer_codes() in the script if your format differs."
        ),
    )
    parser.add_argument(
        "--out-dir", required=True, type=Path,
        metavar="PATH",
        help="Directory where output files will be written (created if absent).",
    )
    parser.add_argument(
        "--genes",
        default="BRCA2,MSH6,PMS2,CHEK2,BRIP1,FANCA",
        metavar="GENE1,GENE2,...",
        help="Comma-separated list of genes to analyse. Default: BRCA2,MSH6,PMS2,CHEK2,BRIP1,FANCA",
    )
    parser.add_argument(
        "--categories",
        default="pLoF,damaging_missense,pLoF_damaging_missense",
        metavar="CAT1,CAT2,...",
        help=(
            "Comma-separated list of annotation categories to analyse. "
            "Default: pLoF,damaging_missense,pLoF_damaging_missense"
        ),
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# STEP 1: Load carrier status from pre-computed TSV files
# ---------------------------------------------------------------------------

def step1_load_carrier_status(genes, categories, carrier_dir):
    print("\n=== STEP 1: Loading carrier status files ===")

    carrier_status = {}   # (gene, category) -> {person_id: 0 or 1}

    for gene in genes:
        for category in categories:
            carrier_path = carrier_dir / f"{gene}_{category}_carrier_status.tsv"
            if not carrier_path.exists():
                print(f"  WARNING: {carrier_path.name} not found — skipping")
                carrier_status[(gene, category)] = {}
                continue

            df = pd.read_csv(carrier_path, sep="\t", dtype={"person_id": str})
            n_carriers = df["carrier"].sum()
            print(f"  {gene} | {category}: {n_carriers} carriers out of {len(df)} people")

            carrier_status[(gene, category)] = dict(
                zip(df["person_id"], df["carrier"])
            )

    return carrier_status


# ---------------------------------------------------------------------------
# STEP 2: Load cancer ICD codes for HCC cases
# ---------------------------------------------------------------------------

def load_cancer_codes(condition_file, case_ids, all_codes_upper):
    """
    Load cancer ICD codes for case individuals.

    Assumes tab-separated file with columns:
        person_id, condition_source_value

    *** Adapt this function if your biobank uses a different format. ***
    Returns: dict {person_id (str) -> set of ICD codes (upper-case str)}
    """
    person_codes = {}
    for chunk in pd.read_csv(
        condition_file, sep="\t",
        usecols=["person_id", "condition_source_value"],
        dtype=str, chunksize=500_000,
    ):
        chunk["person_id"] = chunk["person_id"].astype(str)
        chunk["condition_source_value"] = (
            chunk["condition_source_value"].str.strip().str.upper()
        )
        chunk = chunk[chunk["person_id"].isin(case_ids)]
        chunk = chunk[chunk["condition_source_value"].isin(all_codes_upper)]
        for pid, code in zip(chunk["person_id"], chunk["condition_source_value"]):
            person_codes.setdefault(pid, set()).add(code)
    return person_codes


# ---------------------------------------------------------------------------
# STEP 3: Co-occurrence counts + Fisher's exact test
# ---------------------------------------------------------------------------

def step3_cooccurrence(carrier_status, case_ids, condition_file,
                       cancer_codes, genes, categories, out_dir):
    print("\n=== STEP 2: Loading cancer ICD codes ===")
    all_codes_upper = {c.upper() for codes in cancer_codes.values() for c in codes}
    person_codes = load_cancer_codes(condition_file, case_ids, all_codes_upper)

    print("\n=== STEP 3: Cancer co-occurrence analysis ===")
    rows = []

    for gene in genes:
        for category in categories:
            status = carrier_status.get((gene, category), {})
            if not status:
                continue

            carrier_cases     = {pid for pid, c in status.items() if c == 1}
            non_carrier_cases = case_ids - carrier_cases

            print(f"\n  {gene} | {category}: "
                  f"{len(carrier_cases)} carriers, "
                  f"{len(non_carrier_cases)} non-carriers")

            for cancer, codes in cancer_codes.items():
                codes_upper = {c.upper() for c in codes}

                carrier_dx        = sum(1 for p in carrier_cases
                                        if person_codes.get(p, set()) & codes_upper)
                carrier_no_dx     = len(carrier_cases) - carrier_dx
                non_carrier_dx    = sum(1 for p in non_carrier_cases
                                        if person_codes.get(p, set()) & codes_upper)
                non_carrier_no_dx = len(non_carrier_cases) - non_carrier_dx

                if len(carrier_cases) > 0 and carrier_dx > 0 and non_carrier_dx > 0:
                    table = [[carrier_dx,     carrier_no_dx],
                             [non_carrier_dx, non_carrier_no_dx]]
                    or_val, pval = fisher_exact(table, alternative="two-sided")
                else:
                    or_val, pval = None, None

                rows.append({
                    "Gene":               gene,
                    "Category":           category,
                    "Cancer":             cancer,
                    "Carrier_n":          len(carrier_cases),
                    "Carrier_with_dx":    carrier_dx,
                    "Carrier_no_dx":      carrier_no_dx,
                    "Carrier_pct":        round(100 * carrier_dx / len(carrier_cases), 1)
                                          if carrier_cases else 0,
                    "NonCarrier_n":       len(non_carrier_cases),
                    "NonCarrier_with_dx": non_carrier_dx,
                    "NonCarrier_no_dx":   non_carrier_no_dx,
                    "NonCarrier_pct":     round(100 * non_carrier_dx / len(non_carrier_cases), 1)
                                          if non_carrier_cases else 0,
                    "OddsRatio":          round(or_val, 4) if or_val is not None else "NA",
                    "Fisher_pvalue":      round(pval, 4)   if pval   is not None else "NA",
                })

    results = pd.DataFrame(rows)
    out_path = out_dir / "carrier_cancer_cooccurrence.csv"
    results.to_csv(out_path, index=False)

    print("\n" + "=" * 80)
    print(results.to_string(index=False))
    print(f"\n>>> Summary file saved: {out_path}")
    print(">>> Please send this CSV file back for meta-analysis.")
    return results


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    genes      = [g.strip() for g in args.genes.split(",")]
    categories = [c.strip() for c in args.categories.split(",")]

    # Load HCC cases
    print("Loading HCC cases...")
    case_ids = set(
        pd.read_csv(args.cases_file, sep="\t", header=None, dtype=str)[0]
    )
    print(f"  {len(case_ids):,} cases loaded")

    # Step 1
    carrier_status = step1_load_carrier_status(
        genes       = genes,
        categories  = categories,
        carrier_dir = args.carrier_dir,
    )

    # Steps 2 + 3
    step3_cooccurrence(
        carrier_status = carrier_status,
        case_ids       = case_ids,
        condition_file = args.condition_file,
        cancer_codes   = CANCER_CODES,
        genes          = genes,
        categories     = categories,
        out_dir        = args.out_dir,
    )
