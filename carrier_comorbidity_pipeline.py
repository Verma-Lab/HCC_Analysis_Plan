"""
carrier_comorbidity_pipeline.py

Pipeline to test whether carriers of rare variants in BRCA2, MSH6, PMS2,
CHEK2, BRIP1, or FANCA are enriched for liver disease comorbidities
compared to non-carriers, within a defined case population.

Steps
-----
  1. Parse SAIGE-GENE+ group files -> per-gene x annotation variant lists
  2. Run PLINK2 to extract carrier genotypes (calls plink2 via subprocess)
  3. Build binary carrier status and run Fisher's exact enrichment test
  4. Write summary-level output CSV  

OUTPUT
------
  carrier_liver_disease_enrichment.csv
      One row per gene x annotation category x disease.
      Contains only counts, percentages, odds ratios, and p-values 

USAGE
-----
  python carrier_comorbidity_pipeline.py \\
    --brca2-group-file /path/to/cancer_genes.13.txt \\
    --msh6-group-file  /path/to/cancer_genes.2.txt \\
    --pms2-group-file  /path/to/cancer_genes.7.txt \\
    --chek2-group-file /path/to/cancer_genes.22.txt \\
    --brip1-group-file /path/to/cancer_genes.17.txt \\
    --fanca-group-file /path/to/cancer_genes.16.txt \\
    --bfile            /path/to/plink_bfile_prefix \\
    --covariate-file   /path/to/covariates.csv \\
    --case-column      HCC_Cancer_Free \\
    --condition-file   /path/to/condition_occurrence.txt \\
    --out-dir          /path/to/output_directory

  Run with --help to see all options.

NOTES
-----
  - PLINK2 must be available on your PATH (use --plink2-bin if it is elsewhere).
  - --covariate-file must be a CSV with columns: person_id, <case_column>
    where 1 = case and 0 = control.
  - --condition-file is assumed to be tab-separated with columns:
    [person_id, condition_source_value] containing ICD-9/ICD-10 codes.
    *** If your biobank uses a different format, adapt load_condition_codes(). ***
"""

import argparse
import subprocess
import sys
import math
from collections import defaultdict
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact


# ---------------------------------------------------------------------------
# ICD code definitions 
# ---------------------------------------------------------------------------

DISEASE_CODES = {
    "HBV": [
        "B18.0", "B18.1",
        "070.32", "070.21", "070.22", "070.23", "070.31", "070.33",
    ],
    "HCV": [
        "B18.2", "B19.2", "B17.1",
        "070.54", "070.41", "070.44", "070.5", "070.7",
    ],
    "Alcohol_Related_Liver_Disease": [
        "571.0", "K70.0", "571.1", "K70.1", "571.2", "K70.3", "K70.2",
        "571.3", "K70.4", "K70.40", "K70.41", "K70.9",
        "303.0", "303.9", "F10.229", "F10.20",
    ],
    "MASH_MASLD": [
        "K76.0", "K75.81", "K75.8", "K75.89",
    ],
}

GENES_OF_INTEREST = ["BRCA2", "MSH6", "PMS2", "CHEK2", "BRIP1", "FANCA"]


# ---------------------------------------------------------------------------
# STEP 1: Parse group files -> variant lists + cases keep file
# ---------------------------------------------------------------------------

def normalize_variant_id(var: str) -> str:
    """Convert colon separators to underscores to match BIM format.
    e.g. chr2:12345:A:T -> chr2_12345_A_T
    """
    return var.replace(":", "_")


def step1_parse_group_files(covariate_file, case_column, group_files,
                             genes_of_interest, categories, out_dir,
                             plof_label="pLoF", damaging_label="damaging_missense"):
    print("\n=== STEP 1: Parsing group files ===")

    gene_cat_variants = defaultdict(lambda: defaultdict(list))

    for gene, group_file in group_files.items():
        if not group_file.exists():
            print(f"  ERROR: Group file not found: {group_file}")
            sys.exit(1)

        lines = group_file.read_text().splitlines()
        i = 0
        while i < len(lines) - 1:
            var_parts  = lines[i].split()
            anno_parts = lines[i + 1].split()

            if var_parts[0] == gene and anno_parts[0] == gene:
                variants    = var_parts[1:]
                annotations = anno_parts[1:]
                for var, anno in zip(variants, annotations):
                    if anno in categories:
                        gene_cat_variants[gene][anno].append(
                            normalize_variant_id(var)
                        )
                i += 2
            else:
                i += 1

    # Write per-gene x category variant lists
    variant_files = {}
    for gene in genes_of_interest:
        cats     = gene_cat_variants.get(gene, {})
        plof     = cats.get(plof_label, [])
        dm       = cats.get(damaging_label, [])
        combined = list(dict.fromkeys(plof + dm))  # union, order-preserving

        for label, variants in [
            (plof_label,                              plof),
            (damaging_label,                          dm),
            (f"{plof_label}_{damaging_label}",        combined),
        ]:
            out_path = out_dir / f"{gene}_{label}.txt"
            out_path.write_text("\n".join(variants) + "\n")
            variant_files[(gene, label)] = out_path
            print(f"  {out_path.name}: {len(variants)} variants")

    # Write PLINK keep file — cases only
    cov = pd.read_csv(covariate_file, usecols=["person_id", case_column])
    cov["person_id"] = cov["person_id"].astype(str)
    cases = cov.loc[cov[case_column] == 1, "person_id"]

    keep_path = out_dir / "cases_only.txt"
    keep_path.write_text(
        "\n".join(f"{pid}\t{pid}" for pid in cases) + "\n"
    )
    print(f"  cases_only.txt: {len(cases):,} cases")

    return variant_files, set(cases)


# ---------------------------------------------------------------------------
# STEP 2: PLINK2 extraction
# ---------------------------------------------------------------------------

def step2_plink_extract(plink_bfile, variant_files, out_dir, plink2_bin):
    print("\n=== STEP 2: PLINK2 extraction ===")
    keep_file = out_dir / "cases_only.txt"

    for (gene, category), var_file in variant_files.items():
        out_prefix = out_dir / f"{gene}_{category}_carriers"
        raw_file   = Path(str(out_prefix) + ".raw")

        if raw_file.exists():
            print(f"  {raw_file.name} already exists — skipping")
            continue

        cmd = [
            plink2_bin,
            "--bfile",   str(plink_bfile),
            "--extract", str(var_file),
            "--keep",    str(keep_file),
            "--export",  "A",
            "--out",     str(out_prefix),
        ]
        print(f"  Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"  ERROR running plink2:\n{result.stderr}")
            sys.exit(1)
        print(f"  -> {raw_file.name}")


# ---------------------------------------------------------------------------
# STEP 3: Build carrier status + Fisher's exact enrichment test
# ---------------------------------------------------------------------------

def load_condition_codes(condition_file, case_ids, all_codes_upper):
    """
    Load ICD codes for case individuals.

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


def build_carrier_status(gene, category, out_dir):
    """Returns dict {person_id (str) -> carrier (0 or 1)}, or {} if file missing."""
    raw_path = out_dir / f"{gene}_{category}_carriers.raw"
    if not raw_path.exists():
        print(f"  WARNING: {raw_path.name} not found — skipping")
        return {}

    raw = pd.read_csv(raw_path, sep=r"\s+", dtype=str)
    genotype_cols = raw.columns[6:]
    print(f"  {gene} | {category}: {len(genotype_cols)} variant columns")

    # PLINK2 --export A counts REF alleles by default.
    # REF allele count < 2 means at least one ALT allele -> carrier
    records = []
    for _, row in raw.iterrows():
        pid = str(row.iloc[1])
        is_carrier = 0
        for c in genotype_cols:
            val = row[c]
            if val in ("NA", "nan") or val != val:
                continue
            try:
                if float(val) < 2:
                    is_carrier = 1
                    break
            except ValueError:
                pass
        records.append({"person_id": pid, "carrier": is_carrier})

    df = pd.DataFrame(records)
    n_carriers = df["carrier"].sum()
    print(f"    -> {n_carriers} carriers out of {len(df)} people")

    carrier_path = out_dir / f"{gene}_{category}_carrier_status.tsv"
    df.to_csv(carrier_path, sep="\t", index=False)

    return dict(zip(df["person_id"], df["carrier"]))


def step3_enrichment(genes, categories, case_ids, condition_file,
                     disease_codes, out_dir):
    print("\n=== STEP 3: Carrier enrichment tests ===")

    combined_label  = "_".join(categories)
    all_categories  = list(categories) + [combined_label]
    all_codes_upper = {c.upper() for codes in disease_codes.values() for c in codes}

    print("  Loading condition codes...")
    person_codes = load_condition_codes(condition_file, case_ids, all_codes_upper)

    rows = []

    for gene in genes:
        for category in all_categories:
            status = build_carrier_status(gene, category, out_dir)
            if not status:
                continue

            carrier_cases     = {pid for pid, c in status.items() if c == 1}
            non_carrier_cases = case_ids - carrier_cases

            print(f"\n  {gene} | {category}: "
                  f"{len(carrier_cases)} carriers, "
                  f"{len(non_carrier_cases)} non-carriers")

            for disease, codes in disease_codes.items():
                codes_upper = {c.upper() for c in codes}

                carrier_dx        = sum(1 for p in carrier_cases
                                        if person_codes.get(p, set()) & codes_upper)
                carrier_no_dx     = len(carrier_cases) - carrier_dx
                non_carrier_dx    = sum(1 for p in non_carrier_cases
                                        if person_codes.get(p, set()) & codes_upper)
                non_carrier_no_dx = len(non_carrier_cases) - non_carrier_dx

                table = [[carrier_dx,     carrier_no_dx],
                         [non_carrier_dx, non_carrier_no_dx]]
                or_val, pval = fisher_exact(table, alternative="two-sided")

                rows.append({
                    "Gene":               gene,
                    "Category":           category,
                    "Disease":            disease,
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
                    "OddsRatio":          round(or_val, 4),
                    "Fisher_pvalue":      round(pval, 4),
                })

    results = pd.DataFrame(rows)
    out_path = out_dir / "carrier_liver_disease_enrichment.csv"
    results.to_csv(out_path, index=False)

    print("\n" + "=" * 80)
    print(results.to_string(index=False))
    print(f"\n>>> Summary file saved: {out_path}")
    print(">>> Please send this CSV file back for meta-analysis.")
    return results


# ---------------------------------------------------------------------------
# ARGUMENT PARSING
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Carrier comorbidity pipeline: tests enrichment of liver disease "
            "diagnoses in BRCA2/MSH6/PMS2/CHEK2/BRIP1/FANCA rare variant carriers "
            "vs non-carriers within a case population. Outputs summary-level CSV only."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  python carrier_comorbidity_pipeline.py \\\n"
            "    --brca2-group-file /path/to/cancer_genes.13.txt \\\n"
            "    --msh6-group-file  /path/to/cancer_genes.2.txt \\\n"
            "    --pms2-group-file  /path/to/cancer_genes.7.txt \\\n"
            "    --chek2-group-file /path/to/cancer_genes.22.txt \\\n"
            "    --brip1-group-file /path/to/cancer_genes.17.txt \\\n"
            "    --fanca-group-file /path/to/cancer_genes.16.txt \\\n"
            "    --bfile            /path/to/plink_bfile_prefix \\\n"
            "    --covariate-file   /path/to/covariates.csv \\\n"
            "    --case-column      HCC_Cancer_Free \\\n"
            "    --condition-file   /path/to/condition_occurrence.txt \\\n"
            "    --out-dir          /path/to/output_directory\n"
        ),
    )

    parser.add_argument(
        "--brca2-group-file", required=True, type=Path,
        metavar="PATH",
        help="SAIGE-GENE+ group file for the chromosome containing BRCA2 (chr13).",
    )
    parser.add_argument(
        "--msh6-group-file", required=True, type=Path,
        metavar="PATH",
        help="SAIGE-GENE+ group file for the chromosome containing MSH6 (chr2).",
    )
    parser.add_argument(
        "--pms2-group-file", required=True, type=Path,
        metavar="PATH",
        help="SAIGE-GENE+ group file for the chromosome containing PMS2 (chr7).",
    )
    parser.add_argument(
        "--chek2-group-file", required=True, type=Path,
        metavar="PATH",
        help="SAIGE-GENE+ group file for the chromosome containing CHEK2 (chr22).",
    )
    parser.add_argument(
        "--brip1-group-file", required=True, type=Path,
        metavar="PATH",
        help="SAIGE-GENE+ group file for the chromosome containing BRIP1 (chr17).",
    )
    parser.add_argument(
        "--fanca-group-file", required=True, type=Path,
        metavar="PATH",
        help="SAIGE-GENE+ group file for the chromosome containing FANCA (chr16).",
    )
    parser.add_argument(
        "--bfile", required=True, type=Path,
        metavar="PREFIX",
        help="PLINK bfile prefix (no .bed/.bim/.fam extension).",
    )
    parser.add_argument(
        "--covariate-file", required=True, type=Path,
        metavar="PATH",
        help="CSV file with at least two columns: person_id and <case-column>.",
    )
    parser.add_argument(
        "--case-column", required=True,
        metavar="COLUMN_NAME",
        help="Column name in the covariate file that flags cases (1) vs controls (0).",
    )
    parser.add_argument(
        "--condition-file", required=True, type=Path,
        metavar="PATH",
        help=(
            "Tab-separated file with columns [person_id, condition_source_value] "
            "containing ICD-9/ICD-10 codes. "
            "Adapt load_condition_codes() in the script if your format differs."
        ),
    )
    parser.add_argument(
        "--out-dir", required=True, type=Path,
        metavar="PATH",
        help="Directory where all output files will be written (created if absent).",
    )
    parser.add_argument(
        "--plink2-bin", default="plink2",
        metavar="PATH",
        help="Path to the plink2 binary. Default: 'plink2' (assumes it is on PATH).",
    )
    parser.add_argument(
        "--plof-label", default="pLoF",
        metavar="LABEL",
        help="Annotation label used for pLoF variants in group files (default: pLoF).",
    )
    parser.add_argument(
        "--damaging-label", default="damaging_missense",
        metavar="LABEL",
        help="Annotation label used for damaging missense variants in group files (default: damaging_missense).",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    group_files = {
        "BRCA2": args.brca2_group_file,
        "MSH6":  args.msh6_group_file,
        "PMS2":  args.pms2_group_file,
        "CHEK2": args.chek2_group_file,
        "BRIP1": args.brip1_group_file,
        "FANCA": args.fanca_group_file,
    }

    categories = [args.plof_label, args.damaging_label]

    # Step 1
    variant_files, case_ids = step1_parse_group_files(
        covariate_file    = args.covariate_file,
        case_column       = args.case_column,
        group_files       = group_files,
        genes_of_interest = GENES_OF_INTEREST,
        categories        = categories,
        out_dir           = args.out_dir,
        plof_label        = args.plof_label,
        damaging_label    = args.damaging_label,
    )

    # Step 2
    step2_plink_extract(
        plink_bfile = args.bfile,
        variant_files = variant_files,
        out_dir     = args.out_dir,
        plink2_bin  = args.plink2_bin,
    )

    # Step 3
    step3_enrichment(
        genes          = GENES_OF_INTEREST,
        categories     = categories,
        case_ids       = case_ids,
        condition_file = args.condition_file,
        disease_codes  = DISEASE_CODES,
        out_dir        = args.out_dir,
    )
