"""
age_at_onset_pipeline.py

Pipeline to compare age at HCC diagnosis between rare variant carriers
and non-carriers for BRCA2 and MSH6

Steps
-----
  1. Compute age at HCC diagnosis per case (birth date vs earliest HCC ICD date)
  2. Merge with carrier status files and run Welch's two-sample t-test
     (carrier vs non-carrier) for each gene x annotation category
  3. Write summary-level output CSV 

OUTPUT
------
  {cohort}_age_at_onset_results.csv
      One row per gene x category.
      Contains only n, mean, SD, median, IQR for carriers and non-carriers,
      plus Welch t-statistic and p-value 

USAGE
-----
  python age_at_onset_pipeline.py \\
    --cases-file     /path/to/cases_only.txt \\
    --person-file    /path/to/phenotype_person.txt \\
    --condition-file /path/to/condition_occurrence.txt \\
    --carrier-dir    /path/to/carrier_status_files/ \\
    --out-dir        /path/to/output_directory \\
    --cohort         PMBB

  Run with --help to see all options.

NOTES
-----
  - --cases-file: two-column tab-separated file with no header where the person ID is 
    repeated in both columns (this is the cases_only.txt produced by carrier_comorbidity_pipeline.py)
  - --person-file: tab-separated, must have columns: person_id, birth_datetime
    Date format: YYYY-MM-DD
  - --condition-file: tab-separated, must have columns: person_id,
    condition_source_value, condition_start_date
    Date format: YYYY-MM-DD
  - --carrier-dir: directory containing {GENE}_{CATEGORY}_carrier_status.tsv files
    (produced by carrier_comorbidity_pipeline.py)
  - Statistical test: Welch's t-test (two-sided, unequal variances)
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from scipy.stats import ttest_ind


# ---------------------------------------------------------------------------
# ARGUMENT PARSING
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Age at HCC onset pipeline: computes age at diagnosis then compares "
            "carriers vs non-carriers per gene x category using Welch's t-test. "
            "Outputs summary-level CSV only."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  python age_at_onset_pipeline.py \\\n"
            "    --cases-file     /path/to/cases_only.txt \\\n"
            "    --person-file    /path/to/phenotype_person.txt \\\n"
            "    --condition-file /path/to/condition_occurrence.txt \\\n"
            "    --carrier-dir    /path/to/carrier_status_files/ \\\n"
            "    --out-dir        /path/to/output_directory \\\n"
            "    --cohort         PMBB\n"
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
        "--person-file", required=True, type=Path,
        metavar="PATH",
        help="Tab-separated file with columns: person_id, birth_datetime (YYYY-MM-DD).",
    )
    parser.add_argument(
        "--condition-file", required=True, type=Path,
        metavar="PATH",
        help=(
            "Tab-separated file with columns: person_id, condition_source_value, "
            "condition_start_date (YYYY-MM-DD). Used to find earliest HCC diagnosis date."
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
        "--out-dir", required=True, type=Path,
        metavar="PATH",
        help="Directory where output files will be written (created if absent).",
    )
    parser.add_argument(
        "--cohort", default=None,
        metavar="LABEL",
        help="Cohort label used in output filenames and plot title (e.g. PMBB).",
    )
    parser.add_argument(
        "--genes", default="BRCA2,MSH6",
        metavar="GENE1,GENE2,...",
        help="Comma-separated list of genes to analyse. Default: BRCA2,MSH6",
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
    parser.add_argument(
        "--hcc-codes",
        default="155.0,155.2,C22.0,C22.2,C22.7,C22.8,C22.9",
        metavar="CODE1,CODE2,...",
        help=(
            "Comma-separated HCC ICD codes used to find diagnosis date. "
            "Default: 155.0,155.2,C22.0,C22.2,C22.7,C22.8,C22.9"
        ),
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# STEP 1: Compute age at HCC diagnosis
# ---------------------------------------------------------------------------

def step1_compute_age(cases_file, person_file, condition_file,
                      hcc_codes_upper, out_dir, prefix):
    print("\n=== STEP 1: Computing age at HCC diagnosis ===")

    # Load cases
    hcc_cases = set(
        pd.read_csv(cases_file, sep="\t", header=None, dtype=str)[0]
    )
    print(f"  {len(hcc_cases):,} cases loaded")

    # Load birth dates
    person = pd.read_csv(
        person_file, sep="\t",
        usecols=["person_id", "birth_datetime"],
        dtype={"person_id": str},
    )
    person["birth_datetime"] = pd.to_datetime(
        person["birth_datetime"], format="%Y-%m-%d", errors="coerce"
    )
    person = person[person["person_id"].isin(hcc_cases)]
    print(f"  {len(person):,} cases found in person file")

    # Find earliest HCC diagnosis date per person
    print(f"  Scanning condition file for HCC codes: {', '.join(sorted(hcc_codes_upper))}...")
    hcc_dates = {}
    for chunk in pd.read_csv(
        condition_file, sep="\t",
        usecols=["person_id", "condition_source_value", "condition_start_date"],
        dtype={"person_id": str, "condition_source_value": str},
        chunksize=500_000,
    ):
        chunk["condition_source_value"] = (
            chunk["condition_source_value"].str.strip().str.upper()
        )
        chunk = chunk[
            chunk["person_id"].isin(hcc_cases) &
            chunk["condition_source_value"].isin(hcc_codes_upper)
        ]
        if chunk.empty:
            continue
        chunk["condition_start_date"] = pd.to_datetime(
            chunk["condition_start_date"], format="%Y-%m-%d", errors="coerce"
        )
        chunk = chunk.dropna(subset=["condition_start_date"])
        for pid, date in zip(chunk["person_id"], chunk["condition_start_date"]):
            if pid not in hcc_dates or date < hcc_dates[pid]:
                hcc_dates[pid] = date

    print(f"  HCC date found for {len(hcc_dates):,} people")

    # Compute age
    hcc_date_df = pd.DataFrame(
        list(hcc_dates.items()), columns=["person_id", "hcc_date"]
    )
    merged = person.merge(hcc_date_df, on="person_id", how="inner")
    merged["age_at_hcc_dx"] = (
        (merged["hcc_date"] - merged["birth_datetime"]).dt.days / 365.25
    ).round(2)

    n_missing = len(hcc_cases) - len(merged)
    if n_missing > 0:
        print(f"  WARNING: {n_missing} cases had no date or birth info — excluded")

    print(f"  Age computed for {len(merged):,} people")
    print(f"  Median age: {merged['age_at_hcc_dx'].median():.1f} years  "
          f"Range: {merged['age_at_hcc_dx'].min():.1f}–{merged['age_at_hcc_dx'].max():.1f}")

    # Save intermediate age file (not for sharing — individual-level)
    age_path = out_dir / f"{prefix}age_at_hcc_diagnosis.tsv"
    merged[["person_id", "age_at_hcc_dx"]].to_csv(age_path, sep="\t", index=False)
    print(f"  Intermediate file saved: {age_path.name}")

    return merged[["person_id", "age_at_hcc_dx"]].dropna()


# ---------------------------------------------------------------------------
# STEP 2: Carrier vs non-carrier age comparison (Welch's t-test)
# ---------------------------------------------------------------------------

CATEGORY_LABELS = {
    "pLoF":                   "pLoF",
    "damaging_missense":      "Damaging Missense",
    "pLoF_damaging_missense": "pLoF + Damaging Missense",
}

CARRIER_COLOR    = "#E53935"
NONCARRIER_COLOR = "#90A4AE"


def summary_stats(arr):
    if len(arr) == 0:
        return dict(n=0, mean=None, sd=None, median=None, q1=None, q3=None)
    return dict(
        n      = len(arr),
        mean   = round(float(np.mean(arr)), 1),
        sd     = round(float(np.std(arr, ddof=1)), 1) if len(arr) >= 2 else None,
        median = round(float(np.median(arr)), 1),
        q1     = round(float(np.percentile(arr, 25)), 1),
        q3     = round(float(np.percentile(arr, 75)), 1),
    )


def step2_age_comparison(age_df, carrier_dir, genes, categories,
                          cohort_label, prefix, out_dir):
    print("\n=== STEP 2: Carrier vs non-carrier age comparison ===")

    n_genes = len(genes)
    n_cats  = len(categories)

    fig, axes = plt.subplots(
        n_genes, n_cats,
        figsize=(14, n_genes * 2.8),
        sharey=False,
        squeeze=False,
    )

    rows = []

    for gi, gene in enumerate(genes):
        for ci, category in enumerate(categories):
            ax = axes[gi][ci]

            carrier_path = carrier_dir / f"{gene}_{category}_carrier_status.tsv"
            if not carrier_path.exists():
                print(f"  WARNING: {carrier_path.name} not found — skipping")
                ax.set_visible(False)
                continue

            carrier_df = pd.read_csv(
                carrier_path, sep="\t", dtype={"person_id": str}
            )
            merged = age_df.merge(carrier_df, on="person_id", how="inner")

            carrier_ages    = merged.loc[merged["carrier"] == 1, "age_at_hcc_dx"].values
            noncarrier_ages = merged.loc[merged["carrier"] == 0, "age_at_hcc_dx"].values

            # Welch's t-test (two-sided, unequal variances)
            if len(carrier_ages) >= 2 and len(noncarrier_ages) >= 2:
                t_stat, pval = ttest_ind(carrier_ages, noncarrier_ages, equal_var=False)
            else:
                t_stat, pval = None, None

            cs = summary_stats(carrier_ages)
            ns = summary_stats(noncarrier_ages)

            print(f"  {gene} | {category}: "
                  f"{cs['n']} carriers (median {cs['median']} yrs), "
                  f"{ns['n']} non-carriers (median {ns['median']} yrs), "
                  f"p={round(pval, 4) if pval is not None else 'NA'}")

            rows.append({
                "Gene":                  gene,
                "Category":              category,
                "Carrier_n":             cs["n"],
                "Carrier_mean_age":      cs["mean"],
                "Carrier_sd_age":        cs["sd"],
                "Carrier_median_age":    cs["median"],
                "Carrier_IQR":           f"{cs['q1']}–{cs['q3']}" if cs["q1"] is not None else "NA",
                "NonCarrier_n":          ns["n"],
                "NonCarrier_mean_age":   ns["mean"],
                "NonCarrier_sd_age":     ns["sd"],
                "NonCarrier_median_age": ns["median"],
                "NonCarrier_IQR":        f"{ns['q1']}–{ns['q3']}" if ns["q1"] is not None else "NA",
                "Welch_t":               round(t_stat, 3) if t_stat is not None else "NA",
                "pvalue":                round(pval,   4) if pval   is not None else "NA",
            })

            # Box plot
            data_to_plot, labels, colors = [], [], []
            if len(carrier_ages) > 0:
                data_to_plot.append(carrier_ages)
                labels.append(f"Carrier\n(n={cs['n']})")
                colors.append(CARRIER_COLOR)
            if len(noncarrier_ages) > 0:
                data_to_plot.append(noncarrier_ages)
                labels.append(f"Non-Carrier\n(n={ns['n']})")
                colors.append(NONCARRIER_COLOR)

            bp = ax.boxplot(
                data_to_plot,
                patch_artist=True,
                widths=0.5,
                medianprops=dict(color="black", linewidth=1.5),
            )
            for patch, color in zip(bp["boxes"], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.8)

            ax.set_xticklabels(labels, fontsize=7)
            ax.set_ylabel("Age at HCC diagnosis (years)", fontsize=7)
            ax.tick_params(axis="y", labelsize=7)

            # Significance bracket
            if pval is not None and pval < 0.05 and len(data_to_plot) == 2:
                y_max = max(
                    np.max(data_to_plot[0]) if len(data_to_plot[0]) > 0 else 0,
                    np.max(data_to_plot[1]) if len(data_to_plot[1]) > 0 else 0,
                )
                y_line = y_max * 1.05
                ax.plot([1, 1, 2, 2],
                        [y_line, y_line + y_max * 0.02,
                         y_line + y_max * 0.02, y_line],
                        color="black", linewidth=1)
                ax.text(1.5, y_max * 1.08, "*",
                        ha="center", va="bottom", fontsize=12,
                        fontweight="bold", color="black")
                ax.set_ylim(top=y_max * 1.22)

            title = f"{gene}\n{CATEGORY_LABELS.get(category, category)}"
            if pval is not None:
                title += f"\np={pval:.4f}"
            ax.set_title(title, fontsize=8, fontweight="bold")

    # Legend and layout
    legend_handles = [
        mpatches.Patch(color=CARRIER_COLOR,    label="Carrier"),
        mpatches.Patch(color=NONCARRIER_COLOR, label="Non-Carrier"),
        plt.Line2D([0], [0], marker="*", color="black", linestyle="None",
                   markersize=10, label="p < 0.05"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=3,
               fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(
        f"Age at HCC Diagnosis: Carriers vs Non-Carriers ({cohort_label})\n"
        "Welch's Two-Sample t-Test",
        fontsize=12, fontweight="bold", y=1.01,
    )
    plt.tight_layout()

    plot_path = out_dir / f"{prefix}age_at_onset_boxplot.png"
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    print(f"\n  Box plot saved: {plot_path.name}")

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    genes      = [g.strip() for g in args.genes.split(",")]
    categories = [c.strip() for c in args.categories.split(",")]
    hcc_codes_upper = {c.strip().upper() for c in args.hcc_codes.split(",")}

    cohort_label = args.cohort or "cohort"
    prefix       = f"{args.cohort}_" if args.cohort else ""

    # Step 1
    age_df = step1_compute_age(
        cases_file     = args.cases_file,
        person_file    = args.person_file,
        condition_file = args.condition_file,
        hcc_codes_upper = hcc_codes_upper,
        out_dir        = args.out_dir,
        prefix         = prefix,
    )

    # Step 2
    results = step2_age_comparison(
        age_df       = age_df,
        carrier_dir  = args.carrier_dir,
        genes        = genes,
        categories   = categories,
        cohort_label = cohort_label,
        prefix       = prefix,
        out_dir      = args.out_dir,
    )

    # Save summary CSV
    csv_path = args.out_dir / f"{prefix}age_at_onset_results.csv"
    results.to_csv(csv_path, index=False)

    print("\n" + "=" * 80)
    print(results.to_string(index=False))
    print(f"\n>>> Summary file saved: {csv_path}")
    print(">>> Please send this CSV file back for meta-analysis.")
