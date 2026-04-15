# Sensitivity analysis version of variant_group_summary.py
#
# PURPOSE:
#   Regenerates burden test group files for MSH6 and PMS2 restricting both pLoF and
#   damaging missense categories to ClinVar pathogenic/likely pathogenic (P/LP) variants
#   only. Variants included solely based on computational predictions (LOFTEE, REVEL, CADD, SpliceAI) were excluded
#   Consequence type determines category assignment:
#     - pLoF: ClinVar P/LP + stop_gained, frameshift_variant, splice_donor/acceptor/region, etc.
#     - Damaging missense: ClinVar P/LP + missense_variant, inframe indel, etc. (not already in pLoF)
#   MAF <= 1% in gnomAD and MANE transcript filters are unchanged from the original script.
#
# INPUT FILE:
#   The per-gene VEP annotation file with ClinVar inclusion/exclusion flag columns
#   (clinvar_hcc_inc, clinvar_hcc_excl) already added. This is the same input file
#   used by variant_group_summary.py. The filename does not matter, any tab-separated
#   file with the required columns will work. Required columns: Uploaded_variation,
#   Consequence, MANE_SELECT, CANONICAL, LoF, SpliceAI_pred, REVEL_score, CADD_PHRED,
#   gnomADe_AF, SYMBOL, clinvar_hcc_inc, clinvar_hcc_excl.
#
# ARGUMENTS:
#   --input_file   Path to the merged VEP + ClinVar annotation file for a single gene
#   --gene         Gene name in lowercase (e.g., msh6, pms2)
#   --chr          Chromosome number only, no 'chr' prefix (e.g., 2, 7)
#   --output_dir   Directory where output group file and summary will be written
#
# OUTPUTS:
#   - cancer_genes_sensitivity.{chr}.txt        (group file for Nextflow)
#   - summary_results_sensitivity_{gene}.txt    (variant count summary)
#
# HOW TO RUN:
#   This script must be run twice — once for MSH6 and once for PMS2.
#   Replace /path/to/ with your actual file paths.
#
#   Run 1 — MSH6:
#   python variant_group_summary_sensitivity.py \
#       --input_file /path/to/vep_annot_msh6_chr2_with_curated_list_cols_merging_on_chr_pos_ref_alt.txt \
#       --gene msh6 \
#       --chr 2 \
#       --output_dir /path/to/sensitivity_group_files/
#
#   Run 2 — PMS2:
#   python variant_group_summary_sensitivity.py \
#       --input_file /path/to/vep_annot_pms2_chr7_with_curated_list_cols_merging_on_chr_pos_ref_alt.txt \
#       --gene pms2 \
#       --chr 7 \
#       --output_dir /path/to/sensitivity_group_files/

import pandas as pd
import os
import argparse

def parse_arguments():
    """Set up command line argument parsing."""
    parser = argparse.ArgumentParser(
        description='Process variant file and create ClinVar-restricted sensitivity group file.'
    )

    parser.add_argument(
        '--input_file',
        type=str,
        required=True,
        help='Path to input variant file'
    )

    parser.add_argument(
        '--gene',
        type=str,
        required=True,
        help='Gene name (e.g., msh6, pms2)'
    )

    parser.add_argument(
        '--chr',
        type=str,
        required=True,
        help='Chromosome number'
    )

    parser.add_argument(
        '--output_dir',
        type=str,
        required=True,
        help='Directory where output files will be saved'
    )

    return parser.parse_args()

def process_spliceai(df):
    """Process SpliceAI predictions and return DataFrame with max score."""
    if 'SpliceAI_pred' not in df.columns:
        return df

    splice = df['SpliceAI_pred'].str.split('|', expand=True)
    if splice.empty:
        return df

    splice.columns = ['SYMBOL', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL',
                      'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL']
    splice = splice.mask(splice == 'None')

    score_cols = ['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL',
                  'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL']
    splice[score_cols] = splice[score_cols].astype(float)
    splice['SpliceAI_DS'] = splice[['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL']].max(axis=1)

    return pd.concat([df, splice['SpliceAI_DS']], axis=1)

def calculate_summaries(plof_df, damaging_df, gene, df):
    """Calculate summaries for sensitivity analysis variant categories."""
    summaries = {
        "pLoF": {},
        "damaging_missense": {},
        "pLoF+damaging_missense": {}
    }

    if not plof_df.empty:
        summaries["pLoF"][gene] = {
            "clinvar_hcc_inc": len(plof_df[plof_df['clinvar_hcc_inc'] == 1]['Uploaded_variation'].unique()),
            "VEP_High_Impact": len(plof_df[plof_df['Consequence'].isin([
                "transcript_ablation", "stop_gained", "frameshift_variant",
                "stop_lost", "start_lost", "transcript_amplification",
                "feature_elongation", "feature_truncation"
            ])]),
            "Splice_variants": len(plof_df[plof_df['Consequence'].isin([
                "splice_donor_variant", "splice_acceptor_variant", "splice_region_variant"
            ])]),
            "non_excluded": len(plof_df[plof_df['clinvar_hcc_excl'] != 1]),
            "gnomad combined population AF <= 0.01": len(plof_df[plof_df['gnomADe_AF'] <= 0.01]),
            "gnomad combined population AF == '-'": len(plof_df[plof_df['gnomADe_AF_missing']])
        }
    else:
        print(f"\nWARNING: pLoF DataFrame is empty for {gene}.")
        summaries["pLoF"][gene] = {
            "clinvar_hcc_inc": 0,
            "VEP_High_Impact": 0,
            "Splice_variants": 0,
            "non_excluded": 0,
            "gnomad combined population AF <= 0.01": 0,
            "gnomad combined population AF == '-'": 0
        }

    if not damaging_df.empty:
        summaries["damaging_missense"][gene] = {
            "clinvar_hcc_inc": len(damaging_df[damaging_df['clinvar_hcc_inc'] == 1]['Uploaded_variation'].unique()),
            "Consequence_missense_start_stop_indel_splicing": len(damaging_df[
                damaging_df['Consequence'].isin([
                    "missense_variant", "start_lost", "stop_lost",
                    "inframe_insertion", "inframe_deletion",
                    "splice_region_variant", "splice_donor_variant", "splice_acceptor_variant"
                ])
            ]),
            "gnomad combined population AF <= 0.01": len(damaging_df[damaging_df['gnomADe_AF'] <= 0.01]),
            "gnomad combined population AF == '-'": len(damaging_df[damaging_df['gnomADe_AF_missing']])
        }
    else:
        print(f"\nWARNING: Damaging Missense DataFrame is empty for {gene}.")
        summaries["damaging_missense"][gene] = {
            "clinvar_hcc_inc": 0,
            "Consequence_missense_start_stop_indel_splicing": 0,
            "gnomad combined population AF <= 0.01": 0,
            "gnomad combined population AF == '-'": 0
        }

    if not plof_df.empty or not damaging_df.empty:
        all_variants = pd.concat([
            plof_df['Uploaded_variation'] if not plof_df.empty else pd.Series(),
            damaging_df['Uploaded_variation'] if not damaging_df.empty else pd.Series()
        ]).unique()

        summaries["pLoF+damaging_missense"][gene] = {
            "Total": len(all_variants)
        }

    return summaries

def format_summary_as_dataframe(summaries, gene):
    """Convert summary dictionary to DataFrame format."""
    rows = []

    for metric, value in summaries["pLoF"][gene].items():
        rows.append({'gene': gene, 'category': 'pLoF', 'metric': metric, 'value': value})

    for metric, value in summaries["damaging_missense"][gene].items():
        rows.append({'gene': gene, 'category': 'damaging_missense', 'metric': metric, 'value': value})

    if "pLoF+damaging_missense" in summaries and gene in summaries["pLoF+damaging_missense"]:
        for metric, value in summaries["pLoF+damaging_missense"][gene].items():
            rows.append({'gene': gene, 'category': 'pLoF+damaging_missense', 'metric': metric, 'value': value})

    df = pd.DataFrame(rows)
    df = df.sort_values(['category', 'metric'])
    return df

def process_file(input_file, gene, chromosome, output_dir):
    """Process variant file and create ClinVar-restricted sensitivity group file."""
    try:
        print(f"Reading file: {input_file}")
        df = pd.read_csv(input_file, sep='\t')

        print(f"\nFile loaded. Original row count: {len(df)}")
        print("MANE_SELECT value counts before any processing:")
        print(df['MANE_SELECT'].value_counts().to_string())
        print("\nCANONICAL value counts before any processing:")
        print(df['CANONICAL'].value_counts().to_string())

        # Keep only MANE transcript rows
        df = df[df['MANE_SELECT'] != '-'].copy()
        print(f"Rows with MANE_SELECT: {len(df)}")

        # Process SpliceAI scores (retained for summary reporting, not inclusion)
        df = process_spliceai(df)

        # Convert numeric columns
        df['CADD_PHRED'] = pd.to_numeric(df['CADD_PHRED'], errors='coerce')
        df['REVEL_score'] = pd.to_numeric(df['REVEL_score'], errors='coerce')

        # Flag missing gnomADe_AF values before converting to numeric
        df['gnomADe_AF_missing'] = df['gnomADe_AF'] == '-'
        df['gnomADe_AF'] = pd.to_numeric(df['gnomADe_AF'], errors='coerce')

        # Get unique ClinVar P/LP variants
        clinvar_variants = df[df['clinvar_hcc_inc'] == 1]['Uploaded_variation'].unique()
        clinvar_mask = df['Uploaded_variation'].isin(clinvar_variants)

        print(f"\nClinVar P/LP variants (all consequences, pre-dedup): {clinvar_mask.sum()}")

        # Remove duplicates
        dup_drop_cols = ['Uploaded_variation', 'MANE_SELECT', 'Consequence', 'LoF', 'gnomADe_AF',
                         'clinvar_hcc_inc', 'CADD_PHRED', 'clinvar_hcc_excl', 'SpliceAI_DS', 'REVEL_score']
        df = df.drop_duplicates(subset=dup_drop_cols)
        print(f"After duplicate removal: {len(df)} rows")

        # Recompute clinvar_mask after dedup
        clinvar_variants = df[df['clinvar_hcc_inc'] == 1]['Uploaded_variation'].unique()
        clinvar_mask = df['Uploaded_variation'].isin(clinvar_variants)

        # --- SENSITIVITY pLoF mask ---
        # ClinVar P/LP only, consequence determines category (no LOFTEE/SpliceAI/VEP OR conditions)
        plof_mask = (
            clinvar_mask &
            (df['Consequence'].isin([
                "transcript_ablation", "stop_gained", "frameshift_variant",
                "stop_lost", "start_lost", "transcript_amplification",
                "feature_elongation", "feature_truncation",
                "splice_donor_variant", "splice_acceptor_variant", "splice_region_variant"
            ])) &
            (df['MANE_SELECT'] != '-') &
            ~(df['clinvar_hcc_excl'] == 1) &
            ((df['gnomADe_AF'] <= 0.01) | (df['gnomADe_AF_missing']))
        )

        plof_df = df[plof_mask].copy()
        plof_df['anno'] = 'pLoF'

        # --- SENSITIVITY damaging missense mask ---
        # ClinVar P/LP only, missense-type consequence, not already in pLoF
        damaging_mask = (
            clinvar_mask &
            ~plof_mask &
            (df['Consequence'].isin([
                "missense_variant", "start_lost", "stop_lost",
                "inframe_insertion", "inframe_deletion",
                "splice_region_variant", "splice_donor_variant", "splice_acceptor_variant"
            ])) &
            (df['MANE_SELECT'] != '-') &
            ((df['gnomADe_AF'] <= 0.01) | (df['gnomADe_AF_missing']))
        )

        damaging_df = df[damaging_mask].copy()
        damaging_df['anno'] = 'damaging_missense'

        # Check for and remove duplicate pLoF variant-gene pairs
        plof_count_before = len(plof_df)
        plof_df = plof_df.set_index(['Uploaded_variation', 'SYMBOL'])
        check_dup_plof = plof_df[plof_df.index.duplicated(keep=False)]
        if len(check_dup_plof) > 0:
            print(f"Found {len(check_dup_plof)} duplicate pLoF variant-gene pairs")
            plof_df = plof_df[~plof_df.index.duplicated(keep='first')]
        plof_df = plof_df.reset_index()
        print(f"Removed {plof_count_before - len(plof_df)} duplicate pLoF rows")

        # Check for and remove duplicate damaging missense variant-gene pairs
        damaging_count_before = len(damaging_df)
        damaging_df = damaging_df.set_index(['Uploaded_variation', 'SYMBOL'])
        check_dup_damaging = damaging_df[damaging_df.index.duplicated(keep=False)]
        if len(check_dup_damaging) > 0:
            print(f"Found {len(check_dup_damaging)} duplicate damaging missense variant-gene pairs")
            damaging_df = damaging_df[~damaging_df.index.duplicated(keep='first')]
        damaging_df = damaging_df.reset_index()
        print(f"Removed {damaging_count_before - len(damaging_df)} duplicate damaging missense rows")

        # Create sensitivity group file
        group_file_path = os.path.join(output_dir, f"cancer_genes_sensitivity.{chromosome}.txt")

        combined_df = pd.concat([plof_df, damaging_df], ignore_index=True)

        with open(group_file_path, 'w') as f:
            gene_upper = gene.upper()
            var_line = f"{gene_upper} var {' '.join(combined_df['Uploaded_variation'].astype(str))}"
            anno_line = f"{gene_upper} anno {' '.join(combined_df['anno'].astype(str))}"
            f.write(f"{var_line}\n{anno_line}\n")
        print(f"Created sensitivity group file: {group_file_path}")

        # Calculate summaries
        summaries = calculate_summaries(plof_df, damaging_df, gene, df)
        summary_df = format_summary_as_dataframe(summaries, gene)

        # Count totals from group file
        with open(group_file_path, 'r') as f:
            anno_line = [line for line in f.readlines() if 'anno' in line][0]
            plof_total = anno_line.count("pLoF")
            damaging_total = anno_line.count("damaging_missense")

        total_rows = pd.DataFrame([
            {'gene': gene, 'category': 'pLoF', 'metric': 'Total', 'value': plof_total},
            {'gene': gene, 'category': 'damaging_missense', 'metric': 'Total', 'value': damaging_total}
        ])

        summary_df = pd.concat([summary_df, total_rows], ignore_index=True)

        summary_file_path = os.path.join(output_dir, f"summary_results_sensitivity_{gene}.txt")
        summary_df.to_csv(summary_file_path, sep='\t', index=False)
        print(f"\nSaved summary to: {summary_file_path}")

        print("\nSummary Results:")
        print(summary_df.to_string())

        return True

    except Exception as e:
        print(f"Error processing file: {str(e)}")
        return False

def main():
    args = parse_arguments()
    os.makedirs(args.output_dir, exist_ok=True)

    success = process_file(
        args.input_file,
        args.gene,
        args.chr,
        args.output_dir
    )

    if success:
        print("\nProcessing completed successfully!")
    else:
        print("\nProcessing failed!")

if __name__ == "__main__":
    main()