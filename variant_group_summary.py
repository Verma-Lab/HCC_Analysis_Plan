
#Run Command
# python variant_group_summary.py --input filename --gene gene --chr CHR --output_dir outdir

import pandas as pd
import os
import argparse

def parse_arguments():
    """Set up command line argument parsing."""
    parser = argparse.ArgumentParser(description='Process variant file and create group file.')
    
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
        help='Gene name (e.g., msh6, brca2)'
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
    """Calculate detailed summaries for both variant types."""
    summaries = {
        "pLoF": {},
        "damaging_missense": {}
    }
    
    if not plof_df.empty:
        summaries["pLoF"][gene] = { 
            "clinvar_hcc_inc": len(plof_df[plof_df['clinvar_hcc_inc'] == 1]['Uploaded_variation'].unique()),
            "LoF_HC": len(plof_df[plof_df['LoF'] == "HC"]),
            "Lof_HC_MANE": len(plof_df[(plof_df['LoF'] == "HC") & (plof_df['MANE_SELECT'] != '-')]),
            "Splice_variants": len(plof_df[
                (plof_df['Consequence'].isin(["splice_donor_variant", "splice_acceptor_variant", "splice_region_variant"])) &
                (plof_df['SpliceAI_DS'] >= 0.2)
            ]),
            "Splice_variants_MANE": len(plof_df[
                (plof_df['Consequence'].isin(["splice_donor_variant", "splice_acceptor_variant", "splice_region_variant"])) &
                (plof_df['SpliceAI_DS'] >= 0.2) & 
                (plof_df['MANE_SELECT'] != '-')
            ]),
            "non_excluded": len(plof_df[plof_df['clinvar_hcc_excl'] != 1])
        }
        
    else:
        print(f"\nWARNING: pLoF DataFrame is empty for {gene}. Falling back to full DataFrame.")
        summaries["pLoF"][gene] = {
            "clinvar_hcc_inc": len(df[df['clinvar_hcc_inc'] == 1]['Uploaded_variation'].unique()),
            "LoF_HC": len(df[df['LoF'] == "HC"]),
            "Lof_HC_MANE": len(df[(df['LoF'] == "HC") & (df['MANE_SELECT'] != '-')]),
            "Splice_variants": len(df[
                (df['Consequence'].isin(["splice_donor_variant", "splice_acceptor_variant", "splice_region_variant"])) &
                (df['SpliceAI_DS'] >= 0.2)
            ]),
            "Splice_variants_MANE": len(df[
                (df['Consequence'].isin(["splice_donor_variant", "splice_acceptor_variant", "splice_region_variant"])) &
                (df['SpliceAI_DS'] >= 0.2) & 
                (df['MANE_SELECT'] != '-')
            ]),
            "non_excluded": len(df[df['clinvar_hcc_excl'] != 1])
        }

    # Calculate damaging missense summaries
    if not damaging_df.empty:
        summaries["damaging_missense"][gene] = {
            "Not_plof": len(damaging_df),
            "Consequence_missense_start_stop_indel_splicing": len(damaging_df[
                damaging_df['Consequence'].isin([
                    "missense_variant", "start_lost", "stop_lost",
                    "inframe_insertion", "inframe_deletion",
                    "splice_region_variant", "splice_donor_variant", "splice_acceptor_variant"
                ])
            ]),
            "Revel_score": len(damaging_df[damaging_df['REVEL_score'] >= 0.773]),
            "CADD_PHRED_merge": len(damaging_df[damaging_df['CADD_PHRED_merge'] >= 28.1]),
            "Splice_variants": len(damaging_df[damaging_df['SpliceAI_DS'] >= 0.2]),
            "Splice_variants_MANE": len(damaging_df[
                (damaging_df['SpliceAI_DS'] >= 0.2) & (damaging_df['MANE_SELECT'] != '-')
            ]),
            "LoF_LC": len(damaging_df[damaging_df['LoF'] == "LC"])
        }

    else:
        print(f"\nWARNING: Damaging Missense DataFrame is empty for {gene}. Falling back to full DataFrame.")
        summaries["damaging_missense"][gene] = {
            "Not_plof": len(df[~df['Uploaded_variation'].isin(plof_df['Uploaded_variation'])]),
            "Consequence_missense_start_stop_indel_splicing": len(df[
                ~df['Uploaded_variation'].isin(plof_df['Uploaded_variation']) &
                df['Consequence'].isin([
                    "missense_variant", "start_lost", "stop_lost",
                    "inframe_insertion", "inframe_deletion",
                    "splice_region_variant", "splice_donor_variant", "splice_acceptor_variant"
                ])
            ]),
            "Revel_score": len(df[df['REVEL_score'] >= 0.773]),
            "CADD_PHRED_merge": len(df[df['CADD_PHRED_merge'] >= 28.1]),
            "Splice_variants": len(df[df['SpliceAI_DS'] >= 0.2]),
            "Splice_variants_MANE": len(df[
                (df['SpliceAI_DS'] >= 0.2) & (df['MANE_SELECT'] != '-')
            ]),
            "LoF_LC": len(df[df['LoF'] == "LC"])
        }
    # # pLoF summaries
    # summaries["pLoF"][gene] = {
    #     "clinvar_hcc_inc": len(plof_df[plof_df['clinvar_hcc_inc'] == 1]['Uploaded_variation'].unique()),
    #     "LoF_HC": len(plof_df[plof_df['LoF'] == "HC"]),
    #     "Lof_HC_MANE": len(plof_df[(plof_df['LoF'] == "HC") & (plof_df['MANE_SELECT'] != '-')]),
    #     "Splice_variants": len(plof_df[
    #         (plof_df['Consequence'].isin(["splice_donor_variant", "splice_acceptor_variant", "splice_region_variant"])) &
    #         (plof_df['SpliceAI_DS'] >= 0.2)
    #     ]),
    #     "Splice_variants_MANE": len(plof_df[
    #         (plof_df['Consequence'].isin(["splice_donor_variant", "splice_acceptor_variant", "splice_region_variant"])) &
    #         (plof_df['SpliceAI_DS'] >= 0.2) & 
    #         (plof_df['MANE_SELECT'] != '-')
    #     ]),
    #     "non_excluded": len(plof_df[plof_df['clinvar_hcc_excl'] != 1])
    # }
    
    # Damaging missense summaries
    # summaries["damaging_missense"][gene] = {
    #     "Not_plof": len(damaging_df),
    #     "Consequence_missense_start_stop_indel_splicing": len(damaging_df[
    #         damaging_df['Consequence'].isin([
    #             "missense_variant", "start_lost", "stop_lost",
    #             "inframe_insertion", "inframe_deletion",
    #             "splice_region_variant", "splice_donor_variant", "splice_acceptor_variant"
    #         ])
    #     ]),
    #     "Revel_score": len(damaging_df[damaging_df['REVEL_score'] >= 0.773]),
    #     "CADD_PHRED_merge": len(damaging_df[damaging_df['CADD_PHRED_merge'] >= 28.1]),
    #     "Splice_variants": len(damaging_df[damaging_df['SpliceAI_DS'] >= 0.2]),
    #     "Splice_variants_MANE": len(damaging_df[
    #         (damaging_df['SpliceAI_DS'] >= 0.2) & (damaging_df['MANE_SELECT'] != '-')
    #     ]),
    #     "LoF_LC": len(damaging_df[damaging_df['LoF'] == "LC"])
    # }
    
    return summaries

def format_summary_as_dataframe(summaries, gene):
    """Convert summary dictionary to DataFrame format."""
    rows = []
    
    # Process pLoF summary
    plof_stats = summaries["pLoF"][gene]
    for metric, value in plof_stats.items():
        rows.append({
            'gene': gene,
            'category': 'pLoF',
            'metric': metric,
            'value': value
        })
    
    # Process damaging_missense summary
    damaging_stats = summaries["damaging_missense"][gene]
    for metric, value in damaging_stats.items():
        rows.append({
            'gene': gene,
            'category': 'damaging_missense',
            'metric': metric,
            'value': value
        })
    
    # Create DataFrame and sort by category and metric
    df = pd.DataFrame(rows)
    df = df.sort_values(['category', 'metric'])
    return df

def process_file(input_file, gene, chromosome, output_dir):
    """Process variant file and create group file."""
    try:
        # Read input file
        print(f"Reading file: {input_file}")
        df = pd.read_csv(input_file, sep='\t')
        
        #Remove duplicates in Uploaded_variation column
        print("Before deduplication:", df['Uploaded_variation'].nunique())

        df = df.drop_duplicates(subset=['Uploaded_variation'])
        
        print("After deduplication:", df['Uploaded_variation'].nunique())

        # Process SpliceAI scores
        df = process_spliceai(df)
        
        # Convert numeric columns
        df['CADD_PHRED_merge'] = pd.to_numeric(df['CADD_PHRED_merge'], errors='coerce')
        df['REVEL_score'] = pd.to_numeric(df['REVEL_score'], errors='coerce')
        
        # Get unique ClinVar variants
        clinvar_variants = df[df['clinvar_hcc_inc'] == 1]['Uploaded_variation'].unique()
        clinvar_mask = df['Uploaded_variation'].isin(clinvar_variants)
        
        # Identify pLoF variants
        plof_mask = (
            (clinvar_mask |
             ((df['LoF'] == "HC") & (df['MANE_SELECT'] != '-')) |
             ((df['Consequence'].isin(["splice_donor_variant", "splice_acceptor_variant", "splice_region_variant"])) &
              (df['SpliceAI_DS'] >= 0.2))) &
            (df['MANE_SELECT'] != '-') &
            ~(df['clinvar_hcc_excl'] == 1)
        )
        
        plof_df = df[plof_mask].copy()
        plof_df['anno'] = 'pLoF'
        
        # Identify damaging missense variants
        damaging_mask = (
            ~plof_mask &
            (df['Consequence'].isin([
                "missense_variant", "start_lost", "stop_lost",
                "inframe_insertion", "inframe_deletion",
                "splice_region_variant", "splice_donor_variant", "splice_acceptor_variant"
            ])) &
            ((df['REVEL_score'] >= 0.773) |
             (df['CADD_PHRED_merge'] >= 28.1) |
             ((df['SpliceAI_DS'] >= 0.2) & (df['MANE_SELECT'] != '-')) |
             (df['LoF'] == "LC"))
        )
        
        damaging_df = df[damaging_mask].copy()
        damaging_df['anno'] = 'damaging_missense'
        
        # Calculate detailed summaries
        summaries = calculate_summaries(plof_df, damaging_df, gene, df)
        
        # Convert summaries to DataFrame format
        summary_df = format_summary_as_dataframe(summaries, gene)
        
        # Save summary DataFrame
        summary_file_path = os.path.join(output_dir, f"summary_results_{gene}.txt")
        summary_df.to_csv(summary_file_path, sep='\t', index=False)
        print(f"\nSaved summary to: {summary_file_path}")
        
        # Create group file
        group_file_path = os.path.join(output_dir, f"cancer_genes.{chromosome}.txt")
        
        # Combine variants
        combined_df = pd.concat([plof_df, damaging_df], ignore_index=True)
        
        with open(group_file_path, 'w') as f:
            gene_upper = gene.upper()
            var_line = f"{gene_upper} var {' '.join(combined_df['Uploaded_variation'].astype(str))}"
            anno_line = f"{gene_upper} anno {' '.join(combined_df['anno'].astype(str))}"
            f.write(f"{var_line}\n{anno_line}\n")
        
        print(f"Created group file: {group_file_path}")
        
        
        # Print summary to console
        print("\nSummary Results:")
        print(summary_df.to_string())
        
        return True
        
    except Exception as e:
        print(f"Error processing file: {str(e)}")
        return False

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Process the file
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


#Run in Bash 
# for i in *merging_on_chr_pos_ref_alt.txt; 
# do 
#     echo $i; 
#     gene=$(basename $i | cut -d "_" -f3); 
#     chr=$(basename $i | cut -d"_" -f4 | sed 's/chr//g'); 
#     python variant_group_summary.py 
#     --input ${i} 
#     --gene ${gene}  
#     --chr ${chr} 
#     --output_dir /vep_group_out/ ;
# done