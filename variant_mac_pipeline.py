#!/usr/bin/env python3

"""
Variant Mac Pipeline 

Pipeline to calculate minor allele counts (MAC) for cases and controls from variants in a group file.

EXAMPLE USAGE:
python3 variant_mac_pipeline.py \
  --group-files /home/agaro/verma_shared/projects/HCC/PMBB_6Gene_Burden_Analysis_030625/group_files/cancer_genes.*.txt \
  --phenotype-file /project/verma_shared/projects/HCC/covariates/covariates_hcc_cancer_free_controls_exomepcs_ancestry_column_3.0.csv \
  --bfile /static/PMBB/PMBB-Release-2024-3.0/Exome/pVCF/NF_all_variants/PMBB-Release-2024-3.0_genetic_exome_NF \
  --output-dir /home/agaro/verma_shared/projects/HCC/PMBB_6Gene_Burden_Analysis_030625/ALL_ALL/Variant_Level_Analysis_Test \
  --case-column HCC_Cancer_Free \
  --id-column person_id

This pipeline:
1. Extracts variants from group files and combines them 
2. Creates PLINK keep files for cases (HCC_Cancer_Free=1) and controls (HCC_Cancer_Free=0)
3. Runs PLINK2 genotype counting for cases and controls
4. Calculates Minor Allele Counts (MAC) for both groups
5. Creates final summary: MAC_case_control_summary.txt

"""

import pandas as pd
import argparse
import sys
from pathlib import Path 
import subprocess
import os 
import glob 

class VariantMACPipeline:
    def __init__(self, args):
        self.args = args
        self.output_dir = Path(args.output_dir)
        self.output_dir.mkdir(parents = True, exist_ok= True)

    def extract_variants_from_file(self, file_path):
        """Extract variants from a file after the 'var' keyword."""
        try:
            with open(file_path, "r") as file:
                first_line = file.readline().strip()
                parts = first_line.split()
                
                # Find 'var' and extract everything after it
                if "var" in parts:
                    var_index = parts.index("var")
                    return parts[var_index + 1:]
                else:
                    print(f"Warning: 'var' keyword not found in {file_path}")
                    return []
        except FileNotFoundError:
            print(f"Error: File not found - {file_path}")
            return []
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            return []

    def find_group_files(self):
        """Find group files based on the provided pattern or directory."""
        group_files = []
        
        if self.args.group_files:
            # Specific files provided
            for file_pattern in self.args.group_files:
                if '*' in file_pattern or '?' in file_pattern:
                    # Handle wildcards
                    matched_files = glob.glob(file_pattern)
                    group_files.extend(matched_files)
                else:
                    # Direct file path
                    group_files.append(file_pattern)
        elif self.args.group_dir:
            # Directory provided - find all .txt files
            group_dir = Path(self.args.group_dir)
            if group_dir.exists():
                group_files = list(group_dir.glob("*.txt"))
            else:
                print(f"Error: Group directory not found: {group_dir}")
                return []
        else:
            print("Error: Either --group-files or --group-dir must be specified")
            return []
        
        # Convert to Path objects and validate
        validated_files = []
        for file_path in group_files:
            path_obj = Path(file_path)
            if path_obj.exists():
                validated_files.append(path_obj)
            else:
                print(f"Warning: Group file not found: {file_path}")
        
        return validated_files

    def combine_group_files(self):
        """Extract and combine variants from all group files."""
        print("Step 1: Extracting and combining variants from group files...")
        
        group_files = self.find_group_files()
        if not group_files:
            print("No valid group files found!")
            return None
        
        all_variants = []
        variants_per_file = {}
        
        for file_path in group_files:
            print(f"  Processing: {file_path.name}")
            variants = self.extract_variants_from_file(file_path)
            variants_per_file[file_path.name] = variants
            all_variants.extend(variants)
            print(f"    Found {len(variants)} variants")
        
        # Remove duplicates while preserving order
        unique_variants = list(dict.fromkeys(all_variants))
        
        # Write combined variants file
        output_file = self.output_dir / "all_group_file_variants_combined.txt"
        with open(output_file, "w") as file:
            for variant in unique_variants:
                file.write(variant + "\n")
        
        print(f"  Created: {output_file}")
        print(f"  Total unique variants: {len(unique_variants)} (from {len(all_variants)} total)")
        
        return output_file

    def create_plink_files(self):
        """Create PLINK keep files for cases and controls."""
        print("Step 2: Creating PLINK keep files from phenotype data...")
        
        phenotype_file = self.args.phenotype_file
        case_column = self.args.case_column
        id_column = self.args.id_column
        
        try:
            # Try different separators
            separators = [',', '\t', ' ']
            df = None
            
            for sep in separators:
                try:
                    df = pd.read_csv(phenotype_file, sep=sep, header='infer')
                    if len(df.columns) > 1:  # Successfully parsed with multiple columns
                        print(f"  Successfully loaded with separator: '{sep}'")
                        break
                except:
                    continue
            
            if df is None:
                raise ValueError("Could not parse phenotype file with any separator")
            
            print(f"  Loaded phenotype file with {len(df)} individuals and {len(df.columns)} columns")
            print(f"  Columns: {list(df.columns)}")
            
            # Validate required columns exist
            if case_column not in df.columns:
                print(f"Available columns: {list(df.columns)}")
                raise ValueError(f"Case column '{case_column}' not found in phenotype file")
            if id_column not in df.columns:
                print(f"Available columns: {list(df.columns)}")
                raise ValueError(f"ID column '{id_column}' not found in phenotype file")
            
            # Show case/control distribution
            case_counts = df[case_column].value_counts()
            print(f"  Case/control distribution: {dict(case_counts)}")
            
            # Create cases file
            df_cases = df[df[case_column] == 1]
            df_cases_plink = pd.DataFrame({
                'FID': df_cases[id_column],
                'IID': df_cases[id_column]
            })
            
            cases_file = self.output_dir / "hcc_cases_fid_iid_plink_keep.txt"
            df_cases_plink.to_csv(cases_file, sep='\t', index=False, header=True)
            print(f"  Created cases file: {cases_file} ({len(df_cases_plink)} individuals)")
            
            # Create controls file
            df_controls = df[df[case_column] == 0]
            df_controls_plink = pd.DataFrame({
                'FID': df_controls[id_column],
                'IID': df_controls[id_column]
            })
            
            controls_file = self.output_dir / "cancer_free_controls_fid_iid_plink_keep.txt"
            df_controls_plink.to_csv(controls_file, sep='\t', index=False, header=True)
            print(f"  Created controls file: {controls_file} ({len(df_controls_plink)} individuals)")
            
            return cases_file, controls_file
            
        except Exception as e:
            print(f"Error processing phenotype file: {e}")
            sys.exit(1)

    def check_plink_availability(self):
        """Check if PLINK2 is available in the system."""
        try:
            result = subprocess.run(['plink2', '--version'], 
                                  capture_output=True, text=True, check=True)
            version_info = result.stdout.strip().split('\n')[0]
            print(f"  PLINK2 found: {version_info}")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("  Warning: PLINK2 not found in PATH")
            print("  Please ensure PLINK2 is installed and accessible")
            return False

    def run_plink_commands(self, variants_file, cases_file, controls_file):
        """Run PLINK commands to generate geno-counts."""
        print("Step 3: Running PLINK commands...")
        
        # Check PLINK availability first
        if not self.check_plink_availability():
            print("  PLINK2 not available. Cannot proceed with genotype counting.")
            return None, None
        
        bfile = self.args.bfile
        
        # PLINK command for cases
        cases_output = self.output_dir / "cases_geno_counts"
        cmd_cases = [
            "plink2",
            "--bfile", str(bfile),
            "--extract", str(variants_file),
            "--keep", str(cases_file),
            "--geno-counts",
            "--out", str(cases_output)
        ]
        
        # PLINK command for controls
        controls_output = self.output_dir / "controls_geno_counts"
        cmd_controls = [
            "plink2",
            "--bfile", str(bfile),
            "--extract", str(variants_file),
            "--keep", str(controls_file),
            "--geno-counts",
            "--out", str(controls_output)
        ]
        
        print("  Running PLINK for cases...")
        print(f"    Command: {' '.join(cmd_cases)}")
        
        try:
            result_cases = subprocess.run(cmd_cases, capture_output=True, text=True, check=True)
            print(f"    Cases PLINK completed successfully")
        except subprocess.CalledProcessError as e:
            print(f"    Error running PLINK for cases:")
            print(f"    Return code: {e.returncode}")
            if e.stderr:
                print(f"    STDERR: {e.stderr}")
            return None, None
        except FileNotFoundError:
            print(f"    Error: PLINK2 executable not found")
            return None, None
        
        print("  Running PLINK for controls...")
        print(f"    Command: {' '.join(cmd_controls)}")
        
        try:
            result_controls = subprocess.run(cmd_controls, capture_output=True, text=True, check=True)
            print(f"    Controls PLINK completed successfully")
        except subprocess.CalledProcessError as e:
            print(f"    Error running PLINK for controls:")
            print(f"    Return code: {e.returncode}")
            if e.stderr:
                print(f"    STDERR: {e.stderr}")
            return None, None
        except FileNotFoundError:
            print(f"    Error: PLINK2 executable not found")
            return None, None
        
        # Verify output files were created
        cases_gcount_file = f"{cases_output}.gcount"
        controls_gcount_file = f"{controls_output}.gcount"
        
        if not os.path.exists(cases_gcount_file):
            print(f"    Error: Expected output file not found: {cases_gcount_file}")
            return None, None
        
        if not os.path.exists(controls_gcount_file):
            print(f"    Error: Expected output file not found: {controls_gcount_file}")
            return None, None
        
        print(f"    Generated files:")
        print(f"      Cases: {cases_gcount_file}")
        print(f"      Controls: {controls_gcount_file}")
        
        return cases_gcount_file, controls_gcount_file

    def calculate_mac(self, gcount_file):
        """Calculate MAC from a PLINK .gcount file."""
        try:
            df = pd.read_csv(gcount_file, sep="\t")
            
            # Calculate allele counts
            df['REF_allele_count'] = (2 * df['HOM_REF_CT']) + df['HET_REF_ALT_CTS']
            df['ALT_allele_count'] = df['HET_REF_ALT_CTS'] + (2 * df['TWO_ALT_GENO_CTS'])
            
            # Calculate MAC (minimum of REF and ALT allele counts)
            df['MAC'] = df[['REF_allele_count', 'ALT_allele_count']].min(axis=1)
            
            return df
        except Exception as e:
            print(f"Error processing {gcount_file}: {e}")
            return None

    def process_mac_results(self, cases_gcount_file, controls_gcount_file):
        """Process MAC results and create final summary."""
        print("Step 4: Calculating MAC and creating summary...")
        
        # Calculate MAC for cases and controls
        cases_df = self.calculate_mac(cases_gcount_file)
        if cases_df is None:
            return None
            
        controls_df = self.calculate_mac(controls_gcount_file)
        if controls_df is None:
            return None
        
        # Rename MAC columns
        cases_df['MAC_Case'] = cases_df['MAC']
        controls_df['MAC_Control'] = controls_df['MAC']
        
        # Merge results
        try:
            merged_df = pd.merge(
                cases_df[['#CHROM', 'ID', 'REF', 'ALT', 'MAC_Case']],
                controls_df[['#CHROM', 'ID', 'REF', 'ALT', 'MAC_Control']],
                on=['#CHROM', 'ID', 'REF', 'ALT']
            )
            
            # Save results
            output_file = self.output_dir / "MAC_case_control_summary.txt"
            merged_df.to_csv(output_file, sep='\t', index=False)
            
            print(f"  Created final summary: {output_file}")
            print(f"  Summary contains {len(merged_df)} variants")
            
            # Print some summary statistics
            print("\nSummary Statistics:")
            print(f"  Mean MAC in cases: {merged_df['MAC_Case'].mean():.2f}")
            print(f"  Mean MAC in controls: {merged_df['MAC_Control'].mean():.2f}")
            print(f"  Variants with MAC > 5 in cases: {sum(merged_df['MAC_Case'] > 5)}")
            print(f"  Variants with MAC > 5 in controls: {sum(merged_df['MAC_Control'] > 5)}")
            
            return output_file
            
        except Exception as e:
            print(f"Error merging results: {e}")
            return None

    def run_pipeline(self):
        """Run the complete analysis pipeline."""
        print("Starting Variant MAC Pipeline...")
        print("="*50)
        
        # Step 1: Combine group files
        variants_file = self.combine_group_files()
        if variants_file is None:
            return False
        
        # Step 2: Create PLINK files
        cases_file, controls_file = self.create_plink_files()
        
        if self.args.skip_plink:
            print("\nSkipping PLINK execution as requested")
            print("Files created:")
            print(f"  - {variants_file}")
            print(f"  - {cases_file}")
            print(f"  - {controls_file}")
            return True
        
        # Step 3: Run PLINK commands
        cases_gcount, controls_gcount = self.run_plink_commands(variants_file, cases_file, controls_file)
        
        if cases_gcount is None or controls_gcount is None:
            print("Pipeline failed at PLINK step")
            return False
        
        # Step 4: Process MAC results
        final_output = self.process_mac_results(cases_gcount, controls_gcount)
        
        if final_output is None:
            print("Pipeline failed at MAC processing step")
            return False
        
        print("="*50)
        print("Pipeline completed successfully!")
        print(f"Final output: {final_output}")
        return True


def check_dependencies():
    """Check if required dependencies are available."""
    try:
        import pandas
        print("✓ pandas available")
        return True
    except ImportError:
        print("✗ pandas not found - please install with: pip install pandas")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Variant MAC Pipeline - Extract variants, process phenotypes, run PLINK, and calculate MAC",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with specific group files
  python3 variant_mac_pipeline.py \\
    --group-files cancer_genes.*.txt \\
    --phenotype-file phenotypes.csv \\
    --bfile /path/to/plink_files \\
    --output-dir results/

  # Using a directory of group files
  python3 variant_mac_pipeline.py \\
    --group-dir /path/to/group/files/ \\
    --phenotype-file phenotypes.csv \\
    --bfile /path/to/plink_files \\
    --output-dir results/

  # Test without running PLINK
  python3 variant_mac_pipeline.py \\
    --group-files *.txt \\
    --phenotype-file phenotypes.csv \\
    --bfile /path/to/plink_files \\
    --output-dir results/ \\
    --skip-plink

  # Custom column names
  python3 variant_mac_pipeline.py \\
    --group-files cancer_genes.*.txt \\
    --phenotype-file phenotypes.csv \\
    --bfile /path/to/plink_files \\
    --output-dir results/ \\
    --case-column disease_status \\
    --id-column subject_id
        """
    )
    
    # Group files - either specific files or directory
    group_args = parser.add_mutually_exclusive_group(required=True)
    group_args.add_argument(
        "--group-files", 
        nargs="+", 
        help="Group files containing variants (supports wildcards like *.txt)"
    )
    group_args.add_argument(
        "--group-dir", 
        help="Directory containing group files (will use all .txt files)"
    )
    
    # Required arguments
    parser.add_argument(
        "--phenotype-file", 
        required=True,
        help="Phenotype file (CSV/TSV) with case/control status"
    )
    parser.add_argument(
        "--bfile", 
        required=True,
        help="PLINK binary file prefix (without .bed/.bim/.fam extension)"
    )
    parser.add_argument(
        "--output-dir", 
        required=True,
        help="Output directory for results"
    )
    
    # Optional arguments with defaults
    parser.add_argument(
        "--case-column", 
        default="HCC_Cancer_Free",
        help="Column name for case/control status (default: HCC_Cancer_Free)"
    )
    parser.add_argument(
        "--id-column", 
        default="person_id",
        help="Column name for individual IDs (default: person_id)"
    )
    
    # Flags
    parser.add_argument(
        "--skip-plink", 
        action="store_true",
        help="Skip PLINK execution (for testing)"
    )
    parser.add_argument(
        "--version", 
        action="version", 
        version="Variant MAC Pipeline v1.0"
    )
    
    args = parser.parse_args()
    
    # Show what we're doing
    print("Variant MAC Pipeline")
    print("=" * 50)
    print("Settings:")
    if args.group_files:
        print(f"  Group files: {args.group_files}")
    else:
        print(f"  Group directory: {args.group_dir}")
    print(f"  Phenotype file: {args.phenotype_file}")
    print(f"  PLINK bfile: {args.bfile}")
    print(f"  Output directory: {args.output_dir}")
    print(f"  Case column: {args.case_column}")
    print(f"  ID column: {args.id_column}")
    print(f"  Skip PLINK: {args.skip_plink}")
    print()
    
    # Check dependencies
    if not check_dependencies():
        sys.exit(1)
    
    # Run pipeline
    try:
        pipeline = VariantMACPipeline(args)
        success = pipeline.run_pipeline()
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        print("\nPipeline interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\nUnexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()


