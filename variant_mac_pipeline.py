#!/usr/bin/env python3

"""
Variant Mac Pipeline 

Pipeline to calculate minor allele counts (MAC) and carrier frequencies for cases 
and controls from variants in group files, with gnomAD allele frequencies from VEP 
annotation files and variant annotations derived from group file classifications.

Also includes ancestry-stratified carrier frequnecy tables to address reviewer 
requests for frequency of each variant category (pLoF, damaging_missense, combines)
in cases and controls per ancestry subgroup.

EXAMPLE USAGE:
python3 variant_mac_pipeline.py \
  --group-files /home/agaro/verma_shared/projects/HCC/PMBB_6Gene_Burden_Analysis_030625/group_files/cancer_genes.*.txt \
  --vep-files /home/agaro/verma_shared/projects/HCC/PMBB_6Gene_Burden_Analysis_030625/variant_annotations/vep_annot_msh6_chr_2_releasev3.txt /home/agaro/verma_shared/projects/HCC/PMBB_6Gene_Burden_Analysis_030625/variant_annotations/vep_annot_pms2_chr_7_releasev3.txt \
  --phenotype-file /project/verma_shared/projects/HCC/covariates/covariates_hcc_cancer_free_controls_exomepcs_ancestry_column_3.0.csv \
  --bfile /static/PMBB/PMBB-Release-2024-3.0/Exome/pVCF/NF_all_variants/PMBB-Release-2024-3.0_genetic_exome_NF \
  --output-dir /home/agaro/verma_shared/projects/HCC/PMBB_6Gene_Burden_Analysis_030625/ALL_ALL/Variant_Level_Analysis_Test \
  --case-column HCC_Cancer_Free \
  --id-column person_id \
  --biobank PMBB \
  --ancestry-column ancestry \
  --clinvar-file /path/to/Clinvar_inclusion_list.txt


This pipeline:
1. Extracts variants from group files and combines them 
2. Creates PLINK keep files for cases (HCC_Cancer_Free=1) and controls (HCC_Cancer_Free=0), including ancestry-stratiied keep files for each subgroup
3. Runs PLINK2 genotype counting for cases and controls
4. Calculates Minor Allele Counts (MAC) for both groups
5. Adds gnomAD allele frequency column and variant annotation columns
6. Adds ClinVar clinical significance based on annotation type
7. Creates final summary: MAC_case_control_summary.txt
8. Creates per-gene carrier frequency tables stratified by ancestry and case/controls: (gene_carrier_frequency_table.tsv and gene_carrier_frequency_wide.tsv)

"""

import pandas as pd
import argparse
import sys
from pathlib import Path 
import subprocess
import os 
import glob 
import re
from collections import defaultdict

class VariantMACPipeline:
    def __init__(self, args):
        self.args = args
        self.output_dir = Path(args.output_dir)
        self.output_dir.mkdir(parents = True, exist_ok= True)

        # Store the group file to annotation mapping
        self.group_annotations = {}
        #Store variant -> gene name mapping
        self.variant_to_gene = {}
        # Store VEP data for lookup 
        self.vep_data = {}
        # Store ClinVar data for lookup
        self.clinvar_data = {}
        # Store VEP files by chromosome for CLIN_SIG lookup
        self.vep_files_by_chr = {}

    # ──────────────────────────────────────────────────────────────────────────
    # Group file parsing
    # ──────────────────────────────────────────────────────────────────────────

    def extract_gene_from_filename(self, file_path):
        """Extract gene name from a group file name.

        Handles naming patterns like:
          cancer_genes.BRCA2.txt  ->  BRCA2
          cancer_genes.chr2.txt   ->  MSH6  (via CHR_TO_GENE lookup)
          BRCA2_groupfile.txt     ->  BRCA2_groupfile
        Falls back to the full stem if no dot-separated gene can be identified.
        """
        CHR_TO_GENE = {
            '2':  'MSH6',
            '7':  'PMS2',
            '13': 'BRCA2',
            '16': 'FANCA',
            '17': 'BRIP1',
            '22': 'CHEK2',
        }

        stem = Path(file_path).stem          # e.g. "cancer_genes.chr2"
        parts = stem.split('.')
        if len(parts) >= 2:
            token = parts[-1]                # rightmost dot-separated token
            # If the token looks like a chromosome identifier, map it to a gene name
            chr_match = re.match(r'^chr[_\s]?(\d+|X|Y|MT)$', token, re.IGNORECASE)
            if chr_match:
                chr_num = chr_match.group(1)
                return CHR_TO_GENE.get(chr_num, token)
            if token.isdigit():
                return CHR_TO_GENE.get(token, token)
            return token
        return stem


    def extract_variants_from_file(self, file_path):
        """Extract variants from the group file after the 'var' keyword and map to annotations and gene."""
        try:
            with open(file_path, "r") as file:
                lines = file.readlines()

                if len(lines) < 2: 
                    print(f"Warning: Expected at least 2 lines in {file_path}")
                    return [], {}
                
                # Parse first line (gene name and variants)
                first_line = lines[0].strip()
                first_parts = first_line.split() 

                # Parse second line (annotations)
                second_line = lines[1].strip()
                second_parts = second_line.split()

                variants = []
                variant_annotations = {}

                # Find 'var' and extract everything after it
                if "var" in first_parts:
                    var_index = first_parts.index("var")
                    variants = first_parts[var_index + 1:]
                    
                    # Find 'anno' and extract everything after it
                    if "anno" in second_parts:
                        anno_index = second_parts.index("anno")
                        annotations = second_parts[anno_index + 1:]

                        # Map variants to their annotations and gene 
                        if len(variants) == len(annotations):
                            for variant, annotation in zip(variants, annotations):
                               variant_annotations[variant] = annotation
                            print(f"  Mapped {len(variants)} variants to annotations in {file_path.name}")
                        else:
                            print(f"Warning: Mismatch in {file_path}: {len(variants)} variants vs {len(annotations)} annotations")
                            return [], {}
                    else:
                        print(f"Warning: 'anno' keyword not found in {file_path}")
                        return [], {}
                else:
                    print(f"Warning: 'var' keyword not found in {file_path}")
                    return [], {}
            
                return variants, variant_annotations
            
        except FileNotFoundError:
            print(f"Error: File not found - {file_path}")
            return [], {}
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            return [], {} 

    # ──────────────────────────────────────────────────────────────────────────
    # VEP / ClinVar loading
    # ──────────────────────────────────────────────────────────────────────────

    def load_vep_file(self, vep_file):  
        """Load a single VEP file and add to lookup dictionary."""
        vep_path = Path(vep_file)
        if not vep_path.exists():
            print(f"Warning: VEP file not found: {vep_file}")
            return
        
        # Extract chromosomes from filename for later CLIN_SIG lookup
        chr_match = re.search(r'chr[_\s]?(\d+|X|Y|MT)', vep_path.name, re.IGNORECASE)
        if chr_match:
            chr_num = chr_match.group(1)
            self.vep_files_by_chr[chr_num] = vep_path
            print(f"  Mapped chromosome {chr_num} to {vep_path.name}")

        try:
            # Try different separators
            separators = ['\t', ',']
            vep_df = None
            
            for sep in separators:
                try:
                    vep_df = pd.read_csv(vep_path, sep=sep, low_memory=False)
                    if len(vep_df.columns) > 5:  # Reasonable number of columns for VEP
                        break
                except:
                    continue
            
            if vep_df is None:
                print(f"Warning: Could not parse VEP file: {vep_file}")
                return
            
            print(f"  Loaded {vep_file}: {len(vep_df)} rows, {len(vep_df.columns)} columns")
            
            # Determine which gnomAD column to use based on sequencing type
            gnomad_column = self.get_gnomad_column()
            
            # Handle column names that may have # prefix
            uploaded_var_col = None
            if 'Uploaded_variation' in vep_df.columns:
                uploaded_var_col = 'Uploaded_variation'
            elif '#Uploaded_variation' in vep_df.columns:
                uploaded_var_col = '#Uploaded_variation'

            if uploaded_var_col is None:
                print(f" Warning: Could not find Uploaded_variation column in {vep_file}")
                print(f"  Available columns: {list(vep_df.columns)}")
                return
            
            #Check for gnomAD column 
            if gnomad_column not in vep_df.columns:
                print(f" Warning: Missing column in {vep_file}: {gnomad_column}")
                print(f"  Available columns: {list(vep_df.columns)}")
                # Try to find alternative gnomAD columns
                available_gnomad_cols = [col for col in vep_df.columns if 'gnomad' in col.lower()]
                if available_gnomad_cols: 
                    print(f" Available gnomAD columns: {available_gnomad_cols}")
                return
            
            # Store VEP data with Uploaded_variation as key
            for _, row in vep_df.iterrows():
                variant_id = row[uploaded_var_col]
                self.vep_data[variant_id] = {
                    'gnomAD_AF': row[gnomad_column],  
                    'source_file': vep_path.name,
                    'sequencing_type': self.args.sequencing_type,
                    'CLIN_SIG': row.get('CLIN_SIG', '-')
                }

            print(f"  Added {len(vep_df)} variant annotations from {vep_path.name}")
            
        except Exception as e:
            print(f"Error loading VEP file {vep_file}: {e}")
    
    def load_all_vep_data(self):
        """Load all VEP files and create lookup dictionary."""
        if not self.args.vep_files:
            print("No VEP files provided - skipping VEP annotation")
            return
        
        # Determine which gnomAD column to use based on sequencing type
        gnomad_column = self.get_gnomad_column()
        print(f"Loading VEP files for annotation lookup (using {gnomad_column})...")
        
        for vep_file in self.args.vep_files:
            self.load_vep_file(vep_file)
        
        print(f"Total VEP annotations loaded: {len(self.vep_data)}")
    
    def load_clinvar_data(self):
        """Load ClinVar file for pLoF variant clinical significance lookup."""
        if not self.args.clinvar_file:
            print("No ClinVar file provided - skipping ClinVar annotations")
            return
        
        clinvar_path = Path(self.args.clinvar_file)
        if not clinvar_path.exists():
            print(f"Warning: ClinVar file not found: {self.args.clinvar_file}")
            return
        
        print(f"Loading ClinVar file for clinical significance lookup...")
        
        try:
            # Try different separators
            separators = ['\t', ',']
            clinvar_df = None
            
            for sep in separators:
                try:
                    clinvar_df = pd.read_csv(clinvar_path, sep=sep, low_memory=False)
                    if len(clinvar_df.columns) > 5:  # Reasonable number of columns
                        print(f"  Successfully loaded ClinVar with separator: '{sep}'")
                        break
                except:
                    continue
            
            if clinvar_df is None:
                print(f"Warning: Could not parse ClinVar file: {self.args.clinvar_file}")
                return
            
            print(f"  Loaded ClinVar: {len(clinvar_df)} rows, {len(clinvar_df.columns)} columns")
            
            # Check for required columns
            if 'ID' not in clinvar_df.columns:
                print(f"Warning: 'ID' column not found in ClinVar file")
                print(f"  Available columns: {list(clinvar_df.columns)}")
                return
            
            if 'ClinVar_202407_GRCh38.clinical_significance_ordered' not in clinvar_df.columns:
                print(f"Warning: 'ClinVar_202407_GRCh38.clinical_significance_ordered' column not found")
                print(f"  Available columns: {list(clinvar_df.columns)}")
                return
            
            # Store ClinVar data with ID as key
            for _, row in clinvar_df.iterrows():
                variant_id = row['ID']
                self.clinvar_data[variant_id] = {
                    'clinical_significance': row['ClinVar_202407_GRCh38.clinical_significance_ordered']
                }
            
            print(f"  Added {len(self.clinvar_data)} ClinVar annotations")
            
        except Exception as e:
            print(f"Error loading ClinVar file: {e}")

    def get_gnomad_column(self):
        """Determine which gnomAD column to use based on sequencing type."""
        if self.args.sequencing_type.lower() == 'wes':
            return 'gnomADe_AF'  # Exome allele frequencies
        elif self.args.sequencing_type.lower() == 'wgs':
            return 'gnomADg_AF'  # Genome allele frequencies
        else:
            print(f"Warning: Unknown sequencing type '{self.args.sequencing_type}'. Defaulting to WES (gnomADe_AF)")
            return 'gnomADe_AF'

    # ──────────────────────────────────────────────────────────────────────────
    # Group file discovery and combination
    # ──────────────────────────────────────────────────────────────────────────

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
        """Extract and combine variants from all group files with annotation mapping.

        Also populates self.variant_to_gene so that the frequency table step
        can aggregate by gene name.
        """
        print("Step 1: Extracting and combining variants from group files...")

        group_files = self.find_group_files()
        if not group_files:
            print("No valid group files found!")
            return None

        all_variants = []
        variants_per_file = {}

        for file_path in group_files:
            print(f"  Processing: {file_path.name}")
            variants, variant_annotations = self.extract_variants_from_file(file_path)

            gene_name = self.extract_gene_from_filename(file_path)

            variants_per_file[file_path.name] = {
                'variants': variants,
                'annotations': variant_annotations,
                'gene': gene_name
            }

            for variant, annotation in variant_annotations.items():
                self.group_annotations[variant] = annotation
                self.variant_to_gene[variant] = gene_name   # NEW

            all_variants.extend(variants)

            if variant_annotations:
                anno_counts = {}
                for anno in variant_annotations.values():
                    anno_counts[anno] = anno_counts.get(anno, 0) + 1
                print(f"    Gene: {gene_name} | {len(variants)} variants | annotations: {dict(anno_counts)}")
            else:
                print(f"    Gene: {gene_name} | {len(variants)} variants (no annotations)")

        unique_variants = list(dict.fromkeys(all_variants))

        output_file = self.output_dir / "all_group_file_variants_combined.txt"
        with open(output_file, "w") as file:
            for variant in unique_variants:
                file.write(variant + "\n")

        print(f"  Created: {output_file}")
        print(f"  Total unique variants: {len(unique_variants)} (from {len(all_variants)} total)")

        annotation_file = self.output_dir / "variant_annotation_mapping.txt"
        with open(annotation_file, "w") as file:
            file.write("Variant_ID\tGene\tAnnotation\tSource_File\n")
            for variant, annotation in self.group_annotations.items():
                source_file = "unknown"
                gene = self.variant_to_gene.get(variant, "unknown")
                for fname, fdata in variants_per_file.items():
                    if variant in fdata.get('annotations', {}):
                        source_file = fname
                        break
                file.write(f"{variant}\t{gene}\t{annotation}\t{source_file}\n")
        print(f"  Created annotation mapping: {annotation_file}")

        overall_anno_counts = {}
        for anno in self.group_annotations.values():
            overall_anno_counts[anno] = overall_anno_counts.get(anno, 0) + 1
        print(f"  Overall annotation distribution: {dict(overall_anno_counts)}")

        return output_file
    
    # ──────────────────────────────────────────────────────────────────────────
    # PLINK keep-file creation (overall + ancestry-stratified)
    # ──────────────────────────────────────────────────────────────────────────

    def create_plink_files(self):
        """Create PLINK keep files for cases/controls (overall and per ancestry).

        Returns
        -------
        cases_file, controls_file : Path
            Overall case/control keep files (unchanged from original).
        ancestry_keep_files : dict
            Keys are (ancestry, 'cases'|'controls'), values are (Path, int N).
        """
        print("Step 2: Creating PLINK keep files from phenotype data...")

        phenotype_file = self.args.phenotype_file
        case_column = self.args.case_column
        id_column = self.args.id_column
        ancestry_column = self.args.ancestry_column

        try:
            separators = [',', '\t', ' ']
            df = None

            for sep in separators:
                try:
                    df = pd.read_csv(phenotype_file, sep=sep, header='infer')
                    if len(df.columns) > 1:
                        print(f"  Successfully loaded with separator: '{sep}'")
                        break
                except Exception:
                    continue

            if df is None:
                raise ValueError("Could not parse phenotype file with any separator")

            print(f"  Loaded phenotype file with {len(df)} individuals and {len(df.columns)} columns")
            print(f"  Columns: {list(df.columns)}")

            if case_column not in df.columns:
                raise ValueError(f"Case column '{case_column}' not found in phenotype file")
            if id_column not in df.columns:
                raise ValueError(f"ID column '{id_column}' not found in phenotype file")

            case_counts = df[case_column].value_counts()
            print(f"  Case/control distribution: {dict(case_counts)}")

            # ── Overall cases ────────────────────────────────────────────────
            df_cases = df[df[case_column] == 1]
            df_cases_plink = pd.DataFrame({'FID': df_cases[id_column], 'IID': df_cases[id_column]})
            cases_file = self.output_dir / "hcc_cases_fid_iid_plink_keep.txt"
            df_cases_plink.to_csv(cases_file, sep='\t', index=False, header=True)
            print(f"  Created cases file: {cases_file} ({len(df_cases_plink)} individuals)")

            # ── Overall controls ─────────────────────────────────────────────
            df_controls = df[df[case_column] == 0]
            df_controls_plink = pd.DataFrame({'FID': df_controls[id_column], 'IID': df_controls[id_column]})
            controls_file = self.output_dir / "cancer_free_controls_fid_iid_plink_keep.txt"
            df_controls_plink.to_csv(controls_file, sep='\t', index=False, header=True)
            print(f"  Created controls file: {controls_file} ({len(df_controls_plink)} individuals)")

            # ── Ancestry-stratified keep files ───────────────────────────────
            ancestry_keep_files = {}

            if ancestry_column and ancestry_column in df.columns:
                ancestries = sorted(df[ancestry_column].dropna().unique())
                print(f"  Ancestry groups found: {list(ancestries)}")

                for ancestry in ancestries:
                    for status, label in [(1, 'cases'), (0, 'controls')]:
                        subset = df[(df[ancestry_column] == ancestry) & (df[case_column] == status)]
                        if len(subset) == 0:
                            continue
                        plink_df = pd.DataFrame({'FID': subset[id_column], 'IID': subset[id_column]})
                        fname = self.output_dir / f"{ancestry}_{label}_plink_keep.txt"
                        plink_df.to_csv(fname, sep='\t', index=False, header=True)
                        ancestry_keep_files[(ancestry, label)] = (fname, len(subset))
                        print(f"    Created: {fname.name} ({len(subset)} individuals)")
            else:
                if ancestry_column:
                    print(f"  Warning: Ancestry column '{ancestry_column}' not found — skipping ancestry stratification")
                else:
                    print("  No --ancestry-column provided — skipping ancestry stratification")

            return cases_file, controls_file, ancestry_keep_files

        except Exception as e:
            print(f"Error processing phenotype file: {e}")
            sys.exit(1)

    # ──────────────────────────────────────────────────────────────────────────
    # PLINK execution
    # ──────────────────────────────────────────────────────────────────────────

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

    def _run_single_plink_gcount(self, variants_file, keep_file, output_prefix):
        """Run one plink2 --geno-counts call. Returns .gcount path or None on failure."""
        cmd = [
            "plink2",
            "--bfile", str(self.args.bfile),
            "--extract", str(variants_file),
            "--keep", str(keep_file),
            "--geno-counts",
            "--out", str(output_prefix)
        ]
        print(f"    Command: {' '.join(cmd)}")
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"    PLINK error (return code {e.returncode})")
            if e.stderr:
                print(f"    STDERR: {e.stderr[:500]}")
            return None
        except FileNotFoundError:
            print("    Error: plink2 executable not found")
            return None

        gcount_file = f"{output_prefix}.gcount"
        if not os.path.exists(gcount_file):
            print(f"    Error: Expected output not found: {gcount_file}")
            return None
        return gcount_file

    def run_plink_commands(self, variants_file, cases_file, controls_file, ancestry_keep_files=None):
        """Run PLINK geno-counts for overall cases/controls and all ancestry strata."""
        print("Step 3: Running PLINK commands...")

        if not self.check_plink_availability():
            print("  PLINK2 not available. Cannot proceed.")
            return None, None, {}

        # Overall
        print("  Running PLINK for all cases...")
        cases_gcount = self._run_single_plink_gcount(
            variants_file, cases_file, self.output_dir / "cases_geno_counts"
        )

        print("  Running PLINK for all controls...")
        controls_gcount = self._run_single_plink_gcount(
            variants_file, controls_file, self.output_dir / "controls_geno_counts"
        )

        if cases_gcount is None or controls_gcount is None:
            return None, None, {}

        # Ancestry-stratified
        ancestry_gcount_files = {}
        if ancestry_keep_files:
            print("  Running PLINK for ancestry-stratified groups...")
            for (ancestry, label), (keep_file, n) in ancestry_keep_files.items():
                print(f"    {ancestry} {label}...")
                prefix = self.output_dir / f"{ancestry}_{label}_geno_counts"
                gcount = self._run_single_plink_gcount(variants_file, keep_file, prefix)
                if gcount:
                    ancestry_gcount_files[(ancestry, label)] = (gcount, n)
                else:
                    print(f"    Warning: PLINK failed for {ancestry} {label} — skipping")

        return cases_gcount, controls_gcount, ancestry_gcount_files

    # ──────────────────────────────────────────────────────────────────────────
    # MAC calculation 
    # ──────────────────────────────────────────────────────────────────────────

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

    def _carrier_counts_from_gcount(self, gcount_file):
        """Load a .gcount file and add n_carriers and N columns.

        n_carriers = individuals carrying ≥1 ALT allele (HET + HOM_ALT)
        N          = total genotyped individuals in this group
        """
        try:
            df = pd.read_csv(gcount_file, sep="\t")
            df['n_carriers'] = df['HET_REF_ALT_CTS'] + df['TWO_ALT_GENO_CTS']
            df['N'] = df['HOM_REF_CT'] + df['HET_REF_ALT_CTS'] + df['TWO_ALT_GENO_CTS']
            return df
        except Exception as e:
            print(f"Error loading {gcount_file}: {e}")
            return None

    def _aggregate_carriers_by_gene(self, gcount_file, group_n_override=None):
        """Aggregate carrier counts by (gene, annotation) for one group.

        Parameters
        ----------
        gcount_file : str
            Path to a PLINK .gcount file.
        group_n_override : int or None
            If provided, use this as N instead of the per-row N from gcount
            (useful when PLINK may drop individuals with missing genotypes).

        Returns
        -------
        dict : {gene: {annotation: {'n_carriers': int, 'N': int}}}
        """
        df = self._carrier_counts_from_gcount(gcount_file)
        if df is None:
            return {}

        results = {}

        for _, row in df.iterrows():
            variant_id = row['ID']
            gene = self.variant_to_gene.get(variant_id, 'unknown')
            annotation = self.group_annotations.get(variant_id, 'unknown')
            n_carriers = int(row['n_carriers'])
            N = group_n_override if group_n_override is not None else int(row['N'])

            if gene not in results:
                results[gene] = {}
            if annotation not in results[gene]:
                results[gene][annotation] = {'n_carriers': 0, 'N': N}

            results[gene][annotation]['n_carriers'] += n_carriers
            # Take the maximum N seen (should be constant across variants in same group)
            results[gene][annotation]['N'] = max(results[gene][annotation]['N'], N)

        return results
    
    def create_frequency_summary_table(self, cases_gcount, controls_gcount,
                                       ancestry_gcount_files,
                                       cases_n, controls_n):
        """Create per-gene carrier frequency tables stratified by ancestry.

        Produces two output files:
          gene_carrier_frequency_table.tsv  — long format (one row per gene/group/category)
          gene_carrier_frequency_wide.tsv   — wide format matching Table 1 layout

        Parameters
        ----------
        cases_gcount, controls_gcount : str
            Overall gcount file paths.
        ancestry_gcount_files : dict
            {(ancestry, label): (gcount_path, N)}
        cases_n, controls_n : int
            Total N for overall cases/controls (from phenotype file).
        """
        print("Step 5 (NEW): Building per-gene carrier frequency table...")

        ANNOTATIONS = ['pLoF', 'damaging_missense']

        # ── Collect all groups ───────────────────────────────────────────────
        # Each entry: (group_label, gcount_file, N)
        groups = [
            ('All_Cases',    cases_gcount,    cases_n),
            ('All_Controls', controls_gcount, controls_n),
        ]
        for (ancestry, label), (gcount_path, n) in sorted(ancestry_gcount_files.items()):
            group_label = f"{ancestry}_{'Cases' if label == 'cases' else 'Controls'}"
            groups.append((group_label, gcount_path, n))

        # ── Aggregate per group ──────────────────────────────────────────────
        all_group_data = {}   # group_label -> {gene -> {anno -> {n_carriers, N}}}
        genes_seen = set()

        for group_label, gcount_file, n in groups:
            agg = self._aggregate_carriers_by_gene(gcount_file, group_n_override=n)
            all_group_data[group_label] = agg
            genes_seen.update(agg.keys())

        genes_sorted = sorted(g for g in genes_seen if g != 'unknown')
        if 'unknown' in genes_seen:
            genes_sorted.append('unknown')

        # ── Long-format output ───────────────────────────────────────────────
        long_rows = []
        for group_label, _, _ in groups:
            agg = all_group_data.get(group_label, {})
            for gene in genes_sorted:
                gene_data = agg.get(gene, {})
                for anno in ANNOTATIONS:
                    anno_data = gene_data.get(anno, {'n_carriers': 0, 'N': 0})
                    n_c = anno_data['n_carriers']
                    N   = anno_data['N']
                    freq = (n_c / N) if N > 0 else 'NA'
                    long_rows.append({
                        'Biobank':     self.args.biobank,
                        'Gene':        gene,
                        'Group':       group_label,
                        'Category':    anno,
                        'n_carriers':  n_c,
                        'N':           N,
                        'Frequency':   round(freq, 6) if freq != 'NA' else 'NA'
                    })

                # Combined = pLoF + damaging_missense (sum, matching Table 1 approach)
                plof_n = gene_data.get('pLoF', {'n_carriers': 0})['n_carriers']
                dam_n  = gene_data.get('damaging_missense', {'n_carriers': 0})['n_carriers']
                N_combined = gene_data.get('pLoF', gene_data.get('damaging_missense', {'N': 0})).get('N', 0)
                comb_n = plof_n + dam_n
                freq_comb = (comb_n / N_combined) if N_combined > 0 else 'NA'
                long_rows.append({
                    'Biobank':    self.args.biobank,
                    'Gene':       gene,
                    'Group':      group_label,
                    'Category':   'pLoF_and_damaging',
                    'n_carriers': comb_n,
                    'N':          N_combined,
                    'Frequency':  round(freq_comb, 6) if freq_comb != 'NA' else 'NA'
                })

        long_df = pd.DataFrame(long_rows)
        long_file = self.output_dir / "gene_carrier_frequency_table.tsv"
        long_df.to_csv(long_file, sep='\t', index=False)
        print(f"  Saved long-format frequency table: {long_file}")

        # ── Wide-format output (mirrors Table 1 layout) ──────────────────────
        # Columns: Gene | Group | n_pLoF | pLoF_freq | n_damaging | damaging_freq | n_combined | combined_freq | N
        wide_rows = []
        for group_label, _, _ in groups:
            agg = all_group_data.get(group_label, {})
            for gene in genes_sorted:
                gene_data = agg.get(gene, {})

                plof   = gene_data.get('pLoF',             {'n_carriers': 0, 'N': 0})
                dam    = gene_data.get('damaging_missense', {'n_carriers': 0, 'N': 0})
                N      = max(plof['N'], dam['N'])
                comb_n = plof['n_carriers'] + dam['n_carriers']

                def fmt_freq(n, denom):
                    if denom == 0:
                        return 'NA'
                    return round(n / denom, 6)

                wide_rows.append({
                    'Biobank':        self.args.biobank,
                    'Gene':           gene,
                    'Group':          group_label,
                    'n_pLoF':         plof['n_carriers'],
                    'pLoF_freq':      fmt_freq(plof['n_carriers'], N),
                    'n_damaging':     dam['n_carriers'],
                    'damaging_freq':  fmt_freq(dam['n_carriers'], N),
                    'n_combined':     comb_n,
                    'combined_freq':  fmt_freq(comb_n, N),
                    'N':              N
                })

        wide_df = pd.DataFrame(wide_rows)
        wide_file = self.output_dir / "gene_carrier_frequency_wide.tsv"
        wide_df.to_csv(wide_file, sep='\t', index=False)
        print(f"  Saved wide-format frequency table: {wide_file}")

        # ── Print a preview ──────────────────────────────────────────────────
        print("\nFrequency table preview (All Cases / All Controls):")
        preview = wide_df[wide_df['Group'].isin(['All_Cases', 'All_Controls'])]
        print(preview.to_string(index=False))

        return long_file, wide_file

    # ──────────────────────────────────────────────────────────────────────────
    # ClinVar / clinical significance
    # ──────────────────────────────────────────────────────────────────────────

    def get_clin_sig_from_vep(self, variant_id, chromosome):
        """Get CLIN_SIG from VEP file for damaging_missense variants."""
        # First check if we already have it in memory
        if variant_id in self.vep_data and 'CLIN_SIG' in self.vep_data[variant_id]:
            clin_sig = self.vep_data[variant_id]['CLIN_SIG']
            # Return '-' if it's empty or NA
            if pd.isna(clin_sig) or clin_sig == '' or clin_sig == '-':
                return '-'
            return clin_sig
        
        # If not in memory, return '-'
        return '-'

    def add_clinical_significance(self, merged_df):
        """Add ClinVar_CLIN_SIG column based on annotation type."""
        if not self.args.clinvar_file:
            print("  Skipping ClinVar clinical significance (no ClinVar file provided)")
            return merged_df
        
        print("  Adding ClinVar clinical significance...")
        
        clin_sig_values = []
        plof_not_found = []
        
        for idx, row in merged_df.iterrows():
            variant_id = row['ID']
            annotation = row['Annotation']
            chromosome = str(row['#CHROM'])
            
            if annotation == 'pLoF':
                # Convert variants with ":" separator to "_" separator for ClinVar lookup
                clinvar_variant_id = variant_id.replace(':', '_')
                if clinvar_variant_id in self.clinvar_data:
                    clin_sig = self.clinvar_data[clinvar_variant_id]['clinical_significance']
                    clin_sig_values.append(clin_sig)
                else:
                    clin_sig_values.append('NOT_FOUND_IN_CLINVAR')
                    plof_not_found.append(variant_id)
            
            elif annotation == 'damaging_missense':
                # Look up in VEP CLIN_SIG column
                clin_sig = self.get_clin_sig_from_vep(variant_id, chromosome)
                clin_sig_values.append(clin_sig)
            
            else:
                clin_sig_values.append('UNKNOWN_ANNOTATION')
        
        merged_df['ClinVar_CLIN_SIG'] = clin_sig_values
        
        # Report pLoF variants not found in ClinVar
        if plof_not_found:
            print(f"  WARNING: {len(plof_not_found)} pLoF variants not found in ClinVar file:")
            print(f"           This may indicate ref/alt allele flipping issues")
            not_found_file = self.output_dir / "plof_variants_not_found_in_clinvar.txt"
            with open(not_found_file, 'w') as f:
                f.write("Variant_ID\n")
                for var in plof_not_found:
                    f.write(f"{var}\n")
            print(f"           List saved to: {not_found_file}")
        
        return merged_df

    # ──────────────────────────────────────────────────────────────────────────
    # MAC summary 
    # ──────────────────────────────────────────────────────────────────────────

    def process_mac_results(self, cases_gcount_file, controls_gcount_file):
        """Process MAC results and create final summary with annotations and gnomAD column."""
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
            
            # Add gnomAD column
            gnomad_column_name = f"gnomAD_AF_{self.args.sequencing_type.upper()}"
            print(f"  Adding VEP annotations ({gnomad_column_name})...")
            merged_df[gnomad_column_name] = merged_df['ID'].map(
                lambda x: self.vep_data.get(x, {}).get('gnomAD_AF', 'NA')
            )

            # Add variant annotations
            print("  Adding group-based annotations...")
            merged_df['Annotation'] = merged_df['ID'].map(
                lambda x: self.group_annotations.get(x, 'unknown')
            )
            
            CHR_TO_GENE = {
                '2':  'MSH6',
                '7':  'PMS2',
                '13': 'BRCA2',
                '16': 'FANCA',
                '17': 'BRIP1',
                '22': 'CHEK2',
            }

            merged_df['Gene'] = merged_df['#CHROM'].astype(str).str.replace(
                r'^chr', '', regex=True
            ).map(lambda x: CHR_TO_GENE.get(x, x))


            # Add ClinVar clinical significance
            merged_df = self.add_clinical_significance(merged_df)

            # Reorder columns
            gnomad_column_name = f"gnomAD_AF_{self.args.sequencing_type.upper()}"
            column_order = ['#CHROM', 'ID', 'REF', 'ALT', 'Gene', 'Annotation',
                            'MAC_Case', 'MAC_Control', gnomad_column_name, 'ClinVar_CLIN_SIG']
            merged_df = merged_df[column_order]

            # Save results
            output_file = self.output_dir / "group_file_variants_MAC_case_control_summary.txt"
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
        print("=" * 50)

        self.load_all_vep_data()
        self.load_clinvar_data()

        # Step 1
        variants_file = self.combine_group_files()
        if variants_file is None:
            return False

        # Step 2
        cases_file, controls_file, ancestry_keep_files = self.create_plink_files()

        if self.args.skip_plink:
            print("\nSkipping PLINK execution as requested")
            print("Files created:")
            print(f"  - {variants_file}")
            print(f"  - {cases_file}")
            print(f"  - {controls_file}")
            for (anc, lbl), (fp, n) in ancestry_keep_files.items():
                print(f"  - {fp}  (N={n})")
            return True

        # Step 3
        cases_gcount, controls_gcount, ancestry_gcount_files = self.run_plink_commands(
            variants_file, cases_file, controls_file, ancestry_keep_files
        )

        if cases_gcount is None or controls_gcount is None:
            print("Pipeline failed at PLINK step")
            return False

        # Step 4 — existing MAC summary
        final_output = self.process_mac_results(cases_gcount, controls_gcount)
        if final_output is None:
            print("Pipeline failed at MAC processing step")
            return False

        # Step 5 — NEW: per-gene carrier frequency table
        cases_n    = sum(1 for line in open(cases_file)) - 1     # subtract header
        controls_n = sum(1 for line in open(controls_file)) - 1

        self.create_frequency_summary_table(
            cases_gcount, controls_gcount,
            ancestry_gcount_files,
            cases_n, controls_n
        )

        print("=" * 50)
        print("Pipeline completed successfully!")
        print(f"Variant-level MAC summary : {final_output}")
        print(f"Gene frequency table (long): {self.output_dir / 'gene_carrier_frequency_table.tsv'}")
        print(f"Gene frequency table (wide): {self.output_dir / 'gene_carrier_frequency_wide.tsv'}")
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
        description="Variant MAC Pipeline — MAC, frequencies, and ancestry-stratified carrier tables",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 variant_mac_pipeline.py \\
    --group-files cancer_genes.*.txt \\
    --vep-files vep_gene1.txt vep_gene2.txt \\
    --clinvar-file ClinVar_Inclusion.txt \\
    --sequencing-type WES \\
    --phenotype-file phenotypes.csv \\
    --bfile /path/to/plink_files \\
    --biobank PMBB \\
    --ancestry-column ancestry \\
    --output-dir results/
        """
    )

    group_args = parser.add_mutually_exclusive_group(required=True)
    group_args.add_argument("--group-files", nargs="+",
                            help="Group files containing variants (supports wildcards)")
    group_args.add_argument("--group-dir",
                            help="Directory containing group files (all .txt files)")

    parser.add_argument("--vep-files", nargs="*",
                        help="VEP annotation files (optional)")
    parser.add_argument("--clinvar-file",
                        help="ClinVar file for clinical significance annotation (optional)")
    parser.add_argument("--sequencing-type", choices=['WES', 'WGS', 'wes', 'wgs'],
                        default='WES',
                        help="Sequencing type: WES=gnomADe_AF, WGS=gnomADg_AF (default: WES)")

    parser.add_argument("--phenotype-file", required=True,
                        help="Phenotype file (CSV/TSV) with case/control status")
    parser.add_argument("--bfile", required=True,
                        help="PLINK binary file prefix (no .bed/.bim/.fam)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for results")

    parser.add_argument("--case-column", default="HCC_Cancer_Free",
                        help="Column for case/control status (default: HCC_Cancer_Free)")
    parser.add_argument("--id-column", default="person_id",
                        help="Column for individual IDs (default: person_id)")

    # NEW argument
    parser.add_argument("--biobank", required=True,
                        help="Biobank name to label output rows (e.g. PMBB, AOU, Mayo, MVP)")

    parser.add_argument("--ancestry-column", default="ancestry",
                        help="Column for ancestry labels (e.g. EUR, AFR, AMR) used to "
                             "produce ancestry-stratified frequency tables "
                             "(default: ancestry). Set to empty string to skip.")

    parser.add_argument("--skip-plink", action="store_true",
                        help="Skip PLINK execution (for testing)")
    parser.add_argument("--version", action="version",
                        version="Variant MAC Pipeline v2.0")

    args = parser.parse_args()

    # Normalise: empty string -> None for ancestry column
    if args.ancestry_column == '':
        args.ancestry_column = None

    print("Variant MAC Pipeline v2.0")
    print("=" * 50)
    print("Settings:")
    if args.group_files:
        print(f"  Group files:      {args.group_files}")
    else:
        print(f"  Group directory:  {args.group_dir}")
    if args.vep_files:
        col = 'gnomADe_AF' if args.sequencing_type.upper() == 'WES' else 'gnomADg_AF'
        print(f"  VEP files:        {args.vep_files}")
        print(f"  Sequencing type:  {args.sequencing_type} ({col})")
    else:
        print("  VEP files:        None")
    print(f"  ClinVar file:     {args.clinvar_file or 'None'}")
    print(f"  Phenotype file:   {args.phenotype_file}")
    print(f"  PLINK bfile:      {args.bfile}")
    print(f"  Output directory: {args.output_dir}")
    print(f"  Case column:      {args.case_column}")
    print(f"  ID column:        {args.id_column}")
    print(f"  Biobank:          {args.biobank}")
    print(f"  Ancestry column:  {args.ancestry_column or '(none)'}")
    print(f"  Skip PLINK:       {args.skip_plink}")
    print()

    if not check_dependencies():
        sys.exit(1)

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

