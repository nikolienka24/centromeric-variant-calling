#!/usr/bin/env python3

"""
De Novo Variant Filter (Minimap2/Paftools Version)
Description: Compares variant calls from two generations.
             Identifies mutations present in the Gen1-to-Gen2 comparison
             that were NOT passed down to (or seen in) the Gen3-to-Gen2 comparison.
Usage: python3 filter_variants.py -v1 <vs_gp.vcf> -v2 <vs_gd.vcf> -o <output_dir>
"""

import pandas as pd
import argparse
import os
import sys

def main():
    # --- 1. ARGUMENT PARSING ---
    parser = argparse.ArgumentParser(description='Filter unique mutations by subtracting next-generation variants.')
    parser.add_argument('-v1', '--vs_gen1', required=True, help='Variants from Gen1-vs-Gen2 (e.g., vs_gp.filtered.tsv)')
    parser.add_argument('-v2', '--vs_gen2', required=True, help='Variants from Gen3-vs-Gen2 (e.g., vs_gd.filtered.tsv)')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory to save the final_de_novo.tsv')
    args = parser.parse_args()

    # Ensure output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    try:
        # --- 2. DATA LOADING ---
        print(f"Loading variant files...")
        # Note: Using sep='\t' as Paftools/VCF based TSVs are tab-delimited
        df_gp = pd.read_csv(args.vs_gen1, sep='\t')
        df_gd = pd.read_csv(args.vs_gen2, sep='\t')

        # Validation: Check for the POS column (standard for Paftools VCF-to-TSV)
        if 'POS' not in df_gp.columns or 'POS' not in df_gd.columns:
            print("ERROR: Column 'POS' not found. Ensure the input files are correctly formatted TSVs.")
            sys.exit(1)

        # --- 3. FILTERING LOGIC ---
        # Goal: Keep variants from Gen1-vs-Gen2, subtract any POS found in Gen3-vs-Gen2
        # Using a set for O(1) lookup efficiency
        next_gen_positions = set(df_gd['POS'])

        # Keep rows where the POS is NOT in the next generation's set
        result = df_gp[~df_gp['POS'].isin(next_gen_positions)]

        # --- 4. SAVING RESULTS ---
        output_file = os.path.join(args.output_dir, 'final_de_novo.tsv')
        result.to_csv(output_file, sep='\t', index=False)

        # --- 5. TERMINAL SUMMARY ---
        print("-" * 45)
        print(f"De Novo Filtering Complete")
        print(f"Total variants in Gen1 vs Gen2:  {len(df_gp)}")
        print(f"Variants seen in next gen:      {len(df_gd)}")
        print(f"Final De Novo count:            {len(result)}")
        print(f"Result saved to:                {output_file}")
        print("-" * 45)

    except Exception as e:
        print(f"RUNTIME ERROR: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()