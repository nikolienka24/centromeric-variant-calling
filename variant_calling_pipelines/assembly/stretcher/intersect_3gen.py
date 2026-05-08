#!/usr/bin/env python3

"""
De Novo Variant Filter
Description: Compares two alignments involving Gen2.
             Identifies mutations present in the Gen1-to-Gen2 alignment
             that are NOT present in the Gen3-to-Gen2 alignment.
Usage: python3 filter_denovo.py -gp <gen1.bedpe> -gd <gen2.bedpe> -o <output_dir>
"""

import pandas as pd
import argparse
import os
import sys

def main():
    # --- 1. ARGUMENT PARSING ---
    parser = argparse.ArgumentParser(description='Filter unique mutations between generations.')
    parser.add_argument('-gp', '--input_gp', required=True, help='Path to Generation 1 BEDPE file')
    parser.add_argument('-gd', '--input_gd', required=True, help='Path to Generation 2 BEDPE file')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory to save the final_de_novo.tsv')
    args = parser.parse_args()

    # Ensure output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    try:
        # --- 2. DATA LOADING ---
        print(f"Loading files...")
        gp = pd.read_csv(args.input_gp, sep='\t')
        gd = pd.read_csv(args.input_gd, sep='\t')

        # Check if required column exists
        if 'start1' not in gp.columns or 'start1' not in gd.columns:
            print("ERROR: Column 'start1' not found in one of the input files.")
            sys.exit(1)

        # --- 3. FILTERING LOGIC ---
        # Create a set of unique positions from the granddaughter (GD)
        gd_positions = set(gd['start1'])

        # Keep rows where start1 in GP is NOT present in GD positions
        result = gp[~gp['start1'].isin(gd_positions)]

        # --- 4. SAVING RESULTS ---
        output_file = os.path.join(args.output_dir, 'final_de_novo.tsv')
        result.to_csv(output_file, sep='\t', index=False)

        # --- 5. TERMINAL SUMMARY ---
        print("-" * 40)
        print(f"Processing Complete")
        print(f"Generation 1 total:   {len(gp)}")
        print(f"Unique mutations:     {len(result)}")
        print(f"Result saved to:      {output_file}")
        print("-" * 40)

    except Exception as e:
        print(f"RUNTIME ERROR: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()