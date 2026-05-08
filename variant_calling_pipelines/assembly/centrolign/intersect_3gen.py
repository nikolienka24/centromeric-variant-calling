#!/usr/bin/env python3

"""
Centrolign Mutation Comparison Tool
Description: Compares mutation results from two generations. Identifies
             mutations present in the first generation (GP) that are NOT
             present in the second generation (GD) based on Ref_Start.
Usage: python3 compare_mutations.py -g1 <input_gp.tsv> -g2 <input_gd.tsv> -o <output.tsv>
"""

import pandas as pd
import argparse
import os
import sys

def main():
    # --- 1. ARGUMENT PARSING ---
    parser = argparse.ArgumentParser(description='Compare mutation positions between two generations.')
    parser.add_argument('-g1', '--input_gp', required=True, help='Path to Generation 1 (GP) TSV results')
    parser.add_argument('-g2', '--input_gd', required=True, help='Path to Generation 2 (GD) TSV results')
    parser.add_argument('-o', '--output', required=True, help='Path for the filtered output TSV')
    args = parser.parse_args()

    # --- 2. DATA LOADING ---
    try:
        print(f"Loading files...")
        gp = pd.read_csv(args.input_gp, sep='\t')
        gd = pd.read_csv(args.input_gd, sep='\t')
    except Exception as e:
        print(f"RUNTIME ERROR: Failed to read input files: {e}")
        sys.exit(1)

    # Validate presence of coordinate column
    if 'Ref_Start' not in gp.columns or 'Ref_Start' not in gd.columns:
        print("ERROR: Column 'Ref_Start' not found. Ensure the inputs are Centrolign result TSVs.")
        sys.exit(1)

    # --- 3. FILTERING LOGIC ---
    # Create a set of positions from the daughter/granddaughter (GD) for exclusion
    gd_starts = set(gd['Ref_Start'])

    # Keep rows from the parent (GP) only if Ref_Start is NOT in the exclusion set
    result = gp[~gp['Ref_Start'].isin(gd_starts)]

    # --- 4. SAVING RESULTS ---
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    result.to_csv(args.output, sep='\t', index=False)

    # --- 5. TERMINAL SUMMARY ---
    print("-" * 50)
    print(f"Mutation Comparison Complete")
    print(f"GP Source:          {os.path.basename(args.input_gp)}")
    print(f"GD Source:          {os.path.basename(args.input_gd)}")
    print("-" * 50)
    print(f"Total mutations (GP):        {len(gp)}")
    print(f"Total mutations (GD):        {len(gd)}")
    print(f"Unique mutations (De Novo):  {len(result)}")
    print("-" * 50)
    print(f"Result saved to: {args.output}")

if __name__ == "__main__":
    main()