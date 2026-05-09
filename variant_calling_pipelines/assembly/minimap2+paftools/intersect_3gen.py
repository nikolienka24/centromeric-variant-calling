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


def parse_args():
    """Parses and returns command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Filter unique mutations by subtracting next-generation variants.')
    parser.add_argument('-v1', '--vs_gen1', required=True,
                        help='Variants from Gen1-vs-Gen2 (e.g., vs_gp.filtered.tsv)')
    parser.add_argument('-v2', '--vs_gen2', required=True,
                        help='Variants from Gen3-vs-Gen2 (e.g., vs_gd.filtered.tsv)')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Directory to save the final_de_novo.tsv')
    return parser.parse_args()


def main():
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load input files
    try:
        print("Loading variant files...")
        df_gp = pd.read_csv(args.vs_gen1, sep='\t')
        df_gd = pd.read_csv(args.vs_gen2, sep='\t')
    except Exception as e:
        print(f"ERROR: Failed to read input files: {e}", file=sys.stderr)
        sys.exit(1)

    # Validate presence of coordinate column
    if 'POS' not in df_gp.columns or 'POS' not in df_gd.columns:
        print("ERROR: Column 'POS' not found. Ensure the input files are correctly formatted TSVs.",
              file=sys.stderr)
        sys.exit(1)

    # Keep variants from Gen1-vs-Gen2 whose POS is NOT present in Gen3-vs-Gen2
    next_gen_positions = set(df_gd['POS'])
    result = df_gp[~df_gp['POS'].isin(next_gen_positions)]

    # Save results
    output_file = os.path.join(args.output_dir, 'final_de_novo.tsv')
    result.to_csv(output_file, sep='\t', index=False)

    print("-" * 45)
    print("De Novo Filtering Complete")
    print(f"Total variants in Gen1 vs Gen2:  {len(df_gp)}")
    print(f"Variants seen in next gen:       {len(df_gd)}")
    print(f"Final De Novo count:             {len(result)}")
    print(f"Result saved to:                 {output_file}")
    print("-" * 45)


if __name__ == "__main__":
    main()