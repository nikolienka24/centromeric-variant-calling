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


def parse_args():
    """Parses and returns command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Filter unique mutations between generations.')
    parser.add_argument('-gp', '--input_gp', required=True,
                        help='Path to Generation 1 BEDPE file')
    parser.add_argument('-gd', '--input_gd', required=True,
                        help='Path to Generation 2 BEDPE file')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Directory to save the final_de_novo.tsv')
    return parser.parse_args()


def main():
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load input files
    try:
        print("Loading files...")
        gp = pd.read_csv(args.input_gp, sep='\t')
        gd = pd.read_csv(args.input_gd, sep='\t')
    except Exception as e:
        print(f"ERROR: Failed to read input files: {e}", file=sys.stderr)
        sys.exit(1)

    # Validate presence of coordinate column
    if 'start1' not in gp.columns or 'start1' not in gd.columns:
        print("ERROR: Column 'start1' not found. Ensure the inputs are Stretcher BEDPE files.",
              file=sys.stderr)
        sys.exit(1)

    # Keep rows from GP whose start1 is NOT present in GD
    gd_positions = set(gd['start1'])
    result = gp[~gp['start1'].isin(gd_positions)]

    # Save results
    output_file = os.path.join(args.output_dir, 'final_de_novo.tsv')
    result.to_csv(output_file, sep='\t', index=False)

    print("-" * 40)
    print("Processing Complete")
    print(f"Generation 1 total:   {len(gp)}")
    print(f"Unique mutations:     {len(result)}")
    print(f"Result saved to:      {output_file}")
    print("-" * 40)


if __name__ == "__main__":
    main()