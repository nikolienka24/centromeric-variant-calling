#!/usr/bin/env python3

"""
Variant Comparison Tool
Description: Compares two specific sequences from a transposed matrix
             and identifies positional differences.
"""

import sys
import argparse

def main():
    # --- Argument Parsing ---
    parser = argparse.ArgumentParser(description='Compare two specific rows in a transposed matrix.')
    parser.add_argument('-i', '--input', required=True, help='Input matrix.tsv')
    parser.add_argument('-o', '--output', required=True, help='Output variants.tsv')
    parser.add_argument('-r1', '--row1', required=True, help='Name of the first row (e.g., REF)')
    parser.add_argument('-r2', '--row2', required=True, help='Name of the second row (e.g., MUT)')
    args = parser.parse_args()

    data = {}

    # --- Data Loading ---
    try:
        with open(args.input, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if not parts:
                    continue
                # Key: Row Name, Value: List of sequence values
                data[parts[0]] = parts[1:]

        # Validate existence of rows
        if args.row1 not in data:
            print(f"ERROR: Row '{args.row1}' not found.")
            sys.exit(1)
        if args.row2 not in data:
            print(f"ERROR: Row '{args.row2}' not found.")
            sys.exit(1)

        seq1 = data[args.row1]
        seq2 = data[args.row2]

        # --- Comparison Logic ---
        with open(args.output, 'w') as f_out:
            # Write header
            f_out.write(f"COL_IDX\t{args.row1}\t{args.row2}\n")

            diff_count = 0
            # Iterate through sequences and log mismatches
            for i, (val1, val2) in enumerate(zip(seq1, seq2)):
                if val1 != val2:
                    f_out.write(f"{i}\t{val1}\t{val2}\n")
                    diff_count += 1

        # --- Terminal Output ---
        print("-" * 30)
        print(f"Comparison complete")
        print(f"Reference: {args.row1}")
        print(f"Query:     {args.row2}")
        print(f"Variants:  {diff_count}")
        print("-" * 30)

    except Exception as e:
        print(f"RUNTIME ERROR: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()