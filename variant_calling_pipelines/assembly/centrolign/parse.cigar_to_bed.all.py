#!/usr/bin/env python3

"""
CIGAR to BED Converter
Description: Parses a CIGAR string and converts alignment operations into
             a BED-like format tracking both Reference and Mutation coordinates.
Usage: python3 cigar_to_bed.py -i <input.cigar> -o <output.bed> --chrom <name>
"""

import re
import sys
import argparse


def cigar_to_bed_robust(cigar_path, output_bed, chrom_name):
    # --- 1. DATA LOADING ---
    try:
        with open(cigar_path, 'r') as f:
            cigar = f.read().strip()
    except FileNotFoundError:
        print(f"ERROR: CIGAR file not found: {cigar_path}")
        sys.exit(1)

    # Regex finding all operations: digits followed by a valid CIGAR character
    operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar)

    ref_pos = 1  # 1-based start
    mut_pos = 1

    # --- 2. CIGAR PARSING LOGIC ---
    print(f"Parsing CIGAR for {chrom_name}...")

    with open(output_bed, 'w') as out:
        # Standard BED-style header with additional coordinate tracking
        out.write("#Chrom\tRef_Start\tRef_End\tMut_Start\tMut_End\tOp\tLength\tDescription\n")

        for length, op in operations:
            length = int(length)
            r_start, m_start = ref_pos, mut_pos

            # CASE 1: Consumes both (Alignment / Match / Mismatch)
            if op in ('M', '=', 'X'):
                ref_pos += length
                mut_pos += length
                out.write(
                    f"{chrom_name}\t{r_start}\t{ref_pos - 1}\t{m_start}\t{mut_pos - 1}\t{op}\t{length}\talignment\n")

            # CASE 2: Consumes Reference only (Deletion / Gap)
            elif op in ('D', 'N'):
                ref_pos += length
                out.write(f"{chrom_name}\t{r_start}\t{ref_pos - 1}\t{mut_pos}\t{mut_pos}\t{op}\t{length}\tgap_in_mut\n")

            # CASE 3: Consumes Mutated only (Insertion / Soft-clip)
            elif op in ('I', 'S'):
                mut_pos += length
                out.write(
                    f"{chrom_name}\t{ref_pos}\t{ref_pos}\t{m_start}\t{mut_pos - 1}\t{op}\t{length}\tadded_to_mut\n")

            # CASE 4: Consumes neither (Hard-clip / Padding)
            elif op in ('H', 'P'):
                # These don't move the pointers, but are logged for completeness
                out.write(
                    f"{chrom_name}\t{ref_pos}\t{ref_pos}\t{mut_pos}\t{mut_pos}\t{op}\t{length}\tclipped_or_padded\n")

    print(f"Conversion complete. Results saved to: {output_bed}")


def main():
    # --- 3. ARGUMENT PARSING ---
    parser = argparse.ArgumentParser(description='Convert CIGAR strings to an informative BED file.')
    parser.add_argument('-i', '--input', required=True, help='Input CIGAR text file')
    parser.add_argument('-o', '--output', required=True, help='Output BED file')
    parser.add_argument('--chrom', required=True, help='Chromosome name for the output records')

    args = parser.parse_args()

    cigar_to_bed_robust(
        args.input,
        args.output,
        args.chrom
    )


if __name__ == "__main__":
    main()