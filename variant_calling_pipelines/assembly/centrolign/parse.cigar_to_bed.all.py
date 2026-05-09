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


def parse_args():
    """Parses and returns command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Convert CIGAR strings to an informative BED file.')
    parser.add_argument('-i', '--input', required=True, help='Input CIGAR text file')
    parser.add_argument('-o', '--output', required=True, help='Output BED file')
    parser.add_argument('--chrom', required=True, help='Chromosome name for the output records')
    return parser.parse_args()


def cigar_to_bed(cigar_path, output_bed, chrom_name):
    """
    Parses a CIGAR string file and writes alignment operations to a BED-like file
    tracking both reference and mutation coordinates.
    """
    try:
        with open(cigar_path, 'r') as f:
            cigar = f.read().strip()
    except FileNotFoundError:
        print(f"ERROR: CIGAR file not found: {cigar_path}", file=sys.stderr)
        sys.exit(1)

    operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar)

    if not operations:
        print(f"ERROR: No valid CIGAR operations found in {cigar_path}. "
              f"File may be empty or malformed.", file=sys.stderr)
        sys.exit(1)

    ref_pos = 1  # 1-based start
    mut_pos = 1

    print(f"Parsing CIGAR for {chrom_name}...")

    with open(output_bed, 'w') as out:
        out.write("#Chrom\tRef_Start\tRef_End\tMut_Start\tMut_End\tOp\tLength\tDescription\n")

        for length, op in operations:
            length = int(length)
            r_start, m_start = ref_pos, mut_pos

            if op in ('M', '=', 'X'):
                # Consumes both reference and mutation (alignment / match / mismatch)
                ref_pos += length
                mut_pos += length
                out.write(
                    f"{chrom_name}\t{r_start}\t{ref_pos - 1}\t{m_start}\t{mut_pos - 1}"
                    f"\t{op}\t{length}\talignment\n")

            elif op in ('D', 'N'):
                # Consumes reference only (deletion / gap)
                ref_pos += length
                out.write(
                    f"{chrom_name}\t{r_start}\t{ref_pos - 1}\t{mut_pos}\t{mut_pos}"
                    f"\t{op}\t{length}\tgap_in_mut\n")

            elif op in ('I', 'S'):
                # Consumes mutation only (insertion / soft-clip)
                mut_pos += length
                out.write(
                    f"{chrom_name}\t{ref_pos}\t{ref_pos}\t{m_start}\t{mut_pos - 1}"
                    f"\t{op}\t{length}\tadded_to_mut\n")

            elif op in ('H', 'P'):
                # Consumes neither (hard-clip / padding) — logged for completeness
                out.write(
                    f"{chrom_name}\t{ref_pos}\t{ref_pos}\t{mut_pos}\t{mut_pos}"
                    f"\t{op}\t{length}\tclipped_or_padded\n")

    print(f"Conversion complete. Results saved to: {output_bed}")


def main():
    args = parse_args()
    cigar_to_bed(args.input, args.output, args.chrom)


if __name__ == "__main__":
    main()