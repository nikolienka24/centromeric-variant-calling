#!/usr/bin/env python3

"""
CIGAR to Variant TSV Converter
Description: Parses a CIGAR string and two FASTA sequences to identify
             SNPs, Insertions (INS), and Deletions (DEL) while
             applying genomic offsets.
Usage: python3 cigar_to_tsv.py -r <ref.fa> -m <mut.fa> -c <cigar.txt> \
                               -o <output.tsv> --chrom <name> \
                               --off_ref <int> --off_mut <int>
"""

import re
import sys
import argparse


def load_fasta_seq(path):
    """Loads FASTA sequence, stripping headers and joining lines."""
    try:
        with open(path, 'r') as f:
            return "".join(line.strip() for line in f if not line.startswith(">"))
    except FileNotFoundError:
        print(f"ERROR: File not found: {path}")
        sys.exit(1)


def extract_mutations(ref_path, mut_path, cigar_path, output_file, chrom_name, offset_ref, offset_mut):
    # --- 1. DATA LOADING ---
    print(f"Loading sequences for {chrom_name}...")
    ref_seq = load_fasta_seq(ref_path)
    mut_seq = load_fasta_seq(mut_path)

    try:
        with open(cigar_path, 'r') as f:
            cigar = f.read().strip()
    except FileNotFoundError:
        print(f"ERROR: CIGAR file not found: {cigar_path}")
        sys.exit(1)

    # Extract operations (e.g., 10M, 2D, 5I) using regex
    operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar)

    ref_pos = 0
    mut_pos = 0

    # --- 2. CIGAR PARSING & VARIANT EXTRACTION ---
    with open(output_file, 'w') as out:
        # Header includes Ref_End for easier BED conversion (0-based)
        out.write("Name\tType\tRef_Start\tRef_End\tMut_Pos\tRef_Base\tAlt_Base\tLength\n")

        for length, op in operations:
            length = int(length)

            # SOFT-CLIPPING (S): Advance in query (mut) sequence only
            if op == 'S':
                mut_pos += length
                continue

            # MATCH / MISMATCH (M, =, X): Check for SNPs
            if op in ('M', '=', 'X'):
                r_sub = ref_seq[ref_pos: ref_pos + length]
                m_sub = mut_seq[mut_pos: mut_pos + length]

                for j in range(len(r_sub)):
                    if r_sub[j] != m_sub[j]:
                        # 0-based absolute positions
                        abs_ref_start = ref_pos + j + offset_ref
                        abs_mut_start = mut_pos + j + offset_mut
                        out.write(
                            f"{chrom_name}\tSNP\t{abs_ref_start}\t{abs_ref_start + 1}\t"
                            f"{abs_mut_start}\t{r_sub[j]}\t{m_sub[j]}\t1\n"
                        )

                ref_pos += length
                mut_pos += length

            # DELETION (D): Bases missing in mutated sequence
            elif op == 'D':
                deleted_seq = ref_seq[ref_pos: ref_pos + length]
                abs_ref_start = ref_pos + offset_ref
                abs_mut_start = mut_pos + offset_mut

                out.write(
                    f"{chrom_name}\tDEL\t{abs_ref_start}\t{abs_ref_start + length}\t"
                    f"{abs_mut_start}\t{deleted_seq}\t-\t{length}\n"
                )
                ref_pos += length

            # INSERTION (I): Extra bases in mutated sequence
            elif op == 'I':
                inserted_seq = mut_seq[mut_pos: mut_pos + length]
                abs_ref_start = ref_pos + offset_ref
                abs_mut_start = mut_pos + offset_mut

                # For insertions, Ref_Start == Ref_End (point between bases in BED)
                out.write(
                    f"{chrom_name}\tINS\t{abs_ref_start}\t{abs_ref_start}\t"
                    f"{abs_mut_start}\t-\t{inserted_seq}\t{length}\n"
                )
                mut_pos += length

    print(f"Extraction complete. Results saved to: {output_file}")


def main():
    # --- 3. ARGUMENT PARSING ---
    parser = argparse.ArgumentParser(description='Extract variants from CIGAR alignment.')
    parser.add_argument('-r', '--ref', required=True, help='Reference FASTA')
    parser.add_argument('-m', '--mut', required=True, help='Mutated (Query) FASTA')
    parser.add_argument('-c', '--cigar', required=True, help='Input CIGAR text file')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    parser.add_argument('--chrom', required=True, help='Chromosome name for output')
    parser.add_argument('--off_ref', type=int, default=0, help='Reference genomic offset')
    parser.add_argument('--off_mut', type=int, default=0, help='Mutation genomic offset')

    args = parser.parse_args()

    extract_mutations(
        args.ref,
        args.mut,
        args.cigar,
        args.output,
        args.chrom,
        args.off_ref,
        args.off_mut
    )


if __name__ == "__main__":
    main()