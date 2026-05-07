import random
import pysam
import collections
import argparse
import sys
import os

"""
SYNTHETIC MUTATION GENERATOR FOR GENOMIC REGIONS
-----------------------------------------------
This script generates random mutations (SNPs, Substitutions, Indels) 
within specified BED regions. It ensures that mutations are placed 
only on valid nucleotide bases (A, C, G, T) by checking the reference genome.

The output is a standard VCF v4.2 file.
"""


def get_random_bases(length):
    """Generates a random DNA sequence string of a specified length."""
    return ''.join(random.choice(['A', 'C', 'G', 'T']) for _ in range(length))


def generate_distributed_indels(genome, regions, all_mutations, count, min_len, max_len, max_per_chrom=None):
    """
    Generates random insertions and deletions across specified regions.

    Args:
        genome: pysam.FastaFile object.
        regions: List of (chrom, start, end) tuples.
        all_mutations: List to append generated mutation tuples to.
        count: Total number of indels to generate.
        min_len: Minimum length of the indel.
        max_len: Maximum length of the indel.
        max_per_chrom: Limit of how many indels can occur on a single chromosome.
    """
    mutations_added = 0
    chrom_stats = collections.defaultdict(int)

    while mutations_added < count:
        chrom, start, end = random.choice(regions)

        # Enforce per-chromosome limits if specified
        if max_per_chrom and chrom_stats[chrom] >= max_per_chrom:
            available = [r for r in regions if chrom_stats[r[0]] < max_per_chrom]
            if not available:
                # Increment limit slightly to prevent infinite loop if regions are exhausted
                max_per_chrom += 1
                continue
            chrom, start, end = random.choice(available)

        length = random.randint(min_len, max_len)
        is_insertion = random.random() > 0.5

        if is_insertion:
            pos = random.randint(start + 1, end)
            anchor = genome.fetch(chrom, pos - 1, pos).upper()
            if anchor in ['A', 'C', 'G', 'T']:
                # VCF Insertion: REF is anchor, ALT is anchor + new sequence
                all_mutations.append(
                    (chrom, pos, anchor, anchor + get_random_bases(length), f"SVTYPE=INS;LEN={length}"))
                chrom_stats[chrom] += 1
                mutations_added += 1
        else:
            # Ensure the deletion doesn't exceed region boundaries
            if (end - start) < (length + 1):
                continue
            pos = random.randint(start + 1, end - length)
            ref_check = genome.fetch(chrom, pos - 1, pos + length).upper()
            if len(ref_check) == (length + 1) and all(b in ['A', 'C', 'G', 'T'] for b in ref_check):
                # VCF Deletion: REF is anchor + deleted sequence, ALT is just anchor
                all_mutations.append((chrom, pos, ref_check, ref_check[0], f"SVTYPE=DEL;LEN={length}"))
                chrom_stats[chrom] += 1
                mutations_added += 1


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic mutations in VCF format for specific BED regions.")

    # File Arguments
    parser.add_argument("-r", "--reference", required=True, help="Path to reference FASTA file")
    parser.add_argument("-b", "--bed", required=True, help="Path to input BED file defining target regions")
    parser.add_argument("-o", "--output", required=True, help="Path for the generated output VCF")

    # Mutation Count Arguments
    parser.add_argument("--snps", type=int, default=40, help="Number of SNPs to generate")
    parser.add_argument("--sub2", type=int, default=20, help="Number of 2-bp substitutions to generate")
    parser.add_argument("--short_indels", type=int, default=20, help="Number of short indels (3-50bp) to generate")
    parser.add_argument("--long_indels", type=int, default=20,
                        help="Number of long structural indels (51-4000bp) to generate")

    args = parser.parse_args()

    # Load Genome
    if not os.path.exists(args.reference):
        print(f"Error: Reference file {args.reference} not found.")
        sys.exit(1)
    genome = pysam.FastaFile(args.reference)

    # Parse BED
    regions = []
    with open(args.bed, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.split()
                regions.append((parts[0], int(parts[1]), int(parts[2])))

    all_mutations = []

    # 1. Generate SNPs
    print(f"Generating {args.snps} SNPs...")
    for _ in range(args.snps):
        while True:
            chrom, start, end = random.choice(regions)
            pos = random.randint(start + 1, end)
            ref_base = genome.fetch(chrom, pos - 1, pos).upper()
            if ref_base in ['A', 'C', 'G', 'T']:
                alt_base = random.choice([b for b in ['A', 'C', 'G', 'T'] if b != ref_base])
                all_mutations.append((chrom, pos, ref_base, alt_base, "."))
                break

    # 2. Generate 2-bp Substitutions
    print(f"Generating {args.sub2} 2-bp substitutions...")
    for _ in range(args.sub2):
        while True:
            chrom, start, end = random.choice(regions)
            pos = random.randint(start + 1, end - 1)
            ref_seq = genome.fetch(chrom, pos - 1, pos + 1).upper()
            if len(ref_seq) == 2 and all(b in ['A', 'C', 'G', 'T'] for b in ref_seq):
                while True:
                    alt_seq = get_random_bases(2)
                    if alt_seq[0] != ref_seq[0] and alt_seq[1] != ref_seq[1]:
                        break
                all_mutations.append((chrom, pos, ref_seq, alt_seq, "."))
                break

    # 3. Generate Short Indels (3-50 bp)
    print(f"Generating {args.short_indels} short indels...")
    generate_distributed_indels(genome, regions, all_mutations, args.short_indels, 3, 50)

    # 4. Generate Long Indels (51-4000 bp, Max 1 per chromosome)
    print(f"Generating {args.long_indels} long structural indels...")
    generate_distributed_indels(genome, regions, all_mutations, args.long_indels, 51, 4000, max_per_chrom=1)

    # Sort mutations by chromosome and position for valid VCF output
    all_mutations.sort(key=lambda x: (x[0], x[1]))

    # Write VCF Output
    print(f"Writing results to {args.output}...")
    with open(args.output, 'w') as f:
        f.write("##fileformat=VCFv4.2\n")
        for name, length in zip(genome.references, genome.lengths):
            f.write(f"##contig=<ID={name},length={length}>\n")
        f.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
        f.write("##INFO=<ID=LEN,Number=1,Type=Integer,Description=\"Length of variant\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for m in all_mutations:
            f.write(f"{m[0]}\t{m[1]}\t.\t{m[2]}\t{m[3]}\t999\tPASS\t{m[4]}\n")

    print(f"Done! Successfully generated {len(all_mutations)} mutations.")


if __name__ == "__main__":
    main()