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

VALID_BASES = {'A', 'C', 'G', 'T'}

# Maximum number of attempts to find a valid position before giving up
MAX_ATTEMPTS = 100_000

# INFO field templates
INFO_INS = "SVTYPE=INS;LEN={length}"
INFO_DEL = "SVTYPE=DEL;LEN={length}"


def parse_args():
    """Parses and returns command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate synthetic mutations in VCF format for specific BED regions.")

    # File arguments
    parser.add_argument("-r", "--reference", required=True, help="Path to reference FASTA file")
    parser.add_argument("-b", "--bed", required=True, help="Path to input BED file defining target regions")
    parser.add_argument("-o", "--output", required=True, help="Path for the generated output VCF")

    # Mutation count arguments
    parser.add_argument("--snps", type=int, default=40, help="Number of SNPs to generate")
    parser.add_argument("--sub2", type=int, default=20, help="Number of 2-bp substitutions to generate")
    parser.add_argument("--short_indels", type=int, default=20,
                        help="Number of short indels (3-50 bp) to generate")
    parser.add_argument("--long_indels", type=int, default=20,
                        help="Number of long structural indels (51-4000 bp) to generate")

    return parser.parse_args()


def get_random_bases(length):
    """Generates a random DNA sequence string of a specified length."""
    return ''.join(random.choice(list(VALID_BASES)) for _ in range(length))


def parse_bed(bed_path):
    """
    Parses a BED file and returns a list of (chrom, start, end) tuples.
    Skips empty lines, comments, and track/browser headers.
    """
    regions = []
    with open(bed_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            regions.append((parts[0], int(parts[1]), int(parts[2])))
    return regions


def generate_snps(genome, regions, count):
    """
    Generates random SNPs within the specified regions.

    Args:
        genome: pysam.FastaFile object.
        regions: List of (chrom, start, end) tuples.
        count: Number of SNPs to generate.

    Returns:
        List of mutation tuples and count of successfully generated SNPs.
    """
    mutations = []
    added = 0
    for _ in range(count):
        attempts = 0
        while attempts < MAX_ATTEMPTS:
            attempts += 1
            chrom, start, end = random.choice(regions)
            pos = random.randint(start + 1, end)
            ref_base = genome.fetch(chrom, pos - 1, pos).upper()
            if ref_base in VALID_BASES:
                alt_base = random.choice([b for b in VALID_BASES if b != ref_base])
                mutations.append((chrom, pos, ref_base, alt_base, "."))
                added += 1
                break
        else:
            print(f"WARNING: Could not place all SNPs after {MAX_ATTEMPTS} attempts.", file=sys.stderr)
    return mutations, added


def generate_substitutions(genome, regions, count):
    """
    Generates random 2-bp substitutions within the specified regions.

    Args:
        genome: pysam.FastaFile object.
        regions: List of (chrom, start, end) tuples.
        count: Number of 2-bp substitutions to generate.

    Returns:
        List of mutation tuples and count of successfully generated substitutions.
    """
    mutations = []
    added = 0
    for _ in range(count):
        attempts = 0
        while attempts < MAX_ATTEMPTS:
            attempts += 1
            chrom, start, end = random.choice(regions)
            pos = random.randint(start + 1, end - 1)
            ref_seq = genome.fetch(chrom, pos - 1, pos + 1).upper()
            if len(ref_seq) == 2 and all(b in VALID_BASES for b in ref_seq):
                inner_attempts = 0
                while inner_attempts < MAX_ATTEMPTS:
                    inner_attempts += 1
                    alt_seq = get_random_bases(2)
                    if alt_seq[0] != ref_seq[0] and alt_seq[1] != ref_seq[1]:
                        mutations.append((chrom, pos, ref_seq, alt_seq, "."))
                        added += 1
                        break
                break
        else:
            print(f"WARNING: Could not place all 2-bp substitutions after {MAX_ATTEMPTS} attempts.",
                  file=sys.stderr)
    return mutations, added


def generate_distributed_indels(genome, regions, count, min_len, max_len, max_per_chrom=None):
    """
    Generates random insertions and deletions across specified regions.

    Args:
        genome: pysam.FastaFile object.
        regions: List of (chrom, start, end) tuples.
        count: Total number of indels to generate.
        min_len: Minimum length of the indel.
        max_len: Maximum length of the indel.
        max_per_chrom: Limit of how many indels can occur on a single chromosome.

    Returns:
        List of mutation tuples and count of successfully generated indels.
    """
    mutations = []
    added = 0
    chrom_stats = collections.defaultdict(int)
    attempts = 0

    while added < count:
        attempts += 1
        if attempts > MAX_ATTEMPTS:
            print(f"WARNING: Could only generate {added}/{count} indels "
                  f"(len {min_len}-{max_len} bp) after {MAX_ATTEMPTS} attempts. "
                  f"Regions may be too small or too few.", file=sys.stderr)
            break

        if max_per_chrom:
            available = [r for r in regions if chrom_stats[r[0]] < max_per_chrom]
            if not available:
                print(f"WARNING: Per-chromosome limit reached for all chromosomes. "
                      f"Generated {added}/{count} indels.", file=sys.stderr)
                break
            chrom, start, end = random.choice(available)
        else:
            chrom, start, end = random.choice(regions)

        length = random.randint(min_len, max_len)
        is_insertion = random.random() > 0.5

        if is_insertion:
            pos = random.randint(start + 1, end)
            anchor = genome.fetch(chrom, pos - 1, pos).upper()
            if anchor in VALID_BASES:
                mutations.append(
                    (chrom, pos, anchor, anchor + get_random_bases(length),
                     INFO_INS.format(length=length)))
                chrom_stats[chrom] += 1
                added += 1
        else:
            if (end - start) < (length + 1):
                continue
            pos = random.randint(start + 1, end - length)
            ref_check = genome.fetch(chrom, pos - 1, pos + length).upper()
            if len(ref_check) == (length + 1) and all(b in VALID_BASES for b in ref_check):
                mutations.append(
                    (chrom, pos, ref_check, ref_check[0],
                     INFO_DEL.format(length=length)))
                chrom_stats[chrom] += 1
                added += 1

    return mutations, added


def write_vcf(output_path, mutations, genome):
    """
    Writes sorted mutations to a VCF v4.2 file.

    Args:
        output_path: Path to the output VCF file.
        mutations: List of (chrom, pos, ref, alt, info) tuples.
        genome: pysam.FastaFile object (used for contig headers).
    """
    mutations.sort(key=lambda x: (x[0], x[1]))
    with open(output_path, 'w') as f:
        f.write("##fileformat=VCFv4.2\n")
        for name, length in zip(genome.references, genome.lengths):
            f.write(f"##contig=<ID={name},length={length}>\n")
        f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        f.write('##INFO=<ID=LEN,Number=1,Type=Integer,Description="Length of variant">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for m in mutations:
            f.write(f"{m[0]}\t{m[1]}\t.\t{m[2]}\t{m[3]}\t999\tPASS\t{m[4]}\n")


def main():
    args = parse_args()

    # Validate input files
    if not os.path.exists(args.reference):
        print(f"ERROR: Reference file {args.reference} not found.", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.bed):
        print(f"ERROR: BED file {args.bed} not found.", file=sys.stderr)
        sys.exit(1)

    # Load genome and regions
    genome = pysam.FastaFile(args.reference)
    regions = parse_bed(args.bed)
    if not regions:
        print("ERROR: No valid regions found in BED file.", file=sys.stderr)
        sys.exit(1)
    print(f"Loaded {len(regions)} regions from BED file.")

    # Generate mutations
    print(f"Generating {args.snps} SNPs...")
    snp_muts, snps_added = generate_snps(genome, regions, args.snps)
    print(f"  -> Generated {snps_added}/{args.snps} SNPs.")

    print(f"Generating {args.sub2} 2-bp substitutions...")
    sub_muts, subs_added = generate_substitutions(genome, regions, args.sub2)
    print(f"  -> Generated {subs_added}/{args.sub2} 2-bp substitutions.")

    print(f"Generating {args.short_indels} short indels (3-50 bp)...")
    short_muts, short_added = generate_distributed_indels(
        genome, regions, args.short_indels, 3, 50)
    print(f"  -> Generated {short_added}/{args.short_indels} short indels.")

    print(f"Generating {args.long_indels} long structural indels (51-4000 bp)...")
    long_muts, long_added = generate_distributed_indels(
        genome, regions, args.long_indels, 51, 4000, max_per_chrom=1)
    print(f"  -> Generated {long_added}/{args.long_indels} long indels.")

    # Combine and write
    all_mutations = snp_muts + sub_muts + short_muts + long_muts
    print(f"Writing {len(all_mutations)} mutations to {args.output}...")
    write_vcf(args.output, all_mutations, genome)

    print(f"Done! Successfully wrote {len(all_mutations)} mutations "
          f"(SNPs: {snps_added}, SUBs: {subs_added}, "
          f"short indels: {short_added}, long indels: {long_added}).")


if __name__ == "__main__":
    main()