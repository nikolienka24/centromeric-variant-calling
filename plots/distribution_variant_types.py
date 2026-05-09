import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import argparse
import os
import sys

"""
VARIANT CONCORDANCE VISUALIZER
------------------------------
This script analyzes the output of a consensus variant caller to visualize
how many tools (2 or 3) agree on specific types of genomic variants.

Features:
1. Classifies variants into SNPs, Substitutions, and Indels (Short/Long).
2. Generates a stacked bar chart showing tool concordance.
3. Uses a high-contrast 'Magma' color palette for scientific clarity.
"""


def parse_args():
    """Parses and returns command-line arguments."""
    parser = argparse.ArgumentParser(description="Visualize variant concordance across multiple tools.")
    parser.add_argument("-i", "--input", required=True, help="Path to the combined consensus TSV file")
    parser.add_argument("-o", "--output", required=True, help="Path to save the generated PNG plot")
    return parser.parse_args()


def classify_variant(row):
    """Categorizes variants based on the lengths of REF and ALT sequences."""
    ref_seq = str(row['ref_seq']) if pd.notna(row['ref_seq']) else ''
    alt_seq = str(row['alt_seq']) if pd.notna(row['alt_seq']) else ''

    ref_len = len(ref_seq)
    alt_len = len(alt_seq)

    if ref_len == 0 and alt_len == 0:
        return 'Unknown'

    if ref_len == 1 and alt_len == 1:
        return 'SNP'
    elif ref_len == alt_len:
        return 'Substitution'
    else:
        diff = abs(ref_len - alt_len)
        if diff <= 50:
            return 'Indel (≤50bp)'
        else:
            return 'Indel (>50bp)'


def main():
    args = parse_args()

    # Validate input file
    if not os.path.exists(args.input):
        print(f"ERROR: Input file {args.input} not found.", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(args.input, sep='\t')

    if df.empty:
        print("ERROR: Input file is empty.", file=sys.stderr)
        sys.exit(1)

    required_cols = {'ref_seq', 'alt_seq', 'tools_count'}
    missing = required_cols - set(df.columns)
    if missing:
        print(f"ERROR: Input file is missing required columns: {missing}", file=sys.stderr)
        sys.exit(1)

    # Classify variants and prepare plot data
    print("Classifying variants and calculating concordance...")
    df['variant_type'] = df.apply(classify_variant, axis=1)

    plot_data = df.groupby(['variant_type', 'tools_count']).size().unstack(fill_value=0)

    for level in [2, 3]:
        if level not in plot_data.columns:
            plot_data[level] = 0
    plot_data = plot_data.reindex(columns=[2, 3])

    # Generate stacked bar plot
    # 0.3 is dark purple (2 tools), 0.8 is bright orange (3 tools)
    magma_colors = [cm.magma(0.3), cm.magma(0.8)]

    print("Generating stacked bar plot...")
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_data.plot(kind='bar', stacked=True, ax=ax, color=magma_colors)

    ax.set_title('Variant Concordance by Tool Agreement', fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel('Variant Type', fontsize=12)
    ax.set_ylabel('Count of Variants', fontsize=12)
    ax.legend(title='Tool Agreement', labels=['2 Tools', '3 Tools'], frameon=True)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    ax.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()

    # Save output
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    plt.savefig(args.output, dpi=300, bbox_inches='tight')

    print("-" * 50)
    print("SUCCESS")
    print(f"Total variants analyzed: {len(df)}")
    print(f"Plot saved to: {args.output}")
    print("-" * 50)


if __name__ == "__main__":
    main()