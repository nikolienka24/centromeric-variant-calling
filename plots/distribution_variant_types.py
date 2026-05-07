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


def classify_variant(row):
    """Categorizes variants based on the lengths of REF and ALT sequences."""
    # Handle potential NaN values
    ref_seq = str(row['ref_seq']) if pd.notna(row['ref_seq']) else ''
    alt_seq = str(row['alt_seq']) if pd.notna(row['alt_seq']) else ''

    ref_len = len(ref_seq)
    alt_len = len(alt_seq)

    # Basic safety check for empty data
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
    parser = argparse.ArgumentParser(description="Visualize variant concordance across multiple tools.")
    parser.add_argument("-i", "--input", required=True, help="Path to the combined consensus TSV file")
    parser.add_argument("-o", "--output", required=True, help="Path to save the generated PNG plot")
    args = parser.parse_args()

    # --- 1. Data Loading ---
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} not found.")
        sys.exit(1)

    df = pd.read_csv(args.input, sep='\t')

    # --- 2. Processing & Classification ---
    print("Classifying variants and calculating concordance...")
    df['variant_type'] = df.apply(classify_variant, axis=1)

    # Prepare data for stacked bar plot (Tool counts are columns, types are rows)
    plot_data = df.groupby(['variant_type', 'tools_count']).size().unstack(fill_value=0)

    # Ensure both concordance levels (2 and 3) are represented
    for level in [2, 3]:
        if level not in plot_data.columns:
            plot_data[level] = 0

    plot_data = plot_data.reindex(columns=[2, 3])

    # --- 3. Visualization ---
    # Select discrete colors from the Magma colormap
    # 0.3 is dark purple (2 tools), 0.8 is bright orange (3 tools)
    magma_colors = [cm.magma(0.3), cm.magma(0.8)]

    print(f"Generating stacked bar plot...")
    plt.figure(figsize=(10, 6))
    ax = plot_data.plot(kind='bar', stacked=True, figsize=(10, 6), color=magma_colors)

    # Styling the plot
    plt.title('Variant Concordance by Tool Agreement', fontsize=14, fontweight='bold', pad=15)
    plt.xlabel('Variant Type', fontsize=12)
    plt.ylabel('Count of Variants', fontsize=12)
    plt.legend(title='Tool Agreement', labels=['2 Tools', '3 Tools'], frameon=True)

    plt.xticks(rotation=0)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()

    # --- 4. Export ---
    # Ensure output directory exists
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