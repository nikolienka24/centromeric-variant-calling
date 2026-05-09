import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys

"""
PEDIGREE BENCHMARKING ANALYSIS
------------------------------
This script processes resource usage logs (time and memory) for various
bioinformatics tools. It generates three distinct plots to visualize
performance trade-offs and exports a summary statistics table.

Visualizations:
1. Total Runtime (Sum) per tool.
2. Peak Memory (Max) per tool.
3. Efficiency Scatter (Time vs. Memory correlation).
"""


def parse_args():
    """Parses and returns command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate benchmarking plots and stats for tool comparisons.")
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Directory containing tool .tsv performance files")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Directory to save generated plots")
    parser.add_argument("-s", "--stats_out", required=True,
                        help="Path to save the final summary TSV file")
    return parser.parse_args()


def main():
    args = parse_args()

    # Validate input directory
    if not os.path.isdir(args.input_dir):
        print(f"ERROR: Input directory {args.input_dir} not found.", file=sys.stderr)
        sys.exit(1)

    # Mapping tool names to their expected log filenames
    file_map = {
        'Minimap2/Deepvariant': 'map+deepvariant.tsv',
        'Centrolign_Align': 'centrolign_align.tsv',
        'Minimap2/Paftools': 'minimap2+paftools.tsv',
        'Stretcher': 'stretcher_align.tsv',
        'Centrolign_Pang': 'centrolign_pangenomes.tsv'
    }

    all_data = []

    # Load and normalize data
    print(f"Reading performance logs from: {args.input_dir}")
    for tool_name, file_name in file_map.items():
        file_path = os.path.join(args.input_dir, file_name)

        if not os.path.exists(file_path):
            print(f"Warning: File not found: {file_path}", file=sys.stderr)
            continue

        try:
            temp_df = pd.read_csv(file_path, sep='\s+').dropna(how='all')

            if temp_df.empty:
                print(f"Warning: {file_name} is empty.", file=sys.stderr)
                continue

            # Standardize time column names
            if 'Time_sec' in temp_df.columns:
                temp_df = temp_df.rename(columns={'Time_sec': 'RealTime_sec'})

            if 'RealTime_sec' not in temp_df.columns or 'MaxMem_MB' not in temp_df.columns:
                print(f"Warning: Missing required columns in {file_name}.", file=sys.stderr)
                continue

            temp_df['Tool'] = tool_name
            all_data.append(temp_df[['Tool', 'RealTime_sec', 'MaxMem_MB']])

        except Exception as e:
            print(f"ERROR: Could not read {file_name}: {e}", file=sys.stderr)

    if not all_data:
        print("ERROR: No valid data loaded. Check input directory and file headers.", file=sys.stderr)
        sys.exit(1)

    # Combine and clean
    df = pd.concat(all_data, ignore_index=True)
    df['RealTime_sec'] = pd.to_numeric(df['RealTime_sec'], errors='coerce')
    df['MaxMem_MB'] = pd.to_numeric(df['MaxMem_MB'], errors='coerce')
    df = df.dropna(subset=['RealTime_sec', 'MaxMem_MB'])

    # Statistics
    tool_order = sorted(df['Tool'].unique())
    summary = df.groupby('Tool').agg({
        'RealTime_sec': ['mean', 'sum'],
        'MaxMem_MB': ['mean', 'max']
    }).round(2)
    summary.columns = ['Avg_Time_sec', 'Total_Time_sec', 'Avg_Mem_MB', 'Peak_Mem_MB']
    summary = summary.reset_index()

    os.makedirs(args.output_dir, exist_ok=True)
    sns.set_theme(style="whitegrid")

    # Plot 1: Total runtime (excluding Centrolign_Pang to maintain readable scale)
    summary_no_pang = summary[summary['Tool'] != 'Centrolign_Pang'].copy()
    tool_order_no_pang = [t for t in tool_order if t != 'Centrolign_Pang']

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.barplot(data=summary_no_pang, x='Tool', y='Total_Time_sec',
                palette='magma', order=tool_order_no_pang, ax=ax)
    ax.set_yscale('log')
    ax.set_title('Total Runtime (Excluding Pangenome)', fontweight='bold', fontsize=14)
    ax.set_ylabel('Total RealTime [sec] (log10)')
    ax.tick_params(axis='x', rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, '1_total_runtime.png'), dpi=300)
    plt.close()

    # Plot 2: Peak memory
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.barplot(data=summary, x='Tool', y='Peak_Mem_MB',
                palette='magma', order=tool_order, ax=ax)
    ax.set_yscale('log')
    ax.set_title('Peak Memory Usage per Tool', fontweight='bold', fontsize=14)
    ax.set_ylabel('Peak MaxMem [MB] (log10)')
    ax.tick_params(axis='x', rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, '2_peak_memory.png'), dpi=300)
    plt.close()

    # Plot 3: Efficiency scatter (all observations)
    fig, ax = plt.subplots(figsize=(10, 7))
    sns.scatterplot(data=df, x='RealTime_sec', y='MaxMem_MB', hue='Tool',
                    palette='magma', hue_order=tool_order, s=100, alpha=0.8, ax=ax)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title('Efficiency: Time vs. Memory', fontweight='bold', fontsize=14)
    ax.set_xlabel('RealTime [sec] (log10)')
    ax.set_ylabel('MaxMem [MB] (log10)')
    ax.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, '3_efficiency_scatter.png'), dpi=300)
    plt.close()

    # Export summary statistics
    summary.to_csv(args.stats_out, sep='\t', index=False)

    print("-" * 80)
    print("BENCHMARKING SUMMARY")
    print("-" * 80)
    print(summary.to_string(index=False))
    print("-" * 80)
    print(f"Results saved to: {args.output_dir}")
    print(f"Stats exported to: {args.stats_out}")


if __name__ == "__main__":
    main()