import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys

"""
BENCHMARKING VISUALIZATION DASHBOARD
------------------------------------
This script aggregates performance metrics (Runtime and Memory) from various
bioinformatics tools used in the pedigree analysis pipeline.

It performs the following:
1. Loads performance TSV files for different tools.
2. Unifies column names (e.g., Time_sec vs RealTime_sec).
3. Calculates descriptive statistics (means, sums, and peaks).
4. Generates a multi-panel visualization (Dashboard) with log-scale axes.
5. Exports a summary statistics table.
"""


def parse_args():
    """Parses and returns command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate benchmarking dashboard from tool performance logs.")
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Directory containing tool-specific .tsv performance files")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Directory where plots and summary stats will be saved")
    return parser.parse_args()


def main():
    args = parse_args()

    # Validate input directory
    if not os.path.isdir(args.input_dir):
        print(f"ERROR: Input directory {args.input_dir} not found.", file=sys.stderr)
        sys.exit(1)

    # Mapping tool names to their expected log filenames
    file_map = {
        'Centrolign_Align': 'centrolign_align.tsv',
        'Centrolign_Pang': 'centrolign_pangenomes.tsv',
        'Minimap2/Deepvariant': 'map+deepvariant.tsv',
        'Minimap2/Paftools': 'minimap2+paftools.tsv',
        'Stretcher': 'stretcher_align.tsv'
    }

    all_data = []

    # Load and unify data
    print(f"Loading data from: {args.input_dir}")
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

            # Unify time column naming conventions
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

    # Dashboard visualization
    print(f"Generating plots in: {args.output_dir}")
    os.makedirs(args.output_dir, exist_ok=True)

    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 3, figsize=(22, 7))

    # Plot 1: Total runtime (excluding Centrolign_Pang to preserve scale)
    summary_plot = summary[summary['Tool'] != 'Centrolign_Pang']
    order_plot = [t for t in tool_order if t != 'Centrolign_Pang']

    sns.barplot(ax=axes[0], data=summary_plot, x='Tool', y='Total_Time_sec',
                palette='magma', order=order_plot)
    axes[0].set_yscale('log')
    axes[0].set_title('Total Runtime (Excl. Pangenome)', fontweight='bold')
    axes[0].set_ylabel('Total Time [sec] (log10)')
    axes[0].tick_params(axis='x', rotation=45)

    # Plot 2: Peak memory
    sns.barplot(ax=axes[1], data=summary, x='Tool', y='Peak_Mem_MB',
                palette='magma', order=tool_order)
    axes[1].set_yscale('log')
    axes[1].set_title('Peak Memory Usage', fontweight='bold')
    axes[1].set_ylabel('Max Memory [MB] (log10)')
    axes[1].tick_params(axis='x', rotation=45)

    # Plot 3: Efficiency scatter (all observations)
    sns.scatterplot(ax=axes[2], data=df, x='RealTime_sec', y='MaxMem_MB', hue='Tool',
                    palette='magma', hue_order=tool_order, s=100, alpha=0.6)
    axes[2].set_xscale('log')
    axes[2].set_yscale('log')
    axes[2].set_title('Efficiency: Time vs. Memory', fontweight='bold')
    axes[2].set_xlabel('RealTime [sec] (log10)')
    axes[2].set_ylabel('MaxMem [MB] (log10)')
    axes[2].legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()

    out_img = os.path.join(args.output_dir, 'benchmarking_dashboard.png')
    plt.savefig(out_img, dpi=300, bbox_inches='tight')
    plt.close()

    # Export summary statistics
    stats_file = os.path.join(args.output_dir, 'benchmarking_summary_stats.tsv')
    summary.to_csv(stats_file, sep='\t', index=False)

    print("-" * 50)
    print(f"Success! Dashboard saved to: {out_img}")
    print(f"Summary table saved to: {stats_file}")
    print("-" * 50)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()