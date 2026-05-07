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


def main():
    # --- 1. Argument Parsing ---
    parser = argparse.ArgumentParser(description="Generate benchmarking plots and stats for tool comparisons.")
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing tool .tsv performance files")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save generated plots")
    parser.add_argument("-s", "--stats_out", required=True, help="Path to save the final summary TSV file")

    args = parser.parse_args()

    # Define the mapping of labels to filenames in the input directory
    file_map = {
        'Minimap2/Deepvariant': 'map+deepvariant.tsv',
        'Centrolign_Align': 'centrolign_align.tsv',
        'Minimap2/Paftools': 'minimap2+paftools.tsv',
        'Stretcher': 'stretcher_align.tsv',
        'Centrolign_Pang': 'centrolign_pangenomes.tsv'
    }

    all_data = []

    # --- 2. Data Loading and Normalization ---
    print(f"Reading performance logs from: {args.input_dir}")

    for tool_name, file_name in file_map.items():
        file_path = os.path.join(args.input_dir, file_name)

        if os.path.exists(file_path):
            try:
                # Load file using whitespace/tab separator
                temp_df = pd.read_csv(file_path, sep='\s+').dropna(how='all')

                if not temp_df.empty:
                    # Standardize time column names (Minimap2/Paftools often uses Time_sec)
                    if 'Time_sec' in temp_df.columns:
                        temp_df = temp_df.rename(columns={'Time_sec': 'RealTime_sec'})

                    temp_df['Tool'] = tool_name
                    # Keep only core metrics
                    all_data.append(temp_df[['Tool', 'RealTime_sec', 'MaxMem_MB']])
            except Exception as e:
                print(f"Error reading {file_name}: {e}")
        else:
            print(f"Warning: File not found: {file_path}")

    if not all_data:
        print("Error: No data loaded. Please check your file paths.")
        sys.exit(1)

    # Combine all dataframes and clean types
    df = pd.concat(all_data, ignore_index=True)
    df['RealTime_sec'] = pd.to_numeric(df['RealTime_sec'], errors='coerce')
    df['MaxMem_MB'] = pd.to_numeric(df['MaxMem_MB'], errors='coerce')
    df = df.dropna(subset=['RealTime_sec', 'MaxMem_MB'])

    # --- 3. Statistical Calculations ---
    tool_order = sorted(df['Tool'].unique())
    summary = df.groupby('Tool').agg({
        'RealTime_sec': ['mean', 'sum'],
        'MaxMem_MB': ['mean', 'max']
    }).round(2)

    summary.columns = ['Avg_Time_sec', 'Total_Time_sec', 'Avg_Mem_MB', 'Peak_Mem_MB']
    summary = summary.reset_index()

    # --- 4. Plotting Setup ---
    os.makedirs(args.output_dir, exist_ok=True)
    sns.set_theme(style="whitegrid")

    # PLOT 1: Total Runtime (Sum)
    # Exclude Centrolign_Pang for runtime plotting to maintain a readable scale
    plt.figure(figsize=(8, 6))
    summary_no_pang = summary[summary['Tool'] != 'Centrolign_Pang'].copy()
    tool_order_no_pang = [t for t in tool_order if t != 'Centrolign_Pang']

    ax1 = sns.barplot(data=summary_no_pang, x='Tool', y='Total_Time_sec',
                      palette='magma', order=tool_order_no_pang)
    ax1.set_yscale('log')
    plt.title('Total Runtime (Excluding Pangenome)', fontweight='bold', fontsize=14)
    plt.ylabel('Total RealTime [sec] (log10)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, '1_total_runtime.png'), dpi=300)
    plt.close()

    # PLOT 2: Peak Memory (Max)
    plt.figure(figsize=(8, 6))
    ax2 = sns.barplot(data=summary, x='Tool', y='Peak_Mem_MB',
                      palette='magma', order=tool_order)
    ax2.set_yscale('log')
    plt.title('Peak Memory Usage per Tool', fontweight='bold', fontsize=14)
    plt.ylabel('Peak MaxMem [MB] (log10)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, '2_peak_memory.png'), dpi=300)
    plt.close()

    # PLOT 3: Efficiency Scatter (Time vs. Memory)
    plt.figure(figsize=(10, 7))
    ax3 = sns.scatterplot(data=df, x='RealTime_sec', y='MaxMem_MB', hue='Tool',
                          palette='magma', hue_order=tool_order, s=100, alpha=0.8)
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    plt.title('Efficiency: Time vs. Memory', fontweight='bold', fontsize=14)
    plt.xlabel('RealTime [sec] (log10)')
    plt.ylabel('MaxMem [MB] (log10)')
    plt.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, '3_efficiency_scatter.png'), dpi=300)
    plt.close()

    # --- 5. Export and Logging ---
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