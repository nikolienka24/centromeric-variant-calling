import pandas as pd
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import argparse
import os
import sys

"""
GENOMIC VARIANT CONCORDANCE & HOMOPOLYMER ANALYSIS
--------------------------------------------------
Generates a 3-way Venn diagram showing concordance across three variant callers
(Stretcher, Minimap2, Centrolign). Each region is annotated with regular and
homopolymer variant counts. Chr22 paternal/haplotype2 data is excluded.
"""


def parse_args():
    """Parses and returns command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate a Venn diagram for genomic variant concordance.")
    parser.add_argument("-i", "--intersect", required=True,
                        help="Path to the intersection TSV file (e.g., 2gen.combined.2_of_3.tsv)")
    parser.add_argument("-m", "--master", required=True,
                        help="Path to the master variant TSV file (e.g., master.2gen.tsv)")
    parser.add_argument("-o", "--output", required=True,
                        help="Path where the resulting Venn diagram PNG will be saved")
    return parser.parse_args()


def is_homopolymer(row):
    """Detects if a variant is a single-base repeat expansion/contraction."""
    ref = str(row.get('ref_base', row.get('ref_seq', ''))).upper().strip()
    alt = str(row.get('alt_base', row.get('alt_seq', ''))).upper().strip()

    if not ref or not alt or ref == alt or len(ref) == len(alt):
        return False

    combined = (ref + alt).replace('-', '')
    if not combined:
        return False
    return len(set(combined)) == 1


def split_homo(df):
    """Splits a dataframe count into (non-homopolymer, homopolymer)."""
    if df.empty:
        return 0, 0
    is_h = df.apply(is_homopolymer, axis=1)
    h_count = is_h.sum()
    return len(df) - h_count, h_count


def apply_strict_cleanup(df, label):
    """Filters out chr22 paternal/haplotype2 data specifically."""
    col_name = next((c for c in df.columns if c.lower() in ['folder', 'chrom', '#chrom']), None)
    if col_name:
        vals = df[col_name].astype(str).str.lower()
        to_remove = vals.str.contains('chr22', na=False) & vals.str.contains('paternal|haplotype2', na=False)
        before = len(df)
        df = df[~to_remove].copy()
        print(f"{label}: Removed {before - len(df)} rows.")
    return df


def main():
    args = parse_args()

    # Validate input files
    for path, label in [(args.intersect, "Intersection"), (args.master, "Master")]:
        if not os.path.exists(path):
            print(f"ERROR: {label} file {path} not found.", file=sys.stderr)
            sys.exit(1)

    # Load data
    print(f"Loading data from: {args.intersect} and {args.master}")
    df_p = pd.read_csv(args.intersect, sep='\t')
    df_m = pd.read_csv(args.master, sep='\t')

    if df_p.empty or df_m.empty:
        print("ERROR: One or both input files are empty.", file=sys.stderr)
        sys.exit(1)

    # Regional cleanup
    df_p = apply_strict_cleanup(df_p, "Intersection Data")
    df_m = apply_strict_cleanup(df_m, "Master Data")

    # Position processing
    df_m['ref_pos'] = pd.to_numeric(df_m['ref_pos'], errors='coerce')
    df_m = df_m.dropna(subset=['ref_pos'])

    inter_pos = set()
    for col in ['ref_pos_stretcher', 'ref_pos_minimap', 'ref_pos_centrolign']:
        if col in df_p.columns:
            df_p[col] = pd.to_numeric(df_p[col], errors='coerce')
            inter_pos.update(df_p[df_p[col].notna()][col].astype(int).tolist())

    master_s = df_m[df_m['tool'].str.contains('stretcher', case=False)].copy()
    master_s['adj_pos'] = master_s['ref_pos'].astype(int) + 1

    master_m = df_m[df_m['tool'].str.contains('minimap', case=False)].copy()
    master_m['adj_pos'] = master_m['ref_pos'].astype(int)

    master_c = df_m[df_m['tool'].str.contains('centrolign', case=False)].copy()
    master_c['adj_pos'] = master_c['ref_pos'].astype(int)

    # Venn subset generation
    results = {
        '100': split_homo(master_s[~master_s['adj_pos'].isin(inter_pos)]),
        '010': split_homo(master_m[~master_m['adj_pos'].isin(inter_pos)]),
        '001': split_homo(master_c[~master_c['adj_pos'].isin(inter_pos)]),
        '110': split_homo(df_p[
            df_p['ref_pos_stretcher'].notna() &
            df_p['ref_pos_minimap'].notna() &
            df_p['ref_pos_centrolign'].isna()]),
        '101': split_homo(df_p[
            df_p['ref_pos_stretcher'].notna() &
            df_p['ref_pos_minimap'].isna() &
            df_p['ref_pos_centrolign'].notna()]),
        '011': split_homo(df_p[
            df_p['ref_pos_stretcher'].isna() &
            df_p['ref_pos_minimap'].notna() &
            df_p['ref_pos_centrolign'].notna()]),
        '111': split_homo(df_p[
            df_p['ref_pos_stretcher'].notna() &
            df_p['ref_pos_minimap'].notna() &
            df_p['ref_pos_centrolign'].notna()])
    }

    # Visualization
    fig, ax = plt.subplots(figsize=(10, 8))
    v = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels=('Stretcher', 'Minimap', 'Centrolign'), ax=ax)

    for vid, (norm, homo) in results.items():
        label = v.get_label_by_id(vid)
        if label:
            label.set_text(f"{norm}\n+{homo}")
            label.set_fontweight('bold')

    ax.set_title("Variant Concordance (Master vs Intersection)\nExcluded: chr22 Paternal/Hap2", pad=20)

    # Save output
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Plot saved successfully to: {args.output}")


if __name__ == "__main__":
    main()