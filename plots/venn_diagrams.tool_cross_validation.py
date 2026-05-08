import pandas as pd
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import argparse
import os

"""
GENOMIC VARIANT CONCORDANCE & HOMOPOLYMER ANALYSIS
--------------------------------------------------
Modified to accept command-line arguments for file paths.
"""


def is_homopolymer(row):
    """Detects if a variant is a single-base repeat expansion/contraction."""
    ref = str(row.get('ref_base', row.get('ref_seq', ''))).upper().strip()
    alt = str(row.get('alt_base', row.get('alt_seq', ''))).upper().strip()

    if not ref or not alt or ref == alt or len(ref) == len(alt):
        return False

    combined = (ref + alt).replace('-', '')
    if not combined: return False
    return len(set(combined)) == 1


def split_homo(df):
    """Splits a dataframe count into (non-homopolymer, homopolymer)."""
    if df.empty: return 0, 0
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
    parser = argparse.ArgumentParser(description="Generate a Venn diagram for genomic variant concordance.")

    # Define arguments
    parser.add_argument("-i", "--intersect", required=True,
                        help="Path to the intersection TSV file (e.g., 2gen.combined.2_of_3.tsv)")
    parser.add_argument("-m", "--master", required=True,
                        help="Path to the master variant TSV file (e.g., master.2gen.tsv)")
    parser.add_argument("-o", "--output", required=True, help="Path where the resulting Venn diagram PNG will be saved")

    args = parser.parse_args()

    # --- 1. Data Loading ---
    print(f"Loading data from: {args.intersect} and {args.master}")
    df_p = pd.read_csv(args.intersect, sep='\t')
    df_m = pd.read_csv(args.master, sep='\t')

    # --- 2. Regional Cleanup ---
    df_p = apply_strict_cleanup(df_p, "Intersection Data")
    df_m = apply_strict_cleanup(df_m, "Master Data")

    # --- 3. Position Processing ---
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

    # --- 4. Venn Subset Generation ---
    results = {
        '100': split_homo(master_s[~master_s['adj_pos'].isin(inter_pos)]),
        '010': split_homo(master_m[~master_m['adj_pos'].isin(inter_pos)]),
        '001': split_homo(master_c[~master_c['adj_pos'].isin(inter_pos)]),
        '110': split_homo(df_p[df_p['ref_pos_stretcher'].notna() & df_p['ref_pos_minimap'].notna() & df_p[
            'ref_pos_centrolign'].isna()]),
        '101': split_homo(df_p[df_p['ref_pos_stretcher'].notna() & df_p['ref_pos_minimap'].isna() & df_p[
            'ref_pos_centrolign'].notna()]),
        '011': split_homo(df_p[df_p['ref_pos_stretcher'].isna() & df_p['ref_pos_minimap'].notna() & df_p[
            'ref_pos_centrolign'].notna()]),
        '111': split_homo(df_p[df_p['ref_pos_stretcher'].notna() & df_p['ref_pos_minimap'].notna() & df_p[
            'ref_pos_centrolign'].notna()])
    }

    # --- 5. Visualization ---
    plt.figure(figsize=(10, 8))
    v = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels=('Stretcher', 'Minimap', 'Centrolign'))

    for vid, (norm, homo) in results.items():
        label = v.get_label_by_id(vid)
        if label:
            label.set_text(f"{norm}\n+{homo}")
            label.set_fontweight('bold')

    plt.title("Variant Concordance (Master vs Intersection)\nExcluded: chr22 Paternal/Hap2", pad=20)

    # Ensure output directory exists
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"Plot saved successfully to: {args.output}")


if __name__ == "__main__":
    main()