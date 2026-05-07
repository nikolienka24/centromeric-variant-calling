import pandas as pd
import os
import edlib

"""
CONSENSUS VARIANT CALLER
-----------------------
This script standardizes and compares mutation data from three bioinformatic tools:
Centrolign, Stretcher, and Minimap2. It identifies variants that are consistently 
reported across multiple methodologies to increase the reliability of calls.
"""


def normalize_sequences(ref, alt):
    """
    Standardizes sequences by converting to uppercase and removing common
    leading and trailing characters (prefix/suffix trimming).
    """
    s_ref = str(ref).upper().strip()
    s_alt = str(alt).upper().strip()

    if s_ref == s_alt:
        return s_ref, s_alt

    # Remove identical bases at the beginning
    while len(s_ref) > 0 and len(s_alt) > 0 and s_ref[0] == s_alt[0]:
        s_ref, s_alt = s_ref[1:], s_alt[1:]

    # Remove identical bases at the end
    while len(s_ref) > 0 and len(s_alt) > 0 and s_ref[-1] == s_alt[-1]:
        s_ref, s_alt = s_ref[:-1], s_alt[:-1]

    return (s_ref if s_ref != "" else "-"), (s_alt if s_alt != "" else "-")


def get_similarity(ref1, alt1, ref2, alt2):
    """
    Calculates a similarity score between two variants using sequence length
    ratios and edit distance (Levenshtein) via the edlib library.
    """
    s_ref1, s_alt1 = normalize_sequences(ref1, alt1)
    s_ref2, s_alt2 = normalize_sequences(ref2, alt2)

    len_diff1 = abs(len(s_ref1) - len(s_alt1))
    len_diff2 = abs(len(s_ref2) - len(s_alt2))

    # Avoid division by zero
    max_len_diff = max(len_diff1, len_diff2, 1)
    len_sim = 1.0 - (abs(len_diff1 - len_diff2) / max_len_diff)

    # Proceed to sequence alignment only if length similarity is sufficient
    if len_sim >= 0.50:
        seq1 = s_ref1 if len(s_ref1) > len(s_alt1) else s_alt1
        seq2 = s_ref2 if len(s_ref2) > len(s_alt2) else s_alt2

        # HW mode (Infix) for best partial overlap alignment
        res = edlib.align(seq1, seq2, mode="HW", task="distance")
        overlap_len = min(len(seq1), len(seq2))

        if overlap_len == 0:
            return 0.0

        actual_sim = 1.0 - (res["editDistance"] / overlap_len)
        return actual_sim

    return 0.0


def load_centrolign(path):
    """Parses Centrolign TSV output and renames columns for uniformity."""
    if not os.path.exists(path):
        print(f"WARNING: Centrolign file not found at: {path}")
        return pd.DataFrame(columns=['pos', 'alt_pos', 'ref', 'alt', 'tool'])
    try:
        df = pd.read_csv(path, sep='\t')
        if 'Ref_Pos' in df.columns:
            df = df.rename(columns={'Ref_Pos': 'pos', 'Mut_Pos': 'alt_pos', 'Ref_Base': 'ref', 'Alt_Base': 'alt'})
        df = df.assign(tool='centrolign')
        return df[['pos', 'alt_pos', 'ref', 'alt', 'tool']]
    except Exception as e:
        print(f"ERROR: Failed to read Centrolign data: {e}")
        return pd.DataFrame(columns=['pos', 'alt_pos', 'ref', 'alt', 'tool'])


def load_stretcher(path):
    """Parses EMBOSS Stretcher alignment output and converts to 1-based indexing."""
    if not os.path.exists(path):
        print(f"WARNING: Stretcher file not found at: {path}")
        return pd.DataFrame(columns=['pos', 'alt_pos', 'ref', 'alt', 'tool'])
    try:
        df = pd.read_csv(path, sep='\t')
        if 'start1' in df.columns:
            df = df.rename(columns={'start1': 'pos', 'start2': 'alt_pos', 'sequence1': 'ref', 'sequence2': 'alt'})
            # Adjusting coordinates from 0-based to 1-based
            df['pos'] = df['pos'] + 1
            df['alt_pos'] = df['alt_pos'] + 1
        df = df.assign(tool='stretcher')
        return df[['pos', 'alt_pos', 'ref', 'alt', 'tool']]
    except Exception as e:
        print(f"ERROR: Failed to read Stretcher data: {e}")
        return pd.DataFrame(columns=['pos', 'alt_pos', 'ref', 'alt', 'tool'])


def load_minimap2(path):
    """Parses Minimap2/Paftools VCF output."""
    if not os.path.exists(path):
        print(f"WARNING: Minimap2 file not found at: {path}")
        return pd.DataFrame(columns=['pos', 'alt_pos', 'ref', 'alt', 'tool'])
    try:
        df = pd.read_csv(path, sep='\t', comment='#', header=None,
                         names=['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'sample'])
        df['alt_pos'] = None
        df = df.assign(tool='minimap2+paftools')
        return df[['pos', 'alt_pos', 'ref', 'alt', 'tool']]
    except Exception as e:
        print(f"ERROR: Failed to read Minimap2 data: {e}")
        return pd.DataFrame(columns=['pos', 'alt_pos', 'ref', 'alt', 'tool'])


def find_matches(df_a, df_b, label_a, label_b):
    """
    Identifies variant pairs between two datasets that fall within a defined
    distance (slack) and meet a minimum sequence similarity threshold.
    """
    pairs = []
    if df_a.empty or df_b.empty:
        return pairs

    print(f"\n--- Comparing {label_a} vs {label_b} ---")

    used_b = set()
    for i, row_a in df_a.iterrows():
        s_ref, s_alt = normalize_sequences(row_a['ref'], row_a['alt'])
        variant_len = max(len(s_ref), len(s_alt))

        # Use a larger window for Structural Variants (>50bp)
        current_slack = 20000 if variant_len > 50 else 100

        candidates = df_b[(df_b['pos'] - row_a['pos']).abs() <= current_slack].copy()

        best_sim = -1
        best_match_idx = None

        for j, row_b in candidates.iterrows():
            if j in used_b:
                continue

            sim = get_similarity(row_a['ref'], row_a['alt'], row_b['ref'], row_b['alt'])
            if sim > best_sim:
                best_sim = sim
                best_match_idx = j

        # Acceptance logic
        if best_match_idx is not None and best_sim >= 0.5:
            status = "[OK]" if best_sim >= 0.75 else "[SKIP]"
            print(f"  {status} {row_a['pos']} vs {df_b.loc[best_match_idx, 'pos']} | Similarity: {best_sim:.2f}")

            # Only pairs with similarity >= 0.75 are eligible for final consensus
            if round(best_sim, 2) >= 0.75:
                pairs.append({
                    'idx_a': i, 'pos_a': row_a['pos'], 'alt_pos_a': row_a['alt_pos'],
                    'idx_b': best_match_idx, 'pos_b': df_b.loc[best_match_idx, 'pos'],
                    'alt_pos_b': df_b.loc[best_match_idx, 'alt_pos'],
                    'ref': row_a['ref'], 'alt': row_a['alt'], 'sim': best_sim
                })
                used_b.add(best_match_idx)

    return pairs


# --- Configuration and File Paths ---
home_dir = "/storage/praha5-elixir/projects/bioinf-fi/polakova/BP/__results/"
output_path = os.path.join(home_dir, "validation/3gen/chrX_maternal/")
os.makedirs(output_path, exist_ok=True)

centrolign_input = home_dir + 'centrolign_align/pedigree_3generations/chrX_maternal_3gen/mutations_3gen_final_clean.final.tsv'
minimap_input = home_dir + 'minimap2_paftools/pedigree.3generations/chrX_maternal_3gen/mutations_3gen_final.final.vcf'
stretcher_input = home_dir + "stretcher/pedigree.3generations/PAN027.chrX.maternal/intersected_generations.bedpe"

# --- Data Loading ---
df_stretcher = load_stretcher(stretcher_input)
df_minimap = load_minimap2(minimap_input)
df_centrolign = load_centrolign(centrolign_input)

print(f"Record counts:\nStretcher: {len(df_stretcher)}\nMinimap2: {len(df_minimap)}\nCentrolign: {len(df_centrolign)}")

# --- Pairwise Tool Comparisons ---
pairs_stretcher_minimap = find_matches(df_stretcher, df_minimap, "Stretcher", "Minimap")
pairs_stretcher_centrolign = find_matches(df_stretcher, df_centrolign, "Stretcher", "Centrolign")
pairs_minimap_centrolign = find_matches(df_minimap, df_centrolign, "Minimap", "Centrolign")

final_consensus_list = []
used_s_indices, used_m_indices, used_c_indices = set(), set(), set()
triplets_count = 0

# 1. Logic for Triplets (Intersects all three tools)
for pair in pairs_stretcher_minimap:
    s_idx, m_idx = pair['idx_a'], pair['idx_b']
    match_sc = next((x for x in pairs_stretcher_centrolign if x['idx_a'] == s_idx), None)
    match_mc = next((x for x in pairs_minimap_centrolign if x['idx_a'] == m_idx), None)

    if match_sc and match_mc and match_sc['idx_b'] == match_mc['idx_b']:
        c_idx = match_sc['idx_b']
        final_consensus_list.append({
            'ref_pos_stretcher': pair['pos_a'], 'alt_pos_stretcher': pair['alt_pos_a'],
            'ref_pos_minimap': pair['pos_b'],
            'ref_pos_centrolign': match_sc['pos_b'], 'alt_pos_centrolign': match_sc['alt_pos_b'],
            'ref_seq': pair['ref'], 'alt_seq': pair['alt'], 'tools_count': 3,
            'sim_S_M': round(pair['sim'], 2), 'sim_S_C': round(match_sc['sim'], 2), 'sim_M_C': round(match_mc['sim'], 2)
        })
        used_s_indices.add(s_idx)
        used_m_indices.add(m_idx)
        used_c_indices.add(c_idx)
        triplets_count += 1

# 2. Logic for Doublets (Stretcher + Minimap only)
for pair in pairs_stretcher_minimap:
    if pair['idx_a'] not in used_s_indices and pair['idx_b'] not in used_m_indices:
        final_consensus_list.append({
            'ref_pos_stretcher': pair['pos_a'], 'alt_pos_stretcher': pair['alt_pos_a'],
            'ref_pos_minimap': pair['pos_b'],
            'ref_pos_centrolign': None, 'alt_pos_centrolign': None, 'ref_seq': pair['ref'], 'alt_seq': pair['alt'],
            'tools_count': 2, 'sim_S_M': round(pair['sim'], 2), 'sim_S_C': None, 'sim_M_C': None})
        used_s_indices.add(pair['idx_a'])
        used_m_indices.add(pair['idx_b'])

# 3. Logic for Doublets (Stretcher + Centrolign only)
for pair in pairs_stretcher_centrolign:
    if pair['idx_a'] not in used_s_indices and pair['idx_b'] not in used_c_indices:
        final_consensus_list.append({
            'ref_pos_stretcher': pair['pos_a'], 'alt_pos_stretcher': pair['alt_pos_a'], 'ref_pos_minimap': None,
            'ref_pos_centrolign': pair['pos_b'], 'alt_pos_centrolign': pair['alt_pos_b'], 'ref_seq': pair['ref'],
            'alt_seq': pair['alt'], 'tools_count': 2, 'sim_S_M': None, 'sim_S_C': round(pair['sim'], 2),
            'sim_M_C': None})
        used_s_indices.add(pair['idx_a'])
        used_c_indices.add(pair['idx_b'])

# 4. Logic for Doublets (Minimap + Centrolign only)
for pair in pairs_minimap_centrolign:
    if pair['idx_a'] not in used_m_indices and pair['idx_b'] not in used_c_indices:
        final_consensus_list.append({
            'ref_pos_stretcher': None, 'alt_pos_stretcher': None, 'ref_pos_minimap': pair['pos_a'],
            'ref_pos_centrolign': pair['pos_b'], 'alt_pos_centrolign': pair['alt_pos_b'],
            'ref_seq': pair['ref'], 'alt_seq': pair['alt'], 'tools_count': 2,
            'sim_M_C': round(pair['sim'], 2), 'sim_S_M': None, 'sim_S_C': None})
        used_m_indices.add(pair['idx_a'])
        used_c_indices.add(pair['idx_b'])

# --- Final Processing and Output ---
final_df = pd.DataFrame(final_consensus_list)
if not final_df.empty:
    # Sort variants based on the earliest reference position found among tools
    final_df['sort_pos'] = final_df['ref_pos_stretcher'].fillna(final_df['ref_pos_minimap']).fillna(
        final_df['ref_pos_centrolign'])
    final_df = final_df.sort_values('sort_pos').drop(columns=['sort_pos'])

    # Save variants supported by 2 or more tools
    final_df.to_csv(os.path.join(output_path, 'final_consensus_2_of_3.tsv'), sep='\t', index=False)

    # Save variants supported by exactly all 3 tools
    df_3_of_3 = final_df[final_df['tools_count'] == 3]
    df_3_of_3.to_csv(os.path.join(output_path, 'final_consensus_3_of_3.tsv'), sep='\t', index=False)

print("\n" + "=" * 70)
print(f"VARIANTS SUPPORTED BY ALL 3 TOOLS: {triplets_count}")
print(f"TOTAL UNIQUE CONSENSUS VARIANTS (2+ tools): {len(final_df)}")
print("=" * 70)
print(f"Files saved to: {output_path}")