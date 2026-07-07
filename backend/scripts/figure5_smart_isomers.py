#!/usr/bin/env python3
"""
Figure5 — Sélection intelligente des isomères NMR
====================================================
Lit Figure5_spectral_candidates.csv + les 4 fichiers d'annotation
et génère 3 nouvelles figures :

  1. Figure5_bar_normalized_smart_isomers_spectral
     Barres 100% normalisées par bin de combined_score.
     Catégories : MS unique | MS best isomère NMR | MS close isomère ±5 pts |
                  MS non-match isomère | Spectral analog | No match

  2. Figure5_bar_normalized_smart_isomers_ECFP4
     Même avec ECFP4_candidates.csv

  3. Figure5_nonmatch_isomers_spectral_similarity
     Distribution de similarité spectrale (Tanimoto moyen 1H+13C+HSQC+COSY)
     des isomères non-match vs leur best isomère NMR dans le même groupe de masse.

Logique isomères (par groupe exact_mass, MS confirmés) :
  • best_isomere   = combined_score le plus élevé
  • close_isomere  = combined_score >= best - 5
  • nonmatch       = combined_score < best - 5
"""

import os, csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Configuration ──────────────────────────────────────────────────────────────
FOLDER       = r"C:\Users\seck-s\Documents\_UNC\Projets\Article\RESULTAT"
CSV_SPECTRAL = os.path.join(FOLDER, "Figure5_spectral_candidates.csv")
CSV_ECFP4    = os.path.join(FOLDER, "Figure5_ECFP4_candidates.csv")
ANNOT_FILES  = {
    '1H':   os.path.join(FOLDER, "annotation_1H_latest.csv"),
    '13C':  os.path.join(FOLDER, "annotation_13C_latest.csv"),
    'HSQC': os.path.join(FOLDER, "annotation_HSQC_latest.csv"),
    'COSY': os.path.join(FOLDER, "annotation_COSY_latest.csv"),
}

NMR_THRESHOLD = 5.0   # écart de score (pts)
TOL = {'1H': 0.15, '13C': 2.0, 'HSQC': (3.0, 0.2), 'COSY': (0.2, 0.2)}

BIN_WIDTH  = 5
FIG_DPI    = 300
FIG_WIDTH  = 11.0
FIG_HEIGHT = 4.5

COLOR_MS_UNIQUE = '#6B3A2A'
COLOR_BEST_ISO  = '#2ECC71'
COLOR_CLOSE_ISO = '#F39C12'
COLOR_NONMATCH  = '#E74C3C'
COLOR_ANALOG    = '#007FFF'
COLOR_NOMATCH   = '#154360'

# ── Utilitaires ────────────────────────────────────────────────────────────────

def read_auto(path):
    with open(path, 'r', encoding='utf-8-sig', errors='ignore') as f:
        sample = f.read(2048)
    try:
        sep = csv.Sniffer().sniff(sample, delimiters=[',', ';', '\t']).delimiter
    except csv.Error:
        sep = ','
    return pd.read_csv(path, sep=sep, engine='python')


def parse_1d(s):
    if pd.isna(s) or not str(s).strip():
        return []
    return [float(x) for x in str(s).split('|') if x.strip()]


def parse_2d(s):
    if pd.isna(s) or not str(s).strip():
        return []
    out = []
    for pair in str(s).split('||'):
        parts = pair.split('|')
        if len(parts) == 2:
            try:
                out.append((float(parts[0]), float(parts[1])))
            except ValueError:
                pass
    return list(set(out))


def tanimoto_1d(A, B, tol):
    if not A or not B:
        return 0.0
    A_, B_ = sorted(set(A)), sorted(set(B))
    matched, Bc = 0, B_.copy()
    for a in A_:
        diffs = [abs(b - a) for b in Bc]
        if not diffs:
            break
        i = min(range(len(diffs)), key=lambda i: diffs[i])
        if diffs[i] <= tol:
            matched += 1
            Bc.pop(i)
    union = len(set(A_)) + len(set(B_)) - matched
    return matched / union if union > 0 else 0.0


def tanimoto_2d(A, B, t1, t2):
    if not A or not B:
        return 0.0
    matched, used = 0, set()
    for a in A:
        for i, b in enumerate(B):
            if i not in used and abs(b[0]-a[0]) <= t1 and abs(b[1]-a[1]) <= t2:
                matched += 1
                used.add(i)
                break
    union = len(A) + len(B) - matched
    return matched / union if union > 0 else 0.0


def mean_spectral_tanimoto(shifts_a, shifts_b):
    t1 = tanimoto_1d(shifts_a['1H'],   shifts_b['1H'],   TOL['1H'])
    t2 = tanimoto_1d(shifts_a['13C'],  shifts_b['13C'],  TOL['13C'])
    t3 = tanimoto_2d(shifts_a['HSQC'], shifts_b['HSQC'], TOL['HSQC'][0], TOL['HSQC'][1])
    t4 = tanimoto_2d(shifts_a['COSY'], shifts_b['COSY'], TOL['COSY'][0], TOL['COSY'][1])
    return (t1 + t2 + t3 + t4) / 4


# ── Chargement des shifts depuis les fichiers d'annotation ─────────────────────

def load_shifts():
    """Retourne un dict {mol_name: {'1H': [...], '13C': [...], 'HSQC': [...], 'COSY': [...]}}
    Le nom dans les annotations est Source = mol_name + '_predicted'.
    """
    shifts = {}  # {mol_name: {method: data}}

    for method, path in ANNOT_FILES.items():
        if not os.path.exists(path):
            print(f"  ⚠ Annotation manquante : {path}")
            continue
        df = read_auto(path)
        if 'Source' not in df.columns or 'ppm_molecule' not in df.columns:
            print(f"  ⚠ Colonnes manquantes dans {path}")
            continue
        for _, row in df.iterrows():
            src = str(row['Source'])
            # Enlever le suffixe _predicted
            mol = src.replace('_predicted', '') if src.endswith('_predicted') else src
            if mol not in shifts:
                shifts[mol] = {'1H': [], '13C': [], 'HSQC': [], 'COSY': []}
            if method in ('1H', '13C'):
                shifts[mol][method] = parse_1d(row.get('ppm_molecule', ''))
            else:
                shifts[mol][method] = parse_2d(row.get('ppm_molecule', ''))

    print(f"  ✓ Shifts chargés pour {len(shifts)} molécules")
    return shifts


# ── Classification des isomères ────────────────────────────────────────────────

def classify_isomers(df, shifts, threshold=5.0):
    """Ajoute 'nmr_iso_category' et 'spectral_sim_vs_best' au DataFrame."""
    df = df.copy()
    df['nmr_iso_category'] = df.get('category_spectral', df.get('category_ecfp4', 'No match'))
    df['spectral_sim_vs_best'] = np.nan

    ms_mask = df['found_mass'].str.lower() == 'yes'
    ms_df   = df[ms_mask]

    if 'exact_mass' not in df.columns:
        return df

    for mass_val, group in ms_df.groupby('exact_mass'):
        if len(group) == 1:
            continue  # MS unique — inchangé

        best_score = group['combined_score'].max()
        best_idx   = group['combined_score'].idxmax()
        best_mol   = df.at[best_idx, 'Source']
        shifts_best = shifts.get(best_mol, {'1H':[], '13C':[], 'HSQC':[], 'COSY':[]})

        for idx in group.index:
            score = group.at[idx, 'combined_score']
            gap   = best_score - score

            if idx == best_idx:
                df.at[idx, 'nmr_iso_category'] = 'MS_best_isomer'
            elif gap <= threshold:
                df.at[idx, 'nmr_iso_category'] = 'MS_close_isomer'
            else:
                df.at[idx, 'nmr_iso_category'] = 'MS_nonmatch_isomer'
                nm_mol = df.at[idx, 'Source']
                shifts_nm = shifts.get(nm_mol, {'1H':[], '13C':[], 'HSQC':[], 'COSY':[]})
                sim = mean_spectral_tanimoto(shifts_nm, shifts_best) * 100
                df.at[idx, 'spectral_sim_vs_best'] = sim

    return df


# ── Figures ────────────────────────────────────────────────────────────────────

def make_bins():
    return np.arange(20, 101, BIN_WIDTH)


def plot_normalized_bar(ax, df, bins, title, analog_label):
    bin_edges  = bins
    bin_labels = [f'{int(bins[i])}-{int(bins[i+1])}' for i in range(len(bins)-1)]

    cat_n = df['nmr_iso_category'].value_counts()
    categories = [
        ('MS_unique',          COLOR_MS_UNIQUE, f"MS unique (n={cat_n.get('MS_unique', 0)})"),
        ('MS_best_isomer',     COLOR_BEST_ISO,  f"MS best isomère NMR (n={cat_n.get('MS_best_isomer', 0)})"),
        ('MS_close_isomer',    COLOR_CLOSE_ISO, f"MS isomère proche ±{NMR_THRESHOLD:.0f} pts (n={cat_n.get('MS_close_isomer', 0)})"),
        ('MS_nonmatch_isomer', COLOR_NONMATCH,  f"MS non-match isomère (n={cat_n.get('MS_nonmatch_isomer', 0)})"),
        (analog_label,         COLOR_ANALOG,    f"Analogue {analog_label} (n={cat_n.get(analog_label, 0)})"),
        ('No match',           COLOR_NOMATCH,   f"Sans match (n={cat_n.get('No match', 0)})"),
    ]

    counts_list, colors_list, labels_list = [], [], []
    for cat, color, label in categories:
        subset = df[df['nmr_iso_category'] == cat]['combined_score']
        if len(subset) > 0:
            c, _ = np.histogram(subset, bins=bin_edges)
            counts_list.append(c.astype(float))
            colors_list.append(color)
            labels_list.append(label)

    if not counts_list:
        return

    totals = sum(counts_list)
    x      = np.arange(len(bin_labels))
    bottom = np.zeros(len(bin_labels))

    for c, color, label in zip(counts_list, colors_list, labels_list):
        norm = np.zeros_like(c)
        mask = totals > 0
        norm[mask] = (c[mask] / totals[mask]) * 100
        ax.bar(x, norm, 0.8, bottom=bottom, color=color,
               edgecolor='white', linewidth=0.5, label=label, alpha=0.87)
        bottom += norm

    for i, total in enumerate(totals):
        if total > 0:
            ax.text(i, 102, f'n={int(total)}', ha='center', va='bottom',
                    fontsize=6, fontweight='bold', color='#333333')

    ax.set_xlabel('Combined dereplication score (mean of 4 experiments)', fontsize=11)
    ax.set_ylabel('Percentage (%)', fontsize=11)
    ax.set_ylim(0, 115)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels, rotation=45, ha='right', fontsize=8)
    ax.legend(loc='upper left', fontsize=8.5, frameon=True, edgecolor='#CCCCCC')
    ax.spines['top'].set_visible(False)
    ax.tick_params(labelsize=9)
    ax.axhline(y=100, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.set_title(title, fontsize=10, pad=8)


def plot_spectral_sim(ax, df):
    """Distribution de similarité spectrale (0-100%) des isomères non-match
    vs leur best isomère NMR dans le même groupe de masse.
    """
    nm = df[df['nmr_iso_category'] == 'MS_nonmatch_isomer'].dropna(subset=['spectral_sim_vs_best'])

    if len(nm) == 0:
        ax.text(0.5, 0.5, 'Aucun isomère non-match', transform=ax.transAxes,
                ha='center', va='center', fontsize=12, color='gray')
        ax.set_title('Spectral similarity — non-match isomers vs best NMR isomer', fontsize=10, pad=8)
        return

    sim_bins   = np.arange(0, 105, 10)
    sim_labels = [f'{int(sim_bins[i])}-{int(sim_bins[i+1])}%' for i in range(len(sim_bins)-1)]
    counts, _  = np.histogram(nm['spectral_sim_vs_best'], bins=sim_bins)

    x = np.arange(len(sim_labels))
    ax.bar(x, counts, 0.8, color=COLOR_NONMATCH, edgecolor='white',
           linewidth=0.5, alpha=0.87, label=f'Non-match isomers (n={len(nm)})')

    for i, c in enumerate(counts):
        if c > 0:
            ax.text(i, c + 0.2, str(int(c)), ha='center', va='bottom',
                    fontsize=8, fontweight='bold', color='#333333')

    ax.set_xlabel('Spectral similarity vs best NMR isomer (%)', fontsize=11)
    ax.set_ylabel('Number of isomers', fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(sim_labels, rotation=45, ha='right', fontsize=8)
    ax.legend(loc='upper right', fontsize=9, frameon=True, edgecolor='#CCCCCC')
    ax.spines['top'].set_visible(False)
    ax.tick_params(labelsize=9)
    ax.set_title(
        f'Spectral similarity of non-match isomers (>{NMR_THRESHOLD:.0f} pts below best NMR isomer)',
        fontsize=10, pad=8)


def save_fig(fig, name):
    for ext in ['png', 'tiff', 'svg', 'pdf']:
        path = os.path.join(FOLDER, f'{name}.{ext}')
        fig.savefig(path, dpi=FIG_DPI, bbox_inches='tight', facecolor='white', format=ext)
        print(f"  ✓ {path}")
    plt.close(fig)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    bins = make_bins()

    print("Chargement des shifts spectraux...")
    shifts = load_shifts()

    # --- Figure 1 : Spectral ---
    print("\nLecture Figure5_spectral_candidates.csv...")
    df_sp = read_auto(CSV_SPECTRAL)
    df_sp = classify_isomers(df_sp, shifts, NMR_THRESHOLD)

    print("\nCatégories NMR (spectral) :")
    print(df_sp['nmr_iso_category'].value_counts().to_string())

    fig1, ax1 = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))
    plot_normalized_bar(ax1, df_sp, bins,
        title=f'NMR isomer selection (±{NMR_THRESHOLD:.0f} pts) — Spectral analogs',
        analog_label='Spectral_analog')
    plt.tight_layout()
    print("\nSauvegarde Figure5_bar_normalized_smart_isomers_spectral :")
    save_fig(fig1, 'Figure5_bar_normalized_smart_isomers_spectral')

    # --- Figure 2 : ECFP4 ---
    if os.path.exists(CSV_ECFP4):
        print("\nLecture Figure5_ECFP4_candidates.csv...")
        df_ec = read_auto(CSV_ECFP4)
        df_ec = classify_isomers(df_ec, shifts, NMR_THRESHOLD)

        # Remplacer la catégorie analogue : utiliser category_ecfp4 pour les non-MS
        if 'category_ecfp4' in df_ec.columns:
            non_ms = df_ec['found_mass'].str.lower() != 'yes'
            df_ec.loc[non_ms, 'nmr_iso_category'] = df_ec.loc[non_ms, 'category_ecfp4'].replace(
                {'MS_isomers': 'No match', 'MS_unique': 'No match'})

        print("\nCatégories NMR (ECFP4) :")
        print(df_ec['nmr_iso_category'].value_counts().to_string())

        fig2, ax2 = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))
        plot_normalized_bar(ax2, df_ec, bins,
            title=f'NMR isomer selection (±{NMR_THRESHOLD:.0f} pts) — ECFP4 analogs',
            analog_label='ECFP4_analog')
        plt.tight_layout()
        print("\nSauvegarde Figure5_bar_normalized_smart_isomers_ECFP4 :")
        save_fig(fig2, 'Figure5_bar_normalized_smart_isomers_ECFP4')
    else:
        print(f"\n⚠ {CSV_ECFP4} non trouvé — figure ECFP4 ignorée")

    # --- Figure 3 : Similarité spectrale des non-match ---
    print("\nStats similarité spectrale des non-match isomères :")
    nm = df_sp[df_sp['nmr_iso_category'] == 'MS_nonmatch_isomer'].dropna(subset=['spectral_sim_vs_best'])
    if len(nm) > 0:
        print(nm['spectral_sim_vs_best'].describe().to_string())
    else:
        print("  Aucun non-match isomère avec similarité calculée")

    fig3, ax3 = plt.subplots(figsize=(FIG_WIDTH * 0.7, FIG_HEIGHT))
    plot_spectral_sim(ax3, df_sp)
    plt.tight_layout()
    print("\nSauvegarde Figure5_nonmatch_isomers_spectral_similarity :")
    save_fig(fig3, 'Figure5_nonmatch_isomers_spectral_similarity')

    # --- Export CSV ---
    out_cols = [c for c in ['Source', 'score_1H', 'score_13C', 'score_HSQC', 'score_COSY',
                             'combined_score', 'found_mass', 'exact_mass',
                             'category_spectral', 'nmr_iso_category', 'spectral_sim_vs_best']
                if c in df_sp.columns]
    csv_out = os.path.join(FOLDER, 'Figure5_smart_isomers_candidates.csv')
    df_sp[out_cols].to_csv(csv_out, index=False)
    print(f"\n✓ CSV exporté : {csv_out}")
    print("\n✅ Terminé.")


if __name__ == '__main__':
    main()
