#!/usr/bin/env python3
"""
NMRSCOUT — Pipeline complet unifié
====================================
Script 1 : Construction du réseau moléculaire spectral
Script 2 : Validation MS (recherche d'adducts HRMS)
Script 3 : Génération Figure 5 (distribution 3 couches)

CORRECTION : intersection_coverage pour 1H/13C, coverage_ransy pour HSQC/COSY

Entrées :
  - ID_SMILES_ALL_ASP.csv          (SMILES des candidats)
  - 1H_With_Ransy_and_ECFP4.csv    (annotation 1H)
  - 13_With_Ransy_and_ECFP4.csv    (annotation 13C)
  - HSQC_annotation_2D_v12_mod_with_ECFP4.csv  (annotation HSQC)
  - COSY_annotation_2D_v12_with_ECFP4.csv      (annotation COSY)
  - exact_masses.csv                (masses exactes théoriques)
  - Pos_MZMINE.csv                  (spectres MS mode positif)
  - Neg_MZMINE.csv                  (spectres MS mode négatif)

Sorties :
  - clusters_nodes_structural_first.csv
  - spectral_network_structural_first.graphml
  - metadata_with_found_mass_with_adducts.csv
  - spectral_network_annot_with_adducts.graphml
  - adducts_detection_report.txt
  - Figure5_score_distribution.png / .tiff / .svg / .pdf
"""

import os
import csv
import time
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')  # pour serveurs sans écran
import matplotlib.pyplot as plt
from datetime import datetime
from tqdm import tqdm
from itertools import combinations
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors
from rdkit.DataStructs import ExplicitBitVect, TanimotoSimilarity

# ================================================================
# CONFIGURATION — À MODIFIER SELON VOTRE INSTALLATION
# ================================================================
FOLDER = r"C:\Users\seck-s\Documents\_UNC\Projets\Article\RESULTAT"

# Fichiers d'entrée
EXACT_MASSES_FILE = 'exact_masses.csv'
POS_MS_FILE       = 'Pos_MZMINE.csv'
NEG_MS_FILE       = 'Neg_MZMINE.csv'

ANNOTATION_FILES = [
    ('annotation_1H_latest.csv',    '1H'),
    ('annotation_13C_latest.csv',   '13C'),
    ('annotation_HSQC_latest.csv',  'HSQC'),
    ('annotation_COSY_latest.csv',  'COSY'),
]

# Colonne de score par expérience (CORRECTION)
SCORE_COL = {
    '1H':   'intersection_coverage',
    '13C':  'intersection_coverage',
    'HSQC': 'coverage_ransy',
    'COSY': 'coverage_ransy',
}

# Seuils
TAN_THRESH   = 50.0   # Tanimoto pour clustering structural
HIT_THRESH   = {'1H': 70.0, '13C': 60.0, 'HSQC': 60.0, 'COSY': 50.0}
SPEC_THRESH  = {'1H': 70.0, '13C': 60.0, 'HSQC': 60.0, 'COSY': 50.0}
TOL          = {'1H': (0.1,), '13C': (2.0,), 'HSQC': (3.0, 0.2), 'COSY': (0.2, 0.2)}

# MS validation
PROTON_MASS = 1.007276466621
TOL_PPM     = 5
ADDUCTS_POS = [
    (PROTON_MASS,  '[M+H]+'),
    (21.9820,      '[M+Na-H]+'),
    (37.9555,      '[M+K-H]+'),
]
ADDUCTS_NEG = [
    (-PROTON_MASS, '[M-H]-'),
    (34.9694,      '[M+Cl-H]-'),
    (44.9982,      '[M+FA-H]-'),
    (20.9747,      '[M+Na-2H]-'),
]

# Figure 5 - Tanimoto SPECTRAL
TANIMOTO_SPECTRAL = 0.7    # seuil pour les analogues spectraux (70%)
TOL_PPM_SPECTRAL  = 10     # tolérance MS pour figures spectrales
BIN_WIDTH         = 5

# Tolérances pour Tanimoto spectral (Figure 5)
TOL_FIG5 = {
    '1H':   (0.15,),           # 0.15 ppm
    '13C':  (2.0,),            # 2 ppm
    'HSQC': (3.0, 0.2),        # F1: 3 ppm (13C), F2: 0.2 ppm (1H)
    'COSY': (0.2, 0.2),        # F1: 0.2 ppm, F2: 0.2 ppm
}

# Figure 5 - Tanimoto ECFP4
TANIMOTO_ECFP4    = 0.7    # seuil pour les analogues ECFP4 (70%)
TOL_PPM_ECFP4     = 10     # tolérance MS pour figures ECFP4
FIG_DPI           = 300
FIG_WIDTH         = 7.0    # double colonne Anal. Chem.
FIG_HEIGHT        = 4.0
COLOR_ALL         = '#13456B'
COLOR_ANALOG      = '#007FFF'  # Bleu azur pour analogues (MS non confirmed + Tanimoto)
COLOR_MS          = '#6B3A2A'
COLOR_MS_AMBIG    = '#C0392B'  # Rouge pour MS ambigus (masse partagée)
COLOR_UNMATCHED   = '#154360'  # Bleu foncé pour sans match (MS non confirmed)


# ================================================================
# FONCTIONS UTILITAIRES
# ================================================================

def log_progress(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"\n[{ts}] {msg}")

def read_csv_auto(path, **kwargs):
    with open(path, 'r', encoding=kwargs.get('encoding', 'utf-8'), errors='ignore') as f:
        sample = f.read(2048)
    try:
        sep = csv.Sniffer().sniff(sample, delimiters=[',', ';', '\t']).delimiter
    except csv.Error:
        sep = ','
    return pd.read_csv(path, sep=sep, engine='python', **kwargs)

def parse_1D_signals(s):
    if pd.isna(s) or not s:
        return []
    return [float(x) for x in str(s).split('|') if x]

def parse_2D_signals(s):
    if pd.isna(s) or not s:
        return []
    out = []
    for pair in str(s).split('||'):
        parts = pair.split('|')
        if len(parts) == 2:
            try:
                out.append((float(parts[0]), float(parts[1])))
            except:
                pass
    return list(set(out))

def match_one2one_1D(A, B, tol):
    A_, B_ = sorted(set(A)), sorted(set(B))
    matched, Bc = [], B_.copy()
    for a in A_:
        if not Bc:
            break
        diffs = [abs(b - a) for b in Bc]
        i = min(range(len(diffs)), key=lambda i: diffs[i])
        if diffs[i] <= tol:
            matched.append(Bc.pop(i))
    return matched

def match_one2one_2D(A, B, tol1, tol2):
    A_, B_ = list({tuple(x) for x in A}), list({tuple(x) for x in B})
    matched, used = [], set()
    for a in A_:
        cands = [(i, b) for i, b in enumerate(B_)
                 if i not in used
                 and abs(b[0] - a[0]) <= tol1
                 and abs(b[1] - a[1]) <= tol2]
        if not cands:
            continue
        bi, bb = min(cands, key=lambda ib: max(abs(ib[1][0] - a[0]), abs(ib[1][1] - a[1])))
        matched.append(bb)
        used.add(bi)
    return matched

def spectral_tanimoto_1D(A, B, tol):
    """
    Calcule le Tanimoto spectral 1D (1H ou 13C).
    Tanimoto = |matched| / (|A| + |B| - |matched|)
    Matching one-to-one : chaque shift ne matche qu'une fois.
    """
    if not A or not B:
        return 0.0
    A_set = sorted(set(A))
    B_set = sorted(set(B))
    matched = match_one2one_1D(A_set, B_set, tol)
    n_matched = len(matched)
    union = len(A_set) + len(B_set) - n_matched
    if union == 0:
        return 0.0
    return n_matched / union

def spectral_tanimoto_2D(A, B, tol1, tol2):
    """
    Calcule le Tanimoto spectral 2D (HSQC ou COSY).
    Tanimoto = |matched| / (|A| + |B| - |matched|)
    Matching one-to-one : chaque pic ne matche qu'une fois.
    """
    if not A or not B:
        return 0.0
    A_set = list({tuple(x) for x in A})
    B_set = list({tuple(x) for x in B})
    matched = match_one2one_2D(A_set, B_set, tol1, tol2)
    n_matched = len(matched)
    union = len(A_set) + len(B_set) - n_matched
    if union == 0:
        return 0.0
    return n_matched / union

def tanimoto_from_bitstr(bs1, bs2):
    n = len(bs1)
    v1, v2 = ExplicitBitVect(n), ExplicitBitVect(n)
    for i, ch in enumerate(bs1):
        if ch == '1':
            v1.SetBit(i)
    for i, ch in enumerate(bs2):
        if ch == '1':
            v2.SetBit(i)
    return DataStructs.TanimotoSimilarity(v1, v2) * 100

def bitstr_to_vec(bs):
    if pd.isna(bs) or not bs:
        return None
    bs = str(bs)
    n = len(bs)
    v = ExplicitBitVect(n)
    for i, ch in enumerate(bs):
        if ch == '1':
            v.SetBit(i)
    return v

def spectral_similarity_1D(A, B, tol):
    sA = len(match_one2one_1D(A, B, tol)) / len(A) * 100 if A else 0
    sB = len(match_one2one_1D(B, A, tol)) / len(B) * 100 if B else 0
    return sA, sB, (sA + sB) / 2

def spectral_similarity_2D(A, B, tol1, tol2):
    sA = len(match_one2one_2D(A, B, tol1, tol2)) / len(A) * 100 if A else 0
    sB = len(match_one2one_2D(B, A, tol1, tol2)) / len(B) * 100 if B else 0
    return sA, sB, (sA + sB) / 2


# ================================================================
# PARTIE 1 : CONSTRUCTION DU RÉSEAU SPECTRAL
# ================================================================

def build_spectral_network():
    log_progress("=" * 60)
    log_progress("PARTIE 1 : CONSTRUCTION DU RÉSEAU SPECTRAL")
    log_progress("=" * 60)

    # 1.1) Charger SMILES + masses depuis exact_masses.csv
    # Colonnes attendues: SMILES, Source, exact_mass
    log_progress("📂 Chargement SMILES + masses depuis exact_masses.csv")
    smiles_df = read_csv_auto(os.path.join(FOLDER, EXACT_MASSES_FILE), on_bad_lines='skip')

    # Colonnes connues
    smiles_df['MolName'] = smiles_df['Source']
    smiles_df['Mass'] = smiles_df['exact_mass']

    log_progress(f"  ✓ {len(smiles_df)} molécules chargées")

    # 1.2) Charger annotations
    log_progress("📊 Chargement annotations")
    method_presence = {}
    coverage_map = {lbl: {} for _, lbl in ANNOTATION_FILES}
    frames = {}

    for fname, label in ANNOTATION_FILES:
        log_progress(f"  → {fname}")
        # Forcer sep=';' pour les fichiers d'annotation (évite les problèmes de détection)
        df = pd.read_csv(os.path.join(FOLDER, fname), sep=';', on_bad_lines='skip')
        log_progress(f"     📊 Lignes chargées: {len(df)} (attendu: 3250)")
        # Retirer le suffixe _predicted pour matcher avec exact_masses.csv
        df['MolName'] = df['Source'].str.replace('_predicted', '', regex=False)

        # Merger avec masses
        df = df.merge(smiles_df[['MolName', 'Mass']], on='MolName', how='left')

        # Colonne de score correcte
        col = SCORE_COL[label]
        if col not in df.columns:
            log_progress(f"  ⚠ '{col}' non trouvée pour {label}, fallback 'coverage_ransy'")
            col = 'coverage_ransy' if 'coverage_ransy' in df.columns else None

        seen = set()
        for _, row in df.iterrows():
            m = row['MolName']
            if m not in seen:
                coverage_map[label][m] = (row.get(col, 0) or 0) if col else 0
                seen.add(m)
                method_presence.setdefault(m, set()).add(label)
        frames[label] = df[['ppm_molecule', 'ECFP4', 'MolName', 'Mass']].to_dict('records')
        log_progress(f"     ✓ {len(seen)} molécules (score via '{col}')")

    # 1.3) Clustering structural exhaustif
    log_progress("🏗️ Clustering structural (exhaustif)")
    ecfp_map = {r['MolName']: r['ECFP4'] for r in frames['1H']}
    structural_edges = []
    mols_list = sorted(ecfp_map)
    for i in tqdm(range(len(mols_list)), desc="Tanimoto struct"):
        for j in range(i + 1, len(mols_list)):
            m1, m2 = mols_list[i], mols_list[j]
            t = tanimoto_from_bitstr(ecfp_map[m1], ecfp_map[m2])
            if t >= TAN_THRESH:
                structural_edges.append((m1, m2, t))

    # 1.4) Fusion des doublons exacts - DÉSACTIVÉE
    # Les 3250 molécules originales sont conservées
    log_progress("🔀 Fusion des doublons: DÉSACTIVÉE (toutes les molécules conservées)")

    # 1.5) Construction des clusters
    Gs = nx.Graph()
    Gs.add_edges_from([(u, v) for u, v, _ in structural_edges])
    for m in ecfp_map:
        if m not in Gs:
            Gs.add_node(m)
    struct_map = {
        m: cid + 1
        for cid, comp in enumerate(nx.connected_components(Gs))
        for m in comp
    }
    log_progress(f"   → {len(set(struct_map.values()))} clusters formés")

    # 1.6) Similarité spectrale intra-cluster
    log_progress("🔬 Calcul similarités spectrales intra‑cluster")
    all_specs = []
    for label in frames:
        tol1 = TOL[label][0]
        tol2 = TOL[label][1] if len(TOL[label]) > 1 else None
        clusters = {}
        for r in frames[label]:
            cid = struct_map.get(r['MolName'], 0)
            if cid > 0:
                clusters.setdefault(cid, []).append(r)
        for cid, members in clusters.items():
            if len(members) < 2:
                continue
            for r1, r2 in tqdm(combinations(members, 2),
                               total=len(members) * (len(members) - 1) // 2,
                               desc=f"Spectral {label} clu{cid}"):
                if r1['ECFP4'] == r2['ECFP4'] and r1['Mass'] == r2['Mass']:
                    continue
                if label in ('1H', '13C'):
                    A = parse_1D_signals(r1['ppm_molecule'])
                    B = parse_1D_signals(r2['ppm_molecule'])
                    sA, sB, sM = spectral_similarity_1D(A, B, tol1)
                    nb = len(match_one2one_1D(A, B, tol1))
                else:
                    A = parse_2D_signals(r1['ppm_molecule'])
                    B = parse_2D_signals(r2['ppm_molecule'])
                    sA, sB, sM = spectral_similarity_2D(A, B, tol1, tol2)
                    nb = len(match_one2one_2D(A, B, tol1, tol2))
                if sM < SPEC_THRESH[label]:
                    continue
                t = tanimoto_from_bitstr(r1['ECFP4'], r2['ECFP4'])
                all_specs.append({
                    'Nstar': r1['MolName'], 'Nend': r2['MolName'],
                    'interaction': label, 'Score': round(sM, 2),
                    'Score_A_to_B': round(sA, 2), 'Score_B_to_A': round(sB, 2),
                    'Nb': nb, 'ECFP4_star': r1['ECFP4'], 'ECFP4_end': r2['ECFP4'],
                    'Mass_star': r1['Mass'], 'Mass_end': r2['Mass'], 'tanimoto': t
                })
    df_specs = pd.DataFrame(all_specs)
    log_progress(f"   ✓ {len(df_specs)} liens spectraux")

    # 1.7) Construction du MultiGraph
    log_progress("🌐 Construction du graphe")
    log_progress(f"  🔍 Debug: method_presence a {len(method_presence)} molécules")
    Gf = nx.MultiGraph()
    mols = sorted(method_presence)
    for m in tqdm(mols, desc="Nœuds"):
        s1 = coverage_map['1H'].get(m, 0)
        s2 = coverage_map['13C'].get(m, 0)
        s3 = coverage_map['COSY'].get(m, 0)
        s4 = coverage_map['HSQC'].get(m, 0)
        Gf.add_node(m,
                     struct_cluster=struct_map.get(m, 0),
                     hit_1H=int(s1 >= HIT_THRESH['1H']),
                     hit_13C=int(s2 >= HIT_THRESH['13C']),
                     hit_COSY=int(s3 >= HIT_THRESH['COSY']),
                     hit_HSQC=int(s4 >= HIT_THRESH['HSQC']),
                     score_1H=s1, score_13C=s2, score_COSY=s3, score_HSQC=s4)

    nx.set_node_attributes(Gf, {m: coverage_map['13C'].get(m, 0) for m in Gf.nodes()}, name='score_13C')
    log_progress(f"  🔍 Debug: Gf a {len(Gf.nodes)} nœuds")

    for u, v, t in tqdm(structural_edges, desc="Liens structuraux"):
        Gf.add_edge(u, v, interaction='structural', tanimoto=t)
    for _, r in tqdm(df_specs.iterrows(), total=len(df_specs), desc="Liens spectraux"):
        Gf.add_edge(r['Nstar'], r['Nend'],
                     interaction=r['interaction'], spec_score=r['Score'],
                     spec_score_A_to_B=r['Score_A_to_B'], spec_score_B_to_A=r['Score_B_to_A'],
                     tanimoto=r['tanimoto'], Nb=int(r['Nb']))

    # 1.8) Export
    log_progress("💾 Export réseau")
    nodes_df = pd.DataFrame({
        'Molecule': mols,
        'struct_cluster': [Gf.nodes[m]['struct_cluster'] for m in mols],
        'hit_1H': [Gf.nodes[m]['hit_1H'] for m in mols],
        'hit_13C': [Gf.nodes[m]['hit_13C'] for m in mols],
        'hit_HSQC': [Gf.nodes[m]['hit_HSQC'] for m in mols],
        'hit_COSY': [Gf.nodes[m]['hit_COSY'] for m in mols],
        'score_1H': [Gf.nodes[m]['score_1H'] for m in mols],
        'score_13C': [Gf.nodes[m]['score_13C'] for m in mols],
        'score_HSQC': [Gf.nodes[m]['score_HSQC'] for m in mols],
        'score_COSY': [Gf.nodes[m]['score_COSY'] for m in mols],
        'ECFP4': [ecfp_map.get(m, '') for m in mols],
        'methods_detected': [
            ','.join(sorted([lbl for lbl in HIT_THRESH if Gf.nodes[m][f'hit_{lbl}']]))
            for m in mols
        ]
    })
    nodes_path = os.path.join(FOLDER, 'clusters_nodes_structural_first.csv')
    nodes_df.to_csv(nodes_path, index=False)

    graphml_path = os.path.join(FOLDER, 'spectral_network_structural_first.graphml')
    nx.write_graphml(Gf, graphml_path, prettyprint=False)
    with open(graphml_path, 'r', encoding='utf-8') as f:
        content = f.read()
    with open(graphml_path, 'w', encoding='utf-8') as f:
        f.write(content.replace('attr.type="long"', 'attr.type="int"'))

    # Vérification des scores
    log_progress("📊 Vérification des scores :")
    for label in ['1H', '13C', 'HSQC', 'COSY']:
        vals = [coverage_map[label].get(m, 0) for m in mols]
        non_zero = [v for v in vals if v > 0]
        log_progress(f"  {label} ({SCORE_COL[label]}) : {len(non_zero)} non-zéro / {len(vals)}")
        if non_zero:
            log_progress(f"    Min={min(non_zero):.1f}  Max={max(non_zero):.1f}  Moy={sum(non_zero)/len(non_zero):.1f}")
        else:
            log_progress(f"    ⚠ TOUS LES SCORES SONT À 0 — vérifiez le nom de la colonne !")

    log_progress("✅ Partie 1 terminée")
    return Gf, nodes_df, coverage_map, mols, frames


# ================================================================
# PARTIE 2 : VALIDATION MS (RECHERCHE D'ADDUCTS)
# ================================================================

def run_ms_validation(Gf, nodes_df):
    log_progress("=" * 60)
    log_progress("PARTIE 2 : VALIDATION MS (ADDUCTS)")
    log_progress("=" * 60)

    # 2.1) Chargement des données MS
    meta_df = pd.read_csv(os.path.join(FOLDER, EXACT_MASSES_FILE))
    pos_df = pd.read_csv(os.path.join(FOLDER, POS_MS_FILE))
    neg_df = pd.read_csv(os.path.join(FOLDER, NEG_MS_FILE))

    for df, name in [(pos_df, 'positif'), (neg_df, 'negatif')]:
        if 'Precursor mz' not in df.columns:
            raise ValueError(f"Le tableau {name} doit contenir la colonne 'Precursor mz'")

    log_progress(f"  Masses théoriques : {len(meta_df)}")
    log_progress(f"  Pics MS+ : {len(pos_df)} | MS- : {len(neg_df)}")

    # 2.2) Fonctions de recherche d'adducts
    def find_ppm_difference(mass, spec_df, delta):
        target = mass + delta
        tol = TOL_PPM * target / 1e6
        window = spec_df['Precursor mz'].between(target - tol, target + tol)
        if not window.any():
            return None, None
        candidates = spec_df.loc[window, 'Precursor mz']
        ppm_diffs = (candidates - target) / target * 1e6
        best_idx = ppm_diffs.abs().idxmin()
        return float(ppm_diffs.loc[best_idx]), None

    def find_best_adduct(mass, spec_df, adduct_list):
        best_ppm, best_adduct, best_abs, best_mz = None, None, float('inf'), None
        for delta, name in adduct_list:
            ppm_diff, _ = find_ppm_difference(mass, spec_df, delta)
            if ppm_diff is not None and abs(ppm_diff) < best_abs:
                best_abs = abs(ppm_diff)
                best_ppm = ppm_diff
                best_adduct = name
                best_mz = mass + delta  # m/z théorique de l'adduct
        return best_ppm, best_adduct, best_mz

    # 2.3) Recherche
    log_progress("🔍 Recherche adducts mode positif...")
    pos_results = meta_df['exact_mass'].apply(lambda m: find_best_adduct(m, pos_df, ADDUCTS_POS))
    meta_df['diff_ppm_pos'] = pos_results.apply(lambda x: x[0])
    meta_df['adduct_pos'] = pos_results.apply(lambda x: x[1])
    meta_df['matched_mz_pos'] = pos_results.apply(lambda x: x[2])

    log_progress("🔍 Recherche adducts mode négatif...")
    neg_results = meta_df['exact_mass'].apply(lambda m: find_best_adduct(m, neg_df, ADDUCTS_NEG))
    meta_df['diff_ppm_neg'] = neg_results.apply(lambda x: x[0])
    meta_df['adduct_neg'] = neg_results.apply(lambda x: x[1])
    meta_df['matched_mz_neg'] = neg_results.apply(lambda x: x[2])

    # 2.4) Flags
    meta_df['found_pos'] = meta_df['diff_ppm_pos'].notna().map({True: 'yes', False: 'no'})
    meta_df['found_neg'] = meta_df['diff_ppm_neg'].notna().map({True: 'yes', False: 'no'})
    meta_df['found_mass'] = meta_df.apply(
        lambda r: 'yes' if (r['found_pos'] == 'yes' or r['found_neg'] == 'yes') else 'no', axis=1
    )

    # 2.5) Statistiques
    total_pos = len(meta_df[meta_df['found_pos'] == 'yes'])
    total_neg = len(meta_df[meta_df['found_neg'] == 'yes'])
    total_mass = len(meta_df[meta_df['found_mass'] == 'yes'])

    log_progress(f"\n📊 STATISTIQUES ADDUCTS :")
    log_progress(f"  Positif : {total_pos} détections")
    adduct_pos_counts = meta_df[meta_df['found_pos'] == 'yes']['adduct_pos'].value_counts()
    for adduct, count in adduct_pos_counts.items():
        log_progress(f"    {adduct} : {count}")

    log_progress(f"  Négatif : {total_neg} détections")
    adduct_neg_counts = meta_df[meta_df['found_neg'] == 'yes']['adduct_neg'].value_counts()
    for adduct, count in adduct_neg_counts.items():
        log_progress(f"    {adduct} : {count}")

    log_progress(f"  Total (pos OU neg) : {total_mass} / {len(meta_df)} ({total_mass/len(meta_df)*100:.1f}%)")

    # 2.6) Injection dans le graphe
    log_progress("🔗 Injection dans le graphe...")
    for _, row in meta_df.iterrows():
        node = row['Source']
        if node in Gf.nodes:
            Gf.nodes[node]['found_pos'] = row['found_pos']
            Gf.nodes[node]['diff_ppm_pos'] = row['diff_ppm_pos'] if pd.notna(row['diff_ppm_pos']) else ''
            Gf.nodes[node]['adduct_pos'] = row['adduct_pos'] if pd.notna(row['adduct_pos']) else ''
            Gf.nodes[node]['found_neg'] = row['found_neg']
            Gf.nodes[node]['diff_ppm_neg'] = row['diff_ppm_neg'] if pd.notna(row['diff_ppm_neg']) else ''
            Gf.nodes[node]['adduct_neg'] = row['adduct_neg'] if pd.notna(row['adduct_neg']) else ''
            Gf.nodes[node]['found_mass'] = row['found_mass']

    # 2.7) Merge avec les nœuds
    # Debug: vérifier le matching
    graph_nodes_sample = list(Gf.nodes)[:3]
    meta_source_sample = meta_df['Source'].head(3).tolist()
    log_progress(f"  🔍 Debug: Gf.nodes sample = {graph_nodes_sample}")
    log_progress(f"  🔍 Debug: meta_df['Source'] sample = {meta_source_sample}")

    # Chercher les différences
    graph_nodes_set = set(Gf.nodes)
    meta_sources_set = set(meta_df['Source'])
    missing_from_graph = meta_sources_set - graph_nodes_set
    log_progress(f"  🔍 Debug: {len(missing_from_graph)} Sources de exact_masses.csv absentes du graphe")
    if len(missing_from_graph) > 0 and len(missing_from_graph) <= 5:
        log_progress(f"     Exemples: {list(missing_from_graph)[:5]}")

    meta_filtered = meta_df[meta_df['Source'].isin(Gf.nodes)].copy()
    n_missing = len(meta_df) - len(meta_filtered)
    log_progress(f"  🔍 Debug: meta_filtered a {len(meta_filtered)} lignes (sur {len(meta_df)})")
    if n_missing > 0:
        log_progress(f"  ⚠️ {n_missing} molécules dans exact_masses.csv non présentes dans les annotations")

    clusters = {n: cid + 1 for cid, comp in enumerate(nx.connected_components(Gf)) for n in comp}
    meta_filtered['cluster'] = meta_filtered['Source'].map(clusters)

    def get_hits_str(node):
        hits = []
        for k in ['hit_1H', 'hit_13C', 'hit_COSY', 'hit_HSQC']:
            if Gf.nodes[node].get(k, 0) == 1:
                hits.append(k)
        return ','.join(hits)

    meta_filtered['hits'] = meta_filtered['Source'].map(get_hits_str)

    log_progress(f"  🔍 Debug: nodes_df['Molecule'] sample = {nodes_df['Molecule'].head(3).tolist()}")

    meta_merged = pd.merge(meta_filtered, nodes_df,
                           left_on='Source', right_on='Molecule', how='left')
    log_progress(f"  🔍 Debug: meta_merged a {len(meta_merged)} lignes après merge")
    if 'Molecule' in meta_merged.columns:
        meta_merged = meta_merged.drop(columns=['Molecule'])

    # 2.8) Export
    out_csv = os.path.join(FOLDER, 'metadata_with_found_mass_with_adducts.csv')
    meta_merged.to_csv(out_csv, index=False)
    log_progress(f"  ✓ CSV : {out_csv}")

    graph_out = os.path.join(FOLDER, 'spectral_network_annot_with_adducts.graphml')
    for n in Gf.nodes:
        for k in ['hit_1H', 'hit_13C', 'hit_COSY', 'hit_HSQC']:
            if k in Gf.nodes[n]:
                Gf.nodes[n][k] = int(Gf.nodes[n][k])
    nx.write_graphml(Gf, graph_out)
    with open(graph_out, 'r', encoding='utf-8') as f:
        content = f.read()
    with open(graph_out, 'w', encoding='utf-8') as f:
        f.write(content.replace('attr.type="long"', 'attr.type="int"'))
    log_progress(f"  ✓ GraphML : {graph_out}")

    # 2.9) Rapport
    report_path = os.path.join(FOLDER, 'adducts_detection_report.txt')
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("RAPPORT DE DÉTECTION DES ADDUCTS\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Tolérance : {TOL_PPM} ppm\n")
        f.write(f"Total molécules : {len(meta_df)}\n\n")
        f.write("MODE POSITIF :\n")
        for ad, cnt in adduct_pos_counts.items():
            f.write(f"  {ad} : {cnt} ({cnt/max(total_pos,1)*100:.1f}%)\n")
        f.write(f"  Total : {total_pos}\n\n")
        f.write("MODE NÉGATIF :\n")
        for ad, cnt in adduct_neg_counts.items():
            f.write(f"  {ad} : {cnt} ({cnt/max(total_neg,1)*100:.1f}%)\n")
        f.write(f"  Total : {total_neg}\n\n")
        f.write(f"TOTAL DÉTECTIONS : {total_mass} ({total_mass/len(meta_df)*100:.1f}%)\n")
    log_progress(f"  ✓ Rapport : {report_path}")

    log_progress("✅ Partie 2 terminée")
    return meta_merged


# ================================================================
# PARTIE 3 : FIGURE 5 (DISTRIBUTION 3 COUCHES)
# ================================================================

def generate_figure5(meta_merged, frames):
    log_progress("=" * 60)
    log_progress("PARTIE 3 : FIGURE 5 (DISTRIBUTION 3 COUCHES)")
    log_progress("=" * 60)

    # 3.1) Score combiné = MOYENNE des 4 expériences
    score_cols_fig = []
    for col_name in ['score_1H', 'score_13C', 'score_HSQC', 'score_COSY']:
        # Gérer les suffixes _x ou _y après merge
        if col_name in meta_merged.columns:
            score_cols_fig.append(col_name)
        else:
            alt = [c for c in meta_merged.columns if c.startswith(col_name)]
            if alt:
                meta_merged[col_name] = meta_merged[alt[0]]
                score_cols_fig.append(col_name)
                log_progress(f"  ⚠ '{col_name}' → utilisation de '{alt[0]}'")

    for col in score_cols_fig:
        meta_merged[col] = pd.to_numeric(meta_merged[col], errors='coerce').fillna(0)

    meta_merged['combined_score'] = meta_merged[score_cols_fig].mean(axis=1)

    log_progress(f"\n📊 Score combiné (moyenne des 4) :")
    log_progress(f"  Min={meta_merged['combined_score'].min():.1f}  Max={meta_merged['combined_score'].max():.1f}")
    log_progress(f"  Moy={meta_merged['combined_score'].mean():.1f}  Méd={meta_merged['combined_score'].median():.1f}")

    # 3.2) Identifier les groupes
    mask_ms = meta_merged['found_mass'].astype(str).str.lower() == 'yes'
    df_ms = meta_merged[mask_ms]
    log_progress(f"\n🟤 MS confirmés : {len(df_ms)}")

    # 3.2b) Séparer MS unique vs MS isomères (même masse exacte)
    mass_col = 'exact_mass' if 'exact_mass' in df_ms.columns else 'Mass'
    if mass_col in df_ms.columns:
        # Trouver les masses qui apparaissent plus d'une fois parmi les MS-confirmés
        mass_counts = df_ms[mass_col].value_counts()
        duplicate_masses = mass_counts[mass_counts > 1].index.tolist()

        # Séparer en isomères et uniques
        mask_isomers = df_ms[mass_col].isin(duplicate_masses)
        df_ms_isomers = df_ms[mask_isomers]
        df_ms_unique = df_ms[~mask_isomers]

        log_progress(f"  🟤 MS unique (masse unique) : {len(df_ms_unique)}")
        log_progress(f"  🔴 MS isomères (même masse) : {len(df_ms_isomers)} ({len(duplicate_masses)} masses partagées)")
    else:
        df_ms_unique = df_ms
        df_ms_isomers = pd.DataFrame()
        log_progress(f"  ⚠ Colonne masse non trouvée, tous MS traités comme uniques")

    # 3.2c) Sélection intelligente des isomères par NMR (seuil 5 points)
    # Pour chaque groupe d'isomères (même masse exacte) :
    #   - best : le meilleur combined_score
    #   - close : autres isomères dans les 5 pts du meilleur (score >= best - 5)
    #   - nonmatch : isomères avec score > 5 pts en dessous du meilleur
    NMR_ISOMER_THRESHOLD = 5.0  # seuil en points de score

    best_isomer_indices = set()
    close_isomer_indices = set()
    nonmatch_isomer_indices = set()
    nonmatch_isomer_sim = {}  # {idx: spectral_sim [0-1] vs best isomer du même groupe}

    if mass_col in df_ms.columns and len(df_ms_isomers) > 0:
        for mass_val, group in df_ms.groupby(mass_col):
            if len(group) <= 1:
                continue
            best_score = group['combined_score'].max()
            match_mask = group['combined_score'] >= (best_score - NMR_ISOMER_THRESHOLD)
            match_group = group[match_mask]
            nonmatch_group = group[~match_mask]

            best_idx = group['combined_score'].idxmax()
            best_isomer_indices.add(best_idx)

            for idx in match_group.index:
                if idx != best_idx:
                    close_isomer_indices.add(idx)

            for idx in nonmatch_group.index:
                nonmatch_isomer_indices.add(idx)

        log_progress(f"  ✅ Sélection NMR isomères (seuil ±{NMR_ISOMER_THRESHOLD} pts) :")
        log_progress(f"     Meilleur scoreur par groupe : {len(best_isomer_indices)}")
        log_progress(f"     Isomères proches (dans ±5 pts) : {len(close_isomer_indices)}")
        log_progress(f"     Isomères non-match (>5 pts) : {len(nonmatch_isomer_indices)}")
    else:
        log_progress(f"  ⚠ Pas d'isomères MS ou colonne masse manquante — sélection NMR ignorée")

    # Analogues : Tanimoto SPECTRAL moyen > seuil avec un MS-confirmé, SANS match MS
    log_progress(f"🔶 Identification des analogues (Tanimoto spectral moyen > {TANIMOTO_SPECTRAL})...")
    log_progress(f"   Tolérances : 1H={TOL_FIG5['1H'][0]} ppm, 13C={TOL_FIG5['13C'][0]} ppm")
    log_progress(f"   HSQC=({TOL_FIG5['HSQC'][0]}, {TOL_FIG5['HSQC'][1]}) ppm, COSY=({TOL_FIG5['COSY'][0]}, {TOL_FIG5['COSY'][1]}) ppm")

    # Construire un dictionnaire MolName -> ppm_molecule par méthode
    ppm_by_method = {}
    for method in ['1H', '13C', 'HSQC', 'COSY']:
        ppm_by_method[method] = {}
        if method in frames:
            for rec in frames[method]:
                mol_name = rec.get('MolName', '')
                ppm_str = rec.get('ppm_molecule', '')
                if mol_name:
                    ppm_by_method[method][mol_name] = ppm_str

    # Identifier la colonne MolName dans meta_merged
    mol_col = 'Source' if 'Source' in meta_merged.columns else 'MolName'

    # Pré-parser les shifts des MS-confirmés
    ms_shifts = {}
    for idx, row in df_ms.iterrows():
        mol_name = row.get(mol_col, '')
        ms_shifts[idx] = {
            '1H': parse_1D_signals(ppm_by_method['1H'].get(mol_name, '')),
            '13C': parse_1D_signals(ppm_by_method['13C'].get(mol_name, '')),
            'HSQC': parse_2D_signals(ppm_by_method['HSQC'].get(mol_name, '')),
            'COSY': parse_2D_signals(ppm_by_method['COSY'].get(mol_name, '')),
        }

    analog_indices = set()
    non_ms_df = meta_merged[~mask_ms]

    for idx, row in tqdm(non_ms_df.iterrows(), total=len(non_ms_df), desc="Analogues spectraux"):
        mol_name = row.get(mol_col, '')

        # Parser les shifts de cette molécule
        shifts_A = {
            '1H': parse_1D_signals(ppm_by_method['1H'].get(mol_name, '')),
            '13C': parse_1D_signals(ppm_by_method['13C'].get(mol_name, '')),
            'HSQC': parse_2D_signals(ppm_by_method['HSQC'].get(mol_name, '')),
            'COSY': parse_2D_signals(ppm_by_method['COSY'].get(mol_name, '')),
        }

        # Comparer avec chaque MS-confirmé
        for _, shifts_B in ms_shifts.items():
            # Calculer Tanimoto spectral pour chaque méthode
            tan_1H = spectral_tanimoto_1D(shifts_A['1H'], shifts_B['1H'], TOL_FIG5['1H'][0])
            tan_13C = spectral_tanimoto_1D(shifts_A['13C'], shifts_B['13C'], TOL_FIG5['13C'][0])
            tan_HSQC = spectral_tanimoto_2D(shifts_A['HSQC'], shifts_B['HSQC'],
                                            TOL_FIG5['HSQC'][0], TOL_FIG5['HSQC'][1])
            tan_COSY = spectral_tanimoto_2D(shifts_A['COSY'], shifts_B['COSY'],
                                            TOL_FIG5['COSY'][0], TOL_FIG5['COSY'][1])

            # Moyenne des 4 Tanimoto
            mean_tan = (tan_1H + tan_13C + tan_HSQC + tan_COSY) / 4

            if mean_tan >= TANIMOTO_SPECTRAL:
                analog_indices.add(idx)
                break  # Un seul match suffit

    df_analog_spectral = meta_merged.loc[list(analog_indices)]
    log_progress(f"  ✓ {len(df_analog_spectral)} analogues spectraux (Tanimoto > {TANIMOTO_SPECTRAL}, sans MS)")

    # 3.2e) Analogues ECFP4 (Tanimoto > 0.7 avec un MS-confirmé)
    log_progress(f"🔶 Identification des analogues ECFP4 (Tanimoto > {TANIMOTO_ECFP4})...")
    ecfp_col = 'ECFP4' if 'ECFP4' in meta_merged.columns else None
    if ecfp_col is None:
        ecfp_candidates = [c for c in meta_merged.columns if 'ECFP4' in c.upper()]
        if ecfp_candidates:
            ecfp_col = ecfp_candidates[0]

    analog_ecfp4_indices = set()

    if ecfp_col:
        # Pré-calculer les vecteurs ECFP4 des MS-confirmés
        ms_ecfp_vecs = {}
        for idx, row in df_ms.iterrows():
            v = bitstr_to_vec(row.get(ecfp_col, ''))
            if v is not None:
                ms_ecfp_vecs[idx] = v

        # Pour chaque non-MS, chercher si Tanimoto ECFP4 > seuil avec un MS-confirmé
        non_ms_df = meta_merged[~mask_ms]
        for idx, row in tqdm(non_ms_df.iterrows(), total=len(non_ms_df), desc="Analogues ECFP4"):
            v = bitstr_to_vec(row.get(ecfp_col, ''))
            if v is None:
                continue
            for _, ms_v in ms_ecfp_vecs.items():
                if len(v) == len(ms_v):
                    t = TanimotoSimilarity(v, ms_v)
                    if t >= TANIMOTO_ECFP4:
                        analog_ecfp4_indices.add(idx)
                        break

    df_analog_ecfp4 = meta_merged.loc[list(analog_ecfp4_indices)]
    log_progress(f"  ✓ {len(df_analog_ecfp4)} analogues ECFP4 (Tanimoto > {TANIMOTO_ECFP4}, sans MS)")

    # 3.2e-bis) Calculer la similarité spectrale des isomères non-match vs leur meilleur isomère
    if len(nonmatch_isomer_indices) > 0 and mass_col in df_ms.columns:
        log_progress(f"🔶 Calcul similarité spectrale des isomères non-match vs meilleur isomère...")
        ppm_by_method_iso = {}
        for method in ['1H', '13C', 'HSQC', 'COSY']:
            ppm_by_method_iso[method] = {}
            if method in frames:
                for rec in frames[method]:
                    mol_name = rec.get('MolName', '')
                    ppm_str = rec.get('ppm_molecule', '')
                    if mol_name:
                        ppm_by_method_iso[method][mol_name] = ppm_str

        mol_col_iso = 'Source' if 'Source' in meta_merged.columns else 'MolName'

        for mass_val, group in df_ms.groupby(mass_col):
            if len(group) <= 1:
                continue
            best_idx = group['combined_score'].idxmax()
            best_mol = group.loc[best_idx].get(mol_col_iso, '')
            shifts_best = {
                '1H':   parse_1D_signals(ppm_by_method_iso['1H'].get(best_mol, '')),
                '13C':  parse_1D_signals(ppm_by_method_iso['13C'].get(best_mol, '')),
                'HSQC': parse_2D_signals(ppm_by_method_iso['HSQC'].get(best_mol, '')),
                'COSY': parse_2D_signals(ppm_by_method_iso['COSY'].get(best_mol, '')),
            }
            for idx in group.index:
                if idx not in nonmatch_isomer_indices:
                    continue
                nm_mol = group.loc[idx].get(mol_col_iso, '')
                shifts_nm = {
                    '1H':   parse_1D_signals(ppm_by_method_iso['1H'].get(nm_mol, '')),
                    '13C':  parse_1D_signals(ppm_by_method_iso['13C'].get(nm_mol, '')),
                    'HSQC': parse_2D_signals(ppm_by_method_iso['HSQC'].get(nm_mol, '')),
                    'COSY': parse_2D_signals(ppm_by_method_iso['COSY'].get(nm_mol, '')),
                }
                t1 = spectral_tanimoto_1D(shifts_nm['1H'],   shifts_best['1H'],   TOL_FIG5['1H'][0])
                t2 = spectral_tanimoto_1D(shifts_nm['13C'],  shifts_best['13C'],  TOL_FIG5['13C'][0])
                t3 = spectral_tanimoto_2D(shifts_nm['HSQC'], shifts_best['HSQC'], TOL_FIG5['HSQC'][0], TOL_FIG5['HSQC'][1])
                t4 = spectral_tanimoto_2D(shifts_nm['COSY'], shifts_best['COSY'], TOL_FIG5['COSY'][0], TOL_FIG5['COSY'][1])
                nonmatch_isomer_sim[idx] = (t1 + t2 + t3 + t4) / 4

        log_progress(f"  ✓ Similarité calculée pour {len(nonmatch_isomer_sim)} isomères non-match")

    # DataFrames pour les nouvelles catégories
    df_ms_best_isomer     = df_ms.loc[list(best_isomer_indices)]      if best_isomer_indices      else pd.DataFrame()
    df_ms_close_isomer    = df_ms.loc[list(close_isomer_indices)]     if close_isomer_indices     else pd.DataFrame()
    df_ms_nonmatch_isomer = df_ms.loc[list(nonmatch_isomer_indices)]  if nonmatch_isomer_indices  else pd.DataFrame()

    if len(df_ms_nonmatch_isomer) > 0:
        df_ms_nonmatch_isomer = df_ms_nonmatch_isomer.copy()
        df_ms_nonmatch_isomer['spectral_sim_vs_best'] = [
            nonmatch_isomer_sim.get(i, 0.0) * 100 for i in df_ms_nonmatch_isomer.index
        ]

    # 3.2f) Identifier les non-matchés pour chaque méthode
    # Pour Spectral
    all_matched_spectral = set(df_ms.index) | analog_indices
    unmatched_spectral = set(meta_merged.index) - all_matched_spectral
    df_unmatched_spectral = meta_merged.loc[list(unmatched_spectral)]

    # Pour ECFP4
    all_matched_ecfp4 = set(df_ms.index) | analog_ecfp4_indices
    unmatched_ecfp4 = set(meta_merged.index) - all_matched_ecfp4
    df_unmatched_ecfp4 = meta_merged.loc[list(unmatched_ecfp4)]

    log_progress(f"🔵 Sans match spectral : {len(df_unmatched_spectral)}")
    log_progress(f"🔵 Sans match ECFP4 : {len(df_unmatched_ecfp4)}")
    log_progress(f"📊 Tous les candidats : {len(meta_merged)}")

    # 3.3) Tracé
    log_progress("🎨 Génération de la figure...")

    # Vérifier qu'on a des données valides
    if len(meta_merged) == 0 or meta_merged['combined_score'].isna().all():
        log_progress("⚠️ ERREUR: Aucune donnée valide pour générer la figure")
        log_progress("   Vérifiez que les noms de molécules matchent entre les fichiers")
        return

    # Bins de 20 à 100 par pas de 5
    bins = np.arange(20, 101, BIN_WIDTH)  # 20, 25, 30, ..., 95, 100

    # Fonction helper pour créer un histogramme
    def create_histogram(ax, df_ms_layer, df_analog_layer, df_unmatched_layer,
                         tanimoto_thresh, method_name, with_isomers=False):
        data_to_stack = []
        colors_to_use = []
        labels_to_use = []

        if with_isomers:
            # Avec distinction isomères
            if len(df_ms_unique) > 0:
                data_to_stack.append(df_ms_unique['combined_score'])
                colors_to_use.append(COLOR_MS)
                labels_to_use.append(f'MS confirmed unique (n = {len(df_ms_unique)})')
            if len(df_ms_isomers) > 0:
                data_to_stack.append(df_ms_isomers['combined_score'])
                colors_to_use.append(COLOR_MS_AMBIG)
                labels_to_use.append(f'MS confirmed isomers (n = {len(df_ms_isomers)})')
        else:
            # Sans distinction isomères
            if len(df_ms_layer) > 0:
                data_to_stack.append(df_ms_layer['combined_score'])
                colors_to_use.append(COLOR_MS)
                labels_to_use.append(f'MS confirmed (n = {len(df_ms_layer)})')

        if len(df_analog_layer) > 0:
            data_to_stack.append(df_analog_layer['combined_score'])
            colors_to_use.append(COLOR_ANALOG)
            labels_to_use.append(f'MS non confirmed + {method_name} (n = {len(df_analog_layer)})')

        if len(df_unmatched_layer) > 0:
            data_to_stack.append(df_unmatched_layer['combined_score'])
            colors_to_use.append(COLOR_UNMATCHED)
            labels_to_use.append(f'MS non confirmed (n = {len(df_unmatched_layer)})')

        if len(data_to_stack) > 0:
            # Tracer l'histogramme empilé et récupérer les données
            n_vals, bins_out, _ = ax.hist(data_to_stack, bins=bins, stacked=True,
                    color=colors_to_use, alpha=0.85,
                    edgecolor='white', linewidth=0.5,
                    label=labels_to_use)

            # Ajouter les labels "n" sur chaque barre (total de la barre)
            if len(data_to_stack) == 1:
                totals = n_vals
            else:
                totals = n_vals[-1]

            bin_centers = (bins_out[:-1] + bins_out[1:]) / 2
            for count, x in zip(totals, bin_centers):
                if count > 0:
                    ax.text(x, count + 0.5, f'{int(count)}', ha='center', va='bottom',
                            fontsize=7, fontweight='bold', color='#333333')

        ax.hist(meta_merged['combined_score'], bins=bins, histtype='step',
                color=COLOR_ALL, linewidth=2,
                label=f'All candidates (n = {len(meta_merged)})', zorder=4)

        ax.set_xlabel('Score coverage (mean of 4 experiments)', fontsize=11)
        ax.set_ylabel('Number of candidates', fontsize=11)
        ax.set_xlim(20, 100)
        ax.set_xticks(np.arange(20, 101, 5))
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.18),
                  fontsize=7, frameon=False, ncol=3)
        ax.spines['top'].set_visible(False)
        ax.tick_params(labelsize=9)

    # Fonction helper pour créer un diagramme en barres empilées
    def create_barchart(ax, df_ms_layer, df_analog_layer, df_unmatched_layer,
                        tanimoto_thresh, method_name, with_isomers=False):
        # Calculer les counts par bin pour chaque catégorie
        bin_edges = bins
        bin_labels = [f'{int(bins[i])}-{int(bins[i+1])}' for i in range(len(bins)-1)]
        bar_width = 0.8

        # Préparer les données
        counts_list = []
        colors_list = []
        labels_list = []

        if with_isomers:
            if len(df_ms_unique) > 0:
                counts, _ = np.histogram(df_ms_unique['combined_score'], bins=bin_edges)
                counts_list.append(counts)
                colors_list.append(COLOR_MS)
                labels_list.append(f'MS confirmed unique (n = {len(df_ms_unique)})')
            if len(df_ms_isomers) > 0:
                counts, _ = np.histogram(df_ms_isomers['combined_score'], bins=bin_edges)
                counts_list.append(counts)
                colors_list.append(COLOR_MS_AMBIG)
                labels_list.append(f'MS confirmed isomers (n = {len(df_ms_isomers)})')
        else:
            if len(df_ms_layer) > 0:
                counts, _ = np.histogram(df_ms_layer['combined_score'], bins=bin_edges)
                counts_list.append(counts)
                colors_list.append(COLOR_MS)
                labels_list.append(f'MS confirmed (n = {len(df_ms_layer)})')

        if len(df_analog_layer) > 0:
            counts, _ = np.histogram(df_analog_layer['combined_score'], bins=bin_edges)
            counts_list.append(counts)
            colors_list.append(COLOR_ANALOG)
            labels_list.append(f'MS non confirmed + {method_name} (n = {len(df_analog_layer)})')

        if len(df_unmatched_layer) > 0:
            counts, _ = np.histogram(df_unmatched_layer['combined_score'], bins=bin_edges)
            counts_list.append(counts)
            colors_list.append(COLOR_UNMATCHED)
            labels_list.append(f'MS non confirmed (n = {len(df_unmatched_layer)})')

        # Tracer les barres empilées
        x = np.arange(len(bin_labels))
        bottom = np.zeros(len(bin_labels))

        for counts, color, label in zip(counts_list, colors_list, labels_list):
            ax.bar(x, counts, bar_width, bottom=bottom, color=color,
                   edgecolor='white', linewidth=0.5, label=label, alpha=0.85)
            bottom += counts

        # Ajouter les labels "n" sur chaque barre
        for i, total in enumerate(bottom):
            if total > 0:
                ax.text(i, total + 0.5, f'{int(total)}', ha='center', va='bottom',
                        fontsize=7, fontweight='bold', color='#333333')

        ax.set_xlabel('Score coverage (mean of 4 experiments)', fontsize=11)
        ax.set_ylabel('Number of candidates', fontsize=11)
        ax.set_xticks(x)
        ax.set_xticklabels(bin_labels, rotation=45, ha='right', fontsize=8)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.28),
                  fontsize=7, frameon=False, ncol=3)
        ax.spines['top'].set_visible(False)
        ax.tick_params(labelsize=9)

    # Fonction helper pour créer un diagramme en barres NORMALISÉ À 100%
    def create_barchart_normalized(ax, df_ms_layer, df_analog_layer, df_unmatched_layer,
                                    tanimoto_thresh, method_name, with_isomers=False):
        """Diagramme en barres empilées normalisé à 100%."""
        bin_edges = bins
        bin_labels = [f'{int(bins[i])}-{int(bins[i+1])}' for i in range(len(bins)-1)]
        bar_width = 0.8

        # Préparer les données
        counts_list = []
        colors_list = []
        labels_list = []

        if with_isomers:
            if len(df_ms_unique) > 0:
                counts, _ = np.histogram(df_ms_unique['combined_score'], bins=bin_edges)
                counts_list.append(counts.astype(float))
                colors_list.append(COLOR_MS)
                labels_list.append(f'MS confirmed unique (n = {len(df_ms_unique)})')
            if len(df_ms_isomers) > 0:
                counts, _ = np.histogram(df_ms_isomers['combined_score'], bins=bin_edges)
                counts_list.append(counts.astype(float))
                colors_list.append(COLOR_MS_AMBIG)
                labels_list.append(f'MS confirmed isomers (n = {len(df_ms_isomers)})')
        else:
            if len(df_ms_layer) > 0:
                counts, _ = np.histogram(df_ms_layer['combined_score'], bins=bin_edges)
                counts_list.append(counts.astype(float))
                colors_list.append(COLOR_MS)
                labels_list.append(f'MS confirmed (n = {len(df_ms_layer)})')

        if len(df_analog_layer) > 0:
            counts, _ = np.histogram(df_analog_layer['combined_score'], bins=bin_edges)
            counts_list.append(counts.astype(float))
            colors_list.append(COLOR_ANALOG)
            labels_list.append(f'MS non confirmed + {method_name} (n = {len(df_analog_layer)})')

        if len(df_unmatched_layer) > 0:
            counts, _ = np.histogram(df_unmatched_layer['combined_score'], bins=bin_edges)
            counts_list.append(counts.astype(float))
            colors_list.append(COLOR_UNMATCHED)
            labels_list.append(f'MS non confirmed (n = {len(df_unmatched_layer)})')

        # Calculer les totaux par bin
        totals = np.zeros(len(bin_labels))
        for counts in counts_list:
            totals += counts

        # Normaliser à 100%
        normalized_list = []
        for counts in counts_list:
            normalized = np.zeros_like(counts)
            mask = totals > 0
            normalized[mask] = (counts[mask] / totals[mask]) * 100
            normalized_list.append(normalized)

        # Tracer les barres empilées normalisées
        x = np.arange(len(bin_labels))
        bottom = np.zeros(len(bin_labels))

        for normalized, color, label in zip(normalized_list, colors_list, labels_list):
            ax.bar(x, normalized, bar_width, bottom=bottom, color=color,
                   edgecolor='white', linewidth=0.5, label=label, alpha=0.85)
            bottom += normalized

        # Ajouter les labels "n" (nombre total) sur chaque barre
        for i, total in enumerate(totals):
            if total > 0:
                ax.text(i, 101, f'n={int(total)}', ha='center', va='bottom',
                        fontsize=6, fontweight='bold', color='#333333')

        ax.set_xlabel('Score coverage (mean of 4 experiments)', fontsize=11)
        ax.set_ylabel('Percentage (%)', fontsize=11)
        ax.set_ylim(0, 110)
        ax.set_yticks([0, 20, 40, 60, 80, 100])
        ax.set_xticks(x)
        ax.set_xticklabels(bin_labels, rotation=45, ha='right', fontsize=8)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.18),
                  fontsize=7, frameon=False, ncol=3)
        ax.spines['top'].set_visible(False)
        ax.tick_params(labelsize=9)

    # ================================================================
    # HISTOGRAMMES SPECTRAL
    # ================================================================
    log_progress(f"\n📊 HISTOGRAMMES SPECTRAL (seuil={TANIMOTO_SPECTRAL}, MS={TOL_PPM_SPECTRAL}ppm)")

    # Histogramme 1: Spectral avec isomères
    log_progress("  📊 Histogramme Spectral avec isomères...")
    fig1, ax1 = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))
    create_histogram(ax1, df_ms, df_analog_spectral, df_unmatched_spectral,
                     TANIMOTO_SPECTRAL, "Spectral sim.", with_isomers=True)
    plt.tight_layout()
    for ext, fmt in [('png', 'png'), ('tiff', 'tiff'), ('svg', 'svg'), ('pdf', 'pdf')]:
        out = os.path.join(FOLDER, f'Figure5_hist_spectral_with_isomers.{ext}')
        fig1.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor='white', format=fmt)
        log_progress(f"    ✓ {ext.upper()} : {out}")
    plt.close(fig1)

    # Histogramme 2: Spectral sans isomères
    log_progress("  📊 Histogramme Spectral sans isomères...")
    fig2, ax2 = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))
    create_histogram(ax2, df_ms, df_analog_spectral, df_unmatched_spectral,
                     TANIMOTO_SPECTRAL, "Spectral sim.", with_isomers=False)
    plt.tight_layout()
    for ext, fmt in [('png', 'png'), ('tiff', 'tiff'), ('svg', 'svg'), ('pdf', 'pdf')]:
        out = os.path.join(FOLDER, f'Figure5_hist_spectral_without_isomers.{ext}')
        fig2.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor='white', format=fmt)
        log_progress(f"    ✓ {ext.upper()} : {out}")
    plt.close(fig2)

    # ================================================================
    # BARRES SPECTRAL
    # ================================================================
    log_progress(f"\n📊 BARRES SPECTRAL (seuil={TANIMOTO_SPECTRAL}, MS={TOL_PPM_SPECTRAL}ppm)")

    # Barres 1: Spectral avec isomères
    log_progress("  📊 Barres Spectral avec isomères...")
    fig_b1, ax_b1 = plt.subplots(figsize=(FIG_WIDTH + 2, FIG_HEIGHT))
    create_barchart(ax_b1, df_ms, df_analog_spectral, df_unmatched_spectral,
                    TANIMOTO_SPECTRAL, "Spectral sim.", with_isomers=True)
    plt.tight_layout()
    for ext, fmt in [('png', 'png'), ('tiff', 'tiff'), ('svg', 'svg'), ('pdf', 'pdf')]:
        out = os.path.join(FOLDER, f'Figure5_bar_spectral_with_isomers.{ext}')
        fig_b1.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor='white', format=fmt)
        log_progress(f"    ✓ {ext.upper()} : {out}")
    plt.close(fig_b1)

    # Barres 2: Spectral sans isomères
    log_progress("  📊 Barres Spectral sans isomères...")
    fig_b2, ax_b2 = plt.subplots(figsize=(FIG_WIDTH + 2, FIG_HEIGHT))
    create_barchart(ax_b2, df_ms, df_analog_spectral, df_unmatched_spectral,
                    TANIMOTO_SPECTRAL, "Spectral sim.", with_isomers=False)
    plt.tight_layout()
    for ext, fmt in [('png', 'png'), ('tiff', 'tiff'), ('svg', 'svg'), ('pdf', 'pdf')]:
        out = os.path.join(FOLDER, f'Figure5_bar_spectral_without_isomers.{ext}')
        fig_b2.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor='white', format=fmt)
        log_progress(f"    ✓ {ext.upper()} : {out}")
    plt.close(fig_b2)

    # ================================================================
    # HISTOGRAMMES ECFP4
    # ================================================================
    log_progress(f"\n📊 HISTOGRAMMES ECFP4 (seuil={TANIMOTO_ECFP4}, MS={TOL_PPM_ECFP4}ppm)")

    # Histogramme 3: ECFP4 avec isomères
    log_progress("  📊 Histogramme ECFP4 avec isomères...")
    fig3, ax3 = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))
    create_histogram(ax3, df_ms, df_analog_ecfp4, df_unmatched_ecfp4,
                     TANIMOTO_ECFP4, "ECFP4 Tan.", with_isomers=True)
    plt.tight_layout()
    for ext, fmt in [('png', 'png'), ('tiff', 'tiff'), ('svg', 'svg'), ('pdf', 'pdf')]:
        out = os.path.join(FOLDER, f'Figure5_hist_ECFP4_with_isomers.{ext}')
        fig3.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor='white', format=fmt)
        log_progress(f"    ✓ {ext.upper()} : {out}")
    plt.close(fig3)

    # Histogramme 4: ECFP4 sans isomères
    log_progress("  📊 Histogramme ECFP4 sans isomères...")
    fig4, ax4 = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))
    create_histogram(ax4, df_ms, df_analog_ecfp4, df_unmatched_ecfp4,
                     TANIMOTO_ECFP4, "ECFP4 Tan.", with_isomers=False)
    plt.tight_layout()
    for ext, fmt in [('png', 'png'), ('tiff', 'tiff'), ('svg', 'svg'), ('pdf', 'pdf')]:
        out = os.path.join(FOLDER, f'Figure5_hist_ECFP4_without_isomers.{ext}')
        fig4.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor='white', format=fmt)
        log_progress(f"    ✓ {ext.upper()} : {out}")
    plt.close(fig4)

    # ================================================================
    # BARRES ECFP4
    # ================================================================
    log_progress(f"\n📊 BARRES ECFP4 (seuil={TANIMOTO_ECFP4}, MS={TOL_PPM_ECFP4}ppm)")

    # Barres 3: ECFP4 avec isomères
    log_progress("  📊 Barres ECFP4 avec isomères...")
    fig_b3, ax_b3 = plt.subplots(figsize=(FIG_WIDTH + 2, FIG_HEIGHT))
    create_barchart(ax_b3, df_ms, df_analog_ecfp4, df_unmatched_ecfp4,
                    TANIMOTO_ECFP4, "ECFP4 Tan.", with_isomers=True)
    plt.tight_layout()
    for ext, fmt in [('png', 'png'), ('tiff', 'tiff'), ('svg', 'svg'), ('pdf', 'pdf')]:
        out = os.path.join(FOLDER, f'Figure5_bar_ECFP4_with_isomers.{ext}')
        fig_b3.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor='white', format=fmt)
        log_progress(f"    ✓ {ext.upper()} : {out}")
    plt.close(fig_b3)

    # Barres 4: ECFP4 sans isomères
    log_progress("  📊 Barres ECFP4 sans isomères...")
    fig_b4, ax_b4 = plt.subplots(figsize=(FIG_WIDTH + 2, FIG_HEIGHT))
    create_barchart(ax_b4, df_ms, df_analog_ecfp4, df_unmatched_ecfp4,
                    TANIMOTO_ECFP4, "ECFP4 Tan.", with_isomers=False)
    plt.tight_layout()
    for ext, fmt in [('png', 'png'), ('tiff', 'tiff'), ('svg', 'svg'), ('pdf', 'pdf')]:
        out = os.path.join(FOLDER, f'Figure5_bar_ECFP4_without_isomers.{ext}')
        fig_b4.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor='white', format=fmt)
        log_progress(f"    ✓ {ext.upper()} : {out}")
    plt.close(fig_b4)

    # ================================================================
    # BARRES NORMALISÉES À 100% - SPECTRAL
    # ================================================================
    log_progress(f"\n📊 BARRES NORMALISÉES 100% - SPECTRAL")

    # Barres normalisées Spectral sans isomères
    log_progress("  📊 Barres normalisées Spectral sans isomères...")
    fig_n1, ax_n1 = plt.subplots(figsize=(FIG_WIDTH + 2, FIG_HEIGHT))
    create_barchart_normalized(ax_n1, df_ms, df_analog_spectral, df_unmatched_spectral,
                               TANIMOTO_SPECTRAL, "Spectral sim.", with_isomers=False)
    plt.tight_layout()
    for ext, fmt in [('png', 'png'), ('tiff', 'tiff'), ('svg', 'svg'), ('pdf', 'pdf')]:
        out = os.path.join(FOLDER, f'Figure5_bar_normalized_spectral.{ext}')
        fig_n1.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor='white', format=fmt)
        log_progress(f"    ✓ {ext.upper()} : {out}")
    plt.close(fig_n1)

    # Barres normalisées Spectral avec isomères
    log_progress("  📊 Barres normalisées Spectral avec isomères...")
    fig_n2, ax_n2 = plt.subplots(figsize=(FIG_WIDTH + 2, FIG_HEIGHT))
    create_barchart_normalized(ax_n2, df_ms, df_analog_spectral, df_unmatched_spectral,
                               TANIMOTO_SPECTRAL, "Spectral sim.", with_isomers=True)
    plt.tight_layout()
    for ext, fmt in [('png', 'png'), ('tiff', 'tiff'), ('svg', 'svg'), ('pdf', 'pdf')]:
        out = os.path.join(FOLDER, f'Figure5_bar_normalized_spectral_isomers.{ext}')
        fig_n2.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor='white', format=fmt)
        log_progress(f"    ✓ {ext.upper()} : {out}")
    plt.close(fig_n2)

    # ================================================================
    # BARRES NORMALISÉES À 100% - ECFP4
    # ================================================================
    log_progress(f"\n📊 BARRES NORMALISÉES 100% - ECFP4")

    # Barres normalisées ECFP4 sans isomères
    log_progress("  📊 Barres normalisées ECFP4 sans isomères...")
    fig_n3, ax_n3 = plt.subplots(figsize=(FIG_WIDTH + 2, FIG_HEIGHT))
    create_barchart_normalized(ax_n3, df_ms, df_analog_ecfp4, df_unmatched_ecfp4,
                               TANIMOTO_ECFP4, "ECFP4 Tan.", with_isomers=False)
    plt.tight_layout()
    for ext, fmt in [('png', 'png'), ('tiff', 'tiff'), ('svg', 'svg'), ('pdf', 'pdf')]:
        out = os.path.join(FOLDER, f'Figure5_bar_normalized_ECFP4.{ext}')
        fig_n3.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor='white', format=fmt)
        log_progress(f"    ✓ {ext.upper()} : {out}")
    plt.close(fig_n3)

    # Barres normalisées ECFP4 avec isomères
    log_progress("  📊 Barres normalisées ECFP4 avec isomères...")
    fig_n4, ax_n4 = plt.subplots(figsize=(FIG_WIDTH + 2, FIG_HEIGHT))
    create_barchart_normalized(ax_n4, df_ms, df_analog_ecfp4, df_unmatched_ecfp4,
                               TANIMOTO_ECFP4, "ECFP4 Tan.", with_isomers=True)
    plt.tight_layout()
    for ext, fmt in [('png', 'png'), ('tiff', 'tiff'), ('svg', 'svg'), ('pdf', 'pdf')]:
        out = os.path.join(FOLDER, f'Figure5_bar_normalized_ECFP4_isomers.{ext}')
        fig_n4.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor='white', format=fmt)
        log_progress(f"    ✓ {ext.upper()} : {out}")
    plt.close(fig_n4)

    # ================================================================
    # NOUVELLES FIGURES — SÉLECTION INTELLIGENTE DES ISOMÈRES NMR
    # best isomère = meilleur NMR par groupe de masse
    # close isomère = dans ±5 pts du best (comptés aussi)
    # nonmatch isomère = >5 pts en dessous du best → sim spectrale calculée
    # ================================================================
    COLOR_BEST_ISO     = '#2ECC71'  # vert
    COLOR_CLOSE_ISO    = '#F39C12'  # orange
    COLOR_NONMATCH_ISO = '#E74C3C'  # rouge

    def create_barchart_normalized_smart_isomers(ax, df_best, df_close, _df_nonmatch,
                                                  df_unique, df_analog, df_unmatched,
                                                  method_name):
        """Barres normalisées 100% : MS unique | best isomère NMR | close isomère ±5 pts
           | analogue (méthode) | sans match.
           Les non-match isomères sont exclus ici (figure dédiée ci-dessous).
        """
        bin_edges = bins
        bin_labels = [f'{int(bins[i])}-{int(bins[i+1])}' for i in range(len(bins)-1)]
        bar_width = 0.8
        counts_list, colors_list, labels_list = [], [], []

        if len(df_unique) > 0:
            c, _ = np.histogram(df_unique['combined_score'], bins=bin_edges)
            counts_list.append(c.astype(float)); colors_list.append(COLOR_MS)
            labels_list.append(f'MS unique (n={len(df_unique)})')
        if len(df_best) > 0:
            c, _ = np.histogram(df_best['combined_score'], bins=bin_edges)
            counts_list.append(c.astype(float)); colors_list.append(COLOR_BEST_ISO)
            labels_list.append(f'MS best isomère NMR (n={len(df_best)})')
        if len(df_close) > 0:
            c, _ = np.histogram(df_close['combined_score'], bins=bin_edges)
            counts_list.append(c.astype(float)); colors_list.append(COLOR_CLOSE_ISO)
            labels_list.append(f'MS isomère proche ±5 pts (n={len(df_close)})')
        if len(df_analog) > 0:
            c, _ = np.histogram(df_analog['combined_score'], bins=bin_edges)
            counts_list.append(c.astype(float)); colors_list.append(COLOR_ANALOG)
            labels_list.append(f'Analogue {method_name} (n={len(df_analog)})')
        if len(df_unmatched) > 0:
            c, _ = np.histogram(df_unmatched['combined_score'], bins=bin_edges)
            counts_list.append(c.astype(float)); colors_list.append(COLOR_UNMATCHED)
            labels_list.append(f'Sans match (n={len(df_unmatched)})')

        totals = np.zeros(len(bin_labels))
        for c in counts_list:
            totals += c

        x = np.arange(len(bin_labels))
        bottom = np.zeros(len(bin_labels))
        for c, color, label in zip(counts_list, colors_list, labels_list):
            norm = np.zeros_like(c)
            mask = totals > 0
            norm[mask] = (c[mask] / totals[mask]) * 100
            ax.bar(x, norm, bar_width, bottom=bottom, color=color,
                   edgecolor='white', linewidth=0.5, label=label, alpha=0.85)
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
        ax.legend(loc='upper left', fontsize=9, frameon=True, edgecolor='#CCCCCC')
        ax.spines['top'].set_visible(False)
        ax.tick_params(labelsize=9)
        ax.axhline(y=100, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

    def create_barchart_nonmatch_isomer_similarity(ax, df_nonmatch):
        """Distribution de la similarité spectrale des isomères non-match
        vs leur meilleur isomère NMR dans le même groupe de masse (bins 0-100%).
        """
        if len(df_nonmatch) == 0 or 'spectral_sim_vs_best' not in df_nonmatch.columns:
            ax.text(0.5, 0.5, 'Aucun isomère non-match', transform=ax.transAxes,
                    ha='center', va='center', fontsize=12, color='gray')
            ax.set_xlabel('Spectral similarity vs best NMR isomer (%)', fontsize=11)
            ax.set_ylabel('Number of isomers', fontsize=11)
            return

        sim_bins = np.arange(0, 105, 10)
        sim_labels = [f'{int(sim_bins[i])}-{int(sim_bins[i+1])}%' for i in range(len(sim_bins)-1)]
        counts, _ = np.histogram(df_nonmatch['spectral_sim_vs_best'], bins=sim_bins)

        x = np.arange(len(sim_labels))
        ax.bar(x, counts, 0.8, color=COLOR_NONMATCH_ISO, edgecolor='white',
               linewidth=0.5, alpha=0.85, label=f'Non-match isomers (n={len(df_nonmatch)})')

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
            f'Spectral similarity of non-match isomers (score >{NMR_ISOMER_THRESHOLD} pts below best)',
            fontsize=10, pad=8)

    # Figure smart isomères — Spectral
    if len(best_isomer_indices) > 0 or len(close_isomer_indices) > 0:
        log_progress(f"\n📊 BARRES NORMALISÉES — SÉLECTION INTELLIGENTE ISOMÈRES NMR (Spectral)")
        fig_si1, ax_si1 = plt.subplots(figsize=(FIG_WIDTH + 2, FIG_HEIGHT))
        create_barchart_normalized_smart_isomers(
            ax_si1, df_ms_best_isomer, df_ms_close_isomer, df_ms_nonmatch_isomer,
            df_ms_unique, df_analog_spectral, df_unmatched_spectral, "Spectral sim.")
        ax_si1.set_title(
            f'NMR isomer selection (±{NMR_ISOMER_THRESHOLD} pts) — Spectral analogs (Tanimoto ≥ {TANIMOTO_SPECTRAL})',
            fontsize=10, pad=8)
        plt.tight_layout()
        for ext, fmt in [('png', 'png'), ('tiff', 'tiff'), ('svg', 'svg'), ('pdf', 'pdf')]:
            out = os.path.join(FOLDER, f'Figure5_bar_normalized_smart_isomers_spectral.{ext}')
            fig_si1.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor='white', format=fmt)
            log_progress(f"    ✓ {ext.upper()} : {out}")
        plt.close(fig_si1)

        # Figure smart isomères — ECFP4
        log_progress(f"\n📊 BARRES NORMALISÉES — SÉLECTION INTELLIGENTE ISOMÈRES NMR (ECFP4)")
        fig_si2, ax_si2 = plt.subplots(figsize=(FIG_WIDTH + 2, FIG_HEIGHT))
        create_barchart_normalized_smart_isomers(
            ax_si2, df_ms_best_isomer, df_ms_close_isomer, df_ms_nonmatch_isomer,
            df_ms_unique, df_analog_ecfp4, df_unmatched_ecfp4, "ECFP4 Tan.")
        ax_si2.set_title(
            f'NMR isomer selection (±{NMR_ISOMER_THRESHOLD} pts) — ECFP4 analogs (Tanimoto ≥ {TANIMOTO_ECFP4})',
            fontsize=10, pad=8)
        plt.tight_layout()
        for ext, fmt in [('png', 'png'), ('tiff', 'tiff'), ('svg', 'svg'), ('pdf', 'pdf')]:
            out = os.path.join(FOLDER, f'Figure5_bar_normalized_smart_isomers_ECFP4.{ext}')
            fig_si2.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor='white', format=fmt)
            log_progress(f"    ✓ {ext.upper()} : {out}")
        plt.close(fig_si2)
    else:
        log_progress(f"\n⚠ Pas d'isomères NMR détectés — figures smart isomères ignorées")

    # Figure similarité spectrale des non-match isomers
    if len(df_ms_nonmatch_isomer) > 0:
        log_progress(f"\n📊 DISTRIBUTION SIMILARITÉ SPECTRALE — ISOMÈRES NON-MATCH NMR")
        fig_nm, ax_nm = plt.subplots(figsize=(FIG_WIDTH + 2, FIG_HEIGHT))
        create_barchart_nonmatch_isomer_similarity(ax_nm, df_ms_nonmatch_isomer)
        plt.tight_layout()
        for ext, fmt in [('png', 'png'), ('tiff', 'tiff'), ('svg', 'svg'), ('pdf', 'pdf')]:
            out = os.path.join(FOLDER, f'Figure5_nonmatch_isomers_spectral_similarity.{ext}')
            fig_nm.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor='white', format=fmt)
            log_progress(f"    ✓ {ext.upper()} : {out}")
        plt.close(fig_nm)
    else:
        log_progress(f"\n⚠ Pas d'isomères non-match NMR — figure similarité ignorée")

    # 3.5) Statistiques pour l'article
    log_progress("\n" + "=" * 60)
    log_progress("STATISTIQUES POUR L'ARTICLE")
    log_progress("=" * 60)
    log_progress(f"Score combiné = moyenne(S_1H, S_13C, S_HSQC, S_COSY)")
    log_progress(f"\nTous (n={len(meta_merged)}) : moy={meta_merged['combined_score'].mean():.1f}, méd={meta_merged['combined_score'].median():.1f}")
    if len(df_ms) > 0:
        log_progress(f"MS confirmés (n={len(df_ms)}) : moy={df_ms['combined_score'].mean():.1f}, méd={df_ms['combined_score'].median():.1f}")
    if len(df_analog_spectral) > 0:
        log_progress(f"Analogues spectraux (n={len(df_analog_spectral)}) : moy={df_analog_spectral['combined_score'].mean():.1f}, méd={df_analog_spectral['combined_score'].median():.1f}")
    if len(df_analog_ecfp4) > 0:
        log_progress(f"Analogues ECFP4 (n={len(df_analog_ecfp4)}) : moy={df_analog_ecfp4['combined_score'].mean():.1f}, méd={df_analog_ecfp4['combined_score'].median():.1f}")
    if len(df_unmatched_spectral) > 0:
        log_progress(f"Sans match spectral (n={len(df_unmatched_spectral)}) : moy={df_unmatched_spectral['combined_score'].mean():.1f}, méd={df_unmatched_spectral['combined_score'].median():.1f}")

    # Test statistique
    if len(df_ms) > 5:
        try:
            from scipy.stats import mannwhitneyu
            stat, pval = mannwhitneyu(
                df_ms['combined_score'],
                meta_merged[~mask_ms]['combined_score'],
                alternative='greater'
            )
            log_progress(f"\nMann-Whitney (MS vs non-MS) : U={stat:.0f}, p={pval:.2e}")
            if pval < 0.001:
                log_progress(f"  → Hautement significatif (p < 0.001)")
        except ImportError:
            log_progress("  scipy non installé — test Mann-Whitney ignoré")

    # 3.6) Export CSV détaillé des candidats (2 fichiers : Spectral et ECFP4)
    log_progress("\n💾 Export CSV détaillé des candidats...")

    # Colonnes de base à exporter
    base_cols = ['Source']

    # Scores par méthode
    for col in ['score_1H', 'score_13C', 'score_HSQC', 'score_COSY']:
        if col in meta_merged.columns:
            base_cols.append(col)

    # Score combiné
    base_cols.append('combined_score')

    # MS positif
    if 'found_pos' in meta_merged.columns:
        base_cols.append('found_pos')
    if 'adduct_pos' in meta_merged.columns:
        base_cols.append('adduct_pos')
    if 'matched_mz_pos' in meta_merged.columns:
        base_cols.append('matched_mz_pos')
    if 'diff_ppm_pos' in meta_merged.columns:
        base_cols.append('diff_ppm_pos')

    # MS négatif
    if 'found_neg' in meta_merged.columns:
        base_cols.append('found_neg')
    if 'adduct_neg' in meta_merged.columns:
        base_cols.append('adduct_neg')
    if 'matched_mz_neg' in meta_merged.columns:
        base_cols.append('matched_mz_neg')
    if 'diff_ppm_neg' in meta_merged.columns:
        base_cols.append('diff_ppm_neg')

    # MS global
    if 'found_mass' in meta_merged.columns:
        base_cols.append('found_mass')

    # Masse exacte
    if 'exact_mass' in meta_merged.columns:
        base_cols.append('exact_mass')

    # ========== CSV 1 : SPECTRAL ==========
    meta_merged['category_spectral'] = 'No match'
    meta_merged.loc[df_ms_unique.index, 'category_spectral'] = 'MS_unique'
    meta_merged.loc[df_ms_isomers.index, 'category_spectral'] = 'MS_isomers'
    meta_merged.loc[list(analog_indices), 'category_spectral'] = 'Spectral_analog'

    export_cols_spectral = base_cols + ['category_spectral']
    export_cols_spectral = [c for c in export_cols_spectral if c in meta_merged.columns]

    csv_spectral = os.path.join(FOLDER, 'Figure5_spectral_candidates.csv')
    meta_merged[export_cols_spectral].to_csv(csv_spectral, index=False)
    log_progress(f"  ✓ CSV Spectral : {csv_spectral}")

    # ========== CSV 2 : ECFP4 ==========
    meta_merged['category_ecfp4'] = 'No match'
    meta_merged.loc[df_ms_unique.index, 'category_ecfp4'] = 'MS_unique'
    meta_merged.loc[df_ms_isomers.index, 'category_ecfp4'] = 'MS_isomers'
    meta_merged.loc[list(analog_ecfp4_indices), 'category_ecfp4'] = 'ECFP4_analog'

    export_cols_ecfp4 = base_cols + ['category_ecfp4']
    export_cols_ecfp4 = [c for c in export_cols_ecfp4 if c in meta_merged.columns]

    csv_ecfp4 = os.path.join(FOLDER, 'Figure5_ECFP4_candidates.csv')
    meta_merged[export_cols_ecfp4].to_csv(csv_ecfp4, index=False)
    log_progress(f"  ✓ CSV ECFP4 : {csv_ecfp4}")

    log_progress(f"    {len(meta_merged)} candidats exportés")

    log_progress("✅ Partie 3 terminée")


# ================================================================
# MAIN — EXÉCUTION SÉQUENTIELLE
# ================================================================

if __name__ == '__main__':
    start_time = time.time()

    log_progress("🚀 NMRSCOUT PIPELINE COMPLET")
    log_progress(f"   Dossier : {FOLDER}")

    # Partie 1 : Réseau spectral
    Gf, nodes_df, coverage_map, mols, frames = build_spectral_network()

    # Partie 2 : Validation MS
    meta_merged = run_ms_validation(Gf, nodes_df)

    # Partie 3 : Figure 5
    generate_figure5(meta_merged, frames)

    elapsed = time.time() - start_time
    log_progress(f"\n{'='*60}")
    log_progress(f"🏁 PIPELINE COMPLET TERMINÉ EN {elapsed:.1f}s ({elapsed/60:.1f} min)")
    log_progress(f"{'='*60}")
    log_progress("Fichiers générés :")
    log_progress("  1. clusters_nodes_structural_first.csv")
    log_progress("  2. spectral_network_structural_first.graphml")
    log_progress("  3. metadata_with_found_mass_with_adducts.csv")
    log_progress("  4. spectral_network_annot_with_adducts.graphml")
    log_progress("  5. adducts_detection_report.txt")
    log_progress("  HISTOGRAMMES:")
    log_progress("  6. Figure5_hist_spectral_with_isomers.png/.tiff/.svg/.pdf")
    log_progress("  7. Figure5_hist_spectral_without_isomers.png/.tiff/.svg/.pdf")
    log_progress("  8. Figure5_hist_ECFP4_with_isomers.png/.tiff/.svg/.pdf")
    log_progress("  9. Figure5_hist_ECFP4_without_isomers.png/.tiff/.svg/.pdf")
    log_progress("  DIAGRAMMES EN BARRES:")
    log_progress(" 10. Figure5_bar_spectral_with_isomers.png/.tiff/.svg/.pdf")
    log_progress(" 11. Figure5_bar_spectral_without_isomers.png/.tiff/.svg/.pdf")
    log_progress(" 12. Figure5_bar_ECFP4_with_isomers.png/.tiff/.svg/.pdf")
    log_progress(" 13. Figure5_bar_ECFP4_without_isomers.png/.tiff/.svg/.pdf")
    log_progress("  DIAGRAMMES NORMALISÉS (100%):")
    log_progress(" 14. Figure5_bar_normalized_spectral.png/.tiff/.svg/.pdf")
    log_progress(" 15. Figure5_bar_normalized_spectral_isomers.png/.tiff/.svg/.pdf")
    log_progress(" 16. Figure5_bar_normalized_ECFP4.png/.tiff/.svg/.pdf")
    log_progress(" 17. Figure5_bar_normalized_ECFP4_isomers.png/.tiff/.svg/.pdf")
    log_progress("  CSV:")
    log_progress(" 18. Figure5_spectral_candidates.csv")
    log_progress(" 19. Figure5_ECFP4_candidates.csv")