#!/usr/bin/env python
# -*- coding: utf-8 -*-
# IMPORTANT: Suppress ALL warnings FIRST before any imports
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

# Suppress RDKit logging (valence errors, etc.) BEFORE importing rdkit
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

"""
fingerprint_db_predictor.py
===========================
Prédiction des shifts NMR basée sur les Morgan Fingerprints.
Hiérarchie: Radius 4 (très exact) → Radius 3 (exact) → Radius 2 (fallback) → ML

Usage:
    from fingerprint_db_predictor import FingerprintPredictor

    predictor = FingerprintPredictor()
    result = predictor.predict_atom(mol_with_h, atom_idx, solvent, nucleus)
    # result = (shift, std, count, radius_used) ou None si pas trouvé
"""

import os
import re
import pickle
from collections import defaultdict
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import CalcNumHeteroatoms

# Beyond features pour composite fingerprint (Morgan + features chimiques + stéréochimie)
# CHANGEMENT: Utilisation de beyond_features_stereo.py au lieu de beyond_features.py
#             pour séparer automatiquement les stéréoisomères (R/S et E/Z)
try:
    from beyond_features_stereo import (
        compute_beyond_with_stereo_safe as compute_beyond_discrete_safe,
        ensure_gasteiger_charges
    )
    BEYOND_AVAILABLE = True
except ImportError:
    BEYOND_AVAILABLE = False
    def compute_beyond_discrete_safe(mol, atom_idx, include_chiral=True):
        # Retourne 14 zéros pour correspondre au nouveau format beyond_features_stereo
        # (8 features 2D + 2 features stéréo + 2 features discriminantes + 2 features distance)
        # Étape 1: +nH +hybridization | Étape 3: config relative | Étape 4: +dist_hetero +dist_arom
        return (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    def ensure_gasteiger_charges(mol):
        return False

def is_exchangeable_proton(mol, atom_idx):
    """
    Détecte si un proton est échangeable (OH, NH, SH).
    Les protons échangeables ont des shifts très variables
    selon les conditions (concentration, température, liaisons H, échange solvant)
    et doivent utiliser ML uniquement.
    """
    atom = mol.GetAtomWithIdx(atom_idx)

    # Doit être un hydrogène
    if atom.GetSymbol() != 'H':
        return False

    # Vérifier si attaché à O, N, ou S
    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() in ['O', 'N', 'S']:
            return True

    return False


def classify_carbon_radius_requirement(mol, atom_idx):
    """
    Rayon MINIMUM pour une prédiction fiable.
    Si ce rayon n'est pas disponible → retourner None → ML.

    Principe : ML > R2 faux

    Returns:
        2  → R2 suffit (CH₃, CH₂ alkyl simple, vinyl isolé)
        3  → R3 minimum (conjugué, aromatique fusionné, C=O, benzylic)
        4  → R4 minimum (CHO conjugué, jonctions polycycliques)
        99 → ML obligatoire (gem-diméthyl E/Z, through-space)
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    symbol = atom.GetSymbol()

    # H → classer par le parent
    if symbol == 'H':
        neighbors = list(atom.GetNeighbors())
        if not neighbors:
            return 2
        atom = neighbors[0]
        atom_idx = atom.GetIdx()
        symbol = atom.GetSymbol()

    # Non-carbone → R2 suffit (O, N, S rarement prédits en ¹³C)
    if symbol != 'C':
        return 2

    ring_info = mol.GetRingInfo()

    # ================================================================
    # HELPERS locaux
    # ================================================================
    def has_double_bond_to(a_idx, target_symbol=None):
        a = mol.GetAtomWithIdx(a_idx)
        for n in a.GetNeighbors():
            if target_symbol and n.GetSymbol() != target_symbol:
                continue
            bond = mol.GetBondBetweenAtoms(a_idx, n.GetIdx())
            if bond and bond.GetBondTypeAsDouble() >= 2.0:
                return True
        return False

    def count_rings(a_idx):
        return sum(1 for r in ring_info.AtomRings() if a_idx in r)

    def is_conjugated_sp2(a_idx):
        """Vérifie si C(sp2) est dans un système conjugué étendu."""
        a = mol.GetAtomWithIdx(a_idx)
        if 'SP2' not in str(a.GetHybridization()):
            return False
        for n in a.GetNeighbors():
            if n.GetSymbol() == 'C' and 'SP2' in str(n.GetHybridization()):
                return True
            if n.GetIsAromatic():
                return True
        return False

    def neighbor_symbols():
        return {n.GetSymbol() for n in atom.GetNeighbors()}

    def neighbor_has_hetero_aromatic():
        """Un voisin est aromatique ET dans un hétérocycle."""
        for n in atom.GetNeighbors():
            if n.GetIsAromatic():
                for r in ring_info.AtomRings():
                    if n.GetIdx() in r:
                        ring_atoms = [mol.GetAtomWithIdx(i) for i in r]
                        if any(a.GetSymbol() in ('N', 'O', 'S') for a in ring_atoms):
                            return True
        return False

    # FIX ERREUR 3: Compter les H directement depuis les voisins
    # GetTotalNumHs() peut être incohérent après AddHs sur certaines versions RDKit
    # Solution robuste: compter directement les voisins H
    n_H = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == 'H')
    # Degree = voisins non-H (compatible avec H explicites ou implicites)
    degree = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() != 'H')
    is_arom = atom.GetIsAromatic()
    is_in_ring = atom.IsInRing()
    nsymbols = neighbor_symbols()

    # ================================================================
    # 99 → ML OBLIGATOIRE (aucun fingerprint ne peut capturer)
    # ================================================================

    # gem-Diméthyl avec stéréochimie E/Z adjacente
    # Les deux CH₃ ont le même FP à TOUT rayon → ML seul peut
    # les distinguer via descripteurs 3D
    if n_H == 3 and degree == 1:
        heavy_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != 'H']
        if heavy_neighbors:
            parent = heavy_neighbors[0]
            if parent.GetSymbol() == 'C' and parent.GetTotalNumHs() == 0:
                # C quaternaire avec 2+ méthyles
                ch3_count = sum(
                    1 for n in parent.GetNeighbors()
                    if n.GetSymbol() == 'C' and n.GetTotalNumHs() == 3 and n.GetDegree() == 1
                )
                if ch3_count >= 2:
                    # Vérifier si C=C adjacent (E/Z qui différencie les CH₃)
                    for n in parent.GetNeighbors():
                        if n.GetSymbol() == 'C' and has_double_bond_to(parent.GetIdx(), 'C'):
                            return 99  # ML obligatoire

    # ================================================================
    # R2 SUFFIT — environnement purement local
    # ================================================================

    # CH₃ terminal simple (pas gem-diméthyl, pas adjacent C=C)
    if n_H == 3 and degree == 1:
        heavy_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != 'H']
        if heavy_neighbors:
            parent = heavy_neighbors[0]
            if parent.GetSymbol() == 'C':
                if not has_double_bond_to(parent.GetIdx()):
                    if not parent.GetIsAromatic():
                        return 2

    # CH₂ alkyl simple (pas adjacent hétéroatome, pas C=C, pas cycle, pas allylique)
    if n_H == 2 and not is_arom and not is_in_ring:
        has_hetero = any(s not in ('C', 'H') for s in nsymbols)
        has_db = any(
            mol.GetBondBetweenAtoms(atom_idx, n.GetIdx()).GetBondTypeAsDouble() >= 2.0
            for n in atom.GetNeighbors()
            if mol.GetBondBetweenAtoms(atom_idx, n.GetIdx()) is not None
        )
        # Vérifier si un voisin C a une double liaison (cas allylique)
        has_adjacent_db = any(
            has_double_bond_to(n.GetIdx(), 'C')
            for n in atom.GetNeighbors()
            if n.GetSymbol() == 'C'
        )
        if not has_hetero and not has_db and not has_adjacent_db:
            return 2

    # CH₃ sur hétéroatome (O-CH₃, N-CH₃) — shift très caractéristique
    if n_H == 3 and degree == 1:
        heavy_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != 'H']
        if heavy_neighbors:
            parent = heavy_neighbors[0]
            if parent.GetSymbol() in ('O', 'N', 'S'):
                return 2

    # CH vinyl isolé (=CH₂ terminal, pas conjugué)
    if n_H >= 1 and not is_arom:
        if has_double_bond_to(atom_idx, 'C'):
            if not is_conjugated_sp2(atom_idx):
                return 2  # Vinyl isolé

    # ================================================================
    # R3 NÉCESSAIRE
    # ================================================================

    # Tout C=O (carbonyle) : taille de cycle, conjugaison, H-bond
    if has_double_bond_to(atom_idx, 'O'):
        return 3

    # C aromatique dans système fusionné (indole, quinoline, etc.)
    if is_arom and count_rings(atom_idx) >= 2:
        return 3

    # C aromatique avec substituant mésomère (OH, NH₂, OCH₃, NO₂)
    if is_arom:
        for n in atom.GetNeighbors():
            if n.GetSymbol() in ('O', 'N') and not n.GetIsAromatic():
                return 3

    # C vinylique conjugué (diène, énone, vinyl-aromatique)
    if not is_arom and is_conjugated_sp2(atom_idx):
        return 3

    # CH₂ benzylic/allylique adjacent à hétérocycle
    if n_H >= 1 and not is_arom:
        if neighbor_has_hetero_aromatic():
            return 3

    # CH₂/CH₃ allylique (adjacent à C=C) — effet électronique du π
    if n_H >= 2 and not is_arom and not is_in_ring:
        for n in atom.GetNeighbors():
            if n.GetSymbol() == 'C' and has_double_bond_to(n.GetIdx(), 'C'):
                return 3  # Allylique

    # C dans cycle avec hétéroatome (pipéridine, pyrrolidine, DKP)
    if is_in_ring and not is_arom:
        for r in ring_info.AtomRings():
            if atom_idx in r:
                ring_atoms = [mol.GetAtomWithIdx(i) for i in r]
                if any(a.GetSymbol() in ('N', 'O', 'S') for a in ring_atoms):
                    return 3

    # C aromatique simple (benzène monosubstitué) → R3 par prudence
    # Les effets ortho/meta/para des substituants propagent à 3 liaisons
    if is_arom:
        return 3

    # ================================================================
    # R4 NÉCESSAIRE
    # ================================================================

    # C quaternaire à la jonction de 3+ cycles (spiro, bridgehead)
    if n_H == 0 and count_rings(atom_idx) >= 3:
        return 4

    # CHO conjugué à un aromatique (aldéhyde salicylique, etc.)
    if n_H == 1 and has_double_bond_to(atom_idx, 'O'):
        # C'est un aldéhyde — vérifier conjugaison avec aromatique
        for n in atom.GetNeighbors():
            if n.GetSymbol() == 'C':
                if n.GetIsAromatic():
                    return 4
                # Ou conjugué via C=C vers aromatique
                if is_conjugated_sp2(n.GetIdx()):
                    return 4

    # ================================================================
    # DÉFAUT : R2
    # ================================================================
    return 2


def parse_mol_block_robust(mol_block: str):
    """
    Parse robuste d'un bloc MOL avec valences incorrectes.
    Reconstruit la molécule manuellement à partir des atomes et bonds.
    
    Ignore le header counts (souvent incorrect ou mal formaté) et compte
    réellement les lignes d'atomes et de bonds.
    
    Returns:
        (mol_clean, num_heavy_atoms) ou (None, None) si échec
    """
    lines = mol_block.strip().split('\n')
    
    # Trouver la ligne des compteurs (header V2000)
    header_idx = None
    for i, line in enumerate(lines):
        if 'V2000' in line or 'V3000' in line:
            header_idx = i
            break
    
    if header_idx is None:
        return None, None
    
    # Créer molécule vide
    mol = Chem.RWMol()
    
    # Lire les atomes (ignorer le header, compter réellement)
    # Format SDF V2000 avec positions fixes:
    # Colonnes 1-10: x, 11-20: y, 21-30: z, 31-34: symbole
    atom_lines = []
    idx = header_idx + 1
    while idx < len(lines):
        line = lines[idx]  # NE PAS strip() pour garder les positions
        
        # Vérifier fin du bloc atomes
        if len(line.strip()) == 0 or line.strip().startswith('M '):
            break
        
        try:
            # Extraire selon positions fixes (colonnes 1-based)
            # x: colonnes 1-10, y: 11-20, z: 21-30
            if len(line) < 34:
                # Ligne trop courte, pas un atome
                break
            
            x_str = line[0:10].strip()
            y_str = line[10:20].strip()
            z_str = line[20:30].strip()
            symbol = line[30:34].strip()  # Symbole aux colonnes 31-34
            
            # Vérifier que ce sont bien des coordonnées
            float(x_str)
            float(y_str)
            float(z_str)
            
            # Vérifier que le symbole est valide (1-2 caractères)
            if len(symbol) > 0 and len(symbol) <= 3 and symbol[0].isalpha():
                atom_lines.append((symbol, line))
                idx += 1
            else:
                break
        except:
            # Pas une ligne d'atome
            break
    
    # Créer les atomes
    for symbol, line in atom_lines:
        try:
            atom = Chem.Atom(symbol)
            mol.AddAtom(atom)
        except:
            return None, None
    
    n_atoms = len(atom_lines)
    if n_atoms == 0:
        return None, None
    
    # Lire les bonds (jusqu'à M  END ou M CHG)
    bond_start = header_idx + 1 + n_atoms
    for i in range(bond_start, len(lines)):
        line = lines[i].strip()
        if line.startswith('M '):
            break
        
        try:
            parts = line.split()
            if len(parts) < 3:
                continue
            
            atom1 = int(parts[0]) - 1  # SDF est 1-based
            atom2 = int(parts[1]) - 1
            bond_type = int(parts[2])
            
            # Vérifier que les indices sont valides
            if atom1 < 0 or atom1 >= n_atoms or atom2 < 0 or atom2 >= n_atoms:
                continue
            
            if bond_type == 1:
                bt = Chem.BondType.SINGLE
            elif bond_type == 2:
                bt = Chem.BondType.DOUBLE
            elif bond_type == 3:
                bt = Chem.BondType.TRIPLE
            else:
                bt = Chem.BondType.SINGLE
            
            mol.AddBond(atom1, atom2, bt)
        except:
            pass  # Ignorer les bonds invalides
    
    # Convertir en molécule standard
    try:
        mol_clean = mol.GetMol()
        # Sanitizer pour calculer aromaticité, hybridation, etc.
        Chem.SanitizeMol(mol_clean)
        return mol_clean, n_atoms
    except:
        # Si sanitization échoue, retourner la molécule brute
        try:
            mol_clean = mol.GetMol()
            # Tentative partielle
            mol_clean.UpdatePropertyCache(strict=False)
            return mol_clean, n_atoms
        except:
            return None, None


# =============================================================================
# CONFIGURATION - Chemins relatifs au script
# =============================================================================
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_DATA_DIR = os.path.join(_THIS_DIR, '..', 'ml_dependencies', 'data')

BMRB_DIR = os.path.join(_DATA_DIR, 'sdf_bmrb_rdkit_v2')
NMRSHIFTDB_FILE = os.path.join(_DATA_DIR, 'nmrshiftdb2withsignals.sd')
NMR_SPECTRA_FILE = os.path.join(_DATA_DIR, 'all_nmr_spectra_nmrshiftdb_format.sdf')
NPMRD_SDF_DIR = os.path.join(_DATA_DIR, 'NPMRD')  # Répertoire avec 1607 fichiers SDF NPMRD
EXPERIMENTAL_DATA_DIR = os.path.join(_DATA_DIR, 'with_experimental_data')
DB_CACHE_FILE = os.path.join(_DATA_DIR, 'fingerprint_database.pkl')

# Radius supportés (6=ultra précis → 2=fallback)
RADIUS_ULTRA_EXACT = 6   # Environnement ultra-précis (6 liaisons)
RADIUS_HYPER_EXACT = 5   # Environnement hyper-précis (5 liaisons)
RADIUS_VERY_EXACT = 4    # Environnement très précis (4 liaisons)
RADIUS_EXACT = 3         # Environnement exact (3 liaisons)
RADIUS_FALLBACK = 2      # Fallback (2 liaisons)
SUPPORTED_RADII = [RADIUS_ULTRA_EXACT, RADIUS_HYPER_EXACT, RADIUS_VERY_EXACT, RADIUS_EXACT, RADIUS_FALLBACK]  # Du plus précis au moins précis

# SubstructMatchParameters avec maxMatches=1
from rdkit.Chem import SubstructMatchParameters as _SMP
_SUBSTRUCT_PARAMS = _SMP()
_SUBSTRUCT_PARAMS.maxMatches = 1

def _substruct_match_with_timeout(mol, mol_canon, timeout=1.0):
    """GetSubstructMatch avec timeout (évite les boucles infinies sur macrocycles)."""
    import threading
    result = []
    def _run():
        try:
            result.append(mol.GetSubstructMatch(mol_canon, _SUBSTRUCT_PARAMS))
        except Exception:
            pass
    t = threading.Thread(target=_run, daemon=True)
    t.start()
    t.join(timeout)
    if t.is_alive():
        # Timeout — on retourne un tuple vide, la molécule ira en fingerprint uniquement
        return ()
    return result[0] if result else ()

# FIX ERREUR 8: Seuils de qualité ajustés (0.2 ppm était trop strict pour 1H)
# Les shifts 1H varient naturellement de 0.1-0.3 ppm entre solvants/conditions
MAX_STD_13C = 3.0   # ppm - si std > 3.0 ppm pour 13C, utiliser ML (variabilité naturelle 1-3 ppm)
MAX_STD_1H = 0.35   # ppm - si std > 0.35 ppm pour 1H, utiliser ML (plus réaliste pour variabilité)

# Seuils FALLBACK (accepter même avec std élevée si pas d'autre choix)
FALLBACK_MAX_STD_13C = 5.0   # ppm - seuil absolu pour 13C
FALLBACK_MAX_STD_1H = 0.6    # ppm - seuil absolu pour 1H

# Seuils minimum de sources par radius
# Principe: Plus le radius est élevé, plus l'environnement est spécifique
# → Augmenté à 5 pour forcer fallback quand données insuffisantes
MIN_SAMPLES_R6 = 1       # R6 ultra-spécifique, 1 sample suffit
MIN_SAMPLES_R5 = 1       # R5 hyper-spécifique, 1 sample suffit
MIN_SAMPLES_R4 = 1       # R4 très spécifique, 1 sample suffit
MIN_SAMPLES_R3 = 3       # R3: besoin de 3 pour validation
MIN_SAMPLES_R2 = 3       # R2 fallback: besoin de 3 samples

# Logging
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)  # DEBUG pour développement, WARNING pour production

# =============================================================================
# NORMALISATION DES SOLVANTS
# =============================================================================

# Map pour normaliser les variantes de solvants
SOLVENT_MAP = {
    # CDCl3 variants
    'cdcl3': 'CDCl3',
    'chloroform': 'CDCl3',
    'chloroform-d': 'CDCl3',
    'chloroform-d1': 'CDCl3',
    'deuterochloroform': 'CDCl3',
    'chloroform-d1 (cdcl3)': 'CDCl3',
    'chcl3': 'CDCl3',
    'cdcl3,tms': 'CDCl3',

    # D2O variants
    'd2o': 'D2O',
    'deuterium oxide': 'D2O',
    'water': 'D2O',
    'h2o': 'D2O',
    'deuteriumoxide (d2o)': 'D2O',
    'water (d2o)': 'D2O',
    'water (h2o)': 'D2O',

    # DMSO-d6 variants
    'dmso': 'DMSO-d6',
    'dmso-d6': 'DMSO-d6',
    'deuterated dmso': 'DMSO-d6',
    'dimethylsulphoxide-d6 (dmso-d6, c2d6os)': 'DMSO-d6',
    'dimethyl sulfoxide-d6': 'DMSO-d6',
    'dimco-d6': 'DMSO-d6',
    '(cd3)2so': 'DMSO-d6',
    'dmco-d6': 'DMSO-d6',
    'dimethylsulphoxide-d6 (dmso-d6, c2d6so)': 'DMSO-d6',

    # MeOD/Methanol variants
    'meod': 'MeOD',
    'meod-d4': 'MeOD',
    'methanol-d4': 'MeOD',
    'cd3od': 'MeOD',
    'meoh-d4': 'MeOD',
    'methanol-d4 (cd3od)': 'MeOD',
    'methanol-d4 (meod-d4)': 'MeOD',  # Variant trouvée dans nmrshiftdb2
    'methanol': 'MeOD',
    'meoh': 'MeOD',
    'methanol-d1 (ch3od)': 'MeOD',
    'methanol-d3(cd3oh)': 'MeOD',
    'methanol-d4/d2o (cd3od/d4)': 'MeOD',

    # Acetone-d6 variants
    'acetone-d6': 'Acetone-d6',
    'acetone': 'Acetone-d6',
    'acetone-d6 ((cd3)2co)': 'Acetone-d6',
    'cd3cocd': 'Acetone-d6',
    'acetone-d6 + 1 drop d2o': 'Acetone-d6',

    # Pyridine-d5 variants
    'pyridine-d5': 'Pyridine-d5',
    'pyridin-d5 (c5d5n)': 'Pyridine-d5',

    # Benzene-d6 variants
    'benzene-d6': 'Benzene-d6',
    'c6d6': 'Benzene-d6',
    'benzene-d6 (c6d6)': 'Benzene-d6',
    'benzene': 'Benzene-d6',

    # Other common solvents
    'cd2cl2': 'CD2Cl2',
    'methylenchloride-d2 (cd2cl2)': 'CD2Cl2',
    'thf-d8': 'THF-d8',
    'tetrahydrofuran-d8 (thf-d8, c4d4o)': 'THF-d8',
    'toluene-d8': 'Toluene-d8',
    'toluol-d8': 'Toluene-d8',
    'toluene': 'Toluene-d8',
    'cd3cn': 'CD3CN',
    'acetonitrile-d3(cd3cn)': 'CD3CN',
    'cdcn': 'CD3CN',
    'cd3cn or cdcl3': 'CD3CN',
    'ccl4': 'CCl4',
    'tetrachloro-methane (ccl4)': 'CCl4',
    'ethanol': 'Ethanol',
}

def normalize_solvent(raw: str) -> str:
    """Normalise le nom du solvant."""
    if not raw:
        return 'Unreported'

    solvent = raw.strip()
    if ':' in solvent and solvent.split(':')[0].isdigit():
        solvent = solvent.split(':', 1)[1].strip()

    if ';' in solvent:
        solvent = solvent.split(';')[0].strip()

    key = solvent.lower().strip()

    # Unifier tous les solvants inconnus vers "Unreported"
    unknown_variants = ['unknown', 'unreported', 'not reported', 'not specified',
                        'n/a', 'na', 'none', 'undefined', '']
    if key in unknown_variants:
        return 'Unreported'

    return SOLVENT_MAP.get(key, solvent)


def is_valid_solvent(solvent: str) -> bool:
    """
    Vérifie si un solvant est valide.
    Maintenant accepte "Unreported" pour inclure ces données dans les index fingerprint.
    """
    if not solvent:
        return False

    # Tous les solvants sont valides (y compris "Unreported")
    # Les variantes sont normalisées vers "Unreported" dans normalize_solvent()
    return True


# Reverse mapping: solvant normalisé → variantes connues non-normalisées
# Utilisé pour fusionner les données quand la DB contient des variantes non-normalisées
SOLVENT_VARIANTS = {}
for raw_key, normalized in SOLVENT_MAP.items():
    if normalized not in SOLVENT_VARIANTS:
        SOLVENT_VARIANTS[normalized] = set()
    # Ajouter la clé originale en capitalisation originale (si différente)
    if raw_key != normalized.lower():
        # On stocke la version titre/originale
        SOLVENT_VARIANTS[normalized].add(raw_key.title().replace('-D', '-d').replace(' D', ' d'))

# Ajouter manuellement les variantes exactes trouvées dans la DB
SOLVENT_VARIANTS.setdefault('MeOD', set()).add('Methanol-D4 (MeOD-d4)')


def get_solvent_variants(normalized_solvent: str) -> list:
    """
    Retourne toutes les variantes connues d'un solvant normalisé.
    Utile pour chercher dans une DB qui contient des variantes non-normalisées.
    """
    variants = [normalized_solvent]  # Toujours inclure le normalisé
    if normalized_solvent in SOLVENT_VARIANTS:
        variants.extend(SOLVENT_VARIANTS[normalized_solvent])
    return variants


# =============================================================================
# FINGERPRINT CALCULATION
# =============================================================================
def get_atom_fingerprint_raw(mol, atom_idx: int, radius: int = RADIUS_VERY_EXACT, heavy_only: bool = False) -> int:
    """
    Calcule le Morgan fingerprint brut pour un atome.
    
    Args:
        mol: Molécule RDKit
        atom_idx: Index de l'atome
        radius: Radius Morgan (4=très exact, 3=exact, 2=fallback)
        heavy_only: Si True, calcule le FP sur une molécule sans H explicites
                   (capture mieux la topologie des atomes lourds)
    """
    # Si heavy_only, travailler sur une copie sans H
    if heavy_only:
        mol_for_fp = Chem.RemoveHs(mol)
        # Retrouver l'index correspondant dans la molécule sans H
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() == 'H':
            # Pour un H, utiliser son parent lourd
            neighbors = list(atom.GetNeighbors())
            if neighbors:
                parent_idx = neighbors[0].GetIdx()
                # Compter combien d'atomes lourds avant le parent
                heavy_count = 0
                for i in range(parent_idx):
                    if mol.GetAtomWithIdx(i).GetSymbol() != 'H':
                        heavy_count += 1
                atom_idx_in_no_h = heavy_count
            else:
                return 0
        else:
            # Pour un atome lourd, compter combien d'atomes lourds avant lui
            heavy_count = 0
            for i in range(atom_idx):
                if mol.GetAtomWithIdx(i).GetSymbol() != 'H':
                    heavy_count += 1
            atom_idx_in_no_h = heavy_count

        # FIX ERREUR 5: Vérification de bornes pour éviter crashes
        if atom_idx_in_no_h >= mol_for_fp.GetNumAtoms():
            return 0
    else:
        mol_for_fp = mol
        atom_idx_in_no_h = atom_idx

    info = {}
    AllChem.GetMorganFingerprint(mol_for_fp, radius, bitInfo=info)

    # FIX ERREUR 1-2: Collecter (bit_id, radius_du_bit) pour distinguer
    # les mêmes bits à différentes profondeurs et éviter le XOR destructif
    atom_bits_with_radius = []
    for bit_id, atom_list in info.items():
        for a_idx, rad in atom_list:
            if a_idx == atom_idx_in_no_h:
                # Garder le bit_id ET le radius où il a été trouvé
                atom_bits_with_radius.append((bit_id, rad))

    if not atom_bits_with_radius:
        return 0

    # FIX: Utiliser hash(tuple(sorted(...))) au lieu de XOR
    # Cela préserve l'information d'ordre et évite les annulations
    return hash(tuple(sorted(atom_bits_with_radius)))


def get_atom_fingerprint(mol, atom_idx: int, radius: int = RADIUS_VERY_EXACT, heavy_only: bool = False) -> int:
    """
    Calcule le Morgan fingerprint pour un atome.
    
    STRATÉGIE PAR DÉFAUT (heavy_only=False):
    =========================================
    Inclut les H explicites dans le calcul du fingerprint. Les H contribuent
    à distinguer les environnements (ex: CH3 vs CH2 vs CH).
    
    Tests sur la caféine (3 N-CH3):
      - heavy_only=False: 2/3 N-CH3 distingués (C0 ≠ C9=C13) même avec R6
      - heavy_only=True: 0/3 distingués (tous identiques)
      → heavy_only=False est MEILLEUR
    
    Limitation: Certaines molécules symétriques (caféine, théobromine) ont
    des groupes intrinsèquement indistinguables par fingerprint.
    → Solution: détection d'ambiguïté dans integrated_nmr_prediction.py
    
    Args:
        mol: Molécule RDKit (avec H explicites)
        atom_idx: Index de l'atome
        radius: Radius Morgan (4=très exact, 3=exact, 2=fallback)
        heavy_only: Si True, exclut les H du calcul (généralement moins bon)
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    
    if atom.GetSymbol() == 'H':
        neighbors = list(atom.GetNeighbors())
        if neighbors:
            parent_idx = neighbors[0].GetIdx()
            # Calculer le FP du PARENT (C, N, O, etc.) avec ou sans H
            parent_fp = get_atom_fingerprint_raw(mol, parent_idx, radius, heavy_only=heavy_only)
            # Le FP du H lui-même (toujours avec heavy_only pour cohérence)
            h_fp = get_atom_fingerprint_raw(mol, atom_idx, radius, heavy_only=heavy_only)
            return (radius << 62) | (parent_fp << 32) | h_fp
    
    # Pour les atomes lourds (C, N, O, etc.)
    raw_fp = get_atom_fingerprint_raw(mol, atom_idx, radius, heavy_only=heavy_only)
    return (radius << 62) | raw_fp


def get_atom_env_string(mol, atom_idx: int) -> str:
    """Génère une description lisible de l'environnement atomique avec descripteurs 2D."""
    atom = mol.GetAtomWithIdx(atom_idx)
    env_parts = [atom.GetSymbol()]
    hyb = str(atom.GetHybridization()).replace('SP', 'sp')
    env_parts.append(hyb)
    env_parts.append(f"H{atom.GetTotalNumHs()}")
    neighbors = []
    for n in atom.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(atom_idx, n.GetIdx())
        bond_type = str(bond.GetBondType())[0]
        neighbors.append(f"{n.GetSymbol()}{bond_type}")
    neighbors.sort()
    env_parts.append("_".join(neighbors))
    
    # Ajouter descripteurs 2D pour effets conformationnels
    fast_2d = get_fast_2d_descriptors(mol, atom_idx)
    env_parts.append(fast_2d)
    
    return "|".join(env_parts)


def get_fast_2d_descriptors(mol, atom_idx: int) -> str:
    """
    Descripteurs 2D ultra-rapides pour capturer ~80% des effets conformationnels
    sans calcul 3D (temps: ~0.1 ms par atome).
    
    Inspiré de l'analyse: ces descripteurs capturent les effets stériques/anisotropiques
    majeurs sans la complexité temporelle de la conformation 3D.
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    descriptors = []
    
    # 1. Effet de cycle (anisotropie, contrainte conformationnelle)
    ring_info = "acyclic"
    if atom.IsInRing():
        ring_info_obj = mol.GetRingInfo()
        # Trouver la taille du plus petit cycle
        min_ring_size = 999
        for ring_size in [3, 4, 5, 6, 7, 8, 9, 10]:  # Tailles communes
            if ring_info_obj.IsAtomInRingOfSize(atom_idx, ring_size):
                min_ring_size = min(min_ring_size, ring_size)
        if min_ring_size < 999:
            ring_info = f"ring{min_ring_size}"
        else:
            ring_info = "ring"
    descriptors.append(ring_info)
    
    # 2. Groupes anisotropes à proximité (effet longue distance)
    aniso_groups = 0
    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() in ['O', 'N'] and neighbor.GetTotalNumHs() == 0:
            aniso_groups += 1  # C=O, C=N, etc.
    descriptors.append(f"aniso{aniso_groups}")
    
    # 3. Encombrement stérique (compression van der Waals)
    heavy_neighbors = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() > 1)
    descriptors.append(f"heavy{heavy_neighbors}")
    
    return "|".join(descriptors)


# =============================================================================
# PARALLEL PROCESSING WORKER FUNCTION
# =============================================================================
def _process_mol_block_worker(args):
    """
    Worker function pour traitement parallèle des molécules.
    Doit être au niveau module pour être picklable par multiprocessing.

    Returns:
        Tuple (fp_records, exact_records, composite_records, count) ou None si erreur
        - fp_records: list of (radius, fp, solvent, nucleus, shift, env_str)
        - exact_records: list of (inchikey, solvent, nucleus, canon_idx, shift)
        - composite_records: list of (radius, composite_fp, solvent, nucleus, shift)
        - count: nombre de shifts traités
    """
    from rdkit.Chem import inchi as rdkit_inchi

    mol_block, mol_count = args

    try:
        mol_name = mol_block.split('\n')[0].strip() or f'mol#{mol_count}'

        m_end_idx = mol_block.find('M  END')
        if m_end_idx < 0:
            return None

        mol_text = mol_block[:m_end_idx + 6]

        # Parser robuste
        mol, num_heavy = parse_mol_block_robust(mol_text)
        if mol is None:
            return None

        mol_with_h = Chem.AddHs(mol)

        # Calculer charges Gasteiger pour beyond features
        beyond_available = False
        try:
            if BEYOND_AVAILABLE:
                ensure_gasteiger_charges(mol_with_h)
                beyond_available = True
        except:
            pass

        # Calculer InChIKey et canon_map
        inchikey = None
        canon_map = {}
        try:
            inchikey = rdkit_inchi.MolToInchiKey(mol)
            if inchikey:
                canon_smiles = Chem.MolToSmiles(mol)
                mol_canon_no_h = Chem.MolFromSmiles(canon_smiles)
                if mol_canon_no_h is not None:
                    match_no_h = _substruct_match_with_timeout(mol, mol_canon_no_h, timeout=1.0)
                    for canon_idx, sdf_idx in enumerate(match_no_h):
                        canon_map[sdf_idx] = canon_idx

                    mol_canon_with_h = Chem.AddHs(mol_canon_no_h)
                    match_with_h = _substruct_match_with_timeout(mol_with_h, mol_canon_with_h, timeout=1.0)
                    if match_with_h:
                        for canon_idx, sdf_idx in enumerate(match_with_h):
                            canon_map[sdf_idx] = canon_idx
        except Exception:
            inchikey = None

        # Extraire solvant
        solvent_match = re.search(r'<Solvent>\s*\n([^\n]+)', mol_block)
        solvent_map = {}
        if solvent_match:
            solvent_line = solvent_match.group(1).strip()
            parts = re.findall(r'(\d+):(.+?)(?=\s+\d+:|$)', solvent_line)
            for idx, solv in parts:
                solvent_map[idx] = normalize_solvent(solv.strip())

        fp_records = []
        exact_records = []
        composite_records = []  # NOUVEAU: fingerprints composites
        count = 0

        # Parse 13C and 1H spectra
        for nucleus, pattern in [('13C', r'<Spectrum 13C (\d+)>\s*\n([^\n<]+)'),
                                 ('1H', r'<Spectrum 1H (\d+)>\s*\n([^\n<]+)')]:
            for spec_match in re.finditer(pattern, mol_block):
                spec_idx = spec_match.group(1)
                spec_data = spec_match.group(2).strip()
                solvent = solvent_map.get(spec_idx, 'Unknown')
                solvent_normalized = normalize_solvent(solvent)

                if not is_valid_solvent(solvent_normalized):
                    continue

                for entry in spec_data.split('|'):
                    parts = entry.split(';')
                    if len(parts) >= 3:
                        try:
                            shift = float(parts[0])
                            atom_idx_no_h = int(parts[2])

                            if atom_idx_no_h >= mol.GetNumAtoms():
                                continue

                            # Filtrage outliers
                            if nucleus == '1H':
                                if shift < -5.0 or shift > 14.0:
                                    continue
                            elif nucleus == '13C':
                                if shift < -5.0 or shift > 250.0:
                                    continue

                            if nucleus == '1H':
                                # Trouver les H attachés
                                atom = mol_with_h.GetAtomWithIdx(atom_idx_no_h)
                                if atom.GetSymbol() == 'H':
                                    h_indices = [atom_idx_no_h]
                                else:
                                    h_indices = [n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == 'H']

                                for h_idx in h_indices:
                                    env_str = get_atom_env_string(mol_with_h, h_idx)
                                    for radius in SUPPORTED_RADII:
                                        fp = get_atom_fingerprint(mol_with_h, h_idx, radius)
                                        fp_records.append((radius, fp, solvent_normalized, nucleus, shift, env_str))

                                        # COMPOSITE: Ajouter fingerprint composite si beyond disponible
                                        if beyond_available:
                                            try:
                                                beyond = compute_beyond_discrete_safe(mol_with_h, h_idx, include_chiral=False)
                                                composite_fp = (fp, beyond)  # Stocker le tuple directement (pas hash)
                                                composite_records.append((radius, composite_fp, solvent_normalized, nucleus, shift))
                                            except Exception as e:
                                                # Ne pas ignorer silencieusement - logger l'erreur
                                                logger.debug(f"Beyond features failed for H{h_idx}: {e}")
                                                pass

                                    if inchikey and h_idx in canon_map:
                                        exact_records.append((inchikey, solvent_normalized, nucleus, canon_map[h_idx], shift))
                                    count += 1
                            else:
                                # 13C
                                env_str = get_atom_env_string(mol_with_h, atom_idx_no_h)
                                for radius in SUPPORTED_RADII:
                                    fp = get_atom_fingerprint(mol_with_h, atom_idx_no_h, radius)
                                    fp_records.append((radius, fp, solvent_normalized, nucleus, shift, env_str))

                                    # COMPOSITE: Ajouter fingerprint composite si beyond disponible
                                    if beyond_available:
                                        try:
                                            beyond = compute_beyond_discrete_safe(mol_with_h, atom_idx_no_h, include_chiral=False)
                                            composite_fp = (fp, beyond)
                                            composite_records.append((radius, composite_fp, solvent_normalized, nucleus, shift))
                                        except:
                                            pass

                                if inchikey and atom_idx_no_h in canon_map:
                                    exact_records.append((inchikey, solvent_normalized, nucleus, canon_map[atom_idx_no_h], shift))
                                count += 1
                        except (ValueError, IndexError):
                            pass

        return (fp_records, exact_records, composite_records, count)

    except Exception:
        return None


def _process_nmr_spectra_worker(args):
    """Worker pour NMR Spectra format (légèrement différent de NMRShiftDB).
    Retourne (fp_records, exact_records, composite_records, count)."""
    from rdkit.Chem import inchi as rdkit_inchi

    mol_block, mol_count = args

    try:
        m_end_idx = mol_block.find('M  END')
        if m_end_idx < 0:
            return None

        mol_text = mol_block[:m_end_idx + 6]
        mol, num_heavy = parse_mol_block_robust(mol_text)
        if mol is None:
            return None

        mol_with_h = Chem.AddHs(mol)

        # Calculer charges Gasteiger pour beyond features
        beyond_available = False
        try:
            if BEYOND_AVAILABLE:
                ensure_gasteiger_charges(mol_with_h)
                beyond_available = True
        except:
            pass

        # Calculer InChIKey et canon_map
        inchikey = None
        canon_map = {}
        try:
            inchikey = rdkit_inchi.MolToInchiKey(mol)
            if inchikey:
                canon_smiles = Chem.MolToSmiles(mol)
                mol_canon_no_h = Chem.MolFromSmiles(canon_smiles)
                if mol_canon_no_h is not None:
                    match_no_h = _substruct_match_with_timeout(mol, mol_canon_no_h, timeout=1.0)
                    for canon_idx, sdf_idx in enumerate(match_no_h):
                        canon_map[sdf_idx] = canon_idx
                    mol_canon_with_h = Chem.AddHs(mol_canon_no_h)
                    match_with_h = _substruct_match_with_timeout(mol_with_h, mol_canon_with_h, timeout=1.0)
                    if match_with_h:
                        for canon_idx, sdf_idx in enumerate(match_with_h):
                            canon_map[sdf_idx] = canon_idx
        except Exception:
            inchikey = None

        # Format NMR Spectra: <Solvent> (1)\n...
        solvent_match = re.search(r'<Solvent>\s*\(\d+\)\s*\n([^\n]+)', mol_block)
        solvent_map = {}
        if solvent_match:
            solvent_line = solvent_match.group(1).strip()
            parts = re.findall(r'(\d+):(.+?)(?=\s*\d+:|$)', solvent_line)
            for idx, solv in parts:
                solvent_map[idx] = normalize_solvent(solv.strip())
            if not solvent_map:
                solvent_map['0'] = normalize_solvent(solvent_line)

        fp_records = []
        exact_records = []
        composite_records = []  # NOUVEAU: fingerprints composites
        count = 0

        # Format NMR Spectra: <Spectrum 13C 0> (1)\n...
        for nucleus, pattern in [('13C', r'<Spectrum 13C (\d+)>\s*\(\d+\)\s*\n([^\n<]+)'),
                                 ('1H', r'<Spectrum 1H (\d+)>\s*\(\d+\)\s*\n([^\n<]+)')]:
            for spec_match in re.finditer(pattern, mol_block):
                spec_idx = spec_match.group(1)
                spec_data = spec_match.group(2).strip()
                solvent = solvent_map.get(spec_idx, solvent_map.get('0', 'Unknown'))
                solvent_normalized = normalize_solvent(solvent)

                if not is_valid_solvent(solvent_normalized):
                    continue

                for entry in spec_data.split('|'):
                    if not entry.strip():
                        continue
                    parts = entry.split(';')
                    if len(parts) >= 3:
                        try:
                            shift = float(parts[0])
                            atom_idx = int(parts[2])

                            if atom_idx >= mol_with_h.GetNumAtoms():
                                continue

                            # Filtrage outliers
                            if nucleus == '1H' and (shift < -5.0 or shift > 14.0):
                                continue
                            if nucleus == '13C' and (shift < -5.0 or shift > 250.0):
                                continue

                            env_str = get_atom_env_string(mol_with_h, atom_idx)
                            for radius in SUPPORTED_RADII:
                                fp = get_atom_fingerprint(mol_with_h, atom_idx, radius)
                                fp_records.append((radius, fp, solvent_normalized, nucleus, shift, env_str))

                                # COMPOSITE: Ajouter fingerprint composite si beyond disponible
                                if beyond_available:
                                    try:
                                        beyond = compute_beyond_discrete_safe(mol_with_h, atom_idx, include_chiral=False)
                                        composite_fp = (fp, beyond)
                                        composite_records.append((radius, composite_fp, solvent_normalized, nucleus, shift))
                                    except:
                                        pass

                            if inchikey and atom_idx in canon_map:
                                exact_records.append((inchikey, solvent_normalized, nucleus, canon_map[atom_idx], shift))
                            count += 1
                        except (ValueError, IndexError):
                            pass

        return (fp_records, exact_records, composite_records, count)

    except Exception:
        return None


def _process_npmrd_file_worker(args):
    """Worker pour fichiers NPMRD individuels.
    Retourne (fp_records, exact_records, composite_records, count)."""
    from rdkit.Chem import inchi as rdkit_inchi

    filepath, file_idx = args

    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        m_end_idx = content.find('M  END')
        if m_end_idx < 0:
            return None

        mol_text = content[:m_end_idx + 6]
        mol, num_heavy = parse_mol_block_robust(mol_text)
        if mol is None:
            return None

        mol_with_h = Chem.AddHs(mol)

        # Calculer charges Gasteiger pour beyond features
        beyond_available = False
        try:
            if BEYOND_AVAILABLE:
                ensure_gasteiger_charges(mol_with_h)
                beyond_available = True
        except:
            pass

        # Calculer InChIKey et canon_map
        inchikey = None
        canon_map = {}
        try:
            inchikey = rdkit_inchi.MolToInchiKey(mol)
            if inchikey:
                canon_smiles = Chem.MolToSmiles(mol)
                mol_canon_no_h = Chem.MolFromSmiles(canon_smiles)
                if mol_canon_no_h is not None:
                    match_no_h = _substruct_match_with_timeout(mol, mol_canon_no_h, timeout=1.0)
                    for canon_idx, sdf_idx in enumerate(match_no_h):
                        canon_map[sdf_idx] = canon_idx
                    mol_canon_with_h = Chem.AddHs(mol_canon_no_h)
                    match_with_h = _substruct_match_with_timeout(mol_with_h, mol_canon_with_h, timeout=1.0)
                    if match_with_h:
                        for canon_idx, sdf_idx in enumerate(match_with_h):
                            canon_map[sdf_idx] = canon_idx
        except Exception:
            inchikey = None

        # Parser le solvant
        solvent_match = re.search(r'> <Solvent>\s*\n([^\n]+)', content)
        solvent = normalize_solvent(solvent_match.group(1).strip()) if solvent_match else 'Unknown'

        if not is_valid_solvent(solvent):
            solvent = 'Unknown'

        fp_records = []
        exact_records = []
        composite_records = []  # NOUVEAU: fingerprints composites
        count = 0

        # Parser spectre 13C
        spec_13c_match = re.search(r'> <Spectrum_13C>\s*\n([^\n]+)', content)
        if spec_13c_match:
            spec_data = spec_13c_match.group(1).strip()
            for entry in spec_data.split('|'):
                if not entry.strip():
                    continue
                parts = entry.split(';')
                if len(parts) >= 3:
                    try:
                        shift = float(parts[0])
                        # NPMRD: indices sont 1-based
                        atom_idx = int(parts[2]) - 1
                        if atom_idx < 0 or atom_idx >= mol_with_h.GetNumAtoms():
                            continue
                        if shift < -5.0 or shift > 250.0:
                            continue

                        env_str = get_atom_env_string(mol_with_h, atom_idx)
                        for radius in SUPPORTED_RADII:
                            fp = get_atom_fingerprint(mol_with_h, atom_idx, radius)
                            fp_records.append((radius, fp, solvent, '13C', shift, env_str))

                            # COMPOSITE: Ajouter fingerprint composite
                            if beyond_available:
                                try:
                                    beyond = compute_beyond_discrete_safe(mol_with_h, atom_idx, include_chiral=False)
                                    composite_fp = (fp, beyond)
                                    composite_records.append((radius, composite_fp, solvent, '13C', shift))
                                except:
                                    pass

                        if inchikey and atom_idx in canon_map:
                            exact_records.append((inchikey, solvent, '13C', canon_map[atom_idx], shift))
                        count += 1
                    except (ValueError, IndexError):
                        pass

        # Parser spectre 1H
        spec_1h_match = re.search(r'> <Spectrum_1H>\s*\n([^\n]+)', content)
        if not spec_1h_match:
            spec_1h_match = re.search(r'> <Spectrum_Unknown>\s*\n([^\n]+)', content)

        if spec_1h_match:
            spec_data = spec_1h_match.group(1).strip()
            for entry in spec_data.split('|'):
                if not entry.strip():
                    continue
                parts = entry.split(';')
                if len(parts) >= 3:
                    try:
                        shift = float(parts[0])
                        # NPMRD 1H: indices multiples possibles (parts[2:])
                        for idx_str in parts[2:]:
                            atom_idx = int(idx_str) - 1  # 1-based
                            if atom_idx < 0 or atom_idx >= mol_with_h.GetNumAtoms():
                                continue
                            if shift < -5.0 or shift > 14.0:
                                continue

                            env_str = get_atom_env_string(mol_with_h, atom_idx)
                            for radius in SUPPORTED_RADII:
                                fp = get_atom_fingerprint(mol_with_h, atom_idx, radius)
                                fp_records.append((radius, fp, solvent, '1H', shift, env_str))

                                # COMPOSITE: Ajouter fingerprint composite pour 1H
                                if beyond_available:
                                    try:
                                        beyond = compute_beyond_discrete_safe(mol_with_h, atom_idx, include_chiral=False)
                                        composite_fp = (fp, beyond)
                                        composite_records.append((radius, composite_fp, solvent, '1H', shift))
                                    except:
                                        pass

                            if inchikey and atom_idx in canon_map:
                                exact_records.append((inchikey, solvent, '1H', canon_map[atom_idx], shift))
                            count += 1
                    except (ValueError, IndexError):
                        pass

        return (fp_records, exact_records, composite_records, count)

    except Exception:
        return None


def _process_experimental_file_worker(args):
    """Worker pour fichiers Experimental individuels.
    Retourne (fp_records, exact_records, composite_records, count)."""
    from rdkit.Chem import inchi as rdkit_inchi

    filepath, file_idx = args

    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        m_end_idx = content.find('M  END')
        if m_end_idx < 0:
            return None

        mol_text = content[:m_end_idx + 6]
        mol = Chem.MolFromMolBlock(mol_text)
        if mol is None:
            return None

        num_heavy = mol.GetNumAtoms()
        mol_with_h = Chem.AddHs(mol)

        # Calculer charges Gasteiger pour beyond features
        beyond_available = False
        try:
            if BEYOND_AVAILABLE:
                ensure_gasteiger_charges(mol_with_h)
                beyond_available = True
        except:
            pass

        # Calculer InChIKey et canon_map
        inchikey = None
        canon_map = {}
        try:
            inchikey = rdkit_inchi.MolToInchiKey(mol)
            if inchikey:
                canon_smiles = Chem.MolToSmiles(mol)
                mol_canon_no_h = Chem.MolFromSmiles(canon_smiles)
                if mol_canon_no_h is not None:
                    match_no_h = _substruct_match_with_timeout(mol, mol_canon_no_h, timeout=1.0)
                    for canon_idx, sdf_idx in enumerate(match_no_h):
                        canon_map[sdf_idx] = canon_idx
                    mol_canon_with_h = Chem.AddHs(mol_canon_no_h)
                    match_with_h = _substruct_match_with_timeout(mol_with_h, mol_canon_with_h, timeout=1.0)
                    if match_with_h:
                        for canon_idx, sdf_idx in enumerate(match_with_h):
                            canon_map[sdf_idx] = canon_idx
        except Exception:
            inchikey = None

        # Parser solvant
        solvent_match = re.search(r'<Solvent>\s*\(\d+\)\s*\n([^\n]+)', content)
        solvent = normalize_solvent(solvent_match.group(1).strip()) if solvent_match else 'Unknown'

        fp_records = []
        exact_records = []
        composite_records = []  # NOUVEAU: fingerprints composites
        count = 0

        # Parser spectre 1H
        spec_match = re.search(r'<Spectrum 1H>\s*\(\d+\)\s*\n([^\n]+)', content)
        if spec_match:
            spec_data = spec_match.group(1).strip()
            for entry in spec_data.split('|'):
                if not entry.strip():
                    continue
                parts = entry.split(';')
                if len(parts) >= 3:
                    try:
                        shift = float(parts[0])
                        atom_idx = int(parts[2])
                        if atom_idx >= mol_with_h.GetNumAtoms():
                            continue
                        if shift < -5.0 or shift > 14.0:
                            continue

                        # Gérer Format A (atome lourd) et Format B (H direct)
                        atom = mol_with_h.GetAtomWithIdx(atom_idx)
                        if atom.GetSymbol() == 'H':
                            h_indices = [atom_idx]
                        else:
                            h_indices = [n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == 'H']

                        for h_idx in h_indices:
                            env_str = get_atom_env_string(mol_with_h, h_idx)
                            for radius in SUPPORTED_RADII:
                                fp = get_atom_fingerprint(mol_with_h, h_idx, radius)
                                fp_records.append((radius, fp, solvent, '1H', shift, env_str))

                                # COMPOSITE: Ajouter fingerprint composite pour 1H
                                if beyond_available:
                                    try:
                                        beyond = compute_beyond_discrete_safe(mol_with_h, h_idx, include_chiral=False)
                                        composite_fp = (fp, beyond)
                                        composite_records.append((radius, composite_fp, solvent, '1H', shift))
                                    except:
                                        pass

                            if inchikey and h_idx in canon_map:
                                exact_records.append((inchikey, solvent, '1H', canon_map[h_idx], shift))
                            count += 1
                    except (ValueError, IndexError):
                        pass

        # Parser spectre 13C
        spec_13c_match = re.search(r'<Spectrum 13C>\s*\(\d+\)\s*\n([^\n]+)', content)
        if spec_13c_match:
            spec_data_13c = spec_13c_match.group(1).strip()
            for entry in spec_data_13c.split('|'):
                if not entry.strip():
                    continue
                parts = entry.split(';')
                if len(parts) >= 3:
                    try:
                        shift = float(parts[0])
                        atom_idx = int(parts[2])
                        if atom_idx >= mol_with_h.GetNumAtoms():
                            continue
                        if shift < -5.0 or shift > 250.0:
                            continue

                        env_str = get_atom_env_string(mol_with_h, atom_idx)
                        for radius in SUPPORTED_RADII:
                            fp = get_atom_fingerprint(mol_with_h, atom_idx, radius)
                            fp_records.append((radius, fp, solvent, '13C', shift, env_str))

                            # COMPOSITE: Ajouter fingerprint composite pour 13C
                            if beyond_available:
                                try:
                                    beyond = compute_beyond_discrete_safe(mol_with_h, atom_idx, include_chiral=False)
                                    composite_fp = (fp, beyond)
                                    composite_records.append((radius, composite_fp, solvent, '13C', shift))
                                except:
                                    pass

                        if inchikey and atom_idx in canon_map:
                            exact_records.append((inchikey, solvent, '13C', canon_map[atom_idx], shift))
                        count += 1
                    except (ValueError, IndexError):
                        pass

        return (fp_records, exact_records, composite_records, count)

    except Exception:
        return None


# =============================================================================
# FINGERPRINT PREDICTOR CLASS
# =============================================================================
class FingerprintPredictor:
    """
    Prédicteur basé sur les fingerprints Morgan (radius 2 uniquement).
    Charge automatiquement la base de données au premier appel.
    """

    _instance = None
    _db = None
    _db_composite = None  # Index composite: {radius: {composite_fp: {solvent: {nucleus: [shifts]}}}}
    _db_mtime = None  # Timestamp du fichier DB pour détecter les mises à jour
    _exact_db = None  # {inchikey: {solvent: {nucleus: {canon_atom_idx: shift}}}}

    def __new__(cls):
        """Singleton pattern pour éviter de recharger la DB."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        # Vérifier si la DB doit être (re)chargée
        needs_reload = False

        if FingerprintPredictor._db is None:
            needs_reload = True
        elif os.path.exists(DB_CACHE_FILE):
            # Vérifier si le fichier a été modifié depuis le dernier chargement
            current_mtime = os.path.getmtime(DB_CACHE_FILE)
            if FingerprintPredictor._db_mtime is None or current_mtime > FingerprintPredictor._db_mtime:
                print(f"[FingerprintPredictor] Database file updated, reloading...", flush=True)
                needs_reload = True

        if needs_reload:
            self._load_database()
    
    def _load_database(self):
        """Charge la base de données depuis le cache."""
        if not os.path.exists(DB_CACHE_FILE):
            print(f"[FingerprintPredictor] Cache not found: {DB_CACHE_FILE}")
            print("[FingerprintPredictor] Building database...")
            self._build_database()
            return

        try:
            with open(DB_CACHE_FILE, 'rb') as f:
                data = pickle.load(f)

            if data.get('version') not in ('13.0',):
                print(f"[FingerprintPredictor] Cache version mismatch ({data.get('version')} != 13.0), rebuilding...")
                self._build_database()
                return

            # Charger les index pour chaque radius supporté
            FingerprintPredictor._db = {}
            for radius in SUPPORTED_RADII:
                FingerprintPredictor._db[radius] = data['indexes'].get(radius, {})
            FingerprintPredictor._fp_to_env = data.get('fp_to_env', {})
            FingerprintPredictor._exact_db = data.get('exact_db', {})

            # Charger l'index composite
            FingerprintPredictor._db_composite = {}
            composite_indexes = data.get('indexes_composite', {})
            for radius in SUPPORTED_RADII:
                FingerprintPredictor._db_composite[radius] = composite_indexes.get(radius, {})

            # Vérifier si radius 6 existe, sinon reconstruire
            if RADIUS_ULTRA_EXACT not in data['indexes'] or len(data['indexes'].get(RADIUS_ULTRA_EXACT, {})) == 0:
                print(f"[FingerprintPredictor] Radius {RADIUS_ULTRA_EXACT} not in cache, rebuilding...")
                self._build_database()
                return

            for r in SUPPORTED_RADII:
                count = len(FingerprintPredictor._db.get(r, {}))
                print(f"[FingerprintPredictor] Loaded {count} fingerprints (radius {r})")

            # Stats composite
            comp_total = sum(len(FingerprintPredictor._db_composite.get(r, {})) for r in SUPPORTED_RADII)
            if comp_total > 0:
                print(f"[FingerprintPredictor] Loaded {comp_total} composite fingerprints total")

            # Enregistrer le timestamp du fichier chargé
            FingerprintPredictor._db_mtime = os.path.getmtime(DB_CACHE_FILE)

        except Exception as e:
            print(f"[FingerprintPredictor] Error loading cache: {e}")
            self._build_database()
    
    def _build_database(self):
        """Construit la base de données depuis les sources."""
        import time
        start = time.time()

        # Index par radius: {radius: {fingerprint -> solvent -> nucleus -> [shifts]}}
        db_by_radius = {r: defaultdict(lambda: defaultdict(lambda: defaultdict(list))) for r in SUPPORTED_RADII}
        fp_to_env = {}

        # NOUVEAU: Index composite (Morgan + beyond features)
        db_composite_by_radius = {r: defaultdict(lambda: defaultdict(lambda: defaultdict(list))) for r in SUPPORTED_RADII}

        # Exact-match DB: {inchikey: {solvent: {nucleus: {canon_atom_idx: shift}}}}
        exact_db = {}
        FingerprintPredictor._exact_db = exact_db

        # Charger BMRB (avec composite)
        count_bmrb = self._load_bmrb(db_by_radius, fp_to_env, exact_db, db_composite=db_composite_by_radius)

        # Charger NMRShiftDB (parallèle - MAINTENANT avec composite)
        count_nmrdb = self._load_nmrshiftdb(db_by_radius, fp_to_env, exact_db=exact_db, db_composite=db_composite_by_radius)

        # Charger NMR Spectra (avec composite)
        count_spectra = self._load_nmr_spectra(db_by_radius, fp_to_env, exact_db=exact_db, db_composite=db_composite_by_radius)

        # Charger NPMRD (avec composite)
        count_npmrd = self._load_npmrd(db_by_radius, fp_to_env, exact_db=exact_db, db_composite=db_composite_by_radius)

        # Charger données expérimentales (avec composite - non parallèle)
        count_experimental = self._load_experimental_data(db_by_radius, fp_to_env, exact_db, db_composite=db_composite_by_radius)

        elapsed = time.time() - start
        total = count_bmrb + count_nmrdb + count_spectra + count_npmrd + count_experimental
        print(f"[FingerprintPredictor] Built database: {total} shifts in {elapsed:.1f}s")
        for r in SUPPORTED_RADII:
            print(f"[FingerprintPredictor] Radius {r}: {len(db_by_radius[r])} unique fingerprints")

        # Stats composite
        if BEYOND_AVAILABLE:
            for r in SUPPORTED_RADII:
                comp_count = len(db_composite_by_radius[r])
                if comp_count > 0:
                    print(f"[FingerprintPredictor] Composite R{r}: {comp_count} unique fingerprints")

        FingerprintPredictor._db = {r: dict(db_by_radius[r]) for r in SUPPORTED_RADII}
        FingerprintPredictor._db_composite = {r: dict(db_composite_by_radius[r]) for r in SUPPORTED_RADII}
        FingerprintPredictor._fp_to_env = fp_to_env

        # Sauvegarder le cache (si possible)
        self._save_cache()
    
    def _save_cache(self):
        """Sauvegarde le cache."""
        try:
            # Convertir en dict sérialisable
            indexes_serializable = {}
            stats = {}
            for radius in SUPPORTED_RADII:
                indexes_serializable[radius] = {}
                db = FingerprintPredictor._db.get(radius, {})
                for fp, solvent_dict in db.items():
                    indexes_serializable[radius][fp] = {}
                    for solvent, nuc_dict in solvent_dict.items():
                        indexes_serializable[radius][fp][solvent] = dict(nuc_dict)
                stats[radius] = len(db)

            # Sérialiser l'index composite
            composite_serializable = {}
            composite_stats = {}
            for radius in SUPPORTED_RADII:
                composite_serializable[radius] = {}
                db_comp = (FingerprintPredictor._db_composite or {}).get(radius, {})
                for fp, solvent_dict in db_comp.items():
                    composite_serializable[radius][fp] = {}
                    for solvent, nuc_dict in solvent_dict.items():
                        composite_serializable[radius][fp][solvent] = dict(nuc_dict)
                composite_stats[radius] = len(db_comp)

            data = {
                'version': '13.0',  # Version 13.0: beyond features + étapes 1+3+4 (14 features: +nH +hyb, config relative, +dist_hetero +dist_arom)
                'indexes': indexes_serializable,
                'indexes_composite': composite_serializable,
                'fp_to_env': FingerprintPredictor._fp_to_env,
                'exact_db': FingerprintPredictor._exact_db or {},
                'total_records': sum(stats.values()),
                'stats': {'fps_by_radius': stats, 'composite_by_radius': composite_stats}
            }

            os.makedirs(os.path.dirname(DB_CACHE_FILE), exist_ok=True)
            with open(DB_CACHE_FILE, 'wb') as f:
                pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

            print(f"[FingerprintPredictor] Cache saved to {DB_CACHE_FILE}")
        except Exception as e:
            print(f"[FingerprintPredictor] Could not save cache: {e}")
    
    def _add_record(self, db_by_radius, fp_to_env, mol, atom_idx, shift, solvent, nucleus, env_str=None, db_composite_by_radius=None):
        """
        Ajoute un record indexé par FINGERPRINT, pas par SMILES+index.

        SOLUTION AU PROBLÈME D'ORDRE DES ATOMES:
        ========================================
        Au lieu de stocker: database[smiles][atom_index] = shift
        On stocke: database[fingerprint][solvent][nucleus] = [shifts...]

        Le fingerprint est calculé sur l'ENVIRONNEMENT CHIMIQUE de l'atome,
        donc il sera IDENTIQUE peu importe que l'atome soit à l'index 0 ou 9
        dans le fichier SDF ou SMILES.

        NOUVEAU: Si db_composite_by_radius est fourni, stocke aussi le composite
        fingerprint (Morgan + beyond_features) pour une résolution plus fine.
        """
        if env_str is None:
            env_str = get_atom_env_string(mol, atom_idx)

        # NORMALISATION: Normaliser le nom du solvant pour regrouper les variantes
        solvent_normalized = normalize_solvent(solvent)

        # FILTRAGE: Ignorer les solvants invalides
        if not is_valid_solvent(solvent_normalized):
            return

        # FILTRAGE DES OUTLIERS: Rejeter les shifts aberrants
        if nucleus == '1H':
            if shift < -5.0 or shift > 14.0:
                return  # Rejeter ce shift
        elif nucleus == '13C':
            if shift < -5.0 or shift > 250.0:
                return  # Rejeter ce shift

        # Calculer beyond features une fois (si disponible)
        beyond = None
        if db_composite_by_radius is not None and BEYOND_AVAILABLE:
            try:
                # IMPORTANT: include_chiral=False pour fusionner énantiomères (mêmes shifts)
                beyond = compute_beyond_discrete_safe(mol, atom_idx, include_chiral=False)
            except:
                pass

        # Ajouter pour chaque radius
        for radius in SUPPORTED_RADII:
            fp = get_atom_fingerprint(mol, atom_idx, radius)
            db_by_radius[radius][fp][solvent_normalized][nucleus].append(shift)
            if fp not in fp_to_env:
                fp_to_env[fp] = f"r{radius}|{env_str}"

            # COMPOSITE: Stocker aussi avec beyond features
            # IMPORTANT: Utiliser le tuple directement, pas hash() qui n'est pas stable entre processus!
            if beyond is not None and db_composite_by_radius is not None:
                composite_fp = (fp, beyond)
                db_composite_by_radius[radius][composite_fp][solvent_normalized][nucleus].append(shift)

    def _store_1h_shift(self, mol_with_h, atom_idx, shift, solvent, inchikey,
                        canon_map, exact_db, db, fp_to_env, db_composite_by_radius=None):
        """
        Stocke un shift 1H de manière robuste, peu importe si atom_idx
        pointe vers le H ou vers son atome lourd porteur.

        Cette fonction gère automatiquement deux cas:
        - CAS 1: atom_idx pointe directement vers un H (NPMRD, certains NMRShiftDB)
        - CAS 2: atom_idx pointe vers l'atome lourd porteur (Experimental, BMRB)

        Args:
            mol_with_h: Molécule avec hydrogènes explicites
            atom_idx: Index de l'atome (peut être H ou atome lourd)
            shift: Chemical shift en ppm
            solvent: Solvant normalisé
            inchikey: InChIKey de la molécule (ou None)
            canon_map: Dictionnaire de mapping atom_idx -> canon_idx
            exact_db: Base de données EXACT-MATCH (ou None)
            db: Base de données fingerprint (db_by_radius)
            fp_to_env: Dictionnaire fingerprint -> environnement
            db_composite_by_radius: Base de données composite (Morgan + beyond)

        Returns:
            Nombre de shifts stockés (pour comptage)
        """
        atom = mol_with_h.GetAtomWithIdx(atom_idx)
        count = 0

        if atom.GetSymbol() == 'H':
            # CAS 1: L'index pointe DIRECTEMENT vers le H
            # Utilisé par: NPMRD, certains fichiers NMRShiftDB

            # Stocker pour fingerprint prediction (+ composite si disponible)
            self._add_record(db, fp_to_env, mol_with_h, atom_idx, shift, solvent, '1H',
                           db_composite_by_radius=db_composite_by_radius)
            count += 1

            # Stocker pour EXACT-MATCH
            if inchikey and exact_db is not None and atom_idx in canon_map:
                h_canon_idx = canon_map[atom_idx]
                exact_db.setdefault(inchikey, {}).setdefault(solvent, {}).setdefault('1H', {})[h_canon_idx] = shift

        else:
            # CAS 2: L'index pointe vers l'ATOME LOURD porteur
            # Utilisé par: Experimental Data, certains BMRB, certains NMRShiftDB

            # Trouver TOUS les H attachés à cet atome lourd
            h_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'H']

            for h_atom in h_neighbors:
                h_idx = h_atom.GetIdx()

                # Stocker pour fingerprint prediction (+ composite si disponible)
                self._add_record(db, fp_to_env, mol_with_h, h_idx, shift, solvent, '1H',
                               db_composite_by_radius=db_composite_by_radius)
                count += 1

                # Stocker pour EXACT-MATCH
                if inchikey and exact_db is not None and h_idx in canon_map:
                    h_canon_idx = canon_map[h_idx]
                    exact_db.setdefault(inchikey, {}).setdefault(solvent, {}).setdefault('1H', {})[h_canon_idx] = shift

        return count

    def _load_bmrb(self, db, fp_to_env, exact_db=None, db_composite=None) -> int:
        """
        Charge les données BMRB avec le nouveau format <NMR_Shifts>.

        BMRB est un fichier SDF unique contenant plusieurs molécules.
        Si db_composite est fourni, construit aussi l'index composite.

        Format attendu pour chaque molécule:
        >  <NMR_Shifts>  (1)
        INDEX   ATOM   SHIFT_1H   SHIFT_13C
        0       C      0.000      39.309
        5       H      3.105      0.000
        ...
        """
        from rdkit.Chem import inchi as rdkit_inchi

        if not os.path.exists(BMRB_DIR):
            return 0

        # Chercher le fichier SDF BMRB (un seul fichier multi-molécule)
        sdf_files = [f for f in os.listdir(BMRB_DIR) if f.endswith('.sdf')]
        if not sdf_files:
            print("[FingerprintPredictor] BMRB: No SDF file found")
            return 0

        # Prendre le premier fichier SDF (devrait être unique)
        bmrb_file = os.path.join(BMRB_DIR, sdf_files[0])
        count = 0
        exact_count = 0

        try:
            # Utiliser SDMolSupplier pour lire toutes les molécules du fichier
            suppl = Chem.SDMolSupplier(bmrb_file, removeHs=False, sanitize=True)

            for mol_idx, mol in enumerate(suppl):
                if mol is None:
                    continue

                if mol_idx % 100 == 0:
                    mol_name = mol.GetProp('_Name') if mol.HasProp('_Name') else f'mol#{mol_idx}'
                    print(f"  [BMRB] mol {mol_idx}: {mol_name}", flush=True)

                try:
                    # Extraire le solvant
                    raw_solvent = mol.GetProp('Solvent') if mol.HasProp('Solvent') else 'Unknown'
                    solvent = normalize_solvent(raw_solvent)

                    # Extraire les shifts NMR avec le nouveau format
                    if not mol.HasProp('NMR_Shifts'):
                        continue

                    nmr_shifts_text = mol.GetProp('NMR_Shifts').strip()
                    shift_lines = [line.strip() for line in nmr_shifts_text.split('\n') if line.strip()]

                    # Ajouter les H si nécessaire
                    if mol.GetNumAtoms() > 0 and all(atom.GetSymbol() != 'H' for atom in mol.GetAtoms()):
                        mol = Chem.AddHs(mol)

                    # Précalculer charges Gasteiger pour beyond features (composite)
                    if db_composite is not None and BEYOND_AVAILABLE:
                        ensure_gasteiger_charges(mol)

                    # Calculer InChIKey et canon_map pour exact_db
                    # canon_map: atom_idx_in_mol_with_h → canon_idx
                    inchikey = None
                    canon_map = {}
                    if exact_db is not None:
                        try:
                            mol_no_h = Chem.RemoveHs(mol)
                            inchikey = rdkit_inchi.MolToInchiKey(mol_no_h)
                            if inchikey:
                                # Mapper les atomes lourds
                                canon_smiles = Chem.MolToSmiles(mol_no_h)
                                mol_canon_no_h = Chem.MolFromSmiles(canon_smiles)
                                if mol_canon_no_h is not None:
                                    match_no_h = _substruct_match_with_timeout(mol_no_h, mol_canon_no_h, timeout=1.0)
                                    # match[canon_idx] = sdf_idx dans mol_no_h
                                    noh_to_canon = {sdf_i: canon_i for canon_i, sdf_i in enumerate(match_no_h)}
                                    # Mapper les atomes lourds (indices 0 à num_heavy-1 dans mol et mol_with_h)
                                    for noh_idx, canon_idx in noh_to_canon.items():
                                        canon_map[noh_idx] = canon_idx

                                    # Mapper les hydrogènes en utilisant mol_canon avec H
                                    # ✅ FIX: Utiliser mol_with_h (avec H explicites) pour le matching
                                    # car les indices dans NMR_Shifts référencent mol après AddHs
                                    mol_canon_with_h = Chem.AddHs(mol_canon_no_h)
                                    mol_with_h_for_match = Chem.AddHs(mol)
                                    match_with_h = _substruct_match_with_timeout(mol_with_h_for_match, mol_canon_with_h, timeout=1.0)
                                    if match_with_h:
                                        # match_with_h[canon_idx] = sdf_idx dans mol_with_h
                                        for canon_idx, sdf_idx in enumerate(match_with_h):
                                            canon_map[sdf_idx] = canon_idx
                        except Exception:
                            inchikey = None

                    # Parser chaque ligne: INDEX ATOM SHIFT_1H SHIFT_13C
                    for line in shift_lines:
                        parts = line.split()
                        if len(parts) != 4:
                            continue

                        try:
                            atom_idx = int(parts[0])
                            atom_symbol = parts[1]
                            shift_1h = float(parts[2])
                            shift_13c = float(parts[3])
                        except (ValueError, IndexError):
                            continue

                        # Vérifier que l'index est valide
                        if atom_idx >= mol.GetNumAtoms():
                            continue

                        # Ajouter le shift 1H si non nul
                        if shift_1h != 0.0 and atom_symbol == 'H':
                            # Utiliser la fonction générique qui gère H ou atome lourd
                            added = self._store_1h_shift(mol, atom_idx, shift_1h, solvent,
                                                         inchikey, canon_map, exact_db, db, fp_to_env,
                                                         db_composite_by_radius=db_composite)
                            count += added
                            if added > 0:
                                exact_count += added

                        # Ajouter le shift 13C si non nul
                        if shift_13c != 0.0 and atom_symbol == 'C':
                            self._add_record(db, fp_to_env, mol, atom_idx, shift_13c, solvent, '13C',
                                           db_composite_by_radius=db_composite)
                            count += 1
                            if inchikey and exact_db is not None and atom_idx in canon_map:
                                canon_idx = canon_map[atom_idx]
                                exact_db.setdefault(inchikey, {}).setdefault(solvent, {}).setdefault('13C', {})[canon_idx] = shift_13c
                                exact_count += 1

                except Exception:
                    # Ignorer les erreurs sur une molécule spécifique
                    continue

        except Exception:
            # Erreur fatale sur le fichier
            pass

        print(f"[FingerprintPredictor] BMRB: {count} shifts, {exact_count} exact-match entries ({len(exact_db) if exact_db else 0} molecules total in exact_db)")
        return count
    
    def _load_nmrshiftdb(self, db, fp_to_env, exact_db=None, db_composite=None) -> int:
        """Charge les données NMRShiftDB2 avec traitement PARALLÈLE.
        Supporte maintenant les fingerprints composites (Morgan + beyond)."""
        if not os.path.exists(NMRSHIFTDB_FILE):
            return 0

        def _iter_mol_blocks(filepath):
            current_mol = []
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if line.strip() == '$$$$':
                        if current_mol:
                            yield '\n'.join(current_mol)
                            current_mol = []
                    else:
                        current_mol.append(line.rstrip())

        import time as _time
        _t0 = _time.time()

        # Collecter tous les mol_blocks avec leur index
        print("  [NMRShiftDB] Loading mol blocks...", flush=True)
        mol_blocks = [(block, i+1) for i, block in enumerate(_iter_mol_blocks(NMRSHIFTDB_FILE))]
        total_mols = len(mol_blocks)
        print(f"  [NMRShiftDB] {total_mols} molecules loaded in {_time.time() - _t0:.1f}s", flush=True)

        # Déterminer le nombre de workers (utiliser tous les CPUs disponibles)
        n_workers = max(1, multiprocessing.cpu_count() - 1)
        print(f"  [NMRShiftDB] Starting parallel processing with {n_workers} workers...", flush=True)

        count = 0
        exact_count = 0
        processed = 0
        batch_size = 1000  # Taille de batch pour affichage progression

        # Traitement parallèle
        _t_parallel = _time.time()
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            # Soumettre tous les jobs
            futures = {executor.submit(_process_mol_block_worker, args): args[1] for args in mol_blocks}

            for future in as_completed(futures):
                processed += 1

                if processed % batch_size == 0:
                    elapsed = _time.time() - _t_parallel
                    rate = processed / elapsed if elapsed > 0 else 0
                    eta = (total_mols - processed) / rate if rate > 0 else 0
                    print(f"  [NMRShiftDB] {processed}/{total_mols} | {count} shifts | {elapsed:.0f}s | {rate:.0f} mol/s | ETA {eta/60:.1f}min", flush=True)

                result = future.result()
                if result is None:
                    continue

                # MODIFIÉ: Le worker retourne maintenant 4 éléments (avec composite_records)
                if len(result) == 4:
                    fp_records, exact_records, composite_records, rec_count = result
                else:
                    # Rétrocompatibilité
                    fp_records, exact_records, rec_count = result
                    composite_records = []

                count += rec_count

                # Intégrer les résultats dans la base
                for radius, fp, solvent, nucleus, shift, env_str in fp_records:
                    db[radius][fp][solvent][nucleus].append(shift)
                    if fp not in fp_to_env:
                        fp_to_env[fp] = f"r{radius}|{env_str}"

                # NOUVEAU: Intégrer les fingerprints composites
                if db_composite is not None and composite_records:
                    for radius, composite_fp, solvent, nucleus, shift in composite_records:
                        db_composite[radius][composite_fp][solvent][nucleus].append(shift)

                # Intégrer les exact records
                if exact_db is not None:
                    for inchikey, solvent, nucleus, canon_idx, shift in exact_records:
                        exact_db.setdefault(inchikey, {}).setdefault(solvent, {}).setdefault(nucleus, {})[canon_idx] = shift
                        exact_count += 1

        elapsed_total = _time.time() - _t0
        print(f"[FingerprintPredictor] NMRShiftDB: {count} shifts, {exact_count} exact-match entries ({len(exact_db) if exact_db else 0} molecules) in {elapsed_total:.1f}s ({n_workers} workers)")
        return count
    
    def _load_nmr_spectra(self, db, fp_to_env, exact_db=None, db_composite=None) -> int:
        """Charge les données NMR Spectra format avec traitement PARALLÈLE.
        Supporte maintenant les fingerprints composites (Morgan + beyond)."""
        if not os.path.exists(NMR_SPECTRA_FILE):
            return 0

        def _iter_mol_blocks(filepath):
            current_mol = []
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if line.strip() == '$$$$':
                        if current_mol:
                            yield '\n'.join(current_mol)
                            current_mol = []
                    else:
                        current_mol.append(line.rstrip())

        import time as _time
        _t0 = _time.time()

        # Collecter tous les mol_blocks
        print("  [NMRSpectra] Loading mol blocks...", flush=True)
        mol_blocks = [(block, i+1) for i, block in enumerate(_iter_mol_blocks(NMR_SPECTRA_FILE))]
        total_mols = len(mol_blocks)
        print(f"  [NMRSpectra] {total_mols} molecules loaded in {_time.time() - _t0:.1f}s", flush=True)

        n_workers = max(1, multiprocessing.cpu_count() - 1)
        print(f"  [NMRSpectra] Starting parallel processing with {n_workers} workers...", flush=True)

        count = 0
        exact_count = 0
        processed = 0
        batch_size = 1000

        _t_parallel = _time.time()
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(_process_nmr_spectra_worker, args): args[1] for args in mol_blocks}

            for future in as_completed(futures):
                processed += 1
                if processed % batch_size == 0:
                    elapsed = _time.time() - _t_parallel
                    rate = processed / elapsed if elapsed > 0 else 0
                    eta = (total_mols - processed) / rate if rate > 0 else 0
                    print(f"  [NMRSpectra] {processed}/{total_mols} | {count} shifts | {elapsed:.0f}s | {rate:.0f} mol/s | ETA {eta/60:.1f}min", flush=True)

                result = future.result()
                if result is None:
                    continue

                # MODIFIÉ: Le worker retourne maintenant 4 éléments (avec composite_records)
                if len(result) == 4:
                    fp_records, exact_records, composite_records, rec_count = result
                else:
                    fp_records, exact_records, rec_count = result
                    composite_records = []

                count += rec_count

                for radius, fp, solvent, nucleus, shift, env_str in fp_records:
                    db[radius][fp][solvent][nucleus].append(shift)
                    if fp not in fp_to_env:
                        fp_to_env[fp] = f"r{radius}|{env_str}"

                # NOUVEAU: Intégrer les fingerprints composites
                if db_composite is not None and composite_records:
                    for radius, composite_fp, solvent, nucleus, shift in composite_records:
                        db_composite[radius][composite_fp][solvent][nucleus].append(shift)

                if exact_db is not None:
                    for inchikey, solvent, nucleus, canon_idx, shift in exact_records:
                        exact_db.setdefault(inchikey, {}).setdefault(solvent, {}).setdefault(nucleus, {})[canon_idx] = shift
                        exact_count += 1

        elapsed_total = _time.time() - _t0
        print(f"[FingerprintPredictor] NMRSpectra: {count} shifts, {exact_count} exact-match entries ({len(exact_db) if exact_db else 0} molecules) in {elapsed_total:.1f}s ({n_workers} workers)")
        return count
    
    def _load_npmrd(self, db, fp_to_env, exact_db=None, db_composite=None) -> int:
        """Charge les fichiers SDF NPMRD avec traitement PARALLÈLE."""
        if not os.path.exists(NPMRD_SDF_DIR):
            return 0

        import glob
        import time as _time
        _t0 = _time.time()

        sdf_files = glob.glob(os.path.join(NPMRD_SDF_DIR, '*.sdf'))
        total_files = len(sdf_files)
        print(f"  [NPMRD] {total_files} files found", flush=True)

        n_workers = max(1, multiprocessing.cpu_count() - 1)
        print(f"  [NPMRD] Starting parallel processing with {n_workers} workers...", flush=True)

        count = 0
        exact_count = 0
        composite_count = 0
        processed = 0
        batch_size = 100

        # Préparer les arguments pour le worker
        file_args = [(filepath, i) for i, filepath in enumerate(sdf_files)]

        _t_parallel = _time.time()
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(_process_npmrd_file_worker, args): args[1] for args in file_args}

            for future in as_completed(futures):
                processed += 1
                if processed % batch_size == 0:
                    elapsed = _time.time() - _t_parallel
                    rate = processed / elapsed if elapsed > 0 else 0
                    eta = (total_files - processed) / rate if rate > 0 else 0
                    print(f"  [NPMRD] {processed}/{total_files} | {count} shifts | {elapsed:.0f}s | {rate:.0f} files/s | ETA {eta:.1f}s", flush=True)

                result = future.result()
                if result is None:
                    continue

                # Handle 4-element tuple (with composite) or 3-element tuple (legacy)
                if len(result) == 4:
                    fp_records, exact_records, composite_records, rec_count = result
                else:
                    fp_records, exact_records, rec_count = result
                    composite_records = []
                count += rec_count

                for radius, fp, solvent, nucleus, shift, env_str in fp_records:
                    db[radius][fp][solvent][nucleus].append(shift)
                    if fp not in fp_to_env:
                        fp_to_env[fp] = f"r{radius}|{env_str}"

                if exact_db is not None:
                    for inchikey, solvent, nucleus, canon_idx, shift in exact_records:
                        exact_db.setdefault(inchikey, {}).setdefault(solvent, {}).setdefault(nucleus, {})[canon_idx] = shift
                        exact_count += 1

                # COMPOSITE: Stocker les fingerprints composites
                if db_composite is not None and composite_records:
                    for radius, composite_fp, solvent, nucleus, shift in composite_records:
                        db_composite[radius][composite_fp][solvent][nucleus].append(shift)
                        composite_count += 1

        elapsed_total = _time.time() - _t0
        print(f"[FingerprintPredictor] NPMRD: {count} shifts from {total_files} files, {exact_count} exact-match entries, {composite_count} composite ({len(exact_db) if exact_db else 0} molecules) in {elapsed_total:.1f}s ({n_workers} workers)")
        return count

    def _parse_npmrd_spectrum(self, mol_with_h, spec_data, solvent, nucleus, db, fp_to_env, exact_db=None, inchikey=None, canon_map=None) -> tuple:
        """Parse les données de spectre NPMRD.

        Format:
        - 13C: shift;mult;atom_idx|shift;mult;atom_idx|...
        - 1H: shift;mult;atom_idx;atom_idx;...|shift;mult;atom_idx|...

        Note: Les indices dans NPMRD sont 1-based, il faut soustraire 1 pour RDKit.

        Returns: (count, exact_count)
        """
        if canon_map is None:
            canon_map = {}
        count = 0
        exact_count = 0
        entries = spec_data.split('|')

        for entry in entries:
            if not entry.strip():
                continue

            parts = entry.split(';')
            if len(parts) < 3:
                continue

            try:
                shift = float(parts[0])

                # Les indices d'atomes sont dans parts[2:]
                # IMPORTANT: Les indices dans les SDF NPMRD commencent à 1, mais RDKit utilise 0-based
                # Il faut donc soustraire 1 de chaque indice
                atom_indices = []
                for idx_str in parts[2:]:
                    try:
                        sdf_idx = int(idx_str)
                        # Convertir de 1-based (SDF) à 0-based (RDKit)
                        atom_idx = sdf_idx - 1
                        atom_indices.append(atom_idx)
                    except ValueError:
                        pass

                if not atom_indices:
                    continue

                # Pour 13C: un seul index par shift
                # Pour 1H: potentiellement plusieurs indices (protons équivalents)
                # On utilise seulement le premier pour le fingerprint
                atom_idx = atom_indices[0]

                # Vérifier que l'index est valide
                if atom_idx >= mol_with_h.GetNumAtoms():
                    continue

                # Vérifier que c'est le bon type d'atome
                atom = mol_with_h.GetAtomWithIdx(atom_idx)
                atom_symbol = atom.GetSymbol()

                if nucleus == '13C' and atom_symbol != 'C':
                    # Pour 13C, l'index doit pointer vers un C
                    # Note: parfois l'index pointe vers un H attaché au C
                    # Dans ce cas, utiliser le parent
                    if atom_symbol == 'H':
                        neighbors = list(atom.GetNeighbors())
                        if neighbors and neighbors[0].GetSymbol() == 'C':
                            atom_idx = neighbors[0].GetIdx()
                        else:
                            continue
                    else:
                        continue

                elif nucleus == '1H':
                    # Pour 1H: utiliser la fonction générique qui gère automatiquement
                    # - Cas où atom_idx pointe vers un H directement
                    # - Cas où atom_idx pointe vers le C porteur
                    added = self._store_1h_shift(mol_with_h, atom_idx, shift, solvent,
                                                 inchikey, canon_map, exact_db, db, fp_to_env)
                    count += added
                    if inchikey and exact_db is not None:
                        exact_count += added
                    continue

                # Pour 13C: ajouter le record normalement
                self._add_record(db, fp_to_env, mol_with_h, atom_idx, shift, solvent, nucleus)
                count += 1
                if inchikey and exact_db is not None and atom_idx in canon_map:
                    canon_idx = canon_map[atom_idx]
                    exact_db.setdefault(inchikey, {}).setdefault(solvent, {}).setdefault(nucleus, {})[canon_idx] = shift
                    exact_count += 1

            except (ValueError, IndexError) as e:
                # Ignorer les entrées mal formatées
                pass

        return count, exact_count

    def _load_experimental_data(self, db, fp_to_env, exact_db=None, db_composite=None) -> int:
        """Charge les données expérimentales avec traitement PARALLÈLE."""
        if not os.path.exists(EXPERIMENTAL_DATA_DIR):
            return 0

        import glob
        import time as _time
        _t0 = _time.time()

        sdf_files = glob.glob(os.path.join(EXPERIMENTAL_DATA_DIR, '*.sdf'))
        total_files = len(sdf_files)
        print(f"  [Experimental] {total_files} files found", flush=True)

        if total_files == 0:
            return 0

        n_workers = max(1, multiprocessing.cpu_count() - 1)
        print(f"  [Experimental] Starting parallel processing with {n_workers} workers...", flush=True)

        count = 0
        exact_count = 0
        composite_count = 0
        processed = 0
        batch_size = max(10, total_files // 10)

        file_args = [(filepath, i) for i, filepath in enumerate(sdf_files)]

        _t_parallel = _time.time()
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(_process_experimental_file_worker, args): args[1] for args in file_args}

            for future in as_completed(futures):
                processed += 1
                if processed % batch_size == 0:
                    elapsed = _time.time() - _t_parallel
                    rate = processed / elapsed if elapsed > 0 else 0
                    print(f"  [Experimental] {processed}/{total_files} | {count} shifts | {elapsed:.0f}s | {rate:.0f} files/s", flush=True)

                result = future.result()
                if result is None:
                    continue

                # Handle 4-element tuple (with composite) or 3-element tuple (legacy)
                if len(result) == 4:
                    fp_records, exact_records, composite_records, rec_count = result
                else:
                    fp_records, exact_records, rec_count = result
                    composite_records = []
                count += rec_count

                for radius, fp, solvent, nucleus, shift, env_str in fp_records:
                    db[radius][fp][solvent][nucleus].append(shift)
                    if fp not in fp_to_env:
                        fp_to_env[fp] = f"r{radius}|{env_str}"

                if exact_db is not None:
                    for inchikey, solvent, nucleus, canon_idx, shift in exact_records:
                        exact_db.setdefault(inchikey, {}).setdefault(solvent, {}).setdefault(nucleus, {})[canon_idx] = shift
                        exact_count += 1

                # COMPOSITE: Stocker les fingerprints composites
                if db_composite is not None and composite_records:
                    for radius, composite_fp, solvent, nucleus, shift in composite_records:
                        db_composite[radius][composite_fp][solvent][nucleus].append(shift)
                        composite_count += 1

        elapsed_total = _time.time() - _t0
        print(f"[FingerprintPredictor] Experimental data: {count} shifts, {exact_count} exact-match entries, {composite_count} composite in {elapsed_total:.1f}s ({n_workers} workers)")
        return count

    def _select_best_n(self, shifts: list, n: int = 10) -> list:
        """
        Si len(shifts) > n, sélectionne les n valeurs les plus cohérentes.
        
        Algorithme robuste:
        1. IQR filtering pour éliminer outliers extrêmes
        2. Clustering par densité pour trouver le groupe principal
        3. Prendre les n valeurs les plus proches du centre de densité
        
        Gère correctement les distributions bimodales et les outliers.
        """
        import numpy as np
        
        if len(shifts) <= n:
            return shifts
        
        shifts_array = np.array(shifts)
        
        # Étape 1: IQR filtering pour éliminer outliers extrêmes
        q1, q3 = np.percentile(shifts_array, [25, 75])
        iqr = q3 - q1
        
        if iqr > 0:  # Protection division par zéro
            # Filtrer outliers à ±1.5×IQR
            lower_bound = q1 - 1.5 * iqr
            upper_bound = q3 + 1.5 * iqr
            filtered = shifts_array[(shifts_array >= lower_bound) & (shifts_array <= upper_bound)]
        else:
            filtered = shifts_array
        
        # Si après filtrage on a ≤ n valeurs, retourner tout
        if len(filtered) <= n:
            return filtered.tolist()
        
        # Étape 2: Trouver le centre de densité maximale
        sorted_filtered = np.sort(filtered)
        densities = np.array([np.sum(np.abs(sorted_filtered - x) < np.std(sorted_filtered) * 0.5) 
                             for x in sorted_filtered])
        center_idx = np.argmax(densities)
        center_value = sorted_filtered[center_idx]
        
        # Étape 3: Prendre les n valeurs les plus proches du centre
        distances = np.abs(sorted_filtered - center_value)
        nearest_indices = np.argsort(distances)[:n]
        
        return sorted_filtered[nearest_indices].tolist()

    def _is_bimodal(self, shifts_list: list, nucleus: str) -> bool:
        """
        Détecte si une distribution de shifts est bimodale (deux pics distincts).

        USAGE:
        ======
        Une distribution bimodale indique que deux environnements chimiques
        DIFFÉRENTS ont été fusionnés par erreur dans le même fingerprint.

        Exemples de causes:
        - CH₃ ester et CH₃ alkyle fusionnés (avant étape 4)
        - CH₂ benzylique et CH₂ alkyle fusionnés (avant étape 4)
        - Stéréoisomères E/Z non distingués (avant étape 3)
        - Énantiomères dans solvant chiral (rare)

        Si bimodal détecté → retourner None (fallback ML/SGNN)

        Algorithme:
        -----------
        1. Coefficient de bimodalité de Sarle: b = (skew² + 1) / kurtosis
           - b > 0.555 → bimodal (distribution uniforme)
           - b > 0.59 → fortement bimodal
        2. Validation par gap analysis: chercher un gap significatif
           - Gap > 2×std entre deux groupes

        Args:
            shifts_list: Liste des shifts chimiques (après filtrage IQR)
            nucleus: '1H' ou '13C'

        Returns:
            True si distribution bimodale, False sinon
        """
        import numpy as np

        # DÉSACTIVÉ: Le test bimodal est trop agressif pour les radius spécifiques (R6/R5)
        # Même quand bimodal détecté (ex: CDCl3 24 molécules, coeff=0.840),
        # le fallback vers ML donne des résultats similaires (1.263 vs 1.287 ppm).
        # De plus, un des groupes bimodaux peut être le bon (ex: groupe 2 CDCl3 = 1.553 ppm)
        # Mieux vaut utiliser la moyenne de toutes les données disponibles.
        return False  # ← PATCH: Désactiver complètement le test bimodal

        # Ancien seuil: n>=10 était trop restrictif pour R6 ultra-spécifique
        # if len(shifts_list) < 10:
        #     return False

        shifts_array = np.array(shifts_list)

        # ================================================================
        # MÉTHODE 1: Coefficient de bimodalité de Sarle
        # ================================================================
        # Formule: b = (skew² + 1) / kurtosis
        # b > 0.555 → bimodal (uniforme)
        # b > 0.59 → fortement bimodal

        try:
            from scipy import stats as scipy_stats

            skewness = scipy_stats.skew(shifts_array)
            # kurtosis excess (Fisher=True) : kurtosis - 3
            kurtosis_excess = scipy_stats.kurtosis(shifts_array, fisher=True)
            # Kurtosis de Pearson (non-excess)
            kurtosis_pearson = kurtosis_excess + 3

            # Protection division par zéro
            if kurtosis_pearson <= 0:
                kurtosis_pearson = 0.01

            bimodality_coeff = (skewness**2 + 1) / kurtosis_pearson

            # Seuil: 0.555 (uniforme) mais on utilise 0.59 pour être plus strict
            # (éviter faux positifs sur distributions légèrement irrégulières)
            if bimodality_coeff > 0.59:
                logger.debug(f"Bimodal detected: Sarle coefficient = {bimodality_coeff:.3f} > 0.59")
                return True
        except Exception as e:
            logger.debug(f"Sarle coefficient calculation failed (scipy may not be available): {e}")

        # ================================================================
        # MÉTHODE 2: Gap analysis (validation supplémentaire)
        # ================================================================
        # Chercher un gap significatif entre deux groupes
        # Si gap > 2×std ET chaque groupe a ≥3 valeurs → bimodal

        sorted_shifts = np.sort(shifts_array)
        overall_std = np.std(shifts_array)

        # Seuil de gap basé sur le nucleus
        # 1H: gap > 0.5 ppm ET > 2×std
        # 13C: gap > 5 ppm ET > 2×std
        min_gap_absolute = 0.5 if nucleus == '1H' else 5.0
        min_gap_relative = 2.0 * overall_std
        min_gap = max(min_gap_absolute, min_gap_relative)

        # Chercher le plus grand gap
        gaps = np.diff(sorted_shifts)
        max_gap_idx = np.argmax(gaps)
        max_gap = gaps[max_gap_idx]

        if max_gap > min_gap:
            # Vérifier que chaque groupe a au moins 3 valeurs
            group1_size = max_gap_idx + 1
            group2_size = len(sorted_shifts) - group1_size

            if group1_size >= 3 and group2_size >= 3:
                logger.debug(f"Bimodal detected: Gap = {max_gap:.2f} ppm > {min_gap:.2f} ppm (groups: {group1_size}+{group2_size})")
                return True

        # ================================================================
        # MÉTHODE 3: Dip test (optionnel - si scipy disponible)
        # ================================================================
        # Le dip test de Hartigan est le test statistique formel pour bimodalité
        # Malheureusement, il n'est pas dans scipy standard (besoin de diptest package)
        # On peut l'ajouter si nécessaire

        return False

    def predict_atom(self, mol_with_h, atom_idx: int, solvent: str, nucleus: str) -> Optional[Tuple[float, float, int, int]]:
        """
        Prédit le shift par MATCHING DE FINGERPRINT CHIMIQUE.
        
        🔑 SOLUTION AU PROBLÈME D'ORDRE DES ATOMES:
        ============================================
        Cette méthode cherche par FINGERPRINT (empreinte chimique) au lieu de
        chercher par SMILES+index. Donc peu importe si l'atome est à l'index 0
        dans le SMILES et index 9 dans le SDF : s'ils ont le même environnement
        chimique, ils auront le MÊME fingerprint et trouveront le MÊME shift!
        
        RÈGLE IMPORTANTE: Les protons échangeables (OH, NH, SH) utilisent
        TOUJOURS le modèle ML car leurs shifts varient énormément selon
        les conditions expérimentales.
        
        HIÉRARCHIE avec fallback std (1H: >0.3 ppm, 13C: >3.0 ppm):
          Pour chaque radius R4→R3→R2:
            1. Calcule fingerprint de l'atome query
            2. Cherche ce fingerprint dans la base (solvant exact)
            3. Si std acceptable: retourne
            4. Sinon essaye radius suivant
            5. Si tous radius échouent → ML
        
        Args:
            mol_with_h: Molécule RDKit avec H explicites (PEU IMPORTE L'ORDRE DES ATOMES)
            atom_idx: Index de l'atome dans mol_with_h
            solvent: Solvant (CDCl3, DMSO-d6, etc.)
            nucleus: '1H' ou '13C'
        
        Returns:
            (shift_mean, shift_std, count, radius_used) ou None si pas trouvé
        """
        if FingerprintPredictor._db is None:
            return None

        # RÈGLE: Protons échangeables utilisent toujours ML
        if nucleus == '1H' and is_exchangeable_proton(mol_with_h, atom_idx):
            logger.debug(f"Atom {atom_idx}: Exchangeable proton → ML")
            return None

        solvent_norm = normalize_solvent(solvent)
        import numpy as np

        # EXACT-MATCH: vérifier si la molécule est connue exactement par InChIKey
        if FingerprintPredictor._exact_db:
            try:
                from rdkit.Chem import inchi as rdkit_inchi
                mol_no_h = Chem.RemoveHs(mol_with_h)
                inchikey = rdkit_inchi.MolToInchiKey(mol_no_h)
                if inchikey and inchikey in FingerprintPredictor._exact_db:
                    exact_entry = FingerprintPredictor._exact_db[inchikey]
                    # Trouver le canon_idx de l'atome courant
                    # On passe par SMILES canonique pour obtenir le mapping
                    canon_smiles = Chem.MolToSmiles(mol_no_h)
                    mol_canon_no_h = Chem.MolFromSmiles(canon_smiles)
                    if mol_canon_no_h is not None:
                        atom = mol_with_h.GetAtomWithIdx(atom_idx)
                        canon_idx = None

                        if atom.GetSymbol() != 'H':
                            # Pour les atomes lourds: mapper via mol_no_h
                            match_no_h = _substruct_match_with_timeout(mol_no_h, mol_canon_no_h, timeout=1.0)
                            sdf_to_canon = {sdf_i: canon_i for canon_i, sdf_i in enumerate(match_no_h)}
                            # Compter combien d'atomes lourds sont avant atom_idx dans mol_with_h
                            heavy_idx_in_noh = sum(
                                1 for i in range(atom_idx)
                                if mol_with_h.GetAtomWithIdx(i).GetAtomicNum() > 1
                            )
                            canon_idx = sdf_to_canon.get(heavy_idx_in_noh)
                        else:
                            # Pour les H: mapper via mol_canon_with_h
                            mol_canon_with_h = Chem.AddHs(mol_canon_no_h)
                            match_with_h = _substruct_match_with_timeout(mol_with_h, mol_canon_with_h, timeout=1.0)
                            if match_with_h:
                                # match_with_h[canon_idx] = sdf_idx dans mol_with_h
                                # On veut: sdf_idx → canon_idx
                                sdf_to_canon_wh = {sdf_i: canon_i for canon_i, sdf_i in enumerate(match_with_h)}
                                canon_idx = sdf_to_canon_wh.get(atom_idx)
                        if canon_idx is not None:
                            for sol_key in [solvent_norm] + list(exact_entry.keys()):
                                nuc_dict = exact_entry.get(sol_key, {}).get(nucleus, {})
                                # Keys can be int or str depending on how the data was stored
                                shift = nuc_dict.get(canon_idx) or nuc_dict.get(str(canon_idx))
                                if shift is not None:
                                    logger.debug(f"Atom {atom_idx}: Exact-match InChIKey {inchikey}, shift={shift}")
                                    return shift, 0.0, 99, 0  # count=99, radius=0 = exact match
            except Exception as e:
                logger.debug(f"Exact-match error: {e}")

        # Déterminer le rayon minimum acceptable pour cet atome
        min_radius_required = classify_carbon_radius_requirement(mol_with_h, atom_idx)

        # ML obligatoire pour certains environnements (gem-diméthyl E/Z)
        if min_radius_required == 99:
            logger.debug(f"Atom {atom_idx}: Requires 3D/ML (gem-dimethyl E/Z) → skip fingerprint")
            return None

        # ════════════════════════════════════════════════════════════════════
        # NIVEAU 2: COMPOSITE FINGERPRINT (Morgan + Beyond)
        # ════════════════════════════════════════════════════════════════════
        # Le composite distingue ester/cétone/amide → environnement exact
        # MIN_SAMPLES = 1 car l'environnement est pur (pas de mélange)
        if BEYOND_AVAILABLE and FingerprintPredictor._db_composite:
            try:
                ensure_gasteiger_charges(mol_with_h)
                # IMPORTANT: include_chiral=False pour fusionner énantiomères (mêmes shifts)
                beyond = compute_beyond_discrete_safe(mol_with_h, atom_idx, include_chiral=False)

                for radius in [RADIUS_ULTRA_EXACT, RADIUS_HYPER_EXACT, RADIUS_VERY_EXACT, RADIUS_EXACT, RADIUS_FALLBACK]:
                    if radius < min_radius_required:
                        continue

                    db_comp = FingerprintPredictor._db_composite.get(radius, {})
                    if not db_comp:
                        continue

                    fp = get_atom_fingerprint(mol_with_h, atom_idx, radius)
                    # IMPORTANT: Utiliser le tuple directement, pas hash() qui n'est pas stable entre processus!
                    composite_fp = (fp, beyond)

                    if composite_fp in db_comp:
                        # Hiérarchie de solvants pour composite (identique à Morgan)
                        COMPOSITE_SOLVENT_FALLBACK = {
                            'CDCl3': ['Acetone-d6', 'MeOD', 'DMSO-d6', 'Unreported'],
                            'Acetone-d6': ['CDCl3', 'MeOD', 'DMSO-d6', 'Unreported'],
                            'MeOD': ['Acetone-d6', 'CDCl3', 'DMSO-d6', 'Unreported'],
                            'DMSO-d6': ['MeOD', 'Acetone-d6', 'CDCl3', 'Unreported'],
                            'Unreported': ['CDCl3', 'DMSO-d6', 'MeOD', 'Acetone-d6'],
                            'D2O': []  # Pas de fallback (trop différent)
                        }
                        composite_fallbacks = COMPOSITE_SOLVENT_FALLBACK.get(solvent_norm, [])

                        # Seuils minimum par radius (identiques à Morgan)
                        MIN_SAMPLES_COMPOSITE = {
                            RADIUS_ULTRA_EXACT: MIN_SAMPLES_R6,  # R6 = 1
                            RADIUS_HYPER_EXACT: MIN_SAMPLES_R5,  # R5 = 1
                            RADIUS_VERY_EXACT: MIN_SAMPLES_R4,   # R4 = 2
                            RADIUS_EXACT: MIN_SAMPLES_R3,        # R3 = 3
                            RADIUS_FALLBACK: MIN_SAMPLES_R2      # R2 = 3
                        }
                        min_samples_comp = MIN_SAMPLES_COMPOSITE.get(radius, 3)

                        # Collecter shifts du solvant exact + variantes (fusion runtime)
                        solvent_variants = get_solvent_variants(solvent_norm)
                        exact_shifts = []
                        for sv in solvent_variants:
                            if sv in db_comp[composite_fp]:
                                exact_shifts.extend(db_comp[composite_fp][sv].get(nucleus, []))

                        if len(exact_shifts) >= min_samples_comp:
                            # ÉTAPE 5: Vérifier bimodalité avant de retourner
                            if self._is_bimodal(exact_shifts, nucleus):
                                logger.debug(f"Atom {atom_idx}: COMPOSITE R{radius} {solvent_norm} BIMODAL → fallback ML")
                                continue  # Essayer radius suivant ou fallback ML

                            import numpy as np
                            mean = float(np.mean(exact_shifts))
                            std = float(np.std(exact_shifts)) if len(exact_shifts) > 1 else 0.0
                            logger.debug(f"Atom {atom_idx}: COMPOSITE R{radius} {solvent_norm}, n={len(exact_shifts)}")
                            return mean, std, len(exact_shifts), radius

                        # Solvant exact a trop peu de données - essayer fallbacks HIÉRARCHIQUEMENT
                        # (dans l'ordre de priorité, pas agrégation globale)
                        for try_solvent in composite_fallbacks:
                            # Fusionner les variantes pour ce fallback
                            fallback_variants = get_solvent_variants(try_solvent)
                            fallback_shifts = []
                            for fv in fallback_variants:
                                if fv in db_comp[composite_fp]:
                                    fallback_shifts.extend(db_comp[composite_fp][fv].get(nucleus, []))

                            # Si ce fallback a assez de données, le tester
                            if len(fallback_shifts) >= min_samples_comp:
                                # ÉTAPE 5: Vérifier bimodalité
                                if self._is_bimodal(fallback_shifts, nucleus):
                                    logger.debug(f"Atom {atom_idx}: COMPOSITE R{radius} {try_solvent} BIMODAL → essayer fallback suivant")
                                    continue  # Essayer le fallback suivant

                                # OK - retourner ce fallback
                                import numpy as np
                                mean = float(np.mean(fallback_shifts))
                                std = float(np.std(fallback_shifts)) if len(fallback_shifts) > 1 else 0.0
                                logger.debug(f"Atom {atom_idx}: COMPOSITE R{radius} {try_solvent} (fallback n={len(fallback_shifts)})")
                                return mean, std, len(fallback_shifts), radius

                        # Aucun fallback n'a fonctionné (pas assez de données ou tous bimodaux)
                        # → Continuer vers "dernier recours" (low n) ci-dessous

                        # Dernier recours: utiliser solvant exact même avec peu de données (< min_samples)
                        # Seulement si aucun fallback n'a fonctionné
                        if exact_shifts:
                            # ÉTAPE 5: Vérifier bimodalité (même pour low n)
                            if len(exact_shifts) >= 10 and self._is_bimodal(exact_shifts, nucleus):
                                logger.debug(f"Atom {atom_idx}: COMPOSITE R{radius} {solvent_norm} (low n) BIMODAL → fallback ML")
                                continue  # Essayer radius suivant ou fallback ML

                            import numpy as np
                            mean = float(np.mean(exact_shifts))
                            std = float(np.std(exact_shifts)) if len(exact_shifts) > 1 else 0.0
                            logger.debug(f"Atom {atom_idx}: COMPOSITE R{radius} {solvent_norm} (low n={len(exact_shifts)})")
                            return mean, std, len(exact_shifts), radius

                        # Vraiment aucune donnée composite pour ce FP → essayer radius suivant
            except Exception as e:
                logger.debug(f"Composite lookup error: {e}")

        # ════════════════════════════════════════════════════════════════════
        # NIVEAU 3: MORGAN FINGERPRINT SEUL (fallback)
        # ════════════════════════════════════════════════════════════════════

        # Seuils minimum par radius
        MIN_COUNT_BY_RADIUS = {
            RADIUS_ULTRA_EXACT: MIN_SAMPLES_R6,
            RADIUS_HYPER_EXACT: MIN_SAMPLES_R5,
            RADIUS_VERY_EXACT: MIN_SAMPLES_R4,
            RADIUS_EXACT: MIN_SAMPLES_R3,
            RADIUS_FALLBACK: MIN_SAMPLES_R2
        }
        
        # Seuils de std pour fallback
        MAX_STD = MAX_STD_13C if nucleus == '13C' else MAX_STD_1H
        FALLBACK_STD = FALLBACK_MAX_STD_13C if nucleus == '13C' else FALLBACK_MAX_STD_1H
        MAX_VALUES_FOR_MEAN = 999999  # MODIFIE: Utiliser TOUTES les valeurs (pas de limite à 20)
        
        # Solvants "sans solvant"
        SOLVENTLESS_FALLBACK = ['Unreported', 'unknown', 'Unknown']

        # Hiérarchies de fallback par solvant (basées sur la similarité chimique)
        SOLVENT_FALLBACK_HIERARCHY = {
            'CDCl3': ['Acetone-d6', 'MeOD', 'DMSO-d6', 'Unreported'],
            'Acetone-d6': ['CDCl3', 'MeOD', 'DMSO-d6', 'Unreported'],
            'MeOD': ['Acetone-d6', 'CDCl3', 'DMSO-d6', 'Unreported'],
            'DMSO-d6': ['MeOD', 'Acetone-d6', 'CDCl3', 'Unreported'],
            'Unreported': ['CDCl3', 'DMSO-d6', 'MeOD', 'Acetone-d6'],
            'D2O': []  # Pas de fallback (trop différent chimiquement)
        }
        
        # Pré-calculer tous les fingerprints
        fp_by_radius = {r: get_atom_fingerprint(mol_with_h, atom_idx, r) 
                       for r in SUPPORTED_RADII}
        
        # Stocker le meilleur résultat (même avec std élevée)
        # Format: (mean, std, count, shifts_list) pour permettre vérification bimodale
        best_result = None
        best_radius = None
        best_solvent = None

        # FIX ERREUR 10: AGRÉGATION des solvants similaires au lieu de fallback séquentiel
        # ================================================================================
        # Stratégie: Pour chaque radius R6 → R2:
        #   1. Collecter les données du solvant demandé
        #   2. Si pas assez ou std trop élevée, agréger avec solvants similaires
        #   3. Calcul de moyenne pondérée par le nombre d'échantillons
        #   4. Retourner si std acceptable, sinon passer au radius suivant

        # Construire la liste des solvants compatibles (pour agrégation)
        fallback_solvents = SOLVENT_FALLBACK_HIERARCHY.get(solvent_norm, [])

        # Liste des radii à essayer (du plus précis au moins précis)
        radii_to_try = [RADIUS_ULTRA_EXACT, RADIUS_HYPER_EXACT, RADIUS_VERY_EXACT, RADIUS_EXACT, RADIUS_FALLBACK]

        for radius in radii_to_try:
            # Skip si min_radius_required > radius (ex: aromatiques nécessitent R3+)
            if radius < min_radius_required:
                continue

            db = FingerprintPredictor._db.get(radius, {})
            if not db:
                continue

            fp = fp_by_radius[radius]
            min_count = MIN_COUNT_BY_RADIUS.get(radius, 1)

            if fp not in db:
                continue

            # ÉTAPE 1: Essayer d'abord le solvant exact (+ variantes non-normalisées)
            # Fusion des variantes: ex. MeOD + "Methanol-D4 (MeOD-d4)"
            solvent_variants = get_solvent_variants(solvent_norm)
            exact_shifts = []
            for sv in solvent_variants:
                if sv in db[fp]:
                    exact_shifts.extend(db[fp][sv].get(nucleus, []))

            if len(exact_shifts) >= min_count:
                mean, std, count = self._compute_stats_no_check(exact_shifts, MAX_VALUES_FOR_MEAN, nucleus)

                # ÉTAPE 5: Vérifier bimodalité avant de retourner
                if self._is_bimodal(exact_shifts, nucleus):
                    logger.debug(f"Atom {atom_idx}: Morgan R{radius} {solvent_norm} BIMODAL → skip")
                    # Ne pas retourner, ne pas mémoriser comme best, continuer vers agrégation/radius suivant
                else:
                    # Si std acceptable avec solvant exact, retourner immédiatement
                    if std <= MAX_STD:
                        logger.debug(f"Atom {atom_idx}: {solvent_norm} R{radius}, n={count}, std={std:.2f}")
                        return mean, std, count, radius
                    # Garder comme meilleur candidat si acceptable pour fallback
                    if std <= FALLBACK_STD and (best_result is None or std < best_result[1]):
                        best_result = (mean, std, count, exact_shifts)  # Stocker shifts pour vérification finale
                        best_radius = radius
                        best_solvent = solvent_norm

            # ÉTAPE 2: FALLBACK HIÉRARCHIQUE - Essayer les solvants de fallback UN PAR UN
            # (dans l'ordre de priorité, PAS d'agrégation globale)
            if len(exact_shifts) < min_count or (len(exact_shifts) >= min_count and self._is_bimodal(exact_shifts, nucleus)):
                # Solvant exact insuffisant ou bimodal → essayer les fallbacks
                for solvent_candidate in fallback_solvents:
                    candidate_variants = get_solvent_variants(solvent_candidate)
                    candidate_shifts = []
                    for cv in candidate_variants:
                        if cv in db[fp]:
                            candidate_shifts.extend(db[fp][cv].get(nucleus, []))

                    # Si ce fallback a assez de données, le tester
                    if len(candidate_shifts) >= min_count:
                        mean, std, count = self._compute_stats_no_check(candidate_shifts, MAX_VALUES_FOR_MEAN, nucleus)

                        # ÉTAPE 5: Vérifier bimodalité
                        if self._is_bimodal(candidate_shifts, nucleus):
                            logger.debug(f"Atom {atom_idx}: Morgan R{radius} {solvent_candidate} BIMODAL → essayer fallback suivant")
                            continue  # Essayer le fallback suivant

                        # OK - std acceptable, retourner ce fallback
                        if std <= MAX_STD:
                            logger.debug(f"Atom {atom_idx}: {solvent_candidate} R{radius} (fallback), n={count}, std={std:.2f}")
                            return mean, std, count, radius

                        # Mémoriser comme candidat si acceptable
                        if std <= FALLBACK_STD and (best_result is None or std < best_result[1]):
                            best_result = (mean, std, count, candidate_shifts)
                            best_radius = radius
                            best_solvent = solvent_candidate

        # Si on a un résultat (même avec std élevée), le retourner
        # MAIS vérifier d'abord la bimodalité (ÉTAPE 5)
        if best_result is not None:
            mean, std, count, shifts_list = best_result

            # Vérifier bimodalité sur le meilleur résultat
            if len(shifts_list) >= 10 and self._is_bimodal(shifts_list, nucleus):
                logger.debug(f"Atom {atom_idx}: Best result R{best_radius} from {best_solvent} is BIMODAL → ML")
                return None  # Fallback ML

            logger.debug(f"Atom {atom_idx}: Using best result R{best_radius} from {best_solvent}")
            return (mean, std, count, best_radius)

        # Pas de match → directement ML
        logger.debug(f"Atom {atom_idx}: No match in {solvent_norm} (R6/R5/R4/R3/R2) → ML")
        return None
    
    def _search_any_solvent(self, mol_with_h, atom_idx: int, nucleus: str, 
                          min_count_by_radius: dict, max_values: int, exclude_solventless: bool = False) -> Optional[Tuple[float, float, int, int]]:
        """
        Recherche avec n'importe quel solvant (optionnellement exclure solventless).
        
        OPTIMISATION: Agrège TOUTES les valeurs des solvants concernés au lieu
        de retourner au premier match. Cela donne de meilleures statistiques.
        """
        import numpy as np
        
        SOLVENTLESS_FALLBACK = ['Unreported', 'unknown', 'Unknown']
        
        for radius in SUPPORTED_RADII:
            db = FingerprintPredictor._db.get(radius, {})
            if not db:
                continue
                
            fp = get_atom_fingerprint(mol_with_h, atom_idx, radius)
            
            if fp in db:
                # Agrégation: collecter TOUTES les valeurs des solvants pertinents
                all_shifts = []
                for solvent, nuc_dict in db[fp].items():
                    # Exclure les solvants "sans solvant" si demandé
                    if exclude_solventless and solvent in SOLVENTLESS_FALLBACK:
                        continue
                    
                    shifts = nuc_dict.get(nucleus, [])
                    if shifts:
                        all_shifts.extend(shifts)
                
                # Vérifier si on a assez de données
                if len(all_shifts) >= min_count_by_radius.get(radius, 1):
                    mean, std, count = self._compute_stats(all_shifts, max_values, nucleus)
                    logger.debug(f"Atom {atom_idx}: Any solvent R{radius}, n={count}, std={std:.2f}")
                    return mean, std, count, radius
        
        return None
    
    def _select_best_solventless(self, shifts: list, max_count: int) -> list:
        """
        Sélectionne les meilleures valeurs pour données solventless.
        
        Stratégie robuste avec clustering par densité:
        1. Si n ≤ max_count: prendre toutes
        2. Filtrage IQR pour éliminer outliers extrêmes
        3. Clustering par histogramme pour trouver le groupe principal
        4. Sélectionner max_count valeurs du cluster le plus dense
        """
        import numpy as np
        
        if len(shifts) <= max_count:
            return shifts
        
        shifts_array = np.array(shifts)
        
        # Étape 1: IQR filtering pour éliminer outliers
        q1, q3 = np.percentile(shifts_array, [25, 75])
        iqr = q3 - q1
        
        if iqr > 0:
            lower_bound = q1 - 2.0 * iqr  # Plus permissif pour solventless
            upper_bound = q3 + 2.0 * iqr
            filtered = shifts_array[(shifts_array >= lower_bound) & (shifts_array <= upper_bound)]
        else:
            filtered = shifts_array
        
        if len(filtered) <= max_count:
            return filtered.tolist()
        
        # Étape 2: Clustering par histogramme adaptatif
        n_bins = min(20, max(5, len(filtered) // 10))
        hist, bin_edges = np.histogram(filtered, bins=n_bins)
        
        # Trouver le bin le plus dense
        densest_bin_idx = np.argmax(hist)
        densest_bin_center = (bin_edges[densest_bin_idx] + bin_edges[densest_bin_idx + 1]) / 2
        
        # Étape 3: Sélectionner les max_count valeurs les plus proches du centre dense
        distances = np.abs(filtered - densest_bin_center)
        nearest_indices = np.argsort(distances)[:max_count]
        
        return filtered[nearest_indices].tolist()
    
    def _remove_outliers_iqr(self, shifts_list: list, iqr_factor: float = 1.5) -> list:
        """
        Filtre les outliers majeurs et mineurs avec la méthode IQR (Interquartile Range).
        
        Args:
            shifts_list: Liste des shifts chimiques
            iqr_factor: Facteur multiplicateur pour IQR (1.5 = standard, détecte outliers majeurs et mineurs)
        
        Returns:
            Liste des shifts sans outliers
        """
        import numpy as np
        
        if len(shifts_list) < 4:
            # Pas assez de données pour filtrage statistique
            return shifts_list
        
        shifts_array = np.array(shifts_list)
        
        # Calcul des quartiles
        q1 = np.percentile(shifts_array, 25)
        q3 = np.percentile(shifts_array, 75)
        iqr = q3 - q1
        
        # Définir les bornes
        lower_bound = q1 - iqr_factor * iqr
        upper_bound = q3 + iqr_factor * iqr
        
        # Filtrer les valeurs dans l'intervalle acceptable
        filtered = shifts_array[(shifts_array >= lower_bound) & (shifts_array <= upper_bound)]
        
        return filtered.tolist()
    
    def _filter_coherent_shifts(self, shifts_list: list, nucleus: str) -> list:
        """
        Filtre les shifts pour ne garder que ceux qui contribuent à un std acceptable.
        Utilise une approche par percentiles pour trouver le cluster principal.
        
        Args:
            shifts_list: Liste des shifts après filtrage IQR
            nucleus: '13C' ou '1H'
        
        Returns:
            Liste des shifts cohérents (std acceptable)
        """
        import numpy as np
        
        # Seuils de std maximum
        max_std = MAX_STD_13C if nucleus == '13C' else MAX_STD_1H
        
        # Si déjà peu de données, ne pas filtrer davantage
        if len(shifts_list) <= 3:
            return shifts_list
        
        shifts_array = np.array(shifts_list)
        
        # Si le std actuel est déjà acceptable, retourner tout
        current_std = np.std(shifts_array)
        if current_std <= max_std:
            return shifts_list
        
        # Stratégie: Trouver le percentile range le plus large avec std acceptable
        # On teste des intervalles de percentiles de plus en plus étroits
        best_subset = []
        best_size = 0
        
        # Percentile ranges à tester (du plus large au plus étroit)
        percentile_ranges = [
            (0, 100),    # Tout
            (5, 95),     # Enlever 5% de chaque côté
            (10, 90),    # Enlever 10% de chaque côté
            (15, 85),
            (20, 80),
            (25, 75),    # Quartiles
            (30, 70),
            (35, 65),
            (40, 60)
        ]
        
        for pct_low, pct_high in percentile_ranges:
            low_val = np.percentile(shifts_array, pct_low)
            high_val = np.percentile(shifts_array, pct_high)
            
            # Garder les valeurs dans cet intervalle
            mask = (shifts_array >= low_val) & (shifts_array <= high_val)
            subset = shifts_array[mask]
            
            if len(subset) < 2:
                continue
            
            subset_std = np.std(subset)
            
            if subset_std <= max_std and len(subset) > best_size:
                best_subset = subset.tolist()
                best_size = len(subset)
                break  # On garde le plus large intervalle acceptable
        
        # Si aucun intervalle acceptable, garder les valeurs proches de la médiane
        if not best_subset:
            median = np.median(shifts_array)
            distances = np.abs(shifts_array - median)
            # Garder au moins 3 valeurs ou 25% des données
            n_keep = max(3, len(shifts_array) // 4)
            closest_indices = np.argsort(distances)[:n_keep]
            best_subset = shifts_array[closest_indices].tolist()
        
        return best_subset
    
    def _compute_stats(self, shifts_list: list, max_values: int, nucleus: str = '13C') -> Tuple[float, float, int]:
        """Calcule mean/std avec filtrage automatique des outliers (méthode IQR) et filtrage de cohérence."""
        import numpy as np
        
        original_count = len(shifts_list)
        
        # Étape 1: Filtrage des outliers avec méthode IQR
        cleaned_shifts = self._remove_outliers_iqr(shifts_list)
        
        # Étape 2: Filtrage de cohérence (garder uniquement shifts donnant std acceptable)
        coherent_shifts = self._filter_coherent_shifts(cleaned_shifts, nucleus)
        
        # Étape 3: Sélection des meilleurs si trop nombreux
        if len(coherent_shifts) > max_values:
            selected = self._select_best_n(coherent_shifts, max_values)
        else:
            selected = coherent_shifts
        
        return float(np.mean(selected)), float(np.std(selected)), original_count
    
    def _compute_stats_no_check(self, shifts_list: list, max_values: int, nucleus: str = '13C') -> Tuple[float, float, int]:
        """
        Calcule mean/std SANS validation de seuil std, AVEC filtrage outliers ET cohérence.
        Utilisé pour fallback radius: permet d'évaluer la qualité avant d'accepter/rejeter.
        """
        import numpy as np

        original_count = len(shifts_list)

        # Étape 1: Filtrage des outliers avec méthode IQR
        cleaned_shifts = self._remove_outliers_iqr(shifts_list)

        # Étape 2: Filtrage de cohérence (garder uniquement shifts donnant std acceptable)
        coherent_shifts = self._filter_coherent_shifts(cleaned_shifts, nucleus)

        # Étape 3: Sélection des meilleurs si trop nombreux
        if len(coherent_shifts) > max_values:
            selected = self._select_best_n(coherent_shifts, max_values)
        else:
            selected = coherent_shifts

        return float(np.mean(selected)), float(np.std(selected)), original_count
    
    def _compute_stats_with_original_count(self, selected_shifts: list, original_count: int) -> Tuple[float, float, int]:
        """Calcule stats sur shifts sélectionnés mais retourne le count original."""
        import numpy as np
        return float(np.mean(selected_shifts)), float(np.std(selected_shifts)), original_count

    def _select_best_n(self, shifts_list: list, n: int) -> list:
        """Sélectionne les n valeurs avec std minimal (fenêtre glissante)."""
        import numpy as np
        
        if len(shifts_list) <= n:
            return shifts_list
        
        # Trier les shifts
        sorted_shifts = sorted(shifts_list)
        
        # Chercher la fenêtre de taille n avec le std minimal
        best_std = float('inf')
        best_start = 0
        
        for i in range(len(sorted_shifts) - n + 1):
            window = sorted_shifts[i:i+n]
            window_std = np.std(window)
            if window_std < best_std:
                best_std = window_std
                best_start = i
        
        return sorted_shifts[best_start:best_start+n]
    
    def get_environment(self, mol_with_h, atom_idx: int) -> str:
        """Retourne la description de l'environnement atomique."""
        return get_atom_env_string(mol_with_h, atom_idx)
    
    def is_available(self) -> bool:
        """Vérifie si la base de données est chargée."""
        return FingerprintPredictor._db is not None and len(FingerprintPredictor._db) > 0


# =============================================================================
# CONVENIENCE FUNCTION
# =============================================================================
def get_predictor() -> FingerprintPredictor:
    """Retourne l'instance singleton du prédicteur."""
    return FingerprintPredictor()


# =============================================================================
# TEST
# =============================================================================
if __name__ == '__main__':
    import sys
    
    smiles = "CCO"  # Éthanol
    solvent = "CDCl3"
    
    if len(sys.argv) >= 2:
        smiles = sys.argv[1]
    if len(sys.argv) >= 3:
        solvent = sys.argv[2]
    
    print(f"Testing FingerprintPredictor with: {smiles} in {solvent}")
    print("=" * 60)
    
    predictor = FingerprintPredictor()
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES: {smiles}")
        sys.exit(1)
    
    mol_with_h = Chem.AddHs(mol)
    
    print(f"\nPredictions (multi-radius hierarchy):")
    print("-" * 60)
    
    found_count = 0
    not_found_count = 0
    
    for atom_idx in range(mol_with_h.GetNumAtoms()):
        atom = mol_with_h.GetAtomWithIdx(atom_idx)
        symbol = atom.GetSymbol()
        
        if symbol == 'C':
            nucleus = '13C'
        elif symbol == 'H':
            nucleus = '1H'
        else:
            continue
        
        env = predictor.get_environment(mol_with_h, atom_idx)
        result = predictor.predict_atom(mol_with_h, atom_idx, solvent, nucleus)
        
        if result:
            shift, std, count, radius_used = result
            print(f"  {atom_idx:3d} | {nucleus:3s} | {env:30s} | {shift:7.2f} ± {std:5.2f} (n={count}, R{radius_used})")
            found_count += 1
        else:
            print(f"  {atom_idx:3d} | {nucleus:3s} | {env:30s} | NOT FOUND (use ML fallback)")
            not_found_count += 1
    
    print("-" * 60)
    print(f"Found: {found_count}, Not found (need ML): {not_found_count}")
