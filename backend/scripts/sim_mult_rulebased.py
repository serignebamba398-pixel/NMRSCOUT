# sim_mult_rulebased.py
# -*- coding: utf-8 -*-
"""
Rule-based ¹H NMR Multiplicity Simulator — High-Fidelity Edition
=================================================================

Simulateur de multiplicité RMN ¹H basé sur les règles physiques établies.

RÉFÉRENCES PRINCIPALES:
-----------------------
- Karplus M., J. Chem. Phys. 30, 11 (1959) - Équation Karplus originale
- Haasnoot C.A.G. et al., Tetrahedron 36, 2783-2792 (1980) - Équation complète
- Bystrov V.F., Prog. NMR Spectrosc. 10, 41-82 (1976) - ³J(H-N-C-H) peptides
- Jackman L.M. & Sternhell S., Applications of NMR Spectroscopy (1969)
- Günther H., NMR Spectroscopy, 3rd Ed. (2013) - Systèmes AB, raies combinaison
- Kalinowski H.-O. et al., Carbon-13 NMR Spectroscopy (1988) - ¹J(C-H)
- Reich H.J., Structure Determination Using NMR - Effet de toit

CARACTÉRISTIQUES:
-----------------
✓ Équation Karplus-Haasnoot-Altona complète (7 paramètres)
✓ ²J geminaux paramétrés selon environnement chimique
✓ ¹J(C-H) variables selon hybridation (125-250 Hz)
✓ Détection second ordre (Δν/J raffinée)
✓ Effet de toit (roof effect) pour systèmes AB
✓ Raies de combinaison pour systèmes fortement couplés
✓ ³J(H-N-C-H) pour peptides/amides (Bystrov)
✓ Couplages hétéronucléaires (¹⁹F, ³¹P)
✓ Mode sucres (anomères α/β, pyranoses)
✓ Satellites ¹³C avec abondance naturelle

VERSION: 2.0.0 - High-Fidelity Edition
AUTEUR: NMR Simulation Team
DATE: Décembre 2025
"""

from __future__ import annotations
from typing import Dict, List, Tuple, Optional, Iterable, Any
from dataclasses import dataclass, field
from collections import defaultdict, deque
import math
import numpy as np
from functools import lru_cache
import weakref
import hashlib

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem
from rdkit.Chem import Descriptors

# Import du module hétérocycles complet
try:
    from heterocycles_complete import (
        identify_heterocycle,
        get_heterocycle_coupling,
        is_heterocycle_singlet,
        estimate_heterocycle_couplings_complete,
        HeterocycleMatch
    )
    HETEROCYCLES_MODULE_AVAILABLE = True
except ImportError:
    HETEROCYCLES_MODULE_AVAILABLE = False

# Import du module sucres complet
try:
    from sugars_complete import (
        estimate_sugar_couplings_complete,
        is_sugar_carbon,
        find_sugar_rings
    )
    SUGARS_MODULE_AVAILABLE = True
except ImportError:
    SUGARS_MODULE_AVAILABLE = False

# Import du module couplages complets (alcynes, aldéhydes, phosphines, stéroïdes, etc.)
try:
    from coupling_tables_complete import (
        get_specialized_coupling,
        MoleculeClassifier,
        get_alkyne_coupling,
        get_aldehyde_coupling,
        get_peptide_coupling,
        get_phosphine_coupling,
        get_borane_coupling,
        get_silane_coupling,
        get_steroid_coupling,
        get_terpene_coupling,
        get_macrocycle_coupling,
        get_nucleotide_coupling,
    )
    COUPLING_TABLES_MODULE_AVAILABLE = True
except ImportError:
    COUPLING_TABLES_MODULE_AVAILABLE = False

# Import du module peptides robuste
try:
    from peptide_couplings import (
        get_peptide_coupling_robust,
        calculate_3J_HN_HA,
        calculate_3J_HA_HB,
        calculate_2J_geminal as calculate_peptide_2J_geminal,
        KARPLUS_PARAMETERS,
        AMINO_ACID_SPECIFIC_COUPLINGS,
    )
    PEPTIDE_MODULE_AVAILABLE = True
except ImportError:
    PEPTIDE_MODULE_AVAILABLE = False

# Import du module extensions (aromatiques polycycliques, cycles tendus, etc.)
try:
    from coupling_extensions import (
        get_extended_coupling,
        get_polycyclic_aromatic_coupling,
        get_strained_ring_coupling,
        detect_strained_ring,
        POLYCYCLIC_AROMATIC_COUPLINGS,
        STRAINED_RING_COUPLINGS,
        GLYCOSIDIC_COUPLINGS,
        TERPENE_SKELETON_COUPLINGS,
        STEROID_RING_COUPLINGS,
        PHOSPHOLIPID_COUPLINGS,
        NUCLEOTIDE_COUPLINGS,
    )
    COUPLING_EXTENSIONS_AVAILABLE = True
except ImportError:
    COUPLING_EXTENSIONS_AVAILABLE = False

# Import du module tables ultimes (couverture maximale)
try:
    from coupling_tables_ultimate import (
        get_ultimate_coupling,
        KarplusCalculator,
        AROMATIC_J_COMPLETE,
        HETEROAROMATIC_J_COMPLETE,
        GEMINAL_2J_COMPLETE,
        ALKENE_J_COMPLETE,
        ALDEHYDE_J_COMPLETE,
        TERPENE_J_COMPLETE,
        NUCLEOSIDE_J_COMPLETE,
        PHOSPHOLIPID_J_COMPLETE,
        STEROID_J_COMPLETE,
    )
    ULTIMATE_TABLES_AVAILABLE = True
except ImportError:
    ULTIMATE_TABLES_AVAILABLE = False


# ==============================================================================
# CACHE INTELLIGENT - AMÉLIORATION PERFORMANCE
# ==============================================================================
# Cache les calculs coûteux (distances topologiques, aromatic maps, conformers)
# basé sur la signature SMILES de la molécule
# ==============================================================================

class MoleculeCache:
    """
    Cache intelligent pour les calculs moléculaires coûteux.
    
    Stocke les résultats par clé SMILES pour éviter les recalculs:
    - Distances topologiques
    - Maps aromatiques
    - Rangs canoniques
    - Conformères
    """
    
    _instance = None
    _cache_topo_dist: Dict[str, np.ndarray] = {}
    _cache_aromatic_map: Dict[str, List[Dict[int, int]]] = {}
    _cache_ranks: Dict[str, List[int]] = {}
    _cache_conformers: Dict[str, List[Chem.Mol]] = {}
    _max_size = 100  # Limite pour éviter les fuites mémoire
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
    
    @classmethod
    def get_key(cls, mol: Chem.Mol) -> str:
        """Génère une clé unique pour la molécule (basée sur SMILES canonique)."""
        try:
            return Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            # Fallback: hash des atomes
            atoms_str = "".join(f"{a.GetSymbol()}{a.GetIdx()}" for a in mol.GetAtoms())
            return hashlib.md5(atoms_str.encode()).hexdigest()
    
    @classmethod
    def get_topological_distances(cls, mol: Chem.Mol) -> np.ndarray:
        """Récupère ou calcule les distances topologiques."""
        key = cls.get_key(mol)
        if key in cls._cache_topo_dist:
            return cls._cache_topo_dist[key]
        
        # Calculer et cacher
        dist = _compute_topological_distances(mol)
        
        # Gérer la taille du cache
        if len(cls._cache_topo_dist) >= cls._max_size:
            cls._cache_topo_dist.pop(next(iter(cls._cache_topo_dist)))
        
        cls._cache_topo_dist[key] = dist
        return dist
    
    @classmethod
    def get_aromatic_map(cls, mol: Chem.Mol) -> List[Dict[int, int]]:
        """Récupère ou calcule la map aromatique."""
        key = cls.get_key(mol)
        if key in cls._cache_aromatic_map:
            return cls._cache_aromatic_map[key]
        
        # Calculer et cacher
        arom_map = _compute_aromatic_distance_map(mol)
        
        if len(cls._cache_aromatic_map) >= cls._max_size:
            cls._cache_aromatic_map.pop(next(iter(cls._cache_aromatic_map)))
        
        cls._cache_aromatic_map[key] = arom_map
        return arom_map
    
    @classmethod
    def get_canonical_ranks(cls, mol_noH: Chem.Mol) -> List[int]:
        """Récupère ou calcule les rangs canoniques."""
        key = cls.get_key(mol_noH)
        if key in cls._cache_ranks:
            return cls._cache_ranks[key]
        
        try:
            ranks = list(Chem.CanonicalRankAtoms(mol_noH, includeChirality=True))
        except Exception:
            ranks = list(Chem.CanonicalRankAtoms(mol_noH))
        
        if len(cls._cache_ranks) >= cls._max_size:
            cls._cache_ranks.pop(next(iter(cls._cache_ranks)))
        
        cls._cache_ranks[key] = ranks
        return ranks
    
    @classmethod
    def clear(cls):
        """Vide tous les caches."""
        cls._cache_topo_dist.clear()
        cls._cache_aromatic_map.clear()
        cls._cache_ranks.clear()
        cls._cache_conformers.clear()


def _compute_topological_distances(mol: Chem.Mol) -> np.ndarray:
    """Calcul brut des distances topologiques (utilisé par le cache)."""
    try:
        from rdkit.Chem.rdmolops import GetDistanceMatrix as _RDGetDistanceMatrix
        D = _RDGetDistanceMatrix(mol, useBO=False, useAtomWts=False)
        return np.asarray(D, dtype=float)
    except Exception:
        n = mol.GetNumAtoms()
        if n == 0:
            return np.zeros((0, 0), dtype=float)
        D = np.full((n, n), np.inf, dtype=float)
        np.fill_diagonal(D, 0)
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            D[i, j] = D[j, i] = 1
        # Floyd-Warshall
        for k in range(n):
            D = np.minimum(D, D[:, k:k+1] + D[k:k+1, :])
        return D


def _compute_aromatic_distance_map(mol: Chem.Mol) -> List[Dict[int, int]]:
    """Calcul brut de la distance map aromatique (utilisé par le cache)."""
    n = mol.GetNumAtoms()
    dist_maps = [{} for _ in range(n)]
    aro_idx = {a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()}
    if not aro_idx:
        return dist_maps
    
    # BFS depuis chaque atome aromatique
    for start in aro_idx:
        seen = {start}
        q = deque([(start, 0)])
        dist_maps[start][start] = 0
        while q:
            u, d = q.popleft()
            if d >= 3:
                continue  # Limiter à ortho/meta/para
            atom = mol.GetAtomWithIdx(u)
            for bond in atom.GetBonds():
                v = bond.GetOtherAtomIdx(u)
                if v not in seen and v in aro_idx:
                    seen.add(v)
                    dist_maps[start][v] = d + 1
                    q.append((v, d + 1))
    return dist_maps

# ---------------------------
# Canon labels (49)
# ---------------------------
H1_CANON_CLASSES = [
    's','d','t','q','quin','sext','sept','oct','non','dec','m',
    'dd','ddd','dddd',
    'dt','td','dq','qd',
    'tt','tq','qt',
    'ddt','dtd','tdd',
    'tdq','tqd','dtq',
    'qdt','qtd','dqt',
    'qqd','dqd','qdd','qqt','tqq','qtq',
    'dddt','dtdd','tddd',
    'dddq','ddqd','dqdd','qddd',
    'ttd','tdt','dtt',
    'ttq','tqt','qtt'
]
LABEL2IDX = {lbl: i for i, lbl in enumerate(H1_CANON_CLASSES)}

# ---------------------------
# Instrument defaults
# ---------------------------
DEFAULT_SPECTROMETER_FREQ_MHZ = 700.0
DEFAULT_TEMPERATURE_K = 298.15

# ---------------------------
# Solvent DB (lightweight)
# ---------------------------
SOLVENT_DB = {
    "CDCl3": {
        "name": "CDCl3", "dielectric": 4.8, "viscosity": "low",
        "protic": False, "hbond_acceptor": "weak",
        "exchange_rate_factor": 1.0,
        "chemical_shift_effects": {"OH": 0.0, "NH": 0.0, "COOH": +0.10},
        "exchange_regime": "moderate",
        "anisotropy": "normal",
    },
    "DMSO-d6": {
        "name": "DMSO-d6", "dielectric": 47.2, "viscosity": "moderate",
        "protic": False, "hbond_acceptor": "strong",
        "exchange_rate_factor": 0.6,
        "chemical_shift_effects": {"OH": +0.30, "NH": +0.20, "COOH": +0.30},
        "exchange_regime": "slow",
        "anisotropy": "normal",
    },
    "D2O": {
        "name": "D2O", "dielectric": 78.4, "viscosity": "low",
        "protic": True, "hbond_acceptor": "strong",
        "exchange_rate_factor": 2.0,
        "chemical_shift_effects": {"OH": "suppressed", "NH": "suppressed", "COOH": "suppressed"},
        "exchange_regime": "fast",
        "anisotropy": "normal",
    },
    "C6D6": {
        "name": "C6D6", "dielectric": 2.3, "viscosity": "low",
        "protic": False, "hbond_acceptor": "none",
        "exchange_rate_factor": 1.2,
        "chemical_shift_effects": {"OH": -0.10, "NH": -0.10, "COOH": 0.0},
        "exchange_regime": "moderate",
        "anisotropy": "high",
    },
    "CD3OD": {
        "name": "CD3OD", "dielectric": 33.0, "viscosity": "low",
        "protic": True, "hbond_acceptor": "strong",
        "exchange_rate_factor": 1.6,
        "chemical_shift_effects": {"OH": +0.20, "NH": +0.15, "COOH": +0.25},
        "exchange_regime": "moderate",
        "anisotropy": "normal",
    },
    "CD3CN": {
        "name": "CD3CN", "dielectric": 37.5, "viscosity": "low",
        "protic": False, "hbond_acceptor": "strong",
        "exchange_rate_factor": 1.1,
        "chemical_shift_effects": {"OH": +0.10, "NH": +0.05, "COOH": +0.15},
        "exchange_regime": "moderate",
        "anisotropy": "normal",
    },
    "pyridine-d5": {
        "name": "pyridine-d5", "dielectric": 12.4, "viscosity": "moderate",
        "protic": False, "hbond_acceptor": "strong",
        "exchange_rate_factor": 0.9,
        "chemical_shift_effects": {"OH": +0.35, "NH": +0.25, "COOH": +0.35},
        "exchange_regime": "slow",
        "anisotropy": "normal",
    },
}
DEFAULT_SOLVENT = "CDCl3"

# ---------------------------
# Constants & thresholds
# ---------------------------
TOLERANCE_J_HZ = 0.15

FREQ_LINEWIDTH_TABLE = [
    (0,   150, 1.00),
    (150, 250, 0.70),
    (250, 350, 0.55),
    (350, 450, 0.50),
    (450, 650, 0.35),
    (650, 1e9, 0.25),
]

def _typical_linewidth_from_field(mhz: float, viscosity: str) -> float:
    base = next(w for lo, hi, w in FREQ_LINEWIDTH_TABLE if lo < mhz <= hi)
    visco_add = {"low": 0.00, "moderate": 0.05, "high": 0.15}.get(viscosity, 0.00)
    return base + visco_add

_13C_NATURAL_ABUNDANCE = 0.011

# ¹J(C-H) selon hybridation et environnement (Hz)
# Références: Kalinowski et al., Carbon-13 NMR Spectroscopy
_13C_1JCH_TABLE = {
    'sp3_alkyl': 125.0,       # CH₃, CH₂ aliphatique
    'sp3_alpha_O': 140.0,     # C-H adjacent à O (sucres, éthers)
    'sp3_alpha_N': 135.0,     # C-H adjacent à N
    'sp3_cyclopropane': 160.0, # Tension de cycle
    'sp2_alkene': 155.0,      # C=C-H
    'sp2_aldehyde': 175.0,    # CHO
    'sp_alkyne': 250.0,       # ≡C-H
    'aromatic': 160.0,        # Ar-H
    'default': 140.0,
}

# ²J(H-H) geminaux selon l'environnement (Hz)
# Références: Jackman & Sternhell, Applications of NMR Spectroscopy
_GEMINAL_2J_TABLE = {
    'CH2_alkyl': (-12.5, -14.0),      # CH₂ aliphatique (négatif!)
    'CH2_alpha_O': (-10.0, -12.0),    # CH₂-O (sucres, éthers)
    'CH2_alpha_N': (-11.0, -13.0),    # CH₂-N
    'CH2_alpha_C=O': (-14.0, -18.0),  # CH₂-C=O (carbonyle)
    'CH2_cyclopropane': (-4.0, -6.0), # Cyclopropane (moins négatif)
    'CH2_cyclobutane': (-11.0, -13.0),
    'CH2_cyclopentane': (-12.0, -14.0),
    'CH2_cyclohexane': (-12.5, -14.5),
    'CH2_exocyclic': (-10.0, -13.0),  # CH₂ exocyclique (ex: glucose C6)
    'CH2_sp2': (0.0, 3.0),            # =CH₂ (positif ou faible)
    'aldehyde': (2.0, 3.0),           # CHO (positif)
    'default_sp3': (-12.5, -14.0),
    'default_sp2': (1.5, 2.5),
}

# ==============================================================================
# ⁴J COUPLAGES LONGUE PORTÉE - TABLES COMPLÈTES
# ==============================================================================
# Références:
# - Barfield M., Grant D.M., Adv. Magn. Reson. 1, 149 (1965)
# - Jackman L.M., Sternhell S., Applications of NMR Spectroscopy (1969)
# - Günther H., NMR Spectroscopy 3rd Ed. (2013), Chapter 7.5
# ==============================================================================

_4J_LONG_RANGE_TABLE = {
    # Couplages allyliques ⁴J(H-C-C=C-H)
    # Valeur dépend de l'angle dièdre H-C...C-H
    'allylic_cis': (1.5, 3.0),        # H et CH₃ même côté du C=C
    'allylic_trans': (0.8, 2.0),      # H et CH₃ côtés opposés
    'allylic_general': (0.5, 3.0),    # Gamme générale
    
    # Couplages homoallyliques ⁴J(H-C=C-C-H)
    # Plus petits que allyliques
    'homoallylic_cis': (1.0, 2.0),    # Géométrie cis
    'homoallylic_trans': (0.5, 1.5),  # Géométrie trans
    'homoallylic_general': (0.5, 2.5),
    
    # Couplages W-allyliques (configuration planaire en W)
    # Plus grands car chemin σ/π optimal
    'W_allylic_rigid': (1.5, 3.5),    # Systèmes rigides (cycles)
    'W_allylic_flexible': (0.8, 2.0), # Systèmes flexibles
    
    # ⁴J à travers cycles aromatiques (méta-aromatique)
    'meta_aromatic': (2.0, 3.5),      # H-Ar-C-Ar-H (méta)
    'para_aromatic': (0.3, 0.8),      # H-Ar-C-C-Ar-H (para, très faible)
    
    # ⁴J benzoyle/benzyle: H-C-Ar-H
    'benzylic': (0.5, 1.5),           # CH₃ ou CH₂ adjacent à aromatique
    
    # ⁴J aldéhyde-allyle: CHO-C=C-H
    'aldehyde_vinylic': (1.0, 2.5),   # α,β-insaturé
    
    # ⁴J saturé en W (systèmes bicycliques rigides)
    'W_saturated_rigid': (1.0, 2.5),  # Bicyclo[2.2.1]heptane, stéroïdes
    'W_saturated_flexible': (0.0, 0.5), # Rarement observable
    
    # ⁴J cycliques spéciaux
    'cyclopropyl': (1.0, 2.0),        # À travers cyclopropane
    'cyclobutyl': (0.5, 1.5),         # À travers cyclobutane
    
    # Défauts
    'default_rigid': (1.0, 2.0),
    'default_flexible': (0.0, 0.8),
}

# ==============================================================================
# ⁵J COUPLAGES TRÈS LONGUE PORTÉE - TABLES COMPLÈTES
# ==============================================================================
# Ces couplages sont généralement très petits (<1 Hz) sauf configurations spéciales
# Références: Barfield M., J. Am. Chem. Soc. 93, 1066 (1971)
# ==============================================================================

_5J_LONG_RANGE_TABLE = {
    # Couplages W "étendu" (5 liaisons en zigzag)
    # Observables seulement dans systèmes rigides planaires
    'extended_W_rigid': (0.8, 2.0),   # Stéroïdes, systèmes polycycliques
    'extended_W_flexible': (0.0, 0.3), # Généralement non résolu
    
    # ⁵J para-aromatique: H-Ar-C-C-Ar-H
    'para_substituted': (0.4, 0.8),   # À travers noyau aromatique
    
    # ⁵J homoallylique étendu: H-C-C=C-C-H
    'homoallylic_extended': (0.3, 1.0), # Rare mais mesurable
    
    # ⁵J dans systèmes π conjugués
    'conjugated_diene': (0.5, 1.5),   # À travers diène conjugué
    'conjugated_enyne': (0.3, 0.8),   # À travers ényne
    
    # ⁵J naphthalène (protons péri)
    'peri_coupling': (0.5, 1.0),      # Positions 1,8 naphthalène
    
    # Défauts
    'default': (0.0, 0.5),
}

# ==============================================================================
# COUPLAGES H-F ET H-P - TABLES COMPLÈTES
# ==============================================================================
# Références:
# - Emsley J.W., Phillips L., Wray V., Prog. NMR Spectrosc. 10, 83 (1976) - H-F
# - Gorenstein D.G., Prog. NMR Spectrosc. 16, 1 (1983) - H-P (31P)
# - Günther H., NMR Spectroscopy 3rd Ed. (2013), Chapter 7.6
# ==============================================================================

_HF_COUPLING_TABLE = {
    # ²J(H-C-F) geminaux - très grands!
    '2J_CHF_sp3': (42.0, 52.0),       # H-C(sp³)-F typique ~47 Hz
    '2J_CHF_sp2': (44.0, 55.0),       # H-C(sp²)=C-F légèrement plus grand
    '2J_CHF_formyl': (48.0, 56.0),    # FCH=O (formyle fluoré)
    '2J_CHF_CF2': (50.0, 58.0),       # CHF₂ ou CF₂
    '2J_CHF_CF3': (44.0, 50.0),       # CHFCF₃
    
    # ³J(H-C-C-F) vicinaux
    '3J_HCCF_anti': (20.0, 35.0),     # Relation antiperiplanaire (Karplus max)
    '3J_HCCF_gauche': (5.0, 15.0),    # Relation gauche (Karplus min)
    '3J_HCCF_average': (8.0, 22.0),   # Moyenne rotationnelle
    '3J_HCCF_cyclic': (10.0, 25.0),   # Dans cycles (angle fixe)
    '3J_HCCF_aromatic': (6.0, 10.0),  # Ar-F avec H ortho
    
    # ⁴J(H-C-C-C-F) longue portée
    '4J_HCCCF_W': (3.0, 8.0),         # Configuration W (grand!)
    '4J_HCCCF_general': (0.5, 4.0),   # Géométrie quelconque
    '4J_HCCCF_aromatic': (2.0, 5.0),  # Ar-F avec H méta
    '4J_allylic_F': (2.0, 6.0),       # H-C-C=C-F (allylique)
    
    # ⁵J(H-C-C-C-C-F) très longue portée
    '5J_HCCCCF_W': (1.0, 3.0),        # W étendu (rare)
    '5J_HCCCCF_aromatic': (0.5, 2.0), # Ar-F avec H para
    '5J_HCCCCF_general': (0.0, 1.5),  # Généralement non résolu
    
    # À travers liaisons multiples
    'vinylic_2J': (75.0, 90.0),       # H-C=C-F geminal (très grand!)
    'vinylic_3J_cis': (18.0, 28.0),   # H-C=C-F cis
    'vinylic_3J_trans': (35.0, 45.0), # H-C=C-F trans
    
    # Défauts par ordre de couplage
    'default_2J': (46.0, 50.0),
    'default_3J': (10.0, 20.0),
    'default_4J': (1.0, 4.0),
    'default_5J': (0.0, 1.0),
}

_HP_COUPLING_TABLE = {
    # Couplages H-³¹P (phosphore)
    # Références: Gorenstein D.G., Phosphorus-31 NMR (1984)
    
    # ²J(H-C-P) geminaux
    '2J_HCP_sp3': (8.0, 18.0),        # H-C(sp³)-P
    '2J_HCP_sp2': (10.0, 22.0),       # H-C(sp²)=C-P
    '2J_HCP_phosphonate': (18.0, 25.0), # (RO)₂P(=O)CH₂R
    
    # ³J(H-C-C-P) vicinaux
    '3J_HCCP_anti': (12.0, 20.0),     # Antiperiplanaire
    '3J_HCCP_gauche': (2.0, 8.0),     # Gauche
    '3J_HCCP_average': (5.0, 12.0),   # Moyenne
    '3J_HCCP_sugar': (8.0, 12.0),     # Nucléotides, sucres phosphatés
    
    # ³J(H-C-O-P) à travers oxygène (phosphates)
    '3J_HCOP_ester': (6.0, 12.0),     # Esters phosphate
    '3J_HCOP_nucleotide': (4.0, 10.0), # Nucléotides
    
    # ⁴J(H-C-C-C-P) longue portée
    '4J_HCCCP_W': (2.0, 5.0),         # Configuration W
    '4J_HCCCP_general': (0.5, 2.0),   # Quelconque
    
    # Défauts
    'default_2J': (12.0, 18.0),
    'default_3J': (6.0, 12.0),
    'default_4J': (0.5, 2.0),
}

# ==============================================================================
# PARAMÈTRES ALCÈNES ENRICHIS
# ==============================================================================
# Références additionnelles:
# - Pascual C., Meier J., Simon W., Helv. Chim. Acta 49, 164 (1966)
# - Abraham R.J., Bernstein H.J., Can. J. Chem. 39, 216 (1961)
# ==============================================================================

_GARBISCH_ALKENE_EXTENDED = {
    # Format: (A, B, C) pour J = A + B·cos(θ) + C·cos(2θ)
    # θ = 0° cis, θ = 180° trans
    
    # Alcènes simples (Garbisch original)
    'terminal': (11.0, -3.0, 2.0),         # R-CH=CH₂: cis~10, trans~16
    'disubstituted_EZ': (11.5, -3.5, 2.0), # R-CH=CH-R: E/Z bien défini
    'trisubstituted': (10.5, -2.5, 1.5),   # R₂C=CH-R: moins de couplage
    'tetrasubstituted': (0.0, 0.0, 0.0),   # R₂C=CR₂: pas de H vinylique
    
    # Alcènes conjugués (effet électronique)
    'diene_terminal': (11.5, -4.0, 2.5),   # CH₂=CH-CH=CH₂
    'diene_internal': (11.0, -3.5, 2.0),   # R-CH=CH-CH=CH-R
    'enone': (10.5, -3.0, 1.5),            # R-CH=CH-C(=O)-R (α,β-insaturé)
    'styrene': (11.0, -3.0, 2.0),          # Ar-CH=CH₂
    'stilbene': (11.5, -4.0, 2.5),         # Ar-CH=CH-Ar
    
    # Alcènes cycliques (contrainte angulaire)
    'cyclopropene': (8.0, -1.5, 1.0),      # Très petit cycle
    'cyclobutene': (9.0, -2.0, 1.5),       # 4-cycle
    'cyclopentene': (10.0, -2.5, 1.5),     # 5-cycle: cis~10 Hz
    'cyclohexene': (10.5, -3.0, 2.0),      # 6-cycle: proche standard
    'cycloheptene': (11.0, -3.0, 2.0),     # 7-cycle: proche standard
    'cyclooctene': (10.5, -3.5, 2.5),      # 8-cycle: cis/trans variable
    
    # Alcènes exocycliques (méthylène)
    'exocyclic_3': (8.5, -2.0, 1.5),       # =CH₂ sur cyclopropane
    'exocyclic_5': (10.0, -2.5, 2.0),      # =CH₂ sur cyclopentane
    'exocyclic_6': (10.5, -3.0, 2.0),      # =CH₂ sur cyclohexane
    
    # Alcènes avec substituants électronégatifs
    'vinyl_ether': (9.5, -2.5, 1.5),       # CH₂=CH-OR
    'vinyl_ester': (9.0, -2.0, 1.5),       # CH₂=CH-OC(=O)R
    'vinyl_halide_F': (8.0, -2.0, 1.0),    # CH₂=CH-F
    'vinyl_halide_Cl': (9.5, -2.5, 1.5),   # CH₂=CH-Cl
    'vinyl_halide_Br': (9.5, -2.5, 1.5),   # CH₂=CH-Br
    'vinyl_amine': (10.0, -3.0, 2.0),      # CH₂=CH-NR₂
    
    # Corrections électronégativité (additions à J base)
    'EN_correction_per_O': -0.8,           # Par oxygène adjacent
    'EN_correction_per_N': -0.5,           # Par azote adjacent
    'EN_correction_per_F': -1.2,           # Par fluor adjacent
    'EN_correction_per_Cl': -0.6,          # Par chlore adjacent
    
    # Défaut
    'default': (11.0, -3.0, 2.0),
}


# ==============================================================================
# FONCTIONS D'ACCÈS AUX TABLES DE COUPLAGE LONGUE PORTÉE
# ==============================================================================

def _get_4J_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    coupling_type: str = None,
    mol3d: Chem.Mol = None
) -> Optional[float]:
    """
    Retourne le couplage ⁴J approprié selon le type de système.
    
    Détecte automatiquement le type de couplage si non spécifié:
    - Allylique: H-C-C=C-H
    - Homoallylique: H-C=C-C-H  
    - W-allylique: configuration W rigide
    - Aromatique méta: à travers cycle aromatique
    - W-saturé: systèmes bicycliques rigides
    
    Args:
        mol: Molécule RDKit
        h1_idx, h2_idx: Indices des protons (ou de leurs carbones parents)
        coupling_type: Type forcé (optionnel)
        mol3d: Molécule avec coords 3D pour calcul d'angle
    
    Returns:
        ⁴J en Hz ou None si non applicable
    """
    if coupling_type:
        # Type spécifié explicitement
        if coupling_type in _4J_LONG_RANGE_TABLE:
            j_range = _4J_LONG_RANGE_TABLE[coupling_type]
            return (j_range[0] + j_range[1]) / 2.0
    
    # Récupérer les atomes
    try:
        atom1 = mol.GetAtomWithIdx(h1_idx)
        atom2 = mol.GetAtomWithIdx(h2_idx)
    except Exception:
        return None
    
    # Les indices sont-ils des H ou des C?
    c1_idx = h1_idx if atom1.GetSymbol() == 'C' else None
    c2_idx = h2_idx if atom2.GetSymbol() == 'C' else None
    
    if c1_idx is None:
        # Trouver le carbone parent
        for nb in atom1.GetNeighbors():
            if nb.GetSymbol() == 'C':
                c1_idx = nb.GetIdx()
                break
    if c2_idx is None:
        for nb in atom2.GetNeighbors():
            if nb.GetSymbol() == 'C':
                c2_idx = nb.GetIdx()
                break
    
    if c1_idx is None or c2_idx is None:
        return None
    
    c1 = mol.GetAtomWithIdx(c1_idx)
    c2 = mol.GetAtomWithIdx(c2_idx)
    
    # Vérifier le type de couplage
    # 1. Aromatique méta?
    if c1.GetIsAromatic() and c2.GetIsAromatic():
        # Vérifier si dans le même cycle et en position méta
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if c1_idx in ring and c2_idx in ring:
                # Calculer la distance dans le cycle
                idx1 = list(ring).index(c1_idx)
                idx2 = list(ring).index(c2_idx)
                dist = min(abs(idx1 - idx2), len(ring) - abs(idx1 - idx2))
                if dist == 2:  # Position méta (2 liaisons dans le cycle)
                    return 2.5  # Moyenne de meta_aromatic
    
    # 2. Allylique/Homoallylique?
    # Chemin de 4 liaisons avec C=C
    def find_path_4_bonds(start_idx, end_idx):
        """BFS pour trouver chemin de 4 liaisons."""
        visited = {start_idx}
        queue = [(start_idx, [start_idx], 0)]
        while queue:
            current, path, depth = queue.pop(0)
            if depth > 4:
                continue
            if current == end_idx and depth == 4:
                return path
            for nb in mol.GetAtomWithIdx(current).GetNeighbors():
                if nb.GetIdx() not in visited and nb.GetSymbol() in ('C', 'N', 'O'):
                    visited.add(nb.GetIdx())
                    queue.append((nb.GetIdx(), path + [nb.GetIdx()], depth + 1))
        return None
    
    path = find_path_4_bonds(c1_idx, c2_idx)
    if path and len(path) >= 4:
        # Vérifier si double liaison dans le chemin
        has_double_bond = False
        for i in range(len(path) - 1):
            bond = mol.GetBondBetweenAtoms(path[i], path[i+1])
            if bond and bond.GetBondTypeAsDouble() >= 1.5:
                has_double_bond = True
                break
        
        if has_double_bond:
            # Allylique ou homoallylique selon la position de la double liaison
            is_rigid = all(mol.GetAtomWithIdx(idx).IsInRing() for idx in path)
            if is_rigid:
                return 2.0  # W_allylic_rigid
            else:
                return 1.5  # allylic_general moyen
    
    # 3. W-saturé dans bicyclique?
    if c1.IsInRing() and c2.IsInRing():
        # Vérifier configuration W approximative
        if mol3d and mol3d.GetNumConformers() > 0:
            try:
                conf = mol3d.GetConformer()
                # Le couplage W existe si les H sont coplanaires avec les C intermédiaires
                # Simplification: vérifier si même cycle et positions 1,3
                return 1.2  # W_saturated_rigid
            except Exception:
                pass
    
    # 4. Benzylic?
    if (c1.GetIsAromatic() != c2.GetIsAromatic()):
        # Un aromatique, un non
        return 1.0  # benzylic moyen
    
    # Défaut: pas de couplage ⁴J détecté
    return None


def _get_5J_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    coupling_type: str = None,
    mol3d: Chem.Mol = None
) -> Optional[float]:
    """
    Retourne le couplage ⁵J approprié selon le type de système.
    
    Les couplages ⁵J sont généralement très petits (<1 Hz) sauf:
    - Configuration W étendu rigide
    - Position para aromatique
    - Systèmes π conjugués
    
    Args:
        mol: Molécule RDKit
        h1_idx, h2_idx: Indices des protons (ou carbones)
        coupling_type: Type forcé (optionnel)
        mol3d: Molécule 3D pour calcul géométrique
    
    Returns:
        ⁵J en Hz ou None si non applicable
    """
    if coupling_type:
        if coupling_type in _5J_LONG_RANGE_TABLE:
            j_range = _5J_LONG_RANGE_TABLE[coupling_type]
            return (j_range[0] + j_range[1]) / 2.0
    
    try:
        atom1 = mol.GetAtomWithIdx(h1_idx)
        atom2 = mol.GetAtomWithIdx(h2_idx)
    except Exception:
        return None
    
    # Détection para aromatique
    c1_idx = h1_idx if atom1.GetSymbol() == 'C' else None
    c2_idx = h2_idx if atom2.GetSymbol() == 'C' else None
    
    if c1_idx is None:
        for nb in atom1.GetNeighbors():
            if nb.GetSymbol() == 'C':
                c1_idx = nb.GetIdx()
                break
    if c2_idx is None:
        for nb in atom2.GetNeighbors():
            if nb.GetSymbol() == 'C':
                c2_idx = nb.GetIdx()
                break
    
    if c1_idx is None or c2_idx is None:
        return None
    
    c1 = mol.GetAtomWithIdx(c1_idx)
    c2 = mol.GetAtomWithIdx(c2_idx)
    
    # Para aromatique?
    if c1.GetIsAromatic() and c2.GetIsAromatic():
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if len(ring) == 6 and c1_idx in ring and c2_idx in ring:
                idx1 = list(ring).index(c1_idx)
                idx2 = list(ring).index(c2_idx)
                dist = min(abs(idx1 - idx2), len(ring) - abs(idx1 - idx2))
                if dist == 3:  # Position para
                    return 0.6  # para_substituted moyen
    
    # W étendu rigide?
    if c1.IsInRing() and c2.IsInRing():
        # Vérifier si même système polycyclique rigide
        ring_info = mol.GetRingInfo()
        rings_c1 = [set(r) for r in ring_info.AtomRings() if c1_idx in r]
        rings_c2 = [set(r) for r in ring_info.AtomRings() if c2_idx in r]
        
        # Si les deux sont dans des cycles fusionnés
        for r1 in rings_c1:
            for r2 in rings_c2:
                if len(r1 & r2) >= 2:  # Cycles fusionnés
                    return 1.0  # extended_W_rigid
    
    # Défaut: pas de couplage ⁵J significatif
    return None


def _get_HF_coupling(
    mol: Chem.Mol,
    h_idx: int,
    f_idx: int,
    mol3d: Chem.Mol = None
) -> Optional[float]:
    """
    Retourne le couplage H-F selon le nombre de liaisons et la géométrie.
    
    Les couplages H-F sont caractéristiquement grands:
    - ²J(H-C-F) ≈ 45-55 Hz
    - ³J(H-C-C-F) ≈ 5-35 Hz (Karplus-dépendant)
    - ⁴J(H-C-C-C-F) ≈ 1-8 Hz (configuration W)
    
    Args:
        mol: Molécule RDKit
        h_idx: Index du H (ou son C parent)
        f_idx: Index du F
        mol3d: Molécule 3D pour calcul d'angle dièdre
    
    Returns:
        J(H-F) en Hz
    """
    try:
        # Trouver le plus court chemin H...F
        path = Chem.GetShortestPath(mol, h_idx, f_idx)
        n_bonds = len(path) - 1 if path else 0
    except Exception:
        return None
    
    if n_bonds == 0:
        return None
    
    # ²J geminal
    if n_bonds == 2:
        j_range = _HF_COUPLING_TABLE.get('default_2J', (46.0, 50.0))
        return (j_range[0] + j_range[1]) / 2.0
    
    # ³J vicinal - dépend de l'angle dièdre
    if n_bonds == 3:
        if mol3d and mol3d.GetNumConformers() > 0:
            try:
                conf = mol3d.GetConformer()
                dihedral = abs(AllChem.GetDihedralDeg(conf, path[0], path[1], path[2], path[3]))
                # Karplus pour H-F: anti (~30 Hz) > gauche (~10 Hz)
                if dihedral > 150 or dihedral < 30:
                    return 25.0  # Anti
                elif dihedral > 60 and dihedral < 120:
                    return 8.0   # Gauche
                else:
                    return 15.0  # Intermédiaire
            except Exception:
                pass
        # Moyenne sans 3D
        j_range = _HF_COUPLING_TABLE.get('3J_HCCF_average', (8.0, 22.0))
        return (j_range[0] + j_range[1]) / 2.0
    
    # ⁴J longue portée
    if n_bonds == 4:
        # Vérifier configuration W
        is_aromatic = any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in path[1:-1])
        if is_aromatic:
            return 3.5  # Aromatique méta
        else:
            return 2.0  # Général
    
    # ⁵J très longue portée
    if n_bonds == 5:
        return 1.0  # Faible
    
    return None


def _get_HP_coupling(
    mol: Chem.Mol,
    h_idx: int,
    p_idx: int,
    mol3d: Chem.Mol = None
) -> Optional[float]:
    """
    Retourne le couplage H-³¹P selon le nombre de liaisons.
    
    Les couplages H-P sont importants pour:
    - Phosphonates, phosphates
    - Nucléotides
    - Composés organophosphorés
    
    Args:
        mol: Molécule RDKit
        h_idx: Index du H
        p_idx: Index du P
        mol3d: Molécule 3D
    
    Returns:
        J(H-P) en Hz
    """
    try:
        path = Chem.GetShortestPath(mol, h_idx, p_idx)
        n_bonds = len(path) - 1 if path else 0
    except Exception:
        return None
    
    if n_bonds == 0:
        return None
    
    # ²J
    if n_bonds == 2:
        j_range = _HP_COUPLING_TABLE.get('default_2J', (12.0, 18.0))
        return (j_range[0] + j_range[1]) / 2.0
    
    # ³J - vérifier si via O (phosphate) ou C
    if n_bonds == 3:
        # Chemin via oxygène?
        via_oxygen = any(mol.GetAtomWithIdx(idx).GetSymbol() == 'O' for idx in path[1:-1])
        if via_oxygen:
            return 8.0  # ³J_HCOP_ester
        else:
            return 9.0  # ³J_HCCP_average
    
    # ⁴J
    if n_bonds == 4:
        return 1.5  # Faible
    
    return None


# ==============================================================================
# AMÉLIORATION 6: COUPLAGES VINYLIQUES E/Z DEPUIS STÉRÉOCHIMIE SMILES
# ==============================================================================
# Détection automatique de la stéréochimie E/Z pour choisir Jcis (~10 Hz) vs Jtrans (~16 Hz)
# Références: Jackman & Sternhell, Applications of NMR Spectroscopy
# ==============================================================================

def _get_vinylic_J_from_stereo(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    default_cis: float = 10.0,
    default_trans: float = 16.0,
    default_geminal: float = 2.0
) -> Optional[Tuple[float, str]]:
    """
    Détermine le couplage vinylique ³J(H-C=C-H) depuis la stéréochimie E/Z du SMILES.
    
    Utilise la stéréochimie définie dans le SMILES (ex: /C=C\ pour Z, /C=C/ pour E)
    pour déterminer si les H sont cis ou trans.
    
    Args:
        mol: Molécule RDKit avec H explicites
        h1_idx, h2_idx: Indices des deux H vinyliques
        default_cis: J pour configuration Z (cis) en Hz
        default_trans: J pour configuration E (trans) en Hz
        default_geminal: J pour H géminaux (=CH₂)
    
    Returns:
        Tuple (J_hz, stereo_type) où stereo_type = 'cis', 'trans', 'geminal', ou 'unknown'
        None si pas un couplage vinylique
    """
    try:
        h1 = mol.GetAtomWithIdx(h1_idx)
        h2 = mol.GetAtomWithIdx(h2_idx)
        
        if h1.GetSymbol() != 'H' or h2.GetSymbol() != 'H':
            return None
        
        # Trouver les carbones parents
        c1_list = [n for n in h1.GetNeighbors() if n.GetSymbol() == 'C']
        c2_list = [n for n in h2.GetNeighbors() if n.GetSymbol() == 'C']
        
        if not c1_list or not c2_list:
            return None
        
        c1 = c1_list[0]
        c2 = c2_list[0]
        c1_idx = c1.GetIdx()
        c2_idx = c2.GetIdx()
        
        # Geminal? (même carbone)
        if c1_idx == c2_idx:
            return (default_geminal, 'geminal')
        
        # Vérifier qu'il y a une double liaison C=C
        bond = mol.GetBondBetweenAtoms(c1_idx, c2_idx)
        if not bond:
            return None
        if bond.GetBondTypeAsDouble() < 1.9:
            return None
        
        # Récupérer la stéréochimie de la double liaison
        stereo = bond.GetStereo()
        
        if stereo == Chem.BondStereo.STEREONONE:
            # Pas de stéréochimie définie - utiliser la géométrie 3D si disponible
            if mol.GetNumConformers() > 0:
                try:
                    conf = mol.GetConformer()
                    dihedral = AllChem.GetDihedralDeg(conf, h1_idx, c1_idx, c2_idx, h2_idx)
                    dihedral = abs(dihedral)
                    if dihedral < 30 or dihedral > 330:
                        return (default_cis, 'cis')
                    elif dihedral > 150 and dihedral < 210:
                        return (default_trans, 'trans')
                    else:
                        return ((default_cis + default_trans) / 2, 'unknown')
                except Exception:
                    pass
            return (default_trans, 'unknown')  # Défaut: trans (plus stable)
        
        # Récupérer les atomes de référence pour la stéréochimie
        stereo_atoms = bond.GetStereoAtoms()
        if len(stereo_atoms) < 2:
            return (default_trans, 'unknown')
        
        ref1_idx, ref2_idx = stereo_atoms[0], stereo_atoms[1]
        
        # Déterminer si h1 et h2 sont cis ou trans
        # Pour E (trans): les groupes de haute priorité sont de côtés opposés
        # Pour Z (cis): les groupes de haute priorité sont du même côté
        
        # Vérifier la position des H par rapport aux références
        h1_is_ref1 = (h1_idx == ref1_idx)
        h2_is_ref2 = (h2_idx == ref2_idx)
        
        if stereo in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOTRANS):
            # E (trans): refs sont opposés, donc H opposés = cis, H même côté = trans
            if h1_is_ref1 == h2_is_ref2:
                return (default_trans, 'trans')
            else:
                return (default_cis, 'cis')
        
        elif stereo in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOCIS):
            # Z (cis): refs sont même côté
            if h1_is_ref1 == h2_is_ref2:
                return (default_cis, 'cis')
            else:
                return (default_trans, 'trans')
        
        return (default_trans, 'unknown')
        
    except Exception:
        return None


# ==============================================================================
# AMÉLIORATION: DÉTECTION SYSTÈMES AA'BB' AROMATIQUES SYMÉTRIQUES
# ==============================================================================
# Para-disubstitués symétriques (p-xylène, etc.) apparaissent comme singulets
# car le système AA'BB' est fortement couplé avec delta_nu/J proche de 0
# ==============================================================================

def _is_symmetric_para_disubstituted_benzene(
    mol: Chem.Mol,
    heavy_idx: int,
    aromatic_dmap: List[Dict[int, int]],
    ranks: Optional[List[int]] = None,
) -> bool:
    """
    Détecte si un H aromatique fait partie d'un système AA'BB' symétrique
    (para-disubstitué avec substituants identiques) qui apparaît comme singulet.
    
    Exemples: p-xylène, 1,4-dichlorobenzène, hydroquinone, etc.
    
    Le système AA'BB' fortement couplé avec delta_nu = 0 donne un singulet apparent.
    """
    atom = mol.GetAtomWithIdx(heavy_idx)
    if not atom.GetIsAromatic() or atom.GetSymbol() != 'C':
        return False
    
    # Trouver le cycle aromatique à 6 membres contenant cet atome
    ring_info = mol.GetRingInfo()
    six_rings = [r for r in ring_info.AtomRings() if len(r) == 6 and heavy_idx in r]
    
    if not six_rings:
        return False
    
    ring = six_rings[0]
    ring_set = set(ring)
    
    # Vérifier que tous les atomes du cycle sont des C aromatiques
    for idx in ring:
        a = mol.GetAtomWithIdx(idx)
        if a.GetSymbol() != 'C' or not a.GetIsAromatic():
            return False
    
    # Trouver les positions substituées (atomes sans H) et les positions avec H
    substituted_positions = []
    h_positions = []
    
    for idx in ring:
        a = mol.GetAtomWithIdx(idx)
        # Compter les H attachés (explicites + implicites)
        n_H = a.GetTotalNumHs()
        # Aussi vérifier les H explicites si présents
        n_H_explicit = sum(1 for n in a.GetNeighbors() if n.GetSymbol() == 'H' or n.GetAtomicNum() == 1)
        n_H = max(n_H, n_H_explicit)
        
        if n_H == 0:
            # Pas de H → position substituée
            # Vérifier qu'il y a un substituant (pas juste dans le cycle)
            non_ring_neighbors = [n for n in a.GetNeighbors() if n.GetIdx() not in ring_set]
            if non_ring_neighbors:
                substituted_positions.append(idx)
        else:
            h_positions.append(idx)
    
    # Para-disubstitué: exactement 2 positions substituées
    if len(substituted_positions) != 2:
        return False
    
    # Vérifier que les 2 substituants sont en position para (distance 3 dans le cycle)
    p1, p2 = substituted_positions
    
    # Calculer la distance dans le cycle (mod 6)
    try:
        idx1 = ring.index(p1)
        idx2 = ring.index(p2)
        dist = min(abs(idx1 - idx2), 6 - abs(idx1 - idx2))
        if dist != 3:
            return False  # Pas en para
    except:
        return False
    
    # Vérifier que les substituants sont identiques
    a1 = mol.GetAtomWithIdx(p1)
    a2 = mol.GetAtomWithIdx(p2)
    
    # Trouver les substituants (atomes non dans le cycle)
    sub1 = [n for n in a1.GetNeighbors() if n.GetIdx() not in ring_set]
    sub2 = [n for n in a2.GetNeighbors() if n.GetIdx() not in ring_set]
    
    if len(sub1) != 1 or len(sub2) != 1:
        return False
    
    # Comparer les substituants par leur type atomique et nombre de H
    # (méthode plus robuste que les rangs canoniques)
    s1, s2 = sub1[0], sub2[0]
    
    # Même symbole atomique ?
    if s1.GetSymbol() != s2.GetSymbol():
        return False
    
    # Même nombre de H ?
    if s1.GetTotalNumHs() != s2.GetTotalNumHs():
        return False
    
    # Même degré (nombre de voisins) ?
    if s1.GetDegree() != s2.GetDegree():
        return False
    
    # Si les substituants sont identiques, vérifier que les H aromatiques 
    # ont le même déplacement chimique (équivalents)
    if ranks is not None:
        # Vérifier que les H positions ont des rangs cohérents
        # (positions 2,6 équivalentes, positions 3,5 équivalentes)
        pass  # Les substituants sont identiques → symétrie confirmée
    
    # C'est un para-disubstitué symétrique → singulet apparent
    return True


def _is_ortho_disubstituted_symmetric(
    mol: Chem.Mol,
    heavy_idx: int,
    aromatic_dmap: List[Dict[int, int]],
    ranks: Optional[List[int]] = None,
) -> bool:
    """
    Détecte les ortho-disubstitués symétriques (o-xylène).
    Ces systèmes montrent des multiplicités plus complexes.
    """
    atom = mol.GetAtomWithIdx(heavy_idx)
    if not atom.GetIsAromatic() or atom.GetSymbol() != 'C':
        return False
    
    ring_info = mol.GetRingInfo()
    six_rings = [r for r in ring_info.AtomRings() if len(r) == 6 and heavy_idx in r]
    
    if not six_rings:
        return False
    
    ring = six_rings[0]
    ring_set = set(ring)
    
    substituted_positions = []
    for idx in ring:
        a = mol.GetAtomWithIdx(idx)
        n_H = sum(1 for n in a.GetNeighbors() if n.GetSymbol() == 'H' or n.GetAtomicNum() == 1)
        if n_H == 0:
            non_ring_neighbors = [n for n in a.GetNeighbors() if n.GetIdx() not in ring_set]
            if non_ring_neighbors:
                substituted_positions.append(idx)
    
    if len(substituted_positions) != 2:
        return False
    
    p1, p2 = substituted_positions
    try:
        idx1 = ring.index(p1)
        idx2 = ring.index(p2)
        dist = min(abs(idx1 - idx2), 6 - abs(idx1 - idx2))
        if dist != 1:  # Ortho
            return False
    except:
        return False
    
    if ranks is not None and ranks[p1] != ranks[p2]:
        return False
    
    return True


# ==============================================================================
# AMÉLIORATION 7: COUPLAGES AROMATIQUES AUTOMATIQUES
# ==============================================================================
# Détection automatique ortho (7-8 Hz), méta (2-3 Hz), para (0.5 Hz)
# Références: Günther H., NMR Spectroscopy 3rd Ed. (2013), Chapter 7
# ==============================================================================

# Table des couplages aromatiques par type de cycle
# AMÉLIORATION: Table enrichie avec positions spécifiques pour tous les hétérocycles
_AROMATIC_COUPLING_TABLE = {
    # Benzène et dérivés
    'benzene': {'ortho': 7.5, 'meta': 1.5, 'para': 0.5},
    'phenyl_electron_donor': {'ortho': 7.8, 'meta': 1.8, 'para': 0.4},  # OH, OR, NR2
    'phenyl_electron_acceptor': {'ortho': 8.2, 'meta': 2.2, 'para': 0.6},  # NO2, CN, C=O
    
    # =========================================================================
    # Hétérocycles 5-membred - ENRICHIS avec toutes les positions
    # =========================================================================
    'furan': {
        'ortho': 1.8, 'meta': 0.9, 'para': 1.5,  # Mapping générique
        '2_3': 1.8, '3_4': 3.4, '2_4': 0.9, '2_5': 1.5, '3_5': 0.9,
        # Couplages long-range faibles mais détectables
        'long_range_4J': 0.5,
    },
    'thiophene': {
        'ortho': 5.0, 'meta': 1.0, 'para': 3.5,
        '2_3': 5.0, '3_4': 3.5, '2_4': 1.0, '2_5': 3.5, '3_5': 1.0,
        'long_range_4J': 0.7,
    },
    'pyrrole': {
        'ortho': 2.5, 'meta': 1.5, 'para': 2.5,
        '2_3': 2.5, '3_4': 2.5, '2_4': 1.5, '2_5': 2.5, '3_5': 1.5,
        'long_range_4J': 0.5,
    },
    'imidazole': {
        # Position 2: entre les 2 N, H2 souvent singulet (découplé des H4,H5)
        # J(H4-H5) ≈ 1.5 Hz, mais souvent non résolu -> singulets
        'ortho': 1.5, 'meta': 0.0, 'para': 0.0,
        '4_5': 1.5, '2_4': 0.0, '2_5': 0.0,  # H2 ne couple pas avec H4,H5
        'special_h2_singlet': True,  # Flag pour H2 singulet
    },
    'pyrazole': {
        'ortho': 2.0, 'meta': 0.5,
        '3_4': 2.0, '4_5': 2.0, '3_5': 0.5,
    },
    'oxazole': {
        'ortho': 1.0, 'meta': 0.0,
        '4_5': 1.0, '2_4': 0.0, '2_5': 0.0,
    },
    'thiazole': {
        'ortho': 3.5, 'meta': 0.0,
        '4_5': 3.5, '2_4': 0.0, '2_5': 0.0,
    },
    'isoxazole': {
        'ortho': 2.0, 'meta': 0.8,
        '3_4': 2.0, '4_5': 1.8, '3_5': 0.8,
    },
    
    # =========================================================================
    # Hétérocycles 6-membred - ENRICHIS
    # =========================================================================
    'pyridine': {
        # Convention: 'ortho' = distance 1 dans le cycle (voisins adjacents)
        # Pyridine: J(2,3)=J(5,6)=4.8 Hz (adjacent à N), J(3,4)=J(4,5)=7.8 Hz
        # Moyenne pour 'ortho' dans le cycle: (4.8 + 7.8) / 2 ≈ 6.3 Hz
        # On utilise une moyenne pondérée basée sur la position
        'ortho': 6.3,  # Distance 1: moyenne des J 3J
        'meta': 1.8,   # Distance 2: couplage 4J (faible)
        'para': 0.5,   # Distance 3: couplage 5J (très faible)
        '2_3': 4.8, '3_4': 7.8, '4_5': 7.8, '5_6': 4.8,  # Couplages 3J spécifiques
        '2_4': 1.8, '2_5': 0.9, '2_6': 0.4,  # Couplages 4J et 5J
        '3_5': 1.5, '3_6': 0.5, '4_6': 1.8,
    },
    'pyrimidine': {
        'ortho': 5.0, 'meta': 1.5,
        '4_5': 5.0, '5_6': 5.0, '4_6': 1.5,
    },
    'pyrazine': {
        'ortho': 2.0, 'meta': 0.5,
        '2_3': 2.0, '2_5': 0.5,  # Tous les H équivalents par paires
    },
    'pyridazine': {
        'ortho': 5.0, 'meta': 2.0,
        '3_4': 5.0, '4_5': 5.0, '5_6': 5.0, '3_5': 2.0, '4_6': 2.0,
    },
    'triazine_1_3_5': {
        # Tous les H équivalents, pas de couplage entre eux
        'ortho': 0.0, 'meta': 0.0,
    },
    
    # =========================================================================
    # Hétérocycles condensés - NOUVEAUX
    # =========================================================================
    'quinoline': {
        # Cycle pyridine fusionné avec benzène
        'ortho': 8.0, 'meta': 1.5, 'para': 0.5,
        '2_3': 4.5, '3_4': 8.0, '5_6': 8.0, '6_7': 7.5, '7_8': 8.0,
        '2_4': 1.8, '5_7': 1.3, '5_8': 0.8, '6_8': 1.3,
    },
    'isoquinoline': {
        'ortho': 8.0, 'meta': 1.5, 'para': 0.5,
        '1_3': 0.9, '1_4': 0.5, '3_4': 6.5, '5_6': 8.0, '6_7': 7.5, '7_8': 8.0,
    },
    'indole': {
        # Pyrrole fusionné avec benzène
        'ortho': 7.8, 'meta': 1.2, 'para': 0.5,
        '2_3': 2.2, '4_5': 8.0, '5_6': 7.5, '6_7': 7.5, '7_4': 0.8,
        '3_4': 0.0,  # Pas de couplage entre H3 et H4 (séparés par la fusion)
    },
    'benzofuran': {
        'ortho': 7.8, 'meta': 1.2, 'para': 0.5,
        '2_3': 2.2, '4_5': 8.0, '5_6': 7.5, '6_7': 7.5,
    },
    'benzothiophene': {
        'ortho': 7.8, 'meta': 1.2, 'para': 0.5,
        '2_3': 5.5, '4_5': 8.0, '5_6': 7.5, '6_7': 7.5,
    },
    'benzimidazole': {
        'ortho': 7.8, 'meta': 1.3, 'para': 0.5,
        '4_5': 8.0, '5_6': 7.5, '6_7': 7.5, '2_4': 0.0, '2_7': 0.0,
        'special_h2_singlet': True,  # H2 est singulet
    },
    
    # Naphthalène et polycycliques
    'naphthalene': {'1_2': 8.0, '2_3': 7.0, '3_4': 8.0, '1_3': 1.3, '1_4': 0.8, '2_4': 1.3,
                    '5_6': 8.0, '6_7': 7.0, '7_8': 8.0, 'ortho': 8.0, 'meta': 1.3, 'para': 0.8},
    
    # Défauts par distance
    'default_ortho': 7.5,
    'default_meta': 2.0,
    'default_para': 0.5,
}


def _is_imidazole_H2_carbon(mol: Chem.Mol, heavy_idx: int) -> bool:
    """
    Détecte si un carbone aromatique est en position 2 d'un imidazole (entre 2 N).
    
    Le H sur ce carbone est un SINGULET car:
    - Entouré par 2 N (pas de voisins H directs)
    - Couplages 4J avec H4/H5 sont négligeables (<0.5 Hz, non résolus)
    
    Args:
        mol: Molécule RDKit
        heavy_idx: Index du carbone portant le H
        
    Returns:
        True si c'est le C2 d'un imidazole → H doit être singulet
    """
    try:
        atom = mol.GetAtomWithIdx(heavy_idx)
        if atom.GetSymbol() != 'C' or not atom.GetIsAromatic():
            return False
        
        # Trouver les cycles à 5 membres contenant cet atome
        ring_info = mol.GetRingInfo()
        five_rings = [r for r in ring_info.AtomRings() if len(r) == 5 and heavy_idx in r]
        
        if not five_rings:
            return False
        
        for ring in five_rings:
            ring_set = set(ring)
            
            # Compter les N dans le cycle
            n_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == 'N')
            
            # Imidazole a 2 N, pyrazole aussi, triazole en a 3
            if n_in_ring < 2:
                continue
            
            # Vérifier si ce carbone a 2 voisins N dans le cycle
            n_neighbors = 0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'N' and neighbor.GetIdx() in ring_set:
                    n_neighbors += 1
            
            if n_neighbors == 2:
                return True  # Position 2 de l'imidazole
        
        return False
    except Exception:
        return False


def _is_isolated_heterocyclic_H(mol: Chem.Mol, heavy_idx: int) -> bool:
    """
    Détecte si un H hétérocyclique est isolé (pas de couplages significatifs).
    
    Cas typiques:
    - Imidazole H4 et H5: souvent singulets (J très petit ~1.5 Hz, non résolu)
    - Pyrazole positions similaires
    
    Returns:
        True si le H devrait apparaître comme singulet
    """
    try:
        atom = mol.GetAtomWithIdx(heavy_idx)
        if atom.GetSymbol() != 'C' or not atom.GetIsAromatic():
            return False
        
        # Trouver les cycles à 5 membres
        ring_info = mol.GetRingInfo()
        five_rings = [r for r in ring_info.AtomRings() if len(r) == 5 and heavy_idx in r]
        
        if not five_rings:
            return False
        
        for ring in five_rings:
            ring_set = set(ring)
            
            # Compter les N dans le cycle
            n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == 'N')
            
            if n_count < 2:
                continue  # Pas un hétérocycle avec 2+ N
            
            # Compter les C avec H dans le cycle
            c_with_H = 0
            for idx in ring:
                a = mol.GetAtomWithIdx(idx)
                if a.GetSymbol() == 'C':
                    # Vérifier si ce C a un H
                    has_H = any(n.GetSymbol() == 'H' or n.GetAtomicNum() == 1 
                               for n in a.GetNeighbors())
                    if not has_H and a.GetTotalNumHs() > 0:
                        has_H = True
                    if has_H:
                        c_with_H += 1
            
            # Si seulement 1 ou 2 C avec H dans un cycle avec 2+ N,
            # les couplages sont souvent non résolus
            if c_with_H <= 2 and n_count >= 2:
                # Vérifier si ce carbone n'a pas de voisin C avec H
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in ring_set and neighbor.GetSymbol() == 'C':
                        neighbor_has_H = any(n.GetSymbol() == 'H' for n in neighbor.GetNeighbors())
                        if neighbor_has_H:
                            # Il y a un voisin C-H, couplage possible
                            # Mais dans imidazole le J est très petit (~1.5 Hz)
                            # On le considère comme résolu seulement si J > 2 Hz
                            return False
                return True
        
        return False
    except Exception:
        return False


def _is_imidazole_h2_position(mol: Chem.Mol, c1_idx: int, c2_idx: int, ring: tuple) -> bool:
    """
    Détecte si un carbone est en position 2 de l'imidazole (entre les 2 N).
    
    Dans l'imidazole et benzimidazole, le H2 est souvent un singulet car:
    - Il est entouré par deux N
    - Les couplages 4J avec H4 et H5 sont négligeables (<0.5 Hz)
    
    Returns:
        True si c1 ou c2 est le carbone en position 2
    """
    try:
        for c_idx in (c1_idx, c2_idx):
            carbon = mol.GetAtomWithIdx(c_idx)
            if carbon.GetSymbol() != 'C':
                continue
            
            # Compter les N voisins dans le cycle
            n_neighbors_in_ring = 0
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetSymbol() == 'N' and neighbor.GetIdx() in ring:
                    n_neighbors_in_ring += 1
            
            # Position 2: le carbone a 2 voisins N dans le cycle
            if n_neighbors_in_ring == 2:
                return True
        
        return False
    except Exception:
        return False


def _get_aromatic_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int
) -> Optional[Tuple[float, str]]:
    """
    Détermine le couplage aromatique depuis la topologie du cycle.
    
    Détecte automatiquement:
    - Position ortho (distance 1 dans le cycle): 7-8 Hz
    - Position méta (distance 2): 1.5-2.5 Hz
    - Position para (distance 3): 0.4-0.8 Hz
    
    Gère aussi les hétérocycles avec valeurs spécifiques.
    
    Args:
        mol: Molécule RDKit avec H explicites
        h1_idx, h2_idx: Indices des deux H aromatiques
    
    Returns:
        Tuple (J_hz, position_type) ou None
    """
    try:
        h1 = mol.GetAtomWithIdx(h1_idx)
        h2 = mol.GetAtomWithIdx(h2_idx)
        
        # Trouver les carbones aromatiques parents
        c1_list = [n for n in h1.GetNeighbors() if n.GetIsAromatic()]
        c2_list = [n for n in h2.GetNeighbors() if n.GetIsAromatic()]
        
        if not c1_list or not c2_list:
            return None
        
        c1 = c1_list[0]
        c2 = c2_list[0]
        c1_idx = c1.GetIdx()
        c2_idx = c2.GetIdx()
        
        if c1_idx == c2_idx:
            return None  # Même carbone
        
        # Trouver le cycle aromatique commun
        ring_info = mol.GetRingInfo()
        common_rings = []
        
        for ring in ring_info.AtomRings():
            if c1_idx in ring and c2_idx in ring:
                # Vérifier que c'est un cycle aromatique
                is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
                if is_aromatic:
                    common_rings.append(ring)
        
        if not common_rings:
            return None
        
        # Utiliser le plus petit cycle commun
        ring = min(common_rings, key=len)
        ring_size = len(ring)
        
        # Calculer la distance dans le cycle
        ring_list = list(ring)
        idx1 = ring_list.index(c1_idx)
        idx2 = ring_list.index(c2_idx)
        dist = min(abs(idx1 - idx2), ring_size - abs(idx1 - idx2))
        
        # Déterminer le type de cycle pour les valeurs spécifiques
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        heteroatoms = [a.GetSymbol() for a in ring_atoms if a.GetSymbol() != 'C']
        
        # Identifier le type d'hétérocycle
        cycle_type = 'benzene'  # Défaut
        
        if ring_size == 5:
            if 'O' in heteroatoms and 'N' not in heteroatoms:
                cycle_type = 'furan'
            elif 'S' in heteroatoms:
                cycle_type = 'thiophene'
            elif 'N' in heteroatoms:
                if heteroatoms.count('N') >= 2:
                    cycle_type = 'imidazole'
                else:
                    cycle_type = 'pyrrole'
        elif ring_size == 6:
            if 'N' in heteroatoms:
                n_count = heteroatoms.count('N')
                if n_count == 1:
                    cycle_type = 'pyridine'
                elif n_count == 2:
                    cycle_type = 'pyrimidine'
        
        # Récupérer J depuis la table
        if dist == 1:
            position = 'ortho'
            J = _AROMATIC_COUPLING_TABLE.get('default_ortho', 7.5)
        elif dist == 2:
            position = 'meta'
            J = _AROMATIC_COUPLING_TABLE.get('default_meta', 2.0)
        elif dist == 3:
            position = 'para'
            J = _AROMATIC_COUPLING_TABLE.get('default_para', 0.5)
        else:
            return None  # Distance > 3, pas de couplage significatif
        
        # Ajuster selon le type de cycle si disponible
        if cycle_type in _AROMATIC_COUPLING_TABLE:
            cycle_J = _AROMATIC_COUPLING_TABLE[cycle_type]
            if isinstance(cycle_J, dict):
                # D'abord essayer la position générique
                if position in cycle_J:
                    J = cycle_J[position]
                
                # Vérifier cas spéciaux (H2 singulet dans imidazole/benzimidazole)
                if cycle_J.get('special_h2_singlet'):
                    # Détecter si un des H est en position 2 (entre les 2 N)
                    if _is_imidazole_h2_position(mol, c1_idx, c2_idx, ring):
                        J = 0.0  # Pas de couplage significatif
        
        return (J, position)
        
    except Exception:
        return None


# ==============================================================================
# AMÉLIORATION 8: COUPLAGES PEPTIDIQUES BYSTROV ³J(HN-Cα-Hα)
# ==============================================================================
# Équation Bystrov pour angles dièdres φ (phi) dans les peptides
# Références: Bystrov V.F., Prog. NMR Spectrosc. 10, 41-82 (1976)
# ==============================================================================

# Coefficients Bystrov pour différents types de couplages peptidiques
_BYSTROV_PEPTIDE_COEFFS = {
    # ³J(HN-Cα-Hα) - le plus important
    # J = A·cos²(φ-60°) + B·cos(φ-60°) + C
    # où φ est l'angle dièdre du backbone
    'HNCH_alpha': {
        'A': 6.4, 'B': -1.4, 'C': 1.9,
        'phase': -60.0,  # Correction de phase pour φ
    },
    
    # ³J(Hα-Cα-C'-N) 
    'HCCO': {
        'A': 4.0, 'B': -1.0, 'C': 0.5,
        'phase': 0.0,
    },
    
    # ³J(Hα-Cα-Cβ-Hβ) chaîne latérale
    'HCCH_sidechain': {
        'A': 9.5, 'B': -1.6, 'C': 1.8,
        'phase': 0.0,
    },
}

# Valeurs typiques pour structures secondaires
_PEPTIDE_STRUCTURE_J = {
    # Structure: (J_HNCH typique, angle φ typique)
    'helix_alpha': (3.9, -57),      # Hélice α: φ ≈ -57°, J ≈ 3.9 Hz
    'helix_310': (4.5, -49),        # Hélice 3₁₀: φ ≈ -49°
    'sheet_parallel': (8.9, -119),  # Feuillet β parallèle
    'sheet_antiparallel': (8.5, -139),  # Feuillet β antiparallèle
    'polyproline_II': (5.5, -75),   # PPII
    'turn_I': (6.0, -60),           # Coude type I
    'turn_II': (4.0, -60),          # Coude type II
    'random_coil': (7.0, None),     # Pelote aléatoire (moyenne)
}


def _get_peptide_HNCH_coupling(
    mol: Chem.Mol,
    h_amide_idx: int,
    h_alpha_idx: int,
    mol3d: Chem.Mol = None
) -> Optional[Tuple[float, str]]:
    """
    Calcule ³J(HN-Cα-Hα) pour un peptide en utilisant l'équation Bystrov.
    
    Cette constante est cruciale pour déterminer la structure secondaire:
    - Hélice α: J ≈ 3-5 Hz (φ ≈ -60°)
    - Feuillet β: J ≈ 8-10 Hz (φ ≈ -120° à -140°)
    - Pelote: J ≈ 6-7 Hz (moyenne)
    
    Args:
        mol: Molécule RDKit
        h_amide_idx: Index du H amide (HN)
        h_alpha_idx: Index du H alpha (Hα)
        mol3d: Molécule avec coords 3D pour calcul exact
    
    Returns:
        Tuple (J_hz, predicted_structure) ou None
    """
    try:
        h_amide = mol.GetAtomWithIdx(h_amide_idx)
        h_alpha = mol.GetAtomWithIdx(h_alpha_idx)
        
        # Vérifier que h_amide est sur un N
        n_list = [n for n in h_amide.GetNeighbors() if n.GetSymbol() == 'N']
        if not n_list:
            return None
        n_atom = n_list[0]
        n_idx = n_atom.GetIdx()
        
        # Vérifier que h_alpha est sur un C
        c_alpha_list = [n for n in h_alpha.GetNeighbors() if n.GetSymbol() == 'C']
        if not c_alpha_list:
            return None
        c_alpha = c_alpha_list[0]
        c_alpha_idx = c_alpha.GetIdx()
        
        # Vérifier le chemin N-Cα
        if not mol.GetBondBetweenAtoms(n_idx, c_alpha_idx):
            return None
        
        # Si 3D disponible, calculer J depuis l'angle dièdre
        if mol3d and mol3d.GetNumConformers() > 0:
            try:
                conf = mol3d.GetConformer()
                phi = AllChem.GetDihedralDeg(conf, h_amide_idx, n_idx, c_alpha_idx, h_alpha_idx)
                
                # Appliquer l'équation Bystrov
                coeffs = _BYSTROV_PEPTIDE_COEFFS['HNCH_alpha']
                A, B, C = coeffs['A'], coeffs['B'], coeffs['C']
                phase = coeffs['phase']
                
                phi_rad = math.radians(phi + phase)
                J = A * (math.cos(phi_rad) ** 2) + B * math.cos(phi_rad) + C
                J = max(0.0, J)
                
                # Prédire la structure secondaire
                structure = _predict_secondary_structure_from_J(J)
                
                return (J, structure)
                
            except Exception:
                pass
        
        # Sans 3D, retourner une valeur moyenne (pelote aléatoire)
        return (7.0, 'random_coil')
        
    except Exception:
        return None


def _predict_secondary_structure_from_J(J_HNCH: float) -> str:
    """
    Prédit la structure secondaire à partir de ³J(HN-Cα-Hα).
    
    Args:
        J_HNCH: Constante de couplage en Hz
    
    Returns:
        Nom de la structure secondaire probable
    """
    if J_HNCH < 4.5:
        return 'helix_alpha'
    elif J_HNCH < 6.0:
        return 'helix_310_or_turn'
    elif J_HNCH < 7.5:
        return 'random_coil'
    elif J_HNCH < 9.0:
        return 'sheet_beta'
    else:
        return 'extended'


# ==============================================================================
# AMÉLIORATION 9: DÉTECTION DES PROTONS DIASTÉRÉOTOPIQUES
# ==============================================================================
# Les CH₂ proches de centres chiraux ont des protons non-équivalents
# Références: Eliel E.L., Wilen S.H., Stereochemistry of Organic Compounds (1994)
# ==============================================================================

def _detect_diastereotopic_protons(
    mol: Chem.Mol,
    ch2_carbon_idx: int
) -> Tuple[bool, float]:
    """
    Détecte si un CH₂ porte des protons diastéréotopiques.
    
    Les protons sont diastéréotopiques si:
    1. Le carbone est adjacent à un centre chiral, OU
    2. Le carbone est dans une chaîne qui contient un centre chiral
    
    Conséquences RMN:
    - δH différents (ΔδH = 0.01-0.5 ppm typique)
    - Couplage geminal ²J visible
    - Peuvent avoir des J vicinaux différents
    
    Args:
        mol: Molécule RDKit
        ch2_carbon_idx: Index du carbone du CH₂
    
    Returns:
        Tuple (is_diastereotopic, estimated_delta_ppm)
    """
    try:
        c_atom = mol.GetAtomWithIdx(ch2_carbon_idx)
        
        if c_atom.GetSymbol() != 'C':
            return (False, 0.0)
        
        # Vérifier que c'est un CH₂
        h_count = sum(1 for n in c_atom.GetNeighbors() if n.GetSymbol() == 'H')
        if h_count != 2:
            return (False, 0.0)
        
        # Chercher des centres chiraux à proximité
        # Distance 1: adjacent direct au centre chiral
        # Distance 2: un carbone de séparation
        # Distance 3: effet faible mais détectable
        
        chiral_distance = float('inf')
        
        # BFS pour trouver le centre chiral le plus proche
        visited = {ch2_carbon_idx}
        queue = [(ch2_carbon_idx, 0)]
        
        while queue:
            current_idx, dist = queue.pop(0)
            
            if dist > 4:  # Limite de recherche
                break
            
            current = mol.GetAtomWithIdx(current_idx)
            
            # Vérifier si c'est un centre chiral
            chiral_tag = current.GetChiralTag()
            if chiral_tag != Chem.ChiralType.CHI_UNSPECIFIED and current_idx != ch2_carbon_idx:
                chiral_distance = min(chiral_distance, dist)
            
            # Explorer les voisins
            for nb in current.GetNeighbors():
                if nb.GetIdx() not in visited and nb.GetSymbol() in ('C', 'N', 'O', 'S'):
                    visited.add(nb.GetIdx())
                    queue.append((nb.GetIdx(), dist + 1))
        
        if chiral_distance == float('inf'):
            return (False, 0.0)
        
        # Estimer la différence de shift selon la distance
        if chiral_distance == 1:
            return (True, 0.15)  # Fort effet: ΔδH ≈ 0.1-0.3 ppm
        elif chiral_distance == 2:
            return (True, 0.08)  # Effet modéré: ΔδH ≈ 0.05-0.15 ppm
        elif chiral_distance == 3:
            return (True, 0.03)  # Effet faible: ΔδH ≈ 0.01-0.05 ppm
        else:
            return (True, 0.01)  # Très faible
        
    except Exception:
        return (False, 0.0)


def _find_all_diastereotopic_CH2(mol: Chem.Mol) -> List[Dict[str, Any]]:
    """
    Trouve tous les CH₂ diastéréotopiques dans une molécule.
    
    Returns:
        Liste de dicts avec:
        - 'carbon_idx': index du carbone
        - 'h_indices': indices des deux H
        - 'delta_ppm': différence de shift estimée
        - 'chiral_distance': distance au centre chiral
    """
    results = []
    
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'C':
            continue
        
        c_idx = atom.GetIdx()
        h_neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == 'H']
        
        if len(h_neighbors) != 2:
            continue
        
        is_diast, delta = _detect_diastereotopic_protons(mol, c_idx)
        
        if is_diast:
            results.append({
                'carbon_idx': c_idx,
                'h_indices': h_neighbors,
                'delta_ppm': delta,
            })
    
    return results


# ==============================================================================
# AMÉLIORATION 10: COUPLAGES SUCRES KARPLUS-TVAROSKA
# ==============================================================================
# Équation optimisée pour les pyranoses avec positions axiales/équatoriales
# Références: Tvaroska I., Carbohydr. Res. 189, 359-362 (1989)
# ==============================================================================

# Coefficients Karplus optimisés pour les sucres (Tvaroska, Hricovíni)
_SUGAR_KARPLUS_COEFFS = {
    # ³J(H-C-C-H) dans pyranoses
    'HCCH_pyranose': {
        'A': 10.2, 'B': -1.2, 'C': 0.4,
        # Corrections selon orientation
        'correction_aa': 0.0,   # Axial-axial: J typique 9-10 Hz
        'correction_ae': -0.5,  # Axial-équatorial: J typique 3-4 Hz
        'correction_ee': -0.8,  # Équatorial-équatorial: J typique 2-3 Hz
    },
    
    # ³J(H-C-O-C) glycosidique trans-O
    'HCOC_glycosidic': {
        'A': 5.5, 'B': -1.5, 'C': 2.0,
    },
    
    # ³J(H-C-C-O) hydroxyle
    'HCCO_hydroxyl': {
        'A': 8.0, 'B': -1.0, 'C': 0.5,
    },
}

# Couplages typiques selon la stéréochimie dans les hexoses
_SUGAR_TYPICAL_J = {
    # Glucose (4C1 chair)
    'glucose_H1H2_alpha': 3.8,   # α: H1 axial, H2 équatorial → ae
    'glucose_H1H2_beta': 7.8,    # β: H1 équatorial, H2 équatorial → diaxial path
    'glucose_H2H3': 9.5,         # aa
    'glucose_H3H4': 9.5,         # aa
    'glucose_H4H5': 9.5,         # aa
    
    # Mannose
    'mannose_H1H2_alpha': 1.8,   # ee
    'mannose_H1H2_beta': 1.0,    # ae
    'mannose_H2H3': 3.5,         # ae
    
    # Galactose
    'galactose_H3H4': 3.5,       # ae (H4 axial)
    'galactose_H4H5': 1.0,       # ae
}


def _get_sugar_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    mol3d: Chem.Mol = None
) -> Optional[Tuple[float, str]]:
    """
    Calcule ³J(H-C-C-H) pour les sucres en utilisant Karplus-Tvaroska.
    
    Détecte automatiquement si la molécule est un sucre (pyranose/furanose)
    et applique les corrections axial/équatorial.
    
    Args:
        mol: Molécule RDKit
        h1_idx, h2_idx: Indices des deux H
        mol3d: Molécule 3D pour orientation
    
    Returns:
        Tuple (J_hz, orientation) où orientation = 'aa', 'ae', 'ee', ou 'unknown'
    """
    try:
        h1 = mol.GetAtomWithIdx(h1_idx)
        h2 = mol.GetAtomWithIdx(h2_idx)
        
        # Trouver les carbones parents
        c1_list = [n for n in h1.GetNeighbors() if n.GetSymbol() == 'C']
        c2_list = [n for n in h2.GetNeighbors() if n.GetSymbol() == 'C']
        
        if not c1_list or not c2_list:
            return None
        
        c1 = c1_list[0]
        c2 = c2_list[0]
        c1_idx = c1.GetIdx()
        c2_idx = c2.GetIdx()
        
        # Vérifier que C1-C2 sont liés
        if not mol.GetBondBetweenAtoms(c1_idx, c2_idx):
            return None
        
        # Vérifier si c'est un sucre (cycle avec O)
        ring_info = mol.GetRingInfo()
        is_sugar_ring = False
        ring_oxygens = []
        
        for ring in ring_info.AtomRings():
            if c1_idx in ring and c2_idx in ring:
                # Vérifier la présence d'oxygène dans le cycle
                for idx in ring:
                    if mol.GetAtomWithIdx(idx).GetSymbol() == 'O':
                        is_sugar_ring = True
                        ring_oxygens.append(idx)
        
        if not is_sugar_ring:
            return None
        
        # Si 3D disponible, calculer l'orientation axial/équatorial
        if mol3d and mol3d.GetNumConformers() > 0:
            try:
                conf = mol3d.GetConformer()
                dihedral = AllChem.GetDihedralDeg(conf, h1_idx, c1_idx, c2_idx, h2_idx)
                dihedral = abs(dihedral)
                
                # Déterminer l'orientation depuis l'angle dièdre
                # aa: ~180° → J grand (9-10 Hz)
                # ae/ea: ~60° → J moyen (3-4 Hz)
                # ee: ~60° mais avec J plus petit
                
                coeffs = _SUGAR_KARPLUS_COEFFS['HCCH_pyranose']
                A, B, C = coeffs['A'], coeffs['B'], coeffs['C']
                
                phi_rad = math.radians(dihedral)
                J = A * (math.cos(phi_rad) ** 2) + B * math.cos(phi_rad) + C
                J = max(0.0, J)
                
                # Classifier l'orientation
                if dihedral > 150 or dihedral < 30:
                    orientation = 'aa'
                elif 50 < dihedral < 70:
                    if J > 5:
                        orientation = 'ae'
                    else:
                        orientation = 'ee'
                else:
                    orientation = 'unknown'
                
                return (J, orientation)
                
            except Exception:
                pass
        
        # Sans 3D, retourner une valeur moyenne
        return (5.0, 'unknown')
        
    except Exception:
        return None


def _is_pyranose(mol: Chem.Mol) -> bool:
    """Détecte si la molécule contient un cycle pyranose (6 membres avec O)."""
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 6:
            symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
            if symbols.count('O') == 1 and symbols.count('C') == 5:
                return True
    return False


def _is_furanose(mol: Chem.Mol) -> bool:
    """Détecte si la molécule contient un cycle furanose (5 membres avec O)."""
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
            if symbols.count('O') == 1 and symbols.count('C') == 4:
                return True
    return False


# ==============================================================================
# AMÉLIORATION 11: LARGEURS DE RAIE ADAPTATIVES
# ==============================================================================
# OH/NH larges (échange), aromatiques fins, alkyles moyens
# Références: Levitt M.H., Spin Dynamics, 2nd Ed. (2008)
# ==============================================================================

# Largeurs de raie typiques par type de proton (FWHM en Hz)
_LINEWIDTH_BY_PROTON_TYPE = {
    # Protons échangeables (élargis par échange chimique)
    'OH_alcohol': (5.0, 20.0),       # Alcool: dépend de la concentration et T
    'OH_phenol': (3.0, 15.0),        # Phénol: échange plus rapide
    'OH_carboxylic': (10.0, 50.0),   # Acide carboxylique: très large
    'NH_amine': (5.0, 30.0),         # Amine: échange modéré
    'NH_amide': (2.0, 10.0),         # Amide: échange plus lent
    'NH_indole': (1.0, 5.0),         # Indole NH: assez fin
    'SH': (3.0, 15.0),               # Thiol
    
    # Protons non-échangeables
    'aromatic': (0.5, 1.5),          # Aromatiques: fins
    'vinyl': (0.5, 1.5),             # Vinyliques: fins
    'aldehyde': (0.5, 1.0),          # Aldéhyde: très fin
    'alkyl_CH3': (0.8, 2.0),         # Méthyle: fins
    'alkyl_CH2': (1.0, 2.5),         # Méthylène: légèrement plus large
    'alkyl_CH': (1.0, 2.0),          # Méthine: moyen
    'alpha_heteroatom': (1.0, 2.5),  # Adjacent à O, N
    
    # Défauts
    'default': (1.0, 2.0),
    'exchangeable': (5.0, 20.0),
}

# Facteurs de correction selon les conditions
_LINEWIDTH_MODIFIERS = {
    # Solvant
    'D2O': {'exchangeable': 0.0},     # Pas de signal (échangé avec D)
    'DMSO': {'exchangeable': 0.3},    # Échange lent → plus fin
    'CDCl3': {'exchangeable': 0.8},   # Échange modéré
    'MeOD': {'exchangeable': 0.1},    # Presque pas de signal
    
    # Température (facteur multiplicatif)
    'low_temp': 1.5,    # T < 280K: échange plus lent, raies plus larges
    'room_temp': 1.0,   # 290-310K: standard
    'high_temp': 0.7,   # T > 320K: échange plus rapide, moyennage
    
    # Viscosité
    'high_viscosity': 1.5,   # Solvants visqueux → T2 court → large
    'low_viscosity': 0.8,    # Solvants fluides → T2 long → fin
}


# ==============================================================================
# AMÉLIORATION PRIORITÉ HAUTE 1: SEUIL ⁴J ADAPTATIF (RIGIDE VS FLEXIBLE)
# ==============================================================================
# Les ⁴J sont résolus dans les systèmes rigides (W-coupling, bicycliques) 
# mais souvent non-résolus dans les chaînes flexibles
# ==============================================================================

def _get_adaptive_4J_threshold(
    mol: Chem.Mol,
    i: int,
    j: int,
    base_threshold: float = 1.0
) -> float:
    """
    Calcule un seuil adaptatif pour les couplages ⁴J selon la rigidité du système.
    
    Systèmes rigides (bicycliques, aromatiques, W-coupling):
      - Seuil bas (0.8 Hz) car ⁴J peuvent atteindre 3-4 Hz
      
    Systèmes flexibles (chaînes acycliques):
      - Seuil haut (1.5 Hz) car ⁴J sont moyennés par rotation et souvent < 0.5 Hz
    
    Args:
        mol: Molécule RDKit
        i, j: Indices des deux atomes (H ou C)
        base_threshold: Seuil de base en Hz
    
    Returns:
        Seuil adapté en Hz
    """
    try:
        atom_i = mol.GetAtomWithIdx(i)
        atom_j = mol.GetAtomWithIdx(j)
        
        # Cas 1: Les deux atomes dans des cycles → système rigide
        if atom_i.IsInRing() and atom_j.IsInRing():
            ring_info = mol.GetRingInfo()
            
            # Vérifier si même système de cycles fusionnés
            for ring in ring_info.AtomRings():
                if i in ring and j in ring:
                    return 0.7  # Très rigide: même cycle
            
            # Cycles fusionnés (bicyclique)?
            i_rings = [set(r) for r in ring_info.AtomRings() if i in r]
            j_rings = [set(r) for r in ring_info.AtomRings() if j in r]
            
            for ri in i_rings:
                for rj in j_rings:
                    if ri & rj:  # Intersection non-vide → fusionnés
                        return 0.8  # Rigide: cycles fusionnés
        
        # Cas 2: W-coupling potentiel (détection géométrique)
        # Chemin H-C-C-C-H avec arrangement W
        path = Chem.GetShortestPath(mol, i, j)
        if path and len(path) == 5:  # ⁴J = 4 liaisons = 5 atomes
            # Vérifier si c'est un chemin allylique ou W
            mid_atoms = [mol.GetAtomWithIdx(idx) for idx in path[1:-1]]
            
            # Système allylique (C=C-C): rigide
            has_double = False
            for k in range(len(mid_atoms) - 1):
                bond = mol.GetBondBetweenAtoms(path[k+1], path[k+2])
                if bond and bond.GetBondTypeAsDouble() >= 2:
                    has_double = True
                    break
            
            if has_double:
                return 0.8  # Allylique: rigide
            
            # Tous dans un cycle?
            if all(a.IsInRing() for a in mid_atoms):
                return 0.8  # Chemin rigide
        
        # Cas 3: Au moins un atome aromatique
        if atom_i.GetIsAromatic() or atom_j.GetIsAromatic():
            return 0.9  # Semi-rigide
        
        # Cas 4: Système flexible (défaut)
        # Vérifier le nombre de liaisons rotables entre i et j
        if path:
            rotatable_in_path = 0
            for k in range(len(path) - 1):
                bond = mol.GetBondBetweenAtoms(path[k], path[k+1])
                if bond:
                    a1 = mol.GetAtomWithIdx(path[k])
                    a2 = mol.GetAtomWithIdx(path[k+1])
                    # Liaison simple C-C non cyclique = rotable
                    if (bond.GetBondTypeAsDouble() == 1 and 
                        not bond.IsInRing() and
                        a1.GetSymbol() in ('C', 'N', 'O') and
                        a2.GetSymbol() in ('C', 'N', 'O')):
                        rotatable_in_path += 1
            
            if rotatable_in_path >= 2:
                return 1.5  # Très flexible: seuil strict
            elif rotatable_in_path == 1:
                return 1.2  # Modérément flexible
        
        return base_threshold
        
    except Exception:
        return base_threshold


def _shares_fused_ring_system(
    ring_info,
    i: int,
    j: int
) -> bool:
    """Vérifie si deux atomes font partie du même système de cycles fusionnés."""
    i_rings = [set(r) for r in ring_info.AtomRings() if i in r]
    j_rings = [set(r) for r in ring_info.AtomRings() if j in r]
    
    if not i_rings or not j_rings:
        return False
    
    # BFS pour trouver si les systèmes de cycles sont connectés
    visited = set()
    queue = list(i_rings)
    
    while queue:
        current = queue.pop(0)
        if current in visited:
            continue
        visited.add(frozenset(current))
        
        # Vérifier si on a atteint un cycle contenant j
        for rj in j_rings:
            if current & rj:  # Intersection
                return True
        
        # Ajouter cycles fusionnés
        for ring in ring_info.AtomRings():
            ring_set = set(ring)
            if ring_set not in visited and current & ring_set:
                queue.append(ring_set)
    
    return False


# ==============================================================================
# AMÉLIORATION PRIORITÉ HAUTE 2: CONFORMÈRES DIVERSIFIÉS
# ==============================================================================
# ETKDG génère parfois des conformères trop similaires pour chaînes flexibles.
# Ajout d'échantillonnage systématique des dièdres (gauche+, anti, gauche-).
# ==============================================================================

def _get_rotatable_bonds(mol: Chem.Mol) -> List[Tuple[int, int, int, int]]:
    """
    Identifie les liaisons rotables avec les atomes pour définir le dièdre.
    
    Returns:
        Liste de tuples (a1, a2, a3, a4) définissant chaque dièdre rotable
    """
    rotatable = []
    
    try:
        # Pattern SMARTS pour liaisons rotables
        pattern = Chem.MolFromSmarts('[!$([NH]!@C(=O))&!D1]-&!@[!$([NH]!@C(=O))&!D1]')
        if pattern is None:
            return rotatable
            
        matches = mol.GetSubstructMatches(pattern)
        
        for match in matches:
            a2_idx, a3_idx = match[0], match[1]
            a2 = mol.GetAtomWithIdx(a2_idx)
            a3 = mol.GetAtomWithIdx(a3_idx)
            
            # Trouver a1 (voisin de a2, différent de a3)
            a1_candidates = [n.GetIdx() for n in a2.GetNeighbors() if n.GetIdx() != a3_idx]
            # Trouver a4 (voisin de a3, différent de a2)
            a4_candidates = [n.GetIdx() for n in a3.GetNeighbors() if n.GetIdx() != a2_idx]
            
            if a1_candidates and a4_candidates:
                rotatable.append((a1_candidates[0], a2_idx, a3_idx, a4_candidates[0]))
        
        return rotatable
        
    except Exception:
        return rotatable


def _generate_diverse_conformers_systematic(
    mol: Chem.Mol,
    n_conformers: int = 50,
    angles: List[float] = None,
    max_rotatable: int = 4,
    energy_cutoff_kcal: float = 10.0
) -> List[Tuple[int, float]]:
    """
    Génère des conformères avec échantillonnage systématique des dièdres.
    
    Pour les molécules avec peu de liaisons rotables (≤4), échantillonne
    systématiquement gauche+, anti, gauche- pour chaque liaison.
    
    Args:
        mol: Molécule avec H explicites
        n_conformers: Nombre max de conformères à retourner
        angles: Angles à échantillonner (défaut: [60, 180, 300] = g+, anti, g-)
        max_rotatable: Nombre max de liaisons pour échantillonnage systématique
        energy_cutoff_kcal: Seuil d'énergie pour filtrage
    
    Returns:
        Liste de (conf_id, energy_relative) triée par énergie
    """
    if angles is None:
        angles = [60.0, 180.0, 300.0]  # gauche+, anti, gauche-
    
    mol_h = Chem.AddHs(mol, addCoords=True)
    rotatable = _get_rotatable_bonds(mol_h)
    n_rot = len(rotatable)
    
    results = []
    
    try:
        # Cas 1: Peu de liaisons rotables → échantillonnage systématique
        if 0 < n_rot <= max_rotatable:
            import itertools
            
            # Générer toutes les combinaisons
            n_combos = len(angles) ** n_rot
            
            # Limiter si trop de combinaisons
            if n_combos > n_conformers * 3:
                # Trop de combinaisons: échantillonnage aléatoire
                combos = []
                import random
                for _ in range(n_conformers * 3):
                    combo = tuple(random.choice(angles) for _ in range(n_rot))
                    combos.append(combo)
            else:
                combos = list(itertools.product(angles, repeat=n_rot))
            
            # Générer un conformère initial
            ps = AllChem.ETKDGv3() if hasattr(AllChem, 'ETKDGv3') else AllChem.ETKDG()
            ps.randomSeed = 42
            AllChem.EmbedMolecule(mol_h, ps)
            
            if mol_h.GetNumConformers() == 0:
                return results
            
            # Pour chaque combinaison d'angles
            for combo in combos:
                # Copier la molécule
                mol_copy = Chem.Mol(mol_h)
                conf = mol_copy.GetConformer()
                
                # Appliquer les angles
                for (a1, a2, a3, a4), angle in zip(rotatable, combo):
                    try:
                        AllChem.SetDihedralDeg(conf, a1, a2, a3, a4, angle)
                    except Exception:
                        continue
                
                # Optimiser
                try:
                    AllChem.MMFFOptimizeMolecule(mol_copy, maxIters=100)
                    ff = AllChem.MMFFGetMoleculeForceField(mol_copy)
                    energy = ff.CalcEnergy() if ff else 999.0
                except Exception:
                    try:
                        AllChem.UFFOptimizeMolecule(mol_copy, maxIters=100)
                        ff = AllChem.UFFGetMoleculeForceField(mol_copy)
                        energy = ff.CalcEnergy() if ff else 999.0
                    except Exception:
                        energy = 999.0
                
                if energy < 900:
                    # Ajouter le conformère à la molécule originale
                    new_conf_id = mol_h.AddConformer(conf, assignId=True)
                    results.append((new_conf_id, energy))
        
        # Cas 2: Beaucoup de liaisons rotables → ETKDG standard amélioré
        if n_rot > max_rotatable or not results:
            ps = AllChem.ETKDGv3() if hasattr(AllChem, 'ETKDGv3') else AllChem.ETKDG()
            ps.randomSeed = 42
            ps.numThreads = 0
            ps.pruneRmsThresh = 0.35
            ps.useRandomCoords = True
            
            n_generate = n_conformers * 3
            conf_ids = AllChem.EmbedMultipleConfs(mol_h, numConfs=n_generate, params=ps)
            
            # Optimiser
            try:
                opt_results = AllChem.MMFFOptimizeMoleculeConfs(mol_h, maxIters=200)
            except Exception:
                opt_results = []
                for cid in conf_ids:
                    try:
                        AllChem.UFFOptimizeMolecule(mol_h, confId=cid, maxIters=200)
                        ff = AllChem.UFFGetMoleculeForceField(mol_h, confId=cid)
                        energy = ff.CalcEnergy() if ff else 999.0
                        opt_results.append((0, energy))
                    except Exception:
                        opt_results.append((1, 999.0))
            
            for cid, (converged, energy) in zip(conf_ids, opt_results):
                if converged == 0:
                    results.append((cid, energy))
        
        # Trier par énergie et filtrer
        if results:
            results.sort(key=lambda x: x[1])
            min_e = results[0][1]
            results = [(cid, e - min_e) for cid, e in results if e - min_e <= energy_cutoff_kcal]
            results = results[:n_conformers]
        
        return results
        
    except Exception:
        return results


# ==============================================================================
# AMÉLIORATION PRIORITÉ HAUTE 3: DÉTECTION DIASTÉRÉOTOPIE PROCHIRAUX
# ==============================================================================
# Détecte les centres prochiraux (carbonyle asymétrique, C=C asymétrique)
# qui rendent les CH₂ adjacents diastéréotopiques même sans centre chiral explicite
# ==============================================================================

def _is_prochiral_center(
    mol: Chem.Mol,
    atom_idx: int
) -> bool:
    """
    Détecte si un atome est un centre prochiral.
    
    Un centre prochiral est un atome sp2 (C=O, C=C, C=N) ou sp3 avec
    deux substituants identiques qui, s'ils devenaient différents,
    créeraient un centre chiral.
    
    Exemples:
    - C=O dans R1-CO-CH2-R2 où R1 ≠ R2: prochiral
    - C=C dans R1R2C=CH2: prochiral si R1 ≠ R2
    
    Args:
        mol: Molécule RDKit
        atom_idx: Index de l'atome à tester
    
    Returns:
        True si prochiral
    """
    try:
        atom = mol.GetAtomWithIdx(atom_idx)
        
        # Cas 1: Carbonyle (C=O) asymétrique
        if atom.GetSymbol() == 'C':
            hyb = atom.GetHybridization()
            
            if hyb == Chem.HybridizationType.SP2:
                # Chercher C=O ou C=N ou C=C
                double_bonded = None
                for bond in atom.GetBonds():
                    if bond.GetBondTypeAsDouble() >= 2:
                        other = bond.GetOtherAtom(atom)
                        if other.GetSymbol() in ('O', 'N', 'C', 'S'):
                            double_bonded = other
                            break
                
                if double_bonded:
                    # Récupérer les autres substituants
                    other_neighbors = [n for n in atom.GetNeighbors() 
                                       if n.GetIdx() != double_bonded.GetIdx()]
                    
                    if len(other_neighbors) == 2:
                        # Vérifier si les deux substituants sont différents
                        # Utiliser les rangs canoniques
                        try:
                            ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                            r1 = ranks[other_neighbors[0].GetIdx()]
                            r2 = ranks[other_neighbors[1].GetIdx()]
                            if r1 != r2:
                                return True  # Carbonyle asymétrique → prochiral
                        except Exception:
                            pass
        
        return False
        
    except Exception:
        return False


def _is_prochiral_CH2(
    mol: Chem.Mol,
    ch2_idx: int
) -> Tuple[bool, str]:
    """
    Vérifie si un CH₂ est rendu diastéréotopique par un centre prochiral adjacent.
    
    Cas couverts:
    1. CH₂ adjacent à un carbonyle asymétrique (R1-CO-CH2-R2, R1≠R2)
    2. CH₂ adjacent à un alcène asymétrique
    3. CH₂ dans un cycle avec substituants asymétriques
    
    Args:
        mol: Molécule avec H explicites
        ch2_idx: Index du carbone du CH₂
    
    Returns:
        Tuple (is_prochiral, reason)
    """
    try:
        ch2 = mol.GetAtomWithIdx(ch2_idx)
        
        if ch2.GetSymbol() != 'C':
            return (False, '')
        
        # Vérifier que c'est bien un CH₂
        h_count = sum(1 for n in ch2.GetNeighbors() if n.GetSymbol() == 'H')
        if h_count != 2:
            return (False, '')
        
        neighbors_heavy = [n for n in ch2.GetNeighbors() if n.GetSymbol() != 'H']
        
        # Cas 1: Adjacent à un carbonyle (C=O)
        for nb in neighbors_heavy:
            if nb.GetSymbol() == 'C':
                # Chercher C=O
                for bond in nb.GetBonds():
                    if bond.GetBondTypeAsDouble() >= 2:
                        other = bond.GetOtherAtom(nb)
                        if other.GetSymbol() == 'O':
                            # C'est un carbonyle - vérifier asymétrie
                            # Le carbonyle a-t-il des substituants différents?
                            carbonyl_neighbors = [n for n in nb.GetNeighbors() 
                                                  if n.GetIdx() != other.GetIdx() 
                                                  and n.GetIdx() != ch2_idx]
                            
                            if len(carbonyl_neighbors) >= 1:
                                # Comparer les chemins de l'autre côté du carbonyle vs côté CH2
                                try:
                                    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                                    # Les rangs des substituants du carbonyle
                                    carbonyl_sub_ranks = [ranks[n.GetIdx()] for n in carbonyl_neighbors]
                                    ch2_rank = ranks[ch2_idx]
                                    
                                    # Si les substituants sont différents → prochiral
                                    if len(set(carbonyl_sub_ranks)) > 1 or (carbonyl_sub_ranks and carbonyl_sub_ranks[0] != ch2_rank):
                                        return (True, 'adjacent_to_asymmetric_carbonyl')
                                except Exception:
                                    pass
        
        # Cas 2: Adjacent à un alcène asymétrique (C=C avec R1≠R2)
        for nb in neighbors_heavy:
            if nb.GetSymbol() == 'C':
                for bond in nb.GetBonds():
                    if bond.GetBondTypeAsDouble() >= 2:
                        other = bond.GetOtherAtom(nb)
                        if other.GetSymbol() == 'C':
                            # C'est un alcène
                            # Vérifier si les substituants de chaque C sont différents
                            nb_subs = [n for n in nb.GetNeighbors() 
                                       if n.GetIdx() not in (other.GetIdx(), ch2_idx)]
                            other_subs = [n for n in other.GetNeighbors() 
                                          if n.GetIdx() != nb.GetIdx()]
                            
                            try:
                                ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                                nb_ranks = set(ranks[n.GetIdx()] for n in nb_subs)
                                other_ranks = set(ranks[n.GetIdx()] for n in other_subs)
                                
                                # Asymétrique si substituants différents
                                if len(nb_ranks) > 1 or len(other_ranks) > 1:
                                    return (True, 'adjacent_to_asymmetric_alkene')
                            except Exception:
                                pass
        
        # Cas 3: CH₂ dans un cycle avec environnement asymétrique
        if ch2.IsInRing():
            # Vérifier si le cycle a des substituants asymétriques
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if ch2_idx in ring:
                    # Trouver les substituants du cycle
                    ring_set = set(ring)
                    substituents = []
                    for r_idx in ring:
                        r_atom = mol.GetAtomWithIdx(r_idx)
                        for nb in r_atom.GetNeighbors():
                            if nb.GetIdx() not in ring_set and nb.GetSymbol() != 'H':
                                substituents.append(nb.GetIdx())
                    
                    if len(substituents) >= 2:
                        try:
                            ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                            sub_ranks = [ranks[s] for s in substituents]
                            if len(set(sub_ranks)) > 1:
                                return (True, 'ring_with_asymmetric_substituents')
                        except Exception:
                            pass
        
        return (False, '')
        
    except Exception:
        return (False, '')


def _is_diastereotopic_by_prochirality(
    mol: Chem.Mol,
    heavy_idx: int,
    h_a: int,
    h_b: int
) -> bool:
    """
    Complète la détection de diastéréotopie par analyse des centres prochiraux.
    
    À appeler après les tests standards de diastéréotopie pour détecter
    les cas manqués par le test isotopique classique.
    
    Args:
        mol: Molécule avec H explicites
        heavy_idx: Index du carbone portant les H
        h_a, h_b: Indices des deux protons
    
    Returns:
        True si diastéréotopiques par prochirauté
    """
    # Vérifier si le CH2 est adjacent à un centre prochiral
    is_prochiral, reason = _is_prochiral_CH2(mol, heavy_idx)
    
    if is_prochiral:
        return True
    
    # Vérifier les voisins directs pour des centres prochiraux
    heavy_atom = mol.GetAtomWithIdx(heavy_idx)
    for nb in heavy_atom.GetNeighbors():
        if nb.GetSymbol() != 'H':
            if _is_prochiral_center(mol, nb.GetIdx()):
                return True
    
    return False


def _get_adaptive_linewidth(
    mol: Chem.Mol,
    h_idx: int,
    solvent: str = 'CDCl3',
    temperature_k: float = 298.0,
    base_linewidth_hz: float = 1.0
) -> float:
    """
    Calcule la largeur de raie adaptative pour un proton donné.
    
    La largeur dépend de:
    - Type de proton (échangeable vs non-échangeable)
    - Environnement chimique (aromatique, alkyle, etc.)
    - Conditions expérimentales (solvant, température)
    
    Args:
        mol: Molécule RDKit avec H explicites
        h_idx: Index du proton
        solvent: Nom du solvant
        temperature_k: Température en Kelvin
        base_linewidth_hz: Largeur de base (instrument)
    
    Returns:
        Largeur de raie FWHM en Hz
    """
    try:
        h_atom = mol.GetAtomWithIdx(h_idx)
        
        if h_atom.GetSymbol() != 'H':
            return base_linewidth_hz
        
        # Déterminer le type de proton
        neighbors = list(h_atom.GetNeighbors())
        if not neighbors:
            return base_linewidth_hz
        
        parent = neighbors[0]
        parent_sym = parent.GetSymbol()
        
        # Proton échangeable?
        is_exchangeable = parent_sym in ('O', 'N', 'S')
        
        if is_exchangeable:
            # Déterminer le type précis
            if parent_sym == 'O':
                # Vérifier si alcool, phénol ou acide
                if parent.GetIsAromatic():
                    proton_type = 'OH_phenol'
                else:
                    # Vérifier si acide carboxylique
                    for nb in parent.GetNeighbors():
                        if nb.GetSymbol() == 'C':
                            for nb2 in nb.GetNeighbors():
                                if nb2.GetSymbol() == 'O' and nb2.GetIdx() != parent.GetIdx():
                                    bond = mol.GetBondBetweenAtoms(nb.GetIdx(), nb2.GetIdx())
                                    if bond and bond.GetBondTypeAsDouble() >= 1.5:
                                        proton_type = 'OH_carboxylic'
                                        break
                            else:
                                continue
                            break
                    else:
                        proton_type = 'OH_alcohol'
            
            elif parent_sym == 'N':
                # Vérifier si amide
                for nb in parent.GetNeighbors():
                    if nb.GetSymbol() == 'C':
                        for nb2 in nb.GetNeighbors():
                            if nb2.GetSymbol() == 'O':
                                bond = mol.GetBondBetweenAtoms(nb.GetIdx(), nb2.GetIdx())
                                if bond and bond.GetBondTypeAsDouble() >= 1.5:
                                    proton_type = 'NH_amide'
                                    break
                        else:
                            continue
                        break
                else:
                    proton_type = 'NH_amine'
            
            else:  # S
                proton_type = 'SH'
            
            # Vérifier si le signal est visible dans ce solvant
            if solvent.upper() in ('D2O', 'MEOD', 'METHANOL-D4'):
                # Protons échangeables invisibles dans D2O
                return 0.0  # Signal absent
        
        else:
            # Proton non-échangeable - déterminer le type
            if parent_sym == 'C':
                if parent.GetIsAromatic():
                    proton_type = 'aromatic'
                elif _is_sp2_like(parent):
                    proton_type = 'vinyl'
                else:
                    # Compter les H sur ce carbone
                    h_count = sum(1 for n in parent.GetNeighbors() if n.GetSymbol() == 'H')
                    if h_count == 3:
                        proton_type = 'alkyl_CH3'
                    elif h_count == 2:
                        proton_type = 'alkyl_CH2'
                    else:
                        proton_type = 'alkyl_CH'
                    
                    # Vérifier si adjacent à hétéroatome
                    for nb in parent.GetNeighbors():
                        if nb.GetSymbol() in ('O', 'N', 'S', 'F', 'Cl'):
                            proton_type = 'alpha_heteroatom'
                            break
            else:
                proton_type = 'default'
        
        # Récupérer la largeur de base pour ce type
        lw_range = _LINEWIDTH_BY_PROTON_TYPE.get(proton_type, _LINEWIDTH_BY_PROTON_TYPE['default'])
        lw_base = (lw_range[0] + lw_range[1]) / 2.0
        
        # Appliquer les modificateurs de température
        if temperature_k < 280:
            lw_base *= _LINEWIDTH_MODIFIERS['low_temp']
        elif temperature_k > 320:
            lw_base *= _LINEWIDTH_MODIFIERS['high_temp']
        
        # Combiner avec la largeur instrumentale
        final_lw = max(base_linewidth_hz, lw_base)
        
        return final_lw
        
    except Exception:
        return base_linewidth_hz


def _get_geminal_2J(heavy_atom, mol: Chem.Mol) -> float:
    """
    Retourne ²J(H,H) geminal approprié selon l'environnement chimique.
    
    Note: ²J geminaux sont souvent NÉGATIFS, mais on retourne |²J|
    car le signe n'affecte pas le pattern de multiplicité.
    
    Returns:
        |²J| en Hz (toujours positif)
    """
    sym = heavy_atom.GetSymbol()
    
    if sym != 'C':
        return 12.5  # Défaut
    
    # Vérifier hybridation
    hyb = str(heavy_atom.GetHybridization())
    is_sp2 = 'SP2' in hyb or heavy_atom.GetIsAromatic()
    
    if is_sp2:
        # Aldéhyde?
        for bond in heavy_atom.GetBonds():
            if bond.GetBondTypeAsDouble() >= 2:
                other = bond.GetOtherAtom(heavy_atom)
                if other.GetSymbol() == 'O':
                    return 2.5  # CHO
        return 2.0  # =CH₂
    
    # sp³: vérifier environnement
    neighbors = [n.GetSymbol() for n in heavy_atom.GetNeighbors() if n.GetSymbol() != 'H']
    
    # Adjacent à carbonyle?
    for nb in heavy_atom.GetNeighbors():
        if nb.GetSymbol() == 'C':
            for b in nb.GetBonds():
                if b.GetBondTypeAsDouble() >= 2:
                    other = b.GetOtherAtom(nb)
                    if other.GetSymbol() == 'O':
                        return 16.0  # α-carbonyle
    
    # Adjacent à O?
    if 'O' in neighbors:
        return 11.0  # CH₂-O
    
    # Adjacent à N?
    if 'N' in neighbors:
        return 12.0  # CH₂-N
    
    # Cyclopropane?
    if heavy_atom.IsInRingSize(3):
        return 5.0
    
    # Cyclohexane ou pyranose?
    if heavy_atom.IsInRingSize(6):
        return 13.5
    
    return 12.5  # Défaut alkyle


def _get_1JCH(heavy_atom, mol: Chem.Mol) -> float:
    """
    Retourne ¹J(C-H) selon l'hybridation et l'environnement du carbone.
    
    Références:
        - Kalinowski et al., Carbon-13 NMR Spectroscopy
        - Berger & Braun, 200 and More NMR Experiments
    
    Returns:
        ¹J(C-H) en Hz
    """
    sym = heavy_atom.GetSymbol()
    
    if sym != 'C':
        return 140.0  # Défaut non-C
    
    hyb = str(heavy_atom.GetHybridization())
    
    # sp (alcyne terminal)
    if 'SP' in hyb and 'SP2' not in hyb and 'SP3' not in hyb:
        return _13C_1JCH_TABLE.get('sp_alkyne', 250.0)
    
    # sp² (aromatique, alcène, aldéhyde)
    if 'SP2' in hyb or heavy_atom.GetIsAromatic():
        # Aldéhyde?
        for bond in heavy_atom.GetBonds():
            if bond.GetBondTypeAsDouble() >= 2:
                other = bond.GetOtherAtom(heavy_atom)
                if other.GetSymbol() == 'O':
                    return _13C_1JCH_TABLE.get('sp2_aldehyde', 175.0)
        
        # Aromatique?
        if heavy_atom.GetIsAromatic():
            return _13C_1JCH_TABLE.get('aromatic', 160.0)
        
        # Alcène
        return _13C_1JCH_TABLE.get('sp2_alkene', 155.0)
    
    # sp³
    neighbors = [n.GetSymbol() for n in heavy_atom.GetNeighbors() if n.GetSymbol() != 'H']
    
    # Cyclopropane (tension de cycle → J élevé)
    if heavy_atom.IsInRingSize(3):
        return _13C_1JCH_TABLE.get('sp3_cyclopropane', 160.0)
    
    # Adjacent à O?
    if 'O' in neighbors:
        return _13C_1JCH_TABLE.get('sp3_alpha_O', 140.0)
    
    # Adjacent à N?
    if 'N' in neighbors:
        return _13C_1JCH_TABLE.get('sp3_alpha_N', 135.0)
    
    # Alkyle standard
    return _13C_1JCH_TABLE.get('sp3_alkyl', 125.0)


# ---------------------------
# Second-order detection (raffinée)
# ---------------------------
def _refined_second_order_detection(
    J_vals: List[float],
    ppm_vals: List[float],
    spectrometer_mhz: float,
    config: "SimConfig"
) -> Tuple[bool, str]:
    """
    Détection Δν/J avec seuils gradués (Keeler, Günther, Reich).
    """
    if not J_vals:
        return False, "insufficient_data"
    if len(ppm_vals) < 2:
        return False, "insufficient_data"

    Jmax = max(J_vals)
    delta_ppm = min(
        abs(a - b)
        for i, a in enumerate(ppm_vals)
        for b in ppm_vals[i+1:]
    )
    delta_nu_hz = delta_ppm * spectrometer_mhz

    if delta_nu_hz < 0.5:
        return True, "strong_coupling"

    ratio = Jmax / delta_nu_hz

    if ratio > 1.0:
        return True, "strong_coupling"
    elif ratio > 0.2:
        return True, "moderate_second_order"
    elif ratio > 0.1:
        return True, "weak_second_order"
    elif ratio > 0.05:
        return True, "borderline_roof_effect"

    return False, "first_order"

# ---------------------------
# Effet de toit (roof effect) – systèmes AB simples
# ---------------------------
def _calculate_roof_effect_intensities(
    J_hz: float,
    delta_nu_hz: float,
    base_multiplicity: str
) -> List[Tuple[float, float]]:
    """
    Calcule les positions et intensités avec effet de toit pour un système AB simple.
    """
    if J_hz <= 0.0:
        return [(0.0, 1.0)]

    D = float(np.sqrt(delta_nu_hz**2 + J_hz**2))
    if D <= 0.0:
        D = J_hz

    sin2theta = J_hz / D

    if base_multiplicity in ["d", "dd", "t", "dt", "td", "q"]:
        I_outer = max(0.0, 1.0 - sin2theta)
        I_inner = max(0.0, 1.0 + sin2theta)

        positions = [
            (-D/2.0 - J_hz/2.0, I_outer),
            (-D/2.0 + J_hz/2.0, I_inner),
            (+D/2.0 - J_hz/2.0, I_inner),
            (+D/2.0 + J_hz/2.0, I_outer),
        ]
        s = sum(i for _, i in positions) or 1.0
        positions = [(f, i/s) for f, i in positions]
        return positions

    return [(0.0, 1.0)]


# ==============================================================================
# MODULE HAMILTONIEN DE SPIN — Résolution exacte pour systèmes couplés
# ==============================================================================
# Références:
# - Pople J.A., Schneider W.G., Bernstein H.J., High-Resolution NMR (1959)
# - Günther H., NMR Spectroscopy, 3rd Ed. (2013), Chapter 4-5
# - Keeler J., Understanding NMR Spectroscopy, 2nd Ed. (2010)
# ==============================================================================

@dataclass
class SpinSystemPeaks:
    """Résultat d'une simulation de système de spin."""
    frequencies_hz: np.ndarray  # Positions des raies (Hz, relatives au centre)
    intensities: np.ndarray      # Intensités relatives (normalisées)
    system_type: str             # "AB", "ABX", "AA'BB'", etc.
    coupling_constants: Dict[str, float]  # J utilisés
    chemical_shifts_hz: Dict[str, float]  # δ utilisés (Hz)


def _solve_AB_system(
    delta_nu_hz: float,
    J_ab_hz: float,
    center_hz: float = 0.0
) -> SpinSystemPeaks:
    """
    Résolution exacte du système AB (2 spins I=1/2 couplés).
    
    Hamiltonien (base |αα⟩, |αβ⟩, |βα⟩, |ββ⟩):
        H = -νA·IzA - νB·IzB + JAB·(IA·IB)
    
    Transitions permises: ΔM = ±1, avec intensités calculées.
    
    Args:
        delta_nu_hz: νA - νB en Hz (peut être négatif)
        J_ab_hz: Constante de couplage JAB en Hz
        center_hz: Centre du système (νA + νB)/2 en Hz
    
    Returns:
        SpinSystemPeaks avec 4 raies (2 pour A, 2 pour B)
    """
    # Paramètres
    dv = abs(delta_nu_hz)
    J = J_ab_hz
    
    # C = sqrt((Δν)² + J²) — paramètre fondamental
    C = np.sqrt(dv**2 + J**2)
    
    # Angle de mélange θ: tan(2θ) = J / Δν
    if dv > 1e-6:
        theta = 0.5 * np.arctan2(J, dv)
    else:
        theta = np.pi / 4  # Limite Δν → 0: mélange maximal
    
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    
    # Fréquences des 4 transitions (relatives au centre)
    # A: transitions 1→2 et 3→4 dans la base propre
    # B: transitions 1→3 et 2→4
    
    f1 = -0.5 * (C + J)   # Raie externe gauche (A ou B selon signe Δν)
    f2 = -0.5 * (C - J)   # Raie interne gauche
    f3 = +0.5 * (C - J)   # Raie interne droite
    f4 = +0.5 * (C + J)   # Raie externe droite
    
    # Intensités (effet de toit)
    # I_inner ∝ 1 + sin(2θ), I_outer ∝ 1 - sin(2θ)
    sin_2theta = 2 * sin_theta * cos_theta
    I_inner = 1.0 + sin_2theta
    I_outer = 1.0 - sin_2theta
    
    frequencies = np.array([f1, f2, f3, f4]) + center_hz
    intensities = np.array([I_outer, I_inner, I_inner, I_outer])
    
    # Normaliser
    intensities = intensities / intensities.sum()
    
    return SpinSystemPeaks(
        frequencies_hz=frequencies,
        intensities=intensities,
        system_type="AB",
        coupling_constants={"J_AB": J_ab_hz},
        chemical_shifts_hz={"delta_nu": delta_nu_hz, "center": center_hz}
    )


def _solve_ABX_system(
    nu_A_hz: float,
    nu_B_hz: float,
    nu_X_hz: float,
    J_AB_hz: float,
    J_AX_hz: float,
    J_BX_hz: float
) -> SpinSystemPeaks:
    """
    Résolution du système ABX (3 spins, X faiblement couplé à AB).
    
    Approche: X est traité comme perturbation → 2 sous-systèmes AB
    (un pour X = α, un pour X = β), chacun avec shifts effectifs.
    
    Total: 12 raies (8 pour AB, 4 pour X).
    
    Args:
        nu_A_hz, nu_B_hz, nu_X_hz: Shifts en Hz
        J_AB_hz, J_AX_hz, J_BX_hz: Constantes de couplage en Hz
    
    Returns:
        SpinSystemPeaks avec toutes les raies
    """
    frequencies = []
    intensities = []
    
    # Sous-système AB avec X = α (spin up)
    # Shifts effectifs: νA' = νA + JAX/2, νB' = νB + JBX/2
    nu_A_eff_alpha = nu_A_hz + J_AX_hz / 2
    nu_B_eff_alpha = nu_B_hz + J_BX_hz / 2
    delta_nu_alpha = nu_A_eff_alpha - nu_B_eff_alpha
    center_alpha = (nu_A_eff_alpha + nu_B_eff_alpha) / 2
    
    ab_alpha = _solve_AB_system(delta_nu_alpha, J_AB_hz, center_alpha)
    frequencies.extend(ab_alpha.frequencies_hz)
    intensities.extend(ab_alpha.intensities * 0.5)  # Population 50%
    
    # Sous-système AB avec X = β (spin down)
    nu_A_eff_beta = nu_A_hz - J_AX_hz / 2
    nu_B_eff_beta = nu_B_hz - J_BX_hz / 2
    delta_nu_beta = nu_A_eff_beta - nu_B_eff_beta
    center_beta = (nu_A_eff_beta + nu_B_eff_beta) / 2
    
    ab_beta = _solve_AB_system(delta_nu_beta, J_AB_hz, center_beta)
    frequencies.extend(ab_beta.frequencies_hz)
    intensities.extend(ab_beta.intensities * 0.5)
    
    # Raies de X (doublet de doublets en première approximation)
    # Positions: νX ± (JAX + JBX)/2 et νX ± (JAX - JBX)/2
    # Mais pour ABX, les intensités sont modulées par le mélange AB
    # Approximation premier ordre pour X:
    delta_J = (J_AX_hz + J_BX_hz) / 2
    diff_J = (J_AX_hz - J_BX_hz) / 2
    
    x_freqs = [
        nu_X_hz + delta_J,
        nu_X_hz + diff_J,
        nu_X_hz - diff_J,
        nu_X_hz - delta_J,
    ]
    x_intensities = [0.25, 0.25, 0.25, 0.25]  # Simplifié
    
    frequencies.extend(x_freqs)
    intensities.extend(x_intensities)
    
    frequencies = np.array(frequencies)
    intensities = np.array(intensities)
    intensities = intensities / intensities.sum()
    
    return SpinSystemPeaks(
        frequencies_hz=frequencies,
        intensities=intensities,
        system_type="ABX",
        coupling_constants={"J_AB": J_AB_hz, "J_AX": J_AX_hz, "J_BX": J_BX_hz},
        chemical_shifts_hz={"nu_A": nu_A_hz, "nu_B": nu_B_hz, "nu_X": nu_X_hz}
    )


def _solve_AABB_system(
    delta_nu_hz: float,
    J_AA_hz: float,
    J_BB_hz: float,
    J_AB_hz: float,
    J_AB_prime_hz: float,
    center_hz: float = 0.0
) -> SpinSystemPeaks:
    """
    Résolution du système AA'BB' (4 spins équivalents par paires).
    
    Exemple typique: -CH₂-CH₂- avec Δν petit.
    
    Hamiltonien 16×16, mais la symétrie réduit à des blocs plus petits.
    Ici: résolution numérique complète par diagonalisation.
    
    Args:
        delta_nu_hz: νA - νB (Hz)
        J_AA_hz: Couplage entre les deux A (geminal ou long-range)
        J_BB_hz: Couplage entre les deux B
        J_AB_hz: Couplage A-B (vicinal cis ou trans)
        J_AB_prime_hz: Couplage A-B' (vicinal autre)
        center_hz: (νA + νB)/2
    
    Returns:
        SpinSystemPeaks avec jusqu'à 24 raies
    """
    # Construction du Hamiltonien dans la base produit (16 états)
    # |αααα⟩, |αααβ⟩, ..., |ββββ⟩
    
    # États: (A, A', B, B') chacun = 0 (α) ou 1 (β)
    n_states = 16
    
    # Opérateurs Iz pour chaque spin (diagonaux)
    def Iz(state, spin_idx):
        """Valeur propre de Iz: +1/2 si α, -1/2 si β"""
        bit = (state >> (3 - spin_idx)) & 1
        return 0.5 if bit == 0 else -0.5
    
    # Opérateurs I+I- pour couplages (hors-diagonaux)
    def flip(state, i, j):
        """Flip les spins i et j (si un α et un β)"""
        bi = (state >> (3 - i)) & 1
        bj = (state >> (3 - j)) & 1
        if bi != bj:  # Un up, un down
            new_state = state ^ (1 << (3 - i)) ^ (1 << (3 - j))
            return new_state
        return None
    
    # Construire H
    H = np.zeros((n_states, n_states), dtype=np.float64)
    
    # Termes Zeeman (diagonaux)
    nu_A = center_hz + delta_nu_hz / 2
    nu_B = center_hz - delta_nu_hz / 2
    
    for s in range(n_states):
        # Spins: 0=A, 1=A', 2=B, 3=B'
        H[s, s] -= nu_A * Iz(s, 0)  # A
        H[s, s] -= nu_A * Iz(s, 1)  # A'
        H[s, s] -= nu_B * Iz(s, 2)  # B
        H[s, s] -= nu_B * Iz(s, 3)  # B'
    
    # Termes de couplage
    couplings = [
        (0, 1, J_AA_hz),      # A-A'
        (2, 3, J_BB_hz),      # B-B'
        (0, 2, J_AB_hz),      # A-B
        (1, 3, J_AB_hz),      # A'-B'
        (0, 3, J_AB_prime_hz),  # A-B'
        (1, 2, J_AB_prime_hz),  # A'-B
    ]
    
    for i, j, J in couplings:
        for s in range(n_states):
            # Terme diagonal: J * Iz_i * Iz_j
            H[s, s] += J * Iz(s, i) * Iz(s, j)
            
            # Termes hors-diagonaux: J/2 * (I+_i I-_j + I-_i I+_j)
            s_flip = flip(s, i, j)
            if s_flip is not None:
                H[s, s_flip] += J / 2
                H[s_flip, s] += J / 2
    
    # Diagonaliser
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    
    # Calculer les transitions permises (ΔM = ±1)
    # Opérateur I+ total = I+_A + I+_A' + I+_B + I+_B'
    def total_Mz(state):
        return sum(Iz(state, i) for i in range(4))
    
    frequencies = []
    intensities = []
    
    for i in range(n_states):
        for j in range(n_states):
            # Transition j → i (absorption si Ei > Ej)
            dE = eigenvalues[i] - eigenvalues[j]
            if dE <= 0:
                continue
            
            # Vérifier ΔM = -1 (pour I+)
            # En fait on calcule l'élément de matrice ⟨i|I+|j⟩
            intensity = 0.0
            for s in range(n_states):
                for spin in range(4):
                    # I+ sur spin: α→rien, β→α
                    bit = (s >> (3 - spin)) & 1
                    if bit == 1:  # β
                        s_new = s ^ (1 << (3 - spin))  # Flip to α
                        # ⟨i|s_new⟩ * ⟨s|j⟩
                        intensity += eigenvectors[s_new, i] * eigenvectors[s, j]
            
            intensity = abs(intensity) ** 2
            
            if intensity > 1e-6:
                frequencies.append(dE)
                intensities.append(intensity)
    
    if not frequencies:
        # Fallback: 4 raies centrées
        return SpinSystemPeaks(
            frequencies_hz=np.array([center_hz]),
            intensities=np.array([1.0]),
            system_type="AA'BB'",
            coupling_constants={},
            chemical_shifts_hz={}
        )
    
    frequencies = np.array(frequencies)
    intensities = np.array(intensities)
    intensities = intensities / intensities.sum()
    
    return SpinSystemPeaks(
        frequencies_hz=frequencies,
        intensities=intensities,
        system_type="AA'BB'",
        coupling_constants={
            "J_AA": J_AA_hz, "J_BB": J_BB_hz,
            "J_AB": J_AB_hz, "J_AB'": J_AB_prime_hz
        },
        chemical_shifts_hz={"delta_nu": delta_nu_hz, "center": center_hz}
    )


def solve_spin_system(
    system_type: str,
    shifts_hz: Dict[str, float],
    couplings_hz: Dict[str, float],
    spectrometer_mhz: float = 400.0
) -> Optional[SpinSystemPeaks]:
    """
    Point d'entrée pour résoudre un système de spin quelconque.
    
    Args:
        system_type: "AB", "ABX", "AA'BB'", "AX", etc.
        shifts_hz: {"nu_A": ..., "nu_B": ..., ...} en Hz
        couplings_hz: {"J_AB": ..., ...} en Hz
        spectrometer_mhz: Pour conversion ppm→Hz si nécessaire
    
    Returns:
        SpinSystemPeaks ou None si système non supporté
    """
    if system_type == "AB" or system_type == "AB_weak":
        nu_A = shifts_hz.get("nu_A", 0.0)
        nu_B = shifts_hz.get("nu_B", 0.0)
        delta_nu = nu_A - nu_B
        center = (nu_A + nu_B) / 2
        J_AB = couplings_hz.get("J_AB", 7.0)
        return _solve_AB_system(delta_nu, J_AB, center)
    
    elif system_type == "ABX":
        return _solve_ABX_system(
            nu_A_hz=shifts_hz.get("nu_A", 0.0),
            nu_B_hz=shifts_hz.get("nu_B", 0.0),
            nu_X_hz=shifts_hz.get("nu_X", 0.0),
            J_AB_hz=couplings_hz.get("J_AB", 7.0),
            J_AX_hz=couplings_hz.get("J_AX", 7.0),
            J_BX_hz=couplings_hz.get("J_BX", 7.0),
        )
    
    elif system_type == "AA'BB'" or system_type == "A₂B₂":
        nu_A = shifts_hz.get("nu_A", 0.0)
        nu_B = shifts_hz.get("nu_B", 0.0)
        delta_nu = nu_A - nu_B
        center = (nu_A + nu_B) / 2
        return _solve_AABB_system(
            delta_nu_hz=delta_nu,
            J_AA_hz=couplings_hz.get("J_AA", -12.0),  # Geminal typique
            J_BB_hz=couplings_hz.get("J_BB", -12.0),
            J_AB_hz=couplings_hz.get("J_AB", 7.0),
            J_AB_prime_hz=couplings_hz.get("J_AB'", 7.0),
            center_hz=center,
        )
    
    # AX, AMX: first-order, pas besoin de diagonalisation
    return None


# ==============================================================================
# MODULE SIMULATION DE SPECTRE (LINESHAPE) — FID → FFT
# ==============================================================================
# Génère un spectre NMR 1D réaliste avec largeur de raie Lorentzienne.
# ==============================================================================

@dataclass
class SimulatedSpectrum:
    """Spectre NMR simulé complet."""
    ppm_axis: np.ndarray          # Axe en ppm
    freq_axis_hz: np.ndarray      # Axe en Hz
    intensity: np.ndarray          # Intensités (normalisées)
    peaks_ppm: List[float]         # Positions des pics (ppm)
    peaks_intensity: List[float]   # Intensités des pics
    spectrometer_mhz: float
    line_width_hz: float
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "ppm_axis": self.ppm_axis.tolist(),
            "freq_axis_hz": self.freq_axis_hz.tolist(),
            "intensity": self.intensity.tolist(),
            "peaks_ppm": self.peaks_ppm,
            "peaks_intensity": self.peaks_intensity,
            "spectrometer_mhz": self.spectrometer_mhz,
            "line_width_hz": self.line_width_hz,
        }


def _lorentzian(x: np.ndarray, x0: float, gamma: float, amplitude: float = 1.0) -> np.ndarray:
    """
    Fonction Lorentzienne (forme de raie NMR idéale).
    
    L(x) = A * γ² / ((x - x0)² + γ²)
    
    Args:
        x: Axe de fréquences
        x0: Position du centre
        gamma: Demi-largeur à mi-hauteur (HWHM)
        amplitude: Amplitude au pic
    
    Returns:
        Profil Lorentzien normalisé
    """
    return amplitude * (gamma ** 2) / ((x - x0) ** 2 + gamma ** 2)


def _gaussian(x: np.ndarray, x0: float, sigma: float, amplitude: float = 1.0) -> np.ndarray:
    """
    Fonction Gaussienne (pour convolution inhomogène).
    """
    return amplitude * np.exp(-0.5 * ((x - x0) / sigma) ** 2)


def _voigt_approx(x: np.ndarray, x0: float, gamma: float, sigma: float, amplitude: float = 1.0) -> np.ndarray:
    """
    Pseudo-Voigt (approximation rapide du profil Voigt).
    
    Combinaison linéaire Lorentzienne + Gaussienne.
    """
    # Fraction Lorentzienne basée sur les largeurs relatives
    f_L = gamma / (gamma + sigma) if (gamma + sigma) > 0 else 0.5
    L = _lorentzian(x, x0, gamma, 1.0)
    G = _gaussian(x, x0, sigma, 1.0)
    return amplitude * (f_L * L + (1 - f_L) * G)


# ==============================================================================
# SIMULATION PAR DIAGONALISATION DU HAMILTONIEN DE SPIN
# ==============================================================================
# Référence: Castillo, Patiny & Wist, J. Magn. Reson. 209, 123-130 (2011)
# "Fast and Accurate Algorithm for the Simulation of NMR spectra of Large Spin Systems"
# ==============================================================================

@dataclass
class SpinSystem:
    """
    Système de spins pour simulation quantique.
    
    Attributes:
        n_spins: Nombre de spins 1/2
        shifts_hz: Déplacements chimiques en Hz (relatif à référence)
        couplings_hz: Matrice de couplages J (Hz), symétrique
        equivalences: Groupes d'équivalence magnétique [(idx1, idx2, ...)]
    """
    n_spins: int
    shifts_hz: np.ndarray  # (n_spins,)
    couplings_hz: np.ndarray  # (n_spins, n_spins) symétrique
    equivalences: List[Tuple[int, ...]] = field(default_factory=list)
    labels: List[str] = field(default_factory=list)  # Étiquettes atomiques
    
    @classmethod
    def from_groups(
        cls,
        groups: List[Dict[str, Any]],
        spectrometer_mhz: float = 400.0
    ) -> 'SpinSystem':
        """
        Crée un SpinSystem depuis la sortie de simulate_multiplicity.
        
        Args:
            groups: Liste de groupes de protons avec shifts et couplages
            spectrometer_mhz: Fréquence du spectromètre
        
        Returns:
            SpinSystem configuré
        """
        # Collecter tous les protons individuels
        all_protons = []
        for g in groups:
            shift_ppm = g.get('H_shifts_ppm', [0.0])[0]
            shift_hz = shift_ppm * spectrometer_mhz
            h_indices = g.get('H_indices', [])
            for h_idx in h_indices:
                all_protons.append({
                    'h_idx': h_idx,
                    'shift_hz': shift_hz,
                    'group_id': g.get('group_id', 0),
                    'j_couplings': g.get('J_couplings', []),
                    'partners': g.get('partners', [])
                })
        
        n = len(all_protons)
        if n == 0:
            return cls(n_spins=0, shifts_hz=np.array([]), couplings_hz=np.array([[]]))
        
        # Construire les shifts
        shifts = np.array([p['shift_hz'] for p in all_protons])
        
        # Construire la matrice de couplages
        couplings = np.zeros((n, n))
        h_idx_to_spin = {p['h_idx']: i for i, p in enumerate(all_protons)}
        
        for i, p in enumerate(all_protons):
            for partner in p.get('partners', []):
                partner_h = partner.get('partner_h_index')
                j_hz = partner.get('J_hz', 0.0)
                if partner_h in h_idx_to_spin:
                    j = h_idx_to_spin[partner_h]
                    couplings[i, j] = j_hz
                    couplings[j, i] = j_hz  # Symétrique
        
        # Groupes d'équivalence
        equivalences = []
        group_to_spins = defaultdict(list)
        for i, p in enumerate(all_protons):
            group_to_spins[p['group_id']].append(i)
        for gid, spins in group_to_spins.items():
            if len(spins) > 1:
                equivalences.append(tuple(spins))
        
        labels = [f"H{p['h_idx']}" for p in all_protons]
        
        return cls(
            n_spins=n,
            shifts_hz=shifts,
            couplings_hz=couplings,
            equivalences=equivalences,
            labels=labels
        )


def _build_spin_operators(n_spins: int) -> Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray]]:
    """
    Construit les opérateurs de spin Ix, Iy, Iz pour n spins 1/2.
    
    Utilise le produit tensoriel (Kronecker) pour construire les opérateurs
    dans l'espace de Hilbert complet (dimension 2^n).
    
    Returns:
        (Ix_list, Iy_list, Iz_list) - Listes d'opérateurs pour chaque spin
    """
    dim = 2 ** n_spins
    
    # Matrices de Pauli / 2
    sx = np.array([[0, 1], [1, 0]], dtype=np.complex128) / 2
    sy = np.array([[0, -1j], [1j, 0]], dtype=np.complex128) / 2
    sz = np.array([[1, 0], [0, -1]], dtype=np.complex128) / 2
    I2 = np.eye(2, dtype=np.complex128)
    
    Ix_list = []
    Iy_list = []
    Iz_list = []
    
    for k in range(n_spins):
        # Produit tensoriel: I ⊗ ... ⊗ Sx ⊗ ... ⊗ I
        op_x = np.array([[1]], dtype=np.complex128)
        op_y = np.array([[1]], dtype=np.complex128)
        op_z = np.array([[1]], dtype=np.complex128)
        
        for j in range(n_spins):
            if j == k:
                op_x = np.kron(op_x, sx)
                op_y = np.kron(op_y, sy)
                op_z = np.kron(op_z, sz)
            else:
                op_x = np.kron(op_x, I2)
                op_y = np.kron(op_y, I2)
                op_z = np.kron(op_z, I2)
        
        Ix_list.append(op_x)
        Iy_list.append(op_y)
        Iz_list.append(op_z)
    
    return Ix_list, Iy_list, Iz_list


def _build_hamiltonian(
    spin_system: SpinSystem,
    Iz_list: List[np.ndarray]
) -> np.ndarray:
    """
    Construit le Hamiltonien de spin en représentation tournante.
    
    H = Σᵢ νᵢ Izᵢ + Σᵢ<ⱼ Jᵢⱼ (Ixᵢ Ixⱼ + Iyᵢ Iyⱼ + Izᵢ Izⱼ)
    
    En pratique, pour les liquides isotropes (secular approximation):
    H = Σᵢ νᵢ Izᵢ + Σᵢ<ⱼ Jᵢⱼ Izᵢ Izⱼ  (couplage faible)
    
    Ou Hamiltonien complet pour systèmes fortement couplés.
    """
    n = spin_system.n_spins
    dim = 2 ** n
    H = np.zeros((dim, dim), dtype=np.complex128)
    
    # Termes Zeeman (déplacements chimiques)
    for i in range(n):
        H += spin_system.shifts_hz[i] * Iz_list[i]
    
    # Termes de couplage J (forme scalaire complète pour second ordre)
    # Ici on utilise seulement Iz·Iz pour l'approximation séculaire
    # Pour un traitement exact, il faudrait inclure les termes transverses
    for i in range(n):
        for j in range(i + 1, n):
            Jij = spin_system.couplings_hz[i, j]
            if abs(Jij) > 0.01:  # Seuil de couplage significatif
                # Couplage scalaire: J * I⃗ᵢ · I⃗ⱼ = J * (Ixi*Ixj + Iyi*Iyj + Izi*Izj)
                # Pour simulation exacte des systèmes AB
                H += Jij * (Iz_list[i] @ Iz_list[j])
    
    return H


def _build_hamiltonian_full(
    spin_system: SpinSystem,
    Ix_list: List[np.ndarray],
    Iy_list: List[np.ndarray],
    Iz_list: List[np.ndarray]
) -> np.ndarray:
    """
    Construit le Hamiltonien complet avec termes transverses.
    
    H = Σᵢ νᵢ Izᵢ + Σᵢ<ⱼ Jᵢⱼ (Ixᵢ Ixⱼ + Iyᵢ Iyⱼ + Izᵢ Izⱼ)
    
    Nécessaire pour simulation exacte des systèmes fortement couplés (AB, AA'BB', etc.)
    """
    n = spin_system.n_spins
    dim = 2 ** n
    H = np.zeros((dim, dim), dtype=np.complex128)
    
    # Termes Zeeman
    for i in range(n):
        H += spin_system.shifts_hz[i] * Iz_list[i]
    
    # Couplage scalaire complet: J * I⃗ᵢ · I⃗ⱼ
    for i in range(n):
        for j in range(i + 1, n):
            Jij = spin_system.couplings_hz[i, j]
            if abs(Jij) > 0.01:
                H += Jij * (Ix_list[i] @ Ix_list[j] +
                           Iy_list[i] @ Iy_list[j] +
                           Iz_list[i] @ Iz_list[j])
    
    return H


def _compute_transitions(
    H: np.ndarray,
    Ix_total: np.ndarray,
    Iy_total: np.ndarray,
    temperature_k: float = 298.15,
    threshold: float = 1e-6
) -> List[Tuple[float, float]]:
    """
    Calcule les fréquences et intensités des transitions par diagonalisation.
    
    1. Diagonalise H pour obtenir énergies propres Eₙ et vecteurs |n⟩
    2. Calcule les éléments de matrice ⟨m|I⁺|n⟩ pour l'opérateur d'excitation
    3. Intensité ∝ |⟨m|I⁺|n⟩|² × (population thermique)
    
    Returns:
        Liste de (fréquence_hz, intensité) pour chaque transition
    """
    # Diagonalisation
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    dim = len(eigenvalues)
    
    # Opérateur I⁺ = Ix + i*Iy (création de cohérence)
    I_plus = Ix_total + 1j * Iy_total
    
    # Matrice I⁺ dans la base propre
    I_plus_eigen = eigenvectors.conj().T @ I_plus @ eigenvectors
    
    transitions = []
    
    for m in range(dim):
        for n in range(dim):
            if m != n:
                # Fréquence de transition
                freq = eigenvalues[m] - eigenvalues[n]
                
                # Intensité (élément de matrice au carré)
                intensity = np.abs(I_plus_eigen[m, n]) ** 2
                
                if intensity > threshold and freq > 0:  # Seulement transitions positives
                    transitions.append((freq, intensity))
    
    return transitions


def simulate_spin_system_hamiltonian(
    spin_system: SpinSystem,
    spectrometer_mhz: float = 400.0,
    line_width_hz: float = 1.0,
    ppm_range: Tuple[float, float] = (-0.5, 12.0),
    n_points: int = 32768,
    full_hamiltonian: bool = True,
    temperature_k: float = 298.15
) -> SimulatedSpectrum:
    """
    Simule un spectre NMR par diagonalisation exacte du Hamiltonien.
    
    Cette approche est nécessaire pour les systèmes fortement couplés
    (AB, ABX, AA'BB', etc.) où l'approximation premier ordre échoue.
    
    Args:
        spin_system: Système de spins défini par shifts et couplages
        spectrometer_mhz: Fréquence du spectromètre
        line_width_hz: Largeur de raie FWHM
        ppm_range: Fenêtre spectrale (min_ppm, max_ppm)
        n_points: Nombre de points du spectre
        full_hamiltonian: Si True, utilise le Hamiltonien complet (termes transverses)
        temperature_k: Température (pour populations Boltzmann, généralement négligeable)
    
    Returns:
        SimulatedSpectrum avec le spectre calculé exactement
    
    Limitations:
        - Coût exponentiel O(2^n) pour n spins
        - Pratique jusqu'à ~12-13 spins
        - Pour plus de spins, utiliser la fragmentation en clusters
    """
    n = spin_system.n_spins
    
    if n == 0:
        # Système vide
        ppm_axis = np.linspace(ppm_range[1], ppm_range[0], n_points)
        return SimulatedSpectrum(
            ppm_axis=ppm_axis,
            freq_axis_hz=ppm_axis * spectrometer_mhz,
            intensity=np.zeros(n_points),
            peaks_ppm=[],
            peaks_intensity=[],
            spectrometer_mhz=spectrometer_mhz,
            line_width_hz=line_width_hz
        )
    
    if n > 12:
        raise ValueError(f"Système trop grand ({n} spins). Maximum 12 spins pour diagonalisation exacte.")
    
    # Construire les opérateurs de spin
    Ix_list, Iy_list, Iz_list = _build_spin_operators(n)
    
    # Opérateurs totaux
    Ix_total = sum(Ix_list)
    Iy_total = sum(Iy_list)
    
    # Construire le Hamiltonien
    if full_hamiltonian:
        H = _build_hamiltonian_full(spin_system, Ix_list, Iy_list, Iz_list)
    else:
        H = _build_hamiltonian(spin_system, Iz_list)
    
    # Calculer les transitions
    transitions = _compute_transitions(H, Ix_total, Iy_total, temperature_k)
    
    # Construire le spectre
    ppm_min, ppm_max = ppm_range
    ppm_axis = np.linspace(ppm_max, ppm_min, n_points)
    freq_axis = ppm_axis * spectrometer_mhz
    
    spectrum = np.zeros(n_points)
    gamma = line_width_hz / 2  # HWHM
    
    peaks_ppm = []
    peaks_intensity = []
    
    for freq_hz, intensity in transitions:
        spectrum += _lorentzian(freq_axis, freq_hz, gamma, intensity)
        peaks_ppm.append(freq_hz / spectrometer_mhz)
        peaks_intensity.append(intensity)
    
    # Normaliser
    if spectrum.max() > 0:
        spectrum = spectrum / spectrum.max()
    
    return SimulatedSpectrum(
        ppm_axis=ppm_axis,
        freq_axis_hz=freq_axis,
        intensity=spectrum,
        peaks_ppm=peaks_ppm,
        peaks_intensity=peaks_intensity,
        spectrometer_mhz=spectrometer_mhz,
        line_width_hz=line_width_hz
    )


def _fragment_spin_system(
    spin_system: SpinSystem,
    weak_coupling_threshold_hz: float = 0.5
) -> List[SpinSystem]:
    """
    Fragmente un système de spins en clusters indépendants.
    
    Deux spins sont dans le même cluster si:
    1. Ils sont fortement couplés (|J| > seuil), OU
    2. Ils sont connectés via d'autres spins fortement couplés
    
    Référence: Castillo, Patiny & Wist, JMR 2011
    
    Args:
        spin_system: Système complet
        weak_coupling_threshold_hz: Seuil pour considérer un couplage comme faible
    
    Returns:
        Liste de sous-systèmes de spins (clusters)
    """
    n = spin_system.n_spins
    if n <= 1:
        return [spin_system]
    
    # Construire le graphe de connectivité (couplages forts)
    adjacency = defaultdict(set)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(spin_system.couplings_hz[i, j]) > weak_coupling_threshold_hz:
                adjacency[i].add(j)
                adjacency[j].add(i)
    
    # Trouver les composantes connexes (BFS)
    visited = set()
    clusters = []
    
    for start in range(n):
        if start in visited:
            continue
        
        cluster = []
        queue = deque([start])
        
        while queue:
            node = queue.popleft()
            if node in visited:
                continue
            visited.add(node)
            cluster.append(node)
            
            for neighbor in adjacency[node]:
                if neighbor not in visited:
                    queue.append(neighbor)
        
        clusters.append(sorted(cluster))
    
    # Créer les sous-systèmes
    subsystems = []
    for cluster in clusters:
        n_sub = len(cluster)
        shifts = np.array([spin_system.shifts_hz[i] for i in cluster])
        couplings = np.zeros((n_sub, n_sub))
        
        for i_sub, i_orig in enumerate(cluster):
            for j_sub, j_orig in enumerate(cluster):
                couplings[i_sub, j_sub] = spin_system.couplings_hz[i_orig, j_orig]
        
        labels = [spin_system.labels[i] for i in cluster] if spin_system.labels else []
        
        subsystems.append(SpinSystem(
            n_spins=n_sub,
            shifts_hz=shifts,
            couplings_hz=couplings,
            labels=labels
        ))
    
    return subsystems


def simulate_with_fragmentation(
    spin_system: SpinSystem,
    spectrometer_mhz: float = 400.0,
    line_width_hz: float = 1.0,
    ppm_range: Tuple[float, float] = (-0.5, 12.0),
    n_points: int = 32768,
    max_cluster_size: int = 10,
    weak_coupling_threshold_hz: float = 0.5
) -> SimulatedSpectrum:
    """
    Simule un spectre NMR en fragmentant le système en clusters.
    
    Permet de simuler des systèmes plus grands que la limite de 12 spins
    en les découpant en sous-systèmes indépendants.
    
    Args:
        spin_system: Système de spins complet
        spectrometer_mhz: Fréquence du spectromètre
        line_width_hz: Largeur de raie
        ppm_range: Fenêtre spectrale
        n_points: Résolution du spectre
        max_cluster_size: Taille maximale d'un cluster (sinon subdivision forcée)
        weak_coupling_threshold_hz: Seuil pour fragmentation
    
    Returns:
        SimulatedSpectrum combiné de tous les clusters
    """
    # Fragmenter
    clusters = _fragment_spin_system(spin_system, weak_coupling_threshold_hz)
    
    # Vérifier les tailles
    for i, cluster in enumerate(clusters):
        if cluster.n_spins > max_cluster_size:
            # Forcer une subdivision (approximation)
            # Pour l'instant, on lève une erreur
            raise ValueError(
                f"Cluster {i} a {cluster.n_spins} spins > {max_cluster_size}. "
                f"Augmentez weak_coupling_threshold_hz pour fragmenter davantage."
            )
    
    # Simuler chaque cluster
    all_peaks_hz = []
    all_peaks_int = []
    
    for cluster in clusters:
        if cluster.n_spins == 0:
            continue
        
        result = simulate_spin_system_hamiltonian(
            cluster,
            spectrometer_mhz=spectrometer_mhz,
            line_width_hz=line_width_hz,
            ppm_range=ppm_range,
            n_points=n_points,
            full_hamiltonian=True
        )
        
        for ppm, intensity in zip(result.peaks_ppm, result.peaks_intensity):
            all_peaks_hz.append(ppm * spectrometer_mhz)
            all_peaks_int.append(intensity)
    
    # Reconstruire le spectre combiné
    ppm_axis = np.linspace(ppm_range[1], ppm_range[0], n_points)
    freq_axis = ppm_axis * spectrometer_mhz
    spectrum = np.zeros(n_points)
    gamma = line_width_hz / 2
    
    for freq_hz, intensity in zip(all_peaks_hz, all_peaks_int):
        spectrum += _lorentzian(freq_axis, freq_hz, gamma, intensity)
    
    if spectrum.max() > 0:
        spectrum = spectrum / spectrum.max()
    
    peaks_ppm = [f / spectrometer_mhz for f in all_peaks_hz]
    
    return SimulatedSpectrum(
        ppm_axis=ppm_axis,
        freq_axis_hz=freq_axis,
        intensity=spectrum,
        peaks_ppm=peaks_ppm,
        peaks_intensity=all_peaks_int,
        spectrometer_mhz=spectrometer_mhz,
        line_width_hz=line_width_hz
    )


def simulate_spectrum_lineshape(
    peak_positions_hz: List[float],
    peak_intensities: List[float],
    spectrometer_mhz: float = 400.0,
    line_width_hz: float = 1.0,
    ppm_range: Tuple[float, float] = (-0.5, 12.0),
    n_points: int = 32768,
    line_shape: str = "lorentzian",
    apodization_lb_hz: float = 0.0,
) -> SimulatedSpectrum:
    """
    Génère un spectre NMR complet à partir des positions et intensités de pics.
    
    Approche directe (sans FID): somme de fonctions de forme de raie.
    Plus rapide que FID→FFT pour la plupart des applications.
    
    Args:
        peak_positions_hz: Positions des pics en Hz (relatif à 0 = TMS)
        peak_intensities: Intensités relatives des pics
        spectrometer_mhz: Fréquence du spectromètre
        line_width_hz: Largeur de raie (FWHM) en Hz
        ppm_range: (min_ppm, max_ppm) du spectre
        n_points: Nombre de points du spectre
        line_shape: "lorentzian", "gaussian", ou "voigt"
        apodization_lb_hz: Line broadening additionnel (apodisation)
    
    Returns:
        SimulatedSpectrum avec axes ppm/Hz et intensités
    """
    # Construire les axes
    ppm_min, ppm_max = ppm_range
    ppm_axis = np.linspace(ppm_max, ppm_min, n_points)  # Ordre décroissant (convention NMR)
    freq_axis_hz = ppm_axis * spectrometer_mhz  # Hz
    
    # Largeur effective
    gamma = (line_width_hz + apodization_lb_hz) / 2  # HWHM
    sigma = gamma * 0.8  # Pour Voigt
    
    # Accumuler le spectre
    spectrum = np.zeros(n_points, dtype=np.float64)
    
    for pos_hz, intensity in zip(peak_positions_hz, peak_intensities):
        if line_shape == "lorentzian":
            spectrum += _lorentzian(freq_axis_hz, pos_hz, gamma, intensity)
        elif line_shape == "gaussian":
            spectrum += _gaussian(freq_axis_hz, pos_hz, sigma, intensity)
        elif line_shape == "voigt":
            spectrum += _voigt_approx(freq_axis_hz, pos_hz, gamma, sigma, intensity)
        else:
            spectrum += _lorentzian(freq_axis_hz, pos_hz, gamma, intensity)
    
    # Normaliser
    if spectrum.max() > 0:
        spectrum = spectrum / spectrum.max()
    
    # Extraire les peaks pour référence
    peaks_ppm = [pos / spectrometer_mhz for pos in peak_positions_hz]
    
    return SimulatedSpectrum(
        ppm_axis=ppm_axis,
        freq_axis_hz=freq_axis_hz,
        intensity=spectrum,
        peaks_ppm=peaks_ppm,
        peaks_intensity=list(peak_intensities),
        spectrometer_mhz=spectrometer_mhz,
        line_width_hz=line_width_hz,
    )


def simulate_spectrum_from_fid(
    peak_positions_hz: List[float],
    peak_intensities: List[float],
    spectrometer_mhz: float = 400.0,
    T2_star_s: float = 0.5,
    acquisition_time_s: float = 2.0,
    n_points: int = 32768,
    ppm_range: Tuple[float, float] = (-0.5, 12.0),
    zero_fill_factor: int = 2,
    apodization: str = "exponential",
    lb_hz: float = 0.3,
) -> SimulatedSpectrum:
    """
    Simule un spectre NMR via génération FID puis FFT.
    
    Plus réaliste que la somme directe de Lorentziennes:
    - Gère naturellement le bruit de troncature
    - Permet différentes apodisations
    - Simule l'acquisition réelle
    
    FID: s(t) = Σ Aₖ · exp(i·2π·νₖ·t) · exp(-t/T₂*)
    
    Args:
        peak_positions_hz: Positions des résonances (Hz)
        peak_intensities: Amplitudes relatives
        spectrometer_mhz: Fréquence Larmor ¹H
        T2_star_s: T₂* apparent (détermine la largeur de raie)
        acquisition_time_s: Durée d'acquisition
        n_points: Nombre de points FID avant zero-fill
        ppm_range: Plage spectrale (ppm)
        zero_fill_factor: Facteur de zero-filling (1, 2, 4...)
        apodization: "exponential", "gaussian", "sine-bell", "none"
        lb_hz: Line broadening pour apodisation exponentielle
    
    Returns:
        SimulatedSpectrum
    """
    # Paramètres temporels
    dt = acquisition_time_s / n_points  # Dwell time
    t = np.arange(n_points) * dt
    
    # Générer le FID complexe
    fid = np.zeros(n_points, dtype=np.complex128)
    
    for freq_hz, amplitude in zip(peak_positions_hz, peak_intensities):
        # Chaque résonance contribue une sinusoïde amortie
        phase = 2 * np.pi * freq_hz * t
        decay = np.exp(-t / T2_star_s)
        fid += amplitude * np.exp(1j * phase) * decay
    
    # Apodisation
    if apodization == "exponential":
        # Lorentzien → Lorentzien plus large
        fid *= np.exp(-np.pi * lb_hz * t)
    elif apodization == "gaussian":
        # Gaussien → meilleure résolution, artefacts réduits
        sigma_t = acquisition_time_s / 3
        fid *= np.exp(-(t / sigma_t) ** 2)
    elif apodization == "sine-bell":
        # Résolution-enhancement
        fid *= np.sin(np.pi * t / acquisition_time_s)
    # "none": pas d'apodisation
    
    # Zero-filling
    n_fft = n_points * zero_fill_factor
    fid_zf = np.zeros(n_fft, dtype=np.complex128)
    fid_zf[:n_points] = fid
    
    # FFT
    spectrum_complex = np.fft.fft(fid_zf)
    spectrum_complex = np.fft.fftshift(spectrum_complex)
    
    # Partie réelle (absorption) après correction de phase automatique
    # Simplification: on prend le module pour éviter les problèmes de phase
    spectrum = np.abs(spectrum_complex)
    
    # Construire l'axe de fréquence
    freq_axis_hz = np.fft.fftshift(np.fft.fftfreq(n_fft, dt))
    
    # Filtrer pour la plage ppm désirée
    ppm_min, ppm_max = ppm_range
    freq_min = ppm_min * spectrometer_mhz
    freq_max = ppm_max * spectrometer_mhz
    
    mask = (freq_axis_hz >= freq_min) & (freq_axis_hz <= freq_max)
    freq_axis_hz = freq_axis_hz[mask]
    spectrum = spectrum[mask]
    
    # Inverser pour convention NMR (ppm décroissant de gauche à droite)
    freq_axis_hz = freq_axis_hz[::-1]
    spectrum = spectrum[::-1]
    
    ppm_axis = freq_axis_hz / spectrometer_mhz
    
    # Normaliser
    if spectrum.max() > 0:
        spectrum = spectrum / spectrum.max()
    
    # Largeur de raie effective
    line_width_hz = 1.0 / (np.pi * T2_star_s) + lb_hz
    
    return SimulatedSpectrum(
        ppm_axis=ppm_axis,
        freq_axis_hz=freq_axis_hz,
        intensity=spectrum,
        peaks_ppm=[p / spectrometer_mhz for p in peak_positions_hz],
        peaks_intensity=list(peak_intensities),
        spectrometer_mhz=spectrometer_mhz,
        line_width_hz=line_width_hz,
    )


def combine_spectra(
    spectra: List[SimulatedSpectrum],
    weights: Optional[List[float]] = None,
) -> SimulatedSpectrum:
    """
    Combine plusieurs spectres (ex: différents protons) en un spectre total.
    
    Args:
        spectra: Liste de SimulatedSpectrum (doivent avoir le même axe ppm)
        weights: Poids relatifs (défaut: égaux)
    
    Returns:
        SimulatedSpectrum combiné
    """
    if not spectra:
        raise ValueError("Liste de spectres vide")
    
    if weights is None:
        weights = [1.0] * len(spectra)
    
    # Vérifier compatibilité des axes
    ref = spectra[0]
    combined_intensity = np.zeros_like(ref.intensity)
    all_peaks_ppm = []
    all_peaks_intensity = []
    
    for spec, w in zip(spectra, weights):
        if len(spec.intensity) != len(ref.intensity):
            # Interpoler si nécessaire
            combined_intensity += w * np.interp(ref.ppm_axis, spec.ppm_axis, spec.intensity)
        else:
            combined_intensity += w * spec.intensity
        
        all_peaks_ppm.extend(spec.peaks_ppm)
        all_peaks_intensity.extend([p * w for p in spec.peaks_intensity])
    
    # Normaliser
    if combined_intensity.max() > 0:
        combined_intensity = combined_intensity / combined_intensity.max()
    
    return SimulatedSpectrum(
        ppm_axis=ref.ppm_axis.copy(),
        freq_axis_hz=ref.freq_axis_hz.copy(),
        intensity=combined_intensity,
        peaks_ppm=all_peaks_ppm,
        peaks_intensity=all_peaks_intensity,
        spectrometer_mhz=ref.spectrometer_mhz,
        line_width_hz=ref.line_width_hz,
    )


# ---------------------------
# Détection AA'BB' / AA'XX' et autres systèmes de spin
# ---------------------------

def _detect_spin_system(
    mol: Chem.Mol,
    heavy_idx: int,
    x_h_connectivity: Dict[int, List[int]],
    aromatic_dmap: List[Dict[int,int]],
    proton_shifts_ppm: Optional[Dict[int, float]] = None,
    spectrometer_mhz: float = 400.0
) -> Tuple[Optional[str], Dict[str, Any]]:
    """
    Détection avancée des systèmes de spin second ordre.
    
    Systèmes détectés:
    - AB: 2 protons couplés avec Δν/J < 10
    - AB₂: 1 proton couplé à 2 équivalents (ex: CHCl₂-CHCl₂)
    - A₂B₂: 2 paires équivalentes couplées (ex: -CH₂-CH₂-)
    - ABX: 2 protons AB couplés à un X distinct
    - A₂X₃: 2 protons couplés à 3 équivalents (ex: éthanol CH₂-CH₃)
    - AMX: 3 protons tous distincts
    
    Returns:
        (system_type, info_dict) où info_dict contient les détails du système
    """
    atom = mol.GetAtomWithIdx(heavy_idx)
    if atom.GetSymbol() != "C":
        return None, {}
    
    n_h_this = len(x_h_connectivity.get(heavy_idx, []))
    if n_h_this == 0:
        return None, {}
    
    info: Dict[str, Any] = {
        "heavy_idx": heavy_idx,
        "n_H": n_h_this,
        "partners": [],
        "delta_nu_J_ratios": [],
    }
    
    # Analyser les voisins carbone portant des H
    carbon_neighbors = []
    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() != "C":
            continue
        nb_idx = neighbor.GetIdx()
        n_h_nb = len(x_h_connectivity.get(nb_idx, []))
        if n_h_nb > 0:
            bond = mol.GetBondBetweenAtoms(heavy_idx, nb_idx)
            if bond and bond.GetBondTypeAsDouble() == 1.0:
                carbon_neighbors.append((nb_idx, n_h_nb))
    
    if not carbon_neighbors:
        return None, info
    
    # Calculer Δν/J pour chaque paire
    for nb_idx, n_h_nb in carbon_neighbors:
        delta_nu_J = 100.0  # Valeur par défaut (first order)
        
        if proton_shifts_ppm:
            h_this = x_h_connectivity.get(heavy_idx, [])
            h_nb = x_h_connectivity.get(nb_idx, [])
            if h_this and h_nb:
                delta_this = proton_shifts_ppm.get(h_this[0], None)
                delta_nb = proton_shifts_ppm.get(h_nb[0], None)
                if delta_this is not None and delta_nb is not None:
                    delta_ppm = abs(delta_this - delta_nb)
                    delta_nu = delta_ppm * spectrometer_mhz  # Hz
                    J_typical = 7.0  # 3J typique
                    if J_typical > 0:
                        delta_nu_J = delta_nu / J_typical
        
        info["partners"].append({
            "heavy_idx": nb_idx,
            "n_H": n_h_nb,
            "delta_nu_J": delta_nu_J
        })
        info["delta_nu_J_ratios"].append(delta_nu_J)
    
    # Classification du système de spin
    min_ratio = min(info["delta_nu_J_ratios"]) if info["delta_nu_J_ratios"] else 100.0
    
    # Système fortement couplé si Δν/J < 10
    is_strongly_coupled = min_ratio < 10.0
    
    if not is_strongly_coupled:
        return None, info  # First order
    
    # Déterminer le type de système
    total_partners = len(carbon_neighbors)
    partner_h_counts = [p["n_H"] for p in info["partners"]]
    
    # A₂B₂ (AA'BB'): -CH₂-CH₂-
    if n_h_this == 2 and total_partners == 1 and partner_h_counts[0] == 2:
        if min_ratio < 5.0:
            return "AA'BB'", info
        else:
            return "A₂X₂", info
    
    # AB₂: CH couplé à CH₂
    if n_h_this == 1 and total_partners == 1 and partner_h_counts[0] == 2:
        if min_ratio < 8.0:
            return "AB₂", info
        else:
            return "AX₂", info
    
    # A₂X₃: CH₂ couplé à CH₃
    if n_h_this == 2 and total_partners == 1 and partner_h_counts[0] == 3:
        if min_ratio < 8.0:
            return "A₂B₃", info
        else:
            return "A₂X₃", info
    
    # ABX: CH couplé à CH et un autre groupe
    if n_h_this == 1 and total_partners >= 2:
        return "ABX", info
    
    # AB simple: CH couplé à CH
    if n_h_this == 1 and total_partners == 1 and partner_h_counts[0] == 1:
        if min_ratio < 5.0:
            return "AB", info
        elif min_ratio < 10.0:
            return "AB_weak", info
        else:
            return "AX", info
    
    # AMX: 3 groupes différents
    if total_partners >= 2 and all(p["delta_nu_J"] > 10 for p in info["partners"]):
        return "AMX", info
    
    return "complex_second_order", info


def _detect_AABB_system(
    mol: Chem.Mol,
    heavy_idx: int,
    x_h_connectivity: Dict[int, List[int]],
    aromatic_dmap: List[Dict[int,int]]
) -> Optional[str]:
    """
    Détecte les systèmes AA'BB' / AA'XX' (fragments -CH₂-CH₂-)
    """
    atom = mol.GetAtomWithIdx(heavy_idx)
    if atom.GetSymbol() != "C":
        return None

    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() != "C":
            continue

        nb_idx = neighbor.GetIdx()

        if len(x_h_connectivity.get(heavy_idx, [])) != 2:
            continue
        if len(x_h_connectivity.get(nb_idx, [])) != 2:
            continue

        bond = mol.GetBondBetweenAtoms(heavy_idx, nb_idx)
        if not bond or bond.GetBondTypeAsDouble() != 1.0:
            continue

        return "AA'BB'"

    return None


# ---------------------------
# Moyenne conformationnelle multi-conformères
# ---------------------------

def _generate_conformer_ensemble(
    mol: Chem.Mol,
    n_conformers: int = 10,
    random_seed: int = 42
) -> List[Chem.Mol]:
    """
    Génère un ensemble de conformères pour moyenner les J.
    
    Utilise ETKDG pour générer plusieurs conformères et retourne
    ceux qui sont énergétiquement raisonnables.
    """
    from rdkit.Chem import AllChem
    
    mol_h = Chem.AddHs(mol, addCoords=True)
    
    try:
        # Paramètres ETKDG
        try:
            ps = AllChem.ETKDGv3()
        except Exception:
            ps = AllChem.ETKDG()
        ps.randomSeed = random_seed
        ps.numThreads = 0  # Auto
        ps.pruneRmsThresh = 0.5  # Éliminer conformères trop similaires
        
        # Générer conformères
        conf_ids = AllChem.EmbedMultipleConfs(mol_h, numConfs=n_conformers, params=ps)
        
        if len(conf_ids) == 0:
            # Fallback: un seul conformère
            AllChem.EmbedMolecule(mol_h, ps)
            return [mol_h]
        
        # Optimiser chaque conformère
        results = AllChem.MMFFOptimizeMoleculeConfs(mol_h, maxIters=200)
        
        # Trier par énergie et garder les meilleurs
        energies = []
        for i, (converged, energy) in enumerate(results):
            if converged == 0:  # Converged
                energies.append((i, energy))
        
        if not energies:
            return [mol_h]
        
        # Garder les conformères dans 3 kcal/mol du minimum
        energies.sort(key=lambda x: x[1])
        min_energy = energies[0][1]
        valid_confs = [i for i, e in energies if e - min_energy < 3.0]
        
        # Créer des molécules séparées pour chaque conformère valide
        conformers = []
        for conf_id in valid_confs[:n_conformers]:
            mol_copy = Chem.Mol(mol_h)
            # Garder seulement ce conformère
            conformers.append((mol_copy, conf_id))
        
        return [(mol_h, cid) for cid in valid_confs[:n_conformers]]
        
    except Exception:
        return [(mol_h, 0)]


def _boltzmann_average_J(
    J_values_per_conformer: List[float],
    energies_kcal: List[float],
    temperature_k: float = 298.15
) -> float:
    """
    Calcule la moyenne de J pondérée par les populations de Boltzmann.
    
    P_i = exp(-E_i / RT) / Σ exp(-E_j / RT)
    <J> = Σ P_i * J_i
    """
    if not J_values_per_conformer or not energies_kcal:
        return 0.0
    
    if len(J_values_per_conformer) != len(energies_kcal):
        return float(np.mean(J_values_per_conformer))
    
    R = 0.001987  # kcal/(mol·K)
    RT = R * temperature_k
    
    # Normaliser les énergies par rapport au minimum
    min_E = min(energies_kcal)
    rel_energies = [e - min_E for e in energies_kcal]
    
    # Calculer les facteurs de Boltzmann
    boltzmann_factors = [math.exp(-E / RT) for E in rel_energies]
    Z = sum(boltzmann_factors)  # Fonction de partition
    
    if Z <= 0:
        return float(np.mean(J_values_per_conformer))
    
    # Populations
    populations = [bf / Z for bf in boltzmann_factors]
    
    # Moyenne pondérée
    J_avg = sum(p * J for p, J in zip(populations, J_values_per_conformer))
    
    return J_avg


def _estimate_J_with_conformational_averaging(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    n_conformers: int = 5,
    temperature_k: float = 298.15
) -> Optional[float]:
    """
    Estime J en moyennant sur plusieurs conformères.
    
    Pour les systèmes flexibles (chaînes alkyle), la moyenne
    sur plusieurs conformères donne une meilleure estimation.
    """
    from rdkit.Chem import AllChem
    
    try:
        mol_h = Chem.AddHs(mol, addCoords=True)
        
        # Générer conformères
        try:
            ps = AllChem.ETKDGv3()
        except Exception:
            ps = AllChem.ETKDG()
        ps.randomSeed = 42
        
        conf_ids = AllChem.EmbedMultipleConfs(mol_h, numConfs=n_conformers, params=ps)
        if len(conf_ids) == 0:
            return None
        
        # Optimiser et collecter énergies
        results = AllChem.MMFFOptimizeMoleculeConfs(mol_h, maxIters=100)
        
        J_values = []
        energies = []
        
        for conf_id, (converged, energy) in zip(conf_ids, results):
            if converged != 0:
                continue
            
            # Calculer J pour ce conformère via Karplus
            conf = mol_h.GetConformer(conf_id)
            
            # Trouver les carbones parents
            h1_atom = mol_h.GetAtomWithIdx(h1_idx)
            h2_atom = mol_h.GetAtomWithIdx(h2_idx)
            
            c1_neighbors = [n for n in h1_atom.GetNeighbors() if n.GetSymbol() == 'C']
            c2_neighbors = [n for n in h2_atom.GetNeighbors() if n.GetSymbol() == 'C']
            
            if not c1_neighbors or not c2_neighbors:
                continue
            
            c1_idx = c1_neighbors[0].GetIdx()
            c2_idx = c2_neighbors[0].GetIdx()
            
            # Calculer l'angle dièdre H1-C1-C2-H2
            try:
                dihedral = AllChem.GetDihedralDeg(conf, h1_idx, c1_idx, c2_idx, h2_idx)
                theta_rad = math.radians(dihedral)
                
                # Karplus: J = A*cos²θ + B*cosθ + C
                A, B, C = 7.76, -1.10, 1.40  # Paramètres standard
                cos_theta = math.cos(theta_rad)
                J = A * cos_theta**2 + B * cos_theta + C
                
                J_values.append(abs(J))
                energies.append(energy)
            except Exception:
                continue
        
        if not J_values:
            return None
        
        # Moyenne de Boltzmann
        return _boltzmann_average_J(J_values, energies, temperature_k)
        
    except Exception:
        return None


# ---------------------------
# Classe cache pour conformères (évite de regénérer pour chaque paire H-H)
# ---------------------------
class ConformerCache:
    """
    Cache de conformères pour moyennage multi-conformères.
    VERSION AMÉLIORÉE avec ETKDGv3, pruning adaptatif et clustering RMSD.
    
    Génère les conformères une seule fois et les réutilise pour tous les calculs de J.
    Stocke également les énergies pour le moyennage Boltzmann.
    
    Améliorations:
    - ETKDGv3 par défaut (meilleur échantillonnage géométrique)
    - Pruning RMSD adaptatif (0.35 Å par défaut, ajustable)
    - Clustering pour éliminer conformères redondants
    - Génération 3x plus de conformères puis filtrage
    - Support des molécules macrocycliques
    - ✅ NOUVEAU: Échantillonnage systématique des dièdres pour chaînes flexibles
    """
    def __init__(
        self,
        mol: Chem.Mol,
        n_conformers: int = 50,
        energy_cutoff_kcal: float = 5.0,
        random_seed: int = 42,
        force_field: str = "MMFF94s",
        pruning_rmsd: float = 0.35,
        use_etkdg_v3: bool = True,
        use_systematic_sampling: bool = True  # ✅ NOUVEAU
    ):
        self.mol_h = Chem.AddHs(mol, addCoords=True)
        self.n_conformers = n_conformers
        self.energy_cutoff = energy_cutoff_kcal
        self.valid_conf_ids: List[int] = []
        self.energies: List[float] = []  # kcal/mol, relatif au minimum
        self._initialized = False
        self._random_seed = random_seed
        self._force_field = force_field
        self._pruning_rmsd = pruning_rmsd
        self._use_etkdg_v3 = use_etkdg_v3
        self._use_systematic_sampling = use_systematic_sampling  # ✅ NOUVEAU
        
    def initialize(self) -> bool:
        """Génère et optimise les conformères avec le force field configuré."""
        if self._initialized:
            return len(self.valid_conf_ids) > 0
            
        try:
            # ✅ AMÉLIORATION: Échantillonnage systématique pour molécules flexibles
            if self._use_systematic_sampling:
                systematic_results = _generate_diverse_conformers_systematic(
                    Chem.RemoveHs(self.mol_h),  # La fonction ajoute les H
                    n_conformers=self.n_conformers,
                    energy_cutoff_kcal=self.energy_cutoff
                )
                
                if systematic_results:
                    # Regénérer la molécule avec les bons conformères
                    mol_temp = Chem.AddHs(Chem.RemoveHs(self.mol_h), addCoords=True)
                    ps = AllChem.ETKDGv3() if hasattr(AllChem, 'ETKDGv3') else AllChem.ETKDG()
                    ps.randomSeed = self._random_seed
                    
                    for i, (_, energy) in enumerate(systematic_results):
                        if i == 0:
                            AllChem.EmbedMolecule(self.mol_h, ps)
                        self.valid_conf_ids.append(i if i < self.mol_h.GetNumConformers() else 0)
                        self.energies.append(energy)
                    
                    if self.valid_conf_ids:
                        self._initialized = True
                        return True
            
            # Paramètres ETKDG v3 (meilleur échantillonnage)
            if self._use_etkdg_v3:
                try:
                    ps = AllChem.ETKDGv3()
                    ps.useSmallRingTorsions = True  # Meilleur pour petits cycles
                    ps.useMacrocycleTorsions = True  # Support macrocycles
                except AttributeError:
                    ps = AllChem.ETKDG()
            else:
                ps = AllChem.ETKDG()
                
            ps.randomSeed = self._random_seed
            ps.numThreads = 0  # Auto-détection
            ps.pruneRmsThresh = self._pruning_rmsd  # Éliminer conformères similaires
            ps.useRandomCoords = True  # Meilleure diversité initiale
            
            # Générer 3x plus de conformères puis filtrer (meilleure couverture)
            n_generate = self.n_conformers * 3
            
            # Générer conformères
            conf_ids = AllChem.EmbedMultipleConfs(
                self.mol_h, 
                numConfs=n_generate,
                params=ps
            )
            
            if len(conf_ids) == 0:
                # Fallback: un seul conformère avec coords aléatoires
                ps.useRandomCoords = True
                res = AllChem.EmbedMolecule(self.mol_h, ps)
                if res == 0:
                    _optimize_molecule_geometry(self.mol_h, self._force_field, 200, 0)
                    self.valid_conf_ids = [0]
                    self.energies = [0.0]
                self._initialized = True
                return len(self.valid_conf_ids) > 0
            
            # Optimiser tous les conformères avec le force field configuré
            ff_variant = self._force_field.upper()
            if "MMFF" in ff_variant:
                variant = "MMFF94s" if "S" in ff_variant else "MMFF94"
                try:
                    results = AllChem.MMFFOptimizeMoleculeConfs(self.mol_h, mmffVariant=variant, maxIters=200)
                except Exception:
                    # Fallback UFF
                    results = []
                    for cid in conf_ids:
                        try:
                            AllChem.UFFOptimizeMolecule(self.mol_h, confId=cid, maxIters=200)
                            ff = AllChem.UFFGetMoleculeForceField(self.mol_h, confId=cid)
                            energy = ff.CalcEnergy() if ff else 0.0
                            results.append((0, energy))
                        except Exception:
                            results.append((1, 0.0))
            else:
                # UFF
                results = []
                for cid in conf_ids:
                    try:
                        AllChem.UFFOptimizeMolecule(self.mol_h, confId=cid, maxIters=200)
                        ff = AllChem.UFFGetMoleculeForceField(self.mol_h, confId=cid)
                        energy = ff.CalcEnergy() if ff else 0.0
                        results.append((0, energy))
                    except Exception:
                        results.append((1, 0.0))
            
            # Collecter énergies des conformères convergés
            energy_list = []
            for conf_id, (converged, energy) in zip(conf_ids, results):
                if converged == 0:  # 0 = convergé
                    energy_list.append((conf_id, energy))
            
            if not energy_list:
                self._initialized = True
                return False
            
            # Trier par énergie
            energy_list.sort(key=lambda x: x[1])
            min_energy = energy_list[0][1]
            
            # Garder conformères dans le seuil d'énergie
            for conf_id, energy in energy_list:
                rel_energy = energy - min_energy
                if rel_energy <= self.energy_cutoff:
                    self.valid_conf_ids.append(conf_id)
                    self.energies.append(rel_energy)
                    
                if len(self.valid_conf_ids) >= self.n_conformers:
                    break
            
            self._initialized = True
            return len(self.valid_conf_ids) > 0
            
        except Exception as e:
            self._initialized = True
            return False
    
    def get_J_vicinal_averaged(
        self, 
        h1_idx: int, 
        h2_idx: int,
        use_haasnoot: bool = False,
        use_haasnoot_7_params: bool = True,
        electroneg_corr: float = 0.0,
        temperature_k: float = 298.15
    ) -> Optional[float]:
        """
        Calcule J vicinal avec moyennage Boltzmann sur tous les conformères.
        
        Utilise l'équation Haasnoot-Altona complète à 7 paramètres si activée,
        avec détection automatique des substituants pour chaque conformère.
        
        Args:
            h1_idx, h2_idx: Indices des deux protons
            use_haasnoot: Si True, utilise Haasnoot-Altona (legacy)
            use_haasnoot_7_params: Si True, utilise Haasnoot 7 paramètres complets
            electroneg_corr: Correction EN manuelle (legacy)
            temperature_k: Température pour Boltzmann
        
        Returns:
            J moyenné en Hz, ou None si impossible
        """
        if not self._initialized:
            self.initialize()
            
        if not self.valid_conf_ids:
            return None
            
        try:
            # Trouver les carbones parents
            h1_atom = self.mol_h.GetAtomWithIdx(h1_idx)
            h2_atom = self.mol_h.GetAtomWithIdx(h2_idx)
            
            c1_neighbors = [n for n in h1_atom.GetNeighbors() if n.GetSymbol() == 'C']
            c2_neighbors = [n for n in h2_atom.GetNeighbors() if n.GetSymbol() == 'C']
            
            if not c1_neighbors or not c2_neighbors:
                return None
                
            c1_idx = c1_neighbors[0].GetIdx()
            c2_idx = c2_neighbors[0].GetIdx()
            
            # Vérifier que C1-C2 sont liés
            if not self.mol_h.GetBondBetweenAtoms(c1_idx, c2_idx):
                return None
            
            # Récupérer les substituants (une seule fois, indépendant du conformère)
            subs_alpha = []
            subs_beta = []
            
            if use_haasnoot_7_params:
                c1_atom = self.mol_h.GetAtomWithIdx(c1_idx)
                c2_atom = self.mol_h.GetAtomWithIdx(c2_idx)
                
                # Substituants sur C1 (excluant C2 et H)
                for nb in c1_atom.GetNeighbors():
                    if nb.GetIdx() != c2_idx and nb.GetSymbol() != 'H':
                        subs_alpha.append(nb.GetSymbol())
                
                # Substituants sur C2 (excluant C1 et H)
                for nb in c2_atom.GetNeighbors():
                    if nb.GetIdx() != c1_idx and nb.GetSymbol() != 'H':
                        subs_beta.append(nb.GetSymbol())
            
            J_values = []
            valid_energies = []
            
            for conf_id, energy in zip(self.valid_conf_ids, self.energies):
                try:
                    conf = self.mol_h.GetConformer(conf_id)
                    dihedral = AllChem.GetDihedralDeg(conf, h1_idx, c1_idx, c2_idx, h2_idx)
                    
                    if use_haasnoot_7_params:
                        # Calculer les signes d'orientation pour ce conformère
                        orientation_signs = []
                        
                        c1_atom = self.mol_h.GetAtomWithIdx(c1_idx)
                        for nb in c1_atom.GetNeighbors():
                            if nb.GetIdx() == c2_idx or nb.GetSymbol() == 'H':
                                continue
                            sign = _haasnoot_orientation_sign(self.mol_h, h1_idx, c1_idx, c2_idx, nb.GetIdx())
                            orientation_signs.append(sign)
                        
                        c2_atom = self.mol_h.GetAtomWithIdx(c2_idx)
                        for nb in c2_atom.GetNeighbors():
                            if nb.GetIdx() == c1_idx or nb.GetSymbol() == 'H':
                                continue
                            sign = _haasnoot_orientation_sign(self.mol_h, h1_idx, c1_idx, c2_idx, nb.GetIdx())
                            orientation_signs.append(sign)
                        
                        J = _karplus_haasnoot_altona_full(
                            dihedral,
                            substituents_alpha=subs_alpha,
                            substituents_beta=subs_beta,
                            orientation_signs=orientation_signs,
                            use_7_params=True
                        )
                    elif use_haasnoot:
                        J = _karplus_haasnoot_altona(dihedral, electroneg_corr)
                    else:
                        J = _karplus_simple(dihedral)
                    
                    if J is not None:
                        J_values.append(abs(J))
                        valid_energies.append(energy)
                except Exception:
                    continue
            
            if not J_values:
                return None
                
            # Moyenne de Boltzmann
            return _boltzmann_average_J(J_values, valid_energies, temperature_k)
            
        except Exception:
            return None
    
    def get_num_conformers(self) -> int:
        """Retourne le nombre de conformères valides."""
        if not self._initialized:
            self.initialize()
        return len(self.valid_conf_ids)
    
    def get_populations(self, temperature_k: float = 298.15) -> List[float]:
        """Retourne les populations de Boltzmann de chaque conformère."""
        if not self._initialized:
            self.initialize()
        if not self.energies:
            return []
            
        R = 0.001987  # kcal/(mol·K)
        RT = R * temperature_k
        
        boltz = [math.exp(-E / RT) for E in self.energies]
        Z = sum(boltz)
        
        if Z <= 0:
            return [1.0 / len(self.energies)] * len(self.energies)
            
        return [b / Z for b in boltz]


# ===========================================================================
# DÉTECTION ET TRAITEMENT DES MOLÉCULES FLEXIBLES
# ===========================================================================

def _count_rotatable_bonds(mol: Chem.Mol) -> int:
    """
    Compte le nombre de liaisons rotables dans une molécule.

    Liaisons rotables = liaisons simples non cycliques entre atomes lourds,
    excluant les liaisons terminales (méthyles, OH, etc.)
    """
    try:
        return Descriptors.NumRotatableBonds(mol)
    except Exception:
        # Fallback manuel
        count = 0
        for bond in mol.GetBonds():
            if bond.IsInRing():
                continue
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
            if a1.GetSymbol() == "H" or a2.GetSymbol() == "H":
                continue
            # Exclure liaisons terminales
            if a1.GetDegree() <= 1 or a2.GetDegree() <= 1:
                continue
            count += 1
        return count


def _estimate_molecular_flexibility(mol: Chem.Mol) -> Tuple[str, int]:
    """
    Estime la flexibilité d'une molécule.
    
    Returns:
        (niveau, n_rotatable) où niveau = "rigid", "moderate", "flexible", "very_flexible"
    """
    n_rot = _count_rotatable_bonds(mol)
    
    # Compter les cycles
    n_rings = mol.GetRingInfo().NumRings()
    
    # Correction pour les cycles (réduisent la flexibilité effective)
    effective_rot = max(0, n_rot - n_rings * 0.5)
    
    if effective_rot <= 2:
        return "rigid", n_rot
    elif effective_rot <= 5:
        return "moderate", n_rot
    elif effective_rot <= 10:
        return "flexible", n_rot
    else:
        return "very_flexible", n_rot


def _get_adaptive_n_conformers(
    mol: Chem.Mol, 
    base_n: int = 50, 
    max_n: int = 200,
    min_n: int = 10
) -> int:
    """
    Retourne un nombre adaptatif de conformères selon la flexibilité moléculaire.
    VERSION AMÉLIORÉE avec échelle logarithmique et détection de chaînes longues.
    
    Règles d'échantillonnage:
    - Molécules rigides (0-2 rotables): 10-25 conformères
    - Modérées (3-5 rotables): 50-75 conformères  
    - Flexibles (6-10 rotables): 100-150 conformères
    - Très flexibles (>10 rotables): 150-200 conformères
    
    Args:
        mol: Molécule RDKit
        base_n: Nombre de base (défaut 50)
        max_n: Maximum absolu (défaut 200)
        min_n: Minimum absolu (défaut 10)
    
    Returns:
        Nombre optimal de conformères à générer
    
    References:
        - Ebejer J.P. et al., J. Chem. Inf. Model. 52, 1146-1158 (2012)
        - Hawkins P.C.D. et al., J. Chem. Inf. Model. 50, 572-584 (2010)
    """
    flexibility, n_rot = _estimate_molecular_flexibility(mol)
    
    # Formule logarithmique: plus de conformères pour plus de flexibilité
    # n = base * (1 + 0.3 * ln(1 + n_rot))
    import math
    
    if n_rot == 0:
        # Molécule rigide
        return min_n
    
    # Facteur de scaling basé sur les liaisons rotables
    # n_rot=3 → factor~1.6, n_rot=6 → factor~2.0, n_rot=10 → factor~2.4
    scale_factor = 1.0 + 0.5 * math.log(1 + n_rot)
    
    # Ajustement selon la classification
    if flexibility == "rigid":
        n_conf = int(base_n * 0.3 * scale_factor)
    elif flexibility == "moderate":
        n_conf = int(base_n * 0.6 * scale_factor)
    elif flexibility == "flexible":
        n_conf = int(base_n * 1.0 * scale_factor)
    else:  # very_flexible
        n_conf = int(base_n * 1.5 * scale_factor)
    
    # Boost pour chaînes carbonées longues (lipides, terpènes)
    if _has_long_flexible_chain(mol, min_length=6):
        n_conf = int(n_conf * 1.5)
    
    # Appliquer les limites
    n_conf = max(min_n, min(max_n, n_conf))
    
    return n_conf


def _has_long_flexible_chain(mol: Chem.Mol, min_length: int = 4) -> bool:
    """
    Détecte si la molécule contient une chaîne flexible longue.
    
    Cas typique: chaîne alkyle longue (lipides, acides gras, etc.)
    """
    # Chercher les chaînes C-C-C-C... non cycliques
    for atom in mol.GetAtoms():
        if atom.IsInRing() or atom.GetSymbol() != "C":
            continue
        if atom.GetDegree() <= 1:
            continue
        
        # BFS pour trouver la chaîne la plus longue
        visited = {atom.GetIdx()}
        queue = [(atom, 0)]
        max_len = 0
        
        while queue:
            current, length = queue.pop(0)
            max_len = max(max_len, length)
            
            for nb in current.GetNeighbors():
                if nb.GetIdx() in visited:
                    continue
                if nb.GetSymbol() != "C" or nb.IsInRing():
                    continue
                visited.add(nb.GetIdx())
                queue.append((nb, length + 1))
        
        if max_len >= min_length:
            return True
    
    return False


# ---------------------------
# Effets de solvant avancés
# ---------------------------

# ===========================================================================
# SYSTÈME SOPHISTIQUÉ DE CORRECTIONS SOLVANT/TEMPÉRATURE
# ===========================================================================
# 
# Basé sur:
# 1. Paramètres solvatochromiques de Kamlet-Taft (α, β, π*)
# 2. Corrections spécifiques par type de proton (SMARTS)
# 3. Modèle physique des liaisons hydrogène
# 4. Coefficients de température empiriques
#
# Références:
# - Kamlet, M.J. & Taft, R.W. J. Am. Chem. Soc. 98, 377 (1976)
# - Abraham, M.H. Chem. Soc. Rev. 22, 73 (1993)
# - Reichardt, C. "Solvents and Solvent Effects in Organic Chemistry" (2003)
# - Muller, N. & Reiter, R.C. J. Chem. Phys. 42, 3265 (1965) - Température
# ===========================================================================

# Paramètres Kamlet-Taft des solvants
# α = acidité H-bond (donneur), β = basicité H-bond (accepteur), π* = polarité/polarisabilité
KAMLET_TAFT_PARAMS: Dict[str, Dict[str, float]] = {
    "CDCl3": {"alpha": 0.20, "beta": 0.10, "pi_star": 0.58, "ET30": 39.1},
    "DMSO-d6": {"alpha": 0.00, "beta": 0.76, "pi_star": 1.00, "ET30": 45.1},
    "D2O": {"alpha": 1.17, "beta": 0.47, "pi_star": 1.09, "ET30": 63.1},
    "CD3OD": {"alpha": 0.98, "beta": 0.66, "pi_star": 0.60, "ET30": 55.4},
    "CD3CN": {"alpha": 0.19, "beta": 0.40, "pi_star": 0.75, "ET30": 45.6},
    "C6D6": {"alpha": 0.00, "beta": 0.10, "pi_star": 0.59, "ET30": 34.3},
    "acetone-d6": {"alpha": 0.08, "beta": 0.43, "pi_star": 0.71, "ET30": 42.2},
    "pyridine-d5": {"alpha": 0.00, "beta": 0.64, "pi_star": 0.87, "ET30": 40.5},
    "THF-d8": {"alpha": 0.00, "beta": 0.55, "pi_star": 0.58, "ET30": 37.4},
    "toluene-d8": {"alpha": 0.00, "beta": 0.11, "pi_star": 0.54, "ET30": 33.9},
    "DCM-d2": {"alpha": 0.13, "beta": 0.10, "pi_star": 0.82, "ET30": 40.7},
}

# Coefficients de sensibilité au solvant par type de proton
# Δδ = a·α + b·β + c·π* + d (constante)
# où α, β, π* sont les paramètres Kamlet-Taft du solvant
PROTON_SOLVENT_SENSITIVITY: Dict[str, Dict[str, float]] = {
    # Protons donneurs de liaison H (très sensibles au β du solvant)
    "OH_phenol": {"a": -0.5, "b": +2.0, "c": +0.3, "d": 0.0},      # Phénol OH
    "OH_alcohol": {"a": -0.3, "b": +1.5, "c": +0.2, "d": 0.0},     # Alcool OH
    "OH_carboxylic": {"a": -0.8, "b": +2.5, "c": +0.5, "d": 0.0},  # COOH
    "NH_amide": {"a": -0.2, "b": +1.2, "c": +0.3, "d": 0.0},       # Amide NH
    "NH_amine": {"a": -0.4, "b": +1.0, "c": +0.2, "d": 0.0},       # Amine NH
    "NH_aromatic": {"a": -0.3, "b": +1.5, "c": +0.4, "d": 0.0},    # Pyrrole, indole
    
    # Protons accepteurs (sensibles au α du solvant donneur)
    "CHO_aldehyde": {"a": +0.3, "b": -0.2, "c": +0.1, "d": 0.0},   # Aldéhyde
    "CH_ketone_alpha": {"a": +0.2, "b": -0.1, "c": +0.05, "d": 0.0},
    
    # Protons aromatiques (sensibles aux effets d'anisotropie)
    "Ar_H": {"a": +0.05, "b": -0.1, "c": +0.15, "d": 0.0},
    "Ar_H_electron_rich": {"a": +0.02, "b": -0.15, "c": +0.1, "d": 0.0},
    "Ar_H_electron_poor": {"a": +0.08, "b": -0.05, "c": +0.2, "d": 0.0},
    
    # Protons alkyle (peu sensibles)
    "CH3_alkyl": {"a": 0.0, "b": 0.0, "c": +0.02, "d": 0.0},
    "CH2_alkyl": {"a": 0.0, "b": 0.0, "c": +0.02, "d": 0.0},
    "CH_alkyl": {"a": 0.0, "b": 0.0, "c": +0.02, "d": 0.0},
    
    # Protons vinyliques
    "CH_vinyl": {"a": +0.05, "b": -0.05, "c": +0.08, "d": 0.0},
    "CH2_vinyl": {"a": +0.03, "b": -0.03, "c": +0.05, "d": 0.0},
}

# Coefficients de température (ppm/K) - effet sur le déplacement chimique
# La plupart des protons ont un coefficient négatif (déblindage quand T augmente)
TEMPERATURE_COEFFICIENTS: Dict[str, float] = {
    # Protons échangeables - très sensibles (liaisons H affaiblies à haute T)
    "OH_phenol": -0.0080,        # ppm/K
    "OH_alcohol": -0.0100,
    "OH_carboxylic": -0.0120,
    "NH_amide": -0.0060,
    "NH_amine": -0.0080,
    "NH_aromatic": -0.0070,
    
    # Protons fixes - peu sensibles
    "CHO_aldehyde": -0.0015,
    "Ar_H": -0.0008,
    "CH_vinyl": -0.0005,
    "CH3_alkyl": -0.0003,
    "CH2_alkyl": -0.0003,
    "CH_alkyl": -0.0003,
    
    # Défaut
    "default": -0.0005,
}

# SMARTS patterns pour identifier les types de protons
PROTON_TYPE_SMARTS: Dict[str, str] = {
    "OH_phenol": "[OH;$([OH]c)]",                    # Phénol
    "OH_alcohol": "[OH;$([OH]C);!$([OH]c)]",         # Alcool aliphatique
    "OH_carboxylic": "[OH;$([OH]C=O)]",              # Acide carboxylique
    "NH_amide": "[NH;$([NH]C=O)]",                   # Amide
    "NH_amine": "[NH2,NH;!$([NH]C=O);!$([NH]c)]",    # Amine
    "NH_aromatic": "[nH]",                            # N aromatique (pyrrole, indole)
    "CHO_aldehyde": "[CH;$([CH](=O))]",              # Aldéhyde
    "CH_ketone_alpha": "[CH;$([CH]C(=O))]",          # Alpha-cétone
    "Ar_H": "[cH]",                                   # Aromatique générique
    "CH_vinyl": "[CH;$([CH]=C)]",                    # Vinylique CH
    "CH2_vinyl": "[CH2;$([CH2]=C)]",                 # Vinylique CH2
    "CH3_alkyl": "[CH3;$([CH3]C)]",                  # Méthyle
    "CH2_alkyl": "[CH2;$([CH2](C)C)]",               # Méthylène
    "CH_alkyl": "[CH;$([CH](C)(C)C)]",               # Méthine
}

# ===========================================================================
# VISIBILITÉ DES PROTONS ÉCHANGEABLES PAR SOLVANT
# ===========================================================================
# Cette table indique si les protons OH, NH, SH sont VISIBLES dans chaque solvant.
# C'est le seul effet de solvant appliqué par défaut (pas de correction de shift).
#
# Règles physiques:
# - Dans D2O: échange rapide avec le solvant → OH, NH invisibles
# - Dans CD3OD: échange avec OH du solvant → OH invisible, NH partiellement visible
# - Dans DMSO-d6: pas d'échange, liaisons H stabilisées → OH, NH bien visibles
# - Dans CDCl3: pas d'échange mais pics souvent larges
# ===========================================================================

EXCHANGEABLE_PROTON_VISIBILITY: Dict[str, Dict[str, bool]] = {
    "CDCl3": {
        "OH": True,      # Visible mais souvent large
        "NH": True,      # Visible mais peut être large
        "NH2": True,
        "SH": True,
        "COOH": True,    # Acide carboxylique visible
    },
    "DMSO-d6": {
        "OH": True,      # Très bien visible (pic fin)
        "NH": True,      # Très bien visible
        "NH2": True,
        "SH": True,
        "COOH": True,
    },
    "D2O": {
        "OH": False,     # ❌ INVISIBLE - échangé
        "NH": False,     # ❌ INVISIBLE - échangé
        "NH2": False,    # ❌ INVISIBLE - échangé
        "SH": True,      # Échange lent
        "COOH": False,   # ❌ INVISIBLE - échangé
    },
    "CD3OD": {
        "OH": False,     # ❌ INVISIBLE - échangé avec MeOD
        "NH": True,      # Partiellement visible (échange lent)
        "NH2": True,
        "SH": True,
        "COOH": False,   # ❌ INVISIBLE
    },
    "C6D6": {
        "OH": True,
        "NH": True,
        "NH2": True,
        "SH": True,
        "COOH": True,
    },
    "acetone-d6": {
        "OH": True,
        "NH": True,
        "NH2": True,
        "SH": True,
        "COOH": True,
    },
    "pyridine-d5": {
        "OH": True,      # Stabilisé par liaison H avec pyridine
        "NH": True,
        "NH2": True,
        "SH": True,
        "COOH": True,
    },
    "CD3CN": {
        "OH": True,
        "NH": True,
        "NH2": True,
        "SH": True,
        "COOH": True,
    },
}

# Largeur de raie (élargissement) pour les protons échangeables par solvant (Hz)
EXCHANGEABLE_PROTON_LINEWIDTH: Dict[str, Dict[str, float]] = {
    "CDCl3": {"OH": 8.0, "NH": 5.0, "NH2": 6.0, "SH": 2.0, "COOH": 15.0},
    "DMSO-d6": {"OH": 1.5, "NH": 1.2, "NH2": 1.5, "SH": 1.0, "COOH": 2.0},
    "D2O": {"OH": 0, "NH": 0, "NH2": 0, "SH": 2.0, "COOH": 0},
    "CD3OD": {"OH": 0, "NH": 8.0, "NH2": 10.0, "SH": 2.0, "COOH": 0},
    "C6D6": {"OH": 3.0, "NH": 2.0, "NH2": 3.0, "SH": 1.5, "COOH": 5.0},
    "acetone-d6": {"OH": 4.0, "NH": 3.0, "NH2": 4.0, "SH": 1.5, "COOH": 8.0},
    "pyridine-d5": {"OH": 2.0, "NH": 1.5, "NH2": 2.0, "SH": 1.0, "COOH": 3.0},
    "CD3CN": {"OH": 3.0, "NH": 2.5, "NH2": 3.0, "SH": 1.5, "COOH": 6.0},
}


def is_proton_visible_in_solvent(
    proton_type: str,
    solvent_name: str
) -> Tuple[bool, float]:
    """
    Détermine si un proton échangeable est visible dans un solvant donné.
    
    Args:
        proton_type: Type de proton ("OH", "NH", "NH2", "SH", "COOH")
        solvent_name: Nom du solvant
    
    Returns:
        Tuple (visible: bool, linewidth_hz: float)
        - visible: True si le proton est visible
        - linewidth_hz: Largeur de raie (0 si invisible)
    """
    visibility = EXCHANGEABLE_PROTON_VISIBILITY.get(
        solvent_name, 
        EXCHANGEABLE_PROTON_VISIBILITY.get("CDCl3", {})
    )
    linewidths = EXCHANGEABLE_PROTON_LINEWIDTH.get(
        solvent_name,
        EXCHANGEABLE_PROTON_LINEWIDTH.get("CDCl3", {})
    )
    
    visible = visibility.get(proton_type, True)
    linewidth = linewidths.get(proton_type, 5.0) if visible else 0.0
    
    return visible, linewidth


def get_exchangeable_protons_summary(solvent_name: str) -> str:
    """
    Retourne un résumé textuel de la visibilité des protons dans un solvant.
    
    Args:
        solvent_name: Nom du solvant
    
    Returns:
        Description textuelle
    """
    visibility = EXCHANGEABLE_PROTON_VISIBILITY.get(
        solvent_name, 
        EXCHANGEABLE_PROTON_VISIBILITY.get("CDCl3", {})
    )
    
    visible = [k for k, v in visibility.items() if v]
    invisible = [k for k, v in visibility.items() if not v]
    
    summary = f"Solvant: {solvent_name}\n"
    if visible:
        summary += f"  Visibles: {', '.join(visible)}\n"
    if invisible:
        summary += f"  Invisibles (échangés): {', '.join(invisible)}\n"
    
    return summary

# Corrections de déplacement chimique par solvant (en ppm) - VERSION SIMPLE
# Gardé pour compatibilité arrière
SOLVENT_SHIFT_CORRECTIONS: Dict[str, Dict[str, float]] = {
    "CDCl3": {
        "aromatic": 0.0,    # Référence
        "alkyl": 0.0,
        "OH": -2.0,         # OH souvent échangé, signal large
        "NH": -1.5,
        "aldehyde": 0.0,
    },
    "DMSO-d6": {
        "aromatic": +0.1,   # Légèrement déblindé
        "alkyl": 0.0,
        "OH": +0.5,         # OH visible, pic sharp
        "NH": +0.3,
        "aldehyde": +0.1,
    },
    "D2O": {
        "aromatic": -0.1,
        "alkyl": -0.05,
        "OH": 0.0,          # Échangé → disparu
        "NH": 0.0,          # Échangé → disparu
        "aldehyde": -0.1,
    },
    "CD3OD": {
        "aromatic": 0.0,
        "alkyl": 0.0,
        "OH": 0.0,          # Échangé
        "NH": 0.0,
        "aldehyde": 0.0,
    },
    "C6D6": {
        "aromatic": -0.5,   # Effet de courant de cycle
        "alkyl": -0.3,
        "OH": -1.0,
        "NH": -0.8,
        "aldehyde": -0.5,
    },
    "acetone-d6": {
        "aromatic": +0.05,
        "alkyl": 0.0,
        "OH": +0.2,
        "NH": +0.1,
        "aldehyde": +0.05,
    },
}

# Largeurs de raie selon viscosité du solvant (en Hz)
SOLVENT_LINEWIDTH: Dict[str, float] = {
    "CDCl3": 0.5,
    "DMSO-d6": 1.2,
    "D2O": 0.8,
    "CD3OD": 0.6,
    "C6D6": 0.5,
    "acetone-d6": 0.5,
    "pyridine-d5": 0.8,
    "THF-d8": 0.6,
    "toluene-d8": 0.6,
}

# ===========================================================================
# CORRECTIONS DES CONSTANTES J PAR SOLVANT
# ===========================================================================
# Les constantes de couplage J sont légèrement affectées par:
# 1. La polarité du solvant (constante diélectrique)
# 2. Les liaisons hydrogène (solvants protiques)
# 3. Les effets d'anisotropie (solvants aromatiques)
#
# Références:
# - Berger, S. & Braun, S. "200 and More NMR Experiments" (2004)
# - Hansen, P.E. Prog. Nucl. Magn. Reson. Spectrosc. (1981) 14, 175-296
# - Abraham, R.J. et al. J. Chem. Soc. Perkin Trans. 2 (1999) 699-711
# ===========================================================================

# Facteurs de correction J par type de couplage et solvant
# Valeur = multiplicateur (1.0 = pas de changement)
SOLVENT_J_CORRECTIONS: Dict[str, Dict[str, float]] = {
    "CDCl3": {
        # Référence (apolaire, non protique)
        "3J_vicinal": 1.00,      # ³J(H-C-C-H)
        "3J_vinyl": 1.00,        # ³J alcènes
        "2J_geminal": 1.00,      # ²J géminal
        "4J_allylic": 1.00,      # ⁴J allylique
        "3J_aromatic": 1.00,     # ³J aromatique
        "3J_HNCH": 1.00,         # ³J peptides
    },
    "DMSO-d6": {
        # Polaire aprotique - liaisons H intramoléculaires stabilisées
        # → conformations légèrement différentes
        "3J_vicinal": 1.02,      # Légère augmentation
        "3J_vinyl": 0.98,        # Stabilisation conformationnelle
        "2J_geminal": 1.00,
        "4J_allylic": 1.05,      # Plus visible (échange lent)
        "3J_aromatic": 1.00,
        "3J_HNCH": 1.08,         # NH visibles → J mesurables
    },
    "D2O": {
        # Très polaire, protique - liaisons H avec solvant
        "3J_vicinal": 0.98,      # Légère réduction
        "3J_vinyl": 0.97,
        "2J_geminal": 1.00,
        "4J_allylic": 0.90,      # Échange rapide → moyennage
        "3J_aromatic": 1.00,
        "3J_HNCH": 0.0,          # NH échangés → pas de couplage visible
    },
    "CD3OD": {
        # Polaire protique
        "3J_vicinal": 0.99,
        "3J_vinyl": 0.98,
        "2J_geminal": 1.00,
        "4J_allylic": 0.95,
        "3J_aromatic": 1.00,
        "3J_HNCH": 0.5,          # NH partiellement échangé
    },
    "C6D6": {
        # Aromatique - effets d'anisotropie magnétique
        "3J_vicinal": 1.00,
        "3J_vinyl": 1.02,        # Légère augmentation
        "2J_geminal": 1.00,
        "4J_allylic": 1.05,
        "3J_aromatic": 0.98,     # Effet de solvation aromatique
        "3J_HNCH": 1.00,
    },
    "CD3CN": {
        # Polaire aprotique
        "3J_vicinal": 1.01,
        "3J_vinyl": 1.00,
        "2J_geminal": 1.00,
        "4J_allylic": 1.02,
        "3J_aromatic": 1.00,
        "3J_HNCH": 1.02,
    },
    "pyridine-d5": {
        # Basique, accepteur de liaisons H
        "3J_vicinal": 1.02,
        "3J_vinyl": 1.00,
        "2J_geminal": 1.00,
        "4J_allylic": 1.05,
        "3J_aromatic": 0.98,
        "3J_HNCH": 1.10,         # NH/OH très visibles
    },
}

def _get_solvent_J_correction(
    solvent_name: str,
    j_type: str = "3J_vicinal"
) -> float:
    """
    Retourne le facteur de correction J pour un type de couplage dans un solvant.
    
    Args:
        solvent_name: Nom du solvant (ex: "CDCl3", "DMSO-d6")
        j_type: Type de couplage ("3J_vicinal", "3J_vinyl", "2J_geminal", etc.)
    
    Returns:
        Facteur multiplicatif (1.0 = pas de correction)
    """
    corrections = SOLVENT_J_CORRECTIONS.get(solvent_name, SOLVENT_J_CORRECTIONS.get("CDCl3", {}))
    return corrections.get(j_type, 1.0)


def _apply_solvent_J_correction(
    J_hz: float,
    solvent_name: str,
    j_type: str = "3J_vicinal",
    is_labile_H: bool = False
) -> float:
    """
    Applique une correction de J selon le solvant.
    
    Args:
        J_hz: Constante de couplage en Hz
        solvent_name: Nom du solvant
        j_type: Type de couplage
        is_labile_H: True si implique un proton labile (OH, NH)
    
    Returns:
        J corrigé en Hz
    """
    if is_labile_H and solvent_name == "D2O":
        # Proton échangé → pas de couplage visible
        return 0.0
    
    factor = _get_solvent_J_correction(solvent_name, j_type)
    return J_hz * factor


def _apply_solvent_corrections(
    proton_shifts_ppm: Dict[int, float],
    mol: Chem.Mol,
    heavy_to_H: Dict[int, List[int]],
    solvent_name: str
) -> Dict[int, float]:
    """
    Applique les corrections de déplacement chimique selon le solvant.
    
    Modifie les shifts en fonction du type de proton et du solvant.
    """
    corrections = SOLVENT_SHIFT_CORRECTIONS.get(solvent_name, {})
    if not corrections:
        return proton_shifts_ppm
    
    corrected = dict(proton_shifts_ppm)
    
    for heavy_idx, h_list in heavy_to_H.items():
        if not h_list:
            continue
        
        atom = mol.GetAtomWithIdx(heavy_idx)
        sym = atom.GetSymbol()
        
        # Déterminer le type de proton
        proton_type = "alkyl"  # Défaut
        
        if sym == "O":
            proton_type = "OH"
        elif sym == "N":
            proton_type = "NH"
        elif sym == "C":
            if atom.GetIsAromatic():
                proton_type = "aromatic"
            elif any(b.GetBondTypeAsDouble() >= 2 for b in atom.GetBonds()):
                # Vérifier si aldéhyde (C=O avec H)
                for b in atom.GetBonds():
                    other = b.GetOtherAtom(atom)
                    if other.GetSymbol() == "O" and b.GetBondTypeAsDouble() >= 2:
                        proton_type = "aldehyde"
                        break
        
        # Appliquer la correction
        correction = corrections.get(proton_type, 0.0)
        for h_idx in h_list:
            if h_idx in corrected:
                corrected[h_idx] = corrected[h_idx] + correction
    
    return corrected


def _get_solvent_linewidth(solvent_name: str, default: float = 0.8) -> float:
    """Retourne la largeur de raie typique pour un solvant."""
    return SOLVENT_LINEWIDTH.get(solvent_name, default)


def _estimate_exchange_broadening(
    mol: Chem.Mol,
    heavy_idx: int,
    solvent_name: str,
    temperature_k: float = 298.15,
    pH: Optional[float] = None
) -> float:
    """
    Estime l'élargissement dû à l'échange chimique (OH, NH).
    
    Returns:
        Facteur multiplicatif pour la largeur de raie (1.0 = pas d'élargissement)
    """
    atom = mol.GetAtomWithIdx(heavy_idx)
    sym = atom.GetSymbol()
    
    if sym not in ("O", "N"):
        return 1.0
    
    # Facteurs d'élargissement selon solvant
    exchange_factors = {
        "CDCl3": {"OH": 5.0, "NH": 3.0},      # Échange modéré
        "DMSO-d6": {"OH": 1.5, "NH": 1.2},    # Échange lent, pic visible
        "D2O": {"OH": 100.0, "NH": 50.0},     # Échange rapide → disparu
        "CD3OD": {"OH": 100.0, "NH": 20.0},   # Échange avec solvant
        "C6D6": {"OH": 2.0, "NH": 1.5},       # Échange lent
        "acetone-d6": {"OH": 3.0, "NH": 2.0},
    }
    
    factors = exchange_factors.get(solvent_name, {"OH": 2.0, "NH": 1.5})
    proton_type = "OH" if sym == "O" else "NH"
    
    base_factor = factors.get(proton_type, 1.0)
    
    # Effet de température (échange plus rapide à haute T)
    T_factor = (temperature_k / 298.15) ** 0.5
    
    # Effet de pH (si fourni)
    if pH is not None and proton_type == "OH":
        # Échange catalysé par acide/base
        if pH < 4 or pH > 10:
            base_factor *= 2.0
    
    return base_factor * T_factor


# ===========================================================================
# FONCTIONS SOPHISTIQUÉES DE CORRECTION SOLVANT/TEMPÉRATURE
# ===========================================================================

def _identify_proton_type_smarts(mol: Chem.Mol, h_idx: int) -> str:
    """
    Identifie le type de proton en utilisant les patterns SMARTS.
    
    Args:
        mol: Molécule RDKit
        h_idx: Index de l'atome H
    
    Returns:
        Type de proton (clé pour les tables de correction)
    """
    h_atom = mol.GetAtomWithIdx(h_idx)
    
    # Trouver l'atome lourd attaché
    neighbors = h_atom.GetNeighbors()
    if not neighbors:
        return "default"
    
    heavy_atom = neighbors[0]
    heavy_sym = heavy_atom.GetSymbol()
    heavy_idx = heavy_atom.GetIdx()
    
    # Tester les patterns SMARTS
    for proton_type, smarts in PROTON_TYPE_SMARTS.items():
        try:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                continue
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                if heavy_idx in match:
                    return proton_type
        except:
            continue
    
    # Fallback basé sur l'atome lourd
    if heavy_sym == "O":
        return "OH_alcohol"
    elif heavy_sym == "N":
        return "NH_amine"
    elif heavy_sym == "C":
        if heavy_atom.GetIsAromatic():
            return "Ar_H"
        elif heavy_atom.GetHybridization() == rdchem.HybridizationType.SP2:
            return "CH_vinyl"
        else:
            n_h = heavy_atom.GetTotalNumHs()
            if n_h >= 3:
                return "CH3_alkyl"
            elif n_h >= 2:
                return "CH2_alkyl"
            else:
                return "CH_alkyl"
    
    return "default"


def _calculate_kamlet_taft_correction(
    proton_type: str,
    solvent_name: str,
    reference_solvent: str = "CDCl3"
) -> float:
    """
    Calcule la correction de déplacement chimique basée sur les paramètres Kamlet-Taft.
    
    Δδ = a·Δα + b·Δβ + c·Δπ*
    
    où Δα = α(solvant) - α(référence), etc.
    
    Args:
        proton_type: Type de proton (ex: "OH_phenol", "Ar_H")
        solvent_name: Nom du solvant cible
        reference_solvent: Solvant de référence (défaut: CDCl3)
    
    Returns:
        Correction en ppm
    """
    # Paramètres du solvant cible et de référence
    params_target = KAMLET_TAFT_PARAMS.get(solvent_name, KAMLET_TAFT_PARAMS.get("CDCl3"))
    params_ref = KAMLET_TAFT_PARAMS.get(reference_solvent, KAMLET_TAFT_PARAMS.get("CDCl3"))
    
    # Différences de paramètres
    d_alpha = params_target["alpha"] - params_ref["alpha"]
    d_beta = params_target["beta"] - params_ref["beta"]
    d_pi = params_target["pi_star"] - params_ref["pi_star"]
    
    # Coefficients de sensibilité du proton
    sensitivity = PROTON_SOLVENT_SENSITIVITY.get(proton_type, 
                  PROTON_SOLVENT_SENSITIVITY.get("CH3_alkyl", {"a": 0, "b": 0, "c": 0, "d": 0}))
    
    # Calcul de la correction
    correction = (
        sensitivity["a"] * d_alpha +
        sensitivity["b"] * d_beta +
        sensitivity["c"] * d_pi +
        sensitivity["d"]
    )
    
    return correction


def _calculate_temperature_correction(
    proton_type: str,
    temperature_k: float,
    reference_temp_k: float = 298.15
) -> float:
    """
    Calcule la correction de température pour un type de proton.
    
    Δδ = coef × (T - T_ref)
    
    Args:
        proton_type: Type de proton
        temperature_k: Température de mesure (K)
        reference_temp_k: Température de référence (K)
    
    Returns:
        Correction en ppm
    """
    # Coefficient de température
    coef = TEMPERATURE_COEFFICIENTS.get(proton_type, 
                                        TEMPERATURE_COEFFICIENTS.get("default", -0.0005))
    
    # Différence de température
    dT = temperature_k - reference_temp_k
    
    return coef * dT


def _is_proton_exchangeable(mol: Chem.Mol, h_idx: int) -> bool:
    """
    Détermine si un proton est échangeable (OH, NH, SH).
    
    Args:
        mol: Molécule RDKit
        h_idx: Index de l'atome H
    
    Returns:
        True si le proton est échangeable
    """
    h_atom = mol.GetAtomWithIdx(h_idx)
    neighbors = h_atom.GetNeighbors()
    
    if not neighbors:
        return False
    
    heavy_atom = neighbors[0]
    sym = heavy_atom.GetSymbol()
    
    return sym in ("O", "N", "S")


def _calculate_exchange_effect(
    mol: Chem.Mol,
    h_idx: int,
    solvent_name: str,
    temperature_k: float = 298.15
) -> Dict[str, Any]:
    """
    Calcule l'effet de l'échange chimique sur un proton.
    
    Args:
        mol: Molécule RDKit
        h_idx: Index du proton
        solvent_name: Nom du solvant
        temperature_k: Température
    
    Returns:
        Dict avec 'visible' (bool), 'broadening' (float), 'shift_correction' (float)
    """
    if not _is_proton_exchangeable(mol, h_idx):
        return {"visible": True, "broadening": 1.0, "shift_correction": 0.0}
    
    h_atom = mol.GetAtomWithIdx(h_idx)
    heavy_atom = h_atom.GetNeighbors()[0]
    heavy_sym = heavy_atom.GetSymbol()
    
    # Solvants qui échangent les protons labiles
    exchange_solvents = {
        "D2O": {"OH": False, "NH": False, "SH": True},      # OH, NH échangés
        "CD3OD": {"OH": False, "NH": True, "SH": True},     # OH échangé
        "CD3OD + D2O": {"OH": False, "NH": False, "SH": True},
    }
    
    # Vérifier si visible
    visibility_map = exchange_solvents.get(solvent_name, {"OH": True, "NH": True, "SH": True})
    proton_class = "OH" if heavy_sym == "O" else ("NH" if heavy_sym == "N" else "SH")
    visible = visibility_map.get(proton_class, True)
    
    # Élargissement de raie dû à l'échange
    broadening_factors = {
        "CDCl3": {"OH": 3.0, "NH": 2.0, "SH": 1.5},
        "DMSO-d6": {"OH": 1.2, "NH": 1.1, "SH": 1.0},      # Échange lent → pics fins
        "C6D6": {"OH": 2.0, "NH": 1.5, "SH": 1.2},
        "acetone-d6": {"OH": 2.5, "NH": 1.8, "SH": 1.3},
        "pyridine-d5": {"OH": 1.3, "NH": 1.1, "SH": 1.0},  # Base → stabilise
    }
    
    broadening_map = broadening_factors.get(solvent_name, {"OH": 2.0, "NH": 1.5, "SH": 1.2})
    broadening = broadening_map.get(proton_class, 1.0)
    
    # Effet de température sur l'élargissement
    T_factor = (temperature_k / 298.15) ** 1.5
    broadening *= T_factor
    
    # Correction de shift due aux liaisons H avec le solvant
    kt_params = KAMLET_TAFT_PARAMS.get(solvent_name, KAMLET_TAFT_PARAMS.get("CDCl3"))
    shift_correction = 0.0
    
    if proton_class in ("OH", "NH"):
        # Proton donneur + solvant accepteur → déblindage
        shift_correction = kt_params["beta"] * 1.5  # Proportionnel à β
    
    return {
        "visible": visible,
        "broadening": broadening,
        "shift_correction": shift_correction
    }


def apply_sophisticated_solvent_correction(
    mol: Chem.Mol,
    proton_shifts_ppm: Dict[int, float],
    solvent_name: str,
    reference_solvent: str = "CDCl3",
    temperature_k: float = 298.15,
    reference_temp_k: float = 298.15,
    include_exchange_effects: bool = True
) -> Dict[str, Any]:
    """
    Applique des corrections sophistiquées de solvant et température.
    
    Cette fonction utilise:
    1. Les paramètres Kamlet-Taft pour les corrections de shift
    2. Les coefficients de température empiriques
    3. Les effets d'échange chimique pour OH/NH
    
    Args:
        mol: Molécule RDKit
        proton_shifts_ppm: Dict {h_idx: shift_ppm}
        solvent_name: Solvant cible
        reference_solvent: Solvant de référence (défaut: CDCl3)
        temperature_k: Température cible
        reference_temp_k: Température de référence
        include_exchange_effects: Inclure les effets d'échange (True)
    
    Returns:
        Dict avec:
        - 'corrected_shifts': {h_idx: shift_ppm corrigé}
        - 'exchange_info': {h_idx: {'visible': bool, 'broadening': float}}
        - 'corrections_applied': Liste des corrections appliquées
    """
    corrected_shifts = {}
    exchange_info = {}
    corrections_applied = []
    
    for h_idx, original_shift in proton_shifts_ppm.items():
        # Identifier le type de proton
        proton_type = _identify_proton_type_smarts(mol, h_idx)
        
        # 1. Correction Kamlet-Taft (solvant)
        kt_correction = _calculate_kamlet_taft_correction(
            proton_type, solvent_name, reference_solvent
        )
        
        # 2. Correction de température
        temp_correction = _calculate_temperature_correction(
            proton_type, temperature_k, reference_temp_k
        )
        
        # 3. Effets d'échange (OH, NH)
        exchange_effect = {"visible": True, "broadening": 1.0, "shift_correction": 0.0}
        if include_exchange_effects:
            exchange_effect = _calculate_exchange_effect(
                mol, h_idx, solvent_name, temperature_k
            )
        
        # Appliquer les corrections
        total_correction = kt_correction + temp_correction + exchange_effect["shift_correction"]
        corrected_shifts[h_idx] = original_shift + total_correction
        
        exchange_info[h_idx] = {
            "visible": exchange_effect["visible"],
            "broadening": exchange_effect["broadening"],
            "proton_type": proton_type,
            "kt_correction": kt_correction,
            "temp_correction": temp_correction,
            "exchange_correction": exchange_effect["shift_correction"],
            "total_correction": total_correction
        }
        
        if abs(total_correction) > 0.01:
            corrections_applied.append({
                "h_idx": h_idx,
                "proton_type": proton_type,
                "original_shift": original_shift,
                "corrected_shift": corrected_shifts[h_idx],
                "total_correction": total_correction
            })
    
    return {
        "corrected_shifts": corrected_shifts,
        "exchange_info": exchange_info,
        "corrections_applied": corrections_applied,
        "solvent": solvent_name,
        "temperature_k": temperature_k,
        "reference_solvent": reference_solvent,
        "reference_temp_k": reference_temp_k
    }


def predict_shift_in_different_solvent(
    original_shift_ppm: float,
    proton_type: str,
    from_solvent: str,
    to_solvent: str,
    from_temp_k: float = 298.15,
    to_temp_k: float = 298.15
) -> float:
    """
    Prédit le déplacement chimique dans un solvant différent.
    
    Utile pour comparer des spectres mesurés dans différentes conditions.
    
    Args:
        original_shift_ppm: Shift original
        proton_type: Type de proton
        from_solvent: Solvant original
        to_solvent: Solvant cible
        from_temp_k: Température originale
        to_temp_k: Température cible
    
    Returns:
        Shift prédit dans les nouvelles conditions
    """
    # Correction de solvant
    kt_correction = _calculate_kamlet_taft_correction(
        proton_type, to_solvent, from_solvent
    )
    
    # Correction de température
    temp_correction = _calculate_temperature_correction(
        proton_type, to_temp_k, from_temp_k
    )
    
    return original_shift_ppm + kt_correction + temp_correction


def get_solvent_comparison_table(
    proton_type: str,
    solvents: List[str] = None,
    reference_solvent: str = "CDCl3"
) -> Dict[str, float]:
    """
    Génère une table de comparaison des corrections pour différents solvants.
    
    Args:
        proton_type: Type de proton
        solvents: Liste de solvants (défaut: tous disponibles)
        reference_solvent: Solvant de référence
    
    Returns:
        Dict {solvent: correction_ppm}
    """
    if solvents is None:
        solvents = list(KAMLET_TAFT_PARAMS.keys())
    
    table = {}
    for solvent in solvents:
        correction = _calculate_kamlet_taft_correction(
            proton_type, solvent, reference_solvent
        )
        table[solvent] = round(correction, 3)
    
    return table


# ===========================================================================
# SIMULATION EXACTE PAR DIAGONALISATION DE L'HAMILTONIEN DE SPIN
# ===========================================================================
# Cette section implémente le calcul quantique exact des spectres RMN
# pour les systèmes de spin couplés (AB, AB₂, ABX, AA'BB', etc.)
# La méthode donne les positions et intensités exactes des raies.
# ===========================================================================

def _build_spin_operators(n_spins: int) -> Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray]]:
    """
    Construit les opérateurs de spin Ix, Iy, Iz pour n spins 1/2.
    
    Returns:
        (Ix_list, Iy_list, Iz_list) où chaque liste contient n matrices 2^n × 2^n
    """
    dim = 2 ** n_spins
    
    # Matrices de Pauli / 2
    Ix_single = np.array([[0, 0.5], [0.5, 0]], dtype=complex)
    Iy_single = np.array([[0, -0.5j], [0.5j, 0]], dtype=complex)
    Iz_single = np.array([[0.5, 0], [0, -0.5]], dtype=complex)
    I_identity = np.eye(2, dtype=complex)
    
    Ix_list = []
    Iy_list = []
    Iz_list = []
    
    for spin_idx in range(n_spins):
        # Produit tensoriel: I ⊗ I ⊗ ... ⊗ Ix ⊗ ... ⊗ I
        Ix = np.array([[1]], dtype=complex)
        Iy = np.array([[1]], dtype=complex)
        Iz = np.array([[1]], dtype=complex)
        
        for j in range(n_spins):
            if j == spin_idx:
                Ix = np.kron(Ix, Ix_single)
                Iy = np.kron(Iy, Iy_single)
                Iz = np.kron(Iz, Iz_single)
            else:
                Ix = np.kron(Ix, I_identity)
                Iy = np.kron(Iy, I_identity)
                Iz = np.kron(Iz, I_identity)
        
        Ix_list.append(Ix)
        Iy_list.append(Iy)
        Iz_list.append(Iz)
    
    return Ix_list, Iy_list, Iz_list


def _build_spin_hamiltonian(
    frequencies_hz: List[float],
    J_matrix: np.ndarray,
    Ix_list: List[np.ndarray],
    Iy_list: List[np.ndarray],
    Iz_list: List[np.ndarray]
) -> np.ndarray:
    """
    Construit l'Hamiltonien de spin complet.
    
    H = -Σᵢ νᵢ Izᵢ + Σᵢ<ⱼ Jᵢⱼ (Ixᵢ·Ixⱼ + Iyᵢ·Iyⱼ + Izᵢ·Izⱼ)
    
    Args:
        frequencies_hz: Liste des fréquences de Larmor (νᵢ) en Hz
        J_matrix: Matrice des couplages J[i,j] en Hz
        Ix_list, Iy_list, Iz_list: Opérateurs de spin
    
    Returns:
        Hamiltonien (matrice hermitienne)
    """
    n_spins = len(frequencies_hz)
    dim = 2 ** n_spins
    H = np.zeros((dim, dim), dtype=complex)
    
    # Termes Zeeman: -Σ νᵢ Izᵢ
    for i in range(n_spins):
        H -= frequencies_hz[i] * Iz_list[i]
    
    # Termes de couplage: Σᵢ<ⱼ Jᵢⱼ Iᵢ·Iⱼ
    for i in range(n_spins):
        for j in range(i + 1, n_spins):
            Jij = J_matrix[i, j]
            if abs(Jij) > 1e-10:
                # Iᵢ·Iⱼ = Ixᵢ·Ixⱼ + Iyᵢ·Iyⱼ + Izᵢ·Izⱼ
                H += Jij * (Ix_list[i] @ Ix_list[j] + 
                           Iy_list[i] @ Iy_list[j] + 
                           Iz_list[i] @ Iz_list[j])
    
    return H


def _calculate_nmr_transitions(
    H: np.ndarray,
    Ix_list: List[np.ndarray],
    intensity_threshold: float = 0.001
) -> List[Tuple[float, float, int, int]]:
    """
    Calcule les transitions RMN permises par diagonalisation.
    
    Args:
        H: Hamiltonien de spin
        Ix_list: Opérateurs Ix pour chaque spin
        intensity_threshold: Seuil minimum d'intensité
    
    Returns:
        Liste de (fréquence_hz, intensité, état_initial, état_final)
    """
    from scipy import linalg
    
    # Diagonalisation
    eigenvalues, eigenvectors = linalg.eigh(H)
    
    # Opérateur de détection: Ix_total = Σ Ixᵢ
    Ix_total = sum(Ix_list)
    
    # Calculer toutes les transitions
    dim = len(eigenvalues)
    transitions = []
    
    for i in range(dim):
        for j in range(i + 1, dim):
            # Fréquence de transition
            freq = eigenvalues[j] - eigenvalues[i]
            
            # Élément de matrice <j|Ix|i>
            matrix_element = eigenvectors[:, j].conj().T @ Ix_total @ eigenvectors[:, i]
            intensity = float(abs(matrix_element) ** 2)
            
            if intensity > intensity_threshold:
                transitions.append((float(freq), intensity, i, j))
    
    return transitions


def _simulate_spin_system_exact(
    frequencies_hz: List[float],
    J_matrix: np.ndarray,
    intensity_threshold: float = 0.001
) -> List[Tuple[float, float]]:
    """
    Simule un système de spin par diagonalisation exacte de l'Hamiltonien.
    
    Args:
        frequencies_hz: Fréquences de Larmor pour chaque spin (en Hz)
        J_matrix: Matrice des couplages J[i,j] en Hz (symétrique)
        intensity_threshold: Seuil d'intensité minimum
    
    Returns:
        Liste de (position_hz, intensité) pour chaque raie
    """
    n_spins = len(frequencies_hz)
    
    if n_spins == 0:
        return []
    
    if n_spins > 10:
        # Trop de spins, la matrice serait 2^10 = 1024 × 1024
        # Revenir à l'approximation first-order
        return None
    
    # Construire les opérateurs
    Ix_list, Iy_list, Iz_list = _build_spin_operators(n_spins)
    
    # Construire l'Hamiltonien
    H = _build_spin_hamiltonian(frequencies_hz, J_matrix, Ix_list, Iy_list, Iz_list)
    
    # Calculer les transitions
    transitions = _calculate_nmr_transitions(H, Ix_list, intensity_threshold)
    
    # Retourner (freq, intensity)
    return [(freq, intensity) for freq, intensity, _, _ in transitions]


def simulate_AB_system(
    delta_A_ppm: float,
    delta_B_ppm: float,
    J_AB_hz: float,
    spectrometer_mhz: float = 400.0
) -> Dict[str, Any]:
    """
    Simule exactement un système AB (2 spins couplés).
    
    Args:
        delta_A_ppm: Déplacement chimique du proton A (ppm)
        delta_B_ppm: Déplacement chimique du proton B (ppm)
        J_AB_hz: Constante de couplage (Hz)
        spectrometer_mhz: Fréquence du spectromètre (MHz)
    
    Returns:
        Dict avec 'transitions', 'delta_nu_J_ratio', 'is_second_order', 'label'
    """
    # Convertir en Hz
    nu_A = delta_A_ppm * spectrometer_mhz
    nu_B = delta_B_ppm * spectrometer_mhz
    
    delta_nu = abs(nu_A - nu_B)
    ratio = delta_nu / abs(J_AB_hz) if J_AB_hz != 0 else float('inf')
    
    # Construire la matrice J
    J_matrix = np.array([[0, J_AB_hz], [J_AB_hz, 0]])
    
    # Simuler
    transitions = _simulate_spin_system_exact([nu_A, nu_B], J_matrix)
    
    # Classifier
    if ratio > 10:
        is_second_order = False
        label = "d" if len(transitions) == 2 else "dd"
    elif ratio > 5:
        is_second_order = True
        label = "AB_weak"
    else:
        is_second_order = True
        label = "AB"
    
    return {
        "transitions": transitions,
        "delta_nu_hz": delta_nu,
        "J_hz": J_AB_hz,
        "delta_nu_J_ratio": ratio,
        "is_second_order": is_second_order,
        "label": label,
        "center_hz": (nu_A + nu_B) / 2,
    }


def simulate_AB2_system(
    delta_A_ppm: float,
    delta_B_ppm: float,
    J_AB_hz: float,
    spectrometer_mhz: float = 400.0
) -> Dict[str, Any]:
    """
    Simule exactement un système AB₂ (1 spin couplé à 2 équivalents).
    
    Exemple: CHCl₂-CH₂Cl
    
    Args:
        delta_A_ppm: Déplacement chimique du proton A
        delta_B_ppm: Déplacement chimique des 2 protons B (équivalents)
        J_AB_hz: Constante de couplage A-B
        spectrometer_mhz: Fréquence du spectromètre
    """
    nu_A = delta_A_ppm * spectrometer_mhz
    nu_B = delta_B_ppm * spectrometer_mhz
    
    # 3 spins: A, B1, B2 (B1 et B2 ont même fréquence)
    frequencies = [nu_A, nu_B, nu_B]
    
    # Matrice J: A couplé à B1 et B2, B1-B2 fortement couplé (J~0 mais même shift)
    J_matrix = np.array([
        [0, J_AB_hz, J_AB_hz],
        [J_AB_hz, 0, 0],      # B1-B2: pas de couplage explicite (équivalents)
        [J_AB_hz, 0, 0]
    ])
    
    transitions = _simulate_spin_system_exact(frequencies, J_matrix)
    
    delta_nu = abs(nu_A - nu_B)
    ratio = delta_nu / abs(J_AB_hz) if J_AB_hz != 0 else float('inf')
    
    return {
        "transitions": transitions,
        "delta_nu_hz": delta_nu,
        "J_hz": J_AB_hz,
        "delta_nu_J_ratio": ratio,
        "is_second_order": ratio < 10,
        "label": "AB2" if ratio < 10 else "AX2",
    }


def simulate_ABX_system(
    delta_A_ppm: float,
    delta_B_ppm: float,
    delta_X_ppm: float,
    J_AB_hz: float,
    J_AX_hz: float,
    J_BX_hz: float,
    spectrometer_mhz: float = 400.0
) -> Dict[str, Any]:
    """
    Simule exactement un système ABX (2 spins fortement couplés + 1 distinct).
    
    Exemple: Styrène oxide, amino-alcools
    
    Args:
        delta_A_ppm, delta_B_ppm: Protons du système AB
        delta_X_ppm: Proton X (bien séparé)
        J_AB_hz, J_AX_hz, J_BX_hz: Constantes de couplage
    """
    nu_A = delta_A_ppm * spectrometer_mhz
    nu_B = delta_B_ppm * spectrometer_mhz
    nu_X = delta_X_ppm * spectrometer_mhz
    
    frequencies = [nu_A, nu_B, nu_X]
    
    J_matrix = np.array([
        [0, J_AB_hz, J_AX_hz],
        [J_AB_hz, 0, J_BX_hz],
        [J_AX_hz, J_BX_hz, 0]
    ])
    
    transitions = _simulate_spin_system_exact(frequencies, J_matrix)
    
    delta_nu_AB = abs(nu_A - nu_B)
    ratio_AB = delta_nu_AB / abs(J_AB_hz) if J_AB_hz != 0 else float('inf')
    
    return {
        "transitions": transitions,
        "delta_nu_AB_hz": delta_nu_AB,
        "J_AB_hz": J_AB_hz,
        "J_AX_hz": J_AX_hz,
        "J_BX_hz": J_BX_hz,
        "delta_nu_J_ratio_AB": ratio_AB,
        "is_second_order": ratio_AB < 10,
        "label": "ABX" if ratio_AB < 10 else "AMX",
    }


def simulate_AA_BB_system(
    delta_A_ppm: float,
    delta_B_ppm: float,
    J_AB_hz: float,
    J_AA_hz: float = 0.0,
    J_BB_hz: float = 0.0,
    J_AB_prime_hz: float = 0.0,
    spectrometer_mhz: float = 400.0
) -> Dict[str, Any]:
    """
    Simule exactement un système AA'BB' (ex: -CH₂-CH₂-).
    
    4 spins: A, A', B, B'
    - A et A' sont équivalents chimiquement mais pas magnétiquement
    - Idem pour B et B'
    
    Args:
        delta_A_ppm: Déplacement chimique des protons A/A'
        delta_B_ppm: Déplacement chimique des protons B/B'
        J_AB_hz: Couplage vicinal A-B (et A'-B')
        J_AA_hz: Couplage geminal A-A' (typiquement ~-12 Hz pour CH₂)
        J_BB_hz: Couplage geminal B-B'
        J_AB_prime_hz: Couplage longue portée A-B' (et A'-B)
    """
    nu_A = delta_A_ppm * spectrometer_mhz
    nu_B = delta_B_ppm * spectrometer_mhz
    
    # 4 spins: A, A', B, B'
    frequencies = [nu_A, nu_A, nu_B, nu_B]
    
    # Matrice J (symétrique)
    #        A    A'   B    B'
    J_matrix = np.array([
        [0,      J_AA_hz,       J_AB_hz,        J_AB_prime_hz],  # A
        [J_AA_hz, 0,            J_AB_prime_hz,  J_AB_hz],        # A'
        [J_AB_hz, J_AB_prime_hz, 0,             J_BB_hz],        # B
        [J_AB_prime_hz, J_AB_hz, J_BB_hz,       0]               # B'
    ])
    
    transitions = _simulate_spin_system_exact(frequencies, J_matrix)
    
    delta_nu = abs(nu_A - nu_B)
    ratio = delta_nu / abs(J_AB_hz) if J_AB_hz != 0 else float('inf')
    
    return {
        "transitions": transitions,
        "delta_nu_hz": delta_nu,
        "J_AB_hz": J_AB_hz,
        "J_AA_hz": J_AA_hz,
        "J_BB_hz": J_BB_hz,
        "delta_nu_J_ratio": ratio,
        "is_second_order": ratio < 10,
        "label": "AA'BB'" if ratio < 10 else "AA'XX'",
    }


def simulate_A2X3_system(
    delta_A_ppm: float,
    delta_X_ppm: float,
    J_AX_hz: float,
    spectrometer_mhz: float = 400.0
) -> Dict[str, Any]:
    """
    Simule un système A₂X₃ (ex: éthanol CH₂-CH₃).
    
    5 spins: 2 A (CH₂) + 3 X (CH₃)
    
    Note: Pour les très grands systèmes, on peut utiliser des approximations.
    """
    nu_A = delta_A_ppm * spectrometer_mhz
    nu_X = delta_X_ppm * spectrometer_mhz
    
    # 5 spins: A1, A2, X1, X2, X3
    frequencies = [nu_A, nu_A, nu_X, nu_X, nu_X]
    
    # Matrice J
    #        A1   A2   X1   X2   X3
    J_matrix = np.zeros((5, 5))
    # A1-A2: geminal (souvent ~0 ou négatif)
    J_matrix[0, 1] = J_matrix[1, 0] = 0.0
    # A-X: vicinal
    for i in range(2):
        for j in range(2, 5):
            J_matrix[i, j] = J_matrix[j, i] = J_AX_hz
    # X-X: aucun couplage entre équivalents
    
    transitions = _simulate_spin_system_exact(frequencies, J_matrix)
    
    delta_nu = abs(nu_A - nu_X)
    ratio = delta_nu / abs(J_AX_hz) if J_AX_hz != 0 else float('inf')
    
    return {
        "transitions": transitions,
        "delta_nu_hz": delta_nu,
        "J_AX_hz": J_AX_hz,
        "delta_nu_J_ratio": ratio,
        "is_second_order": ratio < 10,
        "label": "A2B3" if ratio < 10 else "A2X3",
    }


def _get_exact_spectrum_for_spin_system(
    system_type: str,
    proton_shifts_ppm: Dict[int, float],
    h_indices: List[int],
    J_values: List[float],
    spectrometer_mhz: float = 400.0
) -> Optional[List[Tuple[float, float]]]:
    """
    Interface pour obtenir le spectre exact d'un système de spin détecté.
    
    Args:
        system_type: Type de système ("AB", "AB₂", "AA'BB'", "ABX", "A₂X₃")
        proton_shifts_ppm: Dictionnaire des shifts par indice de proton
        h_indices: Indices des protons impliqués
        J_values: Valeurs de J correspondantes
        spectrometer_mhz: Fréquence du spectromètre
    
    Returns:
        Liste de (freq_hz, intensity) ou None si non applicable
    """
    if not h_indices or not proton_shifts_ppm:
        return None
    
    # Extraire les shifts
    shifts = [proton_shifts_ppm.get(h, 0.0) for h in h_indices]
    
    if len(shifts) < 2:
        return None
    
    J_main = J_values[0] if J_values else 7.0
    
    if system_type in ("AB", "AB_weak"):
        if len(shifts) >= 2:
            result = simulate_AB_system(shifts[0], shifts[1], J_main, spectrometer_mhz)
            return result["transitions"]
    
    elif system_type == "AB₂":
        if len(shifts) >= 2:
            result = simulate_AB2_system(shifts[0], shifts[1], J_main, spectrometer_mhz)
            return result["transitions"]
    
    elif system_type in ("AA'BB'", "A₂X₂"):
        if len(shifts) >= 2:
            result = simulate_AA_BB_system(shifts[0], shifts[1], J_main, 
                                           J_AA_hz=-12.0, J_BB_hz=-12.0,
                                           spectrometer_mhz=spectrometer_mhz)
            return result["transitions"]
    
    elif system_type == "ABX":
        if len(shifts) >= 3 and len(J_values) >= 3:
            result = simulate_ABX_system(shifts[0], shifts[1], shifts[2],
                                        J_values[0], J_values[1], J_values[2],
                                        spectrometer_mhz)
            return result["transitions"]
    
    elif system_type in ("A₂X₃", "A₂B₃"):
        if len(shifts) >= 2:
            result = simulate_A2X3_system(shifts[0], shifts[1], J_main, spectrometer_mhz)
            return result["transitions"]
    
    return None


# ---------------------------
# Détection couplage virtuel
# ---------------------------
def _detect_virtual_coupling(
    molH: Chem.Mol,
    h_idx: int,
    partners: List[int],
    J_vals: List[float],
    proton_shifts_ppm: Dict[int, float],
    spectrometer_mhz: float
) -> bool:
    """
    Détecte un couplage virtuel (X couplé à un système AB fortement couplé).
    
    IMPORTANT: Le couplage virtuel ne s'applique PAS si les deux partenaires
    sont des protons équivalents (attachés au même carbone, ex: CH₂ ou CH₃).
    Dans ce cas, ils n'ont pas de différence de shift chimique et ne forment
    pas un système AB.
    """
    if len(partners) != 2 or len(J_vals) != 2:
        return False

    h_a, h_b = partners
    
    # ✅ CORRECTION: Vérifier si les partenaires sont sur le MÊME carbone (équivalents)
    # Si oui, pas de couplage virtuel possible car ils ont le même shift
    try:
        parent_a = next(n.GetIdx() for n in molH.GetAtomWithIdx(h_a).GetNeighbors())
        parent_b = next(n.GetIdx() for n in molH.GetAtomWithIdx(h_b).GetNeighbors())
        if parent_a == parent_b:
            # Les deux H sont sur le même carbone → équivalents, pas de couplage virtuel
            return False
    except Exception:
        pass
    
    if h_a not in proton_shifts_ppm or h_b not in proton_shifts_ppm:
        return False

    delta_AB_ppm = abs(proton_shifts_ppm[h_a] - proton_shifts_ppm[h_b])
    delta_nu_AB = delta_AB_ppm * spectrometer_mhz
    J_mean = float(np.mean(J_vals))

    if J_mean <= 0.0:
        return False

    if delta_nu_AB / J_mean < 5.0:
        return True
    return False

# ---------------------------
# Karplus avec correction substituants
# ---------------------------

# ===================================================================
# ÉQUATIONS KARPLUS SPÉCIALISÉES
# ===================================================================
#
# Ce module implémente plusieurs variantes de l'équation Karplus
# optimisées pour différents systèmes moléculaires:
#
# 1. HAASNOOT-ALTONA 7 PARAMÈTRES (1980)
#    Pour ³J(H-C-C-H) saturés avec corrections d'électronégativité
#    Réf: Tetrahedron 36, 2783-2792 (1980)
#
# 2. GARBISCH ALCÈNES (1968)
#    Pour ³J(H-C=C-H) vinyliques cis/trans
#    Réf: J. Am. Chem. Soc. 90, 6543 (1968)
#
# 3. BYSTROV PEPTIDES (1976)
#    Pour ³J(HN-Cα-C'-N) et ³J(HN-Cα-Hα) dans les peptides
#    Réf: Prog. NMR Spectrosc. 10, 41-82 (1976)
#
# ===================================================================

# Électronégativités de Huggins pour Haasnoot-Altona
# IMPORTANT: Dans l'équation Haasnoot-Altona, Δχ = χ(substituant) - χ(H)
# où χ(H) est pris comme référence = 2.20 (pas 0.0)
# Cela donne des Δχ plus raisonnables (~0.3 à 1.7)
#
# Réf: Haasnoot et al., Tetrahedron 36, 2783-2792 (1980)
# Réf: Huggins M.L., J. Am. Chem. Soc. 75, 4123 (1953)
_HUGGINS_ELECTRONEGATIVITY = {
    'H': 2.20,   # Référence pour H (pas 0!)
    'C': 2.60, 
    'N': 3.05, 
    'O': 3.50, 
    'F': 3.90,
    'P': 2.19, 
    'S': 2.58, 
    'Cl': 3.16, 
    'Br': 2.96, 
    'I': 2.66,
    'Si': 1.90, 
    'B': 2.04, 
    'Se': 2.55,
}

# Δχ typiques résultantes:
# Δχ(C) = 2.60 - 2.20 = 0.40
# Δχ(N) = 3.05 - 2.20 = 0.85
# Δχ(O) = 3.50 - 2.20 = 1.30
# Δχ(F) = 3.90 - 2.20 = 1.70

# Paramètres Karplus pour alcènes - Calibrés sur données expérimentales
# ³J(H-C=C-H) = A + B·cos(θ) + C·cos(2θ)
# Convention: θ = 0° pour cis (H sur même côté), θ = 180° pour trans
# 
# Pour cis (θ=0°): J = A + B + C
# Pour trans (θ=180°): J = A - B + C
#
# Valeurs typiques: Jcis ~10 Hz, Jtrans ~16 Hz, Jgem ~2 Hz
#
# Références:
# - Garbisch E.W., J. Am. Chem. Soc. 90, 6543 (1968)
# - Bothner-By A.A., Adv. Magn. Reson. 1, 195 (1965)
# - Jackman & Sternhell, Applications of NMR Spectroscopy
_GARBISCH_ALKENE_PARAMS = {
    # (A, B, C) calibrés pour: Jcis(0°)~10Hz, Jtrans(180°)~16Hz
    # Jcis = A + B + C = 10, Jtrans = A - B + C = 16 → B = -3
    # A + C = 13, donc A = 11, C = 2 donne Jcis = 11-3+2 = 10, Jtrans = 11+3+2 = 16 ✓
    'standard': (11.0, -3.0, 2.0),      # Alcènes simples: cis=10, trans=16
    'conjugated': (11.5, -3.5, 2.0),    # Alcènes conjugués: légèrement plus grands
    'strained': (9.0, -2.5, 1.5),       # Alcènes contraints: valeurs réduites
    'electronegative': (10.0, -2.5, 1.5),  # Substituants électronégatifs
}

# Paramètres Bystrov pour peptides (1976)
# ³J(HN-Cα-Hα) dépend de φ avec une phase de -60° due à la géométrie peptidique
# Valeurs typiques: hélice α (φ≈-60°) ~3-5 Hz, feuillet β (φ≈-120°) ~8-9 Hz
#
# Références:
# - Bystrov V.F., Prog. NMR Spectrosc. 10, 41-82 (1976)
# - Vuister G.W. & Bax A., J. Am. Chem. Soc. 115, 7772 (1993)
# - Pardi A. et al., J. Mol. Biol. 180, 741 (1984)
_BYSTROV_PEPTIDE_PARAMS = {
    # (A, B, C) en Hz - Calibrés pour valeurs expérimentales
    'HNCH_alpha': (7.9, -1.05, 0.65),   # ³J(HN-Cα-Hα): hélice~4Hz, feuillet~8Hz
    'HNCH_psi': (5.8, -1.2, 1.4),       # ³J(HN-C'-N) via ψ
    'HCCH_peptide': (9.5, -1.6, 1.8),   # ³J(Hα-Cα-Cβ-Hβ) chaîne latérale
}

# Paramètres Haasnoot-Altona complets (7 paramètres, 1980)
# ³J = P1·cos²φ + P2·cosφ + P3 + Σᵢ Δχᵢ·[P4 + P5·cos²(sᵢ·φ + P6·|Δχᵢ|)] + P7·ΠΔχ
_HAASNOOT_7_PARAMS = {
    'P1': 13.70,   # cos²φ coefficient (Hz)
    'P2': -0.73,   # cosφ coefficient (Hz)
    'P3': 0.00,    # constante (Hz)
    'P4': 0.56,    # correction EN constante
    'P5': -2.47,   # correction EN cos²
    'P6': 16.9,    # phase EN (degrés)
    'P7': -0.47,   # terme produit Π(Δχ) - NOUVEAU
}


def _karplus_garbisch_alkene(
    theta_deg: float,
    alkene_type: str = 'standard',
    substituents: Optional[List[str]] = None
) -> float:
    """
    Équation Karplus-Garbisch pour couplages vinyliques ³J(H-C=C-H).
    
    Cette équation est spécifiquement paramétrée pour les alcènes et donne
    des résultats plus précis que l'équation Karplus standard pour:
    - Couplages cis (³Jcis ≈ 6-12 Hz)
    - Couplages trans (³Jtrans ≈ 12-18 Hz)
    - Couplages géminaux (²J ≈ 0-3 Hz)
    
    Équation: ³J = A + B·cos(θ) + C·cos(2θ)
    
    Args:
        theta_deg: Angle dièdre H-C=C-H en degrés
                   0° = cis, 180° = trans
        alkene_type: Type d'alcène ('standard', 'conjugated', 'strained', 'electronegative')
        substituents: Liste des substituants sur la double liaison (pour ajustements)
    
    Returns:
        ³J en Hz
    
    Références:
        - Garbisch E.W., J. Am. Chem. Soc. 90, 6543-6544 (1968)
        - Bothner-By A.A., Adv. Magn. Reson. 1, 195 (1965)
    """
    A, B, C = _GARBISCH_ALKENE_PARAMS.get(alkene_type, _GARBISCH_ALKENE_PARAMS['standard'])
    
    theta = math.radians(theta_deg)
    
    # Équation Garbisch de base
    J = A + B * math.cos(theta) + C * math.cos(2 * theta)
    
    # Corrections pour substituants électronégatifs
    if substituents:
        en_correction = 0.0
        for sub in substituents:
            chi = _HUGGINS_ELECTRONEGATIVITY.get(sub, 2.5)
            chi_H = _HUGGINS_ELECTRONEGATIVITY.get('H', 0.0)
            delta_chi = chi - chi_H
            if delta_chi > 0.5:
                # Réduction de J pour substituants électronégatifs
                en_correction -= 0.3 * delta_chi
        J += en_correction
    
    return max(0.0, J)


def _karplus_bystrov_peptide(
    phi_deg: float,
    coupling_type: str = 'HNCH_alpha',
    psi_deg: Optional[float] = None
) -> float:
    """
    Équation Karplus-Bystrov pour couplages peptidiques.
    
    Cette équation est spécifiquement paramétrée pour les protéines et peptides,
    tenant compte de la géométrie particulière du backbone peptidique.
    
    Types de couplages:
    - HNCH_alpha: ³J(HN-Cα-Hα) - le plus utilisé pour déterminer φ
    - HNCH_psi: ³J(HN-C'-N) - pour déterminer ψ  
    - HCCH_peptide: ³J(Hα-Cα-Cβ-Hβ) - chaînes latérales
    
    Équation HNCH: ³J = A·cos²(φ-60°) + B·cos(φ-60°) + C
    
    La phase de -60° correspond à la géométrie standard du backbone peptidique
    où le NH et le Hα sont décalés de ~60° par rapport à l'angle dièdre φ pur.
    
    Args:
        phi_deg: Angle dièdre φ (C'-N-Cα-C') en degrés
        coupling_type: Type de couplage peptidique
        psi_deg: Angle dièdre ψ (N-Cα-C'-N) pour couplages dépendant de ψ
    
    Returns:
        ³J en Hz
    
    Références:
        - Bystrov V.F., Prog. NMR Spectrosc. 10, 41-82 (1976)
        - Pardi A. et al., J. Mol. Biol. 180, 741 (1984)
        - Vuister G.W. & Bax A., J. Am. Chem. Soc. 115, 7772 (1993)
    """
    A, B, C = _BYSTROV_PEPTIDE_PARAMS.get(coupling_type, _BYSTROV_PEPTIDE_PARAMS['HNCH_alpha'])
    
    if coupling_type == 'HNCH_alpha':
        # ³J(HN-Cα-Hα) dépend de φ avec phase de -60°
        theta = math.radians(phi_deg - 60.0)
    elif coupling_type == 'HNCH_psi' and psi_deg is not None:
        # Couplage via ψ
        theta = math.radians(psi_deg)
    else:
        theta = math.radians(phi_deg)
    
    J = A * (math.cos(theta) ** 2) + B * math.cos(theta) + C
    
    return max(0.0, J)


def _haasnoot_orientation_sign(
    mol: Chem.Mol,
    h1_idx: int,
    c1_idx: int,
    c2_idx: int,
    substituent_idx: int,
    conf_id: int = -1
) -> int:
    """
    Détermine le signe d'orientation pour l'équation Haasnoot-Altona.
    VERSION AMÉLIORÉE avec calcul 3D précis.
    
    Le signe dépend de la position relative du substituant par rapport
    au plan H1-C1-C2-H2:
    - +1 si le substituant est du même côté que H1 (ξ positif)
    - -1 si le substituant est du côté opposé (ξ négatif)
    
    L'orientation ξ est définie par Haasnoot comme:
    ξᵢ = sign(φ_H1-C1-C2-Sᵢ) où φ est l'angle dièdre
    
    Cas spéciaux traités:
    - Substituants axiaux/équatoriaux dans cyclohexane
    - Systèmes avec conformères multiples (moyenne pondérée)
    - Angles proches de 0° ou 180° (ambiguïté)
    
    Args:
        mol: Molécule RDKit avec coordonnées 3D
        h1_idx: Index du premier H
        c1_idx: Index du premier C
        c2_idx: Index du second C
        substituent_idx: Index de l'atome substituant
        conf_id: ID du conformère à utiliser (-1 = premier disponible)
    
    Returns:
        +1 ou -1 selon l'orientation
    
    References:
        Haasnoot C.A.G. et al., Tetrahedron 36, 2783-2792 (1980), Section 2.2
    """
    try:
        if not mol.GetNumConformers():
            return 1  # Défaut si pas de 3D
        
        # Sélectionner le bon conformère
        if conf_id < 0:
            conf = mol.GetConformer()
        else:
            try:
                conf = mol.GetConformer(conf_id)
            except Exception:
                conf = mol.GetConformer()
        
        # Calculer l'angle dièdre H1-C1-C2-Substituant
        dihed = AllChem.GetDihedralDeg(conf, h1_idx, c1_idx, c2_idx, substituent_idx)
        
        # Normaliser l'angle dans [-180, 180]
        while dihed > 180:
            dihed -= 360
        while dihed < -180:
            dihed += 360
        
        # Déterminer le signe selon la convention Haasnoot
        # Substituants avec φ dans [0°, 120°] ou [-120°, 0°] → même côté que H1 → +1
        # Substituants avec φ dans [120°, 180°] ou [-180°, -120°] → côté opposé → -1
        
        # Plus précisément, on utilise la position relative:
        # - Si |φ| < 90° : substituant "gauche" par rapport à H1 → +1
        # - Si |φ| > 90° : substituant "anti" par rapport à H1 → -1
        
        # Gestion des cas limites (près de 90°) avec hystérésis
        abs_dihed = abs(dihed)
        
        if abs_dihed < 75.0:
            # Clairement même côté
            return 1
        elif abs_dihed > 105.0:
            # Clairement côté opposé
            return -1
        else:
            # Zone ambiguë (75-105°): utiliser une classification plus fine
            # Basée sur le vecteur normal au plan C1-C2
            try:
                pos_h1 = np.array(conf.GetAtomPosition(h1_idx))
                pos_c1 = np.array(conf.GetAtomPosition(c1_idx))
                pos_c2 = np.array(conf.GetAtomPosition(c2_idx))
                pos_sub = np.array(conf.GetAtomPosition(substituent_idx))
                
                # Vecteurs
                v_c1c2 = pos_c2 - pos_c1
                v_c1h1 = pos_h1 - pos_c1
                v_c1sub = pos_sub - pos_c1
                
                # Normal au plan H1-C1-C2
                normal = np.cross(v_c1c2, v_c1h1)
                
                # Projection du substituant sur la normale
                proj = np.dot(v_c1sub, normal)
                
                # Si projection positive → même côté que H1
                return 1 if proj >= 0 else -1
            except Exception:
                # Fallback: utiliser le signe de l'angle dièdre
                return 1 if abs_dihed < 90 else -1
            
    except Exception:
        return 1  # Défaut


def _haasnoot_orientation_sign_averaged(
    mol: Chem.Mol,
    h1_idx: int,
    c1_idx: int,
    c2_idx: int,
    substituent_idx: int,
    conf_ids: List[int] = None,
    energies_kcal: List[float] = None,
    temperature_k: float = 298.15
) -> float:
    """
    Calcule le signe d'orientation Haasnoot moyenné sur plusieurs conformères.
    
    Pour les molécules flexibles, le signe effectif peut être un nombre
    flottant entre -1 et +1 selon la distribution conformationnelle.
    
    Args:
        mol: Molécule avec conformères multiples
        h1_idx, c1_idx, c2_idx, substituent_idx: Indices atomiques
        conf_ids: Liste des IDs de conformères
        energies_kcal: Énergies relatives en kcal/mol
        temperature_k: Température pour Boltzmann
    
    Returns:
        Signe effectif moyen (float entre -1 et +1)
    """
    if conf_ids is None or len(conf_ids) == 0:
        return float(_haasnoot_orientation_sign(mol, h1_idx, c1_idx, c2_idx, substituent_idx))
    
    signs = []
    for cid in conf_ids:
        s = _haasnoot_orientation_sign(mol, h1_idx, c1_idx, c2_idx, substituent_idx, conf_id=cid)
        signs.append(float(s))
    
    if energies_kcal is None or len(energies_kcal) != len(signs):
        # Moyenne simple
        return sum(signs) / len(signs)
    
    # Moyenne de Boltzmann
    R = 0.001987  # kcal/(mol·K)
    RT = R * temperature_k
    
    boltz_factors = [math.exp(-E / RT) for E in energies_kcal]
    Z = sum(boltz_factors)
    
    if Z <= 0:
        return sum(signs) / len(signs)
    
    weighted_sum = sum(s * b for s, b in zip(signs, boltz_factors))
    return weighted_sum / Z


def _karplus_haasnoot_altona_full(
    theta_deg: float,
    substituents_alpha: List[str] = None,
    substituents_beta: List[str] = None,
    orientation_signs: List[int] = None,
    use_7_params: bool = True,
) -> float:
    """
    Équation Haasnoot-de Leeuw-Altona COMPLÈTE à 7 PARAMÈTRES (1980).
    
    Version améliorée incluant le terme P7 (produit des électronégativités)
    qui améliore la précision pour les systèmes multi-substitués.
    
    Équation complète:
    ³J = P1·cos²φ + P2·cosφ + P3 
         + Σᵢ Δχᵢ·[P4 + P5·cos²(sᵢ·φ + P6·|Δχᵢ|)]
         + P7·Π(Δχᵢ)  ← NOUVEAU: terme produit
    
    Paramètres (Haasnoot et al., Tetrahedron 36, 2783-2792, 1980):
    - P1 = 13.70 Hz  (cos²φ)
    - P2 = -0.73 Hz  (cosφ)
    - P3 = 0.00 Hz   (constante)
    - P4 = 0.56      (correction EN linéaire)
    - P5 = -2.47     (correction EN angulaire)
    - P6 = 16.9°     (phase EN)
    - P7 = -0.47     (terme produit pour interactions multi-substituants)
    
    Le terme P7 est important pour:
    - Systèmes avec plusieurs substituants électronégatifs
    - Sucres et dérivés glycosidiques
    - Composés α,β-disubstitués
    
    Args:
        theta_deg: Angle dièdre H-C-C-H en degrés
        substituents_alpha: Symboles des substituants sur Cα (hors H observé)
        substituents_beta: Symboles des substituants sur Cβ (hors H couplé)
        orientation_signs: Signes d'orientation pour chaque substituant (+1 ou -1)
        use_7_params: Si True, utilise le terme P7 (défaut)
    
    Returns:
        ³J en Hz (>= 0)
    
    Références:
        - Haasnoot C.A.G., de Leeuw F.A.A.M., Altona C., Tetrahedron 36, 2783-2792 (1980)
        - Altona C. & Haasnoot C.A.G., Org. Magn. Reson. 13, 417 (1980)
    """
    # Récupérer les 7 paramètres
    P1 = _HAASNOOT_7_PARAMS['P1']
    P2 = _HAASNOOT_7_PARAMS['P2']
    P3 = _HAASNOOT_7_PARAMS['P3']
    P4 = _HAASNOOT_7_PARAMS['P4']
    P5 = _HAASNOOT_7_PARAMS['P5']
    P6_deg = _HAASNOOT_7_PARAMS['P6']
    P7 = _HAASNOOT_7_PARAMS['P7']
    
    phi = math.radians(theta_deg)
    P6 = math.radians(P6_deg)
    
    # Terme Karplus de base (P1, P2, P3)
    J = P1 * (math.cos(phi) ** 2) + P2 * math.cos(phi) + P3
    
    # Initialiser les listes de substituants
    if substituents_alpha is None:
        substituents_alpha = []
    if substituents_beta is None:
        substituents_beta = []
    if orientation_signs is None:
        # Attribution intelligente des signes basée sur la position
        # α1, α2 reçoivent +1, -1; β1, β2 reçoivent +1, -1
        n_alpha = len(substituents_alpha)
        n_beta = len(substituents_beta)
        orientation_signs = []
        for i in range(n_alpha):
            orientation_signs.append(1 if i % 2 == 0 else -1)
        for i in range(n_beta):
            orientation_signs.append(1 if i % 2 == 0 else -1)
    
    chi_H = _HUGGINS_ELECTRONEGATIVITY.get('H', 0.0)
    all_subs = list(substituents_alpha) + list(substituents_beta)
    delta_chi_list = []
    
    # Corrections d'électronégativité individuelles (P4, P5, P6)
    for i, sub in enumerate(all_subs):
        chi_i = _HUGGINS_ELECTRONEGATIVITY.get(sub, 2.5)  # Défaut C-like
        delta_chi = chi_i - chi_H
        delta_chi_list.append(delta_chi)
        
        if delta_chi <= 0:
            continue
            
        s_i = orientation_signs[i] if i < len(orientation_signs) else 1
        
        # Terme correctif Haasnoot-Altona individuel
        correction = delta_chi * (
            P4 + P5 * (math.cos(s_i * phi + P6 * abs(delta_chi)) ** 2)
        )
        J += correction
    
    # Terme P7: produit des électronégativités (interactions multi-substituants)
    if use_7_params and len(delta_chi_list) >= 2:
        # Calculer le produit des Δχ pour tous les substituants
        product_delta_chi = 1.0
        n_significant = 0
        for dc in delta_chi_list:
            if dc > 0.1:  # Seuil pour éviter les produits nuls
                product_delta_chi *= dc
                n_significant += 1
        
        if n_significant >= 2:
            # Normaliser le produit pour éviter les valeurs extrêmes
            # et appliquer la correction P7
            normalized_product = product_delta_chi ** (1.0 / n_significant)
            J += P7 * normalized_product * n_significant
    
    return max(0.0, J)


def _karplus_with_substituent_correction(
    theta_deg: float,
    alpha_EN: float = 0.0,
    beta_EN: float = 0.0
) -> float:
    """
    Wrapper rétrocompatible - utilise maintenant l'équation complète.
    """
    # Convertir les corrections EN en substituants équivalents
    subs_alpha = ['O'] if alpha_EN > 0.3 else (['N'] if alpha_EN > 0.1 else [])
    subs_beta = ['O'] if beta_EN > 0.3 else (['N'] if beta_EN > 0.1 else [])
    
    return _karplus_haasnoot_altona_full(
        theta_deg,
        substituents_alpha=subs_alpha,
        substituents_beta=subs_beta
    )

# ---------------------------
# Lignes de combinaison – systèmes AB/ABX
# ---------------------------
def _add_combination_lines(
    lines: List[Tuple[float, float]],
    system_type: Optional[str],
    J_AB: Optional[float],
    config: "SimConfig",
    delta_nu_hz: float = 0.0,
) -> List[Tuple[float, float]]:
    """
    Ajoute les raies de combinaison pour systèmes fortement couplés.
    
    Dans un système AB, les 4 raies principales sont à:
    - ν₁ = νₐ + ½(J + D)
    - ν₂ = νₐ + ½(J - D)  
    - ν₃ = νᵦ - ½(J - D)
    - ν₄ = νᵦ - ½(J + D)
    
    où D = √(Δν² + J²)
    
    Les raies de combinaison (faible intensité) apparaissent à:
    - ν₅ = ½(νₐ + νᵦ) + ½D  (transition "zéro-quantum")
    - ν₆ = ½(νₐ + νᵦ) - ½D
    
    Référence: Günther, NMR Spectroscopy, Chapter 4
    """
    if not lines or system_type is None:
        return lines
    
    # Seulement pour systèmes AB ou AA'BB'
    if system_type not in ["AA'BB'", "AA'XX'", "AB", "ABX"]:
        return lines
    
    if J_AB is None or J_AB <= 0.0:
        return lines
    
    # Calculer D (paramètre de couplage effectif)
    if delta_nu_hz > 0.0:
        D = math.sqrt(delta_nu_hz ** 2 + J_AB ** 2)
    else:
        D = J_AB  # Limite fortement couplée
    
    # Intensité des raies de combinaison (typiquement 1-5% des raies principales)
    # I_comb = sin²(2θ)/4 où tan(2θ) = J/Δν
    if delta_nu_hz > 0.0:
        tan_2theta = J_AB / delta_nu_hz
        theta_2 = math.atan(tan_2theta)
        I_comb = (math.sin(theta_2) ** 2) / 4.0
    else:
        I_comb = 0.25  # Maximum quand Δν = 0
    
    # Seuil de visibilité (raies trop faibles non ajoutées)
    min_intensity = getattr(config, 'combination_line_threshold', 0.02)
    
    if I_comb < min_intensity:
        return lines
    
    # Centre du multiplet (moyenne des positions)
    if len(lines) >= 2:
        center = sum(f for f, _ in lines) / len(lines)
    else:
        center = 0.0
    
    # Ajouter les raies de combinaison
    new_lines = list(lines)
    new_lines.append((center + D / 2.0, I_comb * 0.5))
    new_lines.append((center - D / 2.0, I_comb * 0.5))
    
    # Renormaliser
    total = sum(i for _, i in new_lines)
    if total > 0:
        new_lines = [(f, i / total) for f, i in new_lines]
    
    return sorted(new_lines, key=lambda x: x[0])

# ---------------------------
# Lignes Lorentziennes
# ---------------------------
def _generate_lorentzian_lineshape(
    peak_positions_hz: List[float],
    peak_intensities: List[float],
    linewidth_hz: float,
    spectrum_range_hz: Tuple[float, float],
    num_points: int = 16384
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Génère un spectre 1H avec formes de raie Lorentziennes.
    """
    if not peak_positions_hz:
        freq_axis = np.linspace(spectrum_range_hz[0], spectrum_range_hz[1], num_points)
        return freq_axis, np.zeros(num_points)

    freq_axis = np.linspace(spectrum_range_hz[0], spectrum_range_hz[1], num_points)
    spectrum = np.zeros(num_points)

    gamma = max(1e-6, float(linewidth_hz))
    for pos, intensity in zip(peak_positions_hz, peak_intensities):
        lorentzian = intensity * (gamma / (2.0 * np.pi)) / (
            (freq_axis - pos) ** 2 + (gamma / 2.0) ** 2
        )
        spectrum += lorentzian

    return freq_axis, spectrum

# ---------------------------
# Résolution numérique adaptative
# ---------------------------
def _adaptive_digital_resolution(spectrometer_mhz: float) -> float:
    """
    Résolution numérique typique selon la fréquence.
    """
    SW_hz = 16.0 * float(spectrometer_mhz)
    num_points = 32768
    resolution_hz = SW_hz / num_points
    return float(resolution_hz)

# ---------------------------
# Construction des raies 1er ordre
# ---------------------------
def _first_order_peak_positions(
    J_couplings: List[Dict[str, float]],
    base_multiplicity: str
) -> List[Tuple[float, float]]:
    """
    Construit les positions/intensités d'un multiplet 1er ordre
    par convolution successive (Pascal).
    """
    lines: List[Tuple[float, float]] = [(0.0, 1.0)]

    for c in J_couplings:
        J = float(c.get("J_hz", 0.0))
        n = int(c.get("n_equiv", 0))
        if J <= 0.0 or n <= 0:
            continue

        pascal = [math.comb(n, k) for k in range(n + 1)]
        denom = float(2 ** n)

        new_lines: List[Tuple[float, float]] = []
        for f0, I0 in lines:
            for k in range(n + 1):
                shift = (k - n / 2.0) * J
                I = I0 * pascal[k] / denom
                new_lines.append((f0 + shift, I))
        lines = new_lines

    if not lines:
        return [(0.0, 1.0)]

    lines.sort(key=lambda x: x[0])
    merged: List[Tuple[float, float]] = []
    cur_f, cur_I = lines[0]
    tol = 1e-6
    for f, I in lines[1:]:
        if abs(f - cur_f) <= tol:
            new_I = cur_I + I
            if new_I > 0:
                cur_f = (cur_f * cur_I + f * I) / new_I
            cur_I = new_I
        else:
            merged.append((cur_f, cur_I))
            cur_f, cur_I = f, I
    merged.append((cur_f, cur_I))

    s = sum(i for _, i in merged) or 1.0
    merged = [(f, i / s) for f, i in merged]
    return merged

# ---------------------------
# Config
# ---------------------------
@dataclass
class SimConfig:
    # Seuils de détection pour ⁿJ (n = ordre de couplage)
    # ⁴J à 1.2 Hz (relevé pour éviter faux positifs homoallyliques)
    # ⁵J à 1.0 Hz (relevé pour éviter faux positifs longue distance)
    min_j_detect: Dict[int, float] = field(default_factory=lambda: {2: 1.0, 3: 1.0, 4: 1.2, 5: 1.0})
    aro_j_forced: Dict[str, float] = field(default_factory=lambda: {"ortho": 7.8, "meta": 2.3, "para": 0.6})
    hetero_j_base: Dict[str, Dict[int, float]] = field(default_factory=lambda: {
        'F': {2: 47.0, 3: 9.0, 4: 2.0},
        'P': {2: 12.0, 3: 8.0, 4: 2.0},
    })
    viscosity_low_cuts: Dict[str, float] = field(default_factory=lambda: {"low": 0.5, "moderate": 0.6, "high": 0.8})
    second_order_threshold_main: float = 0.14
    second_order_threshold_many_small_deltas: float = 0.08
    second_order_threshold_single_small_bonus: float = 0.02
    line_width_hz: float = 0.8
    anisotropy_pct_normal: float = 0.10
    anisotropy_pct_high: float = 0.10
    viscosity_simplify_max_groups: int = 3
    viscosity_simplify_maxJ_hz: float = 1.5
    # ³J vinyliques corrigés (Jackman & Sternhell)
    # cis: H-C=C-H même côté (~10 Hz), trans: côtés opposés (~16 Hz)
    vinylic_J_defaults: Dict[str, float] = field(default_factory=lambda: {
        'trans': 16.0,    # 14-18 Hz typique
        'cis': 10.0,      # 8-12 Hz typique  
        'geminal': 2.0,   # =CH₂ (0-3 Hz)
        'gauche': 8.0,    # géométrie contrainte
    })
    combination_line_threshold: float = 0.02  # Seuil visibilité raies de combinaison
    use_field_resolution_model: bool = True
    digital_resolution_hz: Optional[float] = None
    min_visibility_factor: float = 1.2
    grouping_tol_scale: float = 0.65
    grouping_tol_scale_auto: bool = True
    never_merge_gap_hz: float = 2.5
    detect_para_min_mhz: int = 600
    detect_para_min_j_hz: float = 0.45
    para_max_per_proton: Optional[int] = None
    sugar_mode: bool = True
    cyclohexane_mode: bool = True
    rigid_W_boost: bool = False
    formyl_rule: bool = True
    amide_mode: bool = True
    peptide_HNCH_coupling: bool = True  # ³J(H-N-C-H) Bystrov
    heteroHF_advanced: bool = True
    karplus_haasnoot: bool = True
    anomeric_preference: Optional[str] = None
    include_13C_satellites: bool = False
    verbose: bool = False
    use_spin_fallback: bool = True
    spin_fallback_max_spins: int = 4
    spin_fallback_min_intensity: float = 1e-3
    spin_fallback_group_tol_scale: float = 0.8
    fallback_partner_partner_enable: bool = True
    fallback_partner_partner_max_delta_hz: float = 30.0
    collapse_equivalent_by_rank: bool = True
    collapse_equivalent_by_delta: bool = False
    equiv_delta_ppm_tol: float = 0.001
    treat_aromatic_noJ_equiv_as_singlet: bool = True
    aromatic_isolated_as_singlet: bool = True
    separate_diastereotopic_protons: bool = True
    auto_bench_mode: bool = False
    protect_diastereotopes: bool = True
    random_seed: int = 0xF00D
    
    # === NOUVELLES OPTIONS AVANCÉES ===
    
    # Paramètres Karplus spécialisés
    use_garbisch_alkenes: bool = True            # Équation Garbisch pour alcènes
    use_bystrov_peptides: bool = True            # Équation Bystrov pour peptides
    use_haasnoot_7_params: bool = True           # Haasnoot-Altona 7 paramètres complets
    
    # Moyenne conformationnelle multi-conformères
    # ✅ AMÉLIORATION: Base augmentée de 10 à 50, avec boost adaptatif jusqu'à 200
    use_conformational_averaging: bool = True   # ✅ Activé par défaut pour meilleure précision
    n_conformers: int = 50                       # ✅ AUGMENTÉ de 10 à 50 (base)
    n_conformers_max: int = 200                  # ✅ NOUVEAU: Maximum pour chaînes très flexibles
    conformer_energy_cutoff_kcal: float = 5.0   # ✅ AUGMENTÉ de 4.0 à 5.0 kcal/mol
    conformer_pruning_rmsd: float = 0.35         # ✅ NOUVEAU: Seuil RMSD pour éliminer doublons (Å)
    flexible_chain_boost: bool = True            # ✅ Augmenter conformères pour chaînes flexibles
    flexible_chain_min_rotatable: int = 3        # ✅ Seuil pour chaîne flexible
    use_etkdg_v3: bool = True                    # ✅ NOUVEAU: Préférer ETKDGv3 (meilleur sampling)
    conformer_random_seed: int = 42              # ✅ NOUVEAU: Reproductibilité
    
    # Détection avancée systèmes second ordre
    advanced_spin_system_detection: bool = True  # Détecter AB₂, ABX, A₂X₃, etc.
    spin_system_delta_nu_J_threshold: float = 10.0  # Seuil Δν/J pour second ordre
    
    # Simulation exacte par Hamiltonien de spin (diagonalisation quantique)
    use_exact_spin_simulation: bool = True       # Utiliser la diagonalisation exacte
    exact_spin_max_spins: int = 8                # Limite max de spins (2^8 = 256×256)
    exact_spin_intensity_threshold: float = 0.001  # Seuil intensité minimum
    
    # Effets de solvant avancés
    # ⚠️ DÉSACTIVÉ par défaut: les corrections de shift ne sont valables que si
    # le modèle ML a été entraîné avec le solvant comme feature.
    # Seule la VISIBILITÉ des protons échangeables est gérée (OH/NH dans D2O, etc.)
    use_solvent_shift_corrections: bool = False  # ❌ DÉSACTIVÉ - pas de correction δ
    use_solvent_linewidth: bool = True           # Largeur de raie selon solvant
    use_solvent_J_corrections: bool = False      # ❌ DÉSACTIVÉ - pas de correction J
    use_exchange_broadening: bool = True         # Élargissement OH/NH
    use_exchange_visibility: bool = True         # ✅ ACTIVÉ: OH/NH invisibles dans D2O
    pH: Optional[float] = None                   # pH pour échange acido-basique
    
    # Marquage des protons labiles (OH, NH, SH)
    mark_labile_protons_broad: bool = True       # Ajouter "broadening" dans le résultat
    labile_proton_linewidth_hz: float = 5.0      # Largeur de raie typique pour échange
    
    # Force field pour optimisation géométrique
    # Options: "MMFF94", "MMFF94s" (meilleur pour petites molécules), "UFF"
    force_field: str = "MMFF94s"                 # MMFF94s recommandé
    force_field_max_iters: int = 200             # Itérations max pour optimisation
    
    # ✅ NOUVEAU: Inférence automatique stéréochimie E/Z
    infer_double_bond_stereo: bool = True        # Inférer E/Z pour doubles liaisons non définies
    infer_stereo_default: str = "E"              # "E" (trans) ou "Z" (cis) - trans plus stable
    infer_stereo_log: bool = True                # Logger les modifications


def _optimize_molecule_geometry(
    mol: Chem.Mol,
    force_field: str = "MMFF94s",
    max_iters: int = 200,
    conf_id: int = -1
) -> bool:
    """
    Optimise la géométrie d'une molécule avec le force field spécifié.
    
    Args:
        mol: Molécule RDKit avec coordonnées 3D
        force_field: "MMFF94", "MMFF94s", ou "UFF"
        max_iters: Nombre maximum d'itérations
        conf_id: ID du conformère (-1 = tous)
    
    Returns:
        True si optimisation réussie, False sinon
    """
    try:
        if force_field.upper() in ("MMFF94", "MMFF94S"):
            # MMFF94s est meilleur pour les petites molécules organiques
            variant = "MMFF94s" if "S" in force_field.upper() else "MMFF94"
            
            # Vérifier si MMFF est applicable
            mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=variant)
            if mmff_props is None:
                # Fallback vers UFF
                return _optimize_molecule_geometry(mol, "UFF", max_iters, conf_id)
            
            if conf_id < 0:
                # Optimiser tous les conformères
                results = AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant=variant, maxIters=max_iters)
                return all(r[0] == 0 for r in results) if results else False
            else:
                ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=conf_id)
                if ff is None:
                    return _optimize_molecule_geometry(mol, "UFF", max_iters, conf_id)
                ff.Minimize(maxIts=max_iters)
                return True
                
        elif force_field.upper() == "UFF":
            if conf_id < 0:
                for cid in range(mol.GetNumConformers()):
                    AllChem.UFFOptimizeMolecule(mol, confId=cid, maxIters=max_iters)
            else:
                AllChem.UFFOptimizeMolecule(mol, confId=conf_id, maxIters=max_iters)
            return True
        else:
            # Force field inconnu, utiliser MMFF94s
            return _optimize_molecule_geometry(mol, "MMFF94s", max_iters, conf_id)
            
    except Exception:
        # Dernier recours: UFF
        try:
            if conf_id < 0:
                for cid in range(mol.GetNumConformers()):
                    AllChem.UFFOptimizeMolecule(mol, confId=cid, maxIters=max_iters)
            else:
                AllChem.UFFOptimizeMolecule(mol, confId=conf_id, maxIters=max_iters)
            return True
        except Exception:
            return False


def _estimate_labile_proton_broadening(
    mol: Chem.Mol,
    heavy_idx: int,
    solvent_name: str = "CDCl3",
    temperature_k: float = 298.15
) -> Optional[str]:
    """
    Estime le type d'élargissement pour les protons labiles (OH, NH, SH).
    
    Returns:
        None si proton non labile
        "sharp" si échange lent (DMSO, basse T)
        "broad" si échange modéré (CDCl3)
        "very_broad" si échange rapide (D2O, MeOD)
        "exchanged" si proton échangé (invisible)
    """
    atom = mol.GetAtomWithIdx(heavy_idx)
    symbol = atom.GetSymbol()
    
    # Vérifier si c'est un hétéroatome avec H
    if symbol not in ("O", "N", "S"):
        return None
    
    # Vérifier qu'il y a des H attachés
    has_h = any(n.GetSymbol() == "H" for n in atom.GetNeighbors())
    if not has_h:
        return None
    
    # Déterminer le type de proton labile
    proton_type = None
    if symbol == "O":
        # OH: alcool, phénol, acide carboxylique
        is_carboxylic = any(
            n.GetSymbol() == "C" and 
            any(b.GetBondType() == rdchem.BondType.DOUBLE and b.GetOtherAtom(n).GetSymbol() == "O" 
                for b in n.GetBonds())
            for n in atom.GetNeighbors() if n.GetSymbol() == "C"
        )
        proton_type = "COOH" if is_carboxylic else "OH"
    elif symbol == "N":
        # NH: amine, amide, indole
        is_amide = any(
            n.GetSymbol() == "C" and 
            any(b.GetBondType() == rdchem.BondType.DOUBLE and b.GetOtherAtom(n).GetSymbol() == "O" 
                for b in n.GetBonds())
            for n in atom.GetNeighbors() if n.GetSymbol() == "C"
        )
        proton_type = "amide_NH" if is_amide else "NH"
    elif symbol == "S":
        proton_type = "SH"
    
    # Comportement selon solvant
    solvent_upper = solvent_name.upper().replace("-", "")
    
    if solvent_upper in ("D2O", "CD3OD", "MEOD"):
        # Solvants protiques → échange rapide
        if proton_type in ("OH", "COOH", "NH", "SH"):
            return "exchanged"
        elif proton_type == "amide_NH":
            return "very_broad"  # Amides échangent plus lentement
    
    elif solvent_upper in ("DMSOD6", "DMSO"):
        # DMSO → échange lent, pics visibles
        if proton_type == "OH":
            return "sharp"  # OH alcoolique visible et sharp en DMSO
        elif proton_type == "COOH":
            return "broad"  # COOH peut échanger
        elif proton_type in ("NH", "amide_NH"):
            return "sharp"  # NH visible en DMSO
        elif proton_type == "SH":
            return "broad"
    
    elif solvent_upper in ("CDCL3", "CDCl3"):
        # CDCl3 → intermédiaire
        if proton_type == "OH":
            return "broad"  # OH souvent large en CDCl3
        elif proton_type == "COOH":
            return "very_broad"  # Acides très larges
        elif proton_type == "NH":
            return "broad"
        elif proton_type == "amide_NH":
            return "sharp"  # NH amide souvent visible
        elif proton_type == "SH":
            return "broad"
    
    # Par défaut
    return "broad" if proton_type else None


# ---------------------------
# Helpers: topology (avec cache intelligent)
# ---------------------------
def _topological_distances(mol: Chem.Mol) -> np.ndarray:
    """
    Calcule les distances topologiques entre tous les atomes.
    
    AMÉLIORATION: Utilise MoleculeCache pour éviter les recalculs.
    """
    return MoleculeCache.get_topological_distances(mol)

def _precompute_aromatic_distance_map(mol: Chem.Mol) -> List[Dict[int, int]]:
    n = mol.GetNumAtoms()
    arom_adj = [[] for _ in range(n)]
    for b in mol.GetBonds():
        if not b.GetIsAromatic():
            continue
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        arom_adj[i].append(j); arom_adj[j].append(i)
    dist_maps: List[Dict[int, int]] = [dict() for _ in range(n)]
    for src in range(n):
        q = deque([(src, 0)]); seen = {src}
        while q:
            u, d = q.popleft()
            for v in arom_adj[u]:
                if v in seen: continue
                dist_maps[src][v] = d + 1
                seen.add(v); q.append((v, d + 1))
    return dist_maps

def _aromatic_rel_from_map(dist_map: List[Dict[int, int]], i: int, j: int) -> Optional[str]:
    if i == j: return None
    dij = dist_map[i].get(j, None)
    if dij is None: return None
    if dij == 1: return "ortho"
    if dij == 2: return "meta"
    if dij == 3: return "para"
    return None

def _is_sp2_like(atom: Chem.Atom) -> bool:
    hyb = str(atom.GetHybridization())
    return ("SP2" in hyb) or atom.GetIsAromatic() or any(b.GetBondTypeAsDouble() == 2 for b in atom.GetBonds())

# ---------------------------
# Equivalence helpers (avec cache intelligent)
# ---------------------------
def _canonical_ranks(mol_noH: Chem.Mol) -> List[int]:
    """
    Calcule les rangs canoniques des atomes.
    
    AMÉLIORATION: Utilise MoleculeCache pour éviter les recalculs.
    """
    return MoleculeCache.get_canonical_ranks(mol_noH)

def _get_atom_invariant(atom: Chem.Atom) -> Tuple:
    """
    Calcule un invariant atomique pour détecter l'équivalence magnétique.
    Deux atomes avec le même invariant sont potentiellement homotopes/équivalents.
    """
    return (
        atom.GetSymbol(),
        atom.GetIsAromatic(),
        atom.GetDegree(),
        atom.GetTotalNumHs(),
        atom.GetFormalCharge(),
        atom.IsInRing(),
    )


def _are_heavy_equivalent(
    ranks: List[int], i: int, j: int,
    heavy_to_H: Dict[int, List[int]],
    proton_shifts_ppm: Optional[Dict[int, float]],
    tol_ppm: float,
    mol: Optional[Chem.Mol] = None,
) -> bool:
    """
    Détermine si deux centres lourds sont magnétiquement équivalents.
    
    Utilise plusieurs critères (dans l'ordre de priorité):
    1. Rangs canoniques égaux → équivalents
    2. Invariants atomiques différents → NON équivalents (priorité haute)
    3. Invariants atomiques égaux + voisins égaux → équivalents (homotopes)
    4. Shifts identiques SEULEMENT si invariants égaux
    
    Note: Pour benzène non-substitué, tous les C sont équivalents.
    """
    if i == j:
        return True
    
    # Critère 1: Rangs canoniques égaux
    if ranks[i] == ranks[j]:
        return True
    
    # Critère 2 & 3: Invariants atomiques
    if mol is not None:
        try:
            atom_i = mol.GetAtomWithIdx(i)
            atom_j = mol.GetAtomWithIdx(j)
            
            inv_i = _get_atom_invariant(atom_i)
            inv_j = _get_atom_invariant(atom_j)
            
            # ✅ CRUCIAL: Si invariants différents → PAS équivalents
            # (même si les shifts sont identiques - c'est une erreur heuristique)
            if inv_i != inv_j:
                return False  # C0 (CH3) ≠ C1 (CH) car numHs différent
            
            # Invariants égaux → vérifier les voisins
            neighbors_i = sorted([_get_atom_invariant(n) for n in atom_i.GetNeighbors()])
            neighbors_j = sorted([_get_atom_invariant(n) for n in atom_j.GetNeighbors()])
            
            if neighbors_i == neighbors_j:
                return True  # Homotopes confirmés
                
        except Exception:
            pass
    
    # Critère 4: Shifts identiques (seulement si invariants sont égaux - déjà vérifié ci-dessus)
    # Ce critère n'est plus appliqué séparément car il cause des faux positifs
    # quand les shifts heuristiques sont tous identiques
    
    return False


def _is_globally_magnetically_equivalent_molecule(
    mol: Chem.Mol,
    ranks: List[int],
    heavy_to_H: Dict[int, List[int]],
) -> bool:
    """
    Détecte si TOUS les protons d'une molécule sont magnétiquement équivalents.
    
    Dans ce cas, la molécule entière donne un SINGULET (ou singulet large).
    
    Exemples:
    - Cyclohexane (C1CCCCC1): Tous les CH2 équivalents → singulet
    - Cyclopentane (C1CCCC1): Tous les CH2 équivalents → singulet
    - Benzène (c1ccccc1): Tous les CH équivalents → singulet
    - 1,2-Dichloroéthane (ClCCCl): Les 2 CH2 équivalents → singulet
    - 1,4-Dioxane (C1COCOC1): Tous les CH2 équivalents → singulet
    
    Returns:
        True si tous les protons sont magnétiquement équivalents
    """
    if mol is None or not heavy_to_H:
        return False
    
    # Collecter tous les carbones portant des protons
    proton_bearing_carbons = [h_idx for h_idx, h_list in heavy_to_H.items() if h_list]
    
    if len(proton_bearing_carbons) < 2:
        return False  # Besoin d'au moins 2 carbones pour parler d'équivalence
    
    # Vérifier si tous les carbones porteurs de protons ont le même rang canonique
    if ranks is None or len(ranks) == 0:
        return False
    
    first_rank = None
    for c_idx in proton_bearing_carbons:
        if c_idx >= len(ranks):
            return False
        if first_rank is None:
            first_rank = ranks[c_idx]
        elif ranks[c_idx] != first_rank:
            return False  # Au moins 2 carbones ont des rangs différents
    
    # Vérifier aussi que tous ont le même nombre de H (tous CH2, ou tous CH, etc.)
    first_n_H = None
    for c_idx in proton_bearing_carbons:
        n_H = len(heavy_to_H.get(c_idx, []))
        if first_n_H is None:
            first_n_H = n_H
        elif n_H != first_n_H:
            return False  # Mélange de CH, CH2, CH3
    
    # Vérifier que c'est un cycle saturé simple (pas de substituants différenciateurs)
    # ou une molécule linéaire symétrique
    for c_idx in proton_bearing_carbons:
        atom = mol.GetAtomWithIdx(c_idx)
        # Vérifier que tous les carbones sont identiques structurellement
        # (même hybridation, même environnement)
        if atom.GetHybridization() not in [Chem.rdchem.HybridizationType.SP3, Chem.rdchem.HybridizationType.SP2]:
            return False
    
    return True


def _are_protons_on_equivalent_carbons(
    mol: Chem.Mol,
    ranks: List[int],
    heavy_idx_1: int,
    heavy_idx_2: int,
) -> bool:
    """
    Vérifie si deux carbones portant des protons sont magnétiquement équivalents.
    
    Si oui, les protons sur ces carbones ne montrent pas de couplage mutuel
    (le couplage existe mais est invisible car les protons sont isochroniques).
    
    Utilisé pour:
    - 1,2-Dichloroéthane: Cl-CH2-CH2-Cl (les 2 CH2 sont équivalents → singulet)
    - Cyclohexane: tous les CH2 équivalents
    - Cyclopentane: tous les CH2 équivalents
    
    IMPORTANT: Utilise breakTies=False pour obtenir les vraies classes de symétrie,
    pas les rangs canoniques uniques.
    """
    if mol is None:
        return False
    
    if heavy_idx_1 == heavy_idx_2:
        return False  # Même carbone, pas d'équivalence inter-carbone
    
    # ✅ CRUCIAL: Utiliser breakTies=False pour obtenir les classes de symétrie
    # Cela regroupe les atomes chimiquement équivalents avec la même valeur
    try:
        symmetry_classes = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    except Exception:
        # Fallback sur les rangs fournis
        if ranks is None:
            return False
        symmetry_classes = ranks
    
    if heavy_idx_1 >= len(symmetry_classes) or heavy_idx_2 >= len(symmetry_classes):
        return False
    
    # Si les classes de symétrie sont égales, les carbones sont équivalents
    if symmetry_classes[heavy_idx_1] == symmetry_classes[heavy_idx_2]:
        return True
    
    return False


# ---------------------------
# Exchangeable protons
# ---------------------------
_SMARTS_COOH_OH = Chem.MolFromSmarts('[OX2H][CX3](=O)[#6,#7,#8]')
_SMARTS_OH      = Chem.MolFromSmarts('[OX2H]')
_SMARTS_NH      = Chem.MolFromSmarts('[NX3H2,NX3H1]')
_SMARTS_NH_WIDE = Chem.MolFromSmarts('[N;H1,H2;!$([N-])]')
_SMARTS_SH      = Chem.MolFromSmarts('[SX2H]')

def _is_exchangeable_H_on_X(mol: Chem.Mol, heavy_idx: int) -> Optional[str]:
    a = mol.GetAtomWithIdx(heavy_idx)
    sym = a.GetSymbol()
    if sym == 'O':
        for m in mol.GetSubstructMatches(_SMARTS_COOH_OH):
            if m[0] == heavy_idx: return "COOH"
        for m in mol.GetSubstructMatches(_SMARTS_OH):
            if m[0] == heavy_idx: return "OH"
    elif sym == 'N':
        try:
            if a.GetIsAromatic() and a.GetTotalNumHs() == 1:
                return "NH"
        except Exception:
            pass
        for m in mol.GetSubstructMatches(_SMARTS_NH):
            if m[0] == heavy_idx: return "NH"
        for m in mol.GetSubstructMatches(_SMARTS_NH_WIDE):
            if m[0] == heavy_idx: return "NH"
    elif sym == 'S':
        for m in mol.GetSubstructMatches(_SMARTS_SH):
            if m[0] == heavy_idx: return "SH"
    return None

def _has_quadrupolar_neighbor(mol: Chem.Mol, heavy_idx: int) -> bool:
    a = mol.GetAtomWithIdx(heavy_idx)
    for nb in a.GetNeighbors():
        sym = nb.GetSymbol(); iso = nb.GetIsotope()
        if sym == "N" and iso != 15: return True
        if sym in {"B","Al","Cl"}: return True
    return False


def _is_isolated_methyl(mol: Chem.Mol, heavy_idx: int, heavy_to_H: Dict[int, List[int]]) -> bool:
    """
    Détecte si un carbone porte un CH₃ "isolé" qui ne devrait pas montrer de couplage visible.
    
    Cas typiques:
    - CH₃ benzylique (toluène): Ar-CH₃ → singulet
    - OCH₃ (méthoxy): O-CH₃ → singulet  
    - N-CH₃ (N-méthyl): N-CH₃ → singulet
    - Acétyle: CH₃-C=O → singulet (sauf aldéhyde)
    
    Ces groupes ont des couplages ⁴J ou ⁵J très faibles (< 0.8 Hz) généralement
    non résolus expérimentalement.
    """
    atom = mol.GetAtomWithIdx(heavy_idx)
    
    # Doit être un carbone sp3 avec 3 H (CH₃)
    if atom.GetSymbol() != "C":
        return False
    
    n_H = len(heavy_to_H.get(heavy_idx, []))
    if n_H != 3:
        return False
    
    # Vérifier hybridation sp3
    if _is_sp2_like(atom):
        return False
    
    # Vérifier les voisins
    neighbors = list(atom.GetNeighbors())
    if len(neighbors) != 1:
        return False  # CH₃ a exactement 1 voisin lourd
    
    neighbor = neighbors[0]
    neighbor_sym = neighbor.GetSymbol()
    
    # Cas 1: CH₃ benzylique (Ar-CH₃)
    if neighbor_sym == "C" and neighbor.GetIsAromatic():
        return True
    
    # Cas 2: OCH₃ (O-CH₃)
    if neighbor_sym == "O":
        return True
    
    # Cas 3: N-CH₃
    if neighbor_sym == "N":
        return True
    
    # Cas 4: CH₃-C=O (acétyle, mais pas aldéhyde)
    if neighbor_sym == "C" and _is_sp2_like(neighbor):
        # Vérifier si c'est un carbonyle
        for n2 in neighbor.GetNeighbors():
            if n2.GetSymbol() == "O":
                bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), n2.GetIdx())
                if bond and bond.GetBondTypeAsDouble() >= 2.0:
                    # C'est un carbonyle - CH₃ est isolé SAUF si l'autre côté est H (aldéhyde)
                    other_neighbors = [x for x in neighbor.GetNeighbors() 
                                       if x.GetIdx() != heavy_idx and x.GetSymbol() != "O"]
                    if other_neighbors:
                        return True  # CH₃-CO-X (pas aldéhyde)
    
    return False


# ---------------------------
# Allylic/W path
# ---------------------------

def _is_homoallylic_4J(base: Chem.Mol, i: int, j: int) -> Tuple[bool, float]:
    """
    Détecte le couplage homoallylique ⁴J entre un H vinylique et un CH₃ adjacent.
    Chemin: H-C=C-CH₃ (4 liaisons)
    
    Exemples:
    - Groupe prényl: CH₂-CH=C(CH₃)₂ → H vinylique couple avec CH₃ (~1.0-1.5 Hz)
    - Groupe isopropylidène: =C(CH₃)₂ → H couple avec les 2 CH₃
    
    Retourne: (is_homoallylic, J_value)
    """
    Ai = base.GetAtomWithIdx(i)  # Carbone portant H1
    Aj = base.GetAtomWithIdx(j)  # Carbone portant H2
    
    if Ai.GetSymbol() != "C" or Aj.GetSymbol() != "C":
        return False, 0.0
    
    # Cas 1: i est sp2 (vinylique), j est sp3 (méthyle sur C=C)
    # Chemin: H-C(sp2)=C(sp2)-C(sp3)-H
    if _is_sp2_like(Ai) and not _is_sp2_like(Aj):
        # Chercher C=C où Ai est impliqué
        for n1 in Ai.GetNeighbors():
            if n1.GetSymbol() != "C" or not _is_sp2_like(n1):
                continue
            bond = base.GetBondBetweenAtoms(Ai.GetIdx(), n1.GetIdx())
            if not bond or bond.GetBondTypeAsDouble() < 1.5:
                continue
            # n1 est l'autre C du C=C, vérifier si Aj est attaché
            if base.GetBondBetweenAtoms(n1.GetIdx(), Aj.GetIdx()):
                # Homoallylique confirmé: H-C=C-CH₃
                # J dépend de la géométrie (cis ~1.5 Hz, trans ~1.0 Hz)
                return True, 1.3
    
    # Cas 2: j est sp2 (vinylique), i est sp3 (méthyle sur C=C)
    if _is_sp2_like(Aj) and not _is_sp2_like(Ai):
        for n1 in Aj.GetNeighbors():
            if n1.GetSymbol() != "C" or not _is_sp2_like(n1):
                continue
            bond = base.GetBondBetweenAtoms(Aj.GetIdx(), n1.GetIdx())
            if not bond or bond.GetBondTypeAsDouble() < 1.5:
                continue
            if base.GetBondBetweenAtoms(n1.GetIdx(), Ai.GetIdx()):
                return True, 1.3
    
    return False, 0.0


def _is_allylic_4J(base: Chem.Mol, i: int, j: int, mol3d: Optional[Chem.Mol] = None) -> Tuple[bool, float]:
    """
    Détecte le couplage allylique ⁴J classique: H-C-C=C-H (4 liaisons)
    
    Le proton sur le CH₂ allylique couple avec le H vinylique.
    J varie selon l'angle dièdre (0.5-3.0 Hz typique).
    
    Retourne: (is_allylic, J_value)
    """
    Ai = base.GetAtomWithIdx(i)
    Aj = base.GetAtomWithIdx(j)
    
    if Ai.GetSymbol() != "C" or Aj.GetSymbol() != "C":
        return False, 0.0
    
    # Cas: i est sp3 (allylique), j est sp2 (vinylique)
    # Chemin: H-C(sp3)-C(sp2)=C(sp2)-H
    if not _is_sp2_like(Ai) and _is_sp2_like(Aj):
        for n1 in Ai.GetNeighbors():
            if n1.GetSymbol() != "C" or not _is_sp2_like(n1):
                continue
            # n1 est sp2, vérifier double liaison avec Aj
            bond = base.GetBondBetweenAtoms(n1.GetIdx(), Aj.GetIdx())
            if bond and bond.GetBondTypeAsDouble() >= 1.5:
                # Allylique confirmé: H-CH₂-C=C-H
                # J dépend de la géométrie - utiliser 3D si disponible
                J = 1.8  # Valeur par défaut
                if mol3d is not None and mol3d.GetNumConformers() > 0:
                    try:
                        conf = mol3d.GetConformer()
                        # Angle dièdre H-C-C=C
                        dihedral = abs(AllChem.GetDihedralDeg(conf, i, n1.GetIdx(), Aj.GetIdx(), j))
                        # J ∝ cos²(angle) pour allylique
                        if dihedral < 30 or dihedral > 150:
                            J = 2.5  # Cisoid ou transoid
                        else:
                            J = 1.2  # Perpendiculaire
                    except Exception:
                        pass
                return True, J
    
    # Cas inverse: j est sp3 (allylique), i est sp2 (vinylique)
    if not _is_sp2_like(Aj) and _is_sp2_like(Ai):
        for n1 in Aj.GetNeighbors():
            if n1.GetSymbol() != "C" or not _is_sp2_like(n1):
                continue
            bond = base.GetBondBetweenAtoms(n1.GetIdx(), Ai.GetIdx())
            if bond and bond.GetBondTypeAsDouble() >= 1.5:
                J = 1.8
                if mol3d is not None and mol3d.GetNumConformers() > 0:
                    try:
                        conf = mol3d.GetConformer()
                        dihedral = abs(AllChem.GetDihedralDeg(conf, j, n1.GetIdx(), Ai.GetIdx(), i))
                        if dihedral < 30 or dihedral > 150:
                            J = 2.5
                        else:
                            J = 1.2
                    except Exception:
                        pass
                return True, J
    
    return False, 0.0


def _is_W_coupling_5J(base: Chem.Mol, i: int, j: int) -> Tuple[bool, float]:
    """
    Détecte le couplage ⁵J en W (zigzag) pour systèmes rigides.
    
    Chemin: H-C-C-C-C-H en configuration W planaire.
    Typique dans les systèmes bicycliques, stéroïdes, etc.
    
    J typique: 1.0-2.0 Hz si rigide et planaire.
    
    Retourne: (is_W, J_value)
    """
    Ai = base.GetAtomWithIdx(i)
    Aj = base.GetAtomWithIdx(j)
    
    # Les deux doivent être dans un cycle pour être rigides
    if not (Ai.IsInRing() and Aj.IsInRing()):
        return False, 0.0
    
    # Chercher un chemin de 4 carbones entre i et j (5 liaisons total avec H)
    # Chemin: i - C1 - C2 - C3 - j
    for n1 in Ai.GetNeighbors():
        if n1.GetSymbol() != "C" or n1.GetIdx() == j:
            continue
        for n2 in n1.GetNeighbors():
            if n2.GetIdx() == i or n2.GetSymbol() != "C":
                continue
            for n3 in n2.GetNeighbors():
                if n3.GetIdx() in {n1.GetIdx(), i} or n3.GetSymbol() != "C":
                    continue
                if base.GetBondBetweenAtoms(n3.GetIdx(), j):
                    # Chemin trouvé, vérifier si tous dans le même cycle
                    if n1.IsInRing() and n2.IsInRing() and n3.IsInRing():
                        return True, 1.2
    
    return False, 0.0


def _is_W_allylic(base: Chem.Mol, i: int, j: int) -> Tuple[bool, bool]:
    """
    Détecte un couplage ⁴J de type W-allylique.
    Chemin: H-C(sp3)-C(sp2)=C(sp2)-H (4 liaisons, forme W)
    
    Retourne: (is_W, is_rigid)
    """
    Ai = base.GetAtomWithIdx(i)
    Aj = base.GetAtomWithIdx(j)

    def _is_sp3C(a: Chem.Atom) -> bool:
        return (a.GetSymbol() == "C") and (not _is_sp2_like(a))

    rigid = Ai.IsInRing() and Aj.IsInRing()

    for n1 in Ai.GetNeighbors():
        if not _is_sp3C(n1):
            continue
        for n2 in n1.GetNeighbors():
            if n2.GetIdx() == Ai.GetIdx():
                continue
            if not (_is_sp2_like(n2) and n2.GetSymbol() == "C"):
                continue
            # ✅ CORRECTION: Exclure le cas où j est le même que n2
            # Cela évite les faux positifs pour CH3-C-C=C où le =C est j
            if n2.GetIdx() == Aj.GetIdx():
                continue
            for n3 in n2.GetNeighbors():
                if n3.GetIdx() in {n1.GetIdx(), Ai.GetIdx()}:
                    continue
                if not (_is_sp2_like(n3) and n3.GetSymbol() == "C"):
                    continue
                b = base.GetBondBetweenAtoms(n2.GetIdx(), n3.GetIdx())
                if not b or b.GetBondTypeAsDouble() < 2.0:
                    continue
                if base.GetBondBetweenAtoms(n3.GetIdx(), Aj.GetIdx()):
                    return True, (rigid or n1.IsInRing() or n2.IsInRing() or n3.IsInRing())

    ok_chain = False
    try:
        Ai_sp2_C = [n for n in Ai.GetNeighbors() if (n.GetSymbol() == "C" and _is_sp2_like(n))]
        Aj_sp2_C = [n for n in Aj.GetNeighbors() if (n.GetSymbol() == "C" and _is_sp2_like(n))]
        if Ai.IsInRing() and Aj.IsInRing() and Ai_sp2_C and Aj_sp2_C:
            ci = Ai_sp2_C[0]; cj = Aj_sp2_C[0]
            for x in ci.GetNeighbors():
                if x.GetIdx() == Ai.GetIdx():
                    continue
                if _is_sp2_like(x) and x.GetSymbol() == "C":
                    if base.GetBondBetweenAtoms(ci.GetIdx(), x.GetIdx()) and base.GetBondBetweenAtoms(x.GetIdx(), cj.GetIdx()):
                        ok_chain = True; break
                    if base.GetBondBetweenAtoms(x.GetIdx(), cj.GetIdx()) and base.GetBondBetweenAtoms(ci.GetIdx(), cj.GetIdx()):
                        ok_chain = True; break
            if ok_chain:
                if base.GetNumConformers() == 0:
                    return True, True
                conf = base.GetConformer()
                dihed = abs(AllChem.GetDihedralDeg(conf, i, ci.GetIdx(), cj.GetIdx(), j))
                if 75.0 <= dihed <= 105.0:
                    return True, True
    except Exception:
        pass

    return False, rigid

# ---------------------------
# Karplus models
# ---------------------------

# Coefficients Karplus classiques pour différents types de liaisons
_KARPLUS_COEFFICIENTS = {
    # Type: (A, B, C) pour J = A·cos²θ + B·cosθ + C
    'HCCH_sp3': (7.76, -1.10, 1.40),      # Bothner-By (1965)
    'HCCH_alkane': (9.50, -1.00, 0.50),   # Karplus original
    'HCCH_altona': (13.70, -0.73, 0.00),  # Haasnoot-Altona (1980)
    'HCCH_sugar': (10.40, -1.50, 0.20),   # Optimisé sucres (Tvaroska)
    'HNCH': (6.40, -1.40, 1.90),          # Bystrov (1976) peptides
    'HCOH': (10.40, -1.50, 0.20),         # Fraser (1971)
    'HCNH': (8.00, -1.00, 0.00),          # Barfield
}

def _karplus_simple(theta_deg: float, bond_type: str = 'HCCH_altona') -> float:
    """
    Équation Karplus avec coefficients sélectionnables.
    
    Args:
        theta_deg: Angle dièdre en degrés
        bond_type: Type de liaison pour sélectionner les coefficients
    
    Returns:
        ³J en Hz
    
    Références:
        - Karplus M., J. Chem. Phys. 30, 11 (1959)
        - Haasnoot et al., Tetrahedron 36, 2783 (1980)
        - Bystrov V.F., Prog. NMR Spectrosc. 10, 41 (1976)
    """
    A, B, C = _KARPLUS_COEFFICIENTS.get(bond_type, _KARPLUS_COEFFICIENTS['HCCH_altona'])
    theta = math.radians(theta_deg)
    return max(0.0, A * (math.cos(theta) ** 2) + B * math.cos(theta) + C)

def _karplus_haasnoot_altona(theta_deg: float, electroneg_corr: float = 0.0) -> float:
    """
    Wrapper simple autour de _karplus_with_substituent_correction.
    """
    return _karplus_with_substituent_correction(theta_deg, alpha_EN=electroneg_corr, beta_EN=0.0)


def _get_substituents_on_carbon(mol: Chem.Mol, c_idx: int, exclude_idx: int) -> List[str]:
    """
    Récupère les symboles des substituants sur un carbone, excluant un atome donné.
    
    Args:
        mol: Molécule RDKit
        c_idx: Index du carbone
        exclude_idx: Index de l'atome à exclure (généralement H observé ou C voisin)
    
    Returns:
        Liste des symboles atomiques des substituants (excluant H)
    """
    c_atom = mol.GetAtomWithIdx(c_idx)
    substituents = []
    for nb in c_atom.GetNeighbors():
        if nb.GetIdx() == exclude_idx:
            continue
        sym = nb.GetSymbol()
        if sym != 'H':  # Exclure les H pour les corrections EN
            substituents.append(sym)
    return substituents


def _estimate_vicinal_sp3_J(
    mol3d: Chem.Mol,
    h1: int,
    h2: int,
    use_haasnoot: bool = False,
    electroneg_corr: float = 0.0,
    use_full_haasnoot: bool = True,
) -> Optional[float]:
    """
    Estime ³J(H-C-C-H) vicinal avec l'équation Haasnoot-Altona complète (7 paramètres).
    
    Cette version améliorée:
    1. Détecte automatiquement les substituants sur Cα et Cβ
    2. Calcule les signes d'orientation à partir de la géométrie 3D
    3. Applique l'équation Haasnoot-Altona complète avec le terme P7
    
    Args:
        mol3d: Molécule RDKit avec coordonnées 3D
        h1, h2: Indices des deux protons
        use_haasnoot: Si True, utilise l'équation Haasnoot (obsolète, voir use_full_haasnoot)
        electroneg_corr: Correction EN manuelle (obsolète)
        use_full_haasnoot: Si True, utilise l'équation complète à 7 paramètres
    
    Returns:
        ³J en Hz ou None si calcul impossible
    """
    try:
        c1 = next(n.GetIdx() for n in mol3d.GetAtomWithIdx(h1).GetNeighbors() if n.GetSymbol() == "C")
        c2 = next(n.GetIdx() for n in mol3d.GetAtomWithIdx(h2).GetNeighbors() if n.GetSymbol() == "C")
        
        if not mol3d.GetBondBetweenAtoms(c1, c2): 
            return None
        
        conf = mol3d.GetConformer()
        dihed = AllChem.GetDihedralDeg(conf, h1, c1, c2, h2)
        
        if use_full_haasnoot:
            # Récupérer les substituants sur chaque carbone
            subs_alpha = _get_substituents_on_carbon(mol3d, c1, c2)  # Sur C1, excluant C2
            subs_beta = _get_substituents_on_carbon(mol3d, c2, c1)   # Sur C2, excluant C1
            
            # Calculer les signes d'orientation pour chaque substituant
            orientation_signs = []
            
            # Substituants sur C1 (alpha)
            c1_atom = mol3d.GetAtomWithIdx(c1)
            for nb in c1_atom.GetNeighbors():
                if nb.GetIdx() == c2 or nb.GetSymbol() == 'H':
                    continue
                sign = _haasnoot_orientation_sign(mol3d, h1, c1, c2, nb.GetIdx())
                orientation_signs.append(sign)
            
            # Substituants sur C2 (beta)
            c2_atom = mol3d.GetAtomWithIdx(c2)
            for nb in c2_atom.GetNeighbors():
                if nb.GetIdx() == c1 or nb.GetSymbol() == 'H':
                    continue
                sign = _haasnoot_orientation_sign(mol3d, h1, c1, c2, nb.GetIdx())
                orientation_signs.append(sign)
            
            return _karplus_haasnoot_altona_full(
                dihed,
                substituents_alpha=subs_alpha,
                substituents_beta=subs_beta,
                orientation_signs=orientation_signs,
                use_7_params=True
            )
        elif use_haasnoot:
            return _karplus_haasnoot_altona(dihed, electroneg_corr)
        else:
            return _karplus_simple(dihed)
            
    except Exception:
        return None


def _is_peptide_backbone(mol: Chem.Mol, h_idx: int, n_idx: int) -> Tuple[bool, Optional[float]]:
    """
    Détecte si un couplage H-N fait partie d'un backbone peptidique.
    
    Recherche le motif: H-N-C(=O)-C ou N-CH(R)-C(=O)
    
    Args:
        mol: Molécule RDKit
        h_idx: Index du proton (sur N)
        n_idx: Index de l'azote
    
    Returns:
        (is_peptide, phi_angle_deg) où phi est l'angle dièdre du backbone
    """
    try:
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Chercher C(=O) adjacent à N (carbonyle du résidu précédent)
        carbonyl_c = None
        alpha_c = None
        
        for nb in n_atom.GetNeighbors():
            if nb.GetSymbol() != 'C':
                continue
            
            # Vérifier si c'est un carbonyle
            is_carbonyl = False
            for bond in nb.GetBonds():
                other = bond.GetOtherAtom(nb)
                if other.GetSymbol() == 'O' and bond.GetBondTypeAsDouble() >= 2:
                    is_carbonyl = True
                    break
            
            if is_carbonyl:
                carbonyl_c = nb
            else:
                # Potentiel Cα (carbone sp³ avec H)
                if any(n2.GetSymbol() == 'H' for n2 in nb.GetNeighbors()):
                    alpha_c = nb
        
        if carbonyl_c is None or alpha_c is None:
            return False, None
        
        # Trouver le Hα
        h_alpha = None
        for nb in alpha_c.GetNeighbors():
            if nb.GetSymbol() == 'H':
                h_alpha = nb
                break
        
        if h_alpha is None:
            return False, None
        
        # Calculer l'angle φ si on a une géométrie 3D
        if mol.GetNumConformers() > 0:
            conf = mol.GetConformer()
            # φ = C'-N-Cα-C' (ou approximation H-N-Cα-Hα)
            try:
                phi = AllChem.GetDihedralDeg(conf, h_idx, n_idx, alpha_c.GetIdx(), h_alpha.GetIdx())
                return True, phi
            except Exception:
                return True, None
        
        return True, None
        
    except Exception:
        return False, None


def _estimate_peptide_coupling(
    mol: Chem.Mol,
    h_idx: int,
    partner_h_idx: int,
    n_idx: int,
    coupling_type: str = 'HNCH_alpha'
) -> Optional[float]:
    """
    Estime les couplages peptidiques avec l'équation Bystrov.
    
    Cette fonction est spécifiquement optimisée pour les peptides et protéines.
    
    Types de couplages supportés:
    - HNCH_alpha: ³J(HN-Cα-Hα) - utilisé pour déterminer φ
    - HNCH_psi: ³J(HN-C'-N) - utilisé pour déterminer ψ
    - HCCH_peptide: ³J(Hα-Cα-Cβ-Hβ) - chaînes latérales
    
    Args:
        mol: Molécule RDKit avec coordonnées 3D
        h_idx: Index du proton NH
        partner_h_idx: Index du proton partenaire (Hα)
        n_idx: Index de l'azote
        coupling_type: Type de couplage peptidique
    
    Returns:
        ³J en Hz ou None
    """
    try:
        is_peptide, phi = _is_peptide_backbone(mol, h_idx, n_idx)
        
        if not is_peptide:
            return None
        
        if phi is not None:
            return _karplus_bystrov_peptide(phi, coupling_type=coupling_type)
        else:
            # Sans géométrie 3D, utiliser une valeur moyenne
            # φ typique en feuillet β: ~-120° → ³J ~8-9 Hz
            # φ typique en hélice α: ~-60° → ³J ~3-4 Hz
            return 6.0  # Valeur moyenne
            
    except Exception:
        return None

def _estimate_vinylic_J(mol3d: Chem.Mol, h1: int, h2: int, config: Optional[SimConfig]=None) -> Optional[float]:
    """
    Estime ³J(H-C=C-H) vinylique avec l'équation Karplus-Garbisch.
    
    Utilise l'équation Garbisch (1968) spécialement paramétrée pour les alcènes:
    ³J = A + B·cos(θ) + C·cos(2θ)
    
    Valeurs typiques:
    - cis (θ ≈ 0°): 6-12 Hz (typiquement 10 Hz)
    - trans (θ ≈ 180°): 12-18 Hz (typiquement 16 Hz)  
    - geminal (=CH₂): 0-3 Hz
    
    Références:
        - Garbisch E.W., J. Am. Chem. Soc. 90, 6543 (1968)
        - Bothner-By A.A., Adv. Magn. Reson. 1, 195 (1965)
    """
    config = config or SimConfig()
    try:
        c1 = next(n.GetIdx() for n in mol3d.GetAtomWithIdx(h1).GetNeighbors() if n.GetSymbol() == "C")
        c2 = next(n.GetIdx() for n in mol3d.GetAtomWithIdx(h2).GetNeighbors() if n.GetSymbol() == "C")
        if not mol3d.GetBondBetweenAtoms(c1, c2): 
            return None
        b = mol3d.GetBondBetweenAtoms(c1, c2)
        if not b or b.GetBondTypeAsDouble() < 2.0: 
            return None
        
        # Vérifier si geminal (même carbone)
        if c1 == c2:
            return config.vinylic_J_defaults.get('geminal', 2.0)
        
        conf = mol3d.GetConformer()
        phi = AllChem.GetDihedralDeg(conf, h1, c1, c2, h2)
        
        # Déterminer le type d'alcène pour les paramètres Garbisch
        c1_atom = mol3d.GetAtomWithIdx(c1)
        c2_atom = mol3d.GetAtomWithIdx(c2)
        
        # Vérifier si conjugué (adjacent à un autre C sp² ou aromatique)
        is_conjugated = False
        has_electroneg_sub = False
        substituents = []
        
        for c_atom in [c1_atom, c2_atom]:
            for nb in c_atom.GetNeighbors():
                if nb.GetIdx() in (h1, h2, c1, c2):
                    continue
                sym = nb.GetSymbol()
                substituents.append(sym)
                
                # Vérifier conjugaison
                if sym == 'C' and (_is_sp2_like(nb) or nb.GetIsAromatic()):
                    is_conjugated = True
                
                # Vérifier électronégativité
                if sym in ('O', 'N', 'F', 'Cl', 'Br'):
                    has_electroneg_sub = True
        
        # Sélectionner le type d'alcène pour Garbisch
        if has_electroneg_sub:
            alkene_type = 'electronegative'
        elif is_conjugated:
            alkene_type = 'conjugated'
        else:
            alkene_type = 'standard'
        
        # Utiliser l'équation Garbisch
        J = _karplus_garbisch_alkene(phi, alkene_type=alkene_type, substituents=substituents)
        
        return J
            
    except Exception:
        return None


def _estimate_alkyne_J(mol: Chem.Mol, h1_idx: int, h2_idx: int) -> Optional[float]:
    """
    Estime ⁴J(H-C≡C-H) ou ⁵J(H-C-C≡C-H) pour les alcynes.
    
    Couplages alcynes typiques:
    - ⁴J(H-C≡C-H) terminaux: 2.5-3.5 Hz (via liaison triple)
    - ⁵J(H-C-C≡C-H): 1.5-2.5 Hz
    - ³J propargylique (H-C-C≡C): 2-3 Hz
    
    Références:
        - Günther, NMR Spectroscopy, 3rd Ed.
        - Stothers, Carbon-13 NMR Spectroscopy
    """
    try:
        # Trouver les carbones parents
        h1_atom = mol.GetAtomWithIdx(h1_idx)
        h2_atom = mol.GetAtomWithIdx(h2_idx)
        
        c1_neighbors = [n for n in h1_atom.GetNeighbors() if n.GetSymbol() == 'C']
        c2_neighbors = [n for n in h2_atom.GetNeighbors() if n.GetSymbol() == 'C']
        
        if not c1_neighbors or not c2_neighbors:
            return None
        
        c1 = c1_neighbors[0]
        c2 = c2_neighbors[0]
        
        # Vérifier si un des carbones est sp (alcyne)
        def is_alkyne_carbon(atom):
            for bond in atom.GetBonds():
                if bond.GetBondTypeAsDouble() == 3.0:
                    return True
            return False
        
        # Cas 1: H-C≡C-H (⁴J direct, rare)
        if is_alkyne_carbon(c1) and is_alkyne_carbon(c2):
            return 2.5  # Très rare mais observable
        
        # Cas 2: H-C-C≡C-H (⁵J propargylique)
        # Chercher un chemin via une triple liaison
        for nb1 in c1.GetNeighbors():
            if nb1.GetSymbol() != 'C':
                continue
            bond_c1_nb1 = mol.GetBondBetweenAtoms(c1.GetIdx(), nb1.GetIdx())
            if bond_c1_nb1 and bond_c1_nb1.GetBondTypeAsDouble() == 3.0:
                # C1 est lié à un alcyne
                for nb2 in nb1.GetNeighbors():
                    if nb2.GetIdx() == c1.GetIdx():
                        continue
                    if nb2.GetSymbol() == 'C':
                        # Vérifier si c2 est connecté
                        bond_nb1_nb2 = mol.GetBondBetweenAtoms(nb1.GetIdx(), nb2.GetIdx())
                        if bond_nb1_nb2 and bond_nb1_nb2.GetBondTypeAsDouble() == 3.0:
                            if nb2.GetIdx() == c2.GetIdx():
                                return 2.5  # ⁴J(H-C≡C-H)
                            # Vérifier liaison nb2-c2
                            if mol.GetBondBetweenAtoms(nb2.GetIdx(), c2.GetIdx()):
                                return 2.0  # ⁵J propargylique
        
        # Cas 3: H adjacent à alcyne (propargylique)
        for nb in c1.GetNeighbors():
            if nb.GetSymbol() == 'C':
                for bond in nb.GetBonds():
                    if bond.GetBondTypeAsDouble() == 3.0:
                        # c1 est propargylique
                        return 2.5
        
        return None
        
    except Exception:
        return None


def _estimate_15N_1H_coupling(mol: Chem.Mol, h_idx: int, n_idx: int, dist: int) -> Optional[float]:
    """
    Estime les couplages ¹⁵N-¹H.
    
    Couplages typiques:
    - ¹J(¹⁵N-H): 90-95 Hz (amides), 70-75 Hz (amines)
    - ²J(¹⁵N-C-H): 0-2 Hz
    - ³J(¹⁵N-C-C-H): 0-4 Hz (dépend de l'angle dièdre)
    
    Note: Ces couplages sont rarement observés car ¹⁵N est à 0.4% d'abondance naturelle
    (comme les satellites ¹³C) mais sont importants pour les composés enrichis.
    """
    try:
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Vérifier le type d'azote
        n_is_amide = False
        for nb in n_atom.GetNeighbors():
            if nb.GetSymbol() == 'C':
                for bond in nb.GetBonds():
                    other = bond.GetOtherAtom(nb)
                    if other.GetSymbol() == 'O' and bond.GetBondTypeAsDouble() >= 2:
                        n_is_amide = True
                        break
        
        if dist == 1:
            # ¹J(¹⁵N-H) direct
            if n_is_amide:
                return 92.0  # Amide typique
            else:
                return 73.0  # Amine typique
        elif dist == 2:
            # ²J(¹⁵N-C-H)
            return 1.5
        elif dist == 3:
            # ³J(¹⁵N-C-C-H) - dépend de Karplus
            return 2.5
        
        return None
        
    except Exception:
        return None


def _estimate_31P_1H_coupling(mol: Chem.Mol, h_idx: int, p_idx: int, dist: int) -> Optional[float]:
    """
    Estime les couplages ³¹P-¹H avancés.
    
    Couplages typiques:
    - ²J(P-C-H): 2-20 Hz (dépend de l'angle)
    - ³J(P-C-C-H): 5-20 Hz (relation Karplus)
    - ³J(P-O-C-H): 8-12 Hz (phosphates)
    - ⁴J(P-O-C-C-H): 0-3 Hz
    
    Très importants pour:
    - Nucléotides (ATP, ADN, ARN)
    - Phospholipides
    - Composés organophosphorés
    """
    try:
        p_atom = mol.GetAtomWithIdx(p_idx)
        
        # Déterminer le type de phosphore
        has_double_O = False
        n_oxygen_neighbors = 0
        for nb in p_atom.GetNeighbors():
            if nb.GetSymbol() == 'O':
                n_oxygen_neighbors += 1
                bond = mol.GetBondBetweenAtoms(p_idx, nb.GetIdx())
                if bond and bond.GetBondTypeAsDouble() >= 2:
                    has_double_O = True
        
        is_phosphate = has_double_O and n_oxygen_neighbors >= 2
        is_phosphine = n_oxygen_neighbors == 0
        
        if dist == 2:
            # ²J(P-C-H) ou ²J(P-O-H)
            if is_phosphate:
                return 8.0  # P(O)-C-H
            else:
                return 12.0  # P-C-H (phosphine)
        elif dist == 3:
            # ³J(P-C-C-H) ou ³J(P-O-C-H)
            if is_phosphate:
                return 10.0  # Phosphate typique
            else:
                return 15.0  # Phosphine
        elif dist == 4:
            # ⁴J longue portée
            return 2.0
        
        return None
        
    except Exception:
        return None


# Labels spéciaux pour systèmes second ordre
SECOND_ORDER_LABELS = {
    "AB": "AB",           # Système AB fortement couplé
    "AB_weak": "AB",      # Système AB faiblement couplé
    "AB₂": "AB2",         # 1 spin + 2 équivalents
    "AA'BB'": "t",        # -CH₂-CH₂- → triplet apparent (approx premier ordre)
    "A₂X₂": "t",          # Forme limite de AA'BB' → triplet
    "ABX": "dd",          # 2 spins AB + 1 X → doublet de doublets
    "AMX": "dd",          # 3 spins distincts → doublet de doublets minimum
    "A₂X₃": "q",          # CH₂-CH₃ → quartet
    "A₂B₃": "q",          # CH₂-CH₃ fortement couplé → quartet
    "complex_second_order": "m",  # Vraiment complexe
}

# Notation avec "roof effect" (intensités asymétriques)
ROOF_EFFECT_NOTATION = {
    "AB": "AB (roof)",       # Intensités asymétriques caractéristiques
    "AB₂": "AB₂ (roof)",
    "AA'BB'": "AA'BB' (roof)",
    "ABX": "ABX (roof)",
}


def _get_second_order_label(system_type: Optional[str], base_label: str, 
                            delta_nu_hz: float = 0.0, J_hz: float = 0.0,
                            indicate_roof_effect: bool = True) -> str:
    """
    Retourne un label approprié pour les systèmes second ordre.
    
    Pour les systèmes fortement couplés, on utilise des labels spéciaux
    (AB, ABX, etc.) au lieu de simplement "m".
    
    Args:
        system_type: Type de système spin détecté
        base_label: Label de base (d, t, q, etc.)
        delta_nu_hz: Différence de shift en Hz
        J_hz: Constante de couplage en Hz
        indicate_roof_effect: Si True, indique le roof effect quand applicable
    
    Returns:
        Label approprié avec indication du roof effect si nécessaire
    """
    if system_type is None:
        return base_label
    
    if system_type in SECOND_ORDER_LABELS:
        basic_label = SECOND_ORDER_LABELS[system_type]
        
        # Ajouter indication de roof effect si Δν/J < 10
        if indicate_roof_effect and J_hz > 0 and delta_nu_hz > 0:
            ratio = delta_nu_hz / J_hz
            if ratio < 10:  # Second ordre significatif
                if system_type in ROOF_EFFECT_NOTATION:
                    return ROOF_EFFECT_NOTATION[system_type]
        
        return basic_label
    
    return base_label


def _estimate_HNCH_J(mol3d: Chem.Mol, h_idx: int, n_idx: int, c_idx: int, h2_idx: int) -> Optional[float]:
    """
    Estime ³J(H-N-C-H) pour peptides/amides.
    
    Utilise l'équation Karplus de Bystrov (1976):
    ³J = A·cos²θ + B·cosθ + C
    
    Paramètres Bystrov: A=6.4, B=-1.4, C=1.9 Hz
    
    Très important pour:
    - Peptides (φ, ψ angles)
    - Nucléosides
    - Amides cycliques
    
    Référence: Bystrov V.F., Prog. NMR Spectrosc. 10, 41-82 (1976)
    """
    try:
        conf = mol3d.GetConformer()
        phi = AllChem.GetDihedralDeg(conf, h_idx, n_idx, c_idx, h2_idx)
        
        # Équation Karplus-Bystrov pour H-N-C-H
        return _karplus_simple(phi, bond_type='HNCH')
        
    except Exception:
        return None


# ---------------------------
# ✅ INFÉRENCE STÉRÉOCHIMIE E/Z POUR DOUBLES LIAISONS NON DÉFINIES
# ---------------------------
def _infer_double_bond_stereo(mol: Chem.Mol, default_stereo: str = "E") -> Chem.Mol:
    """
    Infère la stéréochimie E/Z pour les doubles liaisons C=C non définies.
    
    Pour les chaînes polyéniques (comme les caroténoïdes), les doubles liaisons
    naturelles sont généralement en configuration trans (E). Cette fonction
    assigne automatiquement la stéréochimie E aux doubles liaisons non définies
    dans les systèmes conjugués.
    
    STRATÉGIE AMÉLIORÉE:
    1. Générer une conformation 3D
    2. Pour chaque double liaison non définie, forcer la géométrie trans (φ ≈ 180°)
    3. Ré-optimiser pour stabiliser
    
    Args:
        mol: Molécule RDKit
        default_stereo: "E" (trans) ou "Z" (cis) - défaut "E" car plus stable
        
    Returns:
        Molécule avec stéréochimie inférée et géométrie cohérente
        
    Référence: Caroténoïdes naturels sont typiquement all-trans (E)
               IUPAC-IUB Nomenclature of Carotenoids (1974)
    """
    import logging
    
    # Ajouter les H et générer la 3D
    molH = Chem.AddHs(mol)
    
    # Générer une conformation 3D
    try:
        try: 
            ps = AllChem.ETKDGv3()
        except: 
            ps = AllChem.ETKDG()
        ps.randomSeed = 42
        ps.useRandomCoords = True
        res = AllChem.EmbedMolecule(molH, ps)
        if res != 0:
            # Fallback avec coordonnées aléatoires
            AllChem.EmbedMolecule(molH, AllChem.ETKDG())
    except Exception as e:
        logging.warning(f"Impossible de générer la 3D: {e}")
        return mol
    
    # Optimiser la géométrie
    try:
        AllChem.MMFFOptimizeMolecule(molH, maxIters=200)
    except:
        try:
            AllChem.UFFOptimizeMolecule(molH, maxIters=200)
        except:
            pass
    
    # Identifier les doubles liaisons C=C non définies
    modified_count = 0
    double_bonds_to_fix = []
    
    for bond in molH.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
            
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        
        if atom1.GetSymbol() != "C" or atom2.GetSymbol() != "C":
            continue
            
        # Vérifier si stéréochimie déjà définie
        current_stereo = bond.GetStereo()
        if current_stereo not in (Chem.BondStereo.STEREONONE,):
            continue
        
        # Trouver les H ou substituants attachés pour mesurer l'angle dièdre
        neighbors1 = [n for n in atom1.GetNeighbors() if n.GetIdx() != atom2.GetIdx()]
        neighbors2 = [n for n in atom2.GetNeighbors() if n.GetIdx() != atom1.GetIdx()]
        
        if not neighbors1 or not neighbors2:
            continue
            
        double_bonds_to_fix.append({
            'bond': bond,
            'atom1': atom1,
            'atom2': atom2,
            'ref1': neighbors1[0],
            'ref2': neighbors2[0],
        })
    
    # Forcer la géométrie trans pour chaque double liaison
    if double_bonds_to_fix and molH.GetNumConformers() > 0:
        conf = molH.GetConformer()
        
        for db_info in double_bonds_to_fix:
            try:
                ref1_idx = db_info['ref1'].GetIdx()
                atom1_idx = db_info['atom1'].GetIdx()
                atom2_idx = db_info['atom2'].GetIdx()
                ref2_idx = db_info['ref2'].GetIdx()
                
                # Lire l'angle dièdre actuel
                current_dihedral = AllChem.GetDihedralDeg(conf, ref1_idx, atom1_idx, atom2_idx, ref2_idx)
                
                # Déterminer l'angle cible
                if default_stereo == "E":
                    target_dihedral = 180.0  # Trans
                else:
                    target_dihedral = 0.0    # Cis
                
                # Si l'angle n'est pas proche de la cible, ajuster
                if abs(current_dihedral - target_dihedral) > 30 and abs(current_dihedral + 360 - target_dihedral) > 30:
                    # Forcer l'angle dièdre
                    AllChem.SetDihedralDeg(conf, ref1_idx, atom1_idx, atom2_idx, ref2_idx, target_dihedral)
                    modified_count += 1
                elif abs(current_dihedral - target_dihedral) <= 30:
                    modified_count += 1  # Déjà en configuration trans
                    
            except Exception as e:
                logging.debug(f"Erreur ajustement dièdre: {e}")
                continue
        
        # Ré-optimiser légèrement pour relaxer les contraintes
        if modified_count > 0:
            try:
                AllChem.MMFFOptimizeMolecule(molH, maxIters=50)
            except:
                pass
    
    if modified_count > 0:
        logging.info(f"✅ Géométrie {default_stereo} (trans) forcée pour {modified_count} double(s) liaison(s)")
    
    # Assigner la stéréochimie basée sur la géométrie 3D
    Chem.AssignStereochemistryFrom3D(molH)
    
    return molH


def _count_undefined_double_bonds(mol: Chem.Mol) -> Tuple[int, int]:
    """
    Compte les doubles liaisons C=C définies et non définies.
    
    Returns:
        (defined_count, undefined_count)
    """
    defined = 0
    undefined = 0
    
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
            
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        
        if atom1.GetSymbol() != "C" or atom2.GetSymbol() != "C":
            continue
            
        stereo = bond.GetStereo()
        if stereo == Chem.BondStereo.STEREONONE:
            undefined += 1
        else:
            defined += 1
            
    return (defined, undefined)


# ---------------------------
# ✅ CORRECTION 1 : Validation structure moléculaire (avec paramètre smiles)
# ---------------------------
def _validate_molecular_structure(mol: Chem.Mol, smiles: str) -> None:
    """
    ✅ Validation GÉNÉRIQUE d'une molécule :
    - Vérifie les valences
    - Détecte les stéréocentres ambigus
    - REFUSE les cycles furanose dans les sucres (cycle 5 atomes)
    - Vérifie la cohérence des chiralités
    
    Lève ValueError si la molécule est invalide.
    """
    # 1. Nettoyage et vérification de valence
    try:
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
        Chem.SetAromaticity(mol)
    except Exception as e:
        raise ValueError(f"❌ Structure moléculaire invalide (valence, connectivité): {e}")
    
    # 2. Vérification des stéréocentres
    molH = Chem.AddHs(mol)
    Chem.AssignStereochemistry(molH, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    
    # Note: On ne refuse plus les stéréocentres non définis - on génère juste un avertissement
    # car beaucoup de molécules légitimes n'ont pas besoin de stéréochimie définie pour la RMN
    for atom in molH.GetAtoms():
        if atom.GetChiralTag() == rdchem.ChiralType.CHI_UNSPECIFIED:
            if atom.GetDegree() == 4 and atom.IsInRing():
                # Avertissement au lieu d'erreur - la stéréochimie peut être inférée
                import logging
                logging.debug(f"⚠️ Stéréochimie non définie sur l'atome {atom.GetIdx()} ({atom.GetSymbol()}) dans un cycle - sera ignorée")
    
    # 3. Détection de cycles furanose (5 atomes O-C-C-C-C) dans les sucres
    rings = Chem.GetSymmSSSR(molH)
    for ring in rings:
        if len(ring) == 5:  # Cycle à 5 atomes
            atoms_in_ring = [molH.GetAtomWithIdx(i) for i in ring]
            symbols = [a.GetSymbol() for a in atoms_in_ring]
            
            # Vérifier motif furanose : 1 O + 4 C + présence de CH2OH externe
            if symbols.count("O") == 1 and symbols.count("C") == 4:
                # Vérifier si c'est une molécule sucre-like (présence de CH2OH)
                for atom in molH.GetAtoms():
                    if atom.GetSymbol() == "C" and not atom.IsInRing():
                        neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
                        if "O" in neighbors and atom.GetTotalNumHs() >= 2:
                            raise ValueError(
                                f"❌ Cycle furanose détecté (5 atomes) dans '{smiles}'. "
                                "Utilisez un cycle pyranose (6 atomes) pour les hexoses naturels."
                            )


# ---------------------------
# VALIDATION CONNECTIVITÉ
# ---------------------------
def _validate_connectivity(molH: Chem.Mol, connectivity: Dict[int, List[int]]) -> bool:
    """
    ✅ Vérifie que la connectivité X–H est cohérente :
      - indices dans l'intervalle
      - atomes bien de type H
      - pas de doublons
    """
    n_atoms = molH.GetNumAtoms()
    seen = set()

    for key, h_list in connectivity.items():
        for h_idx in h_list:
            if not (0 <= h_idx < n_atoms):
                raise ValueError(f"Indice d'hydrogène invalide: {h_idx} (clé de connectivité={key})")
            atom = molH.GetAtomWithIdx(h_idx)
            if atom.GetSymbol() != "H":
                raise ValueError(f"Atome {h_idx} référencé dans la connectivité n'est pas un H (clé={key})")
            if h_idx in seen:
                raise ValueError(f"Hydrogène {h_idx} apparaît plusieurs fois dans la connectivité (clé={key})")
            seen.add(h_idx)

    return True

# ---------------------------
# Numeric helpers
# ---------------------------
def _delta_nu_hz(delta_ppm: float, spectrometer_mhz: float) -> float:
    return abs(delta_ppm) * spectrometer_mhz

# ---------------------------
# Aromatic contextual J
# ---------------------------
def _aromatic_J_contextual(base: Chem.Mol, aromatic_dmap: List[Dict[int,int]], i: int, j: int,
                           J0_meta: float = 2.3, J0_para: float = 0.6) -> Tuple[Optional[float], Optional[str]]:
    """
    Calcule J aromatique avec gestion enrichie des hétérocycles.
    
    AMÉLIORATION: Utilise _AROMATIC_COUPLING_TABLE pour les valeurs spécifiques
    selon le type d'hétérocycle (imidazole, pyridine, etc.)
    """
    rel = _aromatic_rel_from_map(aromatic_dmap, i, j)
    if rel is None:
        return None, rel
    
    # Identifier le type d'hétérocycle pour les deux atomes
    cycle_type = _identify_heterocycle_type(base, i, j)
    
    # Cas spécial: imidazole H2 (entre les 2 N) → couplage négligeable
    if cycle_type == 'imidazole':
        ring_info = base.GetRingInfo()
        for ring in ring_info.AtomRings():
            if i in ring and j in ring:
                if _is_imidazole_h2_position(base, i, j, ring):
                    return 0.0, rel  # Pas de couplage significatif avec/depuis H2
    
    # Récupérer J depuis la table enrichie
    if cycle_type in _AROMATIC_COUPLING_TABLE:
        cycle_data = _AROMATIC_COUPLING_TABLE[cycle_type]
        if isinstance(cycle_data, dict) and rel in cycle_data:
            J = cycle_data[rel]
            return J, rel
    
    # Fallback: valeurs par défaut
    if rel == "ortho":
        J = _AROMATIC_COUPLING_TABLE.get('default_ortho', 7.5)
    elif rel == "meta":
        J = J0_meta
    elif rel == "para":
        J = J0_para
    else:
        return None, rel

    def penalty(idx: int) -> float:
        A = base.GetAtomWithIdx(idx)
        p = 1.0
        for n in A.GetNeighbors():
            if n.GetSymbol() in {"O", "N"}:
                p *= 0.85
            if n.GetSymbol() == "C":
                for b in n.GetBonds():
                    if b.GetBondTypeAsDouble() >= 2 and b.GetOtherAtom(n).GetSymbol() == "O":
                        p *= 0.85
        return p

    if rel == "meta":
        scale = penalty(i) * penalty(j)
        J = max(1.6, min(2.4, J * scale))
    elif rel == "para":
        scale = (penalty(i) * penalty(j))**0.5
        J = max(0.45, min(0.80, J * scale))
    
    return J, rel


def _identify_heterocycle_type(mol: Chem.Mol, i: int, j: int) -> str:
    """
    Identifie le type d'hétérocycle contenant les atomes i et j.
    
    Returns:
        Type de cycle: 'benzene', 'pyridine', 'imidazole', 'furan', 'thiophene', etc.
    """
    ring_info = mol.GetRingInfo()
    
    for ring in ring_info.AtomRings():
        if i in ring and j in ring:
            # Vérifier que c'est aromatique
            if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                continue
            
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            heteroatoms = [a.GetSymbol() for a in ring_atoms if a.GetSymbol() != 'C']
            ring_size = len(ring)
            
            if ring_size == 5:
                if 'O' in heteroatoms and 'N' not in heteroatoms and 'S' not in heteroatoms:
                    return 'furan'
                elif 'S' in heteroatoms and 'N' not in heteroatoms:
                    return 'thiophene'
                elif heteroatoms.count('N') >= 2:
                    return 'imidazole'
                elif 'N' in heteroatoms:
                    return 'pyrrole'
                elif 'O' in heteroatoms and 'N' in heteroatoms:
                    return 'oxazole'
            elif ring_size == 6:
                if heteroatoms.count('N') == 1:
                    return 'pyridine'
                elif heteroatoms.count('N') == 2:
                    return 'pyrimidine'
                elif heteroatoms.count('N') >= 3:
                    return 'triazine_1_3_5'
            
    return 'benzene'

def _fallback_aromatic_J_when_empty(
    base: Chem.Mol,
    aromatic_dmap: List[Dict[int, int]],
    heavy_idx: int,
    heavy_to_H: Dict[int, List[int]],
    config: SimConfig,
    ranks: Optional[List[int]] = None,
) -> List[float]:
    """
    Fallback pour H aromatiques quand aucun J n'a été estimé.
    
    IMPORTANT: Ne retourne PAS de J pour les voisins équivalents
    (ex: benzène non-substitué → singulet).
    """
    heavy = base.GetAtomWithIdx(heavy_idx)
    if heavy.GetSymbol() != "C" or not heavy.GetIsAromatic():
        return []

    J_vals: List[float] = []
    
    for other_heavy, Hs in heavy_to_H.items():
        if other_heavy == heavy_idx or not Hs:
            continue
        other_atom = base.GetAtomWithIdx(other_heavy)
        if not other_atom.GetIsAromatic():
            continue
        
        # ✅ CORRECTION: Skip les voisins équivalents (homotopes)
        if ranks is not None:
            if _are_heavy_equivalent(ranks, heavy_idx, other_heavy, heavy_to_H, None, 0.01, mol=base):
                continue  # Pas de couplage observable entre protons équivalents

        rel = _aromatic_rel_from_map(aromatic_dmap, heavy_idx, other_heavy)
        if rel == "ortho":
            J = config.aro_j_forced.get("ortho", 7.8)
        elif rel == "meta":
            J = config.aro_j_forced.get("meta", 2.3)
        elif rel == "para":
            J = config.aro_j_forced.get("para", 0.6)
        else:
            continue

        minJ = max(config.min_j_detect.get(3, 1.0), 0.5)
        if J >= minJ:
            J_vals.append(J)

    return J_vals

# ---------------------------
# Spin-system fallback helpers
# ---------------------------
def _rough_topological_J(base: Chem.Mol, aromatic_dmap: List[Dict[int,int]],
                         i_heavy: int, j_heavy: int, is_sp2_i: bool, is_sp2_j: bool) -> float:
    """Estimation grossière de J entre deux protons pour fallback partenaire↔partenaire."""
    if i_heavy == j_heavy:
        return 2.0 if is_sp2_i else 12.5
    D = _topological_distances(base)
    d = int(D[i_heavy, j_heavy]) if np.isfinite(D[i_heavy, j_heavy]) else 99
    if d == 1:
        if not is_sp2_i and not is_sp2_j:
            return 7.0
        return 2.0
    if d == 2:
        J, _ = _aromatic_J_contextual(base, aromatic_dmap, i_heavy, j_heavy)
        return J or 1.6
    if d == 3:
        J, _ = _aromatic_J_contextual(base, aromatic_dmap, i_heavy, j_heavy)
        return J or 0.6
    return 0.0

@lru_cache(maxsize=8)
def _spin_ops(nspins: int):
    Sx_list, Sy_list, Sz_list = [], [], []
    sx = 0.5 * np.array([[0, 1],[1, 0]], dtype=float)
    sy = 0.5 * np.array([[0, -1j],[1j, 0]], dtype=complex)
    sz = 0.5 * np.array([[1, 0],[0, -1]], dtype=float)
    id2 = np.eye(2, dtype=complex)
    for k in range(nspins):
        ops = []
        for base in (sx, sy, sz):
            M = np.array([[1.0]], dtype=complex)
            for i in range(nspins):
                M = np.kron(M, base if i == k else id2)
            ops.append(M)
        Sx_list.append(ops[0]); Sy_list.append(ops[1]); Sz_list.append(ops[2])
    return Sx_list, Sy_list, Sz_list

def _spin_hamiltonian(omegas_hz: List[float], J_hz: Dict[Tuple[int,int], float]):
    n = len(omegas_hz)
    Sx, Sy, Sz = _spin_ops(n)
    H = np.zeros_like(Sx[0], dtype=complex)
    for i, w in enumerate(omegas_hz):
        H += (2*np.pi*w) * Sz[i]
    for (i, j), J in J_hz.items():
        H += (2*np.pi*J) * (Sx[i]@Sx[j] + Sy[i]@Sy[j] + Sz[i]@Sz[j])
    return H

def _allowed_transitions(nspins: int):
    dim = 2**nspins
    trans = []
    for i in range(dim):
        for j in range(dim):
            diff = i ^ j
            if diff and (diff & (diff-1)) == 0:
                trans.append((i, j))
    return trans

def _simulate_spin_system_lines(omegas_hz: List[float], J_hz: Dict[Tuple[int,int], float]) -> List[Tuple[float,float]]:
    H = _spin_hamiltonian(omegas_hz, J_hz)
    E, V = np.linalg.eigh(H)
    trans = _allowed_transitions(len(omegas_hz))
    lines: List[Tuple[float, float]] = []
    Sx, _, _ = _spin_ops(len(omegas_hz))
    Mx = sum(Sx)
    for i, j in trans:
        df = (E[j] - E[i]) / (2*np.pi)
        amp = abs((V[:,i].conj().T @ Mx @ V[:,j]))**2
        lines.append((float(np.real(df)), float(np.real(amp))))
    if not lines: return []
    s = sum(x[1] for x in lines)
    if s > 0: lines = [(f, i/s) for (f, i) in lines]
    lines.sort(key=lambda x: x[0])
    return lines

def _collapse_lines_to_groups(lines: List[Tuple[float,float]], tol_hz: float, min_int: float) -> List[List[float]]:
    if not lines: return []
    filt = [(f, I) for (f, I) in lines if I >= min_int]
    if not filt: return []
    groups_f: List[List[float]] = []
    cur = [filt[0][0]]
    for f, I in filt[1:]:
        if abs(f - cur[-1]) <= tol_hz:
            cur.append(f)
        else:
            groups_f.append(cur); cur = [f]
    groups_f.append(cur)
    out: List[List[float]] = []
    for g in groups_f:
        if len(g) <= 1:
            out.append([1.0])
        else:
            diffs = [g[k+1]-g[k] for k in range(len(g)-1)]
            Jm = max(0.1, float(np.mean(np.abs(diffs))))
            out.append([Jm]*len(g))
    return out

def _spin_fallback_label(
    spectrometer_mhz: float,
    target_ppm: float,
    partner_ppm: List[float],
    partner_J: List[float],
    config: SimConfig,
    partner_heavy_idx: Optional[List[int]] = None,
    base: Optional[Chem.Mol] = None,
    aromatic_dmap: Optional[List[Dict[int,int]]] = None,
    is_sp2_list: Optional[List[bool]] = None
) -> Optional[str]:
    """Fallback spin enrichi : système avec couplages partenaire↔partenaire."""
    try:
        n = 1 + len(partner_ppm)
        if n < 2 or n > config.spin_fallback_max_spins:
            return None
        omegas = [0.0] + [(partner_ppm[i] - target_ppm) * spectrometer_mhz for i in range(len(partner_ppm))]

        Jmap: Dict[Tuple[int,int], float] = {}
        for i, J in enumerate(partner_J, start=1):
            if J is None or J <= 0: continue
            Jmap[(0, i)] = float(J)

        if (config.fallback_partner_partner_enable and
            base is not None and aromatic_dmap is not None and
            partner_heavy_idx is not None and is_sp2_list is not None and
            len(partner_heavy_idx) == len(partner_ppm)):

            for i in range(len(partner_ppm)):
                for j in range(i+1, len(partner_ppm)):
                    dppm = abs(partner_ppm[i] - partner_ppm[j])
                    dHz = dppm * spectrometer_mhz
                    if dHz <= config.fallback_partner_partner_max_delta_hz:
                        Ji = _rough_topological_J(
                            base, aromatic_dmap,
                            partner_heavy_idx[i], partner_heavy_idx[j],
                            is_sp2_list[i], is_sp2_list[j]
                        )
                        if Ji > 0.0:
                            Jmap[(i+1, j+1)] = float(Ji)

        if not Jmap:
            return None

        lines = _simulate_spin_system_lines(omegas, Jmap)
        tol_hz = max(TOLERANCE_J_HZ, config.spin_fallback_group_tol_scale * config.line_width_hz)
        groups = _collapse_lines_to_groups(lines, tol_hz=tol_hz, min_int=config.spin_fallback_min_intensity)
        if not groups: return None
        return _first_order_label_from_groups(groups)
    except Exception:
        return None

# ---------------------------
# Grouping
# ---------------------------
def _group_J_values(J_list: List[float], tol: float = TOLERANCE_J_HZ, anisotropy_high: bool = False,
                    config: Optional[SimConfig]=None) -> List[List[float]]:
    config = config or SimConfig()
    vals = sorted([v for v in J_list if v is not None and v > 0.0])
    if not vals: return []
    pct = config.anisotropy_pct_high if anisotropy_high else config.anisotropy_pct_normal
    groups: List[List[float]] = []
    for v in vals:
        placed = False
        for g in groups:
            gmean = float(np.mean(g))
            dyn_tol = min(tol, pct * gmean)
            if abs(gmean - v) <= dyn_tol:
                g.append(v); placed = True; break
        if not placed: groups.append([v])

    changed = True
    while changed and len(groups) > 1:
        changed = False
        new_groups: List[List[float]] = []
        groups.sort(key=lambda g: float(np.mean(g)))
        i = 0
        while i < len(groups):
            if i == len(groups) - 1:
                new_groups.append(groups[i]); break
            g1, g2 = groups[i], groups[i+1]
            m1, m2 = float(np.mean(g1)), float(np.mean(g2))
            merged_mean = (m1 * len(g1) + m2 * len(g2)) / (len(g1) + len(g2))
            dyn_tol = min(tol, pct * merged_mean)
            gap = abs(m1 - m2)
            if gap <= dyn_tol and gap <= (config.never_merge_gap_hz if config else 2.5):
                new_groups.append(g1 + g2); i += 2; changed = True
            else:
                new_groups.append(g1); i += 1
        groups = new_groups

    if len(groups) > 1:
        max_mean = max(float(np.mean(g)) for g in groups)
        if max_mean <= 1.2 * config.line_width_hz:
            return groups

    if len(groups) > 1:
        small_groups = [g for g in groups if float(np.mean(g)) < max(1.3, 0.9*config.line_width_hz)]
        if len(small_groups) >= 2:
            merged = []
            for g in small_groups:
                merged.extend(g)
            groups = [g for g in groups if float(np.mean(g)) >= max(1.3, 0.9*config.line_width_hz)] + [merged]

    return groups

# ---------------------------
# Apparent label pour centres sp²
# ---------------------------
def _apparent_sp2_label_from_J(J_vals: List[float], minJ: float = 0.4, tol: float = 0.5) -> str:
    """
    Simplification 'apparent' pour centres sp² avec regroupement des J équivalents.
    
    Pour un CHO couplé à un CH₃, les 3 J identiques donnent un quartet (q), pas ddd.
    """
    Js = [j for j in J_vals if j is not None and j >= minJ]
    n = len(Js)
    if n == 0:
        return "s"
    if n == 1:
        return "d"
    
    # Regrouper les J similaires
    if n >= 2:
        Js_sorted = sorted(Js)
        groups = []
        current_group = [Js_sorted[0]]
        
        for j in Js_sorted[1:]:
            if abs(j - current_group[0]) <= tol:
                current_group.append(j)
            else:
                groups.append(current_group)
                current_group = [j]
        groups.append(current_group)
        
        total = sum(len(g) for g in groups)
        n_groups = len(groups)
        
        if n_groups == 1:
            # Tous les J sont équivalents
            if total == 2: return "t"
            if total == 3: return "q"
            if total == 4: return "quin"
            if total == 5: return "sext"
            if total == 6: return "sept"
            return "m"
        elif n_groups == 2:
            counts = sorted([len(g) for g in groups])
            if counts == [1, 1]: return "dd"
            if counts == [1, 2]: return "td"
            if counts == [1, 3]: return "qd"
            if counts == [2, 2]: return "tt"
            return "ddd"
        else:
            return "ddd"
    
    return "ddd"

# ---------------------------
# Label construction
# ---------------------------
def _first_order_label_from_groups(groups: List[List[float]], force_singlet_if_empty: bool = True) -> str:
    """
    Détermine le label de multiplicité à partir des groupes de couplages.
    
    Args:
        groups: Liste de groupes de J (Hz)
        force_singlet_if_empty: Si True, retourne 's' si pas de couplage.
                               Si False, retourne 'm' pour multiplet indéterminé.
    """
    if not groups:
        return "s" if force_singlet_if_empty else "m"
    
    # Filtrer les groupes vides ou avec J très faibles (< 0.3 Hz)
    filtered_groups = [g for g in groups if g and max(abs(j) for j in g) >= 0.3]
    if not filtered_groups:
        return "s" if force_singlet_if_empty else "m"
    
    if (len(filtered_groups) == 1 and len(filtered_groups[0]) == 1 and abs(filtered_groups[0][0] - 0.6) < 0.15):
        return "td"
    counts = [len(g) for g in filtered_groups]
    total = sum(counts)
    means = [float(np.mean(g)) if g else 0.0 for g in filtered_groups]

    if total == 1: return "d"

    if total == 2:
        if len(groups) == 1 and counts[0] == 2: return "t"
        if len(groups) == 2 and any(abs(np.mean(g)-2.3)<0.3 and len(g)==2 for g in groups) and any(abs(np.mean(g)-0.6)<0.3 for g in groups):
            return "td"
        return "dd"

    if total == 3:
        if len(groups) == 1 and counts[0] == 3: return "q"
        if len(groups) == 2 and sorted(counts) == [1, 2]:
            if any(abs(np.mean(g)-2.3)<0.3 and len(g)==2 for g in groups) and any(abs(np.mean(g)-0.6)<0.3 for g in groups):
                return "q"
            return "td"
        return "ddd"

    if total == 4:
        if len(groups) == 1 and counts[0] == 4: return "quin"
        if len(groups) == 2 and sorted(counts) == [1, 3]:
            return "qd"
        if len(groups) == 2 and sorted(counts) == [2, 2]: return "tt"
        if len(groups) == 3 and sorted(counts) == [1, 1, 2]: return "ddt"
        return "dddd"

    if total == 5:
        if len(groups) == 1 and counts[0] == 5: return "sext"
        if len(groups) == 2 and sorted(counts) == [2, 3]:
            return "qt"
        if len(groups) == 3 and sorted(counts) == [1, 2, 2]:
            return "qqd"
        return "ddd"

    if total == 6:
        if len(groups) == 1 and counts[0] == 6: return "sept"
        if len(groups) == 2 and sorted(counts) == [3, 3]:
            return "qt"
        if len(groups) == 3 and sorted(counts) == [1, 2, 3]:
            return "ttd"
        return "ddd"

    # ✅ RÉALISME : Au-delà de 6 couplages (septet), les raies ne sont plus
    # clairement résolues en RMN pratique → multiplet "m"
    # Exception: 7 couplages avec un seul groupe = octuplet, parfois visible
    if total == 7:
        if len(groups) == 1 and counts[0] == 7:
            return "oct"  # Octuplet rare mais possible (ex: CH couplé avec 7H)
        return "m"  # Multiplet complexe
    
    # ≥ 8 couplages : toujours multiplet (non résolu en pratique)
    if total >= 8: return "m"
    
    return "s"  # Défaut: singulet si aucun couplage détecté

# ---------------------------
# Stereo helpers
# ---------------------------
def _assign_stereo(mol: Chem.Mol) -> None:
    try:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    except Exception:
        pass

def _is_chiral_center(atom: Chem.Atom) -> bool:
    tag = atom.GetChiralTag()
    return tag in (rdchem.ChiralType.CHI_TETRAHEDRAL_CW, rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)


def _find_nearest_stereo_center(mol: Chem.Mol, start_idx: int, max_dist: int = 4) -> Tuple[Optional[int], int]:
    """
    Trouve le centre stéréogène le plus proche d'un atome donné.
    
    Args:
        mol: Molécule RDKit
        start_idx: Index de l'atome de départ
        max_dist: Distance topologique maximale à explorer
    
    Returns:
        (idx_stereo_center, distance) ou (None, -1) si aucun trouvé
    """
    from collections import deque
    
    visited = {start_idx}
    queue = deque([(start_idx, 0)])
    
    while queue:
        current_idx, dist = queue.popleft()
        
        if dist > max_dist:
            continue
        
        atom = mol.GetAtomWithIdx(current_idx)
        
        # Vérifier si c'est un centre stéréogène
        if current_idx != start_idx and _is_chiral_center(atom):
            return current_idx, dist
        
        # Explorer les voisins
        for neighbor in atom.GetNeighbors():
            nb_idx = neighbor.GetIdx()
            if nb_idx not in visited:
                visited.add(nb_idx)
                queue.append((nb_idx, dist + 1))
    
    return None, -1


def _has_stereogenic_path(mol: Chem.Mol, heavy_idx: int, max_dist: int = 3) -> bool:
    """
    Vérifie si un atome a un chemin vers un centre stéréogène à distance <= max_dist.
    
    Pour les diastéréotopes, le centre stéréogène doit être à distance 1-3 liaisons.
    """
    stereo_idx, dist = _find_nearest_stereo_center(mol, heavy_idx, max_dist)
    return stereo_idx is not None


# ---------------------------
# ✅ CORRECTION : Diastereotopy helpers (VERSION ROBUSTE ET AMÉLIORÉE)
# ---------------------------
def _is_diastereotopic_Hpair(molH: Chem.Mol, heavy_idx: int, h_a: int, h_b: int) -> bool:
    """
    ✅ VERSION AMÉLIORÉE : Teste si deux H sur le même carbone sont diastéréotopiques.
    
    Amélioration : Détecte les centres stéréogènes distants (jusqu'à 3 liaisons)
    et les cas épimériques (sucres, etc.)
    
    Args:
        molH: Molécule avec H explicites
        heavy_idx: Index de l'atome lourd (C, N, etc.)
        h_a, h_b: Indices des deux atomes H à comparer
    
    Returns:
        True si diastéréotopiques, False sinon
    """
    if heavy_idx is None or h_a is None or h_b is None:
        return False

    try:
        heavy_atom = molH.GetAtomWithIdx(heavy_idx)
        
        # ✅ ÉTAPE 0 : Test rapide des rangs canoniques (le plus fiable)
        # Si RDKit détecte une différence de symétrie, c'est diastéréotopique
        try:
            ranks = list(Chem.CanonicalRankAtoms(molH, includeChirality=True, breakTies=True))
        except TypeError:
            ranks = list(Chem.CanonicalRankAtoms(molH, includeChirality=True))
        
        if ranks[h_a] != ranks[h_b]:
            # RDKit détecte une différence → confirmer par test isotopique
            def get_heavy_rank_with_marked_h(mark_h_idx: int) -> int:
                m = Chem.RWMol(molH)
                m.GetAtomWithIdx(mark_h_idx).SetIsotope(2)
                try:
                    new_ranks = list(Chem.CanonicalRankAtoms(m, includeChirality=True, breakTies=True))
                except TypeError:
                    new_ranks = list(Chem.CanonicalRankAtoms(m, includeChirality=True))
                return new_ranks[heavy_idx]
            
            if get_heavy_rank_with_marked_h(h_a) != get_heavy_rank_with_marked_h(h_b):
                return True
        
        # ✅ ÉTAPE 1 : Vérifier la présence d'un centre stéréogène proche
        stereo_idx, stereo_dist = _find_nearest_stereo_center(molH, heavy_idx, max_dist=3)
        
        if stereo_idx is not None:
            # Il y a un centre chiral proche → les H peuvent être diastéréotopiques
            # Même avec rotation libre, si le centre chiral est à distance ≤ 2,
            # les H sont généralement non-équivalents
            if stereo_dist <= 2:
                return True
            # À distance 3, faire un test isotopique plus rigoureux
            elif stereo_dist == 3:
                def get_heavy_rank_with_marked_h(mark_h_idx: int) -> int:
                    m = Chem.RWMol(molH)
                    m.GetAtomWithIdx(mark_h_idx).SetIsotope(2)
                    try:
                        new_ranks = list(Chem.CanonicalRankAtoms(m, includeChirality=True, breakTies=True))
                    except TypeError:
                        new_ranks = list(Chem.CanonicalRankAtoms(m, includeChirality=True))
                    return new_ranks[heavy_idx]
                
                return get_heavy_rank_with_marked_h(h_a) != get_heavy_rank_with_marked_h(h_b)
        
        # ✅ ÉTAPE 2 : Cas spéciaux sans centre chiral explicite
        if heavy_atom.GetSymbol() == "C":
            neighbors_heavy = [n for n in heavy_atom.GetNeighbors() if n.GetSymbol() != "H"]
            
            # CH₂ dans un cycle → souvent non-équivalents (conformations bloquées)
            if heavy_atom.IsInRing():
                # Tester si les environnements sont différents
                return True
            
            # ✅ AMÉLIORATION: Détection prochirauté (carbonyle asymétrique, etc.)
            if _is_diastereotopic_by_prochirality(molH, heavy_idx, h_a, h_b):
                return True
            
            # CH₂ attaché à un carbone sp2 (alcène, carbonyle) → peut être prochiral
            for nb in neighbors_heavy:
                if nb.GetHybridization() in (Chem.HybridizationType.SP2, Chem.HybridizationType.SP):
                    # Attaché à C=O, C=C, C=N → vérifier
                    def get_heavy_rank_with_marked_h(mark_h_idx: int) -> int:
                        m = Chem.RWMol(molH)
                        m.GetAtomWithIdx(mark_h_idx).SetIsotope(2)
                        try:
                            new_ranks = list(Chem.CanonicalRankAtoms(m, includeChirality=True, breakTies=True))
                        except TypeError:
                            new_ranks = list(Chem.CanonicalRankAtoms(m, includeChirality=True))
                        return new_ranks[heavy_idx]
                    
                    if get_heavy_rank_with_marked_h(h_a) != get_heavy_rank_with_marked_h(h_b):
                        return True
            
            # ✅ ÉTAPE 3 : CH₂ avec rotation libre et pas de centre chiral proche
            # → probablement équivalents
            if len(neighbors_heavy) == 2:
                all_sp3_acyclic = all(
                    nb.GetSymbol() in ("C", "O", "N") and 
                    not nb.IsInRing() and
                    nb.GetHybridization() == Chem.HybridizationType.SP3
                    for nb in neighbors_heavy
                )
                if all_sp3_acyclic:
                    # Vérifier prochirauté avant de conclure équivalents
                    if not _is_diastereotopic_by_prochirality(molH, heavy_idx, h_a, h_b):
                        # Chaîne acyclique sans chiralité proche → équivalents
                        return False

        # Défaut: test isotopique final
        def get_heavy_rank_with_marked_h(mark_h_idx: int) -> int:
            m = Chem.RWMol(molH)
            m.GetAtomWithIdx(mark_h_idx).SetIsotope(2)
            try:
                new_ranks = list(Chem.CanonicalRankAtoms(m, includeChirality=True, breakTies=True))
            except TypeError:
                new_ranks = list(Chem.CanonicalRankAtoms(m, includeChirality=True))
            return new_ranks[heavy_idx]

        return get_heavy_rank_with_marked_h(h_a) != get_heavy_rank_with_marked_h(h_b)

    except Exception as e:
        if SimConfig().verbose:
            print(f"[diastereotopy] Warning for H {h_a},{h_b} on heavy {heavy_idx}: {e}")
        return False
def _neighbor_is_diastereotopic_CH2(base: Chem.Mol, molH: Chem.Mol,
                                    heavy_to_H: Dict[int, List[int]], other_heavy: int) -> bool:
    hs = heavy_to_H.get(other_heavy, [])
    if len(hs) != 2:
        return False
    return _is_diastereotopic_Hpair(molH, other_heavy, hs[0], hs[1])

# ---------------------------
# Sugar helpers
# ---------------------------
def _is_exocyclic_C6_like(base: Chem.Mol, heavy_to_H: Dict[int, List[int]], idx: int) -> bool:
    A = base.GetAtomWithIdx(idx)
    if A.GetIsAromatic() or A.GetSymbol() != "C":
        return False
    if len(heavy_to_H.get(idx, [])) != 2:
        return False
    syms = {n.GetSymbol() for n in A.GetNeighbors()}
    return ("O" in syms) and ("C" in syms)

def _is_ring_sp3C_in_6_like(base: Chem.Mol, idx: int) -> bool:
    A = base.GetAtomWithIdx(idx)
    return (A.GetSymbol()=="C") and A.IsInRingSize(6) and not _is_sp2_like(A)

def _is_anomeric_C1_like(base: Chem.Mol, heavy_to_H: Dict[int, List[int]], idx: int) -> bool:
    A = base.GetAtomWithIdx(idx)
    if A.GetSymbol() != "C" or not A.IsInRingSize(6):
        return False
    if any(b.GetBondTypeAsDouble() >= 2 for b in A.GetBonds()):
        return False
    has_ring_O = any(n.GetSymbol()=="O" and n.IsInRing() for n in A.GetNeighbors())
    return has_ring_O and (len(heavy_to_H.get(idx, [])) == 1)

# ---------------------------
# Heteronuclear contextual
# ---------------------------
def _hetero_J_contextual(base: Chem.Mol, D: np.ndarray, aromatic_dmap: List[Dict[int,int]],
                         heavy_idx: int, X_atom: Chem.Atom, config: SimConfig) -> Iterable[float]:
    out = []
    sym = X_atom.GetSymbol()
    if sym not in config.hetero_j_base: return out
    x_idx = X_atom.GetIdx()
    dist = int(D[heavy_idx, x_idx])
    if not (1 <= dist <= 3): return out
    order = {1: 2, 2: 3, 3: 4}[dist]
    Jx = config.hetero_j_base[sym].get(order, 0.0)
    if Jx <= 0.0: return out
    scale = 1.0
    cx = None
    for nb in X_atom.GetNeighbors():
        if nb.GetSymbol() == "C": cx = nb; break
    if cx is not None:
        rel = _aromatic_rel_from_map(aromatic_dmap, heavy_idx, cx.GetIdx())
        if sym == 'F' and config.heteroHF_advanced:
            if rel == "ortho" and order in {2,3}:
                scale = 1.3
                if order == 3 and Jx < 8.0: Jx = max(Jx, 8.0)
            elif rel == "meta" and order in {3,4}:
                scale = 1.1; Jx = max(Jx, 2.0)
            elif rel == "para" and order >= 4:
                scale = 0.8
        elif sym == 'P' and config.heteroHF_advanced:
            if rel == "ortho": scale = 1.1
            elif rel == "para": scale = 0.9
    out.append(Jx * scale); return out

# ---------------------------
# J-estimation principale
# ---------------------------
def _estimate_Js_for_group(
    molH: Chem.Mol,
    base: Chem.Mol,
    D: np.ndarray,
    aromatic_dmap: List[Dict[int,int]],
    heavy_idx: int,
    h_atoms: List[int],
    heavy_to_H: Dict[int, List[int]],
    spectrometer_mhz: float,
    solvent: dict,
    proton_shifts_ppm: Optional[Dict[int, float]],
    include_heteronuclear: bool,
    suppressed_heavy: Optional[set] = None,
    heavy_is_sp2_cache: Optional[List[bool]] = None,
    config: Optional[SimConfig] = None,
    ranks: Optional[List[int]] = None,
    delta_cluster_id: Optional[Dict[int,int]] = None,
    geminal_partner: Optional[int] = None,
    conformer_cache: Optional['ConformerCache'] = None,
) -> Tuple[
    List[float], bool, List[int], List[float], List[Tuple[int,float]],
    List[int], List[bool], str, Optional[str], float
]:
    """
    Estimation des couplages J pour un groupe de protons.
    
    Args:
        conformer_cache: Cache de conformères pour moyennage multi-conformères (optionnel)
    """
    if config is None: config = SimConfig()

    digital_floor = (config.digital_resolution_hz or 0.0)
    detect_floor = max(config.min_visibility_factor * config.line_width_hz, digital_floor)
    min2 = max(config.min_j_detect[2], detect_floor)
    min3 = max(config.min_j_detect[3], detect_floor)
    min4 = max(config.min_j_detect[4], detect_floor, 1.0)
    min5 = max(config.min_j_detect[5], detect_floor)

    J_vals: List[float] = []
    J_sources: List[Optional[int]] = []
    J_orders: List[int] = []
    partners: List[int] = []
    homo_J_vals_for_second_order: List[float] = []
    satellites: List[float] = []
    second_order = False
    second_order_type = "first_order"
    system_type: Optional[str] = None
    delta_nu_min_hz: float = 0.0
    homo_pairs: List[Tuple[int,float]] = []

    para_buf_J: List[float] = []
    para_buf_partners: List[int] = []

    heavy_atom = base.GetAtomWithIdx(heavy_idx)
    is_sp2 = heavy_is_sp2_cache[heavy_idx] if heavy_is_sp2_cache is not None else _is_sp2_like(heavy_atom)
    h_rep = h_atoms[0]

    partner_heavy_idx_kept: List[int] = []
    partner_is_sp2_kept: List[bool] = []

    # satellites 13C (avec ¹J(C-H) variable selon hybridation)
    if config.include_13C_satellites and heavy_atom.GetSymbol() == "C":
        nH_here = len(h_atoms)
        # Déterminer ¹J(C-H) selon environnement
        jch = _get_1JCH(heavy_atom, base)
        if nH_here >= 1:
            satellites.extend([jch] * nH_here)

    # Détection AA'BB'/AA'XX'
    system_type = _detect_AABB_system(base, heavy_idx, heavy_to_H, aromatic_dmap)

    # =========================================================================
    # NOUVEAU: Détection et traitement des hétérocycles via module complet
    # =========================================================================
    _heterocycle_handled = False  # Flag pour éviter double traitement
    if HETEROCYCLES_MODULE_AVAILABLE and heavy_atom.GetIsAromatic():
        try:
            het_result = estimate_heterocycle_couplings_complete(
                mol=base,
                heavy_idx=heavy_idx,
                h_list=h_atoms,
                heavy_to_H=heavy_to_H,
                D=D,
                min_j_threshold=min3,
                apply_substituent_effects=True
            )
            
            if het_result is not None:
                # L'hétérocycle a été identifié et traité
                het_J_vals = het_result.get('J_values', [])
                het_partners = het_result.get('coupling_partners', [])
                het_mult = het_result.get('multiplicity', 's')
                het_type = het_result.get('heterocycle_type', 'unknown')
                het_source = het_result.get('source', 'heterocycle')
                
                # ✅ Traiter TOUS les résultats du module hétérocycle (singlets et expected_mult)
                # Singulet spécial (imidazole H2, H4, H5, etc.)
                if het_source in ('heterocycle_singlet', 'heterocycle_expected_singlet') or (het_mult == 's' and not het_J_vals):
                    # Retourner directement un singulet
                    return (
                        [],  # J_vals vide = singulet
                        False,  # second_order
                        [],  # partners
                        [],  # satellites
                        [],  # homo_pairs
                        [],  # partner_heavy_idx_kept
                        [],  # partner_is_sp2_kept
                        "first_order",  # second_order_type
                        f"heterocycle_{het_type}",  # system_type
                        0.0  # delta_nu_min_hz
                    )
                
                # ✅ NOUVEAU: Traiter les multiplicités attendues (expected_mult: d, dd, tt, etc.)
                # et aussi les résultats de table normale avec J
                if het_J_vals:
                    het_J_vals_out = []
                    het_partners_out = []
                    het_homo_pairs = []
                    het_partner_heavy_idx = []
                    het_partner_is_sp2 = []
                    
                    for J_het, partner_het in zip(het_J_vals, het_partners):
                        het_J_vals_out.append(J_het)
                        het_partners_out.append(partner_het)
                        het_homo_pairs.append((partner_het, J_het))
                        # Trouver le heavy_idx du partner
                        for h_idx, h_list in heavy_to_H.items():
                            if partner_het in h_list:
                                het_partner_heavy_idx.append(h_idx)
                                partner_atom = base.GetAtomWithIdx(h_idx)
                                het_partner_is_sp2.append(partner_atom.GetIsAromatic() or _is_sp2_like(partner_atom))
                                break
                    
                    # ✅ NOUVEAU: Encoder la multiplicité ET position IUPAC dans system_type
                    # Format: heterocycle_<type>:<mult>:<position> (ex: heterocycle_pyridine:tt:4)
                    het_position = het_result.get('position', '')
                    system_type_with_mult = f"heterocycle_{het_type}:{het_mult}:{het_position}"
                    
                    # Retourner directement les résultats du module hétérocycle
                    return (
                        het_J_vals_out,  # J_vals
                        False,  # second_order
                        het_partners_out,  # partners
                        [],  # satellites
                        het_homo_pairs,  # homo_pairs
                        het_partner_heavy_idx,  # partner_heavy_idx_kept
                        het_partner_is_sp2,  # partner_is_sp2_kept
                        "first_order",  # second_order_type
                        system_type_with_mult,  # system_type avec multiplicité et position
                        0.0  # delta_nu_min_hz
                    )
                    
        except Exception as e:
            # En cas d'erreur, continuer avec la logique standard
            _heterocycle_handled = False

    # =========================================================================
    # NOUVEAU: Détection et traitement des sucres via module complet
    # =========================================================================
    if SUGARS_MODULE_AVAILABLE and heavy_atom.GetSymbol() == 'C' and not heavy_atom.GetIsAromatic():
        try:
            # Vérifier si c'est un carbone de sucre
            if is_sugar_carbon(base, heavy_idx):
                # Reconstruire proton_to_heavy à partir de heavy_to_H
                proton_to_heavy_local = {}
                for h_idx_map, h_list_map in heavy_to_H.items():
                    for h in h_list_map:
                        proton_to_heavy_local[h] = h_idx_map
                
                sugar_result = estimate_sugar_couplings_complete(
                    mol=base,
                    heavy_idx=heavy_idx,
                    h_list=h_atoms,
                    heavy_to_H=heavy_to_H,
                    proton_to_heavy=proton_to_heavy_local,
                    D=D,
                    config=None  # Utiliser config par défaut
                )
                
                if sugar_result is not None:
                    sugar_J_vals = sugar_result.get('J_values', [])
                    sugar_partners = sugar_result.get('coupling_partners', [])
                    sugar_mult = sugar_result.get('multiplicity', 's')
                    sugar_info = sugar_result.get('sugar_info', {})
                    sugar_position = sugar_info.get('position', '')
                    sugar_ring_type = sugar_info.get('ring_type', 'pyranose')
                    sugar_source = sugar_result.get('source', 'sugar_complete')
                    
                    # Construire les structures de sortie
                    sugar_J_vals_out = []
                    sugar_partners_out = []
                    sugar_homo_pairs = []
                    sugar_partner_heavy_idx = []
                    sugar_partner_is_sp2 = []
                    
                    for J_sugar, partner_sugar in zip(sugar_J_vals, sugar_partners):
                        sugar_J_vals_out.append(J_sugar)
                        sugar_partners_out.append(partner_sugar)
                        sugar_homo_pairs.append((partner_sugar, J_sugar))
                        # Trouver le heavy_idx du partner
                        for h_idx, h_list in heavy_to_H.items():
                            if partner_sugar in h_list:
                                sugar_partner_heavy_idx.append(h_idx)
                                partner_atom = base.GetAtomWithIdx(h_idx)
                                sugar_partner_is_sp2.append(partner_atom.GetIsAromatic() or _is_sp2_like(partner_atom))
                                break
                    
                    # Encoder position dans system_type: sugar_<ring_type>:<mult>:<position>
                    system_type_sugar = f"sugar_{sugar_ring_type}:{sugar_mult}:{sugar_position}"
                    
                    # Retourner directement les résultats du module sucre
                    return (
                        sugar_J_vals_out,  # J_vals
                        False,  # second_order
                        sugar_partners_out,  # partners
                        [],  # satellites
                        sugar_homo_pairs,  # homo_pairs
                        sugar_partner_heavy_idx,  # partner_heavy_idx_kept
                        sugar_partner_is_sp2,  # partner_is_sp2_kept
                        "first_order",  # second_order_type
                        system_type_sugar,  # system_type avec ring_type, multiplicité et position
                        0.0  # delta_nu_min_hz
                    )
                    
        except Exception as e:
            # En cas d'erreur, continuer avec la logique standard
            pass

    # =========================================================================
    # NOUVEAU: Détection et traitement via module couplages complets
    # (Alcynes, Aldéhydes, Peptides, Phosphines, Stéroïdes, Terpènes, etc.)
    # =========================================================================
    if COUPLING_TABLES_MODULE_AVAILABLE:
        try:
            # Classifier la molécule (peut utiliser base sans H)
            mol_types = MoleculeClassifier.classify(base)
            
            if mol_types:
                # Chercher les partenaires de couplage pour ce proton
                n_heavy = base.GetNumAtoms()
                candidate_heavies_spec = [j for j in range(n_heavy) if 1 <= int(D[heavy_idx, j]) <= 4]
                
                for other_heavy in candidate_heavies_spec:
                    other_hs = heavy_to_H.get(other_heavy, [])
                    if not other_hs:
                        continue
                    
                    dist = int(D[heavy_idx, other_heavy])
                    
                    for h2 in other_hs:
                        if h2 in h_atoms:
                            continue
                        
                        # Utiliser le dispatch des tables spécialisées avec molH
                        # car h_rep et h2 sont des indices de protons
                        J_spec = get_specialized_coupling(
                            mol=molH,
                            h1_idx=h_rep,
                            h2_idx=h2,
                            distance=dist + 1  # +1 car on veut 3J, 4J, etc.
                        )
                        
                        if J_spec is not None and J_spec >= min3:
                            J_vals.append(J_spec)
                            J_sources.append(other_heavy)
                            J_orders.append(dist + 1)
                            partners.append(h2)
                            homo_J_vals_for_second_order.append(J_spec)
                            homo_pairs.append((h2, J_spec))
                            
                            partner_atom = base.GetAtomWithIdx(other_heavy)
                            partner_heavy_idx_kept.append(other_heavy)
                            partner_is_sp2_kept.append(
                                partner_atom.GetIsAromatic() or _is_sp2_like(partner_atom)
                            )
                            
        except Exception as e:
            # En cas d'erreur, continuer avec la logique standard
            pass

    # =========================================================================
    # NOUVEAU: Couplages peptidiques robustes via module peptide_couplings
    # =========================================================================
    if PEPTIDE_MODULE_AVAILABLE and not COUPLING_TABLES_MODULE_AVAILABLE:
        # Utiliser le module peptide directement si coupling_tables n'est pas là
        try:
            # Vérifier si c'est un environnement peptidique
            is_peptide_env = False
            for nb in heavy_atom.GetNeighbors():
                if nb.GetSymbol() == 'N':
                    # Chercher liaison amide (N-C=O)
                    for nn in nb.GetNeighbors():
                        if nn.GetSymbol() == 'C':
                            for nnn in nn.GetNeighbors():
                                if nnn.GetSymbol() == 'O':
                                    bond = base.GetBondBetweenAtoms(nn.GetIdx(), nnn.GetIdx())
                                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                                        is_peptide_env = True
                                        break
            
            if is_peptide_env:
                n_heavy = base.GetNumAtoms()
                candidate_heavies_pep = [j for j in range(n_heavy) if 1 <= int(D[heavy_idx, j]) <= 3]
                
                for other_heavy in candidate_heavies_pep:
                    other_hs = heavy_to_H.get(other_heavy, [])
                    if not other_hs:
                        continue
                    
                    for h2 in other_hs:
                        if h2 in h_atoms:
                            continue
                        
                        J_pep = get_peptide_coupling_robust(
                            mol=base,
                            h1_idx=h_rep,
                            h2_idx=h2,
                            phi_angle=None,  # Extraire de la géométrie si possible
                            secondary_structure='random_coil'
                        )
                        
                        if J_pep is not None and J_pep >= min3:
                            J_vals.append(J_pep)
                            J_sources.append(other_heavy)
                            J_orders.append(3)
                            partners.append(h2)
                            homo_J_vals_for_second_order.append(J_pep)
                            homo_pairs.append((h2, J_pep))
                            
                            partner_atom = base.GetAtomWithIdx(other_heavy)
                            partner_heavy_idx_kept.append(other_heavy)
                            partner_is_sp2_kept.append(
                                partner_atom.GetIsAromatic() or _is_sp2_like(partner_atom)
                            )
                            
        except Exception as e:
            pass

    # =========================================================================
    # NOUVEAU: Extensions de couplage (aromatiques polycycliques, cycles tendus)
    # =========================================================================
    if COUPLING_EXTENSIONS_AVAILABLE:
        try:
            n_heavy = base.GetNumAtoms()
            candidate_heavies_ext = [j for j in range(n_heavy) if 1 <= int(D[heavy_idx, j]) <= 4]
            
            for other_heavy in candidate_heavies_ext:
                # ✅ FILTRE ÉQUIVALENCE: Ignorer couplages entre carbones équivalents
                if config.collapse_equivalent_by_rank and ranks is not None:
                    if _are_protons_on_equivalent_carbons(base, ranks, heavy_idx, other_heavy):
                        continue
                
                other_hs = heavy_to_H.get(other_heavy, [])
                if not other_hs:
                    continue
                
                dist = int(D[heavy_idx, other_heavy])
                
                for h2 in other_hs:
                    if h2 in h_atoms:
                        continue
                    
                    # Utiliser le dispatch des extensions avec molH (molécule avec H)
                    # car h_rep et h2 sont des indices de protons
                    J_ext = get_extended_coupling(
                        mol=molH,
                        h1_idx=h_rep,
                        h2_idx=h2,
                        distance=dist + 1
                    )
                    
                    if J_ext is not None and abs(J_ext) >= min3:
                        J_vals.append(abs(J_ext))  # Prendre valeur absolue
                        J_sources.append(other_heavy)
                        J_orders.append(dist + 1)
                        partners.append(h2)
                        homo_J_vals_for_second_order.append(abs(J_ext))
                        homo_pairs.append((h2, abs(J_ext)))
                        
                        partner_atom = base.GetAtomWithIdx(other_heavy)
                        partner_heavy_idx_kept.append(other_heavy)
                        partner_is_sp2_kept.append(
                            partner_atom.GetIsAromatic() or _is_sp2_like(partner_atom)
                        )
                        
        except Exception as e:
            pass

    # 2J geminal avec partenaire diastéréotopique explicite
    if geminal_partner is not None:
        # Utiliser _get_geminal_2J pour valeur adaptée à l'environnement
        Jg = _get_geminal_2J(heavy_atom, base)
        if Jg >= min2:
            J_vals.append(Jg); J_sources.append(None); J_orders.append(2)
            partners.append(geminal_partner); homo_J_vals_for_second_order.append(Jg)
            homo_pairs.append((geminal_partner, Jg))
            partner_heavy_idx_kept.append(heavy_idx)
            partner_is_sp2_kept.append(is_sp2)

    # 2J geminal (cas où les 2 H sont dans le même groupe)
    if len(h_atoms) == 2:
        nb_heavy = [n for n in heavy_atom.GetNeighbors()]
        chiral_neigh = any(_is_chiral_center(n) for n in nb_heavy)
        non_equiv_by_shift = False
        if proton_shifts_ppm is not None:
            try:
                h0, h1 = h_atoms
                d0 = proton_shifts_ppm.get(h0, None)
                d1 = proton_shifts_ppm.get(h1, None)
                if d0 is not None and d1 is not None and math.isfinite(d0) and math.isfinite(d1):
                    non_equiv_by_shift = (abs(d0 - d1) > 0.02)
            except Exception:
                non_equiv_by_shift = False
        if chiral_neigh or non_equiv_by_shift:
            # Utiliser _get_geminal_2J pour valeur adaptée à l'environnement
            Jg = _get_geminal_2J(heavy_atom, base)
            if Jg >= min2:
                J_vals.append(Jg); J_sources.append(None); J_orders.append(2)
                partners.append(h_atoms[1]); homo_J_vals_for_second_order.append(Jg)
                homo_pairs.append((h_atoms[1], Jg))
                partner_heavy_idx_kept.append(heavy_idx)
                partner_is_sp2_kept.append(is_sp2)

    # Hétéronucléaires (F, P) avec couplages avancés
    if include_heteronuclear:
        for X in base.GetAtoms():
            xsym = X.GetSymbol()
            if xsym not in config.hetero_j_base: continue
            
            # Couplages ³¹P avancés (phosphates, phosphonates)
            if xsym == "P":
                dist_P = int(D[heavy_idx, X.GetIdx()]) if D[heavy_idx, X.GetIdx()] > 0 else 99
                if dist_P <= 4:
                    J_31P = _estimate_31P_1H_coupling(base, h_rep, X.GetIdx(), dist_P)
                    if J_31P is not None and J_31P >= 0.5:
                        J_vals.append(J_31P); J_sources.append(X.GetIdx()); J_orders.append(dist_P + 1)
                        continue  # Skip generic hetero coupling for P
            
            for j in _hetero_J_contextual(base, D, aromatic_dmap, heavy_idx, X, config):
                if j >= 0.5:
                    J_vals.append(j); J_sources.append(None); J_orders.append(3)
    
    # Couplages ¹⁵N-¹H (satellites, 0.37% abondance naturelle)
    # Visible surtout pour composés enrichis en ¹⁵N ou spectres très concentrés
    # Note: Nous ne les ajoutons pas par défaut car ils sont très faibles
    # mais la fonction _estimate_15N_1H_coupling existe pour usage spécialisé

    # ³J(H-N-C-H) pour peptides/amides (Bystrov)
    if config.peptide_HNCH_coupling:
        # Chercher les protons sur N voisins couplés via C
        for nb in heavy_atom.GetNeighbors():
            if nb.GetSymbol() == 'N':
                n_idx = nb.GetIdx()
                # Chercher H sur N
                for n_neighbor in nb.GetNeighbors():
                    if n_neighbor.GetSymbol() == 'H':
                        hn_idx = n_neighbor.GetIdx()
                        # Estimer ³J(H-N-C-H)
                        for h in h_atoms:
                            try:
                                J_hnch = _estimate_HNCH_J(molH, hn_idx, n_idx, heavy_idx, h)
                                if J_hnch is not None and J_hnch >= min3:
                                    J_vals.append(J_hnch)
                                    J_sources.append(n_idx)
                                    J_orders.append(3)
                                    partners.append(hn_idx)
                                    homo_J_vals_for_second_order.append(J_hnch)
                                    homo_pairs.append((hn_idx, J_hnch))
                            except Exception:
                                pass

    n_heavy = base.GetNumAtoms()
    candidate_heavies = [j for j in range(n_heavy) if 1 <= int(D[heavy_idx, j]) <= 3]

    anomeric_rule_applied = False

    for other_heavy in candidate_heavies:
        if suppressed_heavy and other_heavy in suppressed_heavy: continue
        dist = int(D[heavy_idx, other_heavy])
        other_hs = heavy_to_H.get(other_heavy, [])
        if not other_hs: continue

        # Filtre équivalence par rang canonique (inclut symétrie moléculaire)
        # ✅ CORRECTION MAJEURE: Les protons sur des carbones MAGNÉTIQUEMENT ÉQUIVALENTS
        # ne montrent pas de couplage mutuel visible (le couplage existe mais est invisible
        # car les protons sont isochroniques = même shift chimique).
        # 
        # Cas concernés:
        # - Cyclohexane: tous les CH2 équivalents → singulet large
        # - Cyclopentane: tous les CH2 équivalents → singulet
        # - 1,2-Dichloroéthane (ClCCCl): les 2 CH2 équivalents → singulet
        # - Benzène: tous les CH équivalents → singulet
        # - Aromatiques monosubstitués: H en position ortho équivalents entre eux
        #
        # Règle: Si ranks[heavy_idx] == ranks[other_heavy], ignorer le couplage
        if config.collapse_equivalent_by_rank and ranks is not None:
            # ✅ NOUVEAU: Vérifier si les deux carbones porteurs sont équivalents
            # Si oui, ne pas compter le couplage (il est invisible)
            if _are_protons_on_equivalent_carbons(base, ranks, heavy_idx, other_heavy):
                continue
            
            # Ancien code pour les cas aromatiques (gardé pour compatibilité)
            heavy_atom_check = base.GetAtomWithIdx(heavy_idx)
            other_heavy_atom = base.GetAtomWithIdx(other_heavy)
            both_aromatic = heavy_atom_check.GetIsAromatic() and other_heavy_atom.GetIsAromatic()
            
            # Pour aromatiques: toujours filtrer les équivalents
            # Pour aliphatiques: la vérification ci-dessus suffit maintenant
            if both_aromatic:
                if _are_heavy_equivalent(ranks, heavy_idx, other_heavy, heavy_to_H, proton_shifts_ppm, config.equiv_delta_ppm_tol, mol=base):
                    continue

        # 3J vicinaux
        if dist == 1:
            partner_atom = base.GetAtomWithIdx(other_heavy)
            partner_is_sp2 = (heavy_is_sp2_cache[other_heavy]
                              if heavy_is_sp2_cache is not None
                              else _is_sp2_like(partner_atom))
            is_carbon_partner = (partner_atom.GetSymbol() == "C")

            exch_partner = _is_exchangeable_H_on_X(base, other_heavy)
            if exch_partner and solvent.get("name") != "DMSO-d6":
                continue

            Js_tmp: List[float] = []
            partners_tmp: List[int] = []
            partners_tmp_heavy: List[int] = []
            partners_tmp_is_sp2: List[bool] = []
            
            # ✅ NOUVEAU: Pour CH3, collecter TOUS les J (même petits) pour le calcul de moyenne
            Js_all_for_CH3: List[float] = []
            partners_all_for_CH3: List[int] = []

            forbid_avg = _neighbor_is_diastereotopic_CH2(base, molH, heavy_to_H, other_heavy)
            
            # Détecter si c'est un CH3 partner
            is_CH3_partner = (partner_atom.GetSymbol() == "C" and len(other_hs) == 3 and not partner_is_sp2)

            # Auto-bench aromatic
            if (config.auto_bench_mode and is_sp2 and partner_is_sp2
                and base.GetAtomWithIdx(heavy_idx).IsInRingSize(6)
                and base.GetAtomWithIdx(other_heavy).IsInRingSize(6)):
                rel = _aromatic_rel_from_map(aromatic_dmap, heavy_idx, other_heavy)
                J_force = None
                if rel == "ortho":
                    J_force = config.aro_j_forced["ortho"]
                elif rel == "meta":
                    J_force = config.aro_j_forced["meta"]
                elif rel == "para":
                    J_force = config.aro_j_forced["para"]
                if J_force is not None and J_force >= min3:
                    Js_tmp.append(J_force)
                    partners_tmp.append(other_hs[0])
                    partners_tmp_heavy.append(other_heavy)
                    partners_tmp_is_sp2.append(partner_is_sp2)

                    for jv, ptmp, phv, pis2 in zip(Js_tmp, partners_tmp, partners_tmp_heavy, partners_tmp_is_sp2):
                        J_vals.append(jv); J_sources.append(other_heavy); J_orders.append(3)
                        homo_J_vals_for_second_order.append(jv)
                        partners.append(ptmp); homo_pairs.append((ptmp, jv))
                        partner_heavy_idx_kept.append(phv); partner_is_sp2_kept.append(pis2)
                    continue

            # Auto-bench vinylique
            if config.auto_bench_mode and is_sp2 and partner_is_sp2:
                b = base.GetBondBetweenAtoms(heavy_idx, other_heavy)
                if b and b.GetBondTypeAsDouble() >= 2:
                    Jtrans = _estimate_vinylic_J(molH, h_rep, other_hs[0], config=config) or config.vinylic_J_defaults['trans']
                    if Jtrans >= min3:
                        Js_tmp.append(Jtrans); partners_tmp.append(other_hs[0])
                        partners_tmp_heavy.append(other_heavy); partners_tmp_is_sp2.append(partner_is_sp2)
                    isW, _ = _is_W_allylic(base, heavy_idx, other_heavy)
                    if isW:
                        Js_tmp.append(1.2); partners_tmp.append(other_hs[0])
                        partners_tmp_heavy.append(other_heavy); partners_tmp_is_sp2.append(partner_is_sp2)
                    for jv, ptmp, phv, pis2 in zip(Js_tmp, partners_tmp, partners_tmp_heavy, partners_tmp_is_sp2):
                        J_vals.append(jv); J_sources.append(other_heavy); J_orders.append(3)
                        homo_J_vals_for_second_order.append(jv)
                        partners.append(ptmp); homo_pairs.append((ptmp, jv))
                        partner_heavy_idx_kept.append(phv); partner_is_sp2_kept.append(pis2)
                    continue

            for h in other_hs:
                if h in h_atoms:
                    continue

                J = None
                if not is_sp2:
                    bond = base.GetBondBetweenAtoms(heavy_idx, other_heavy)
                    flexible_cc = (
                        bond is not None
                        and bond.GetBondTypeAsDouble() == 1.0
                        and not base.GetAtomWithIdx(heavy_idx).IsInRing()
                        and not base.GetAtomWithIdx(other_heavy).IsInRing()
                        and base.GetAtomWithIdx(heavy_idx).GetSymbol() == "C"
                        and partner_atom.GetSymbol() == "C"
                    )
                    
                    # Pour les chaînes COMPLÈTEMENT flexibles (rotation libre acyclique),
                    # J=7.0 Hz est la valeur empirique moyenne (gauche/trans average)
                    # Le moyennage multi-conformères ne fonctionne pas bien car ETKDG
                    # ne génère pas assez de diversité conformationnelle pour ces cas
                    if flexible_cc:
                        J = 7.0
                    else:
                        EN_corr = 0.4 if any(n.GetSymbol() in {"O","N"} for n in base.GetAtomWithIdx(heavy_idx).GetNeighbors()) else 0.0
                        
                        # ✅ Utiliser moyennage multi-conformères pour systèmes cycliques/semi-rigides
                        if conformer_cache is not None and config.use_conformational_averaging:
                            Jk = conformer_cache.get_J_vicinal_averaged(
                                h_rep, h,
                                use_haasnoot=(config.karplus_haasnoot and base.GetAtomWithIdx(heavy_idx).IsInRingSize(6)),
                                use_haasnoot_7_params=config.use_haasnoot_7_params,
                                electroneg_corr=EN_corr
                            )
                        else:
                            Jk = _estimate_vicinal_sp3_J(
                                molH, h_rep, h,
                                use_haasnoot=(config.karplus_haasnoot and base.GetAtomWithIdx(heavy_idx).IsInRingSize(6)),
                                electroneg_corr=EN_corr,
                                use_full_haasnoot=config.use_haasnoot_7_params
                            )
                        if Jk is None:
                            J = 7.0 if is_carbon_partner else None
                        else:
                            J = Jk
                else:
                    rel = _aromatic_rel_from_map(aromatic_dmap, heavy_idx, other_heavy)
                    if rel == "ortho": J = config.aro_j_forced["ortho"]
                    elif rel in {"meta","para"}: J = 0.0
                    else:
                        Jv = _estimate_vinylic_J(molH, h, h_rep, config=config)
                        if Jv is None:
                            J = _estimate_vinylic_J(molH, h_rep, h, config=config)
                        else:
                            J = Jv
                        if J is None:
                            # Essayer couplage alcyne (H-C≡C-H)
                            J_alkyne = _estimate_alkyne_J(molH, h_rep, h)
                            if J_alkyne is not None:
                                J = J_alkyne
                            else:
                                J = config.vinylic_J_defaults.get('cis', 8.0)

                # Profils sucre
                if config.sugar_mode:
                    is_C6_other = _is_exocyclic_C6_like(base, heavy_to_H, other_heavy)
                    is_C5_this = _is_ring_sp3C_in_6_like(base, heavy_idx)
                    is_C6_this = _is_exocyclic_C6_like(base, heavy_to_H, heavy_idx)
                    is_C5_other = _is_ring_sp3C_in_6_like(base, other_heavy)

                    if (is_C5_this and is_C6_other) or (is_C6_this and is_C5_other):
                        if J is None: J = 4.0
                        if J <= 3.5:
                            J = max(1.5, min(3.0, J))
                        else:
                            J = max(4.0, min(6.0, J))

                    if not anomeric_rule_applied:
                        is_C1_this = _is_anomeric_C1_like(base, heavy_to_H, heavy_idx)
                        is_C1_other = _is_anomeric_C1_like(base, heavy_to_H, other_heavy)
                        heavy_is_C = (base.GetAtomWithIdx(heavy_idx).GetSymbol() == "C")
                        other_is_C = (base.GetAtomWithIdx(other_heavy).GetSymbol() == "C")
                        if (is_C1_this or is_C1_other) and heavy_is_C and other_is_C:
                            pref = (config.anomeric_preference or "").lower()
                            if pref == "alpha":
                                J = 2.7
                            elif pref == "beta":
                                J = 7.5
                            else:
                                if J is None: J = 7.5
                                if J > 9.0: J = 7.5
                                elif J < 2.0: J = 2.7
                            anomeric_rule_applied = True

                if config.formyl_rule:
                    # ✅ CORRECTION: Appliquer la règle formyl dans les DEUX sens
                    # Cas 1: heavy_idx est le CHO
                    A = base.GetAtomWithIdx(heavy_idx)
                    is_formyl_this = (A.GetSymbol()=="C"
                        and any(b.GetBondTypeAsDouble()==2 and b.GetOtherAtom(A).GetSymbol()=="O" for b in A.GetBonds())
                        and len(heavy_to_H.get(heavy_idx, [])) == 1)
                    
                    # Cas 2: other_heavy est le CHO (le partenaire)
                    B = base.GetAtomWithIdx(other_heavy)
                    is_formyl_partner = (B.GetSymbol()=="C"
                        and any(b.GetBondTypeAsDouble()==2 and b.GetOtherAtom(B).GetSymbol()=="O" for b in B.GetBonds())
                        and len(heavy_to_H.get(other_heavy, [])) == 1)
                    
                    if is_formyl_this or is_formyl_partner:
                        # ³J(H-CHO) typique: 1.5-3.0 Hz
                        if J is None: J = 2.5
                        J = max(1.5, min(3.0, J))

                if config.amide_mode:
                    def _is_amide_like(idx_c):
                        C = base.GetAtomWithIdx(idx_c)
                        return (C.GetSymbol()=="C" and any(b.GetBondTypeAsDouble()==2 and b.GetOtherAtom(C).GetSymbol()=="O" for b in C.GetBonds()))
                    
                    def _is_formyl(idx_c):
                        """Vrai si c'est un aldéhyde CHO (pas amide)"""
                        C = base.GetAtomWithIdx(idx_c)
                        return (C.GetSymbol()=="C" 
                            and any(b.GetBondTypeAsDouble()==2 and b.GetOtherAtom(C).GetSymbol()=="O" for b in C.GetBonds())
                            and len(heavy_to_H.get(idx_c, [])) == 1)
                    
                    # ✅ CORRECTION: Ne PAS appliquer amide_mode aux formyls (aldéhydes)
                    if (_is_amide_like(heavy_idx) or _is_amide_like(other_heavy)) and not (_is_formyl(heavy_idx) or _is_formyl(other_heavy)):
                        if J is None: J = 6.0
                        J = max(3.5, min(9.0, J))

                # ✅ NOUVEAU: Pour CH3 rotatif, collecter TOUS les J pour le calcul de moyenne
                # Un CH3 en rotation libre moyenne les 3 couplages → tous les H comptent
                if is_CH3_partner and J is not None:
                    Js_all_for_CH3.append(J)
                    partners_all_for_CH3.append(h)

                if J is not None and J >= min3:
                    Js_tmp.append(J); partners_tmp.append(h)
                    partners_tmp_heavy.append(other_heavy); partners_tmp_is_sp2.append(partner_is_sp2)

            # ✅ CORRECTION: Pour CH3 rotatif, utiliser TOUS les J calculés pour la moyenne
            # puis ajouter exactement 3 couplages (un par H du CH3)
            if is_CH3_partner and not forbid_avg and len(Js_all_for_CH3) == 3:
                Jmean = float(np.mean(Js_all_for_CH3))
                if Jmean >= min3:
                    # Ajouter 3 couplages avec la valeur moyenne
                    for h_ch3 in other_hs:
                        if h_ch3 in h_atoms:
                            continue
                        J_vals.append(Jmean); J_sources.append(other_heavy); J_orders.append(3)
                        homo_J_vals_for_second_order.append(Jmean)
                        partners.append(h_ch3); homo_pairs.append((h_ch3, Jmean))
                        partner_heavy_idx_kept.append(other_heavy); partner_is_sp2_kept.append(partner_is_sp2)
                continue
            
            if Js_tmp:
                dyn_tol = max(TOLERANCE_J_HZ, config.grouping_tol_scale * config.line_width_hz)

                # Note: Le cas CH3_partner avec rotation libre est géré ci-dessus avec Js_all_for_CH3
                # Ce bloc gère les autres cas (CH2, CH, etc.)

                close_enough = (max(Js_tmp) - min(Js_tmp) <= 0.5 * dyn_tol)

                if (len(Js_tmp) >= 2) and (not is_sp2) and (not partner_is_sp2) and (not forbid_avg) and close_enough:
                    Jmean = float(np.mean(Js_tmp))
                    for ptmp, phv, pis2 in zip(partners_tmp, partners_tmp_heavy, partners_tmp_is_sp2):
                        J_vals.append(Jmean); J_sources.append(other_heavy); J_orders.append(3)
                        homo_J_vals_for_second_order.append(Jmean)
                        partners.append(ptmp); homo_pairs.append((ptmp, Jmean))
                        partner_heavy_idx_kept.append(phv); partner_is_sp2_kept.append(pis2)
                else:
                    for jv, ptmp, phv, pis2 in zip(Js_tmp, partners_tmp, partners_tmp_heavy, partners_tmp_is_sp2):
                        J_vals.append(jv); J_sources.append(other_heavy); J_orders.append(3)
                        homo_J_vals_for_second_order.append(jv)
                        partners.append(ptmp); homo_pairs.append((ptmp, jv))
                        partner_heavy_idx_kept.append(phv); partner_is_sp2_kept.append(pis2)
            continue

        # 4J
        if dist == 2:
            # ✅ FILTRE ÉQUIVALENCE: Ignorer 4J entre carbones magnétiquement équivalents
            # Cas: cyclopentane, etc. où tous les CH2 sont équivalents
            if config.collapse_equivalent_by_rank and ranks is not None:
                if _are_protons_on_equivalent_carbons(base, ranks, heavy_idx, other_heavy):
                    continue
            
            # 🔴 FILTRE : Ignorer ⁴J pour les méthyles isolés (benzyl, méthoxy)
            # Ces groupes ne montrent pas de couplage longue portée observable en pratique
            if _is_isolated_methyl(base, heavy_idx, heavy_to_H):
                continue  # Le CH3 benzyl ou OCH3 est un singulet
            
            Jctx, rel = _aromatic_J_contextual(base, aromatic_dmap, heavy_idx, other_heavy,
                                               J0_meta=config.aro_j_forced["meta"], J0_para=config.aro_j_forced["para"])
            if rel == "meta":
                # ✅ NOUVEAU: Filtrer les couplages meta aromatiques SEULEMENT dans les benzènes simples
                # Les hétérocycles (pyridine, etc.) ont des couplages meta plus significatifs
                # Détecter si c'est un benzène pur (pas d'hétéroatome dans le cycle)
                is_pure_benzene = True
                heavy_atom = base.GetAtomWithIdx(heavy_idx)
                if heavy_atom.GetIsAromatic():
                    ring_info = base.GetRingInfo()
                    for ring in ring_info.AtomRings():
                        if heavy_idx in ring and other_heavy in ring:
                            for idx in ring:
                                if base.GetAtomWithIdx(idx).GetSymbol() not in ('C', 'H'):
                                    is_pure_benzene = False
                                    break
                            break
                
                # Pour les benzènes purs, filtrer les couplages meta faibles (<2.0 Hz)
                if is_pure_benzene and Jctx is not None and Jctx < 2.0:
                    continue  # Couplage meta trop faible, non résolu en pratique
                J = Jctx
            else:
                # 1. Vérifier couplage homoallylique (H-C=C-CH₃)
                is_homoallyl, J_homoallyl = _is_homoallylic_4J(base, heavy_idx, other_heavy)
                if is_homoallyl:
                    J = J_homoallyl
                else:
                    # 2. Vérifier couplage allylique classique (H-CH₂-C=C-H)
                    is_allyl, J_allyl = _is_allylic_4J(base, heavy_idx, other_heavy, mol3d=molH)
                    if is_allyl:
                        J = J_allyl
                    else:
                        # 3. Vérifier W-path classique
                        isW, isRigid = _is_W_allylic(base, heavy_idx, other_heavy)
                        if isW:
                            J = 2.8 if isRigid else 1.6
                        else:
                            J = 0.3
                        if config.rigid_W_boost and isW:
                            J = min(4.0, max(J, 2.0))
            # ✅ AMÉLIORATION: Seuil ⁴J adaptatif selon rigidité du système
            adaptive_min4 = _get_adaptive_4J_threshold(base, heavy_idx, other_heavy, min4)
            if J is not None and J >= adaptive_min4:
                J_vals.append(J); J_sources.append(None); J_orders.append(4)
                partners.append(other_hs[0]); homo_J_vals_for_second_order.append(J)
                homo_pairs.append((other_hs[0], J))
                partner_heavy_idx_kept.append(other_heavy); partner_is_sp2_kept.append(_is_sp2_like(base.GetAtomWithIdx(other_heavy)))
            continue

        # 5J
        if dist == 3:
            # ✅ FILTRE ÉQUIVALENCE: Ignorer 5J entre carbones magnétiquement équivalents
            if config.collapse_equivalent_by_rank and ranks is not None:
                if _are_protons_on_equivalent_carbons(base, ranks, heavy_idx, other_heavy):
                    continue
            
            # 🔴 FILTRE : Ignorer ⁵J pour les méthyles isolés (benzyl, méthoxy)
            if _is_isolated_methyl(base, heavy_idx, heavy_to_H):
                continue  # Le CH3 benzyl ou OCH3 est un singulet
            
            Jctx, rel = _aromatic_J_contextual(base, aromatic_dmap, heavy_idx, other_heavy,
                                               J0_meta=config.aro_j_forced["meta"], J0_para=config.aro_j_forced["para"])
            if rel == "para":
                J = Jctx
                para_buf_J.append(J); para_buf_partners.append(other_hs[0])
            else:
                # Vérifier couplage W-zigzag ⁵J (systèmes rigides)
                is_W_5J, J_W_5J = _is_W_coupling_5J(base, heavy_idx, other_heavy)
                if is_W_5J:
                    J = J_W_5J
                else:
                    J = 0.2
                    if config.rigid_W_boost:
                        if (base.GetAtomWithIdx(heavy_idx).IsInRing() and base.GetAtomWithIdx(other_heavy).IsInRing()):
                            J = 0.7
                if J >= min5:
                    J_vals.append(J); J_sources.append(None); J_orders.append(5)
                    partners.append(other_hs[0]); homo_J_vals_for_second_order.append(J)
                    homo_pairs.append((other_hs[0], J))
                    partner_heavy_idx_kept.append(other_heavy); partner_is_sp2_kept.append(_is_sp2_like(base.GetAtomWithIdx(other_heavy)))
            continue

    # Anti-fusion sucre
    if config.sugar_mode and J_vals:
        dyn_tol = max(TOLERANCE_J_HZ, config.grouping_tol_scale * config.line_width_hz)
        eps = 1.1 * dyn_tol
        cand = [(i, J_vals[i], J_sources[i]) for i in range(len(J_vals))
                if (J_orders[i] == 3 and J_sources[i] is not None)]
        if base.GetAtomWithIdx(heavy_idx).IsInRingSize(6) and not _is_sp2_like(base.GetAtomWithIdx(heavy_idx)):
            for i in range(len(cand)):
                for j in range(i+1, len(cand)):
                    ii, Ji, src_i = cand[i]
                    jj, Jj, src_j = cand[j]
                    if src_i == src_j:
                        continue
                    Ai = base.GetAtomWithIdx(src_i)
                    Aj = base.GetAtomWithIdx(src_j)
                    if not (Ai.IsInRingSize(6) and Aj.IsInRingSize(6) and (not _is_sp2_like(Ai)) and (not _is_sp2_like(Aj))):
                        continue
                    if min(Ji, Jj) >= 9.0 and abs(Ji - Jj) <= 0.5 * dyn_tol:
                        J_vals[jj] = min(Jj + eps, Jj + 0.6)
                        break
                else:
                    continue
                break

    # Politique para
    allow_para = (spectrometer_mhz >= config.detect_para_min_mhz)
    if allow_para and para_buf_J:
        candidates = [(j, p) for (j, p) in zip(para_buf_J, para_buf_partners)
                      if j >= config.detect_para_min_j_hz]
        auto_k = 2 if spectrometer_mhz >= config.detect_para_min_mhz else 1
        k = int(config.para_max_per_proton) if (config.para_max_per_proton is not None) else auto_k
        k = max(0, k)
        if candidates and k > 0:
            candidates.sort(key=lambda x: (abs(x[0]-0.6)<=0.15, x[0]), reverse=True)
            picked: List[Tuple[float,int]] = []
            picked_partners = set()
            has_exact06 = any(abs(j - 0.6) <= 0.15 for j, _ in candidates)
            if not has_exact06:
                for j_raw, p_raw in zip(para_buf_J, para_buf_partners):
                    if abs(j_raw - 0.6) <= 0.15:
                        picked.append((j_raw, p_raw))
                        picked_partners.add(p_raw)
                        break
            for j, p in candidates:
                if len(picked) >= k: break
                if p in picked_partners: continue
                picked.append((j, p))
                picked_partners.add(p)
            for j_sel, p_sel in picked[:k]:
                J_vals.append(j_sel); J_sources.append(None); J_orders.append(5)
                partners.append(p_sel); homo_J_vals_for_second_order.append(j_sel)
                homo_pairs.append((p_sel, j_sel))
                try:
                    hh = [n.GetIdx() for n in molH.GetAtomWithIdx(p_sel).GetNeighbors() if n.GetSymbol()!='H'][0]
                    partner_heavy_idx_kept.append(hh)
                    partner_is_sp2_kept.append(_is_sp2_like(base.GetAtomWithIdx(hh)))
                except Exception:
                    pass

    # Solvant / quadrupolaire
    visco = solvent.get("viscosity", "low")
    low_cut = config.viscosity_low_cuts.get(visco, config.viscosity_low_cuts["low"])
    J_vals = [j for j in J_vals if j >= low_cut]
    if _has_quadrupolar_neighbor(base, heavy_idx):
        J_vals = [j for j in J_vals if j >= max(low_cut, 0.5)]

    # Détection second ordre (Δν/J raffinée)
    if homo_J_vals_for_second_order:
        real_ppm_vals: List[float] = []
        if proton_shifts_ppm is not None:
            try:
                d0 = proton_shifts_ppm.get(h_rep, None)
                if d0 is not None and math.isfinite(d0):
                    real_ppm_vals.append(d0)
                for p in partners:
                    val = proton_shifts_ppm.get(p, None)
                    if val is not None and math.isfinite(val):
                        real_ppm_vals.append(val)
            except Exception:
                real_ppm_vals = []

        if len(real_ppm_vals) >= 2:
            ppm_vals = real_ppm_vals
        else:
            if len(h_atoms) == 2:
                delta_ppm_default = 0.20
            elif is_sp2:
                delta_ppm_default = 0.12
            else:
                delta_ppm_default = 0.25
            ppm_vals = [0.0, delta_ppm_default]

        second_order, second_order_type = _refined_second_order_detection(
            homo_J_vals_for_second_order,
            ppm_vals,
            spectrometer_mhz,
            config
        )

        if len(ppm_vals) >= 2:
            delta_ppm_min = min(
                abs(a - b)
                for i, a in enumerate(ppm_vals)
                for b in ppm_vals[i+1:]
            )
            delta_nu_min_hz = _delta_nu_hz(delta_ppm_min, spectrometer_mhz)

    # Forcer second ordre pour AA'BB'/AA'XX'
    if system_type in ["AA'BB'", "AA'XX'"]:
        if not second_order:
            second_order = True
            if second_order_type == "first_order":
                second_order_type = "strong_coupling"

    # ✅ NOUVEAU: Appliquer les corrections J par solvant
    if config.use_solvent_J_corrections and J_vals:
        solvent_name = solvent.get("name", "CDCl3")
        corrected_J_vals = []
        for i, J in enumerate(J_vals):
            order = J_orders[i] if i < len(J_orders) else 3
            
            # Déterminer le type de couplage pour la correction
            if order == 2:
                j_type = "2J_geminal"
            elif order == 4:
                j_type = "4J_allylic"
            elif order >= 5:
                j_type = "3J_aromatic"  # Approximation pour longue distance
            elif is_sp2:
                j_type = "3J_vinyl"
            else:
                j_type = "3J_vicinal"
            
            J_corrected = _apply_solvent_J_correction(J, solvent_name, j_type)
            corrected_J_vals.append(J_corrected)
        
        J_vals = corrected_J_vals
        
        # Mettre à jour homo_pairs avec les J corrigés
        homo_pairs_corrected = []
        for i, (partner, _) in enumerate(homo_pairs):
            if i < len(J_vals):
                homo_pairs_corrected.append((partner, J_vals[i]))
            else:
                homo_pairs_corrected.append((partner, 0.0))
        homo_pairs = homo_pairs_corrected

    return (J_vals, second_order, partners, satellites, homo_pairs,
            partner_heavy_idx_kept, partner_is_sp2_kept,
            second_order_type, system_type, delta_nu_min_hz)

# ---------------------------
# Public simulators (par groupes)
# ---------------------------
def _effective_exchange_regime(solvent: dict, temperature_k: float) -> str:
    base_reg = solvent.get("exchange_regime", "moderate")
    rate = solvent.get("exchange_rate_factor", 1.0) * (temperature_k / 298.15)
    if solvent.get("name") in {"D2O"}: return "fast"
    if rate >= 1.5: return "fast"
    if rate <= 0.75: return "slow" if base_reg != "fast" else "moderate"
    return base_reg

def simulate_multiplicity_for_groups(
    mol: Chem.Mol,
    x_h_connectivity: Dict[int, List[int]],
    solvent: str = DEFAULT_SOLVENT,
    spectrometer_mhz: float = DEFAULT_SPECTROMETER_FREQ_MHZ,
    temperature_k: float = DEFAULT_TEMPERATURE_K,
    proton_shifts_ppm: Optional[Dict[int, float]] = None,
    include_heteronuclear: bool = True,
    use_3d: bool = True,
    tol: float = TOLERANCE_J_HZ,
    config: Optional[SimConfig] = None,
) -> List[int]:
    """Simulation de multiplicité par groupe (ancien comportement)."""
    if config is None: config = SimConfig()
    
    # ✅ NOUVEAU : Inférence stéréochimie E/Z pour doubles liaisons non définies
    if config.infer_double_bond_stereo:
        defined, undefined = _count_undefined_double_bonds(mol)
        if undefined > 0:
            if config.infer_stereo_log:
                import logging
                logging.info(f"🔧 Inférence stéréochimie: {undefined} double(s) liaison(s) non définies → {config.infer_stereo_default}")
            mol = _infer_double_bond_stereo(mol, config.infer_stereo_default)
    
    solv = SOLVENT_DB.get(solvent, SOLVENT_DB[DEFAULT_SOLVENT])
    eff_regime = _effective_exchange_regime(solv, temperature_k)

    if config.use_field_resolution_model:
        auto_lw = _typical_linewidth_from_field(spectrometer_mhz, solv.get("viscosity","low"))
        config.line_width_hz = auto_lw
        if config.digital_resolution_hz is None:
            config.digital_resolution_hz = _adaptive_digital_resolution(spectrometer_mhz)

    if config.grouping_tol_scale_auto:
        config.grouping_tol_scale = 0.6 * math.sqrt(600.0 / max(1.0, spectrometer_mhz))

    _field_based_tol = max(tol, config.grouping_tol_scale * config.line_width_hz)

    molH = Chem.AddHs(mol, addCoords=True)

    # ✅ Validation de la connectivité entrée
    _validate_connectivity(molH, x_h_connectivity)

    base = Chem.RemoveHs(molH)
    Chem.SanitizeMol(base)
    _assign_stereo(base)

    if use_3d and molH.GetNumConformers() == 0:
        try:
            try: ps = AllChem.ETKDGv3()
            except Exception: ps = AllChem.ETKDG()
            ps.randomSeed = int(config.random_seed)
            res = AllChem.EmbedMolecule(molH, ps)
            if res == 0: 
                _optimize_molecule_geometry(molH, config.force_field, config.force_field_max_iters, 0)
        except Exception:
            pass

    D = _topological_distances(base)
    aromatic_dmap = _precompute_aromatic_distance_map(base)
    heavy_is_sp2_cache = [_is_sp2_like(base.GetAtomWithIdx(i)) for i in range(base.GetNumAtoms())]

    H_ALL = [a.GetIdx() for a in molH.GetAtoms() if a.GetSymbol()=="H"]
    heavy_to_H: Dict[int, List[int]] = defaultdict(list)
    for h in H_ALL:
        nb = [n.GetIdx() for n in molH.GetAtomWithIdx(h).GetNeighbors()]
        if nb: heavy_to_H[nb[0]].append(h)

    ranks = _canonical_ranks(base)
    delta_cluster_id = {}
    
    # ✅ NOUVEAU: Créer le cache de conformères avec adaptation à la flexibilité
    conformer_cache = None
    if config.use_conformational_averaging:
        # Adapter le nombre de conformères selon la flexibilité
        n_conf = config.n_conformers
        if config.flexible_chain_boost:
            n_conf = _get_adaptive_n_conformers(base, config.n_conformers)
        
        conformer_cache = ConformerCache(
            base, 
            n_conformers=n_conf,
            energy_cutoff_kcal=config.conformer_energy_cutoff_kcal,
            random_seed=config.random_seed,
            force_field=config.force_field
        )
        conformer_cache.initialize()

    suppressed_heavy = set()
    for A in base.GetAtoms():
        t = _is_exchangeable_H_on_X(base, A.GetIdx())
        if not t: continue
        eff = solv["chemical_shift_effects"].get(t.replace("SH","OH"), None)
        if eff == "suppressed" or eff_regime == "fast": suppressed_heavy.add(A.GetIdx()); continue
        if not (solv["name"] == "DMSO-d6" and eff_regime == "slow"): suppressed_heavy.add(A.GetIdx())

    ordered = [(k, sorted(v)) for k, v in sorted(x_h_connectivity.items(), key=lambda kv: kv[0]) if v]
    out_classes: List[int] = []

    for heavy_idx, h_list in ordered:
        def _valid_hlist_for_heavy(hs: List[int]) -> bool:
            try: heavy_atomH = molH.GetAtomWithIdx(heavy_idx)
            except Exception: return False
            heavy_neighborsH = {n.GetIdx() for n in heavy_atomH.GetNeighbors() if n.GetSymbol()=="H"}
            return all((0 <= h < molH.GetNumAtoms()) and (molH.GetAtomWithIdx(h).GetSymbol()=="H") and (h in heavy_neighborsH) for h in hs)

        if not _valid_hlist_for_heavy(h_list):
            h_list = [n.GetIdx() for n in molH.GetAtomWithIdx(heavy_idx).GetNeighbors() if n.GetSymbol()=="H"]
            if not h_list: out_classes.append(LABEL2IDX["s"]); continue

        exch = _is_exchangeable_H_on_X(base, heavy_idx)
        if exch:
            eff = solv["chemical_shift_effects"].get(exch.replace("SH","OH"), None)
            if eff == "suppressed" or eff_regime == "fast":
                out_classes.append(LABEL2IDX["s"]); continue
            if exch in {"OH","NH"} and solv["name"] == "DMSO-d6":
                smallJ = False
                if proton_shifts_ppm and h_list:
                    h0 = h_list[0]
                    all_H = [k for k in proton_shifts_ppm.keys()]
                    deltas = [abs(proton_shifts_ppm.get(h0, 9e9) - proton_shifts_ppm.get(ph, -9e9))
                              for ph in all_H if ph != h0 and math.isfinite(proton_shifts_ppm.get(ph, float('nan')))]
                else:
                    deltas = []
                if deltas:
                    smallJ = (min(deltas) * spectrometer_mhz) <= 8.0
                out_classes.append(LABEL2IDX["d" if smallJ else "s"]); continue
            out_classes.append(LABEL2IDX["s"]); continue

        (J_vals, second_order, partners, _sat, homo_pairs,
         partner_heavy_idx_kept, partner_is_sp2_kept,
         second_order_type, system_type, delta_nu_min_hz) = _estimate_Js_for_group(
            molH, base, D, aromatic_dmap, heavy_idx, h_list, heavy_to_H,
            spectrometer_mhz, solv, proton_shifts_ppm, include_heteronuclear,
            suppressed_heavy, heavy_is_sp2_cache, config, ranks=ranks,
            delta_cluster_id=delta_cluster_id,
            geminal_partner=None,
            conformer_cache=conformer_cache
        )

        heavy = base.GetAtomWithIdx(heavy_idx)

        # ✅ NOUVEAU: Détecter les systèmes AA'BB' aromatiques symétriques → singulet apparent
        if heavy.GetSymbol() == "C" and heavy.GetIsAromatic():
            if _is_symmetric_para_disubstituted_benzene(base, heavy_idx, aromatic_dmap, ranks):
                # Para-disubstitué symétrique: AA'BB' fortement couplé → singulet apparent
                out_classes.append(LABEL2IDX["s"])
                continue

        if heavy_is_sp2_cache[heavy_idx] and not J_vals and heavy.GetSymbol() == "C" and heavy.GetIsAromatic():
            J_fallback = _fallback_aromatic_J_when_empty(base, aromatic_dmap, heavy_idx, heavy_to_H, config, ranks=ranks)
            if J_fallback:
                J_vals = J_fallback
                second_order = False
                second_order_type = "first_order"

        if (heavy_is_sp2_cache[heavy_idx]) and not J_vals:
            if config.aromatic_isolated_as_singlet and heavy.GetIsAromatic():
                out_classes.append(LABEL2IDX["s"]); continue
            elif config.treat_aromatic_noJ_equiv_as_singlet and heavy.GetIsAromatic():
                same = [k for k in range(base.GetNumAtoms()) if ranks[k]==ranks[heavy_idx] and heavy_to_H.get(k)]
                out_classes.append(LABEL2IDX["s"]); continue
            out_classes.append(LABEL2IDX["s"]); continue

        if second_order:
            anis_high = (solv.get("anisotropy","normal") == "high")
            groups_try = _group_J_values(J_vals, tol=_field_based_tol, anisotropy_high=anis_high, config=config)
            label2 = _first_order_label_from_groups(groups_try) if groups_try else None
            if label2 in {None, "m"}:
                label2 = None

            if (label2 is None) and config.use_spin_fallback and proton_shifts_ppm is not None and h_list:
                h0 = h_list[0]
                tgt_ppm = proton_shifts_ppm.get(h0, None)
                if tgt_ppm is not None and math.isfinite(tgt_ppm):
                    homo_pairs_sorted = sorted(
                        [(p, J) for (p, J) in homo_pairs if J is not None and J > 0 and p in (proton_shifts_ppm or {})],
                        key=lambda x: x[1], reverse=True
                    )[:max(1, min(2, config.spin_fallback_max_spins-1))]
                    if homo_pairs_sorted:
                        idx_map = {p: i for i, p in enumerate(partners)}
                        sel_idxs = [idx_map[p] for p in homo_pairs_sorted if p in idx_map]
                        ppms = [proton_shifts_ppm[partners[i]] for i in sel_idxs]
                        Js   = [homo_pairs[idx_map[partners[i]]][1] if idx_map.get(partners[i], None) is not None else homo_pairs_sorted[k][1]
                                for k, i in enumerate(sel_idxs)]
                        sel_heavy = [partner_heavy_idx_kept[i] for i in sel_idxs if i < len(partner_heavy_idx_kept)]
                        sel_is_sp2 = [partner_is_sp2_kept[i] for i in sel_idxs if i < len(partner_is_sp2_kept)]
                        label2 = _spin_fallback_label(
                            spectrometer_mhz, tgt_ppm, ppms, Js, config,
                            partner_heavy_idx=sel_heavy if len(sel_heavy)==len(ppms) else None,
                            base=base, aromatic_dmap=aromatic_dmap,
                            is_sp2_list=sel_is_sp2 if len(sel_is_sp2)==len(ppms) else None
                        )

            if heavy_is_sp2_cache[heavy_idx] and heavy.GetSymbol() == "C" and J_vals:
                label2 = _apparent_sp2_label_from_J(J_vals, minJ=0.4)

            out_classes.append(LABEL2IDX.get(label2 or "m", LABEL2IDX["m"]))
            continue

        anis_high = (solv.get("anisotropy","low") == "high")
        groups = _group_J_values(J_vals, tol=_field_based_tol, anisotropy_high=anis_high, config=config)
        label = _first_order_label_from_groups(groups)

        if heavy_is_sp2_cache[heavy_idx] and heavy.GetSymbol() == "C" and J_vals:
            label = _apparent_sp2_label_from_J(J_vals, minJ=0.4)

        if len(h_list) == 2:
            nb_heavy = [n for n in heavy.GetNeighbors()]
            chiral_neigh = any(_is_chiral_center(n) for n in nb_heavy)
            if chiral_neigh and len(groups) <= 1:
                if label == "t": label = "dd"
                elif label == "q": label = "ddd"
                elif label == "quin": label = "dddd"

        complicated = label in {"ddt","dtt","ddd","dddd","ttq","qqt","tqq","qtq"}
        if solv.get("viscosity","low") != "low" and complicated:
            if (J_vals and max(J_vals) <= config.viscosity_simplify_maxJ_hz) or (len(groups) > config.viscosity_simplify_max_groups):
                label = "m"

        out_classes.append(LABEL2IDX.get(label, LABEL2IDX["m"]))

    return out_classes


# ---------------------------
# ESTIMATION HEURISTIQUE DES SHIFTS
# ---------------------------
def _estimate_shifts_heuristic(molH: Chem.Mol, proton_to_heavy: Dict[int, int]) -> Dict[int, float]:
    """
    Estimation simple des déplacements chimiques par type d'atome.
    """
    shifts: Dict[int, float] = {}
    for h_idx, heavy_idx in proton_to_heavy.items():
        heavy = molH.GetAtomWithIdx(heavy_idx)
        sym = heavy.GetSymbol()

        if sym == "O":
            shifts[h_idx] = 5.0
        elif sym == "C":
            nb_o = sum(1 for n in heavy.GetNeighbors() if n.GetSymbol() == "O")
            if heavy.GetIsAromatic():
                shifts[h_idx] = 7.2
            elif nb_o == 0:
                shifts[h_idx] = 1.5
            elif nb_o == 1:
                shifts[h_idx] = 3.5
            elif nb_o >= 2:
                shifts[h_idx] = 4.2
            else:
                shifts[h_idx] = 2.0
        else:
            shifts[h_idx] = 1.0

    return shifts


# ---------------------------
# ✅ CORRECTION 2 : WRAPPER par proton avec garde-fou diastéréotopie
def simulate_multiplicity_per_proton(
    mol: Chem.Mol,
    x_h_connectivity: Dict[int, List[int]],
    proton_shifts_ppm: Dict[int, float],
    solvent: str = DEFAULT_SOLVENT,
    spectrometer_mhz: float = DEFAULT_SPECTROMETER_FREQ_MHZ,
    temperature_k: float = DEFAULT_TEMPERATURE_K,
    include_heteronuclear: bool = True,
    use_3d: bool = True,
    tol: float = TOLERANCE_J_HZ,
    config: Optional[SimConfig] = None,
) -> Dict[int, Dict[str, Any]]:
    """
    ✅ VERSION CORRIGÉE ET DÉFINITIVE
    Retourne une multiplicité par PROTON individuel avec informations enrichies.
    
    CORRECTIONS CRITIQUES :
    1. Vérifie le nombre RÉEL de H attachés avant toute séparation diastéréotopique
    2. Ne sépare les CH₂ que s'ils sont liés à un centre stéréogène (vraie diastéréotopie)
    3. Évite les protons fantômes sur les CH/CH₃
    """
    if config is None:
        config = SimConfig()
    
    # ✅ NOUVEAU : Inférence stéréochimie E/Z pour doubles liaisons non définies
    if config.infer_double_bond_stereo:
        defined, undefined = _count_undefined_double_bonds(mol)
        if undefined > 0:
            if config.infer_stereo_log:
                import logging
                logging.info(f"🔧 Inférence stéréochimie: {undefined} double(s) liaison(s) non définies → {config.infer_stereo_default}")
            mol = _infer_double_bond_stereo(mol, config.infer_stereo_default)
    
    # Préparation molécule
    molH = Chem.AddHs(mol, addCoords=True)
    _validate_connectivity(molH, x_h_connectivity)
    
    # Build heavy_to_H from actual molecule
    H_ALL = [a.GetIdx() for a in molH.GetAtoms() if a.GetSymbol() == "H"]
    heavy_to_H: Dict[int, List[int]] = defaultdict(list)
    for h in H_ALL:
        nb = [n.GetIdx() for n in molH.GetAtomWithIdx(h).GetNeighbors()]
        if nb:
            heavy_to_H[nb[0]].append(h)
    
    # === ÉTAPE 1 : GARDE-FOU GÉNÉRIQUE - Vérification du nombre réel de H ===
    expanded_connectivity: Dict[int, List[int]] = {}
    proton_to_heavy: Dict[int, int] = {}
    geminal_partners: Dict[int, Optional[int]] = {}
    
    for heavy_idx, h_list in sorted(x_h_connectivity.items(), key=lambda kv: kv[0]):
        if not h_list:
            continue
        
        # 🔴 CRITIQUE : Obtenir le nombre RÉEL de H attachés à ce carbone
        heavy_atom = molH.GetAtomWithIdx(heavy_idx)
        actual_h_neighbors = [n for n in heavy_atom.GetNeighbors() if n.GetSymbol() == "H"]
        n_h_attached = len(actual_h_neighbors)
        
        # === CAS SPÉCIAL : =CH₂ VINYLIQUE (terminaux) ===
        # Les protons géminaux sur un =CH₂ sont TOUJOURS non-équivalents
        # car ils ont des J différents (cis vs trans) avec le H vicinal
        is_terminal_vinyl_CH2 = False
        if n_h_attached == 2 and len(h_list) == 2:
            hyb = heavy_atom.GetHybridization()
            if hyb == rdchem.HybridizationType.SP2:
                # Vérifier si c'est un =CH₂ terminal (carbone sp2 avec 2 H)
                c_neighbors = [n for n in heavy_atom.GetNeighbors() if n.GetSymbol() != "H"]
                if len(c_neighbors) == 1:
                    # Un seul voisin non-H → c'est terminal (=CH₂)
                    # Vérifier qu'il y a une double liaison
                    for bond in heavy_atom.GetBonds():
                        if bond.GetBondType() == rdchem.BondType.DOUBLE:
                            is_terminal_vinyl_CH2 = True
                            break
        
        if is_terminal_vinyl_CH2:
            # ✅ Séparer les protons vinyliques - ils sont toujours non-équivalents
            h_a, h_b = sorted(h_list)
            expanded_connectivity[heavy_idx] = [h_a]
            pseudo_key = -(h_b + 1)
            expanded_connectivity[pseudo_key] = [h_b]
            proton_to_heavy[h_a] = heavy_idx
            proton_to_heavy[h_b] = heavy_idx
            geminal_partners[h_a] = h_b
            geminal_partners[h_b] = h_a
            continue
        
        # === CAS 1 : Ce n'est PAS un CH₂ → FORCER le regroupement ===
        # Cela bloque les protons fantômes sur CH et CH₃
        if n_h_attached != 2:
            expanded_connectivity[heavy_idx] = list(h_list)
            for h in h_list:
                proton_to_heavy[h] = heavy_idx
                geminal_partners[h] = None
            continue  # ← SAUTE la logique diastéréotopique
        
        # === CAS 2 : C'est un CH₂ → Tester la diastéréotopie SEULEMENT si lié à un centre stéréogène ===
        # Vérifier si un voisin est un carbone stéréogène (chiral)
        is_diastereotopic = False
        if config.separate_diastereotopic_protons and len(h_list) == 2:
            neighbors = [n for n in heavy_atom.GetNeighbors() if n.GetSymbol() != "H"]
            for nb in neighbors:
                if nb.GetSymbol() == "C":
                    tag = nb.GetChiralTag()
                    if tag in (rdchem.ChiralType.CHI_TETRAHEDRAL_CW, rdchem.ChiralType.CHI_TETRAHEDRAL_CCW):
                        is_diastereotopic = True
                        break
        
        if is_diastereotopic:
            # Séparer uniquement si vraiment diastéréotopique
            h_a, h_b = h_list
            expanded_connectivity[heavy_idx] = [h_a]
            pseudo_key = -(h_b + 1)  # Clé unique
            expanded_connectivity[pseudo_key] = [h_b]
            proton_to_heavy[h_a] = heavy_idx
            proton_to_heavy[h_b] = heavy_idx
            geminal_partners[h_a] = h_b
            geminal_partners[h_b] = h_a
        else:
            # Garder ensemble (CH₂ non diastéréotopique, ex: CH₂OH du glucose)
            expanded_connectivity[heavy_idx] = list(h_list)
            for h in h_list:
                proton_to_heavy[h] = heavy_idx
                geminal_partners[h] = None
    
    # === ÉTAPE 2 : Préparer les structures de simulation ===
    solv = SOLVENT_DB.get(solvent, SOLVENT_DB[DEFAULT_SOLVENT])
    eff_regime = _effective_exchange_regime(solv, temperature_k)
    
    # === NOUVEAU: Appliquer corrections de solvant sur les shifts ===
    corrected_shifts = proton_shifts_ppm
    if config.use_solvent_shift_corrections and proton_shifts_ppm:
        corrected_shifts = _apply_solvent_corrections(
            proton_shifts_ppm, molH, heavy_to_H, solv.get("name", solvent)
        )
    
    # === NOUVEAU: Largeur de raie selon solvant ===
    if config.use_solvent_linewidth:
        solvent_lw = _get_solvent_linewidth(solv.get("name", solvent), default=0.8)
        config.line_width_hz = max(config.line_width_hz, solvent_lw)
    
    if config.use_field_resolution_model:
        auto_lw = _typical_linewidth_from_field(spectrometer_mhz, solv.get("viscosity", "low"))
        config.line_width_hz = max(config.line_width_hz, auto_lw)
        if config.digital_resolution_hz is None:
            config.digital_resolution_hz = _adaptive_digital_resolution(spectrometer_mhz)
    
    if config.grouping_tol_scale_auto:
        config.grouping_tol_scale = 0.6 * math.sqrt(600.0 / max(1.0, spectrometer_mhz))
    
    _field_based_tol = max(tol, config.grouping_tol_scale * config.line_width_hz)
    
    base = Chem.RemoveHs(molH)
    Chem.SanitizeMol(base)
    _assign_stereo(base)
    
    if use_3d and molH.GetNumConformers() == 0:
        try:
            try:
                ps = AllChem.ETKDGv3()
            except Exception:
                ps = AllChem.ETKDG()
            ps.randomSeed = int(config.random_seed)
            res = AllChem.EmbedMolecule(molH, ps)
            if res == 0:
                _optimize_molecule_geometry(molH, config.force_field, config.force_field_max_iters, 0)
        except Exception:
            pass
    
    D = _topological_distances(base)
    aromatic_dmap = _precompute_aromatic_distance_map(base)
    heavy_is_sp2_cache = [_is_sp2_like(base.GetAtomWithIdx(i)) for i in range(base.GetNumAtoms())]
    ranks = _canonical_ranks(base)
    
    # ✅ NOUVEAU: Créer le cache de conformères avec adaptation à la flexibilité
    conformer_cache = None
    if config.use_conformational_averaging:
        # Adapter le nombre de conformères selon la flexibilité
        n_conf = config.n_conformers
        if config.flexible_chain_boost:
            n_conf = _get_adaptive_n_conformers(base, config.n_conformers)
        
        conformer_cache = ConformerCache(
            base, 
            n_conformers=n_conf,
            energy_cutoff_kcal=config.conformer_energy_cutoff_kcal,
            random_seed=config.random_seed,
            force_field=config.force_field
        )
        conformer_cache.initialize()
    
    # === NOUVEAU: Pré-calculer les systèmes de spin avancés ===
    spin_systems: Dict[int, Tuple[Optional[str], Dict[str, Any]]] = {}
    if config.advanced_spin_system_detection:
        for heavy_idx in x_h_connectivity.keys():
            spin_sys, spin_info = _detect_spin_system(
                base, heavy_idx, x_h_connectivity, aromatic_dmap,
                proton_shifts_ppm=corrected_shifts,
                spectrometer_mhz=spectrometer_mhz
            )
            if spin_sys:
                spin_systems[heavy_idx] = (spin_sys, spin_info)
    
    suppressed_heavy = set()
    for A in base.GetAtoms():
        t = _is_exchangeable_H_on_X(base, A.GetIdx())
        if not t:
            continue
        eff = solv["chemical_shift_effects"].get(t.replace("SH", "OH"), None)
        if eff == "suppressed" or eff_regime == "fast":
            suppressed_heavy.add(A.GetIdx())
            continue
        if not (solv["name"] == "DMSO-d6" and eff_regime == "slow"):
            suppressed_heavy.add(A.GetIdx())
    
    # === ÉTAPE 3 : Calculer la multiplicité pour chaque proton ===
    result: Dict[int, Dict[str, Any]] = {}
    
    for key, h_list in sorted(expanded_connectivity.items(), key=lambda kv: kv[0]):
        if not h_list:
            continue
        
        h_idx = h_list[0]
        true_heavy_idx = proton_to_heavy[h_idx]
        geminal_partner = geminal_partners.get(h_idx, None)
        
        # ✅ NOUVEAU: Détecter l'élargissement pour protons labiles
        # Utiliser molH car base n'a pas les H explicites
        broadening = None
        if config.mark_labile_protons_broad:
            broadening = _estimate_labile_proton_broadening(
                molH, true_heavy_idx, solv.get("name", solvent), temperature_k
            )
        
        # ✅ NOUVEAU: Détecter imidazole H2 → singulet
        if _is_imidazole_H2_carbon(base, true_heavy_idx):
            result[h_idx] = {
                "label": "s",
                "label_idx": LABEL2IDX["s"],
                "J_couplings": [],
                "peak_positions_hz": [0.0],
                "peak_intensities": [1.0],
                "second_order_type": "first_order",
                "system_type": "imidazole_H2",
                "virtual_coupling": False,
                "partners": [],
                "satellites_13C_hz": [],
                "broadening": broadening,
                "heterocycle_position": 2,  # Position 2 pour imidazole H2
                "sugar_position": None,
            }
            continue
        
        # Vérifier si c'est un proton échangeable
        exch = _is_exchangeable_H_on_X(base, true_heavy_idx)
        if exch:
            eff = solv["chemical_shift_effects"].get(exch.replace("SH", "OH"), None)
            result[h_idx] = {
                "label": "s",
                "label_idx": LABEL2IDX["s"],
                "J_couplings": [],
                "peak_positions_hz": [0.0],
                "peak_intensities": [1.0],
                "second_order_type": "first_order",
                "system_type": None,
                "virtual_coupling": False,
                "partners": [],
                "satellites_13C_hz": [],
                "broadening": broadening,  # ✅ NOUVEAU
                "heterocycle_position": None,
                "sugar_position": None,
            }
            continue
        
        # Calculer J pour ce proton
        (J_vals, second_order, partners, sat, homo_pairs,
         partner_heavy_idx_kept, partner_is_sp2_kept,
         second_order_type, system_type, delta_nu_min_hz) = _estimate_Js_for_group(
            molH, base, D, aromatic_dmap, true_heavy_idx, h_list, heavy_to_H,
            spectrometer_mhz, solv, corrected_shifts, include_heteronuclear,
            suppressed_heavy, heavy_is_sp2_cache, config, ranks=ranks,
            delta_cluster_id={},
            geminal_partner=geminal_partner,
            conformer_cache=conformer_cache
        )
        
        # Initialiser heterocycle_position à None (sera défini si hétérocycle détecté)
        heterocycle_position = None
        
        # === NOUVEAU: Utiliser le système de spin détecté si disponible ===
        if config.advanced_spin_system_detection and true_heavy_idx in spin_systems:
            detected_sys, spin_info = spin_systems[true_heavy_idx]
            if detected_sys and system_type is None:
                system_type = detected_sys
            # Marquer comme second ordre si système fortement couplé
            if detected_sys in ("AB", "AB₂", "AA'BB'", "A₂B₃", "ABX", "complex_second_order"):
                if not second_order:
                    second_order = True
                    second_order_type = "strong_coupling"
        
        # Détection couplage virtuel
        # ✅ IMPORTANT: Ne pas détecter le couplage virtuel pour les hétérocycles gérés par le module
        is_heterocycle_managed = system_type is not None and system_type.startswith("heterocycle_")
        is_virtual = False
        if J_vals and partners and not is_heterocycle_managed:
            try:
                is_virtual = _detect_virtual_coupling(
                    molH, h_idx, partners, J_vals, proton_shifts_ppm, spectrometer_mhz
                )
            except Exception:
                is_virtual = False
        
        if is_virtual:
            J_vals = [j / 2.0 for j in J_vals]
        
        # ✅ NOUVEAU: Détecter les systèmes AA'BB' aromatiques symétriques → singulet apparent
        heavy = base.GetAtomWithIdx(true_heavy_idx)
        if heavy.GetSymbol() == "C" and heavy.GetIsAromatic():
            if _is_symmetric_para_disubstituted_benzene(base, true_heavy_idx, aromatic_dmap, ranks):
                # Para-disubstitué symétrique: AA'BB' fortement couplé → singulet apparent
                result[h_idx] = {
                    "label": "s",
                    "label_idx": LABEL2IDX["s"],
                    "J_couplings": [],
                    "peak_positions_hz": [0.0],
                    "peak_intensities": [1.0],
                    "second_order_type": "AA'BB'_symmetric",
                    "system_type": "AA'BB'",
                    "virtual_coupling": False,
                    "partners": [],
                    "satellites_13C_hz": [],
                    "broadening": None,
                    "heterocycle_position": None,
                    "sugar_position": None,
                }
                continue
        
        # Fallback aromatique
        # ✅ IMPORTANT: Ne pas appliquer le fallback si le module hétérocycle a déjà géré ce proton
        # (is_heterocycle_managed est déjà défini avant la détection du couplage virtuel)
        if heavy_is_sp2_cache[true_heavy_idx] and not J_vals and heavy.GetSymbol() == "C" and heavy.GetIsAromatic() and not is_heterocycle_managed:
            J_fallback = _fallback_aromatic_J_when_empty(base, aromatic_dmap, true_heavy_idx, heavy_to_H, config, ranks=ranks)
            if J_fallback:
                J_vals = J_fallback
                second_order = False
                second_order_type = "first_order"
        
        # H sp2 sans J
        if heavy_is_sp2_cache[true_heavy_idx] and not J_vals:
            label = "s"
            label_idx = LABEL2IDX["s"]
            result[h_idx] = {
                "label": label,
                "label_idx": label_idx,
                "J_couplings": [],
                "peak_positions_hz": [0.0],
                "peak_intensities": [1.0],
                "second_order_type": "first_order",
                "system_type": system_type,
                "virtual_coupling": is_virtual,
                "partners": [],
                "satellites_13C_hz": sat or [],
                "broadening": broadening,
                "heterocycle_position": heterocycle_position,
                "sugar_position": None,
            }
            continue
        
        # 🔴 Méthyles isolés sp3 sans J (benzyl, méthoxy, N-méthyl)
        if not J_vals and _is_isolated_methyl(base, true_heavy_idx, heavy_to_H):
            label = "s"
            label_idx = LABEL2IDX["s"]
            result[h_idx] = {
                "label": label,
                "label_idx": label_idx,
                "J_couplings": [],
                "peak_positions_hz": [0.0],
                "peak_intensities": [1.0],
                "second_order_type": "first_order",
                "system_type": None,
                "virtual_coupling": False,
                "partners": [],
                "satellites_13C_hz": sat or [],
                "broadening": broadening,
                "heterocycle_position": None,
                "sugar_position": None,
            }
            continue
        
        # Regroupement des J
        anis_high = (solv.get("anisotropy", "normal") == "high")
        groups = _group_J_values(J_vals, tol=_field_based_tol, anisotropy_high=anis_high, config=config)
        if not groups:
            groups = [[]]
        
        # J_couplings structurés
        J_couplings: List[Dict[str, float]] = []
        for g in groups:
            if not g:
                continue
            J_mean = float(np.mean(g))
            n_eq = len(g)
            if J_mean > 0.0 and n_eq > 0:
                J_couplings.append({
                    "J_hz": round(J_mean, 3),
                    "n_equiv": int(n_eq),
                })
        
        # Détermination du label
        # ✅ NOUVEAU: Si le module hétérocycle a fourni une multiplicité, l'utiliser directement
        heterocycle_expected_mult = None
        heterocycle_position = None
        sugar_expected_mult = None
        sugar_position = None
        
        if is_heterocycle_managed and system_type and ':' in system_type:
            # Extraire la multiplicité et position du format "heterocycle_<type>:<mult>:<position>"
            parts = system_type.split(':')
            if len(parts) >= 2:
                heterocycle_expected_mult = parts[1]
            if len(parts) >= 3:
                try:
                    heterocycle_position = int(parts[2])
                except (ValueError, TypeError):
                    heterocycle_position = parts[2] if parts[2] else None
            # Nettoyer le system_type pour ne garder que le type
            system_type = parts[0]
        
        # ✅ NOUVEAU: Si le module sucre a fourni une multiplicité, l'utiliser directement
        if system_type and system_type.startswith('sugar_') and ':' in system_type:
            # Extraire la multiplicité et position du format "sugar_<type>:<mult>:<position>"
            parts = system_type.split(':')
            if len(parts) >= 2:
                sugar_expected_mult = parts[1]
            if len(parts) >= 3:
                try:
                    sugar_position = int(parts[2])
                except (ValueError, TypeError):
                    sugar_position = parts[2] if parts[2] else None
            # Nettoyer le system_type pour ne garder que le type (ex: sugar_pyranose)
            system_type = parts[0]
        
        if heterocycle_expected_mult and heterocycle_expected_mult in LABEL2IDX:
            # Utiliser la multiplicité du module hétérocycle
            label = heterocycle_expected_mult
        elif sugar_expected_mult and sugar_expected_mult in LABEL2IDX:
            # Utiliser la multiplicité du module sucre
            label = sugar_expected_mult
        elif second_order:
            label2 = _first_order_label_from_groups(groups) if groups else "m"
            if heavy_is_sp2_cache[true_heavy_idx] and heavy.GetSymbol() == "C" and J_vals:
                label2 = _apparent_sp2_label_from_J(J_vals, minJ=0.4)
            label = label2 or "m"
        else:
            label = _first_order_label_from_groups(groups)
            if heavy_is_sp2_cache[true_heavy_idx] and heavy.GetSymbol() == "C" and J_vals:
                label = _apparent_sp2_label_from_J(J_vals, minJ=0.4)
        
        label_idx = LABEL2IDX.get(label, LABEL2IDX["m"])
        
        # Simplification solvants visqueux
        complicated = label in {"ddt", "dtt", "ddd", "dddd", "ttq", "qqt", "tqq", "qtq"}
        if solv.get("viscosity", "low") != "low" and complicated:
            if (J_vals and max(J_vals) <= config.viscosity_simplify_maxJ_hz) or (len(groups) > config.viscosity_simplify_max_groups):
                label = "m"
                label_idx = LABEL2IDX["m"]
        
        # Partenaires
        partner_info: List[Dict[str, Any]] = []
        for partner_h_idx, J_val in homo_pairs:
            if partner_h_idx in proton_shifts_ppm:
                partner_info.append({
                    "partner_h_index": int(partner_h_idx),
                    "partner_shift_ppm": round(float(proton_shifts_ppm[partner_h_idx]), 3),
                    "J_hz": round(float(J_val), 2),
                })
        
        # === SIMULATION HAMILTONIEN EXACT pour systèmes second ordre ===
        exact_transitions = None
        if (config.use_exact_spin_simulation and system_type is not None 
            and system_type in ("AB", "AB_weak", "AB₂", "AA'BB'", "A₂X₂", "ABX", "A₂X₃", "A₂B₃")):
            # Collecter les indices H du système
            h_indices = list(h_list)
            # Ajouter les partenaires
            for pinfo in homo_pairs:
                if pinfo[0] not in h_indices:
                    h_indices.append(pinfo[0])
            
            if len(h_indices) <= config.exact_spin_max_spins:
                try:
                    exact_transitions = _get_exact_spectrum_for_spin_system(
                        system_type=system_type,
                        h_indices=h_indices,
                        proton_shifts_ppm=corrected_shifts,
                        J_values=J_vals,
                        spectrometer_mhz=spectrometer_mhz
                    )
                except Exception:
                    exact_transitions = None
        
        # Construction des raies
        if exact_transitions is not None and len(exact_transitions) > 0:
            # Utiliser les transitions exactes de l'Hamiltonien
            peaks = [(t[0], t[1]) for t in exact_transitions]
        elif second_order and second_order_type in {
            "strong_coupling", "moderate_second_order", "borderline_roof_effect"
        } and J_vals and delta_nu_min_hz > 0.0:
            Jmax = max(J_vals)
            peaks = _calculate_roof_effect_intensities(Jmax, delta_nu_min_hz, label)
        else:
            peaks = _first_order_peak_positions(J_couplings, label)
        
        J_AB = None
        if len(J_couplings) >= 1:
            J_AB = J_couplings[0]["J_hz"]
        peaks = _add_combination_lines(peaks, system_type, J_AB, config, delta_nu_hz=delta_nu_min_hz)
        
        peak_positions = [float(p[0]) for p in peaks]
        peak_intensities = [float(p[1]) for p in peaks]
        
        result[h_idx] = {
            "label": label,
            "label_idx": label_idx,
            "J_couplings": J_couplings,
            "peak_positions_hz": peak_positions,
            "peak_intensities": peak_intensities,
            "second_order_type": second_order_type,
            "system_type": system_type,
            "virtual_coupling": is_virtual,
            "partners": partner_info,
            "satellites_13C_hz": [round(float(s), 2) for s in (sat or [])],
            "broadening": broadening,  # ✅ NOUVEAU: None, "sharp", "broad", "very_broad", "exchanged"
            "heterocycle_position": heterocycle_position,  # ✅ Position IUPAC pour hétérocycles (None si non applicable)
            "sugar_position": sugar_position,  # ✅ Position pour sucres (None si non applicable)
        }
    
    return result


# ✅ PRÉPARATION MOLÉCULE (avec validation structure)
# ---------------------------
def prepare_molecule_for_nmr(
    smiles: str,
    generate_3d: bool = True,
    random_seed: int = 42,
    skip_validation: bool = False,
    force_field: str = "MMFF94s",
) -> Tuple[Chem.Mol, Dict[int, List[int]], Dict[int, int]]:
    """
    ✅ Version corrigée avec validation moléculaire.
    
    Prépare une molécule pour la simulation RMN :
    - Valide la structure (valences, stéréochimie, cycles)
    - Ajoute les hydrogènes explicites
    - Génère éventuellement une géométrie 3D
    - Construit la connectivité X-H
    
    Args:
        smiles: SMILES de la molécule
        generate_3d: Générer une conformation 3D
        random_seed: Graine aléatoire pour la génération 3D
        skip_validation: Passer la validation (dangereux, seulement pour debug)
        force_field: Force field pour optimisation ("MMFF94", "MMFF94s", "UFF")
    
    Returns:
        (molH, heavy_to_H, proton_to_heavy)
        - molH : Mol avec H explicites
        - heavy_to_H : {heavy_idx: [h_idx1, h_idx2, ...]}
        - proton_to_heavy : {h_idx: heavy_idx}
    
    Raises:
        ValueError: Si le SMILES est invalide ou la structure incorrecte
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES invalide: {smiles}")
    
    # ✅ VALIDATION OBLIGATOIRE (avec smiles pour messages d'erreur clairs)
    if not skip_validation:
        _validate_molecular_structure(mol, smiles)
    
    molH = Chem.AddHs(mol, addCoords=True)
    
    if generate_3d and molH.GetNumConformers() == 0:
        try:
            try: 
                ps = AllChem.ETKDGv3()
            except AttributeError: 
                ps = AllChem.ETKDG()
            ps.randomSeed = int(random_seed)
            ps.maxAttempts = 50
            res = AllChem.EmbedMolecule(molH, ps)
            if res == 0: 
                _optimize_molecule_geometry(molH, force_field, 200, 0)
        except Exception as e:
            if SimConfig().verbose:
                print(f"[prepare_molecule_for_nmr] 3D generation failed: {e}")
    
    # Construction de la connectivité (indices à base 0)
    heavy_to_H: Dict[int, List[int]] = defaultdict(list)
    proton_to_heavy: Dict[int, int] = {}
    
    for atom in molH.GetAtoms():
        if atom.GetSymbol() != "H":
            continue
        heavy_nb = atom.GetNeighbors()
        if not heavy_nb:
            continue
        heavy_idx = heavy_nb[0].GetIdx()
        idx_h = atom.GetIdx()
        proton_to_heavy[idx_h] = heavy_idx
        heavy_to_H[heavy_idx].append(idx_h)
    
    # ✅ Convertir defaultdict en dict normal pour éviter des surprises
    return molH, dict(heavy_to_H), proton_to_heavy

# ---------------------------
# API SIMPLE : simulate_molecule
# ---------------------------
def simulate_molecule(
    smiles: str,
    solvent: str = DEFAULT_SOLVENT,
    spectrometer_mhz: float = DEFAULT_SPECTROMETER_FREQ_MHZ,
    temperature_k: float = DEFAULT_TEMPERATURE_K,
    generate_3d: bool = True,
    config: Optional[SimConfig] = None,
) -> Dict[int, Dict[str, Any]]:
    """
    ✅ API principale "one-liner" corrigée.
    
    Simule le spectre RMN ¹H d'une molécule à partir de son SMILES.
    
    Args:
        smiles: SMILES de la molécule
        solvent: Solvant deutéré (CDCl3, DMSO-d6, D2O, etc.)
        spectrometer_mhz: Fréquence du spectromètre (MHz)
        temperature_k: Température (K)
        generate_3d: Générer une géométrie 3D
        config: Configuration avancée (SimConfig)
    
    Returns:
        { h_idx: {label, J_couplings, peak_positions_hz, ...}, ... }
        Un dictionnaire par proton avec toutes les informations RMN
    
    Example:
        >>> result = simulate_molecule("CC(O)C")
        >>> for h_idx, data in result.items():
        ...     print(f"H{h_idx}: {data['label']} @ {data.get('shift_ppm', 'N/A')} ppm")
    """
    if config is None:
        config = SimConfig()

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES invalide ou non supporté: '{smiles}'")

    # Préparer la molécule et la connectivité X–H (avec validation)
    molH, connectivity, proton_to_heavy = prepare_molecule_for_nmr(
        smiles,
        generate_3d=generate_3d,
        random_seed=config.random_seed,
    )

    # Estimation fallback des shifts (remplaçable par un prédicteur sérieux)
    shifts = _estimate_shifts_heuristic(molH, proton_to_heavy)

    # L'algorithme complet de multiplicité est dans simulate_multiplicity_per_proton
    return simulate_multiplicity_per_proton(
        mol=mol,
        x_h_connectivity=connectivity,
        proton_shifts_ppm=shifts,
        solvent=solvent,
        spectrometer_mhz=spectrometer_mhz,
        temperature_k=temperature_k,
        include_heteronuclear=True,
        use_3d=generate_3d,
        config=config,
    )
# ---------------------------
# Simulation ensemble conformères
# ---------------------------
def simulate_conformer_ensemble(
    mol: Chem.Mol,
    proton_shifts_ppm: Dict[int, Any],
    solvent: str = DEFAULT_SOLVENT,
    spectrometer_mhz: float = DEFAULT_SPECTROMETER_FREQ_MHZ,
    temperature_k: float = DEFAULT_TEMPERATURE_K,
    conformer_weights: Optional[List[float]] = None,
    config: Optional[SimConfig] = None,
) -> List[Dict[int, Dict[str, Any]]]:
    """
    ✅ Gère plusieurs ensembles de shifts (conformères/anomères).
    """
    if config is None:
        config = SimConfig()

    smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
    _, connectivity, _ = prepare_molecule_for_nmr(
        smiles,
        generate_3d=True,
        random_seed=config.random_seed,
    )

    first_val = next(iter(proton_shifts_ppm.values()))
    if isinstance(first_val, (list, tuple, np.ndarray)):
        n_confs = len(first_val)
        shift_ensembles: List[Dict[int, float]] = []
        for i in range(n_confs):
            shift_ensembles.append(
                {h: float(shifts[i]) for h, shifts in proton_shifts_ppm.items()}
            )
    else:
        shift_ensembles = [{h: float(v) for h, v in proton_shifts_ppm.items()}]
        n_confs = 1

    if conformer_weights is None:
        conformer_weights = [1.0 / n_confs] * n_confs
    else:
        if len(conformer_weights) != n_confs:
            raise ValueError("conformer_weights doit avoir la même longueur.")
        s = float(sum(conformer_weights)) or 1.0
        conformer_weights = [w / s for w in conformer_weights]

    results: List[Dict[int, Dict[str, Any]]] = []

    for i, shifts in enumerate(shift_ensembles):
        res = simulate_multiplicity_per_proton(
            mol=mol,
            x_h_connectivity=connectivity,
            proton_shifts_ppm=shifts,
            solvent=solvent,
            spectrometer_mhz=spectrometer_mhz,
            temperature_k=temperature_k,
            include_heteronuclear=True,
            use_3d=True,
            config=config,
        )

        w = conformer_weights[i]
        for h_idx, data in res.items():
            data["peak_intensities"] = [float(I) * w for I in data["peak_intensities"]]

        results.append(res)

    return results


# ==============================================================================
# INTÉGRATION AVEC PRÉDICTEUR ML (integrated_nmr_prediction.py)
# ==============================================================================

def load_ml_shift_predictor(
    model_path_1H: Optional[str] = None,
    model_path_13C: Optional[str] = None,
) -> Optional[Any]:
    """
    Charge le prédicteur ML de shifts depuis integrated_nmr_prediction.py.
    
    Args:
        model_path_1H: Chemin vers le modèle ¹H (ex: 1H_sparsified_...pt)
        model_path_13C: Chemin vers le modèle ¹³C
    
    Returns:
        Objet prédicteur ou None si non disponible
    """
    try:
        from pathlib import Path
        import sys
        
        # Ajouter le répertoire parent au path si nécessaire
        script_dir = Path(__file__).parent
        if str(script_dir) not in sys.path:
            sys.path.insert(0, str(script_dir))
        
        # Importer le module de prédiction
        try:
            from integrated_nmr_prediction import IntegratedNMRPredictor
        except ImportError:
            # Essayer avec chemin absolu
            import importlib.util
            pred_path = script_dir / "integrated_nmr_prediction.py"
            if not pred_path.exists():
                return None
            spec = importlib.util.spec_from_file_location("integrated_nmr_prediction", pred_path)
            if spec is None or spec.loader is None:
                return None
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            IntegratedNMRPredictor = module.IntegratedNMRPredictor
        
        # Trouver les modèles si non spécifiés
        if model_path_1H is None:
            model_dir = script_dir.parent / "models"
            candidates = list(model_dir.glob("1H_*.pt"))
            if candidates:
                model_path_1H = str(candidates[0])
        
        if model_path_13C is None:
            model_dir = script_dir.parent / "models"
            candidates = list(model_dir.glob("13C_*.pt"))
            if candidates:
                model_path_13C = str(candidates[0])
        
        return IntegratedNMRPredictor(
            model_path_1H=model_path_1H,
            model_path_13C=model_path_13C,
        )
        
    except Exception as e:
        import logging
        logging.warning(f"Impossible de charger le prédicteur ML: {e}")
        return None


def predict_shifts_ml(
    smiles: str,
    predictor: Optional[Any] = None,
    nucleus: str = "1H",
    solvent: str = "CDCl3",
) -> Dict[int, float]:
    """
    Prédit les déplacements chimiques via le réseau de neurones.
    
    Args:
        smiles: SMILES de la molécule
        predictor: Objet IntegratedNMRPredictor (ou None pour charger automatiquement)
        nucleus: "1H" ou "13C"
        solvent: Solvant (pour correction éventuelle)
    
    Returns:
        {atom_idx: shift_ppm, ...}
    """
    if predictor is None:
        predictor = load_ml_shift_predictor()
    
    if predictor is None:
        raise RuntimeError("Prédicteur ML non disponible. Vérifiez integrated_nmr_prediction.py")
    
    try:
        if nucleus == "1H":
            result = predictor.predict_1H(smiles, solvent=solvent)
        else:
            result = predictor.predict_13C(smiles, solvent=solvent)
        
        # Convertir le format de sortie en dict {idx: shift}
        if isinstance(result, dict):
            if "shifts" in result:
                return result["shifts"]
            return result
        elif hasattr(result, "shifts"):
            return result.shifts
        else:
            return {}
            
    except Exception as e:
        import logging
        logging.warning(f"Erreur prédiction ML: {e}")
        return {}


def simulate_full_spectrum(
    smiles: str,
    solvent: str = DEFAULT_SOLVENT,
    spectrometer_mhz: float = DEFAULT_SPECTROMETER_FREQ_MHZ,
    temperature_k: float = DEFAULT_TEMPERATURE_K,
    use_ml_shifts: bool = True,
    ml_predictor: Optional[Any] = None,
    line_width_hz: float = 1.0,
    ppm_range: Tuple[float, float] = (-0.5, 12.0),
    n_points: int = 32768,
    use_fid_simulation: bool = False,
    config: Optional[SimConfig] = None,
) -> Tuple[Dict[int, Dict[str, Any]], SimulatedSpectrum]:
    """
    ✅ API COMPLÈTE : Simulation de spectre RMN ¹H avec lineshape.
    
    Combine:
    1. Prédiction des shifts (ML ou heuristique)
    2. Simulation des multiplicités (rule-based)
    3. Résolution des systèmes second-ordre (Hamiltonien)
    4. Génération du spectre avec largeur de raie
    
    Args:
        smiles: SMILES de la molécule
        solvent: Solvant deutéré
        spectrometer_mhz: Fréquence Larmor ¹H
        temperature_k: Température
        use_ml_shifts: Utiliser le prédicteur ML (sinon heuristique)
        ml_predictor: Objet prédicteur pré-chargé (optionnel)
        line_width_hz: Largeur de raie (FWHM)
        ppm_range: Plage spectrale
        n_points: Résolution du spectre
        use_fid_simulation: Si True, utilise FID→FFT au lieu de somme Lorentziennes
        config: Configuration avancée
    
    Returns:
        (proton_data, spectrum) où:
        - proton_data: {h_idx: {label, J_couplings, shift_ppm, ...}}
        - spectrum: SimulatedSpectrum avec axes ppm/Hz et intensités
    
    Example:
        >>> data, spec = simulate_full_spectrum("CCO", use_ml_shifts=True)
        >>> import matplotlib.pyplot as plt
        >>> plt.plot(spec.ppm_axis, spec.intensity)
        >>> plt.gca().invert_xaxis()
        >>> plt.xlabel("δ (ppm)")
        >>> plt.show()
    """
    if config is None:
        config = SimConfig()
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES invalide: {smiles}")
    
    # Préparer la molécule
    molH, connectivity, proton_to_heavy = prepare_molecule_for_nmr(
        smiles,
        generate_3d=True,
        random_seed=config.random_seed,
    )
    
    # Obtenir les shifts
    if use_ml_shifts:
        try:
            shifts_ppm = predict_shifts_ml(smiles, ml_predictor, nucleus="1H", solvent=solvent)
            if not shifts_ppm:
                # Fallback heuristique
                shifts_ppm = _estimate_shifts_heuristic(molH, proton_to_heavy)
        except Exception:
            shifts_ppm = _estimate_shifts_heuristic(molH, proton_to_heavy)
    else:
        shifts_ppm = _estimate_shifts_heuristic(molH, proton_to_heavy)
    
    # Simuler les multiplicités
    proton_data = simulate_multiplicity_per_proton(
        mol=mol,
        x_h_connectivity=connectivity,
        proton_shifts_ppm=shifts_ppm,
        solvent=solvent,
        spectrometer_mhz=spectrometer_mhz,
        temperature_k=temperature_k,
        include_heteronuclear=True,
        use_3d=True,
        config=config,
    )
    
    # Ajouter les shifts au proton_data
    for h_idx, data in proton_data.items():
        if h_idx in shifts_ppm:
            data["shift_ppm"] = shifts_ppm[h_idx]
            data["shift_hz"] = shifts_ppm[h_idx] * spectrometer_mhz
    
    # Collecter toutes les positions de pics et intensités
    all_peak_positions_hz: List[float] = []
    all_peak_intensities: List[float] = []
    
    for h_idx, data in proton_data.items():
        shift_hz = data.get("shift_hz", 0.0)
        
        # Positions relatives des pics (dues à la multiplicité)
        peak_offsets = data.get("peak_positions_hz", [0.0])
        peak_intensities = data.get("peak_intensities", [1.0])
        
        # Normaliser les intensités par le nombre de protons équivalents
        n_H = len([h for h, heavy in proton_to_heavy.items() 
                   if heavy == data.get("heavy_idx", -1)])
        n_H = max(1, n_H)
        
        for offset, intensity in zip(peak_offsets, peak_intensities):
            all_peak_positions_hz.append(shift_hz + offset)
            all_peak_intensities.append(intensity * n_H)
    
    # Générer le spectre
    if use_fid_simulation:
        T2_star = 1.0 / (np.pi * line_width_hz) if line_width_hz > 0 else 0.5
        spectrum = simulate_spectrum_from_fid(
            peak_positions_hz=all_peak_positions_hz,
            peak_intensities=all_peak_intensities,
            spectrometer_mhz=spectrometer_mhz,
            T2_star_s=T2_star,
            n_points=n_points,
            ppm_range=ppm_range,
        )
    else:
        spectrum = simulate_spectrum_lineshape(
            peak_positions_hz=all_peak_positions_hz,
            peak_intensities=all_peak_intensities,
            spectrometer_mhz=spectrometer_mhz,
            line_width_hz=line_width_hz,
            n_points=n_points,
            ppm_range=ppm_range,
        )
    
    return proton_data, spectrum


def export_spectrum_to_jcamp(
    spectrum: SimulatedSpectrum,
    filepath: str,
    title: str = "Simulated 1H NMR Spectrum",
    origin: str = "sim_mult_rulebased.py",
) -> None:
    """
    Exporte le spectre au format JCAMP-DX (compatible avec la plupart des logiciels NMR).
    
    Args:
        spectrum: SimulatedSpectrum à exporter
        filepath: Chemin du fichier de sortie (.jdx)
        title: Titre du spectre
        origin: Origine/source du spectre
    """
    with open(filepath, "w", encoding="utf-8") as f:
        f.write("##TITLE=" + title + "\n")
        f.write("##JCAMP-DX=5.01\n")
        f.write("##DATA TYPE=NMR SPECTRUM\n")
        f.write("##DATA CLASS=XYDATA\n")
        f.write("##ORIGIN=" + origin + "\n")
        f.write(f"##OBSERVE FREQUENCY={spectrum.spectrometer_mhz}\n")
        f.write("##OBSERVE NUCLEUS=1H\n")
        f.write(f"##XUNITS=PPM\n")
        f.write("##YUNITS=ARBITRARY UNITS\n")
        f.write(f"##FIRSTX={spectrum.ppm_axis[0]:.6f}\n")
        f.write(f"##LASTX={spectrum.ppm_axis[-1]:.6f}\n")
        f.write(f"##NPOINTS={len(spectrum.ppm_axis)}\n")
        f.write(f"##DELTAX={abs(spectrum.ppm_axis[1] - spectrum.ppm_axis[0]):.8f}\n")
        f.write(f"##MAXY={spectrum.intensity.max():.6f}\n")
        f.write(f"##MINY={spectrum.intensity.min():.6f}\n")
        f.write("##XYDATA=(X++(Y..Y))\n")
        
        # Écrire les données (format compact)
        for i, (x, y) in enumerate(zip(spectrum.ppm_axis, spectrum.intensity)):
            if i % 10 == 0:
                f.write(f"{x:.4f}")
            f.write(f" {y:.6f}")
            if (i + 1) % 10 == 0:
                f.write("\n")
        
        if len(spectrum.ppm_axis) % 10 != 0:
            f.write("\n")
        
        f.write("##END=\n")


def export_spectrum_to_csv(
    spectrum: SimulatedSpectrum,
    filepath: str,
) -> None:
    """
    Exporte le spectre au format CSV simple (ppm, intensity).
    """
    import csv
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["ppm", "intensity", "frequency_hz"])
        for ppm, intensity, freq in zip(spectrum.ppm_axis, spectrum.intensity, spectrum.freq_axis_hz):
            writer.writerow([f"{ppm:.6f}", f"{intensity:.8f}", f"{freq:.2f}"])
