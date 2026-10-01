# -*- coding: utf-8 -*-

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

DEBUG_LOGGING = False

try:
    from heterocycles_complete import (
        identify_heterocycle,
        get_heterocycle_coupling,
        is_heterocycle_singlet,
        estimate_heterocycle_couplings_complete,
        identify_saturated_heterocycle,
        get_saturated_heterocycle_restrictions,
        estimate_saturated_heterocycle_couplings,
        should_differentiate_diastereotopic_J,
        HeterocycleMatch
    )
    HETEROCYCLES_MODULE_AVAILABLE = True
except (ImportError, ModuleNotFoundError) as e:

    try:
        from .heterocycles_complete import (
            identify_heterocycle,
            get_heterocycle_coupling,
            is_heterocycle_singlet,
            estimate_heterocycle_couplings_complete,
            identify_saturated_heterocycle,
            get_saturated_heterocycle_restrictions,
            estimate_saturated_heterocycle_couplings,
            should_differentiate_diastereotopic_J,
            HeterocycleMatch
        )
        HETEROCYCLES_MODULE_AVAILABLE = True
    except (ImportError, ModuleNotFoundError):
        HETEROCYCLES_MODULE_AVAILABLE = False
except Exception as e:
    print(f"[ERROR] Unexpected error importing heterocycles_complete: {e}")
    import traceback
    traceback.print_exc()
    HETEROCYCLES_MODULE_AVAILABLE = False

try:
    from sugars_complete import (
        estimate_sugar_couplings_complete,
        is_sugar_carbon,
        find_sugar_rings
    )
    SUGARS_MODULE_AVAILABLE = True
except ImportError:
    SUGARS_MODULE_AVAILABLE = False

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

        get_glycosidic_long_range_coupling,
        get_nucleoside_chi_coupling,
    )
    COUPLING_TABLES_MODULE_AVAILABLE = True
except ImportError:
    COUPLING_TABLES_MODULE_AVAILABLE = False

try:
    from peptide_couplings import (
        get_peptide_coupling_robust,
        calculate_3J_HN_HA,
        calculate_3J_HA_HB,
        calculate_2J_geminal as calculate_peptide_2J_geminal,
        KARPLUS_PARAMETERS,
        AMINO_ACID_SPECIFIC_COUPLINGS,

        rotamer_averaged_J,
        get_diastereotopic_J_pair,
        ROTAMER_POPULATIONS,
        AminoAcid,
    )
    PEPTIDE_MODULE_AVAILABLE = True
    ROTAMER_AVERAGING_AVAILABLE = True
except ImportError:
    PEPTIDE_MODULE_AVAILABLE = False
    ROTAMER_AVERAGING_AVAILABLE = False

try:
    from coupling_extensions import (
        get_extended_coupling,
        get_polycyclic_aromatic_coupling,
        get_strained_ring_coupling,
        detect_strained_ring,

        detect_exotic_carbocycle,
        get_exotic_carbocycle_coupling,

        get_aromatic_F_coupling,
        AROMATIC_FLUORINE_COUPLINGS,
        POLYCYCLIC_AROMATIC_COUPLINGS,
        STRAINED_RING_COUPLINGS,
        GLYCOSIDIC_COUPLINGS,
        TERPENE_SKELETON_COUPLINGS,
        STEROID_RING_COUPLINGS,
        PHOSPHOLIPID_COUPLINGS,
        NUCLEOTIDE_COUPLINGS,

        ANNULENE_COUPLINGS,
        TROPYLIUM_COUPLINGS,
        ANTIAROMATIC_COUPLINGS,
        CURVED_AROMATIC_COUPLINGS,
        EXOTIC_CAGE_COUPLINGS,

        estimate_steroid_couplings_differentiated,
        estimate_terpene_couplings_differentiated,
        estimate_macrocycle_couplings_differentiated,
        STEROID_COUPLINGS_DIFFERENTIATED,
        TERPENE_COUPLINGS_DIFFERENTIATED,
        MACROCYCLE_COUPLINGS_DIFFERENTIATED,
    )
    COUPLING_EXTENSIONS_AVAILABLE = True
except ImportError:
    COUPLING_EXTENSIONS_AVAILABLE = False

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

try:
    from nmr_plugins.dynamics import apply_dynamics_to_coupling, clear_dynamics_cache, should_merge_diastereotopic_CH2
    DYNAMICS_MODULE_AVAILABLE = True
except ImportError:
    try:
        from .nmr_plugins.dynamics import apply_dynamics_to_coupling, clear_dynamics_cache, should_merge_diastereotopic_CH2
        DYNAMICS_MODULE_AVAILABLE = True
    except ImportError:
        DYNAMICS_MODULE_AVAILABLE = False

        def should_merge_diastereotopic_CH2(mol, heavy_idx, h_a, h_b, temperature_k=298.15, spectrometer_mhz=700.0):

            return (False, "Module dynamics non disponible", {})

try:
    from diastereotopic_delta_table import (
        get_diastereotopic_delta,
        get_diastereotopic_delta_simple,
        get_solvent_multiplier,
        get_solvent_additive_shift,
        CH2_CHIRAL_DELTA_BY_DISTANCE,
        VINYL_DELTA_BY_SUBSTITUENT,
        CYCLIC_DELTA,
        AMINO_ACID_BETA_DELTA,
        ANISOTROPIC_SYSTEMS,
        DEFAULT_DELTA,

        BENZYLIC_DELTA,
        ALLYLIC_DELTA,
        EPOXIDE_DELTA,
        AZIRIDINE_DELTA,
        THIIRANE_DELTA,
        SULFOXIDE_SULFONE_DELTA,
        PHOSPHORUS_DELTA,
        N_OXIDE_DELTA,
        SILANE_DELTA,
        BORON_DELTA,
        NUCLEOSIDE_DELTA,
        CROWN_ETHER_DELTA,
        CYCLODEXTRIN_DELTA,
        ROTAMER_DELTA,
        METAL_COMPLEX_DELTA,

        CH3_DIASTEREOTOPIC_DELTA,
        get_ch3_diastereotopic_delta,
    )
    DIASTEREOTOPIC_TABLE_AVAILABLE = True
except ImportError:
    DIASTEREOTOPIC_TABLE_AVAILABLE = False

    CH2_CHIRAL_DELTA_BY_DISTANCE = {
        1: (0.12, 0.20, "distance 1"),
        2: (0.10, 0.15, "distance 2"),
        3: (0.08, 0.12, "distance 3"),
        4: (0.05, 0.08, "distance 4"),
    }
    DEFAULT_DELTA = {
        "ch2_near_chiral_1": 0.12,
        "ch2_near_chiral_2": 0.10,
        "ch2_near_chiral_3": 0.08,
        "flexible_chain": 0.08,
        "unknown": 0.10,
    }

class MoleculeCache:

    _instance = None
    _cache_topo_dist: Dict[str, np.ndarray] = {}
    _cache_aromatic_map: Dict[str, List[Dict[int, int]]] = {}
    _cache_ranks: Dict[str, List[int]] = {}
    _cache_conformers: Dict[str, List[Chem.Mol]] = {}
    _max_size = 100

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    @classmethod
    def get_key(cls, mol: Chem.Mol) -> str:

        try:
            return Chem.MolToSmiles(mol, canonical=True)
        except Exception:

            atoms_str = "".join(f"{a.GetSymbol()}{a.GetIdx()}" for a in mol.GetAtoms())
            return hashlib.md5(atoms_str.encode()).hexdigest()

    @classmethod
    def get_topological_distances(cls, mol: Chem.Mol) -> np.ndarray:

        key = cls.get_key(mol)
        if key in cls._cache_topo_dist:
            return cls._cache_topo_dist[key]

        dist = _compute_topological_distances(mol)

        if len(cls._cache_topo_dist) >= cls._max_size:
            cls._cache_topo_dist.pop(next(iter(cls._cache_topo_dist)))

        cls._cache_topo_dist[key] = dist
        return dist

    @classmethod
    def get_aromatic_map(cls, mol: Chem.Mol) -> List[Dict[int, int]]:

        key = cls.get_key(mol)
        if key in cls._cache_aromatic_map:
            return cls._cache_aromatic_map[key]

        arom_map = _compute_aromatic_distance_map(mol)

        if len(cls._cache_aromatic_map) >= cls._max_size:
            cls._cache_aromatic_map.pop(next(iter(cls._cache_aromatic_map)))

        cls._cache_aromatic_map[key] = arom_map
        return arom_map

    @classmethod
    def get_canonical_ranks(cls, mol_noH: Chem.Mol) -> List[int]:

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

        cls._cache_topo_dist.clear()
        cls._cache_aromatic_map.clear()
        cls._cache_ranks.clear()
        cls._cache_conformers.clear()

def _compute_topological_distances(mol: Chem.Mol) -> np.ndarray:

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

        for k in range(n):
            D = np.minimum(D, D[:, k:k+1] + D[k:k+1, :])
        return D

def _compute_aromatic_distance_map(mol: Chem.Mol) -> List[Dict[int, int]]:

    n = mol.GetNumAtoms()
    dist_maps = [{} for _ in range(n)]
    aro_idx = {a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()}
    if not aro_idx:
        return dist_maps

    for start in aro_idx:
        seen = {start}
        q = deque([(start, 0)])
        dist_maps[start][start] = 0
        while q:
            u, d = q.popleft()
            if d >= 3:
                continue
            atom = mol.GetAtomWithIdx(u)
            for bond in atom.GetBonds():
                v = bond.GetOtherAtomIdx(u)
                if v not in seen and v in aro_idx:
                    seen.add(v)
                    dist_maps[start][v] = d + 1
                    q.append((v, d + 1))
    return dist_maps

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

DEFAULT_SPECTROMETER_FREQ_MHZ = 700.0
DEFAULT_TEMPERATURE_K = 298.15

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

TOLERANCE_J_HZ = 0.15

FREQ_LINEWIDTH_TABLE = [

    (0,   150, 1.80),
    (150, 250, 1.50),
    (250, 350, 1.30),
    (350, 450, 1.20),
    (450, 600, 1.00),
    (600, 800, 0.90),
    (800, 1e9, 0.80),
]

def _typical_linewidth_from_field(mhz: float, viscosity: str) -> float:
    base = next(w for lo, hi, w in FREQ_LINEWIDTH_TABLE if lo < mhz <= hi)
    visco_add = {"low": 0.00, "moderate": 0.05, "high": 0.15}.get(viscosity, 0.00)
    return base + visco_add

_13C_NATURAL_ABUNDANCE = 0.011

_13C_1JCH_TABLE = {
    'sp3_alkyl': 125.0,
    'sp3_alpha_O': 140.0,
    'sp3_alpha_N': 135.0,
    'sp3_cyclopropane': 160.0,
    'sp2_alkene': 155.0,
    'sp2_aldehyde': 175.0,
    'sp_alkyne': 250.0,
    'aromatic': 160.0,
    'default': 140.0,
}

_GEMINAL_2J_TABLE = {

    'CH2_alkyl': (-12.5, -14.0),
    'default_sp3': (-12.5, -14.0),

    'CH2_alpha_O': (-10.0, -11.0),
    'CH2_alpha_N': (-12.5, -13.5),
    'CH2_alpha_S': (-10.5, -11.5),
    'CH2_alpha_P': (-8.5, -9.5),
    'CH2_acetal': (-5.0, -7.0),

    'CH2_cyclopropane': (-4.0, -5.0),
    'CH2_epoxide': (-5.0, -6.0),
    'CH2_aziridine': (-5.5, -6.5),
    'CH2_cyclobutane': (-11.5, -12.5),
    'CH2_cyclopentane': (-12.5, -13.5),
    'CH2_cyclohexane': (-12.5, -13.0),
    'CH2_exocyclic': (-10.0, -13.0),

    'CH2_benzyl': (-13.5, -14.5),
    'CH2_indolyl': (-14.0, -15.0),
    'CH2_allyl': (-14.0, -15.0),
    'CH2_benzyl_allyl': (-14.5, -15.5),
    'CH2_alpha_C=O': (-15.5, -17.0),
    'CH2_di_carbonyl': (-17.0, -18.5),
    'CH2_alpha_nitro': (-17.5, -18.5),

    'CH2_Ar_C=O': (-16.0, -17.5),
    'CH2_Ar_C=C': (-15.0, -16.0),
    'CH2_N_C=O': (-16.5, -17.5),
    'CH2_O_C=O': (-13.5, -14.5),

    'CH2_sp2': (1.5, 2.5),
    'CH2_allene': (7.0, 8.5),
    'aldehyde': (2.0, 3.0),
    'default_sp2': (1.5, 2.5),
}

_GEMINAL_2J_INCREMENTS = {

    'alpha_O': +2.5,
    'alpha_N': -0.5,
    'alpha_S': +1.5,
    'alpha_P': +3.5,

    'aromatic': +1.5,
    'vinyl': -1.25,
    'carbonyl': -1.5,
    'nitro': -5.5,

    'beta_carbonyl': -2.5,

    'ring_3': +8.0,
    'ring_4': +0.5,
    'ring_5': -0.5,
    'ring_6': 0.0,
}

_4J_LONG_RANGE_TABLE = {

    'allylic_cis': (1.5, 3.0),
    'allylic_trans': (0.8, 2.0),
    'allylic_general': (0.5, 3.0),

    'homoallylic_cis': (1.0, 2.0),
    'homoallylic_trans': (0.5, 1.5),
    'homoallylic_general': (0.5, 2.5),

    'W_allylic_rigid': (1.5, 3.5),
    'W_allylic_flexible': (0.8, 2.0),

    'meta_aromatic': (2.0, 3.5),
    'para_aromatic': (0.3, 0.8),

    'benzylic': (0.5, 1.5),

    'aldehyde_vinylic': (1.0, 2.5),

    'W_saturated_rigid': (1.0, 2.5),
    'W_saturated_flexible': (0.0, 0.5),

    'cyclopropyl': (1.0, 2.0),
    'cyclobutyl': (0.5, 1.5),

    'default_rigid': (1.0, 2.0),
    'default_flexible': (0.0, 0.8),
}

_5J_LONG_RANGE_TABLE = {

    'extended_W_rigid': (0.8, 2.0),
    'extended_W_flexible': (0.0, 0.3),

    'para_substituted': (0.4, 0.8),

    'homoallylic_extended': (0.3, 1.0),

    'conjugated_diene': (0.5, 1.5),
    'conjugated_enyne': (0.3, 0.8),

    'peri_coupling': (0.5, 1.0),

    'default': (0.0, 0.5),
}

_HF_COUPLING_TABLE = {

    '2J_CHF_sp3': (42.0, 52.0),
    '2J_CHF_sp2': (44.0, 55.0),
    '2J_CHF_formyl': (48.0, 56.0),
    '2J_CHF_CF2': (50.0, 58.0),
    '2J_CHF_CF3': (44.0, 50.0),

    '3J_HCCF_anti': (20.0, 35.0),
    '3J_HCCF_gauche': (5.0, 15.0),
    '3J_HCCF_average': (8.0, 22.0),
    '3J_HCCF_cyclic': (10.0, 25.0),
    '3J_HCCF_aromatic': (6.0, 10.0),

    '4J_HCCCF_W': (3.0, 8.0),
    '4J_HCCCF_general': (0.5, 4.0),
    '4J_HCCCF_aromatic': (2.0, 5.0),
    '4J_allylic_F': (2.0, 6.0),

    '5J_HCCCCF_W': (1.0, 3.0),
    '5J_HCCCCF_aromatic': (0.5, 2.0),
    '5J_HCCCCF_general': (0.0, 1.5),

    'vinylic_2J': (75.0, 90.0),
    'vinylic_3J_cis': (18.0, 28.0),
    'vinylic_3J_trans': (35.0, 45.0),

    'default_2J': (46.0, 50.0),
    'default_3J': (10.0, 20.0),
    'default_4J': (1.0, 4.0),
    'default_5J': (0.0, 1.0),
}

_HP_COUPLING_TABLE = {

    '2J_HCP_sp3': (8.0, 18.0),
    '2J_HCP_sp2': (10.0, 22.0),
    '2J_HCP_phosphonate': (18.0, 25.0),

    '3J_HCCP_anti': (12.0, 20.0),
    '3J_HCCP_gauche': (2.0, 8.0),
    '3J_HCCP_average': (5.0, 12.0),
    '3J_HCCP_sugar': (8.0, 12.0),

    '3J_HCOP_ester': (6.0, 12.0),
    '3J_HCOP_nucleotide': (4.0, 10.0),

    '4J_HCCCP_W': (2.0, 5.0),
    '4J_HCCCP_general': (0.5, 2.0),

    'default_2J': (12.0, 18.0),
    'default_3J': (6.0, 12.0),
    'default_4J': (0.5, 2.0),
}

_GARBISCH_ALKENE_EXTENDED = {

    'terminal': (11.0, -3.0, 2.0),
    'disubstituted_EZ': (11.5, -3.5, 2.0),
    'trisubstituted': (10.5, -2.5, 1.5),
    'tetrasubstituted': (0.0, 0.0, 0.0),

    'diene_terminal': (11.5, -4.0, 2.5),
    'diene_internal': (11.0, -3.5, 2.0),
    'enone': (10.5, -3.0, 1.5),
    'styrene': (11.0, -3.0, 2.0),
    'stilbene': (11.5, -4.0, 2.5),

    'cyclopropene': (8.0, -1.5, 1.0),
    'cyclobutene': (9.0, -2.0, 1.5),
    'cyclopentene': (10.0, -2.5, 1.5),
    'cyclohexene': (10.5, -3.0, 2.0),
    'cycloheptene': (11.0, -3.0, 2.0),
    'cyclooctene': (10.5, -3.5, 2.5),

    'exocyclic_3': (8.5, -2.0, 1.5),
    'exocyclic_5': (10.0, -2.5, 2.0),
    'exocyclic_6': (10.5, -3.0, 2.0),

    'vinyl_ether': (9.5, -2.5, 1.5),
    'vinyl_ester': (9.0, -2.0, 1.5),
    'vinyl_halide_F': (8.0, -2.0, 1.0),
    'vinyl_halide_Cl': (9.5, -2.5, 1.5),
    'vinyl_halide_Br': (9.5, -2.5, 1.5),
    'vinyl_amine': (10.0, -3.0, 2.0),

    'EN_correction_per_O': -0.8,
    'EN_correction_per_N': -0.5,
    'EN_correction_per_F': -1.2,
    'EN_correction_per_Cl': -0.6,

    'default': (11.0, -3.0, 2.0),
}

_KARPLUS_EXTENDED_PARAMS = {

    'standard_HCCH': {
        'A': 7.76, 'B': -1.10, 'C': 1.40, 'D': 0.0, 'E': 0.0,
        'reference': 'Karplus M., JACS 85, 2870 (1963)'
    },

    'haasnoot_altona': {
        'P1': 13.70, 'P2': -0.73, 'P3': 0.0, 'P4': 0.56, 'P5': -2.47, 'P6': 16.90, 'P7': 0.14,
        'lambda_i': [0.0, 0.0, 0.0, 0.0],
        'reference': 'Haasnoot-de Leeuw-Altona, Tetrahedron 36, 2783 (1980)'
    },

    'vinylic': {
        'A': 11.0, 'B': -3.0, 'C': 2.0, 'D': 0.0, 'E': 0.0,
        'cis': 10.5, 'trans': 17.0, 'geminal': 2.0,
        'reference': 'Garbisch E.W., JACS 86, 5561 (1964)'
    },

    'furanose_altona': {
        'H1H2': {'A': 9.67, 'B': -1.95, 'C': 0.97, 'cos_coeff': 0.462},
        'H2H3': {'A': 9.67, 'B': -1.95, 'C': 0.97, 'cos_coeff': 0.0},
        'H3H4': {'A': 9.67, 'B': -1.95, 'C': 0.97, 'cos_coeff': -0.462},
        'reference': 'Altona-Sundaralingam, JACS 94, 8205 (1972)'
    },

    'peptide_phi': {
        'A': 6.98, 'B': -1.38, 'C': 1.72, 'D': -0.11, 'E': 0.0,
        'phase': -60.0,
        'reference': 'Vuister-Bax, JACS 115, 7772 (1993)'
    },

    'peptide_chi1': {
        'A': 9.50, 'B': -1.68, 'C': 1.80, 'D': 0.0, 'E': 0.0,
        'reference': 'DeMarco-Llinas-Wuthrich, Biopolymers 17, 617 (1978)'
    },

    'nucleoside_pseudorot': {
        'A': 9.3, 'B': -0.7, 'C': 0.3, 'D': 0.0, 'E': 0.0,
        'P_north': 18.0,
        'P_south': 162.0,
        'reference': 'Altona-Sundaralingam, JACS 95, 2333 (1973)'
    },

    'HCCF': {
        'A': 25.0, 'B': -3.0, 'C': 5.0, 'D': 0.0, 'E': 0.0,
        'reference': 'Abraham-Pachler, Mol Phys 7, 165 (1963)'
    },

    'HCCP': {
        'A': 15.0, 'B': -2.0, 'C': 3.0, 'D': 0.0, 'E': 0.0,
        'reference': 'Gorenstein, Phosphorus-31 NMR (1984)'
    },
}

def _get_4J_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    coupling_type: str = None,
    mol3d: Chem.Mol = None
) -> Optional[float]:

    if coupling_type:

        if coupling_type in _4J_LONG_RANGE_TABLE:
            j_range = _4J_LONG_RANGE_TABLE[coupling_type]
            return (j_range[0] + j_range[1]) / 2.0

    try:
        atom1 = mol.GetAtomWithIdx(h1_idx)
        atom2 = mol.GetAtomWithIdx(h2_idx)
    except Exception:
        return None

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

    if c1.GetIsAromatic() and c2.GetIsAromatic():

        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if c1_idx in ring and c2_idx in ring:

                idx1 = list(ring).index(c1_idx)
                idx2 = list(ring).index(c2_idx)
                dist = min(abs(idx1 - idx2), len(ring) - abs(idx1 - idx2))
                if dist == 2:
                    return 2.5

    def find_path_4_bonds(start_idx, end_idx):

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

        has_double_bond = False
        for i in range(len(path) - 1):
            bond = mol.GetBondBetweenAtoms(path[i], path[i+1])
            if bond and bond.GetBondTypeAsDouble() >= 1.5:
                has_double_bond = True
                break

        if has_double_bond:

            is_rigid = all(mol.GetAtomWithIdx(idx).IsInRing() for idx in path)
            if is_rigid:
                return 2.0
            else:
                return 1.5

    if c1.IsInRing() and c2.IsInRing():

        if mol3d and mol3d.GetNumConformers() > 0:
            try:
                conf = mol3d.GetConformer()

                return 1.2
            except Exception:
                pass

    if (c1.GetIsAromatic() != c2.GetIsAromatic()):

        return 1.0

    return None

def _get_5J_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    coupling_type: str = None,
    mol3d: Chem.Mol = None
) -> Optional[float]:

    if coupling_type:
        if coupling_type in _5J_LONG_RANGE_TABLE:
            j_range = _5J_LONG_RANGE_TABLE[coupling_type]
            return (j_range[0] + j_range[1]) / 2.0

    try:
        atom1 = mol.GetAtomWithIdx(h1_idx)
        atom2 = mol.GetAtomWithIdx(h2_idx)
    except Exception:
        return None

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

    if c1.GetIsAromatic() and c2.GetIsAromatic():
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if len(ring) == 6 and c1_idx in ring and c2_idx in ring:
                idx1 = list(ring).index(c1_idx)
                idx2 = list(ring).index(c2_idx)
                dist = min(abs(idx1 - idx2), len(ring) - abs(idx1 - idx2))
                if dist == 3:
                    return 0.6

    if c1.IsInRing() and c2.IsInRing():

        ring_info = mol.GetRingInfo()
        rings_c1 = [set(r) for r in ring_info.AtomRings() if c1_idx in r]
        rings_c2 = [set(r) for r in ring_info.AtomRings() if c2_idx in r]

        for r1 in rings_c1:
            for r2 in rings_c2:
                if len(r1 & r2) >= 2:
                    return 1.0

    return None

def _get_HF_coupling(
    mol: Chem.Mol,
    h_idx: int,
    f_idx: int,
    mol3d: Chem.Mol = None
) -> Optional[float]:

    try:

        path = Chem.GetShortestPath(mol, h_idx, f_idx)
        n_bonds = len(path) - 1 if path else 0
    except Exception:
        return None

    if n_bonds == 0:
        return None

    if n_bonds == 2:
        j_range = _HF_COUPLING_TABLE.get('default_2J', (46.0, 50.0))
        return (j_range[0] + j_range[1]) / 2.0

    if n_bonds == 3:
        if mol3d and mol3d.GetNumConformers() > 0:
            try:
                conf = mol3d.GetConformer()
                dihedral = abs(AllChem.GetDihedralDeg(conf, path[0], path[1], path[2], path[3]))

                if dihedral > 150 or dihedral < 30:
                    return 25.0
                elif dihedral > 60 and dihedral < 120:
                    return 8.0
                else:
                    return 15.0
            except Exception:
                pass

        j_range = _HF_COUPLING_TABLE.get('3J_HCCF_average', (8.0, 22.0))
        return (j_range[0] + j_range[1]) / 2.0

    if n_bonds == 4:

        is_aromatic = any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in path[1:-1])
        if is_aromatic:
            return 3.5
        else:
            return 2.0

    if n_bonds == 5:
        return 1.0

    return None

def _get_HP_coupling(
    mol: Chem.Mol,
    h_idx: int,
    p_idx: int,
    mol3d: Chem.Mol = None
) -> Optional[float]:

    try:
        path = Chem.GetShortestPath(mol, h_idx, p_idx)
        n_bonds = len(path) - 1 if path else 0
    except Exception:
        return None

    if n_bonds == 0:
        return None

    if n_bonds == 2:
        j_range = _HP_COUPLING_TABLE.get('default_2J', (12.0, 18.0))
        return (j_range[0] + j_range[1]) / 2.0

    if n_bonds == 3:

        via_oxygen = any(mol.GetAtomWithIdx(idx).GetSymbol() == 'O' for idx in path[1:-1])
        if via_oxygen:
            return 8.0
        else:
            return 9.0

    if n_bonds == 4:
        return 1.5

    return None

def _get_vinylic_J_from_stereo(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    default_cis: float = 10.0,
    default_trans: float = 16.0,
    default_geminal: float = 2.0
) -> Optional[Tuple[float, str]]:

    try:
        h1 = mol.GetAtomWithIdx(h1_idx)
        h2 = mol.GetAtomWithIdx(h2_idx)

        if h1.GetSymbol() != 'H' or h2.GetSymbol() != 'H':
            return None

        c1_list = [n for n in h1.GetNeighbors() if n.GetSymbol() == 'C']
        c2_list = [n for n in h2.GetNeighbors() if n.GetSymbol() == 'C']

        if not c1_list or not c2_list:
            return None

        c1 = c1_list[0]
        c2 = c2_list[0]
        c1_idx = c1.GetIdx()
        c2_idx = c2.GetIdx()

        if c1_idx == c2_idx:
            return (default_geminal, 'geminal')

        bond = mol.GetBondBetweenAtoms(c1_idx, c2_idx)
        if not bond:
            return None
        if bond.GetBondTypeAsDouble() < 1.9:
            return None

        stereo = bond.GetStereo()

        if stereo == Chem.BondStereo.STEREONONE:

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
            return (default_trans, 'unknown')

        stereo_atoms = bond.GetStereoAtoms()
        if len(stereo_atoms) < 2:
            return (default_trans, 'unknown')

        ref1_idx, ref2_idx = stereo_atoms[0], stereo_atoms[1]

        h1_is_ref1 = (h1_idx == ref1_idx)
        h2_is_ref2 = (h2_idx == ref2_idx)

        if stereo in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOTRANS):

            if h1_is_ref1 == h2_is_ref2:
                return (default_trans, 'trans')
            else:
                return (default_cis, 'cis')

        elif stereo in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOCIS):

            if h1_is_ref1 == h2_is_ref2:
                return (default_cis, 'cis')
            else:
                return (default_trans, 'trans')

        return (default_trans, 'unknown')

    except Exception:
        return None

def _is_symmetric_para_disubstituted_benzene(
    mol: Chem.Mol,
    heavy_idx: int,
    aromatic_dmap: List[Dict[int, int]],
    ranks: Optional[List[int]] = None,
) -> bool:

    atom = mol.GetAtomWithIdx(heavy_idx)
    if not atom.GetIsAromatic() or atom.GetSymbol() != 'C':
        return False

    ring_info = mol.GetRingInfo()
    six_rings = [r for r in ring_info.AtomRings() if len(r) == 6 and heavy_idx in r]

    if not six_rings:
        return False

    ring = six_rings[0]
    ring_set = set(ring)

    for idx in ring:
        a = mol.GetAtomWithIdx(idx)
        if a.GetSymbol() != 'C' or not a.GetIsAromatic():
            return False

    substituted_positions = []
    h_positions = []

    for idx in ring:
        a = mol.GetAtomWithIdx(idx)

        n_H = a.GetTotalNumHs()

        n_H_explicit = sum(1 for n in a.GetNeighbors() if n.GetSymbol() == 'H' or n.GetAtomicNum() == 1)
        n_H = max(n_H, n_H_explicit)

        if n_H == 0:

            non_ring_neighbors = [n for n in a.GetNeighbors() if n.GetIdx() not in ring_set]
            if non_ring_neighbors:
                substituted_positions.append(idx)
        else:
            h_positions.append(idx)

    if len(substituted_positions) != 2:
        return False

    p1, p2 = substituted_positions

    try:
        idx1 = ring.index(p1)
        idx2 = ring.index(p2)
        dist = min(abs(idx1 - idx2), 6 - abs(idx1 - idx2))
        if dist != 3:
            return False
    except:
        return False

    a1 = mol.GetAtomWithIdx(p1)
    a2 = mol.GetAtomWithIdx(p2)

    sub1 = [n for n in a1.GetNeighbors() if n.GetIdx() not in ring_set]
    sub2 = [n for n in a2.GetNeighbors() if n.GetIdx() not in ring_set]

    if len(sub1) != 1 or len(sub2) != 1:
        return False

    s1, s2 = sub1[0], sub2[0]

    if s1.GetSymbol() != s2.GetSymbol():
        return False

    if s1.GetTotalNumHs() != s2.GetTotalNumHs():
        return False

    if s1.GetDegree() != s2.GetDegree():
        return False

    if ranks is not None:

        pass

    return True

def _is_ortho_disubstituted_symmetric(
    mol: Chem.Mol,
    heavy_idx: int,
    aromatic_dmap: List[Dict[int, int]],
    ranks: Optional[List[int]] = None,
) -> bool:

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
        if dist != 1:
            return False
    except:
        return False

    if ranks is not None and ranks[p1] != ranks[p2]:
        return False

    return True

_AROMATIC_COUPLING_TABLE = {

    'benzene': {'ortho': 7.5, 'meta': 1.5, 'para': 0.5},
    'phenyl_electron_donor': {'ortho': 7.8, 'meta': 1.8, 'para': 0.4},
    'phenyl_electron_acceptor': {'ortho': 8.2, 'meta': 2.2, 'para': 0.6},

    'furan': {
        'ortho': 1.8, 'meta': 0.9, 'para': 1.5,
        '2_3': 1.8, '3_4': 3.4, '2_4': 0.9, '2_5': 1.5, '3_5': 0.9,

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

        'ortho': 1.5, 'meta': 0.0, 'para': 0.0,
        '4_5': 1.5, '2_4': 0.0, '2_5': 0.0,
        'special_h2_singlet': True,
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

    'pyridine': {

        'ortho': 6.3,
        'meta': 1.8,
        'para': 0.5,
        '2_3': 4.8, '3_4': 7.8, '4_5': 7.8, '5_6': 4.8,
        '2_4': 1.8, '2_5': 0.9, '2_6': 0.4,
        '3_5': 1.5, '3_6': 0.5, '4_6': 1.8,
    },
    'pyrimidine': {
        'ortho': 5.0, 'meta': 1.5,
        '4_5': 5.0, '5_6': 5.0, '4_6': 1.5,
    },
    'pyrazine': {
        'ortho': 2.0, 'meta': 0.5,
        '2_3': 2.0, '2_5': 0.5,
    },
    'pyridazine': {
        'ortho': 5.0, 'meta': 2.0,
        '3_4': 5.0, '4_5': 5.0, '5_6': 5.0, '3_5': 2.0, '4_6': 2.0,
    },
    'triazine_1_3_5': {

        'ortho': 0.0, 'meta': 0.0,
    },

    'quinoline': {

        'ortho': 8.0, 'meta': 1.5, 'para': 0.5,
        '2_3': 4.5, '3_4': 8.0, '5_6': 8.0, '6_7': 7.5, '7_8': 8.0,
        '2_4': 1.8, '5_7': 1.3, '5_8': 0.8, '6_8': 1.3,
    },
    'isoquinoline': {
        'ortho': 8.0, 'meta': 1.5, 'para': 0.5,
        '1_3': 0.9, '1_4': 0.5, '3_4': 6.5, '5_6': 8.0, '6_7': 7.5, '7_8': 8.0,
    },
    'indole': {

        'ortho': 7.8, 'meta': 1.2, 'para': 0.5,
        '2_3': 2.2, '4_5': 8.0, '5_6': 7.5, '6_7': 7.5, '7_4': 0.8,
        '3_4': 0.0,
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
        'special_h2_singlet': True,
    },

    'naphthalene': {'1_2': 8.0, '2_3': 7.0, '3_4': 8.0, '1_3': 1.3, '1_4': 0.8, '2_4': 1.3,
                    '5_6': 8.0, '6_7': 7.0, '7_8': 8.0, 'ortho': 8.0, 'meta': 1.3, 'para': 0.8},

    'default_ortho': 7.5,
    'default_meta': 2.0,
    'default_para': 0.5,
}

def _is_imidazole_H2_carbon(mol: Chem.Mol, heavy_idx: int) -> bool:

    try:
        atom = mol.GetAtomWithIdx(heavy_idx)
        if atom.GetSymbol() != 'C' or not atom.GetIsAromatic():
            return False

        ring_info = mol.GetRingInfo()
        five_rings = [r for r in ring_info.AtomRings() if len(r) == 5 and heavy_idx in r]

        if not five_rings:
            return False

        for ring in five_rings:
            ring_set = set(ring)

            n_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == 'N')

            if n_in_ring < 2:
                continue

            n_neighbors = 0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'N' and neighbor.GetIdx() in ring_set:
                    n_neighbors += 1

            if n_neighbors == 2:
                return True

        return False
    except Exception:
        return False

def _is_isolated_heterocyclic_H(mol: Chem.Mol, heavy_idx: int) -> bool:

    try:
        atom = mol.GetAtomWithIdx(heavy_idx)
        if atom.GetSymbol() != 'C' or not atom.GetIsAromatic():
            return False

        ring_info = mol.GetRingInfo()
        five_rings = [r for r in ring_info.AtomRings() if len(r) == 5 and heavy_idx in r]

        if not five_rings:
            return False

        for ring in five_rings:
            ring_set = set(ring)

            n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == 'N')

            if n_count < 2:
                continue

            c_with_H = 0
            for idx in ring:
                a = mol.GetAtomWithIdx(idx)
                if a.GetSymbol() == 'C':

                    has_H = any(n.GetSymbol() == 'H' or n.GetAtomicNum() == 1
                               for n in a.GetNeighbors())
                    if not has_H and a.GetTotalNumHs() > 0:
                        has_H = True
                    if has_H:
                        c_with_H += 1

            if c_with_H <= 2 and n_count >= 2:

                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in ring_set and neighbor.GetSymbol() == 'C':
                        neighbor_has_H = any(n.GetSymbol() == 'H' for n in neighbor.GetNeighbors())
                        if neighbor_has_H:

                            return False
                return True

        return False
    except Exception:
        return False

def _is_imidazole_h2_position(mol: Chem.Mol, c1_idx: int, c2_idx: int, ring: tuple) -> bool:

    try:
        for c_idx in (c1_idx, c2_idx):
            carbon = mol.GetAtomWithIdx(c_idx)
            if carbon.GetSymbol() != 'C':
                continue

            n_neighbors_in_ring = 0
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetSymbol() == 'N' and neighbor.GetIdx() in ring:
                    n_neighbors_in_ring += 1

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

    try:
        h1 = mol.GetAtomWithIdx(h1_idx)
        h2 = mol.GetAtomWithIdx(h2_idx)

        c1_list = [n for n in h1.GetNeighbors() if n.GetIsAromatic()]
        c2_list = [n for n in h2.GetNeighbors() if n.GetIsAromatic()]

        if not c1_list or not c2_list:
            return None

        c1 = c1_list[0]
        c2 = c2_list[0]
        c1_idx = c1.GetIdx()
        c2_idx = c2.GetIdx()

        if c1_idx == c2_idx:
            return None

        ring_info = mol.GetRingInfo()
        common_rings = []

        for ring in ring_info.AtomRings():
            if c1_idx in ring and c2_idx in ring:

                is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
                if is_aromatic:
                    common_rings.append(ring)

        if not common_rings:
            return None

        ring = min(common_rings, key=len)
        ring_size = len(ring)

        ring_list = list(ring)
        idx1 = ring_list.index(c1_idx)
        idx2 = ring_list.index(c2_idx)
        dist = min(abs(idx1 - idx2), ring_size - abs(idx1 - idx2))

        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        heteroatoms = [a.GetSymbol() for a in ring_atoms if a.GetSymbol() != 'C']

        cycle_type = 'benzene'

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
            return None

        if cycle_type in _AROMATIC_COUPLING_TABLE:
            cycle_J = _AROMATIC_COUPLING_TABLE[cycle_type]
            if isinstance(cycle_J, dict):

                if position in cycle_J:
                    J = cycle_J[position]

                if cycle_J.get('special_h2_singlet'):

                    if _is_imidazole_h2_position(mol, c1_idx, c2_idx, ring):
                        J = 0.0

        return (J, position)

    except Exception:
        return None

_BYSTROV_PEPTIDE_COEFFS = {

    'HNCH_alpha': {
        'A': 6.4, 'B': -1.4, 'C': 1.9,
        'phase': -60.0,
    },

    'HCCO': {
        'A': 4.0, 'B': -1.0, 'C': 0.5,
        'phase': 0.0,
    },

    'HCCH_sidechain': {
        'A': 9.5, 'B': -1.6, 'C': 1.8,
        'phase': 0.0,
    },
}

_PEPTIDE_STRUCTURE_J = {

    'helix_alpha': (3.9, -57),
    'helix_310': (4.5, -49),
    'sheet_parallel': (8.9, -119),
    'sheet_antiparallel': (8.5, -139),
    'polyproline_II': (5.5, -75),
    'turn_I': (6.0, -60),
    'turn_II': (4.0, -60),
    'random_coil': (7.0, None),
}

def _get_peptide_HNCH_coupling(
    mol: Chem.Mol,
    h_amide_idx: int,
    h_alpha_idx: int,
    mol3d: Chem.Mol = None
) -> Optional[Tuple[float, str]]:

    try:
        h_amide = mol.GetAtomWithIdx(h_amide_idx)
        h_alpha = mol.GetAtomWithIdx(h_alpha_idx)

        n_list = [n for n in h_amide.GetNeighbors() if n.GetSymbol() == 'N']
        if not n_list:
            return None
        n_atom = n_list[0]
        n_idx = n_atom.GetIdx()

        c_alpha_list = [n for n in h_alpha.GetNeighbors() if n.GetSymbol() == 'C']
        if not c_alpha_list:
            return None
        c_alpha = c_alpha_list[0]
        c_alpha_idx = c_alpha.GetIdx()

        if not mol.GetBondBetweenAtoms(n_idx, c_alpha_idx):
            return None

        if mol3d and mol3d.GetNumConformers() > 0:
            try:
                conf = mol3d.GetConformer()
                phi = AllChem.GetDihedralDeg(conf, h_amide_idx, n_idx, c_alpha_idx, h_alpha_idx)

                coeffs = _BYSTROV_PEPTIDE_COEFFS['HNCH_alpha']
                A, B, C = coeffs['A'], coeffs['B'], coeffs['C']
                phase = coeffs['phase']

                phi_rad = math.radians(phi + phase)
                J = A * (math.cos(phi_rad) ** 2) + B * math.cos(phi_rad) + C
                J = max(0.0, J)

                structure = _predict_secondary_structure_from_J(J)

                return (J, structure)

            except Exception:
                pass

        return (7.0, 'random_coil')

    except Exception:
        return None

def _predict_secondary_structure_from_J(J_HNCH: float) -> str:

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

def _detect_diastereotopic_protons(
    mol: Chem.Mol,
    ch2_carbon_idx: int,
    chiral_indices: Optional[set] = None
) -> Tuple[bool, float]:

    try:
        c_atom = mol.GetAtomWithIdx(ch2_carbon_idx)

        if c_atom.GetSymbol() != 'C':
            return (False, 0.0)

        h_count = sum(1 for n in c_atom.GetNeighbors() if n.GetSymbol() == 'H')
        if h_count != 2:
            return (False, 0.0)

        if chiral_indices is None:
            try:

                base_mol = Chem.RemoveHs(mol)
                chiral_centers_list = Chem.FindMolChiralCenters(base_mol, includeUnassigned=True)
                chiral_indices = {idx for idx, _ in chiral_centers_list}
            except Exception:
                chiral_indices = set()

        chiral_distance = float('inf')

        visited = {ch2_carbon_idx}
        queue = [(ch2_carbon_idx, 0)]

        while queue:
            current_idx, dist = queue.pop(0)

            if dist > 4:
                break

            current = mol.GetAtomWithIdx(current_idx)

            is_chiral = False
            chiral_tag = current.GetChiralTag()
            if chiral_tag != Chem.ChiralType.CHI_UNSPECIFIED:
                is_chiral = True

            if current_idx in chiral_indices:
                is_chiral = True

            if is_chiral and current_idx != ch2_carbon_idx:
                chiral_distance = min(chiral_distance, dist)

            for nb in current.GetNeighbors():
                if nb.GetIdx() not in visited and nb.GetSymbol() in ('C', 'N', 'O', 'S'):
                    visited.add(nb.GetIdx())
                    queue.append((nb.GetIdx(), dist + 1))

        if chiral_distance == float('inf'):
            return (False, 0.0)

        if chiral_distance in CH2_CHIRAL_DELTA_BY_DISTANCE:
            delta_ppm = CH2_CHIRAL_DELTA_BY_DISTANCE[chiral_distance][0]
        elif chiral_distance <= 4:

            delta_ppm = DEFAULT_DELTA.get(f"ch2_near_chiral_{min(chiral_distance, 3)}", 0.05)
        else:
            delta_ppm = 0.02

        return (True, delta_ppm)

    except Exception:
        return (False, 0.0)

_LINEWIDTH_BY_PROTON_TYPE = {

    'OH_alcohol': (5.0, 20.0),
    'OH_phenol': (3.0, 15.0),
    'OH_carboxylic': (10.0, 50.0),
    'NH_amine': (5.0, 30.0),
    'NH_amide': (2.0, 10.0),
    'NH_indole': (1.0, 5.0),
    'SH': (3.0, 15.0),

    'aromatic': (0.5, 1.5),
    'vinyl': (0.5, 1.5),
    'aldehyde': (0.5, 1.0),
    'alkyl_CH3': (0.8, 2.0),
    'alkyl_CH2': (1.0, 2.5),
    'alkyl_CH': (1.0, 2.0),
    'alpha_heteroatom': (1.0, 2.5),

    'default': (1.0, 2.0),
    'exchangeable': (5.0, 20.0),
}

_LINEWIDTH_MODIFIERS = {

    'D2O': {'exchangeable': 0.0},
    'DMSO': {'exchangeable': 0.3},
    'CDCl3': {'exchangeable': 0.8},
    'MeOD': {'exchangeable': 0.1},

    'low_temp': 1.5,
    'room_temp': 1.0,
    'high_temp': 0.7,

    'high_viscosity': 1.5,
    'low_viscosity': 0.8,
}

def _get_adaptive_4J_threshold(
    mol: Chem.Mol,
    i: int,
    j: int,
    base_threshold: float = 1.0
) -> float:

    try:
        atom_i = mol.GetAtomWithIdx(i)
        atom_j = mol.GetAtomWithIdx(j)

        if atom_i.IsInRing() and atom_j.IsInRing():
            ring_info = mol.GetRingInfo()

            for ring in ring_info.AtomRings():
                if i in ring and j in ring:
                    return 0.7

            i_rings = [set(r) for r in ring_info.AtomRings() if i in r]
            j_rings = [set(r) for r in ring_info.AtomRings() if j in r]

            for ri in i_rings:
                for rj in j_rings:
                    if ri & rj:
                        return 0.8

        path = Chem.GetShortestPath(mol, i, j)
        if path and len(path) == 5:

            mid_atoms = [mol.GetAtomWithIdx(idx) for idx in path[1:-1]]

            has_double = False
            for k in range(len(mid_atoms) - 1):
                bond = mol.GetBondBetweenAtoms(path[k+1], path[k+2])
                if bond and bond.GetBondTypeAsDouble() >= 2:
                    has_double = True
                    break

            if has_double:
                return 0.8

            if all(a.IsInRing() for a in mid_atoms):
                return 0.8

        if atom_i.GetIsAromatic() or atom_j.GetIsAromatic():
            return 0.9

        if path:
            rotatable_in_path = 0
            for k in range(len(path) - 1):
                bond = mol.GetBondBetweenAtoms(path[k], path[k+1])
                if bond:
                    a1 = mol.GetAtomWithIdx(path[k])
                    a2 = mol.GetAtomWithIdx(path[k+1])

                    if (bond.GetBondTypeAsDouble() == 1 and
                        not bond.IsInRing() and
                        a1.GetSymbol() in ('C', 'N', 'O') and
                        a2.GetSymbol() in ('C', 'N', 'O')):
                        rotatable_in_path += 1

            if rotatable_in_path >= 2:
                return 1.5
            elif rotatable_in_path == 1:
                return 1.2

        return base_threshold

    except Exception:
        return base_threshold

def _get_rotatable_bonds(mol: Chem.Mol) -> List[Tuple[int, int, int, int]]:

    rotatable = []

    try:

        pattern = Chem.MolFromSmarts('[!$([NH]!@C(=O))&!D1]-&!@[!$([NH]!@C(=O))&!D1]')
        if pattern is None:
            return rotatable

        matches = mol.GetSubstructMatches(pattern)

        for match in matches:
            a2_idx, a3_idx = match[0], match[1]
            a2 = mol.GetAtomWithIdx(a2_idx)
            a3 = mol.GetAtomWithIdx(a3_idx)

            a1_candidates = [n.GetIdx() for n in a2.GetNeighbors() if n.GetIdx() != a3_idx]

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

    if angles is None:
        angles = [60.0, 180.0, 300.0]

    mol_h = Chem.AddHs(mol, addCoords=True)
    rotatable = _get_rotatable_bonds(mol_h)
    n_rot = len(rotatable)

    results = []

    try:

        if 0 < n_rot <= max_rotatable:
            import itertools

            n_combos = len(angles) ** n_rot

            if n_combos > n_conformers * 3:

                combos = []
                import random
                for _ in range(n_conformers * 3):
                    combo = tuple(random.choice(angles) for _ in range(n_rot))
                    combos.append(combo)
            else:
                combos = list(itertools.product(angles, repeat=n_rot))

            ps = AllChem.ETKDGv3() if hasattr(AllChem, 'ETKDGv3') else AllChem.ETKDG()
            ps.randomSeed = 42
            AllChem.EmbedMolecule(mol_h, ps)

            if mol_h.GetNumConformers() == 0:
                return results

            for combo in combos:

                mol_copy = Chem.Mol(mol_h)
                conf = mol_copy.GetConformer()

                for (a1, a2, a3, a4), angle in zip(rotatable, combo):
                    try:
                        AllChem.SetDihedralDeg(conf, a1, a2, a3, a4, angle)
                    except Exception:
                        continue

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

                    new_conf_id = mol_h.AddConformer(conf, assignId=True)
                    results.append((new_conf_id, energy))

        if n_rot > max_rotatable or not results:
            ps = AllChem.ETKDGv3() if hasattr(AllChem, 'ETKDGv3') else AllChem.ETKDG()
            ps.randomSeed = 42
            ps.numThreads = 0
            ps.pruneRmsThresh = 0.35
            ps.useRandomCoords = True

            n_generate = n_conformers * 3
            conf_ids = AllChem.EmbedMultipleConfs(mol_h, numConfs=n_generate, params=ps)

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

        if results:
            results.sort(key=lambda x: x[1])
            min_e = results[0][1]
            results = [(cid, e - min_e) for cid, e in results if e - min_e <= energy_cutoff_kcal]
            results = results[:n_conformers]

        return results

    except Exception:
        return results

def _is_prochiral_center(
    mol: Chem.Mol,
    atom_idx: int
) -> bool:

    try:
        atom = mol.GetAtomWithIdx(atom_idx)

        if atom.GetSymbol() == 'C':
            hyb = atom.GetHybridization()

            if hyb == Chem.HybridizationType.SP2:

                double_bonded = None
                for bond in atom.GetBonds():
                    if bond.GetBondTypeAsDouble() >= 2:
                        other = bond.GetOtherAtom(atom)
                        if other.GetSymbol() in ('O', 'N', 'C', 'S'):
                            double_bonded = other
                            break

                if double_bonded:

                    other_neighbors = [n for n in atom.GetNeighbors()
                                       if n.GetIdx() != double_bonded.GetIdx()]

                    if len(other_neighbors) == 2:

                        try:
                            ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                            r1 = ranks[other_neighbors[0].GetIdx()]
                            r2 = ranks[other_neighbors[1].GetIdx()]
                            if r1 != r2:
                                return True
                        except Exception:
                            pass

        return False

    except Exception:
        return False

def _is_prochiral_CH2(
    mol: Chem.Mol,
    ch2_idx: int
) -> Tuple[bool, str]:

    try:
        ch2 = mol.GetAtomWithIdx(ch2_idx)

        if ch2.GetSymbol() != 'C':
            return (False, '')

        h_count = sum(1 for n in ch2.GetNeighbors() if n.GetSymbol() == 'H')
        if h_count != 2:
            return (False, '')

        neighbors_heavy = [n for n in ch2.GetNeighbors() if n.GetSymbol() != 'H']

        for nb in neighbors_heavy:
            if nb.GetSymbol() == 'C':

                for bond in nb.GetBonds():
                    if bond.GetBondTypeAsDouble() >= 2:
                        other = bond.GetOtherAtom(nb)
                        if other.GetSymbol() == 'O':

                            carbonyl_neighbors = [n for n in nb.GetNeighbors()
                                                  if n.GetIdx() != other.GetIdx()
                                                  and n.GetIdx() != ch2_idx]

                            if len(carbonyl_neighbors) >= 1:

                                try:
                                    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))

                                    carbonyl_sub_ranks = [ranks[n.GetIdx()] for n in carbonyl_neighbors]
                                    ch2_rank = ranks[ch2_idx]

                                    if len(set(carbonyl_sub_ranks)) > 1 or (carbonyl_sub_ranks and carbonyl_sub_ranks[0] != ch2_rank):
                                        return (True, 'adjacent_to_asymmetric_carbonyl')
                                except Exception:
                                    pass

        for nb in neighbors_heavy:
            if nb.GetSymbol() == 'C':
                for bond in nb.GetBonds():
                    if bond.GetBondTypeAsDouble() >= 2:
                        other = bond.GetOtherAtom(nb)
                        if other.GetSymbol() == 'C':

                            nb_subs = [n for n in nb.GetNeighbors()
                                       if n.GetIdx() not in (other.GetIdx(), ch2_idx)]
                            other_subs = [n for n in other.GetNeighbors()
                                          if n.GetIdx() != nb.GetIdx()]

                            try:
                                ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                                nb_ranks = set(ranks[n.GetIdx()] for n in nb_subs)
                                other_ranks = set(ranks[n.GetIdx()] for n in other_subs)

                                if len(nb_ranks) > 1 or len(other_ranks) > 1:
                                    return (True, 'adjacent_to_asymmetric_alkene')
                            except Exception:
                                pass

        if ch2.IsInRing():

            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if ch2_idx in ring:

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

        for nb in neighbors_heavy:
            if nb.GetIsAromatic():

                ar_neighbors = [n for n in nb.GetNeighbors()
                               if n.GetIdx() != ch2_idx and n.GetSymbol() != 'H']
                if len(ar_neighbors) >= 2:
                    try:
                        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                        ar_sub_ranks = [ranks[n.GetIdx()] for n in ar_neighbors]

                        if len(set(ar_sub_ranks)) > 1:
                            return (True, 'benzylic_asymmetric_aromatic')
                    except Exception:
                        pass

                if nb.IsInRing():
                    ring_info = mol.GetRingInfo()
                    nb_rings = []
                    for ring in ring_info.AtomRings():
                        if nb.GetIdx() in ring:
                            nb_rings.append(ring)

                    if len(nb_rings) >= 2:
                        return (True, 'benzylic_fused_ring')

        for nb in neighbors_heavy:
            if nb.GetSymbol() in ('O', 'N', 'S'):
                hetero_neighbors = [n for n in nb.GetNeighbors()
                                   if n.GetIdx() != ch2_idx and n.GetSymbol() != 'H']

                if hetero_neighbors:

                    try:
                        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                        hetero_sub_ranks = [ranks[n.GetIdx()] for n in hetero_neighbors]

                        if len(hetero_sub_ranks) >= 2 and len(set(hetero_sub_ranks)) > 1:
                            return (True, f'alpha_to_asymmetric_{nb.GetSymbol()}')

                        for hetero_nb in hetero_neighbors:
                            if hetero_nb.GetSymbol() == 'C':

                                base_mol = Chem.RemoveHs(mol)
                                chiral_list = Chem.FindMolChiralCenters(base_mol, includeUnassigned=True)
                                chiral_indices = {idx for idx, _ in chiral_list}
                                if hetero_nb.GetIdx() in chiral_indices:
                                    return (True, f'alpha_to_{nb.GetSymbol()}_chiral')
                    except Exception:
                        pass

        if ch2.IsInRing():
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if ch2_idx in ring:
                    ring_size = len(ring)

                    if ring_size in (3, 4):

                        ring_set = set(ring)
                        has_substituent = False
                        has_heteroatom = False

                        for r_idx in ring:
                            r_atom = mol.GetAtomWithIdx(r_idx)

                            if r_atom.GetSymbol() in ('O', 'N', 'S'):
                                has_heteroatom = True

                            for nb_atom in r_atom.GetNeighbors():
                                if nb_atom.GetIdx() not in ring_set and nb_atom.GetSymbol() != 'H':
                                    has_substituent = True

                        if has_heteroatom:
                            if ring_size == 3:
                                return (True, 'strained_ring_3_hetero')
                            else:
                                return (True, 'strained_ring_4_hetero')

                        if has_substituent:
                            return (True, f'strained_ring_{ring_size}_substituted')

        if ch2.IsInRing():
            ring_info = mol.GetRingInfo()
            rings_containing_ch2 = []
            for ring in ring_info.AtomRings():
                if ch2_idx in ring:
                    rings_containing_ch2.append(set(ring))

            if len(rings_containing_ch2) >= 2:
                return (True, 'bicyclic_bridge')

            if len(rings_containing_ch2) == 1:
                my_ring = rings_containing_ch2[0]

                for r_idx in my_ring:
                    if r_idx != ch2_idx:
                        for other_ring in ring_info.AtomRings():
                            other_set = set(other_ring)
                            if r_idx in other_set and other_set != my_ring:

                                return (True, 'fused_or_bridged_ring')

        if ch2.IsInRing():
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if ch2_idx in ring:
                    ring_size = len(ring)
                    ring_set = set(ring)

                    if ring_size in (5, 6, 7, 8):
                        substituent_count = 0
                        substituent_positions = []

                        for r_idx in ring:
                            r_atom = mol.GetAtomWithIdx(r_idx)
                            for nb_atom in r_atom.GetNeighbors():
                                if nb_atom.GetIdx() not in ring_set and nb_atom.GetSymbol() != 'H':
                                    substituent_count += 1
                                    substituent_positions.append(r_idx)

                        if substituent_count >= 1:

                            if substituent_positions:

                                sub_pos = substituent_positions[0]
                                ring_list = list(ring)
                                try:
                                    my_pos = ring_list.index(ch2_idx)
                                    sub_idx_in_ring = ring_list.index(sub_pos)
                                    dist_in_ring = min(
                                        abs(my_pos - sub_idx_in_ring),
                                        ring_size - abs(my_pos - sub_idx_in_ring)
                                    )

                                    if dist_in_ring <= 3:
                                        return (True, f'ring_{ring_size}_monosubstituted')
                                except ValueError:
                                    pass

        for nb in neighbors_heavy:
            if nb.GetSymbol() == 'C' and nb.GetIsAromatic():

                ar_neighbors = [n for n in nb.GetNeighbors()
                               if n.GetIsAromatic() and n.GetIdx() != ch2_idx]
                if ar_neighbors:

                    for ar_nb in ar_neighbors:
                        ar_ring_count = 0
                        for ring in mol.GetRingInfo().AtomRings():
                            if ar_nb.GetIdx() in ring:
                                ar_ring_count += 1

                        if ar_ring_count >= 1:

                            ortho_subs = []
                            for ar_n in nb.GetNeighbors():
                                if ar_n.GetIsAromatic() and ar_n.GetIdx() != ar_nb.GetIdx():
                                    for ar_nn in ar_n.GetNeighbors():
                                        if not ar_nn.GetIsAromatic() and ar_nn.GetSymbol() != 'H':
                                            ortho_subs.append(ar_nn.GetIdx())
                            if len(ortho_subs) >= 2:
                                return (True, 'axial_chirality_biaryl')

        for nb in neighbors_heavy:
            if nb.GetSymbol() == 'P':

                p_neighbors = [n for n in nb.GetNeighbors() if n.GetIdx() != ch2_idx]
                if len(p_neighbors) >= 3:
                    try:
                        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                        p_sub_ranks = [ranks[n.GetIdx()] for n in p_neighbors]

                        if len(set(p_sub_ranks)) >= 3:
                            return (True, 'alpha_to_chiral_phosphorus')

                        if len(set(p_sub_ranks)) >= 2:
                            return (True, 'alpha_to_asymmetric_phosphorus')
                    except Exception:
                        pass

        for nb in neighbors_heavy:
            if nb.GetSymbol() == 'C':
                for bond in nb.GetBonds():
                    if bond.GetBondTypeAsDouble() >= 2:
                        other = bond.GetOtherAtom(nb)
                        if other.GetSymbol() == 'C':

                            stereo = bond.GetStereo()
                            if stereo in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ,
                                         Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOTRANS):
                                return (True, 'adjacent_to_EZ_alkene')

        if ch2.IsInRing():
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if ch2_idx in ring:
                    ring_set = set(ring)

                    has_carbonyl_in_ring = False
                    has_N_in_ring = False
                    has_O_in_ring = False

                    for r_idx in ring:
                        r_atom = mol.GetAtomWithIdx(r_idx)
                        if r_atom.GetSymbol() == 'N':
                            has_N_in_ring = True
                        if r_atom.GetSymbol() == 'O' and not any(
                            b.GetBondTypeAsDouble() >= 2 for b in r_atom.GetBonds()
                        ):
                            has_O_in_ring = True
                        if r_atom.GetSymbol() == 'C':
                            for bond in r_atom.GetBonds():
                                if bond.GetBondTypeAsDouble() >= 2:
                                    other_atom = bond.GetOtherAtom(r_atom)
                                    if other_atom.GetSymbol() == 'O':
                                        has_carbonyl_in_ring = True

                    if has_carbonyl_in_ring and has_N_in_ring:
                        return (True, 'lactam_ring')

                    if has_carbonyl_in_ring and has_O_in_ring:
                        return (True, 'lactone_ring')

        chiral_neighbor_count = 0
        try:
            base_mol = Chem.RemoveHs(mol)
            chiral_list = Chem.FindMolChiralCenters(base_mol, includeUnassigned=True)
            chiral_indices = {idx for idx, _ in chiral_list}

            for nb in neighbors_heavy:
                if nb.GetIdx() in chiral_indices:
                    chiral_neighbor_count += 1

            if chiral_neighbor_count >= 2:
                return (True, 'between_two_chiral_centers')
        except Exception:
            pass

        if ch2.IsInRing():
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if ch2_idx in ring:
                    ring_size = len(ring)

                    if ring_size >= 12:
                        return (True, 'macrocycle')

                    if 8 <= ring_size <= 11:
                        ring_set = set(ring)
                        has_hetero = any(
                            mol.GetAtomWithIdx(r).GetSymbol() in ('O', 'N', 'S')
                            for r in ring
                        )
                        if has_hetero:
                            return (True, 'medium_ring_hetero')

        for nb in neighbors_heavy:
            if nb.GetSymbol() == 'S':

                s_neighbors = [n for n in nb.GetNeighbors() if n.GetIdx() != ch2_idx]

                has_sulfoxide = any(
                    bond.GetBondTypeAsDouble() >= 2 and bond.GetOtherAtom(nb).GetSymbol() == 'O'
                    for bond in nb.GetBonds()
                )
                if has_sulfoxide:
                    return (True, 'alpha_to_chiral_sulfoxide')

                if nb.GetFormalCharge() > 0:
                    return (True, 'alpha_to_sulfonium')

                if len(s_neighbors) >= 2:
                    try:
                        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                        s_sub_ranks = [ranks[n.GetIdx()] for n in s_neighbors]
                        if len(set(s_sub_ranks)) >= 2:
                            return (True, 'alpha_to_asymmetric_sulfur')
                    except Exception:
                        pass

        for nb in neighbors_heavy:
            if nb.GetSymbol() == 'Si':
                si_neighbors = [n for n in nb.GetNeighbors() if n.GetIdx() != ch2_idx]
                if len(si_neighbors) >= 3:
                    try:
                        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                        si_sub_ranks = [ranks[n.GetIdx()] for n in si_neighbors]
                        if len(set(si_sub_ranks)) >= 3:
                            return (True, 'alpha_to_chiral_silicon')
                    except Exception:
                        pass

        if ch2.IsInRing():
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if ch2_idx in ring:
                    ring_set = set(ring)

                    for r_idx in ring:
                        r_atom = mol.GetAtomWithIdx(r_idx)
                        if r_atom.GetSymbol() == 'C':
                            rings_for_atom = sum(1 for r in ring_info.AtomRings() if r_idx in r)
                            if rings_for_atom >= 2:

                                for other_ring in ring_info.AtomRings():
                                    if r_idx in other_ring and set(other_ring) != ring_set:
                                        common = ring_set & set(other_ring)
                                        if len(common) == 1:
                                            return (True, 'spiro_ring')

        for nb in neighbors_heavy:
            if nb.GetSymbol() == 'C':

                double_bond_count = sum(
                    1 for bond in nb.GetBonds() if bond.GetBondTypeAsDouble() >= 2
                )
                if double_bond_count >= 2:

                    return (True, 'adjacent_to_allene')

        for nb in neighbors_heavy:
            if nb.GetSymbol() == 'N' and nb.GetFormalCharge() > 0:
                n_neighbors = [n for n in nb.GetNeighbors() if n.GetIdx() != ch2_idx]
                if len(n_neighbors) >= 3:
                    try:
                        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                        n_sub_ranks = [ranks[n.GetIdx()] for n in n_neighbors]
                        if len(set(n_sub_ranks)) >= 3:
                            return (True, 'alpha_to_chiral_quaternary_N')
                        if len(set(n_sub_ranks)) >= 2:
                            return (True, 'alpha_to_asymmetric_quaternary_N')
                    except Exception:
                        pass

        for nb in neighbors_heavy:
            if nb.GetSymbol() == 'B':
                b_neighbors = [n for n in nb.GetNeighbors() if n.GetIdx() != ch2_idx]
                if len(b_neighbors) >= 2:
                    try:
                        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                        b_sub_ranks = [ranks[n.GetIdx()] for n in b_neighbors]
                        if len(set(b_sub_ranks)) >= 2:
                            return (True, 'alpha_to_asymmetric_boron')
                    except Exception:
                        pass

        o_or_s_neighbors = [nb for nb in neighbors_heavy if nb.GetSymbol() in ('O', 'S')]
        if len(o_or_s_neighbors) >= 2:

            try:
                ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                o_s_ranks = [ranks[nb.GetIdx()] for nb in o_or_s_neighbors]
                if len(set(o_s_ranks)) >= 2:
                    return (True, 'acetal_asymmetric')

                for o_s_atom in o_or_s_neighbors:
                    o_s_other_neighbors = [n for n in o_s_atom.GetNeighbors()
                                          if n.GetIdx() != ch2_idx]
                    if o_s_other_neighbors:
                        other_ranks = [ranks[n.GetIdx()] for n in o_s_other_neighbors]
                        if len(set(other_ranks)) > 1:
                            return (True, 'acetal_asymmetric_substituents')
            except Exception:
                pass

        if ch2.IsInRing():
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if ch2_idx in ring:
                    ring_size = len(ring)

                    if ring_size == 5:
                        ring_set = set(ring)
                        hetero_in_ring = [
                            mol.GetAtomWithIdx(r) for r in ring
                            if mol.GetAtomWithIdx(r).GetSymbol() in ('O', 'N', 'S')
                        ]
                        if hetero_in_ring:

                            return (True, 'saturated_5_ring_hetero')

        for nb in neighbors_heavy:
            if nb.GetSymbol() == 'C':

                nb_heavy_neighbors = [n for n in nb.GetNeighbors() if n.GetSymbol() != 'H']

                if len(nb_heavy_neighbors) == 4:

                    has_double_bond = any(
                        bond.GetBondTypeAsDouble() >= 2
                        for bond in nb.GetBonds()
                    )

                    if not has_double_bond:

                        try:
                            ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                            quat_sub_ranks = [ranks[n.GetIdx()] for n in nb_heavy_neighbors]
                            unique_ranks = set(quat_sub_ranks)

                            if len(unique_ranks) >= 3:
                                return (True, 'adjacent_to_prochiral_quaternary_C')

                            if len(unique_ranks) == 2:
                                ch2_rank = ranks[ch2_idx]
                                other_ranks = [r for r in quat_sub_ranks if r != ch2_rank]

                                if len(set(other_ranks)) >= 2 or (other_ranks and other_ranks[0] != ch2_rank):
                                    return (True, 'adjacent_to_prochiral_quaternary_C')
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

    is_prochiral, reason = _is_prochiral_CH2(mol, heavy_idx)

    if is_prochiral:
        return True

    heavy_atom = mol.GetAtomWithIdx(heavy_idx)
    for nb in heavy_atom.GetNeighbors():
        if nb.GetSymbol() != 'H':
            if _is_prochiral_center(mol, nb.GetIdx()):
                return True

    return False

def _get_geminal_2J(heavy_atom, mol: Chem.Mol) -> float:

    sym = heavy_atom.GetSymbol()

    if sym != 'C':
        return 12.5

    hyb = str(heavy_atom.GetHybridization())
    is_sp2 = 'SP2' in hyb or heavy_atom.GetIsAromatic()

    if is_sp2:

        for bond in heavy_atom.GetBonds():
            if bond.GetBondTypeAsDouble() >= 2:
                other = bond.GetOtherAtom(heavy_atom)
                if other.GetSymbol() == 'O':
                    return 2.5

        n_double_bonds = sum(1 for b in heavy_atom.GetBonds() if b.GetBondTypeAsDouble() >= 2)
        if n_double_bonds >= 2:
            return 7.5
        return 2.0

    _GEMINAL_2J_DIRECT_TABLE = {

        'CH2_indolyl':     14.5,
        'CH2_benzyl':      14.0,
        'CH2_imidazolyl':  14.5,

        'CH2_allyl':       14.5,
        'CH2_allyl_aryl':  15.0,

        'CH2_alpha_CO':    16.5,
        'CH2_di_carbonyl': 17.5,

        'CH2_alpha_N':     13.0,
        'CH2_alpha_O':     10.0,
        'CH2_alpha_S':     11.0,
        'CH2_acetal':       6.0,

        'CH2_cyclopropane': 4.5,
        'CH2_epoxide':      5.5,
        'CH2_aziridine':    6.0,
    }

    neighbors = [n for n in heavy_atom.GetNeighbors() if n.GetSymbol() != 'H']
    neighbor_symbols = [n.GetSymbol() for n in neighbors]

    ch2_type = None

    if neighbor_symbols.count('O') >= 2:
        ch2_type = 'CH2_acetal'

    elif heavy_atom.IsInRingSize(3):
        ring_atoms = []
        for ring in mol.GetRingInfo().AtomRings():
            if len(ring) == 3 and heavy_atom.GetIdx() in ring:
                ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                break
        if any(a.GetSymbol() == 'O' for a in ring_atoms):
            ch2_type = 'CH2_epoxide'
        elif any(a.GetSymbol() == 'N' for a in ring_atoms):
            ch2_type = 'CH2_aziridine'
        else:
            ch2_type = 'CH2_cyclopropane'

    else:
        has_aromatic = False
        has_indole = False
        has_imidazole = False
        has_vinyl = False
        has_carbonyl = False
        n_carbonyl = 0
        has_allyl_on_aryl = False

        for nb in neighbors:
            if nb.GetSymbol() != 'C':
                continue

            if nb.GetIsAromatic():
                has_aromatic = True

                ri = mol.GetRingInfo()
                for ring in ri.AtomRings():
                    if nb.GetIdx() in ring:

                        ring_N_count = sum(1 for i in ring
                                           if mol.GetAtomWithIdx(i).GetSymbol() == 'N'
                                           and mol.GetAtomWithIdx(i).GetIsAromatic())
                        if ring_N_count >= 1:

                            if len(ring) == 5 and ring_N_count == 1:
                                has_indole = True

                            elif len(ring) == 5 and ring_N_count >= 2:
                                has_imidazole = True

            elif 'SP2' in str(nb.GetHybridization()):
                for b in nb.GetBonds():
                    if b.GetBondTypeAsDouble() >= 2:
                        other = b.GetOtherAtom(nb)
                        if other.GetSymbol() == 'O':
                            has_carbonyl = True
                            n_carbonyl += 1
                        elif other.GetSymbol() == 'C':
                            has_vinyl = True

                            for nb3 in other.GetNeighbors():
                                if nb3.GetIsAromatic():
                                    has_allyl_on_aryl = True
                        break

        if has_indole:
            ch2_type = 'CH2_indolyl'
        elif has_imidazole:
            ch2_type = 'CH2_imidazolyl'
        elif has_aromatic:
            ch2_type = 'CH2_benzyl'
        elif n_carbonyl >= 2:
            ch2_type = 'CH2_di_carbonyl'
        elif has_carbonyl:
            ch2_type = 'CH2_alpha_CO'
        elif has_allyl_on_aryl:
            ch2_type = 'CH2_allyl_aryl'
        elif has_vinyl:
            ch2_type = 'CH2_allyl'
        elif 'O' in neighbor_symbols:
            ch2_type = 'CH2_alpha_O'
        elif 'N' in neighbor_symbols:
            ch2_type = 'CH2_alpha_N'
        elif 'S' in neighbor_symbols:
            ch2_type = 'CH2_alpha_S'

    if ch2_type and ch2_type in _GEMINAL_2J_DIRECT_TABLE:
        return _GEMINAL_2J_DIRECT_TABLE[ch2_type]

    base_2J = -12.5
    increment = 0.0

    neighbors = [n for n in heavy_atom.GetNeighbors() if n.GetSymbol() != 'H']
    neighbor_symbols = [n.GetSymbol() for n in neighbors]

    if heavy_atom.IsInRingSize(3):

        ring_atoms = []
        for ring in mol.GetRingInfo().AtomRings():
            if len(ring) == 3 and heavy_atom.GetIdx() in ring:
                ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                break
        has_O_in_ring = any(a.GetSymbol() == 'O' for a in ring_atoms)
        has_N_in_ring = any(a.GetSymbol() == 'N' for a in ring_atoms)

        if has_O_in_ring:
            return 5.5
        elif has_N_in_ring:
            return 6.0
        else:
            return 4.5

    if heavy_atom.IsInRingSize(4) and not heavy_atom.IsInRingSize(5):
        base_2J = -12.0

    if heavy_atom.IsInRingSize(5) and not heavy_atom.IsInRingSize(6):
        base_2J = -13.0

    n_O_neighbors = neighbor_symbols.count('O')
    n_N_neighbors = neighbor_symbols.count('N')
    n_S_neighbors = neighbor_symbols.count('S')
    n_P_neighbors = neighbor_symbols.count('P')

    if n_O_neighbors >= 2:
        return 6.0

    if n_O_neighbors > 0:
        increment += _GEMINAL_2J_INCREMENTS.get('alpha_O', +2.5) * n_O_neighbors
    if n_N_neighbors > 0:
        increment += _GEMINAL_2J_INCREMENTS.get('alpha_N', -0.5) * n_N_neighbors
    if n_S_neighbors > 0:
        increment += _GEMINAL_2J_INCREMENTS.get('alpha_S', +1.5) * n_S_neighbors
    if n_P_neighbors > 0:
        increment += _GEMINAL_2J_INCREMENTS.get('alpha_P', +3.5) * n_P_neighbors

    n_aromatic = 0
    n_vinyl = 0
    n_carbonyl = 0
    n_nitro = 0

    for nb in neighbors:
        if nb.GetSymbol() == 'N':

            n_O_on_N = sum(1 for nn in nb.GetNeighbors() if nn.GetSymbol() == 'O')
            if n_O_on_N >= 2:
                n_nitro += 1
                continue

        if nb.GetSymbol() != 'C':
            continue

        nb_hyb = str(nb.GetHybridization())
        nb_is_sp2 = 'SP2' in nb_hyb or nb.GetIsAromatic()

        if nb.GetIsAromatic():
            n_aromatic += 1
        elif nb_is_sp2:

            for b in nb.GetBonds():
                if b.GetBondTypeAsDouble() >= 2:
                    other = b.GetOtherAtom(nb)
                    if other.GetSymbol() == 'O':
                        n_carbonyl += 1
                        break
                    elif other.GetSymbol() == 'C':
                        n_vinyl += 1
                        break

    if n_aromatic > 0:
        increment += _GEMINAL_2J_INCREMENTS.get('aromatic', +1.5) * n_aromatic
    if n_vinyl > 0:
        increment += _GEMINAL_2J_INCREMENTS.get('vinyl', -1.25) * n_vinyl
    if n_carbonyl > 0:
        increment += _GEMINAL_2J_INCREMENTS.get('carbonyl', -1.5) * n_carbonyl
    if n_nitro > 0:
        increment += _GEMINAL_2J_INCREMENTS.get('nitro', -5.5) * n_nitro

    if n_carbonyl == 0:
        n_beta_carbonyl = 0
        for nb in neighbors:
            if nb.GetSymbol() != 'C':
                continue

            for nb2 in nb.GetNeighbors():
                if nb2.GetIdx() == heavy_atom.GetIdx():
                    continue
                if nb2.GetSymbol() != 'C':
                    continue
                nb2_hyb = str(nb2.GetHybridization())
                if 'SP2' in nb2_hyb:

                    for b in nb2.GetBonds():
                        if b.GetBondTypeAsDouble() >= 2:
                            other = b.GetOtherAtom(nb2)
                            if other.GetSymbol() == 'O':
                                n_beta_carbonyl += 1
                                break
                    if n_beta_carbonyl > 0:
                        break
            if n_beta_carbonyl > 0:
                break

        if n_beta_carbonyl > 0:

            if n_aromatic > 0:

                pass
            elif n_vinyl > 0 or n_carbonyl > 0:

                increment += _GEMINAL_2J_INCREMENTS.get('beta_carbonyl', -2.5) * 0.4
            else:

                increment += _GEMINAL_2J_INCREMENTS.get('beta_carbonyl', -2.5)

    final_2J = base_2J + increment

    final_2J = max(-19.0, min(-4.0, final_2J))

    return abs(round(final_2J, 1))

def _get_1JCH(heavy_atom, mol: Chem.Mol) -> float:

    sym = heavy_atom.GetSymbol()

    if sym != 'C':
        return 140.0

    hyb = str(heavy_atom.GetHybridization())

    if 'SP' in hyb and 'SP2' not in hyb and 'SP3' not in hyb:
        return _13C_1JCH_TABLE.get('sp_alkyne', 250.0)

    if 'SP2' in hyb or heavy_atom.GetIsAromatic():

        for bond in heavy_atom.GetBonds():
            if bond.GetBondTypeAsDouble() >= 2:
                other = bond.GetOtherAtom(heavy_atom)
                if other.GetSymbol() == 'O':
                    return _13C_1JCH_TABLE.get('sp2_aldehyde', 175.0)

        if heavy_atom.GetIsAromatic():
            return _13C_1JCH_TABLE.get('aromatic', 160.0)

        return _13C_1JCH_TABLE.get('sp2_alkene', 155.0)

    neighbors = [n.GetSymbol() for n in heavy_atom.GetNeighbors() if n.GetSymbol() != 'H']

    if heavy_atom.IsInRingSize(3):
        return _13C_1JCH_TABLE.get('sp3_cyclopropane', 160.0)

    if 'O' in neighbors:
        return _13C_1JCH_TABLE.get('sp3_alpha_O', 140.0)

    if 'N' in neighbors:
        return _13C_1JCH_TABLE.get('sp3_alpha_N', 135.0)

    return _13C_1JCH_TABLE.get('sp3_alkyl', 125.0)

def _refined_second_order_detection(
    J_vals: List[float],
    ppm_vals: List[float],
    spectrometer_mhz: float,
    config: "SimConfig"
) -> Tuple[bool, str]:

    if not J_vals:
        return False, "insufficient_data"
    if len(ppm_vals) < 2:
        return False, "insufficient_data"

    Jmax = max(J_vals)

    if len(ppm_vals) >= 2:
        obs_shift = ppm_vals[0]
        partner_shifts = set(round(p, 4) for p in ppm_vals[1:])
        if partner_shifts:
            delta_ppm = min(abs(obs_shift - ps) for ps in partner_shifts)
        else:
            delta_ppm = 0.0
    else:
        delta_ppm = 0.0

    delta_nu_hz = delta_ppm * spectrometer_mhz

    if delta_ppm < 0.001:

        if len(set(round(p, 3) for p in ppm_vals)) == 1:

            return False, "first_order"
        else:

            return True, "strong_coupling"

    if delta_nu_hz < 0.5:
        return True, "strong_coupling"

    ratio = Jmax / delta_nu_hz

    if ratio > 1.0:
        return True, "strong_coupling"
    elif ratio > 0.25:
        return True, "moderate_second_order"
    elif ratio > 0.1:
        return True, "weak_second_order"
    else:
        return False, "first_order"

def _calculate_roof_effect_intensities(
    J_hz: float,
    delta_nu_hz: float,
    base_multiplicity: str
) -> List[Tuple[float, float]]:

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

@dataclass
class SpinSystemPeaks:

    frequencies_hz: np.ndarray
    intensities: np.ndarray
    system_type: str
    coupling_constants: Dict[str, float]
    chemical_shifts_hz: Dict[str, float]

def _solve_AB_system(
    delta_nu_hz: float,
    J_ab_hz: float,
    center_hz: float = 0.0
) -> SpinSystemPeaks:

    dv = abs(delta_nu_hz)
    J = J_ab_hz

    C = np.sqrt(dv**2 + J**2)

    if dv > 1e-6:
        theta = 0.5 * np.arctan2(J, dv)
    else:
        theta = np.pi / 4

    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    f1 = -0.5 * (C + J)
    f2 = -0.5 * (C - J)
    f3 = +0.5 * (C - J)
    f4 = +0.5 * (C + J)

    sin_2theta = 2 * sin_theta * cos_theta
    I_inner = 1.0 + sin_2theta
    I_outer = 1.0 - sin_2theta

    frequencies = np.array([f1, f2, f3, f4]) + center_hz
    intensities = np.array([I_outer, I_inner, I_inner, I_outer])

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

    frequencies = []
    intensities = []

    nu_A_eff_alpha = nu_A_hz + J_AX_hz / 2
    nu_B_eff_alpha = nu_B_hz + J_BX_hz / 2
    delta_nu_alpha = nu_A_eff_alpha - nu_B_eff_alpha
    center_alpha = (nu_A_eff_alpha + nu_B_eff_alpha) / 2

    ab_alpha = _solve_AB_system(delta_nu_alpha, J_AB_hz, center_alpha)
    frequencies.extend(ab_alpha.frequencies_hz)
    intensities.extend(ab_alpha.intensities * 0.5)

    nu_A_eff_beta = nu_A_hz - J_AX_hz / 2
    nu_B_eff_beta = nu_B_hz - J_BX_hz / 2
    delta_nu_beta = nu_A_eff_beta - nu_B_eff_beta
    center_beta = (nu_A_eff_beta + nu_B_eff_beta) / 2

    ab_beta = _solve_AB_system(delta_nu_beta, J_AB_hz, center_beta)
    frequencies.extend(ab_beta.frequencies_hz)
    intensities.extend(ab_beta.intensities * 0.5)

    delta_J = (J_AX_hz + J_BX_hz) / 2
    diff_J = (J_AX_hz - J_BX_hz) / 2

    x_freqs = [
        nu_X_hz + delta_J,
        nu_X_hz + diff_J,
        nu_X_hz - diff_J,
        nu_X_hz - delta_J,
    ]
    x_intensities = [0.25, 0.25, 0.25, 0.25]

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

    n_states = 16

    def Iz(state, spin_idx):

        bit = (state >> (3 - spin_idx)) & 1
        return 0.5 if bit == 0 else -0.5

    def flip(state, i, j):

        bi = (state >> (3 - i)) & 1
        bj = (state >> (3 - j)) & 1
        if bi != bj:
            new_state = state ^ (1 << (3 - i)) ^ (1 << (3 - j))
            return new_state
        return None

    H = np.zeros((n_states, n_states), dtype=np.float64)

    nu_A = center_hz + delta_nu_hz / 2
    nu_B = center_hz - delta_nu_hz / 2

    for s in range(n_states):

        H[s, s] -= nu_A * Iz(s, 0)
        H[s, s] -= nu_A * Iz(s, 1)
        H[s, s] -= nu_B * Iz(s, 2)
        H[s, s] -= nu_B * Iz(s, 3)

    couplings = [
        (0, 1, J_AA_hz),
        (2, 3, J_BB_hz),
        (0, 2, J_AB_hz),
        (1, 3, J_AB_hz),
        (0, 3, J_AB_prime_hz),
        (1, 2, J_AB_prime_hz),
    ]

    for i, j, J in couplings:
        for s in range(n_states):

            H[s, s] += J * Iz(s, i) * Iz(s, j)

            s_flip = flip(s, i, j)
            if s_flip is not None:
                H[s, s_flip] += J / 2
                H[s_flip, s] += J / 2

    eigenvalues, eigenvectors = np.linalg.eigh(H)

    def total_Mz(state):
        return sum(Iz(state, i) for i in range(4))

    frequencies = []
    intensities = []

    for i in range(n_states):
        for j in range(n_states):

            dE = eigenvalues[i] - eigenvalues[j]
            if dE <= 0:
                continue

            intensity = 0.0
            for s in range(n_states):
                for spin in range(4):

                    bit = (s >> (3 - spin)) & 1
                    if bit == 1:
                        s_new = s ^ (1 << (3 - spin))

                        intensity += eigenvectors[s_new, i] * eigenvectors[s, j]

            intensity = abs(intensity) ** 2

            if intensity > 1e-6:
                frequencies.append(dE)
                intensities.append(intensity)

    if not frequencies:

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

@dataclass
class SimulatedSpectrum:

    ppm_axis: np.ndarray
    freq_axis_hz: np.ndarray
    intensity: np.ndarray
    peaks_ppm: List[float]
    peaks_intensity: List[float]
    spectrometer_mhz: float
    line_width_hz: float

    @property
    def spectrum(self) -> np.ndarray:

        return self.intensity

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ppm_axis": self.ppm_axis.tolist(),
            "freq_axis_hz": self.freq_axis_hz.tolist(),
            "intensity": self.intensity.tolist(),
            "spectrum": self.intensity.tolist(),
            "peaks_ppm": self.peaks_ppm,
            "peaks_intensity": self.peaks_intensity,
            "spectrometer_mhz": self.spectrometer_mhz,
            "line_width_hz": self.line_width_hz,
        }

def _lorentzian(x: np.ndarray, x0: float, gamma: float, amplitude: float = 1.0) -> np.ndarray:

    return amplitude * (gamma ** 2) / ((x - x0) ** 2 + gamma ** 2)

def _gaussian(x: np.ndarray, x0: float, sigma: float, amplitude: float = 1.0) -> np.ndarray:

    return amplitude * np.exp(-0.5 * ((x - x0) / sigma) ** 2)

def _voigt_approx(x: np.ndarray, x0: float, gamma: float, sigma: float, amplitude: float = 1.0) -> np.ndarray:

    f_L = gamma / (gamma + sigma) if (gamma + sigma) > 0 else 0.5
    L = _lorentzian(x, x0, gamma, 1.0)
    G = _gaussian(x, x0, sigma, 1.0)
    return amplitude * (f_L * L + (1 - f_L) * G)

@dataclass
class SpinSystem:

    n_spins: int
    shifts_hz: np.ndarray
    couplings_hz: np.ndarray
    equivalences: List[Tuple[int, ...]] = field(default_factory=list)
    labels: List[str] = field(default_factory=list)

    @classmethod
    def from_groups(
        cls,
        groups: List[Dict[str, Any]],
        spectrometer_mhz: float = 400.0
    ) -> 'SpinSystem':

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

        shifts = np.array([p['shift_hz'] for p in all_protons])

        couplings = np.zeros((n, n))
        h_idx_to_spin = {p['h_idx']: i for i, p in enumerate(all_protons)}

        for i, p in enumerate(all_protons):
            for partner in p.get('partners', []):
                partner_h = partner.get('partner_h_index')
                j_hz = partner.get('J_hz', 0.0)
                if partner_h in h_idx_to_spin:
                    j = h_idx_to_spin[partner_h]
                    couplings[i, j] = j_hz
                    couplings[j, i] = j_hz

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

    dim = 2 ** n_spins

    sx = np.array([[0, 1], [1, 0]], dtype=np.complex128) / 2
    sy = np.array([[0, -1j], [1j, 0]], dtype=np.complex128) / 2
    sz = np.array([[1, 0], [0, -1]], dtype=np.complex128) / 2
    I2 = np.eye(2, dtype=np.complex128)

    Ix_list = []
    Iy_list = []
    Iz_list = []

    for k in range(n_spins):

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

    n = spin_system.n_spins
    dim = 2 ** n
    H = np.zeros((dim, dim), dtype=np.complex128)

    for i in range(n):
        H += spin_system.shifts_hz[i] * Iz_list[i]

    for i in range(n):
        for j in range(i + 1, n):
            Jij = spin_system.couplings_hz[i, j]
            if abs(Jij) > 0.01:

                H += Jij * (Iz_list[i] @ Iz_list[j])

    return H

def _compute_transitions(
    H: np.ndarray,
    Ix_total: np.ndarray,
    Iy_total: np.ndarray,
    temperature_k: float = 298.15,
    threshold: float = 1e-6
) -> List[Tuple[float, float]]:

    from scipy import linalg

    eigenvalues, eigenvectors = linalg.eigh(H)
    dim = len(eigenvalues)

    transitions = []

    for i in range(dim):
        for j in range(i + 1, dim):

            freq = eigenvalues[j] - eigenvalues[i]

            if freq > 0:

                matrix_element = eigenvectors[:, j].conj().T @ Ix_total @ eigenvectors[:, i]
                intensity = float(abs(matrix_element) ** 2)

                if intensity > threshold:
                    transitions.append((float(freq), intensity))

    return transitions

def _build_spin_hamiltonian(
    frequencies_hz: List[float],
    J_matrix: np.ndarray,
    Ix_list: List[np.ndarray],
    Iy_list: List[np.ndarray],
    Iz_list: List[np.ndarray]
) -> np.ndarray:

    n_spins = len(frequencies_hz)
    dim = 2 ** n_spins
    H = np.zeros((dim, dim), dtype=complex)

    for i in range(n_spins):
        H -= frequencies_hz[i] * Iz_list[i]

    for i in range(n_spins):
        for j in range(i + 1, n_spins):
            Jij = J_matrix[i, j]
            if abs(Jij) > 1e-10:

                H += Jij * (Ix_list[i] @ Ix_list[j] +
                           Iy_list[i] @ Iy_list[j] +
                           Iz_list[i] @ Iz_list[j])

    return H

def _calculate_nmr_transitions(
    H: np.ndarray,
    Ix_list: List[np.ndarray],
    intensity_threshold: float = 0.001
) -> List[Tuple[float, float, int, int]]:

    from scipy import linalg

    eigenvalues, eigenvectors = linalg.eigh(H)

    Ix_total = sum(Ix_list)

    dim = len(eigenvalues)
    transitions = []

    for i in range(dim):
        for j in range(i + 1, dim):

            freq = eigenvalues[j] - eigenvalues[i]

            matrix_element = eigenvectors[:, j].conj().T @ Ix_total @ eigenvectors[:, i]
            intensity = float(abs(matrix_element) ** 2)

            if intensity > intensity_threshold:
                transitions.append((float(freq), intensity, i, j))

    return transitions

def _simulate_spin_system_exact(
    frequencies_hz: List[float],
    J_matrix: np.ndarray,
    intensity_threshold: float = 0.001,
    use_fragmentation: bool = True,
    fragmentation_threshold: float = 0.05
) -> List[Tuple[float, float]]:

    n_spins = len(frequencies_hz)

    if n_spins == 0:
        return []

    if n_spins > 10:

        if use_fragmentation:

            import sys
            print(f"[INFO] Système de {n_spins} spins > 10, tentative de fragmentation", file=sys.stderr)
            result = _simulate_large_spin_system(
                frequencies_hz, J_matrix,
                threshold_ratio=fragmentation_threshold,
                intensity_threshold=intensity_threshold
            )
            if result is not None:
                return result
            else:
                print(f"[WARN] Fragmentation échouée, fallback approximation premier ordre", file=sys.stderr)

        return None

    Ix_list, Iy_list, Iz_list = _build_spin_operators(n_spins)

    H = _build_spin_hamiltonian(frequencies_hz, J_matrix, Ix_list, Iy_list, Iz_list)

    transitions = _calculate_nmr_transitions(H, Ix_list, intensity_threshold)

    return [(freq, intensity) for freq, intensity, _, _ in transitions]

def _first_order_multiplet(center_hz: float, j_couplings: List[float]) -> Tuple[List[float], List[float]]:

    positions = [center_hz]
    intensities = [1.0]

    for J in j_couplings:
        new_positions = []
        new_intensities = []

        for pos, intensity in zip(positions, intensities):

            new_positions.append(pos - J / 2.0)
            new_intensities.append(intensity / 2.0)
            new_positions.append(pos + J / 2.0)
            new_intensities.append(intensity / 2.0)

        positions = new_positions
        intensities = new_intensities

    return positions, intensities

def _build_coupling_graph(
    frequencies_hz: List[float],
    J_matrix: np.ndarray,
    threshold_ratio: float = 0.05
):

    try:
        import networkx as nx
    except ImportError:

        return None

    n = len(frequencies_hz)
    G = nx.Graph()
    G.add_nodes_from(range(n))

    for i in range(n):
        for j in range(i + 1, n):
            Jij = abs(J_matrix[i, j])
            if Jij < 0.01:
                continue

            delta_nu = abs(frequencies_hz[i] - frequencies_hz[j])

            if delta_nu < 0.001:
                G.add_edge(i, j, weight=Jij)
                continue

            ratio = Jij / delta_nu
            if ratio > threshold_ratio:
                G.add_edge(i, j, weight=Jij)

    return G

def _find_spin_fragments(
    frequencies_hz: List[float],
    J_matrix: np.ndarray,
    threshold_ratio: float = 0.05
) -> List[List[int]]:

    G = _build_coupling_graph(frequencies_hz, J_matrix, threshold_ratio)

    if G is None:

        return [list(range(len(frequencies_hz)))]

    try:
        import networkx as nx
        fragments = list(nx.connected_components(G))
        fragments = [sorted(list(f)) for f in fragments]
        return fragments
    except Exception:

        return [list(range(len(frequencies_hz)))]

def _simulate_fragment(
    fragment_indices: List[int],
    frequencies_hz: List[float],
    J_matrix: np.ndarray,
    intensity_threshold: float = 0.001
) -> List[Tuple[float, float]]:

    n = len(fragment_indices)

    if n == 0:
        return []

    if n > 10:

        import sys
        print(f"[WARN] Fragment de {n} spins trop grand pour simulation exacte", file=sys.stderr)
        return None

    freqs_fragment = [frequencies_hz[i] for i in fragment_indices]

    J_fragment = np.zeros((n, n))
    for local_i, global_i in enumerate(fragment_indices):
        for local_j, global_j in enumerate(fragment_indices):
            J_fragment[local_i, local_j] = J_matrix[global_i, global_j]

    return _simulate_spin_system_exact(freqs_fragment, J_fragment, intensity_threshold)

def _simulate_large_spin_system(
    frequencies_hz: List[float],
    J_matrix: np.ndarray,
    threshold_ratio: float = 0.05,
    intensity_threshold: float = 0.001
) -> List[Tuple[float, float]]:

    n = len(frequencies_hz)

    if n == 0:
        return []

    if n <= 10:
        return _simulate_spin_system_exact(frequencies_hz, J_matrix, intensity_threshold)

    import sys
    print(f"[INFO] Fragmentation: {n} spins, seuil J/Δν = {threshold_ratio}", file=sys.stderr)

    fragments = _find_spin_fragments(frequencies_hz, J_matrix, threshold_ratio)

    print(f"[INFO] → {len(fragments)} fragments détectés", file=sys.stderr)
    for i, frag in enumerate(fragments):
        print(f"[INFO]   Fragment {i+1}: {len(frag)} spins {frag}", file=sys.stderr)

    max_fragment_size = max(len(f) for f in fragments)
    if max_fragment_size > 10:
        print(f"[WARN] Fragment maximum de {max_fragment_size} spins > 10, fragmentation insuffisante", file=sys.stderr)
        print(f"[WARN] → Augmenter threshold_ratio ou utiliser approximation premier ordre", file=sys.stderr)
        return None

    all_transitions = []

    for i, frag in enumerate(fragments):
        transitions = _simulate_fragment(frag, frequencies_hz, J_matrix, intensity_threshold)

        if transitions is None:

            print(f"[ERROR] Fragment {i+1} échec de simulation", file=sys.stderr)
            return None

        if transitions:
            all_transitions.extend(transitions)
            print(f"[INFO]   Fragment {i+1}: {len(transitions)} transitions", file=sys.stderr)

    all_transitions.sort(key=lambda x: x[0])

    print(f"[INFO] Total: {len(all_transitions)} transitions combinées", file=sys.stderr)
    return all_transitions

def simulate_with_fragmentation(
    shifts_hz: np.ndarray,
    j_matrix: np.ndarray,
    freq_mhz: float = 400.0,
    linewidth_hz: float = 1.0,
    ppm_axis: np.ndarray = None,
    line_shape: str = "lorentzian",
    enable_dynamics: bool = False,
    ph: float = 7.4,
) -> 'SimulatedSpectrum':

    n_spins = len(shifts_hz)

    if n_spins == 0:
        return None

    import sys
    print(f"[INFO] SpinSystem HAMILTONIAN: {n_spins} spins, line_shape={line_shape}, linewidth={linewidth_hz} Hz", file=sys.stderr)

    spin_system = SpinSystem(
        shifts_hz=shifts_hz,
        couplings_hz=j_matrix,
        n_spins=n_spins
    )

    if ppm_axis is not None:
        ppm_max = float(np.max(ppm_axis))
        ppm_min = float(np.min(ppm_axis))
        n_points = len(ppm_axis)
    else:

        shifts_ppm = shifts_hz / freq_mhz
        ppm_min = max(-1.0, np.min(shifts_ppm) - 1.0)
        ppm_max = min(15.0, np.max(shifts_ppm) + 1.0)
        n_points = 32768

    try:

        print(f"[INFO] Using Version2 exact Hamiltonian for {n_spins} spins", file=sys.stderr)

        frequencies_hz = list(shifts_hz)

        transitions = _simulate_spin_system_exact(
            frequencies_hz=frequencies_hz,
            J_matrix=j_matrix,
            intensity_threshold=0.001,
            use_fragmentation=True,
            fragmentation_threshold=0.05
        )

        if transitions is None:

            print(f"[WARN] Exact simulation returned None, using first-order approximation", file=sys.stderr)

            all_peaks_hz = []
            all_intensities = []

            for i in range(n_spins):

                couplings_i = []
                for j in range(n_spins):
                    if i != j and abs(j_matrix[i, j]) > 0.1:
                        couplings_i.append(j_matrix[i, j])

                if len(couplings_i) == 0:

                    all_peaks_hz.append(float(shifts_hz[i]))
                    all_intensities.append(1.0)
                else:

                    positions, intensities = _first_order_multiplet(float(shifts_hz[i]), couplings_i)
                    all_peaks_hz.extend(positions)
                    all_intensities.extend(intensities)

            result = simulate_spectrum_lineshape(
                peak_positions_hz=all_peaks_hz,
                peak_intensities=all_intensities,
                spectrometer_mhz=freq_mhz,
                line_width_hz=linewidth_hz,
                ppm_range=(ppm_min, ppm_max),
                n_points=n_points,
                line_shape=line_shape
            )
        else:

            peak_positions_hz = [freq for freq, intensity in transitions]
            peak_intensities = [intensity for freq, intensity in transitions]

            result = simulate_spectrum_lineshape(
                peak_positions_hz=peak_positions_hz,
                peak_intensities=peak_intensities,
                spectrometer_mhz=freq_mhz,
                line_width_hz=linewidth_hz,
                ppm_range=(ppm_min, ppm_max),
                n_points=n_points,
                line_shape=line_shape
            )

        return result
    except Exception as e:
        import sys
        print(f"[ERROR] Hamiltonian simulation failed: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        return None

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

    ppm_min, ppm_max = ppm_range
    ppm_axis = np.linspace(ppm_max, ppm_min, n_points)
    freq_axis_hz = ppm_axis * spectrometer_mhz

    gamma = (line_width_hz + apodization_lb_hz) / 2
    sigma = gamma * 0.8

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

    if spectrum.max() > 0:
        spectrum = spectrum / spectrum.max()

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

    dt = acquisition_time_s / n_points
    t = np.arange(n_points) * dt

    fid = np.zeros(n_points, dtype=np.complex128)

    for freq_hz, amplitude in zip(peak_positions_hz, peak_intensities):

        phase = 2 * np.pi * freq_hz * t
        decay = np.exp(-t / T2_star_s)
        fid += amplitude * np.exp(1j * phase) * decay

    if apodization == "exponential":

        fid *= np.exp(-np.pi * lb_hz * t)
    elif apodization == "gaussian":

        sigma_t = acquisition_time_s / 3
        fid *= np.exp(-(t / sigma_t) ** 2)
    elif apodization == "sine-bell":

        fid *= np.sin(np.pi * t / acquisition_time_s)

    n_fft = n_points * zero_fill_factor
    fid_zf = np.zeros(n_fft, dtype=np.complex128)
    fid_zf[:n_points] = fid

    spectrum_complex = np.fft.fft(fid_zf)
    spectrum_complex = np.fft.fftshift(spectrum_complex)

    spectrum = np.abs(spectrum_complex)

    freq_axis_hz = np.fft.fftshift(np.fft.fftfreq(n_fft, dt))

    ppm_min, ppm_max = ppm_range
    freq_min = ppm_min * spectrometer_mhz
    freq_max = ppm_max * spectrometer_mhz

    mask = (freq_axis_hz >= freq_min) & (freq_axis_hz <= freq_max)
    freq_axis_hz = freq_axis_hz[mask]
    spectrum = spectrum[mask]

    freq_axis_hz = freq_axis_hz[::-1]
    spectrum = spectrum[::-1]

    ppm_axis = freq_axis_hz / spectrometer_mhz

    if spectrum.max() > 0:
        spectrum = spectrum / spectrum.max()

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

def _detect_amino_acid_type(mol: Chem.Mol, ch2_heavy_idx: int) -> Optional[str]:

    try:
        atom = mol.GetAtomWithIdx(ch2_heavy_idx)
        if atom.GetSymbol() != 'C':
            return None

        hyb = str(atom.GetHybridization())
        if 'SP2' in hyb or atom.GetIsAromatic():
            return None

        if atom.IsInRing():
            return None

        has_true_ca_neighbor = False
        for nb in atom.GetNeighbors():
            if nb.GetSymbol() != 'C':
                continue
            if nb.GetIsAromatic() or nb.IsInRing():
                continue
            if 'SP3' not in str(nb.GetHybridization()):
                continue

            has_amine = False
            has_amide_N = False
            has_carboxyl = False
            has_amide_CO = False

            for nn in nb.GetNeighbors():
                if nn.GetSymbol() == 'N' and not nn.GetIsAromatic():

                    is_amide = False
                    for nnn in nn.GetNeighbors():
                        if nnn.GetSymbol() == 'C' and nnn.GetIdx() != nb.GetIdx():
                            for nnnn in nnn.GetNeighbors():
                                if nnnn.GetSymbol() == 'O':
                                    bond = mol.GetBondBetweenAtoms(nnn.GetIdx(), nnnn.GetIdx())
                                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                                        is_amide = True
                                        break
                    if is_amide:
                        has_amide_N = True
                    else:
                        has_amine = True

                elif nn.GetSymbol() == 'C' and nn.GetIdx() != ch2_heavy_idx:

                    has_double_O = False
                    has_single_O = False
                    has_N_on_CO = False

                    for nnn in nn.GetNeighbors():
                        if nnn.GetSymbol() == 'O':
                            bond = mol.GetBondBetweenAtoms(nn.GetIdx(), nnn.GetIdx())
                            if bond:
                                if bond.GetBondType() == Chem.BondType.DOUBLE:
                                    has_double_O = True
                                elif bond.GetBondType() == Chem.BondType.SINGLE:
                                    has_single_O = True
                        elif nnn.GetSymbol() == 'N':
                            has_N_on_CO = True

                    if has_double_O and has_single_O:
                        has_carboxyl = True

                    elif has_double_O and has_N_on_CO:
                        has_amide_CO = True

            is_free_amino_acid = has_amine and has_carboxyl
            is_peptide_context = has_amide_N and has_amide_CO

            if is_free_amino_acid or is_peptide_context:
                has_true_ca_neighbor = True
                break

        if not has_true_ca_neighbor:
            return None

        neighbors = list(atom.GetNeighbors())
        neighbor_symbols = [n.GetSymbol() for n in neighbors]

        for nb in neighbors:
            if nb.GetSymbol() == 'C' and nb.GetIsAromatic():

                ring_info = mol.GetRingInfo()
                for ring in ring_info.AtomRings():
                    if nb.GetIdx() in ring:
                        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                        ring_symbols = [a.GetSymbol() for a in ring_atoms]

                        if len(ring) == 5 and ring_symbols.count('N') == 1:
                            return 'TRP'

                        if len(ring) == 5 and ring_symbols.count('N') == 2:
                            return 'HIS'

                        if len(ring) == 6 and all(a.GetIsAromatic() for a in ring_atoms):
                            for ring_idx in ring:
                                ra = mol.GetAtomWithIdx(ring_idx)
                                for rn in ra.GetNeighbors():
                                    if rn.GetSymbol() == 'O' and not rn.IsInRing():
                                        return 'TYR'
                            return 'PHE'

                return 'PHE'

        if 'O' in neighbor_symbols:
            for nb in neighbors:
                if nb.GetSymbol() == 'O':

                    h_count = sum(1 for n in nb.GetNeighbors() if n.GetSymbol() == 'H')
                    if h_count == 0:
                        h_count = nb.GetTotalNumHs()
                    if h_count >= 1:
                        return 'SER'

        if 'S' in neighbor_symbols:
            for nb in neighbors:
                if nb.GetSymbol() == 'S':
                    return 'CYS'

        for nb in neighbors:
            if nb.GetSymbol() == 'C':

                for nn in nb.GetNeighbors():
                    if nn.GetSymbol() == 'O':
                        bond = mol.GetBondBetweenAtoms(nb.GetIdx(), nn.GetIdx())
                        if bond and bond.GetBondType() == Chem.BondType.DOUBLE:

                            for nnn in nb.GetNeighbors():
                                if nnn.GetSymbol() == 'N':
                                    return 'ASN'
                                if nnn.GetSymbol() == 'O' and nnn.GetIdx() != nn.GetIdx():
                                    return 'ASP'

        def count_hydrogens(atom):

            h_explicit = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == 'H')
            return h_explicit if h_explicit > 0 else atom.GetTotalNumHs()

        for nb in neighbors:
            if nb.GetSymbol() == 'C' and not nb.GetIsAromatic():

                ch3_count = sum(1 for nnn in nb.GetNeighbors()
                              if nnn.GetSymbol() == 'C'
                              and count_hydrogens(nnn) == 3)
                if ch3_count >= 2:
                    return 'LEU'

                chain_length = 0
                current = nb
                visited = {ch2_heavy_idx, nb.GetIdx()}

                while chain_length < 6:
                    chain_length += 1
                    next_c = None
                    for nn in current.GetNeighbors():
                        if nn.GetIdx() in visited:
                            continue
                        if nn.GetSymbol() == 'N':

                            if chain_length >= 3:

                                is_guanidinium = False
                                for n_neighbor in nn.GetNeighbors():
                                    if n_neighbor.GetSymbol() == 'C':

                                        n_count = sum(1 for x in n_neighbor.GetNeighbors()
                                                     if x.GetSymbol() == 'N')
                                        if n_count >= 3:
                                            is_guanidinium = True
                                            break
                                if is_guanidinium:
                                    return 'ARG'
                                else:
                                    return 'LYS'
                        if nn.GetSymbol() == 'S':
                            return 'MET'
                        if nn.GetSymbol() == 'C':

                            ch3_count = sum(1 for nnn in nn.GetNeighbors()
                                          if nnn.GetSymbol() == 'C'
                                          and count_hydrogens(nnn) == 3)
                            if ch3_count >= 2:
                                return 'LEU'
                            next_c = nn
                            visited.add(nn.GetIdx())

                    if next_c is None:
                        break
                    current = next_c

        return None

    except Exception:
        return None

def _is_free_amino_acid(mol: Chem.Mol, ca_idx: int) -> bool:

    try:
        ca = mol.GetAtomWithIdx(ca_idx)

        has_free_n_terminus = False
        has_free_c_terminus = False

        for neighbor in ca.GetNeighbors():

            if neighbor.GetSymbol() == 'N':

                n_h_count = sum(1 for n in neighbor.GetNeighbors() if n.GetSymbol() == 'H')

                n_h_count += neighbor.GetTotalNumHs()

                is_amide_n = False
                for nn in neighbor.GetNeighbors():
                    if nn.GetSymbol() == 'C' and nn.GetIdx() != ca_idx:

                        for nnn in nn.GetNeighbors():
                            if nnn.GetSymbol() == 'O':
                                bond = mol.GetBondBetweenAtoms(nn.GetIdx(), nnn.GetIdx())
                                if bond and bond.GetBondTypeAsDouble() == 2.0:
                                    is_amide_n = True
                                    break

                if n_h_count >= 2 and not is_amide_n:
                    has_free_n_terminus = True

            if neighbor.GetSymbol() == 'C':
                c_idx = neighbor.GetIdx()
                has_carbonyl_o = False
                has_hydroxyl_o = False
                is_amide_c = False

                for cn in neighbor.GetNeighbors():
                    if cn.GetSymbol() == 'O':
                        bond = mol.GetBondBetweenAtoms(c_idx, cn.GetIdx())
                        if bond:
                            if bond.GetBondTypeAsDouble() == 2.0:
                                has_carbonyl_o = True
                            else:
                                has_hydroxyl_o = True
                    elif cn.GetSymbol() == 'N' and cn.GetIdx() != ca_idx:

                        is_amide_c = True

                if has_carbonyl_o and not is_amide_c:
                    has_free_c_terminus = True

        return has_free_n_terminus and has_free_c_terminus

    except Exception:
        return False

def _get_rotamer_averaged_J_for_peptide(
    mol: Chem.Mol,
    h_alpha_idx: int,
    h_beta_idx: int,
    ch2_heavy_idx: int,
    is_hb2: bool = True
) -> Optional[Tuple[float, str]]:

    if not ROTAMER_AVERAGING_AVAILABLE:
        return None

    try:

        aa_type = _detect_amino_acid_type(mol, ch2_heavy_idx)

        if aa_type is None:
            return None

        ca_idx = None
        cb_atom = mol.GetAtomWithIdx(ch2_heavy_idx)
        for neighbor in cb_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C':

                for nn in neighbor.GetNeighbors():
                    if nn.GetSymbol() == 'N':
                        ca_idx = neighbor.GetIdx()
                        break
            if ca_idx is not None:
                break

        if ca_idx is not None and _is_free_amino_acid(mol, ca_idx):
            secondary_structure = 'free_amino_acid'
        else:
            secondary_structure = 'random_coil'

        proton = 'HB2' if is_hb2 else 'HB3'
        J_avg, source = rotamer_averaged_J(
            amino_acid=aa_type,
            coupling_type='HA_HB',
            secondary_structure=secondary_structure,
            proton=proton
        )

        if J_avg > 0 and 'fallback' not in source:
            return (J_avg, source)

        return None

    except Exception:
        return None

def _detect_spin_system(
    mol: Chem.Mol,
    heavy_idx: int,
    x_h_connectivity: Dict[int, List[int]],
    aromatic_dmap: List[Dict[int,int]],
    proton_shifts_ppm: Optional[Dict[int, float]] = None,
    spectrometer_mhz: float = 400.0
) -> Tuple[Optional[str], Dict[str, Any]]:

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

    for nb_idx, n_h_nb in carbon_neighbors:
        delta_nu_J = 100.0

        if proton_shifts_ppm:
            h_this = x_h_connectivity.get(heavy_idx, [])
            h_nb = x_h_connectivity.get(nb_idx, [])
            if h_this and h_nb:
                delta_this = proton_shifts_ppm.get(h_this[0], None)
                delta_nb = proton_shifts_ppm.get(h_nb[0], None)
                if delta_this is not None and delta_nb is not None:
                    delta_ppm = abs(delta_this - delta_nb)
                    delta_nu = delta_ppm * spectrometer_mhz
                    J_typical = 7.0
                    if J_typical > 0:
                        delta_nu_J = delta_nu / J_typical

        info["partners"].append({
            "heavy_idx": nb_idx,
            "n_H": n_h_nb,
            "delta_nu_J": delta_nu_J
        })
        info["delta_nu_J_ratios"].append(delta_nu_J)

    min_ratio = min(info["delta_nu_J_ratios"]) if info["delta_nu_J_ratios"] else 100.0

    is_strongly_coupled = min_ratio < 10.0

    if not is_strongly_coupled:
        return None, info

    total_partners = len(carbon_neighbors)
    partner_h_counts = [p["n_H"] for p in info["partners"]]

    if n_h_this == 2 and total_partners == 1 and partner_h_counts[0] == 2:
        if min_ratio < 5.0:
            return "AA'BB'", info
        else:
            return "A₂X₂", info

    if n_h_this == 1 and total_partners == 1 and partner_h_counts[0] == 2:
        if min_ratio < 8.0:
            return "AB₂", info
        else:
            return "AX₂", info

    if n_h_this == 2 and total_partners == 1 and partner_h_counts[0] == 3:
        if min_ratio < 8.0:
            return "A₂B₃", info
        else:
            return "A₂X₃", info

    if n_h_this == 1 and total_partners >= 2:
        return "ABX", info

    if n_h_this == 1 and total_partners == 1 and partner_h_counts[0] == 1:
        if min_ratio < 5.0:
            return "AB", info
        elif min_ratio < 10.0:
            return "AB_weak", info
        else:
            return "AX", info

    if total_partners >= 2 and all(p["delta_nu_J"] > 10 for p in info["partners"]):
        return "AMX", info

    return "complex_second_order", info

def _detect_AABB_system(
    mol: Chem.Mol,
    heavy_idx: int,
    x_h_connectivity: Dict[int, List[int]],
    aromatic_dmap: List[Dict[int,int]],
    proton_shifts_ppm: Optional[Dict[int, float]] = None,
    spectrometer_mhz: float = 700.0,
    J_typical_hz: float = 7.0
) -> Optional[str]:

    atom = mol.GetAtomWithIdx(heavy_idx)
    if atom.GetSymbol() != "C":
        return None

    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() != "C":
            continue

        nb_idx = neighbor.GetIdx()

        h_list_this = x_h_connectivity.get(heavy_idx, [])
        h_list_neighbor = x_h_connectivity.get(nb_idx, [])

        if len(h_list_this) != 2:
            continue
        if len(h_list_neighbor) != 2:
            continue

        bond = mol.GetBondBetweenAtoms(heavy_idx, nb_idx)
        if not bond or bond.GetBondTypeAsDouble() != 1.0:
            continue

        if proton_shifts_ppm is not None:

            shifts_this = [proton_shifts_ppm.get(h, None) for h in h_list_this]
            shifts_neighbor = [proton_shifts_ppm.get(h, None) for h in h_list_neighbor]

            if all(s is not None for s in shifts_this) and all(s is not None for s in shifts_neighbor):
                avg_shift_this = sum(shifts_this) / len(shifts_this)
                avg_shift_neighbor = sum(shifts_neighbor) / len(shifts_neighbor)

                delta_ppm = abs(avg_shift_this - avg_shift_neighbor)
                delta_nu_hz = delta_ppm * spectrometer_mhz

                ratio = delta_nu_hz / J_typical_hz if J_typical_hz > 0 else float('inf')

                if ratio > 10.0:

                    return None

        return "AA'BB'"

    return None

def _check_tables_coverage(mol: Chem.Mol) -> Dict[str, Any]:

    covered_types = []
    uncovered_types = []

    n_alkene = 0
    n_aromatic = 0
    n_saturated_ring = 0
    n_ester_lactone = 0
    n_alcohol_ether = 0
    n_flexible_sp3_sp3 = 0

    ring_info = mol.GetRingInfo()

    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()

            if a1.GetSymbol() == 'C' and a2.GetSymbol() == 'C':
                n_alkene += 1
            elif (a1.GetSymbol() == 'C' and a2.GetSymbol() == 'O') or \
                 (a1.GetSymbol() == 'O' and a2.GetSymbol() == 'C'):
                n_ester_lactone += 1

    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            n_aromatic += 1

    for ring in ring_info.AtomRings():
        is_aromatic = any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        if not is_aromatic:
            n_saturated_ring += 1

    ester_pattern = Chem.MolFromSmarts('[OX2][CX3](=[OX1])')
    if mol.HasSubstructMatch(ester_pattern):
        n_ester_lactone += len(mol.GetSubstructMatches(ester_pattern))

    alcohol_pattern = Chem.MolFromSmarts('[CX4][OX2H]')
    if mol.HasSubstructMatch(alcohol_pattern):
        n_alcohol_ether += len(mol.GetSubstructMatches(alcohol_pattern))

    ether_pattern = Chem.MolFromSmarts('[CX4][OX2][CX4]')
    if mol.HasSubstructMatch(ether_pattern):
        n_alcohol_ether += len(mol.GetSubstructMatches(ether_pattern))

    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        if bond.IsInRing():
            continue

        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()

        if a1.GetSymbol() in ('O', 'N') or a2.GetSymbol() in ('O', 'N'):
            continue

        if a1.GetSymbol() == 'C' and a2.GetSymbol() == 'C':
            hyb1 = a1.GetHybridization()
            hyb2 = a2.GetHybridization()

            if hyb1 == Chem.HybridizationType.SP3 and hyb2 == Chem.HybridizationType.SP2:
                continue
            if hyb1 == Chem.HybridizationType.SP2 and hyb2 == Chem.HybridizationType.SP3:
                continue

            if hyb1 == Chem.HybridizationType.SP3 and hyb2 == Chem.HybridizationType.SP3:

                n_h1 = a1.GetTotalNumHs()
                n_h2 = a2.GetTotalNumHs()

                if n_h1 == 3 or n_h2 == 3:
                    continue

                n_flexible_sp3_sp3 += 1

    if n_alkene > 0:
        covered_types.append(f"alkènes ({n_alkene})")
    if n_aromatic > 0:
        covered_types.append(f"aromatiques ({n_aromatic} atomes)")
    if n_saturated_ring > 0:
        covered_types.append(f"cycles saturés ({n_saturated_ring})")
    if n_ester_lactone > 0:
        covered_types.append(f"esters/lactones ({n_ester_lactone})")
    if n_alcohol_ether > 0:
        covered_types.append(f"alcools/éthers ({n_alcohol_ether})")

    if n_flexible_sp3_sp3 > 0:

        covered_types.append(f"chaînes sp3-sp3 flexibles ({n_flexible_sp3_sp3}) → J=7.0 Hz fallback")

    all_covered = True

    return {
        'all_covered': all_covered,
        'covered_types': covered_types,
        'uncovered_types': uncovered_types,
        'n_flexible_sp3_sp3': n_flexible_sp3_sp3
    }

def _molecule_needs_conformers(mol: Chem.Mol) -> Tuple[bool, str]:

    from rdkit.Chem import Descriptors, Lipinski

    mol_h = Chem.AddHs(mol)
    ring_info = mol.GetRingInfo()

    tables_coverage = _check_tables_coverage(mol)
    if tables_coverage['all_covered']:
        return (False, f"Tables couvrent tous les J-couplings: {tables_coverage['covered_types']}")

    for ring in ring_info.AtomRings():
        ring_size = len(ring)
        if ring_size >= 8:

            n_ch2_in_ring = 0
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 2:
                    n_ch2_in_ring += 1

            return (False, f"Macrocycle ({ring_size} membres, {n_ch2_in_ring} CH2) - table J moyens utilisée (pas de conformères)")

    n_saturated_rings = 0
    n_semi_rigid_rings = 0
    n_flexible_sp3_sp3 = 0
    n_allylic_sp3_sp2 = 0
    n_diastereotopic_ch2 = 0

    reasons_need = []
    reasons_skip = []

    for ring in ring_info.AtomRings():
        ring_size = len(ring)
        if ring_size < 5 or ring_size > 7:
            continue

        is_saturated = True
        is_rigid = False
        has_heteroatom = False
        has_sp3_carbon = False

        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetIsAromatic():
                is_saturated = False
                is_rigid = True
                break
            if atom.GetSymbol() in ('N', 'O', 'S'):
                has_heteroatom = True
            if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP3:
                has_sp3_carbon = True

        if is_saturated:
            for i, atom_idx in enumerate(ring):
                for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(atom_idx, neighbor.GetIdx())
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:

                        is_rigid = True
                        break
                if is_rigid:
                    break

        if is_saturated and not is_rigid and has_sp3_carbon:
            n_saturated_rings += 1
            if ring_size in (5, 6, 7):
                n_semi_rigid_rings += 1
                ring_type = "pipéridine/tétrahydropyrane" if has_heteroatom else "cyclohexane"
                reasons_need.append(f"Cycle saturé flexible {ring_size}-membres ({ring_type})")

    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        if bond.IsInRing():
            continue

        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()

        if a1.GetSymbol() != 'C' or a2.GetSymbol() != 'C':
            continue

        hyb1 = a1.GetHybridization()
        hyb2 = a2.GetHybridization()

        if (hyb1 == Chem.HybridizationType.SP3 and hyb2 == Chem.HybridizationType.SP2) or \
           (hyb1 == Chem.HybridizationType.SP2 and hyb2 == Chem.HybridizationType.SP3):
            n_allylic_sp3_sp2 += 1

        elif hyb1 == Chem.HybridizationType.SP3 and hyb2 == Chem.HybridizationType.SP3:

            n_flexible_sp3_sp3 += 1

            has_nearby_chiral = False
            for neighbor in a1.GetNeighbors():
                if neighbor.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                    has_nearby_chiral = True
                    break
            for neighbor in a2.GetNeighbors():
                if neighbor.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                    has_nearby_chiral = True
                    break

            if has_nearby_chiral:
                n_diastereotopic_ch2 += 1

    n_stereocenters = sum(1 for a in mol.GetAtoms()
                         if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED)
    n_fused_rings = len([r for r in ring_info.AtomRings() if len(r) >= 5])

    is_steroid_like = (
        (n_fused_rings >= 4 and n_stereocenters >= 5) or
        (n_fused_rings >= 3 and n_stereocenters >= 5) or
        (n_fused_rings >= 2 and n_stereocenters >= 5)
    )

    if is_steroid_like:

        n_rotatable = Lipinski.NumRotatableBonds(mol)
        flexible_fraction = n_rotatable / max(1, mol.GetNumHeavyAtoms())

        if flexible_fraction < 0.20:
            return (False, f"Système stéroïde/polycyclique verrouillé ({n_fused_rings} cycles, {n_stereocenters} stéréocentres, {n_rotatable} rotatives) - J fixes par stéréochimie")
        else:

            return (False, f"Stéroïde avec chaîne latérale ({n_fused_rings} cycles, {n_stereocenters} stéréocentres) - squelette rigide, 1-10 conformères suffisants")

    has_fused_lactam = False
    fused_lactam_info = ""

    lactam_patterns = [

        ("[NX3]1[CX3](=O)[CX4,NX3][NX3][CX3](=O)[CX4,NX3]1", "DKP (dicétopipérazine)"),

        ("[NX3]1[CX3](=O)[CX4][CX4][CX4]1", "pyrrolidinone"),

        ("[NX3]1[CX3](=O)[CX4][CX4]1", "β-lactame"),

        ("[OX2]1[CX3](=O)[NX3][CX4][CX4]1", "oxazolidinone"),

        ("[NX3]1[CX3](=O)[NX3][CX4][CX4]1", "imidazolidinone"),

        ("[NX3]1[CX3](=O)@[CX4]@[CX4]@[CX4,NX3,OX2]1", "lactame-5 fusionné"),

        ("[NX3]1[CX3](=O)@[CX4]@[CX4]@[CX4]@[CX4,NX3,OX2]1", "lactame-6 fusionné"),
    ]

    for smarts, name in lactam_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):

            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                match_set = set(match)
                for ring in ring_info.AtomRings():
                    ring_set = set(ring)

                    shared = match_set & ring_set
                    if len(shared) >= 2 and ring_set != match_set:
                        has_fused_lactam = True
                        fused_lactam_info = name
                        break
                if has_fused_lactam:
                    break
        if has_fused_lactam:
            break

    if has_fused_lactam:
        return (False, f"Lactame fusionné détecté ({fused_lactam_info}) → système rigide, J fixes par tables")

    if n_semi_rigid_rings > 0 and not is_steroid_like:

        return (False, f"Cycles saturés ({n_semi_rigid_rings}x) → J tabulés par taille de cycle (ax/eq ou moyenne)")

    n_rotatable = Lipinski.NumRotatableBonds(mol)
    n_heavy = mol.GetNumHeavyAtoms()
    rotatable_ratio = n_rotatable / max(1, n_heavy)

    is_highly_flexible = (n_rotatable >= 6 and rotatable_ratio > 0.25)

    has_mostly_allylic = (n_allylic_sp3_sp2 >= 2 and n_allylic_sp3_sp2 >= n_flexible_sp3_sp3 * 0.5)

    if n_flexible_sp3_sp3 > 0:
        return (False, f"Chaînes C(sp3)-C(sp3) flexibles ({n_flexible_sp3_sp3}x) → J moyen 7.0 Hz (MMFF biaisé, multi-conformer désactivé)")

    if n_allylic_sp3_sp2 > 0 and n_flexible_sp3_sp3 <= 2:
        return (False, f"Liaisons C(sp3)-C(sp2) à rotation libre ({n_allylic_sp3_sp2}x) → J=7.0 Hz fixe")

    n_aromatic = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
    n_total = mol.GetNumHeavyAtoms()
    if n_aromatic > n_total * 0.6:
        return (False, f"Molécule principalement aromatique ({n_aromatic}/{n_total} atomes)")

    if n_total < 15 and n_saturated_rings == 0:
        return (False, f"Petite molécule sans cycles saturés ({n_total} atomes)")

    if n_flexible_sp3_sp3 > 0 and n_diastereotopic_ch2 == 0:
        return (False, f"Chaînes flexibles sans stéréocentres → J moyen 7.0 Hz")

    return (False, "Aucune structure nécessitant un échantillonnage conformationnel détectée")

def _has_amide_in_ring(mol: Chem.Mol, ring_atoms: Tuple[int, ...]) -> bool:

    ring_set = set(ring_atoms)

    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)

        if atom.GetSymbol() == 'N':

            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() in ring_set:

                    for n2 in neighbor.GetNeighbors():
                        if n2.GetSymbol() == 'O':
                            bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), n2.GetIdx())
                            if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                                return True

    return False

def _boltzmann_average_J(
    J_values_per_conformer: List[float],
    energies_kcal: List[float],
    temperature_k: float = 298.15
) -> float:

    if not J_values_per_conformer or not energies_kcal:
        return 0.0

    if len(J_values_per_conformer) != len(energies_kcal):
        return float(np.mean(J_values_per_conformer))

    R = 0.001987
    RT = R * temperature_k

    min_E = min(energies_kcal)
    rel_energies = [e - min_E for e in energies_kcal]

    boltzmann_factors = [math.exp(-E / RT) for E in rel_energies]
    Z = sum(boltzmann_factors)

    if Z <= 0:
        return float(np.mean(J_values_per_conformer))

    populations = [bf / Z for bf in boltzmann_factors]

    J_avg = sum(p * J for p, J in zip(populations, J_values_per_conformer))

    return J_avg

class ConformerCache:

    def __init__(
        self,
        mol: Chem.Mol,
        n_conformers: int = 50,
        energy_cutoff_kcal: float = 5.0,
        random_seed: int = 42,
        force_field: str = "MMFF94s",
        pruning_rmsd: float = 0.35,
        use_etkdg_v3: bool = True,
        use_systematic_sampling: bool = True
    ):
        self.mol_h = Chem.AddHs(mol, addCoords=True)
        self.n_conformers = n_conformers
        self.energy_cutoff = energy_cutoff_kcal
        self.valid_conf_ids: List[int] = []
        self.energies: List[float] = []
        self._initialized = False
        self._random_seed = random_seed
        self._force_field = force_field
        self._pruning_rmsd = pruning_rmsd
        self._use_etkdg_v3 = use_etkdg_v3
        self._use_systematic_sampling = use_systematic_sampling

    def initialize(self) -> bool:

        if self._initialized:
            return len(self.valid_conf_ids) > 0

        try:

            if self._use_systematic_sampling:
                systematic_results = _generate_diverse_conformers_systematic(
                    Chem.RemoveHs(self.mol_h),
                    n_conformers=self.n_conformers,
                    energy_cutoff_kcal=self.energy_cutoff
                )

                if systematic_results:

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

            if self._use_etkdg_v3:
                try:
                    ps = AllChem.ETKDGv3()
                    ps.useSmallRingTorsions = True
                    ps.useMacrocycleTorsions = True
                except AttributeError:
                    ps = AllChem.ETKDG()
            else:
                ps = AllChem.ETKDG()

            ps.randomSeed = self._random_seed
            ps.numThreads = 0
            ps.pruneRmsThresh = self._pruning_rmsd
            ps.useRandomCoords = True

            n_generate = self.n_conformers * 5

            conf_ids = AllChem.EmbedMultipleConfs(
                self.mol_h,
                numConfs=n_generate,
                params=ps
            )

            if len(conf_ids) == 0:

                ps.useRandomCoords = True
                res = AllChem.EmbedMolecule(self.mol_h, ps)
                if res == 0:
                    _optimize_molecule_geometry(self.mol_h, self._force_field, 200, 0)
                    self.valid_conf_ids = [0]
                    self.energies = [0.0]
                self._initialized = True
                return len(self.valid_conf_ids) > 0

            ff_variant = self._force_field.upper()
            if "MMFF" in ff_variant:
                variant = "MMFF94s" if "S" in ff_variant else "MMFF94"
                try:
                    results = AllChem.MMFFOptimizeMoleculeConfs(self.mol_h, mmffVariant=variant, maxIters=200)
                except Exception:

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

                results = []
                for cid in conf_ids:
                    try:
                        AllChem.UFFOptimizeMolecule(self.mol_h, confId=cid, maxIters=200)
                        ff = AllChem.UFFGetMoleculeForceField(self.mol_h, confId=cid)
                        energy = ff.CalcEnergy() if ff else 0.0
                        results.append((0, energy))
                    except Exception:
                        results.append((1, 0.0))

            energy_list = []
            for conf_id, (converged, energy) in zip(conf_ids, results):
                if converged == 0:
                    energy_list.append((conf_id, energy))

            if not energy_list:
                self._initialized = True
                return False

            energy_list.sort(key=lambda x: x[1])
            min_energy = energy_list[0][1]

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

        if not self._initialized:
            self.initialize()

        if not self.valid_conf_ids:
            return None

        try:

            h1_atom = self.mol_h.GetAtomWithIdx(h1_idx)
            h2_atom = self.mol_h.GetAtomWithIdx(h2_idx)

            c1_neighbors = [n for n in h1_atom.GetNeighbors() if n.GetSymbol() == 'C']
            c2_neighbors = [n for n in h2_atom.GetNeighbors() if n.GetSymbol() == 'C']

            if not c1_neighbors or not c2_neighbors:
                return None

            c1_idx = c1_neighbors[0].GetIdx()
            c2_idx = c2_neighbors[0].GetIdx()

            if not self.mol_h.GetBondBetweenAtoms(c1_idx, c2_idx):
                return None

            subs_alpha = []
            subs_beta = []

            if use_haasnoot_7_params:
                c1_atom = self.mol_h.GetAtomWithIdx(c1_idx)
                c2_atom = self.mol_h.GetAtomWithIdx(c2_idx)

                for nb in c1_atom.GetNeighbors():
                    if nb.GetIdx() != c2_idx and nb.GetSymbol() != 'H':
                        subs_alpha.append(nb.GetSymbol())

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

            return _boltzmann_average_J(J_values, valid_energies, temperature_k)

        except Exception:
            return None

    def get_num_conformers(self) -> int:

        if not self._initialized:
            self.initialize()
        return len(self.valid_conf_ids)

    def get_populations(self, temperature_k: float = 298.15) -> List[float]:

        if not self._initialized:
            self.initialize()
        if not self.energies:
            return []

        R = 0.001987
        RT = R * temperature_k

        boltz = [math.exp(-E / RT) for E in self.energies]
        Z = sum(boltz)

        if Z <= 0:
            return [1.0 / len(self.energies)] * len(self.energies)

        return [b / Z for b in boltz]

def _count_rotatable_bonds(mol: Chem.Mol) -> int:

    try:
        return Descriptors.NumRotatableBonds(mol)
    except Exception:

        count = 0
        for bond in mol.GetBonds():
            if bond.IsInRing():
                continue
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
            if a1.GetSymbol() == "H" or a2.GetSymbol() == "H":
                continue

            if a1.GetDegree() <= 1 or a2.GetDegree() <= 1:
                continue
            count += 1
        return count

def _estimate_molecular_flexibility(mol: Chem.Mol) -> Tuple[str, int]:

    n_rot = _count_rotatable_bonds(mol)

    n_rings = mol.GetRingInfo().NumRings()

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

    flexibility, n_rot = _estimate_molecular_flexibility(mol)

    import math

    if n_rot == 0:

        return min_n

    scale_factor = 1.0 + 0.5 * math.log(1 + n_rot)

    if flexibility == "rigid":
        n_conf = int(base_n * 0.3 * scale_factor)
    elif flexibility == "moderate":
        n_conf = int(base_n * 0.6 * scale_factor)
    elif flexibility == "flexible":
        n_conf = int(base_n * 1.0 * scale_factor)
    else:
        n_conf = int(base_n * 1.5 * scale_factor)

    if _has_long_flexible_chain(mol, min_length=6):
        n_conf = int(n_conf * 1.5)

    n_conf = max(min_n, min(max_n, n_conf))

    return n_conf

def _has_long_flexible_chain(mol: Chem.Mol, min_length: int = 4) -> bool:

    for atom in mol.GetAtoms():
        if atom.IsInRing() or atom.GetSymbol() != "C":
            continue
        if atom.GetDegree() <= 1:
            continue

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

PROTON_SOLVENT_SENSITIVITY: Dict[str, Dict[str, float]] = {

    "OH_phenol": {"a": -0.5, "b": +2.0, "c": +0.3, "d": 0.0},
    "OH_alcohol": {"a": -0.3, "b": +1.5, "c": +0.2, "d": 0.0},
    "OH_carboxylic": {"a": -0.8, "b": +2.5, "c": +0.5, "d": 0.0},
    "NH_amide": {"a": -0.2, "b": +1.2, "c": +0.3, "d": 0.0},
    "NH_amine": {"a": -0.4, "b": +1.0, "c": +0.2, "d": 0.0},
    "NH_aromatic": {"a": -0.3, "b": +1.5, "c": +0.4, "d": 0.0},

    "CHO_aldehyde": {"a": +0.3, "b": -0.2, "c": +0.1, "d": 0.0},
    "CH_ketone_alpha": {"a": +0.2, "b": -0.1, "c": +0.05, "d": 0.0},

    "Ar_H": {"a": +0.05, "b": -0.1, "c": +0.15, "d": 0.0},
    "Ar_H_electron_rich": {"a": +0.02, "b": -0.15, "c": +0.1, "d": 0.0},
    "Ar_H_electron_poor": {"a": +0.08, "b": -0.05, "c": +0.2, "d": 0.0},

    "CH3_alkyl": {"a": 0.0, "b": 0.0, "c": +0.02, "d": 0.0},
    "CH2_alkyl": {"a": 0.0, "b": 0.0, "c": +0.02, "d": 0.0},
    "CH_alkyl": {"a": 0.0, "b": 0.0, "c": +0.02, "d": 0.0},

    "CH_vinyl": {"a": +0.05, "b": -0.05, "c": +0.08, "d": 0.0},
    "CH2_vinyl": {"a": +0.03, "b": -0.03, "c": +0.05, "d": 0.0},
}

TEMPERATURE_COEFFICIENTS: Dict[str, float] = {

    "OH_phenol": -0.0080,
    "OH_alcohol": -0.0100,
    "OH_carboxylic": -0.0120,
    "NH_amide": -0.0060,
    "NH_amine": -0.0080,
    "NH_aromatic": -0.0070,

    "CHO_aldehyde": -0.0015,
    "Ar_H": -0.0008,
    "CH_vinyl": -0.0005,
    "CH3_alkyl": -0.0003,
    "CH2_alkyl": -0.0003,
    "CH_alkyl": -0.0003,

    "default": -0.0005,
}

PROTON_TYPE_SMARTS: Dict[str, str] = {
    "OH_phenol": "[OH;$([OH]c)]",
    "OH_alcohol": "[OH;$([OH]C);!$([OH]c)]",
    "OH_carboxylic": "[OH;$([OH]C=O)]",
    "NH_amide": "[NH;$([NH]C=O)]",
    "NH_amine": "[NH2,NH;!$([NH]C=O);!$([NH]c)]",
    "NH_aromatic": "[nH]",
    "CHO_aldehyde": "[CH;$([CH](=O))]",
    "CH_ketone_alpha": "[CH;$([CH]C(=O))]",
    "Ar_H": "[cH]",
    "CH_vinyl": "[CH;$([CH]=C)]",
    "CH2_vinyl": "[CH2;$([CH2]=C)]",
    "CH3_alkyl": "[CH3;$([CH3]C)]",
    "CH2_alkyl": "[CH2;$([CH2](C)C)]",
    "CH_alkyl": "[CH;$([CH](C)(C)C)]",
}

SOLVENT_SHIFT_CORRECTIONS: Dict[str, Dict[str, float]] = {
    "CDCl3": {
        "aromatic": 0.0,
        "alkyl": 0.0,
        "OH": -2.0,
        "NH": -1.5,
        "aldehyde": 0.0,
    },
    "DMSO-d6": {
        "aromatic": +0.1,
        "alkyl": 0.0,
        "OH": +0.5,
        "NH": +0.3,
        "aldehyde": +0.1,
    },
    "D2O": {
        "aromatic": -0.1,
        "alkyl": -0.05,
        "OH": 0.0,
        "NH": 0.0,
        "aldehyde": -0.1,
    },
    "CD3OD": {
        "aromatic": 0.0,
        "alkyl": 0.0,
        "OH": 0.0,
        "NH": 0.0,
        "aldehyde": 0.0,
    },
    "C6D6": {
        "aromatic": -0.5,
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

SOLVENT_J_CORRECTIONS: Dict[str, Dict[str, float]] = {
    "CDCl3": {

        "3J_vicinal": 1.00,
        "3J_vinyl": 1.00,
        "2J_geminal": 1.00,
        "4J_allylic": 1.00,
        "3J_aromatic": 1.00,
        "3J_HNCH": 1.00,
    },
    "DMSO-d6": {

        "3J_vicinal": 1.02,
        "3J_vinyl": 0.98,
        "2J_geminal": 1.00,
        "4J_allylic": 1.05,
        "3J_aromatic": 1.00,
        "3J_HNCH": 1.08,
    },
    "D2O": {

        "3J_vicinal": 0.98,
        "3J_vinyl": 0.97,
        "2J_geminal": 1.00,
        "4J_allylic": 0.90,
        "3J_aromatic": 1.00,
        "3J_HNCH": 0.0,
    },
    "CD3OD": {

        "3J_vicinal": 0.99,
        "3J_vinyl": 0.98,
        "2J_geminal": 1.00,
        "4J_allylic": 0.95,
        "3J_aromatic": 1.00,
        "3J_HNCH": 0.5,
    },
    "C6D6": {

        "3J_vicinal": 1.00,
        "3J_vinyl": 1.02,
        "2J_geminal": 1.00,
        "4J_allylic": 1.05,
        "3J_aromatic": 0.98,
        "3J_HNCH": 1.00,
    },
    "CD3CN": {

        "3J_vicinal": 1.01,
        "3J_vinyl": 1.00,
        "2J_geminal": 1.00,
        "4J_allylic": 1.02,
        "3J_aromatic": 1.00,
        "3J_HNCH": 1.02,
    },
    "pyridine-d5": {

        "3J_vicinal": 1.02,
        "3J_vinyl": 1.00,
        "2J_geminal": 1.00,
        "4J_allylic": 1.05,
        "3J_aromatic": 0.98,
        "3J_HNCH": 1.10,
    },
}

def _get_solvent_J_correction(
    solvent_name: str,
    j_type: str = "3J_vicinal"
) -> float:

    corrections = SOLVENT_J_CORRECTIONS.get(solvent_name, SOLVENT_J_CORRECTIONS.get("CDCl3", {}))
    return corrections.get(j_type, 1.0)

def _apply_solvent_J_correction(
    J_hz: float,
    solvent_name: str,
    j_type: str = "3J_vicinal",
    is_labile_H: bool = False
) -> float:

    if is_labile_H and solvent_name == "D2O":

        return 0.0

    factor = _get_solvent_J_correction(solvent_name, j_type)
    return J_hz * factor

def _apply_solvent_corrections(
    proton_shifts_ppm: Dict[int, float],
    mol: Chem.Mol,
    heavy_to_H: Dict[int, List[int]],
    solvent_name: str
) -> Dict[int, float]:

    corrections = SOLVENT_SHIFT_CORRECTIONS.get(solvent_name, {})
    if not corrections:
        return proton_shifts_ppm

    corrected = dict(proton_shifts_ppm)

    for heavy_idx, h_list in heavy_to_H.items():
        if not h_list:
            continue

        atom = mol.GetAtomWithIdx(heavy_idx)
        sym = atom.GetSymbol()

        proton_type = "alkyl"

        if sym == "O":
            proton_type = "OH"
        elif sym == "N":
            proton_type = "NH"
        elif sym == "C":
            if atom.GetIsAromatic():
                proton_type = "aromatic"
            elif any(b.GetBondTypeAsDouble() >= 2 for b in atom.GetBonds()):

                for b in atom.GetBonds():
                    other = b.GetOtherAtom(atom)
                    if other.GetSymbol() == "O" and b.GetBondTypeAsDouble() >= 2:
                        proton_type = "aldehyde"
                        break

        correction = corrections.get(proton_type, 0.0)
        for h_idx in h_list:
            if h_idx in corrected:
                corrected[h_idx] = corrected[h_idx] + correction

    return corrected

def _get_solvent_linewidth(solvent_name: str, default: float = 0.8) -> float:

    return SOLVENT_LINEWIDTH.get(solvent_name, default)

def _identify_proton_type_smarts(mol: Chem.Mol, h_idx: int) -> str:

    h_atom = mol.GetAtomWithIdx(h_idx)

    neighbors = h_atom.GetNeighbors()
    if not neighbors:
        return "default"

    heavy_atom = neighbors[0]
    heavy_sym = heavy_atom.GetSymbol()
    heavy_idx = heavy_atom.GetIdx()

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

    params_target = KAMLET_TAFT_PARAMS.get(solvent_name, KAMLET_TAFT_PARAMS.get("CDCl3"))
    params_ref = KAMLET_TAFT_PARAMS.get(reference_solvent, KAMLET_TAFT_PARAMS.get("CDCl3"))

    d_alpha = params_target["alpha"] - params_ref["alpha"]
    d_beta = params_target["beta"] - params_ref["beta"]
    d_pi = params_target["pi_star"] - params_ref["pi_star"]

    sensitivity = PROTON_SOLVENT_SENSITIVITY.get(proton_type,
                  PROTON_SOLVENT_SENSITIVITY.get("CH3_alkyl", {"a": 0, "b": 0, "c": 0, "d": 0}))

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

    coef = TEMPERATURE_COEFFICIENTS.get(proton_type,
                                        TEMPERATURE_COEFFICIENTS.get("default", -0.0005))

    dT = temperature_k - reference_temp_k

    return coef * dT

def _is_proton_exchangeable(mol: Chem.Mol, h_idx: int) -> bool:

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

    if not _is_proton_exchangeable(mol, h_idx):
        return {"visible": True, "broadening": 1.0, "shift_correction": 0.0}

    h_atom = mol.GetAtomWithIdx(h_idx)
    heavy_atom = h_atom.GetNeighbors()[0]
    heavy_sym = heavy_atom.GetSymbol()

    exchange_solvents = {
        "D2O": {"OH": False, "NH": False, "SH": True},
        "CD3OD": {"OH": False, "NH": True, "SH": True},
        "CD3OD + D2O": {"OH": False, "NH": False, "SH": True},
    }

    visibility_map = exchange_solvents.get(solvent_name, {"OH": True, "NH": True, "SH": True})
    proton_class = "OH" if heavy_sym == "O" else ("NH" if heavy_sym == "N" else "SH")
    visible = visibility_map.get(proton_class, True)

    broadening_factors = {
        "CDCl3": {"OH": 3.0, "NH": 2.0, "SH": 1.5},
        "DMSO-d6": {"OH": 1.2, "NH": 1.1, "SH": 1.0},
        "C6D6": {"OH": 2.0, "NH": 1.5, "SH": 1.2},
        "acetone-d6": {"OH": 2.5, "NH": 1.8, "SH": 1.3},
        "pyridine-d5": {"OH": 1.3, "NH": 1.1, "SH": 1.0},
    }

    broadening_map = broadening_factors.get(solvent_name, {"OH": 2.0, "NH": 1.5, "SH": 1.2})
    broadening = broadening_map.get(proton_class, 1.0)

    T_factor = (temperature_k / 298.15) ** 1.5
    broadening *= T_factor

    kt_params = KAMLET_TAFT_PARAMS.get(solvent_name, KAMLET_TAFT_PARAMS.get("CDCl3"))
    shift_correction = 0.0

    if proton_class in ("OH", "NH"):

        shift_correction = kt_params["beta"] * 1.5

    return {
        "visible": visible,
        "broadening": broadening,
        "shift_correction": shift_correction
    }

def _build_spin_operators(n_spins: int) -> Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray]]:

    dim = 2 ** n_spins

    Ix_single = np.array([[0, 0.5], [0.5, 0]], dtype=complex)
    Iy_single = np.array([[0, -0.5j], [0.5j, 0]], dtype=complex)
    Iz_single = np.array([[0.5, 0], [0, -0.5]], dtype=complex)
    I_identity = np.eye(2, dtype=complex)

    Ix_list = []
    Iy_list = []
    Iz_list = []

    for spin_idx in range(n_spins):

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

    n_spins = len(frequencies_hz)
    dim = 2 ** n_spins
    H = np.zeros((dim, dim), dtype=complex)

    for i in range(n_spins):
        H -= frequencies_hz[i] * Iz_list[i]

    for i in range(n_spins):
        for j in range(i + 1, n_spins):
            Jij = J_matrix[i, j]
            if abs(Jij) > 1e-10:

                H += Jij * (Ix_list[i] @ Ix_list[j] +
                           Iy_list[i] @ Iy_list[j] +
                           Iz_list[i] @ Iz_list[j])

    return H

def _calculate_nmr_transitions(
    H: np.ndarray,
    Ix_list: List[np.ndarray],
    intensity_threshold: float = 0.001
) -> List[Tuple[float, float, int, int]]:

    from scipy import linalg

    eigenvalues, eigenvectors = linalg.eigh(H)

    Ix_total = sum(Ix_list)

    dim = len(eigenvalues)
    transitions = []

    for i in range(dim):
        for j in range(i + 1, dim):

            freq = eigenvalues[j] - eigenvalues[i]

            matrix_element = eigenvectors[:, j].conj().T @ Ix_total @ eigenvectors[:, i]
            intensity = float(abs(matrix_element) ** 2)

            if intensity > intensity_threshold:
                transitions.append((float(freq), intensity, i, j))

    return transitions

def _simulate_spin_system_exact(
    frequencies_hz: List[float],
    J_matrix: np.ndarray,
    intensity_threshold: float = 0.001,
    use_fragmentation: bool = True,
    fragmentation_threshold: float = 0.05
) -> List[Tuple[float, float]]:

    n_spins = len(frequencies_hz)

    if n_spins == 0:
        return []

    if n_spins > 10:

        if use_fragmentation:

            import sys
            print(f"[INFO] Système de {n_spins} spins > 10, tentative de fragmentation", file=sys.stderr)
            result = _simulate_large_spin_system(
                frequencies_hz, J_matrix,
                threshold_ratio=fragmentation_threshold,
                intensity_threshold=intensity_threshold
            )
            if result is not None:
                return result
            else:
                print(f"[WARN] Fragmentation échouée, fallback approximation premier ordre", file=sys.stderr)

        return None

    Ix_list, Iy_list, Iz_list = _build_spin_operators(n_spins)

    H = _build_spin_hamiltonian(frequencies_hz, J_matrix, Ix_list, Iy_list, Iz_list)

    transitions = _calculate_nmr_transitions(H, Ix_list, intensity_threshold)

    return [(freq, intensity) for freq, intensity, _, _ in transitions]

def simulate_AB_system(
    delta_A_ppm: float,
    delta_B_ppm: float,
    J_AB_hz: float,
    spectrometer_mhz: float = 400.0
) -> Dict[str, Any]:

    nu_A = delta_A_ppm * spectrometer_mhz
    nu_B = delta_B_ppm * spectrometer_mhz

    delta_nu = abs(nu_A - nu_B)
    ratio = delta_nu / abs(J_AB_hz) if J_AB_hz != 0 else float('inf')

    J_matrix = np.array([[0, J_AB_hz], [J_AB_hz, 0]])

    transitions = _simulate_spin_system_exact([nu_A, nu_B], J_matrix)

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

    nu_A = delta_A_ppm * spectrometer_mhz
    nu_B = delta_B_ppm * spectrometer_mhz

    frequencies = [nu_A, nu_B, nu_B]

    J_matrix = np.array([
        [0, J_AB_hz, J_AB_hz],
        [J_AB_hz, 0, 0],
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

    nu_A = delta_A_ppm * spectrometer_mhz
    nu_B = delta_B_ppm * spectrometer_mhz

    frequencies = [nu_A, nu_A, nu_B, nu_B]

    J_matrix = np.array([
        [0,      J_AA_hz,       J_AB_hz,        J_AB_prime_hz],
        [J_AA_hz, 0,            J_AB_prime_hz,  J_AB_hz],
        [J_AB_hz, J_AB_prime_hz, 0,             J_BB_hz],
        [J_AB_prime_hz, J_AB_hz, J_BB_hz,       0]
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

    nu_A = delta_A_ppm * spectrometer_mhz
    nu_X = delta_X_ppm * spectrometer_mhz

    frequencies = [nu_A, nu_A, nu_X, nu_X, nu_X]

    J_matrix = np.zeros((5, 5))

    J_matrix[0, 1] = J_matrix[1, 0] = 0.0

    for i in range(2):
        for j in range(2, 5):
            J_matrix[i, j] = J_matrix[j, i] = J_AX_hz

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

    if not h_indices or not proton_shifts_ppm:
        return None

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

def _detect_virtual_coupling(
    molH: Chem.Mol,
    h_idx: int,
    partners: List[int],
    J_vals: List[float],
    proton_shifts_ppm: Dict[int, float],
    spectrometer_mhz: float
) -> bool:

    if len(partners) != 2 or len(J_vals) != 2:
        return False

    h_a, h_b = partners

    try:
        parent_a = next(n.GetIdx() for n in molH.GetAtomWithIdx(h_a).GetNeighbors())
        parent_b = next(n.GetIdx() for n in molH.GetAtomWithIdx(h_b).GetNeighbors())
        if parent_a == parent_b:

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

_HUGGINS_ELECTRONEGATIVITY = {
    'H': 2.20,
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

_GARBISCH_ALKENE_PARAMS = {

    'standard': (11.0, -3.0, 2.0),
    'conjugated': (11.5, -3.5, 2.0),
    'strained': (9.0, -2.5, 1.5),
    'electronegative': (10.0, -2.5, 1.5),
}

_BYSTROV_PEPTIDE_PARAMS = {

    'HNCH_alpha': (7.9, -1.05, 0.65),
    'HNCH_psi': (5.8, -1.2, 1.4),
    'HCCH_peptide': (9.5, -1.6, 1.8),
}

_HAASNOOT_7_PARAMS = {
    'P1': 13.70,
    'P2': -0.73,
    'P3': 0.00,
    'P4': 0.56,
    'P5': -2.47,
    'P6': 16.9,
    'P7': -0.47,
}

def _karplus_garbisch_alkene(
    theta_deg: float,
    alkene_type: str = 'standard',
    substituents: Optional[List[str]] = None
) -> float:

    A, B, C = _GARBISCH_ALKENE_PARAMS.get(alkene_type, _GARBISCH_ALKENE_PARAMS['standard'])

    theta = math.radians(theta_deg)

    J = A + B * math.cos(theta) + C * math.cos(2 * theta)

    if substituents:
        en_correction = 0.0
        for sub in substituents:
            chi = _HUGGINS_ELECTRONEGATIVITY.get(sub, 2.5)
            chi_H = _HUGGINS_ELECTRONEGATIVITY.get('H', 0.0)
            delta_chi = chi - chi_H
            if delta_chi > 0.5:

                en_correction -= 0.3 * delta_chi
        J += en_correction

    return max(0.0, J)

def _karplus_bystrov_peptide(
    phi_deg: float,
    coupling_type: str = 'HNCH_alpha',
    psi_deg: Optional[float] = None
) -> float:

    A, B, C = _BYSTROV_PEPTIDE_PARAMS.get(coupling_type, _BYSTROV_PEPTIDE_PARAMS['HNCH_alpha'])

    if coupling_type == 'HNCH_alpha':

        theta = math.radians(phi_deg - 60.0)
    elif coupling_type == 'HNCH_psi' and psi_deg is not None:

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

    try:
        if not mol.GetNumConformers():
            return 1

        if conf_id < 0:
            conf = mol.GetConformer()
        else:
            try:
                conf = mol.GetConformer(conf_id)
            except Exception:
                conf = mol.GetConformer()

        dihed = AllChem.GetDihedralDeg(conf, h1_idx, c1_idx, c2_idx, substituent_idx)

        while dihed > 180:
            dihed -= 360
        while dihed < -180:
            dihed += 360

        abs_dihed = abs(dihed)

        if abs_dihed < 75.0:

            return 1
        elif abs_dihed > 105.0:

            return -1
        else:

            try:
                pos_h1 = np.array(conf.GetAtomPosition(h1_idx))
                pos_c1 = np.array(conf.GetAtomPosition(c1_idx))
                pos_c2 = np.array(conf.GetAtomPosition(c2_idx))
                pos_sub = np.array(conf.GetAtomPosition(substituent_idx))

                v_c1c2 = pos_c2 - pos_c1
                v_c1h1 = pos_h1 - pos_c1
                v_c1sub = pos_sub - pos_c1

                normal = np.cross(v_c1c2, v_c1h1)

                proj = np.dot(v_c1sub, normal)

                return 1 if proj >= 0 else -1
            except Exception:

                return 1 if abs_dihed < 90 else -1

    except Exception:
        return 1

def _karplus_haasnoot_altona_full(
    theta_deg: float,
    substituents_alpha: List[str] = None,
    substituents_beta: List[str] = None,
    orientation_signs: List[int] = None,
    use_7_params: bool = True,
) -> float:

    P1 = _HAASNOOT_7_PARAMS['P1']
    P2 = _HAASNOOT_7_PARAMS['P2']
    P3 = _HAASNOOT_7_PARAMS['P3']
    P4 = _HAASNOOT_7_PARAMS['P4']
    P5 = _HAASNOOT_7_PARAMS['P5']
    P6_deg = _HAASNOOT_7_PARAMS['P6']
    P7 = _HAASNOOT_7_PARAMS['P7']

    phi = math.radians(theta_deg)
    P6 = math.radians(P6_deg)

    J = P1 * (math.cos(phi) ** 2) + P2 * math.cos(phi) + P3

    if substituents_alpha is None:
        substituents_alpha = []
    if substituents_beta is None:
        substituents_beta = []
    if orientation_signs is None:

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

    for i, sub in enumerate(all_subs):
        chi_i = _HUGGINS_ELECTRONEGATIVITY.get(sub, 2.5)
        delta_chi = chi_i - chi_H
        delta_chi_list.append(delta_chi)

        if delta_chi <= 0:
            continue

        s_i = orientation_signs[i] if i < len(orientation_signs) else 1

        correction = delta_chi * (
            P4 + P5 * (math.cos(s_i * phi + P6 * abs(delta_chi)) ** 2)
        )
        J += correction

    if use_7_params and len(delta_chi_list) >= 2:

        product_delta_chi = 1.0
        n_significant = 0
        for dc in delta_chi_list:
            if dc > 0.1:
                product_delta_chi *= dc
                n_significant += 1

        if n_significant >= 2:

            normalized_product = product_delta_chi ** (1.0 / n_significant)
            J += P7 * normalized_product * n_significant

    return max(0.0, J)

def _karplus_with_substituent_correction(
    theta_deg: float,
    alpha_EN: float = 0.0,
    beta_EN: float = 0.0
) -> float:

    subs_alpha = ['O'] if alpha_EN > 0.3 else (['N'] if alpha_EN > 0.1 else [])
    subs_beta = ['O'] if beta_EN > 0.3 else (['N'] if beta_EN > 0.1 else [])

    return _karplus_haasnoot_altona_full(
        theta_deg,
        substituents_alpha=subs_alpha,
        substituents_beta=subs_beta
    )

def _add_combination_lines(
    lines: List[Tuple[float, float]],
    system_type: Optional[str],
    J_AB: Optional[float],
    config: "SimConfig",
    delta_nu_hz: float = 0.0,
) -> List[Tuple[float, float]]:

    if not lines or system_type is None:
        return lines

    if system_type not in ["AA'BB'", "AA'XX'", "AB", "ABX"]:
        return lines

    if J_AB is None or J_AB <= 0.0:
        return lines

    if delta_nu_hz > 0.0:
        D = math.sqrt(delta_nu_hz ** 2 + J_AB ** 2)
    else:
        D = J_AB

    if delta_nu_hz > 0.0:
        tan_2theta = J_AB / delta_nu_hz
        theta_2 = math.atan(tan_2theta)
        I_comb = (math.sin(theta_2) ** 2) / 4.0
    else:
        I_comb = 0.25

    min_intensity = getattr(config, 'combination_line_threshold', 0.02)

    if I_comb < min_intensity:
        return lines

    if len(lines) >= 2:
        center = sum(f for f, _ in lines) / len(lines)
    else:
        center = 0.0

    new_lines = list(lines)
    new_lines.append((center + D / 2.0, I_comb * 0.5))
    new_lines.append((center - D / 2.0, I_comb * 0.5))

    total = sum(i for _, i in new_lines)
    if total > 0:
        new_lines = [(f, i / total) for f, i in new_lines]

    return sorted(new_lines, key=lambda x: x[0])

def _adaptive_digital_resolution(spectrometer_mhz: float) -> float:

    SW_hz = 16.0 * float(spectrometer_mhz)
    num_points = 32768
    resolution_hz = SW_hz / num_points
    return float(resolution_hz)

def _first_order_peak_positions(
    J_couplings: List[Dict[str, float]],
    base_multiplicity: str
) -> List[Tuple[float, float]]:

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

TYPICAL_LINEWIDTH_HZ = 0.8

FREQUENCY_RESOLUTION_THRESHOLDS: Dict[int, Dict[str, float]] = {

    200: {
        "min_resolvable": 1.5,
        "min_broadened": 1.0,
        "4J_visible": 2.0,
        "5J_visible": 999.0,
        "delta_J_for_dd": 1.5,
    },

    300: {
        "min_resolvable": 1.2,
        "min_broadened": 0.8,
        "4J_visible": 1.5,
        "5J_visible": 999.0,
        "delta_J_for_dd": 1.2,
    },

    400: {
        "min_resolvable": 1.0,
        "min_broadened": 0.8,
        "4J_visible": 1.0,
        "5J_visible": 2.0,
        "delta_J_for_dd": 1.0,
    },

    500: {
        "min_resolvable": 0.8,
        "min_broadened": 0.6,
        "4J_visible": 0.8,
        "5J_visible": 1.5,
        "delta_J_for_dd": 0.8,
    },

    600: {
        "min_resolvable": 0.8,
        "min_broadened": 0.5,
        "4J_visible": 0.6,
        "5J_visible": 1.0,
        "delta_J_for_dd": 0.5,
    },

    700: {
        "min_resolvable": 0.8,
        "min_broadened": 0.5,
        "4J_visible": 0.5,
        "5J_visible": 0.8,
        "delta_J_for_dd": 0.5,
    },

    800: {
        "min_resolvable": 0.8,
        "min_broadened": 0.4,
        "4J_visible": 0.4,
        "5J_visible": 0.6,
        "delta_J_for_dd": 0.4,
    },

    900: {
        "min_resolvable": 0.8,
        "min_broadened": 0.4,
        "4J_visible": 0.4,
        "5J_visible": 0.5,
        "delta_J_for_dd": 0.4,
    },

    1000: {
        "min_resolvable": 0.8,
        "min_broadened": 0.3,
        "4J_visible": 0.3,
        "5J_visible": 0.5,
        "delta_J_for_dd": 0.3,
    },

    1200: {
        "min_resolvable": 0.8,
        "min_broadened": 0.3,
        "4J_visible": 0.3,
        "5J_visible": 0.4,
        "delta_J_for_dd": 0.3,
    },
}

def get_resolution_thresholds(frequency_mhz: float) -> Dict[str, float]:

    freqs = sorted(FREQUENCY_RESOLUTION_THRESHOLDS.keys())

    if frequency_mhz <= freqs[0]:
        return FREQUENCY_RESOLUTION_THRESHOLDS[freqs[0]]
    if frequency_mhz >= freqs[-1]:
        return FREQUENCY_RESOLUTION_THRESHOLDS[freqs[-1]]

    lower_freq = freqs[0]
    upper_freq = freqs[-1]
    for i, f in enumerate(freqs):
        if f <= frequency_mhz:
            lower_freq = f
        if f >= frequency_mhz:
            upper_freq = f
            break

    if lower_freq == upper_freq:
        return FREQUENCY_RESOLUTION_THRESHOLDS[lower_freq]

    lower_vals = FREQUENCY_RESOLUTION_THRESHOLDS[lower_freq]
    upper_vals = FREQUENCY_RESOLUTION_THRESHOLDS[upper_freq]

    ratio = (frequency_mhz - lower_freq) / (upper_freq - lower_freq)

    result = {}
    for key in lower_vals:
        result[key] = lower_vals[key] + ratio * (upper_vals[key] - lower_vals[key])

    return result

def filter_J_by_frequency(
    J_values: List[float],
    frequency_mhz: float,
    coupling_orders: Optional[List[int]] = None,
    config: Optional["SimConfig"] = None,
    return_broadening_info: bool = False
) -> Union[Tuple[List[float], List[int]], Tuple[List[float], List[int], List[bool]]]:

    if not J_values:
        if return_broadening_info:
            return [], [], []
        return [], []

    thresholds = get_resolution_thresholds(frequency_mhz)
    min_resolvable = thresholds["min_resolvable"]
    min_broadened = thresholds.get("min_broadened", min_resolvable * 0.7)
    min_4J = thresholds["4J_visible"]
    min_5J = thresholds["5J_visible"]

    if config is not None:
        if hasattr(config, 'frequency_min_resolvable_override') and config.frequency_min_resolvable_override:
            min_resolvable = config.frequency_min_resolvable_override
            min_broadened = min_resolvable * 0.7

    filtered_J = []
    kept_indices = []
    is_broadened = []

    for i, J in enumerate(J_values):

        order = coupling_orders[i] if coupling_orders and i < len(coupling_orders) else 3

        if order == 2:
            threshold_net = min_resolvable
            threshold_broad = min_broadened
        elif order == 3:
            threshold_net = min_resolvable
            threshold_broad = min_broadened
        elif order == 4:
            threshold_net = min_4J
            threshold_broad = min_4J * 0.7
        elif order >= 5:
            threshold_net = min_5J
            threshold_broad = min_5J * 0.8
        else:
            threshold_net = min_resolvable
            threshold_broad = min_broadened

        J_abs = abs(J)

        if J_abs >= threshold_broad:
            filtered_J.append(J)
            kept_indices.append(i)

            is_broadened.append(J_abs < threshold_net)

    if return_broadening_info:
        return filtered_J, kept_indices, is_broadened
    return filtered_J, kept_indices

def compute_multiplicity_from_J(
    J_values: List[float],
    frequency_mhz: float,
    coupling_orders: Optional[List[int]] = None,
    config: Optional["SimConfig"] = None
) -> str:

    resolved_J, _ = filter_J_by_frequency(
        J_values, frequency_mhz, coupling_orders, config
    )

    if not resolved_J:
        return 's'

    thresholds = get_resolution_thresholds(frequency_mhz)
    delta_J_min = thresholds["delta_J_for_dd"]

    resolved_J = sorted(resolved_J, reverse=True)
    n_J = len(resolved_J)

    if n_J == 1:
        return 'd'

    groups = []
    current_group = [resolved_J[0]]

    for i in range(1, n_J):
        if abs(resolved_J[i] - current_group[-1]) < delta_J_min:
            current_group.append(resolved_J[i])
        else:
            groups.append(current_group)
            current_group = [resolved_J[i]]
    groups.append(current_group)

    n_groups = len(groups)
    group_sizes = [len(g) for g in groups]

    if n_groups == 1:

        n_equiv = group_sizes[0]
        if n_equiv == 1:
            return 'd'
        elif n_equiv == 2:
            return 't'
        elif n_equiv == 3:
            return 'q'
        elif n_equiv >= 4:

            return 'm'
        else:
            return 'm'

    elif n_groups == 2:

        s1, s2 = group_sizes[0], group_sizes[1]

        if s1 == 1 and s2 == 1:
            return 'dd'
        elif s1 == 1 and s2 == 2:
            return 'dt'
        elif s1 == 2 and s2 == 1:
            return 'td'
        elif s1 == 1 and s2 == 3:
            return 'dq'
        elif s1 == 3 and s2 == 1:
            return 'qd'
        elif s1 == 2 and s2 == 2:
            return 'tt'
        elif s1 == 2 and s2 == 3:
            return 'tq'
        elif s1 == 3 and s2 == 2:
            return 'qt'
        else:
            return 'm'

    elif n_groups == 3:

        s1, s2, s3 = group_sizes[0], group_sizes[1], group_sizes[2]

        if s1 == 1 and s2 == 1 and s3 == 1:
            return 'ddd'
        else:

            return 'm'

    elif n_groups == 4:
        if all(s == 1 for s in group_sizes):
            return 'dddd'
        else:
            return 'm'

    else:
        return 'm'

@dataclass
class SimConfig:

    use_frequency_dependent_resolution: bool = True

    frequency_min_resolvable_override: Optional[float] = None

    frequency_strict_mode: bool = False

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

    vinylic_J_defaults: Dict[str, float] = field(default_factory=lambda: {
        'trans': 16.0,
        'cis': 10.0,
        'geminal': 2.0,
        'gauche': 8.0,
    })
    combination_line_threshold: float = 0.02
    use_field_resolution_model: bool = True
    digital_resolution_hz: Optional[float] = None
    min_visibility_factor: float = 1.2
    grouping_tol_scale: float = 0.65
    grouping_tol_scale_auto: bool = True

    never_merge_gap_hz: float = 2.5

    separate_geminal_vicinal_J: bool = True

    detect_para_min_mhz: int = 600
    detect_para_min_j_hz: float = 0.45
    para_max_per_proton: Optional[int] = None
    sugar_mode: bool = True
    cyclohexane_mode: bool = True
    rigid_W_boost: bool = False
    formyl_rule: bool = True
    amide_mode: bool = True
    peptide_HNCH_coupling: bool = True
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

    use_garbisch_alkenes: bool = True
    use_bystrov_peptides: bool = True
    use_haasnoot_7_params: bool = True

    use_conformational_averaging: bool = True
    n_conformers: int = 200
    n_conformers_max: int = 300
    conformer_energy_cutoff_kcal: float = 10.0
    conformer_pruning_rmsd: float = 0.25
    flexible_chain_boost: bool = True
    flexible_chain_min_rotatable: int = 4
    use_etkdg_v3: bool = True
    conformer_random_seed: int = 42

    LARGE_MOLECULE_THRESHOLD: int = 40

    auto_disable_averaging_above_atoms: int = 40
    fast_mode_for_large_molecules: bool = True

    skip_3d_for_large_molecules: bool = True
    skip_3d_above_heavy_atoms: int = 50
    skip_3d_above_fused_rings: int = 5

    ultra_fast_mode: bool = True
    ultra_fast_above_heavy_atoms: int = 40
    ultra_fast_default_vicinal_J: float = 7.0

    smart_conformer_mode: bool = True

    advanced_spin_system_detection: bool = True
    spin_system_delta_nu_J_threshold: float = 10.0

    use_exact_spin_simulation: bool = True
    exact_spin_max_spins: int = 8
    exact_spin_intensity_threshold: float = 0.001

    use_spin_fragmentation: bool = True
    fragmentation_threshold_ratio: float = 0.05

    use_solvent_shift_corrections: bool = False
    use_solvent_linewidth: bool = True
    use_solvent_J_corrections: bool = False
    use_exchange_broadening: bool = True
    use_exchange_visibility: bool = True
    pH: Optional[float] = None

    mark_labile_protons_broad: bool = True
    labile_proton_linewidth_hz: float = 5.0

    force_field: str = "MMFF94s"
    force_field_max_iters: int = 200

    infer_double_bond_stereo: bool = True
    infer_stereo_default: str = "E"
    infer_stereo_log: bool = True

def _optimize_molecule_geometry(
    mol: Chem.Mol,
    force_field: str = "MMFF94s",
    max_iters: int = 200,
    conf_id: int = -1
) -> bool:

    try:
        if force_field.upper() in ("MMFF94", "MMFF94S"):

            variant = "MMFF94s" if "S" in force_field.upper() else "MMFF94"

            mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=variant)
            if mmff_props is None:

                return _optimize_molecule_geometry(mol, "UFF", max_iters, conf_id)

            if conf_id < 0:

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

            return _optimize_molecule_geometry(mol, "MMFF94s", max_iters, conf_id)

    except Exception:

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

    atom = mol.GetAtomWithIdx(heavy_idx)
    symbol = atom.GetSymbol()

    if symbol not in ("O", "N", "S"):
        return None

    has_h = any(n.GetSymbol() == "H" for n in atom.GetNeighbors())
    if not has_h:
        return None

    proton_type = None
    if symbol == "O":

        is_carboxylic = any(
            n.GetSymbol() == "C" and
            any(b.GetBondType() == rdchem.BondType.DOUBLE and b.GetOtherAtom(n).GetSymbol() == "O"
                for b in n.GetBonds())
            for n in atom.GetNeighbors() if n.GetSymbol() == "C"
        )
        proton_type = "COOH" if is_carboxylic else "OH"
    elif symbol == "N":

        is_amide = any(
            n.GetSymbol() == "C" and
            any(b.GetBondType() == rdchem.BondType.DOUBLE and b.GetOtherAtom(n).GetSymbol() == "O"
                for b in n.GetBonds())
            for n in atom.GetNeighbors() if n.GetSymbol() == "C"
        )
        proton_type = "amide_NH" if is_amide else "NH"
    elif symbol == "S":
        proton_type = "SH"

    solvent_upper = solvent_name.upper().replace("-", "")

    if solvent_upper in ("D2O", "CD3OD", "MEOD"):

        if proton_type in ("OH", "COOH", "NH", "SH"):
            return "exchanged"
        elif proton_type == "amide_NH":
            return "very_broad"

    elif solvent_upper in ("DMSOD6", "DMSO"):

        if proton_type == "OH":
            return "sharp"
        elif proton_type == "COOH":
            return "broad"
        elif proton_type in ("NH", "amide_NH"):
            return "sharp"
        elif proton_type == "SH":
            return "broad"

    elif solvent_upper in ("CDCL3", "CDCl3"):

        if proton_type == "OH":
            return "broad"
        elif proton_type == "COOH":
            return "very_broad"
        elif proton_type == "NH":
            return "broad"
        elif proton_type == "amide_NH":
            return "sharp"
        elif proton_type == "SH":
            return "broad"

    return "broad" if proton_type else None

def _topological_distances(mol: Chem.Mol) -> np.ndarray:

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

def _canonical_ranks(mol_noH: Chem.Mol) -> List[int]:

    return MoleculeCache.get_canonical_ranks(mol_noH)

def _get_atom_invariant(atom: Chem.Atom) -> Tuple:

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

    if i == j:
        return True

    if mol is not None:
        try:
            symmetry_classes = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
            if i < len(symmetry_classes) and j < len(symmetry_classes):

                if symmetry_classes[i] == symmetry_classes[j]:
                    return True

                else:
                    return False
        except Exception:
            pass

    if ranks is not None and i < len(ranks) and j < len(ranks):
        if ranks[i] == ranks[j]:
            return True

    return False

def _are_protons_on_equivalent_carbons(
    mol: Chem.Mol,
    ranks: List[int],
    heavy_idx_1: int,
    heavy_idx_2: int,
) -> bool:

    if mol is None:
        return False

    if heavy_idx_1 == heavy_idx_2:
        return False

    try:
        symmetry_classes = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    except Exception:

        if ranks is None:
            return False
        symmetry_classes = ranks

    if heavy_idx_1 >= len(symmetry_classes) or heavy_idx_2 >= len(symmetry_classes):
        return False

    if symmetry_classes[heavy_idx_1] == symmetry_classes[heavy_idx_2]:
        return True

    return False

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

    atom = mol.GetAtomWithIdx(heavy_idx)

    if atom.GetSymbol() != "C":
        return False

    n_H = len(heavy_to_H.get(heavy_idx, []))
    if n_H != 3:
        return False

    if _is_sp2_like(atom):
        return False

    neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != "H"]
    if len(neighbors) != 1:
        return False

    neighbor = neighbors[0]
    neighbor_sym = neighbor.GetSymbol()

    if neighbor_sym == "C" and neighbor.GetIsAromatic():
        return True

    if neighbor_sym == "O":
        return True

    if neighbor_sym == "N":
        return True

    if neighbor_sym == "C" and _is_sp2_like(neighbor):

        for n2 in neighbor.GetNeighbors():
            if n2.GetSymbol() == "O":
                bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), n2.GetIdx())
                if bond and bond.GetBondTypeAsDouble() >= 2.0:

                    other_neighbors = [x for x in neighbor.GetNeighbors()
                                       if x.GetIdx() != heavy_idx and x.GetSymbol() not in ("O", "H")]
                    if other_neighbors:
                        return True

    if neighbor_sym == "C" and not _is_sp2_like(neighbor):

        neighbor_idx = neighbor.GetIdx()
        n_H_on_neighbor = len(heavy_to_H.get(neighbor_idx, []))
        if n_H_on_neighbor == 0:

            n_methyl_neighbors = 0
            for n2 in neighbor.GetNeighbors():
                if n2.GetSymbol() == "C" and not _is_sp2_like(n2):
                    n2_idx = n2.GetIdx()
                    n_H_on_n2 = len(heavy_to_H.get(n2_idx, []))

                    n_heavy_neighbors = len([x for x in n2.GetNeighbors() if x.GetSymbol() != "H"])
                    if n_H_on_n2 == 3 and n_heavy_neighbors == 1:
                        n_methyl_neighbors += 1
            if n_methyl_neighbors >= 2:

                return True

    return False

def _is_benzylic_CH2_4J_coupling(
    mol: Chem.Mol,
    heavy_idx1: int,
    heavy_idx2: int,
    heavy_to_H: Dict[int, List[int]]
) -> bool:

    atom1 = mol.GetAtomWithIdx(heavy_idx1)
    atom2 = mol.GetAtomWithIdx(heavy_idx2)

    benzylic_c = None
    aromatic_c = None

    if (atom1.GetSymbol() == "C" and not atom1.GetIsAromatic() and
        atom2.GetSymbol() == "C" and atom2.GetIsAromatic()):
        benzylic_c = atom1
        aromatic_c = atom2

    elif (atom2.GetSymbol() == "C" and not atom2.GetIsAromatic() and
          atom1.GetSymbol() == "C" and atom1.GetIsAromatic()):
        benzylic_c = atom2
        aromatic_c = atom1
    else:
        return False

    if _is_sp2_like(benzylic_c):
        return False

    n_H_benzylic = len(heavy_to_H.get(benzylic_c.GetIdx(), []))
    if n_H_benzylic == 0:
        return False

    benzylic_attached_to_aromatic = False
    for neighbor in benzylic_c.GetNeighbors():
        if neighbor.GetIsAromatic():
            benzylic_attached_to_aromatic = True
            break

    if not benzylic_attached_to_aromatic:
        return False

    n_H_aromatic = len(heavy_to_H.get(aromatic_c.GetIdx(), []))
    if n_H_aromatic == 0:
        return False

    return True

def _is_aldehyde_CHO(mol: Chem.Mol, heavy_idx: int, heavy_to_H: Dict[int, List[int]]) -> bool:

    atom = mol.GetAtomWithIdx(heavy_idx)

    if atom.GetSymbol() != "C":
        return False

    n_H = len(heavy_to_H.get(heavy_idx, []))
    if n_H != 1:
        return False

    if not _is_sp2_like(atom):
        return False

    has_carbonyl_O = False
    has_C_neighbor = False

    neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != "H"]

    for neighbor in neighbors:
        if neighbor.GetSymbol() == "O":
            bond = mol.GetBondBetweenAtoms(heavy_idx, neighbor.GetIdx())
            if bond and bond.GetBondTypeAsDouble() >= 2.0:
                has_carbonyl_O = True
        elif neighbor.GetSymbol() == "C":
            has_C_neighbor = True

    return has_carbonyl_O and has_C_neighbor

def _is_aldehyde_4J_coupling(
    mol: Chem.Mol,
    heavy_idx: int,
    other_heavy: int,
    heavy_to_H: Dict[int, List[int]],
    D: Any
) -> bool:

    dist = int(D[heavy_idx, other_heavy]) if D is not None else 99

    if dist != 2:
        return False

    if _is_aldehyde_CHO(mol, heavy_idx, heavy_to_H):

        other_atom = mol.GetAtomWithIdx(other_heavy)
        if other_atom.GetSymbol() == "C" and not _is_sp2_like(other_atom):
            return True

    if _is_aldehyde_CHO(mol, other_heavy, heavy_to_H):
        atom = mol.GetAtomWithIdx(heavy_idx)
        if atom.GetSymbol() == "C" and not _is_sp2_like(atom):
            return True

    return False

def _is_formamide_CHO(mol: Chem.Mol, heavy_idx: int, heavy_to_H: Dict[int, List[int]]) -> bool:

    atom = mol.GetAtomWithIdx(heavy_idx)

    if atom.GetSymbol() != "C":
        return False

    n_H = len(heavy_to_H.get(heavy_idx, []))
    if n_H != 1:
        return False

    if not _is_sp2_like(atom):
        return False

    has_carbonyl_O = False
    has_N_neighbor = False

    neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != "H"]

    for neighbor in neighbors:
        if neighbor.GetSymbol() == "O":

            bond = mol.GetBondBetweenAtoms(heavy_idx, neighbor.GetIdx())
            if bond and bond.GetBondTypeAsDouble() >= 2.0:
                has_carbonyl_O = True
        elif neighbor.GetSymbol() == "N":
            has_N_neighbor = True

    return has_carbonyl_O and has_N_neighbor

def _is_NCH3_coupling_to_formamide_CHO(
    mol: Chem.Mol,
    heavy_idx: int,
    other_heavy: int,
    heavy_to_H: Dict[int, List[int]]
) -> bool:

    if _is_formamide_CHO(mol, heavy_idx, heavy_to_H):

        other_atom = mol.GetAtomWithIdx(other_heavy)
        if other_atom.GetSymbol() != "C":
            return False

        n_H_other = len(heavy_to_H.get(other_heavy, []))
        if n_H_other != 3:
            return False

        for neighbor in other_atom.GetNeighbors():
            if neighbor.GetSymbol() == "N":
                return True

        return False

    if _is_formamide_CHO(mol, other_heavy, heavy_to_H):
        atom = mol.GetAtomWithIdx(heavy_idx)
        if atom.GetSymbol() != "C":
            return False

        n_H = len(heavy_to_H.get(heavy_idx, []))
        if n_H != 3:
            return False

        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == "N":
                return True

    return False

def _is_geminal_methyl(base: Chem.Mol, methyl_c_idx: int, heavy_to_H: Dict[int, List[int]]) -> bool:

    n_H = len(heavy_to_H.get(methyl_c_idx, []))
    if n_H != 3:
        return False

    methyl_atom = base.GetAtomWithIdx(methyl_c_idx)
    if methyl_atom.GetSymbol() != "C":
        return False

    for neighbor in methyl_atom.GetNeighbors():
        if neighbor.GetSymbol() != "C":
            continue
        parent_idx = neighbor.GetIdx()

        n_methyl_siblings = 0
        for sibling in neighbor.GetNeighbors():
            if sibling.GetIdx() == methyl_c_idx:
                continue
            if sibling.GetSymbol() != "C":
                continue

            n_H_sibling = len(heavy_to_H.get(sibling.GetIdx(), []))
            if n_H_sibling == 3:
                n_methyl_siblings += 1

        if n_methyl_siblings >= 1:
            return True

    return False

def _is_homoallylic_4J(base: Chem.Mol, i: int, j: int, mol3d: Optional[Chem.Mol] = None) -> Tuple[bool, float]:

    Ai = base.GetAtomWithIdx(i)
    Aj = base.GetAtomWithIdx(j)

    if Ai.GetSymbol() != "C" or Aj.GetSymbol() != "C":
        return False, 0.0

    def _get_geminal_methyl_J(base: Chem.Mol, vinyl_C_idx: int, double_bond_C_idx: int,
                               methyl_C_idx: int, mol3d: Optional[Chem.Mol] = None) -> float:

        double_bond_atom = base.GetAtomWithIdx(double_bond_C_idx)

        methyls_on_same_C = []
        for neigh in double_bond_atom.GetNeighbors():
            if neigh.GetSymbol() == "C" and not _is_sp2_like(neigh):

                n_H = neigh.GetTotalNumHs()
                if n_H == 3:
                    methyls_on_same_C.append(neigh.GetIdx())

        if len(methyls_on_same_C) >= 2:

            methyls_on_same_C.sort()

            if mol3d is not None and mol3d.GetNumConformers() > 0:
                try:
                    conf = mol3d.GetConformer()

                    vinyl_atom_3d = mol3d.GetAtomWithIdx(vinyl_C_idx)
                    vinyl_H_idx = None
                    for nn in vinyl_atom_3d.GetNeighbors():
                        if nn.GetSymbol() == "H":
                            vinyl_H_idx = nn.GetIdx()
                            break

                    if vinyl_H_idx is not None:
                        dihedral = abs(AllChem.GetDihedralDeg(
                            conf, vinyl_H_idx, vinyl_C_idx, double_bond_C_idx, methyl_C_idx))

                        if dihedral < 90:
                            return 1.50
                        else:
                            return 1.25
                except Exception:
                    pass

            if methyl_C_idx == methyls_on_same_C[0]:
                return 1.50
            else:
                return 1.25

        return 1.3

    if _is_sp2_like(Ai) and not _is_sp2_like(Aj):

        for n1 in Ai.GetNeighbors():
            if n1.GetSymbol() != "C" or not _is_sp2_like(n1):
                continue
            bond = base.GetBondBetweenAtoms(Ai.GetIdx(), n1.GetIdx())
            if not bond or bond.GetBondTypeAsDouble() < 1.5:
                continue

            if base.GetBondBetweenAtoms(n1.GetIdx(), Aj.GetIdx()):

                J = _get_geminal_methyl_J(base, Ai.GetIdx(), n1.GetIdx(), Aj.GetIdx(), mol3d)
                return True, J

    if _is_sp2_like(Aj) and not _is_sp2_like(Ai):
        for n1 in Aj.GetNeighbors():
            if n1.GetSymbol() != "C" or not _is_sp2_like(n1):
                continue
            bond = base.GetBondBetweenAtoms(Aj.GetIdx(), n1.GetIdx())
            if not bond or bond.GetBondTypeAsDouble() < 1.5:
                continue
            if base.GetBondBetweenAtoms(n1.GetIdx(), Ai.GetIdx()):
                J = _get_geminal_methyl_J(base, Aj.GetIdx(), n1.GetIdx(), Ai.GetIdx(), mol3d)
                return True, J

    return False, 0.0

def _is_allylic_4J(base: Chem.Mol, i: int, j: int, mol3d: Optional[Chem.Mol] = None) -> Tuple[bool, float]:

    Ai = base.GetAtomWithIdx(i)
    Aj = base.GetAtomWithIdx(j)

    if Ai.GetSymbol() != "C" or Aj.GetSymbol() != "C":
        return False, 0.0

    if not _is_sp2_like(Ai) and _is_sp2_like(Aj):
        for n1 in Ai.GetNeighbors():
            if n1.GetSymbol() != "C" or not _is_sp2_like(n1):
                continue

            bond = base.GetBondBetweenAtoms(n1.GetIdx(), Aj.GetIdx())
            if bond and bond.GetBondTypeAsDouble() >= 1.5:

                J = 1.8
                if mol3d is not None and mol3d.GetNumConformers() > 0:
                    try:
                        conf = mol3d.GetConformer()

                        dihedral = abs(AllChem.GetDihedralDeg(conf, i, n1.GetIdx(), Aj.GetIdx(), j))

                        if dihedral < 30 or dihedral > 150:
                            J = 2.5
                        else:
                            J = 1.2
                    except Exception:
                        pass
                return True, J

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

    Ai = base.GetAtomWithIdx(i)
    Aj = base.GetAtomWithIdx(j)

    if not (Ai.IsInRing() and Aj.IsInRing()):
        return False, 0.0

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

                    if n1.IsInRing() and n2.IsInRing() and n3.IsInRing():
                        return True, 1.2

    return False, 0.0

def _is_W_allylic(base: Chem.Mol, i: int, j: int) -> Tuple[bool, bool]:

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

_KARPLUS_COEFFICIENTS = {

    'HCCH_sp3': (7.76, -1.10, 1.40),
    'HCCH_alkane': (9.50, -1.00, 0.50),
    'HCCH_altona': (13.70, -0.73, 0.00),
    'HCCH_sugar': (10.40, -1.50, 0.20),
    'HNCH': (6.40, -1.40, 1.90),
    'HCOH': (10.40, -1.50, 0.20),
    'HCNH': (8.00, -1.00, 0.00),
}

def _karplus_simple(theta_deg: float, bond_type: str = 'HCCH_altona') -> float:

    A, B, C = _KARPLUS_COEFFICIENTS.get(bond_type, _KARPLUS_COEFFICIENTS['HCCH_altona'])
    theta = math.radians(theta_deg)
    return max(0.0, A * (math.cos(theta) ** 2) + B * math.cos(theta) + C)

def _karplus_haasnoot_altona(theta_deg: float, electroneg_corr: float = 0.0) -> float:

    return _karplus_with_substituent_correction(theta_deg, alpha_EN=electroneg_corr, beta_EN=0.0)

def _get_substituents_on_carbon(mol: Chem.Mol, c_idx: int, exclude_idx: int) -> List[str]:

    c_atom = mol.GetAtomWithIdx(c_idx)
    substituents = []
    for nb in c_atom.GetNeighbors():
        if nb.GetIdx() == exclude_idx:
            continue
        sym = nb.GetSymbol()
        if sym != 'H':
            substituents.append(sym)
    return substituents

def _get_flexible_chain_fallback_J(
    mol: Chem.Mol,
    heavy_idx1: int,
    heavy_idx2: int
) -> Tuple[float, str]:

    heteroatoms_1 = []
    charges_1 = []
    for nb in mol.GetAtomWithIdx(heavy_idx1).GetNeighbors():
        if nb.GetIdx() == heavy_idx2:
            continue
        sym = nb.GetSymbol()
        if sym in ('O', 'N', 'F', 'Cl', 'S'):
            heteroatoms_1.append(sym)
            charges_1.append(nb.GetFormalCharge())

    heteroatoms_2 = []
    charges_2 = []
    for nb in mol.GetAtomWithIdx(heavy_idx2).GetNeighbors():
        if nb.GetIdx() == heavy_idx1:
            continue
        sym = nb.GetSymbol()
        if sym in ('O', 'N', 'F', 'Cl', 'S'):
            heteroatoms_2.append(sym)
            charges_2.append(nb.GetFormalCharge())

    has_O_side1 = 'O' in heteroatoms_1
    has_O_side2 = 'O' in heteroatoms_2
    has_F_side1 = 'F' in heteroatoms_1
    has_F_side2 = 'F' in heteroatoms_2

    if (has_O_side1 and has_O_side2) or (has_F_side1 and has_F_side2):
        return (5.0, "gauche_effect_double")

    if has_O_side1 or has_O_side2 or has_F_side1 or has_F_side2:

        if has_F_side1 or has_F_side2:
            return (5.5, "gauche_effect_F")
        else:
            return (6.5, "gauche_effect_O")

    has_cation = any(c > 0 for c in charges_1 + charges_2)
    if has_cation:
        return (7.5, "extended_chain_cation")

    return (7.0, "fallback_average")

def _estimate_vicinal_sp3_J(
    mol3d: Chem.Mol,
    h1: int,
    h2: int,
    use_haasnoot: bool = False,
    electroneg_corr: float = 0.0,
    use_full_haasnoot: bool = True,
) -> Optional[float]:

    try:
        c1 = next(n.GetIdx() for n in mol3d.GetAtomWithIdx(h1).GetNeighbors() if n.GetSymbol() == "C")
        c2 = next(n.GetIdx() for n in mol3d.GetAtomWithIdx(h2).GetNeighbors() if n.GetSymbol() == "C")

        if not mol3d.GetBondBetweenAtoms(c1, c2):
            return None

        conf = mol3d.GetConformer()
        dihed = AllChem.GetDihedralDeg(conf, h1, c1, c2, h2)

        if use_full_haasnoot:

            subs_alpha = _get_substituents_on_carbon(mol3d, c1, c2)
            subs_beta = _get_substituents_on_carbon(mol3d, c2, c1)

            orientation_signs = []

            c1_atom = mol3d.GetAtomWithIdx(c1)
            for nb in c1_atom.GetNeighbors():
                if nb.GetIdx() == c2 or nb.GetSymbol() == 'H':
                    continue
                sign = _haasnoot_orientation_sign(mol3d, h1, c1, c2, nb.GetIdx())
                orientation_signs.append(sign)

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

    try:
        n_atom = mol.GetAtomWithIdx(n_idx)

        carbonyl_c = None
        alpha_c = None

        for nb in n_atom.GetNeighbors():
            if nb.GetSymbol() != 'C':
                continue

            is_carbonyl = False
            for bond in nb.GetBonds():
                other = bond.GetOtherAtom(nb)
                if other.GetSymbol() == 'O' and bond.GetBondTypeAsDouble() >= 2:
                    is_carbonyl = True
                    break

            if is_carbonyl:
                carbonyl_c = nb
            else:

                if any(n2.GetSymbol() == 'H' for n2 in nb.GetNeighbors()):
                    alpha_c = nb

        if carbonyl_c is None or alpha_c is None:
            return False, None

        h_alpha = None
        for nb in alpha_c.GetNeighbors():
            if nb.GetSymbol() == 'H':
                h_alpha = nb
                break

        if h_alpha is None:
            return False, None

        if mol.GetNumConformers() > 0:
            conf = mol.GetConformer()

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

    try:
        is_peptide, phi = _is_peptide_backbone(mol, h_idx, n_idx)

        if not is_peptide:
            return None

        if phi is not None:
            return _karplus_bystrov_peptide(phi, coupling_type=coupling_type)
        else:

            return 6.0

    except Exception:
        return None

def _estimate_vinylic_J(mol3d: Chem.Mol, h1: int, h2: int, config: Optional[SimConfig]=None) -> Optional[float]:

    config = config or SimConfig()
    try:
        c1 = next(n.GetIdx() for n in mol3d.GetAtomWithIdx(h1).GetNeighbors() if n.GetSymbol() == "C")
        c2 = next(n.GetIdx() for n in mol3d.GetAtomWithIdx(h2).GetNeighbors() if n.GetSymbol() == "C")
        if not mol3d.GetBondBetweenAtoms(c1, c2):
            return None
        b = mol3d.GetBondBetweenAtoms(c1, c2)
        if not b or b.GetBondTypeAsDouble() < 2.0:
            return None

        if c1 == c2:
            return config.vinylic_J_defaults.get('geminal', 2.0)

        conf = mol3d.GetConformer()
        phi = AllChem.GetDihedralDeg(conf, h1, c1, c2, h2)

        c1_atom = mol3d.GetAtomWithIdx(c1)
        c2_atom = mol3d.GetAtomWithIdx(c2)

        is_conjugated = False
        has_electroneg_sub = False
        substituents = []

        for c_atom in [c1_atom, c2_atom]:
            for nb in c_atom.GetNeighbors():
                if nb.GetIdx() in (h1, h2, c1, c2):
                    continue
                sym = nb.GetSymbol()
                substituents.append(sym)

                if sym == 'C' and (_is_sp2_like(nb) or nb.GetIsAromatic()):
                    is_conjugated = True

                if sym in ('O', 'N', 'F', 'Cl', 'Br'):
                    has_electroneg_sub = True

        if has_electroneg_sub:
            alkene_type = 'electronegative'
        elif is_conjugated:
            alkene_type = 'conjugated'
        else:
            alkene_type = 'standard'

        J = _karplus_garbisch_alkene(phi, alkene_type=alkene_type, substituents=substituents)

        return J

    except Exception:
        return None

def _estimate_alkyne_J(mol: Chem.Mol, h1_idx: int, h2_idx: int) -> Optional[float]:

    try:

        if COUPLING_TABLES_MODULE_AVAILABLE:

            try:
                path = Chem.GetShortestPath(mol, h1_idx, h2_idx)
                if path and len(path) >= 2:
                    J_table = get_alkyne_coupling(mol, h1_idx, h2_idx, list(path))
                    if J_table is not None:
                        return J_table
            except Exception:
                pass

        h1_atom = mol.GetAtomWithIdx(h1_idx)
        h2_atom = mol.GetAtomWithIdx(h2_idx)

        c1_neighbors = [n for n in h1_atom.GetNeighbors() if n.GetSymbol() == 'C']
        c2_neighbors = [n for n in h2_atom.GetNeighbors() if n.GetSymbol() == 'C']

        if not c1_neighbors or not c2_neighbors:
            return None

        c1 = c1_neighbors[0]
        c2 = c2_neighbors[0]

        def is_alkyne_carbon(atom):
            for bond in atom.GetBonds():
                if bond.GetBondTypeAsDouble() == 3.0:
                    return True
            return False

        if is_alkyne_carbon(c1) and is_alkyne_carbon(c2):
            return 2.5

        for nb1 in c1.GetNeighbors():
            if nb1.GetSymbol() != 'C':
                continue
            bond_c1_nb1 = mol.GetBondBetweenAtoms(c1.GetIdx(), nb1.GetIdx())
            if bond_c1_nb1 and bond_c1_nb1.GetBondTypeAsDouble() == 3.0:

                for nb2 in nb1.GetNeighbors():
                    if nb2.GetIdx() == c1.GetIdx():
                        continue
                    if nb2.GetSymbol() == 'C':

                        bond_nb1_nb2 = mol.GetBondBetweenAtoms(nb1.GetIdx(), nb2.GetIdx())
                        if bond_nb1_nb2 and bond_nb1_nb2.GetBondTypeAsDouble() == 3.0:
                            if nb2.GetIdx() == c2.GetIdx():
                                return 2.5

                            if mol.GetBondBetweenAtoms(nb2.GetIdx(), c2.GetIdx()):
                                return 2.0

        for nb in c1.GetNeighbors():
            if nb.GetSymbol() == 'C':
                for bond in nb.GetBonds():
                    if bond.GetBondTypeAsDouble() == 3.0:

                        return 2.5

        return None

    except Exception:
        return None

def _estimate_31P_1H_coupling(mol: Chem.Mol, h_idx: int, p_idx: int, dist: int) -> Optional[float]:

    try:
        p_atom = mol.GetAtomWithIdx(p_idx)

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

            if is_phosphate:
                return 8.0
            else:
                return 12.0
        elif dist == 3:

            if is_phosphate:
                return 10.0
            else:
                return 15.0
        elif dist == 4:

            return 2.0

        return None

    except Exception:
        return None

def _estimate_HNCH_J(mol3d: Chem.Mol, h_idx: int, n_idx: int, c_idx: int, h2_idx: int) -> Optional[float]:

    try:
        conf = mol3d.GetConformer()
        phi = AllChem.GetDihedralDeg(conf, h_idx, n_idx, c_idx, h2_idx)

        return _karplus_simple(phi, bond_type='HNCH')

    except Exception:
        return None

def _infer_double_bond_stereo(mol: Chem.Mol, default_stereo: str = "E") -> Chem.Mol:

    import logging

    molH = Chem.AddHs(mol)

    n_heavy = molH.GetNumHeavyAtoms()
    if n_heavy > 50:
        logging.info(f"[SKIP 3D] Molecule trop grande ({n_heavy} atomes), pas d'EmbedMolecule")
        return mol

    try:
        try:
            ps = AllChem.ETKDGv3()
        except:
            ps = AllChem.ETKDG()
        ps.randomSeed = 42
        ps.useRandomCoords = True
        res = AllChem.EmbedMolecule(molH, ps)
        if res != 0:

            AllChem.EmbedMolecule(molH, AllChem.ETKDG())
    except Exception as e:
        logging.warning(f"Impossible de générer la 3D: {e}")
        return mol

    try:
        AllChem.MMFFOptimizeMolecule(molH, maxIters=200)
    except:
        try:
            AllChem.UFFOptimizeMolecule(molH, maxIters=200)
        except:
            pass

    modified_count = 0
    double_bonds_to_fix = []

    for bond in molH.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue

        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        if atom1.GetSymbol() != "C" or atom2.GetSymbol() != "C":
            continue

        current_stereo = bond.GetStereo()
        if current_stereo not in (Chem.BondStereo.STEREONONE,):
            continue

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

    if double_bonds_to_fix and molH.GetNumConformers() > 0:
        conf = molH.GetConformer()

        for db_info in double_bonds_to_fix:
            try:
                ref1_idx = db_info['ref1'].GetIdx()
                atom1_idx = db_info['atom1'].GetIdx()
                atom2_idx = db_info['atom2'].GetIdx()
                ref2_idx = db_info['ref2'].GetIdx()

                current_dihedral = AllChem.GetDihedralDeg(conf, ref1_idx, atom1_idx, atom2_idx, ref2_idx)

                if default_stereo == "E":
                    target_dihedral = 180.0
                else:
                    target_dihedral = 0.0

                if abs(current_dihedral - target_dihedral) > 30 and abs(current_dihedral + 360 - target_dihedral) > 30:

                    AllChem.SetDihedralDeg(conf, ref1_idx, atom1_idx, atom2_idx, ref2_idx, target_dihedral)
                    modified_count += 1
                elif abs(current_dihedral - target_dihedral) <= 30:
                    modified_count += 1

            except Exception as e:
                logging.debug(f"Erreur ajustement dièdre: {e}")
                continue

        if modified_count > 0:
            try:
                AllChem.MMFFOptimizeMolecule(molH, maxIters=50)
            except:
                pass

    if modified_count > 0:
        logging.info(f" Géométrie {default_stereo} (trans) forcée pour {modified_count} double(s) liaison(s)")

    Chem.AssignStereochemistryFrom3D(molH)

    return molH

def _count_undefined_double_bonds(mol: Chem.Mol) -> Tuple[int, int]:

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

def _validate_molecular_structure(mol: Chem.Mol, smiles: str) -> None:

    try:
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
        Chem.SetAromaticity(mol)
    except Exception as e:
        raise ValueError(f" Structure moléculaire invalide (valence, connectivité): {e}")

    molH = Chem.AddHs(mol)
    Chem.AssignStereochemistry(molH, cleanIt=True, force=True, flagPossibleStereoCenters=True)

    for atom in molH.GetAtoms():
        if atom.GetChiralTag() == rdchem.ChiralType.CHI_UNSPECIFIED:
            if atom.GetDegree() == 4 and atom.IsInRing():

                import logging
                logging.debug(f" Stéréochimie non définie sur l'atome {atom.GetIdx()} ({atom.GetSymbol()}) dans un cycle - sera ignorée")

    rings = Chem.GetSymmSSSR(molH)
    for ring in rings:
        if len(ring) == 5:
            atoms_in_ring = [molH.GetAtomWithIdx(i) for i in ring]
            symbols = [a.GetSymbol() for a in atoms_in_ring]

            if symbols.count("O") == 1 and symbols.count("C") == 4:

                for atom in molH.GetAtoms():
                    if atom.GetSymbol() == "C" and not atom.IsInRing():
                        neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
                        if "O" in neighbors and atom.GetTotalNumHs() >= 2:
                            raise ValueError(
                                f" Cycle furanose détecté (5 atomes) dans '{smiles}'. "
                                "Utilisez un cycle pyranose (6 atomes) pour les hexoses naturels."
                            )

def _validate_connectivity(molH: Chem.Mol, connectivity: Dict[int, List[int]]) -> bool:

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

def _delta_nu_hz(delta_ppm: float, spectrometer_mhz: float) -> float:
    return abs(delta_ppm) * spectrometer_mhz

def _aromatic_J_contextual(base: Chem.Mol, aromatic_dmap: List[Dict[int,int]], i: int, j: int,
                           J0_meta: float = 2.3, J0_para: float = 0.6) -> Tuple[Optional[float], Optional[str]]:

    rel = _aromatic_rel_from_map(aromatic_dmap, i, j)
    if rel is None:
        return None, rel

    cycle_type = _identify_heterocycle_type(base, i, j)

    if cycle_type == 'imidazole':
        ring_info = base.GetRingInfo()
        for ring in ring_info.AtomRings():
            if i in ring and j in ring:
                if _is_imidazole_h2_position(base, i, j, ring):
                    return 0.0, rel

    if cycle_type in _AROMATIC_COUPLING_TABLE:
        cycle_data = _AROMATIC_COUPLING_TABLE[cycle_type]
        if isinstance(cycle_data, dict) and rel in cycle_data:
            J = cycle_data[rel]
            return J, rel

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

    ring_info = mol.GetRingInfo()

    for ring in ring_info.AtomRings():
        if i in ring and j in ring:

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

        if ranks is not None:
            if _are_heavy_equivalent(ranks, heavy_idx, other_heavy, heavy_to_H, None, 0.01, mol=base):
                continue

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

def _rough_topological_J(base: Chem.Mol, aromatic_dmap: List[Dict[int,int]],
                         i_heavy: int, j_heavy: int, is_sp2_i: bool, is_sp2_j: bool) -> float:

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

def _group_J_values_with_orders(
    J_list: List[float],
    J_orders: List[int],
    tol: float = TOLERANCE_J_HZ,
    anisotropy_high: bool = False,
    config: Optional[SimConfig] = None,
    separate_geminal_vicinal: bool = True
) -> List[List[float]]:

    config = config or SimConfig()

    if not J_list or len(J_list) != len(J_orders):
        return _group_J_values(J_list, tol, anisotropy_high, config)

    if not separate_geminal_vicinal:
        return _group_J_values(J_list, tol, anisotropy_high, config)

    geminal_J = [(j, i) for i, (j, order) in enumerate(zip(J_list, J_orders)) if order == 2 and j is not None and j > 0]
    vicinal_J = [(j, i) for i, (j, order) in enumerate(zip(J_list, J_orders)) if order == 3 and j is not None and j > 0]
    other_J = [(j, i) for i, (j, order) in enumerate(zip(J_list, J_orders)) if order not in (2, 3) and j is not None and j > 0]

    def group_category(j_list):
        if not j_list:
            return []
        vals = [j for j, _ in j_list]
        return _group_J_values(vals, tol, anisotropy_high, config)

    geminal_groups = group_category(geminal_J)
    vicinal_groups = group_category(vicinal_J)
    other_groups = group_category(other_J)

    all_groups = geminal_groups + vicinal_groups + other_groups

    all_groups.sort(key=lambda g: float(np.mean(g)) if g else 0, reverse=True)

    return all_groups

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

def _apparent_sp2_label_from_J(J_vals: List[float], minJ: float = 0.4, tol: float = 0.5) -> str:

    Js = [j for j in J_vals if j is not None and j >= minJ]
    n = len(Js)
    if n == 0:
        return "s"
    if n == 1:
        return "d"

    if n >= 2:
        Js_sorted = sorted(Js)
        groups = []
        current_group = [Js_sorted[0]]

        for j in Js_sorted[1:]:

            ref_J = current_group[0]
            if ref_J < 2.0:
                effective_tol = 0.15
            else:
                effective_tol = tol

            if abs(j - ref_J) <= effective_tol:
                current_group.append(j)
            else:
                groups.append(current_group)
                current_group = [j]
        groups.append(current_group)

        total = sum(len(g) for g in groups)
        n_groups = len(groups)

        if n_groups == 1:

            if total == 2: return "t"
            if total == 3: return "q"
            if total >= 4:

                return "m"
            return "m"
        elif n_groups == 2:
            counts = sorted([len(g) for g in groups])
            if counts == [1, 1]: return "dd"
            if counts == [1, 2]: return "td"
            if counts == [1, 3]: return "qd"
            if counts == [2, 2]: return "tt"
            if counts == [2, 3]: return "tq"
            return "ddd"
        elif n_groups == 3:
            counts = sorted([len(g) for g in groups])

            if counts == [1, 1, 2]: return "ddt"
            if counts == [1, 1, 1]: return "ddd"
            if counts == [1, 2, 2]: return "tdt"
            return "m"
        else:
            return "m"

    return "ddd"

def _first_order_label_from_groups(groups: List[List[float]], force_singlet_if_empty: bool = True) -> str:

    if not groups:
        return "s" if force_singlet_if_empty else "m"

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
        if len(groups) == 1 and counts[0] == 4:

            return "m"
        if len(groups) == 2 and sorted(counts) == [1, 3]:
            return "qd"
        if len(groups) == 2 and sorted(counts) == [2, 2]: return "tt"
        if len(groups) == 3 and sorted(counts) == [1, 1, 2]:

            return "m"
        return "dddd"

    if total == 5:
        if len(groups) == 1 and counts[0] == 5:

            return "m"
        if len(groups) == 2 and sorted(counts) == [2, 3]:
            return "qt"
        if len(groups) == 3 and sorted(counts) == [1, 2, 2]:

            return "m"
        return "m"

    if total == 6:
        if len(groups) == 1 and counts[0] == 6:

            return "m"
        if len(groups) == 2 and sorted(counts) == [3, 3]:
            return "qt"
        if len(groups) == 3 and sorted(counts) == [1, 2, 3]:

            return "m"
        return "m"

    if total == 7:

        return "m"

    if total == 8:

        return "m"

    if total >= 9: return "m"

    return "s"

def _assign_stereo(mol: Chem.Mol) -> None:
    try:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    except Exception:
        pass

def _is_chiral_center(atom: Chem.Atom) -> bool:
    tag = atom.GetChiralTag()
    return tag in (rdchem.ChiralType.CHI_TETRAHEDRAL_CW, rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)

def _find_nearest_stereo_center(mol: Chem.Mol, start_idx: int, max_dist: int = 4) -> Tuple[Optional[int], int]:

    from collections import deque

    visited = {start_idx}
    queue = deque([(start_idx, 0)])

    while queue:
        current_idx, dist = queue.popleft()

        if dist > max_dist:
            continue

        atom = mol.GetAtomWithIdx(current_idx)

        if current_idx != start_idx and _is_chiral_center(atom):
            return current_idx, dist

        for neighbor in atom.GetNeighbors():
            nb_idx = neighbor.GetIdx()
            if nb_idx not in visited:
                visited.add(nb_idx)
                queue.append((nb_idx, dist + 1))

    return None, -1

def _has_stereogenic_path(mol: Chem.Mol, heavy_idx: int, max_dist: int = 3) -> bool:

    stereo_idx, dist = _find_nearest_stereo_center(mol, heavy_idx, max_dist)
    return stereo_idx is not None

def _is_diastereotopic_Hpair(molH: Chem.Mol, heavy_idx: int, h_a: int, h_b: int) -> bool:

    if heavy_idx is None or h_a is None or h_b is None:
        return False

    try:
        heavy_atom = molH.GetAtomWithIdx(heavy_idx)

        try:
            ranks = list(Chem.CanonicalRankAtoms(molH, includeChirality=True, breakTies=True))
        except TypeError:
            ranks = list(Chem.CanonicalRankAtoms(molH, includeChirality=True))

        if ranks[h_a] != ranks[h_b]:

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

        stereo_idx, stereo_dist = _find_nearest_stereo_center(molH, heavy_idx, max_dist=3)

        if stereo_idx is not None:

            if stereo_dist <= 2:
                return True

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

        if heavy_atom.GetSymbol() == "C":
            neighbors_heavy = [n for n in heavy_atom.GetNeighbors() if n.GetSymbol() != "H"]

            if heavy_atom.IsInRing():

                return True

            if _is_diastereotopic_by_prochirality(molH, heavy_idx, h_a, h_b):
                return True

            for nb in neighbors_heavy:
                if nb.GetHybridization() in (Chem.HybridizationType.SP2, Chem.HybridizationType.SP):

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

            if len(neighbors_heavy) == 2:
                all_sp3_acyclic = all(
                    nb.GetSymbol() in ("C", "O", "N") and
                    not nb.IsInRing() and
                    nb.GetHybridization() == Chem.HybridizationType.SP3
                    for nb in neighbors_heavy
                )
                if all_sp3_acyclic:

                    if not _is_diastereotopic_by_prochirality(molH, heavy_idx, h_a, h_b):

                        return False

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
    ultra_fast_mode: bool = False,
    diastereotopic_cache: Optional[Dict[int, Tuple[bool, float]]] = None,
) -> Tuple[
    List[float], bool, List[int], List[float], List[Tuple[int,float]],
    List[int], List[bool], str, Optional[str], float
]:

    if config is None: config = SimConfig()

    digital_floor = (config.digital_resolution_hz or 0.0)
    detect_floor = max(config.min_visibility_factor * config.line_width_hz, digital_floor)
    min2 = max(config.min_j_detect[2], detect_floor)
    min3 = max(config.min_j_detect[3], detect_floor)
    min4 = max(config.min_j_detect[4], detect_floor, 1.0)
    min5 = max(config.min_j_detect[5], detect_floor)

    J_vals: List[float] = []
    J_sources: List[Optional[int]] = []
    J_source_types: List[str] = []
    J_orders: List[int] = []
    partners: List[int] = []
    homo_J_vals_for_second_order: List[float] = []
    satellites: List[float] = []
    second_order = False
    second_order_type = "first_order"
    system_type: Optional[str] = None
    delta_nu_min_hz: float = 0.0

    homo_pairs: List[Tuple[int, float, str]] = []

    para_buf_J: List[float] = []
    para_buf_partners: List[int] = []

    if ultra_fast_mode:
        heavy_atom = base.GetAtomWithIdx(heavy_idx)
        is_sp2 = heavy_is_sp2_cache[heavy_idx] if heavy_is_sp2_cache else _is_sp2_like(heavy_atom)
        h_rep = h_atoms[0] if h_atoms else None

        n_heavy = base.GetNumAtoms()
        candidate_heavies = [j for j in range(n_heavy) if 1 <= int(D[heavy_idx, j]) <= 2]

        J_fast_vals = []
        partners_fast = []

        for other_heavy in candidate_heavies:
            if suppressed_heavy and other_heavy in suppressed_heavy:
                continue

            other_hs = heavy_to_H.get(other_heavy, [])
            if not other_hs:
                continue

            dist = int(D[heavy_idx, other_heavy])
            partner_atom = base.GetAtomWithIdx(other_heavy)
            partner_is_sp2 = heavy_is_sp2_cache[other_heavy] if heavy_is_sp2_cache else _is_sp2_like(partner_atom)

            if is_sp2 and partner_is_sp2 and heavy_atom.GetIsAromatic() and partner_atom.GetIsAromatic():
                rel = _aromatic_rel_from_map(aromatic_dmap, heavy_idx, other_heavy)
                if rel == "ortho":
                    J_fast_vals.append(config.aro_j_forced["ortho"])
                    partners_fast.append(other_hs[0])
                elif rel == "meta":
                    J_fast_vals.append(config.aro_j_forced["meta"])
                    partners_fast.append(other_hs[0])

            elif dist == 1 and not (is_sp2 and partner_is_sp2):
                for h in other_hs:
                    if h not in h_atoms:
                        J_fast_vals.append(config.ultra_fast_default_vicinal_J)
                        partners_fast.append(h)

        if geminal_partner is not None:
            J_fast_vals.append(12.0)
            partners_fast.append(geminal_partner)

        second_order_type_fast = "first_order"
        is_second_order = False

        if proton_shifts_ppm and h_rep is not None and J_fast_vals:

            my_shift_ppm = proton_shifts_ppm.get(h_rep)
            if my_shift_ppm is not None:

                min_ratio = float('inf')
                for partner_h, j_val in zip(partners_fast, J_fast_vals):

                    if partner_h in h_atoms:
                        continue

                    if j_val < 4.0:
                        continue
                    partner_shift = proton_shifts_ppm.get(partner_h)
                    if partner_shift is not None:
                        delta_nu_hz = abs(my_shift_ppm - partner_shift) * spectrometer_mhz
                        ratio = delta_nu_hz / j_val if j_val > 0 else float('inf')
                        if ratio < min_ratio:
                            min_ratio = ratio

                if min_ratio < float('inf'):
                    if min_ratio < 3.0:

                        second_order_type_fast = "strong_coupling"
                        is_second_order = True
                    elif min_ratio < 8.0:

                        second_order_type_fast = "moderate_second_order"
                        is_second_order = True

        return (
            J_fast_vals,
            is_second_order,
            partners_fast,
            [],
            [(p, j, "fallback_average") for p, j in zip(partners_fast, J_fast_vals)],
            [],
            [],
            second_order_type_fast,
            "ultra_fast",
            0.0,
            {}
        )

    processed_partners: set = set()

    heavy_atom = base.GetAtomWithIdx(heavy_idx)
    is_sp2 = heavy_is_sp2_cache[heavy_idx] if heavy_is_sp2_cache is not None else _is_sp2_like(heavy_atom)
    h_rep = h_atoms[0]

    partner_heavy_idx_kept: List[int] = []
    partner_is_sp2_kept: List[bool] = []

    if config.include_13C_satellites and heavy_atom.GetSymbol() == "C":
        nH_here = len(h_atoms)

        jch = _get_1JCH(heavy_atom, base)
        if nH_here >= 1:
            satellites.extend([jch] * nH_here)

    system_type = _detect_AABB_system(
        base, heavy_idx, heavy_to_H, aromatic_dmap,
        proton_shifts_ppm=proton_shifts_ppm,
        spectrometer_mhz=spectrometer_mhz,
        J_typical_hz=7.0
    )

    _heterocycle_handled = False
    if HETEROCYCLES_MODULE_AVAILABLE and heavy_atom.GetIsAromatic():
        try:
            het_result = estimate_heterocycle_couplings_complete(
                mol=molH,
                heavy_idx=heavy_idx,
                h_list=h_atoms,
                heavy_to_H=heavy_to_H,
                D=D,
                min_j_threshold=min3,
                apply_substituent_effects=True,
                spectrometer_mhz=spectrometer_mhz
            )

            if het_result is not None:

                het_J_vals = het_result.get('J_values', [])
                het_partners = het_result.get('coupling_partners', [])
                het_mult = het_result.get('multiplicity', 's')
                het_type = het_result.get('heterocycle_type', 'unknown')
                het_source = het_result.get('source', 'heterocycle')

                if het_source in ('heterocycle_singlet', 'heterocycle_expected_singlet') or (het_mult == 's' and not het_J_vals):

                    return (
                        [],
                        False,
                        [],
                        [],
                        [],
                        [],
                        [],
                        "first_order",
                        f"heterocycle_{het_type}",
                        0.0,
                        {}
                    )

                if het_J_vals:
                    het_J_vals_out = []
                    het_partners_out = []
                    het_homo_pairs: List[Tuple[int, float, str]] = []
                    het_partner_heavy_idx = []
                    het_partner_is_sp2 = []

                    for J_het, partner_het in zip(het_J_vals, het_partners):
                        het_J_vals_out.append(J_het)
                        het_partners_out.append(partner_het)
                        het_homo_pairs.append((partner_het, J_het, "heterocycle"))

                        for h_idx, h_list in heavy_to_H.items():
                            if partner_het in h_list:
                                het_partner_heavy_idx.append(h_idx)
                                partner_atom = base.GetAtomWithIdx(h_idx)
                                het_partner_is_sp2.append(partner_atom.GetIsAromatic() or _is_sp2_like(partner_atom))
                                break

                    het_position = het_result.get('position', '')
                    system_type_with_mult = f"heterocycle_{het_type}:{het_mult}:{het_position}"

                    return (
                        het_J_vals_out,
                        False,
                        het_partners_out,
                        [],
                        het_homo_pairs,
                        het_partner_heavy_idx,
                        het_partner_is_sp2,
                        "first_order",
                        system_type_with_mult,
                        0.0,
                        {}
                    )

        except Exception as e:

            _heterocycle_handled = False

    _saturated_heterocycle_handled = False

    _is_in_saturated_ring = False
    if HETEROCYCLES_MODULE_AVAILABLE and not heavy_atom.GetIsAromatic():
        try:
            _sat_type_check = identify_saturated_heterocycle(molH, heavy_idx)

            if _sat_type_check in ('diketopiperazine', 'lactam_5', 'lactam_6', 'saturated_heterocycle_fused',
                                   'diketopiperazine_adjacent', 'lactam_5_adjacent', 'lactam_6_adjacent',
                                   'dkp_pyrrolidine_fused', 'pyrrolidine_fused', 'pyrrolidine'):
                _is_in_saturated_ring = True
        except Exception:
            pass

    use_heterocycle_tables = _is_in_saturated_ring or (not config.use_conformational_averaging)

    if HETEROCYCLES_MODULE_AVAILABLE and not heavy_atom.GetIsAromatic() and not _heterocycle_handled and use_heterocycle_tables:
        try:
            sat_result = estimate_saturated_heterocycle_couplings(
                mol=molH,
                heavy_idx=heavy_idx,
                H_indices=h_atoms,
                heavy_to_H_mapping=heavy_to_H,
                topological_distances=D,
                spectrometer_mhz=spectrometer_mhz
            )

            if sat_result is not None:
                sat_J_vals = sat_result.get('J_values', [])
                sat_partners = sat_result.get('coupling_partners', [])
                sat_mult = sat_result.get('multiplicity', 's')
                sat_type = sat_result.get('saturated_type', 'unknown')
                max_dist = sat_result.get('max_distance', 3)

                sat_J_vals_out = []
                sat_partners_out = []
                sat_homo_pairs: List[Tuple[int, float, str]] = []
                sat_partner_heavy_idx = []
                sat_partner_is_sp2 = []

                partner_to_heavy = {}
                for p_sat in sat_partners:
                    p_atom = molH.GetAtomWithIdx(p_sat)
                    for nb in p_atom.GetNeighbors():
                        if nb.GetSymbol() != 'H':
                            partner_to_heavy[p_sat] = nb.GetIdx()
                            break

                heavy_to_partners_sat: Dict[int, List[int]] = {}
                for p_sat, h_idx_sat in partner_to_heavy.items():
                    heavy_to_partners_sat.setdefault(h_idx_sat, []).append(p_sat)

                heavy_atom_H = molH.GetAtomWithIdx(heavy_idx)
                current_h_siblings = [n.GetIdx() for n in heavy_atom_H.GetNeighbors() if n.GetSymbol() == 'H']
                current_h_siblings_sorted = sorted(current_h_siblings)

                for J_sat, partner_sat in zip(sat_J_vals, sat_partners):
                    if J_sat >= min3:
                        partner_heavy = partner_to_heavy.get(partner_sat)

                        J_final = J_sat
                        source_final = "saturated_heterocycle"

                        if len(current_h_siblings) == 2 and partner_heavy is not None:
                            try:
                                differentiate = should_differentiate_diastereotopic_J(
                                    molH, heavy_idx, partner_heavy, sat_type
                                )
                            except Exception:
                                differentiate = False

                            if differentiate:

                                if len(h_atoms) == 1:
                                    my_h_idx = h_atoms[0]

                                    if my_h_idx == current_h_siblings_sorted[0]:
                                        J_final = 6.5
                                    else:
                                        J_final = 10.0
                                    source_final = "saturated_heterocycle"

                        sat_J_vals_out.append(J_final)
                        sat_partners_out.append(partner_sat)

                        sat_homo_pairs.append((partner_sat, J_final, source_final))

                if sat_J_vals_out:
                    J_vals.extend(sat_J_vals_out)
                    J_sources.extend([None] * len(sat_J_vals_out))
                    J_orders.extend([3] * len(sat_J_vals_out))
                    partners.extend(sat_partners_out)

                    homo_J_vals_for_second_order.extend([hp[1] for hp in sat_homo_pairs])

                    if geminal_partner is not None:

                        Jg = _get_geminal_2J(heavy_atom, base)

                        if Jg >= min2:

                            geminal_shifts_differ = False
                            if proton_shifts_ppm is not None:
                                try:
                                    h_rep = h_atoms[0] if len(h_atoms) == 1 else h_atoms[0]
                                    d_rep = proton_shifts_ppm.get(h_rep, None)
                                    d_gem = proton_shifts_ppm.get(geminal_partner, None)
                                    if d_rep is not None and d_gem is not None:
                                        if math.isfinite(d_rep) and math.isfinite(d_gem):
                                            geminal_shifts_differ = abs(d_rep - d_gem) > 0.005
                                except Exception:
                                    pass

                            sat_partners_out.append(geminal_partner)
                            sat_homo_pairs.append((geminal_partner, Jg, "geminal_calc"))

                            if geminal_shifts_differ:
                                sat_J_vals_out.append(Jg)

                    _saturated_heterocycle_handled = True

                    return (
                        sat_J_vals_out,
                        False,
                        sat_partners_out,
                        [],
                        sat_homo_pairs,
                        [],
                        [],
                        "first_order",
                        f"saturated_{sat_type}",
                        0.0,
                        {}
                    )

        except Exception:

            _saturated_heterocycle_handled = False

    use_sugar_tables = not config.use_conformational_averaging
    if SUGARS_MODULE_AVAILABLE and heavy_atom.GetSymbol() == 'C' and not heavy_atom.GetIsAromatic() and use_sugar_tables:
        try:

            if is_sugar_carbon(base, heavy_idx):

                proton_to_heavy_local = {}
                for h_idx_map, h_list_map in heavy_to_H.items():
                    for h in h_list_map:
                        proton_to_heavy_local[h] = h_idx_map

                sugar_result = estimate_sugar_couplings_complete(
                    mol=molH,
                    heavy_idx=heavy_idx,
                    h_list=h_atoms,
                    heavy_to_H=heavy_to_H,
                    proton_to_heavy=proton_to_heavy_local,
                    D=D,
                    config=None,
                    spectrometer_mhz=spectrometer_mhz
                )

                if sugar_result is not None:
                    sugar_J_vals = sugar_result.get('J_values', [])
                    sugar_partners = sugar_result.get('coupling_partners', [])
                    sugar_mult = sugar_result.get('multiplicity', 's')
                    sugar_info = sugar_result.get('sugar_info', {})
                    sugar_position = sugar_info.get('position', '')
                    sugar_ring_type = sugar_info.get('ring_type', 'pyranose')
                    sugar_source = sugar_result.get('source', 'sugar_complete')

                    sugar_J_vals_out = []
                    sugar_partners_out = []
                    sugar_homo_pairs: List[Tuple[int, float, str]] = []
                    sugar_partner_heavy_idx = []
                    sugar_partner_is_sp2 = []

                    for J_sugar, partner_sugar in zip(sugar_J_vals, sugar_partners):
                        sugar_J_vals_out.append(J_sugar)
                        sugar_partners_out.append(partner_sugar)
                        sugar_homo_pairs.append((partner_sugar, J_sugar, "sugar"))

                        for h_idx, h_list in heavy_to_H.items():
                            if partner_sugar in h_list:
                                sugar_partner_heavy_idx.append(h_idx)
                                partner_atom = base.GetAtomWithIdx(h_idx)
                                sugar_partner_is_sp2.append(partner_atom.GetIsAromatic() or _is_sp2_like(partner_atom))
                                break

                    system_type_sugar = f"sugar_{sugar_ring_type}:{sugar_mult}:{sugar_position}"

                    return (
                        sugar_J_vals_out,
                        False,
                        sugar_partners_out,
                        [],
                        sugar_homo_pairs,
                        sugar_partner_heavy_idx,
                        sugar_partner_is_sp2,
                        "first_order",
                        system_type_sugar,
                        0.0,
                        {}
                    )

        except Exception as e:

            pass

    if COUPLING_EXTENSIONS_AVAILABLE and heavy_atom.GetSymbol() == 'C' and not heavy_atom.GetIsAromatic():

        proton_to_heavy_ext = {}
        for h_idx_map, h_list_map in heavy_to_H.items():
            for h in h_list_map:
                proton_to_heavy_ext[h] = h_idx_map

        all_h_on_heavy = heavy_to_H.get(heavy_idx, h_atoms)

        try:
            steroid_result = estimate_steroid_couplings_differentiated(
                mol=molH,
                heavy_idx=heavy_idx,
                h_list=all_h_on_heavy,
                heavy_to_H=heavy_to_H,
                proton_to_heavy=proton_to_heavy_ext,
            )

            if steroid_result is not None:
                steroid_type = steroid_result.get('system_type', 'steroid')
                steroid_couplings = steroid_result.get('couplings', {})

                h_rep_couplings = steroid_couplings.get(h_rep, [])

                if h_rep_couplings:
                    steroid_J_vals_out = []
                    steroid_partners_out = []
                    steroid_homo_pairs: List[Tuple[int, float, str]] = []

                    for partner_h, J_val, source in h_rep_couplings:
                        if J_val >= min3:
                            steroid_J_vals_out.append(J_val)
                            steroid_partners_out.append(partner_h)
                            steroid_homo_pairs.append((partner_h, J_val, source))

                    if steroid_J_vals_out:

                        if geminal_partner is not None:
                            Jg = _get_geminal_2J(heavy_atom, base)
                            if Jg >= min2:
                                steroid_partners_out.append(geminal_partner)
                                steroid_homo_pairs.append((geminal_partner, Jg, "geminal_calc"))

                                if proton_shifts_ppm is not None:
                                    try:
                                        d_rep = proton_shifts_ppm.get(h_rep, None)
                                        d_gem = proton_shifts_ppm.get(geminal_partner, None)
                                        if d_rep and d_gem and abs(d_rep - d_gem) > 0.005:
                                            steroid_J_vals_out.append(Jg)
                                    except Exception:
                                        pass

                        return (
                            steroid_J_vals_out,
                            False,
                            steroid_partners_out,
                            [],
                            steroid_homo_pairs,
                            [],
                            [],
                            "first_order",
                            steroid_type,
                            0.0,
                            {}
                        )
        except Exception:
            pass

        try:
            terpene_result = estimate_terpene_couplings_differentiated(
                mol=molH,
                heavy_idx=heavy_idx,
                h_list=all_h_on_heavy,
                heavy_to_H=heavy_to_H,
                proton_to_heavy=proton_to_heavy_ext,
            )

            if terpene_result is not None:
                terpene_type = terpene_result.get('system_type', 'terpene')
                terpene_couplings = terpene_result.get('couplings', {})

                h_rep_couplings = terpene_couplings.get(h_rep, [])

                if h_rep_couplings:
                    terpene_J_vals_out = []
                    terpene_partners_out = []
                    terpene_homo_pairs: List[Tuple[int, float, str]] = []

                    for partner_h, J_val, source in h_rep_couplings:
                        if J_val >= min3:
                            terpene_J_vals_out.append(J_val)
                            terpene_partners_out.append(partner_h)
                            terpene_homo_pairs.append((partner_h, J_val, source))

                    if terpene_J_vals_out:

                        if geminal_partner is not None:
                            Jg = _get_geminal_2J(heavy_atom, base)
                            if Jg >= min2:
                                terpene_partners_out.append(geminal_partner)
                                terpene_homo_pairs.append((geminal_partner, Jg, "geminal_calc"))
                                if proton_shifts_ppm is not None:
                                    try:
                                        d_rep = proton_shifts_ppm.get(h_rep, None)
                                        d_gem = proton_shifts_ppm.get(geminal_partner, None)
                                        if d_rep and d_gem and abs(d_rep - d_gem) > 0.005:
                                            terpene_J_vals_out.append(Jg)
                                    except Exception:
                                        pass

                        return (
                            terpene_J_vals_out,
                            False,
                            terpene_partners_out,
                            [],
                            terpene_homo_pairs,
                            [],
                            [],
                            "first_order",
                            terpene_type,
                            0.0,
                            {}
                        )
        except Exception:
            pass

        try:
            macrocycle_result = estimate_macrocycle_couplings_differentiated(
                mol=molH,
                heavy_idx=heavy_idx,
                h_list=all_h_on_heavy,
                heavy_to_H=heavy_to_H,
                proton_to_heavy=proton_to_heavy_ext,
            )

            if macrocycle_result is not None:
                macro_type = macrocycle_result.get('system_type', 'macrocycle')
                macro_couplings = macrocycle_result.get('couplings', {})

                h_rep_couplings = macro_couplings.get(h_rep, [])

                if h_rep_couplings:
                    macro_J_vals_out = []
                    macro_partners_out = []
                    macro_homo_pairs: List[Tuple[int, float, str]] = []

                    for partner_h, J_val, source in h_rep_couplings:
                        if J_val >= min3:
                            macro_J_vals_out.append(J_val)
                            macro_partners_out.append(partner_h)
                            macro_homo_pairs.append((partner_h, J_val, source))

                    if macro_J_vals_out:

                        if geminal_partner is not None:
                            Jg = _get_geminal_2J(heavy_atom, base)
                            if Jg >= min2:
                                macro_partners_out.append(geminal_partner)
                                macro_homo_pairs.append((geminal_partner, Jg, "geminal_calc"))
                                if proton_shifts_ppm is not None:
                                    try:
                                        d_rep = proton_shifts_ppm.get(h_rep, None)
                                        d_gem = proton_shifts_ppm.get(geminal_partner, None)
                                        if d_rep and d_gem and abs(d_rep - d_gem) > 0.005:
                                            macro_J_vals_out.append(Jg)
                                    except Exception:
                                        pass

                        return (
                            macro_J_vals_out,
                            False,
                            macro_partners_out,
                            [],
                            macro_homo_pairs,
                            [],
                            [],
                            "first_order",
                            macro_type,
                            0.0,
                            {}
                        )
        except Exception:
            pass

        try:
            if len(all_h_on_heavy) == 2 and heavy_atom.IsInRing():
                ri = molH.GetRingInfo()

                in_ring_6 = False
                for ring in ri.AtomRings():
                    if heavy_idx in ring and len(ring) == 6:
                        in_ring_6 = True
                        break

                if in_ring_6:

                    J_large_default = 11.0
                    J_small_default = 4.0

                    generic_J_vals_out = []
                    generic_partners_out = []
                    generic_homo_pairs: List[Tuple[int, float, str]] = []

                    for neighbor in heavy_atom.GetNeighbors():
                        if neighbor.GetSymbol() != 'C':
                            continue
                        neighbor_idx = neighbor.GetIdx()

                        in_same_ring_6 = False
                        for ring in ri.AtomRings():
                            if heavy_idx in ring and neighbor_idx in ring and len(ring) == 6:
                                in_same_ring_6 = True
                                break

                        if not in_same_ring_6:
                            continue

                        partner_hs = heavy_to_H.get(neighbor_idx, [])
                        if not partner_hs:
                            continue

                        differentiate = should_differentiate_diastereotopic_J(
                            molH, heavy_idx, neighbor_idx, 'generic_ring6'
                        )

                        for partner_h in partner_hs:
                            if differentiate:

                                if h_rep == all_h_on_heavy[0]:
                                    J = J_large_default
                                else:
                                    J = J_small_default
                            else:
                                J = 7.0

                            if J >= min3:
                                generic_J_vals_out.append(J)
                                generic_partners_out.append(partner_h)
                                generic_homo_pairs.append((partner_h, J, 'ring6_generic'))

                    if generic_J_vals_out:

                        if geminal_partner is not None:
                            Jg = _get_geminal_2J(heavy_atom, base)
                            if Jg >= min2:
                                generic_partners_out.append(geminal_partner)
                                generic_homo_pairs.append((geminal_partner, Jg, "geminal_calc"))
                                if proton_shifts_ppm is not None:
                                    try:
                                        d_rep = proton_shifts_ppm.get(h_rep, None)
                                        d_gem = proton_shifts_ppm.get(geminal_partner, None)
                                        if d_rep and d_gem and abs(d_rep - d_gem) > 0.005:
                                            generic_J_vals_out.append(Jg)
                                    except Exception:
                                        pass

                        return (
                            generic_J_vals_out,
                            False,
                            generic_partners_out,
                            [],
                            generic_homo_pairs,
                            [],
                            [],
                            "first_order",
                            "ring6_generic",
                            0.0,
                            {}
                        )
        except Exception:
            pass

    if COUPLING_TABLES_MODULE_AVAILABLE:
        try:

            mol_types = MoleculeClassifier.classify(base)

            if mol_types:

                n_heavy = base.GetNumAtoms()
                candidate_heavies_spec = [j for j in range(n_heavy) if 1 <= int(D[heavy_idx, j]) <= 3]

                for other_heavy in candidate_heavies_spec:
                    other_hs = heavy_to_H.get(other_heavy, [])
                    if not other_hs:
                        continue

                    exch_spec = _is_exchangeable_H_on_X(base, other_heavy)
                    if exch_spec and solvent.get("name") != "DMSO-d6":
                        continue

                    dist = int(D[heavy_idx, other_heavy])

                    for h2 in other_hs:
                        if h2 in h_atoms:
                            continue

                        if h2 in processed_partners:
                            continue

                        J_spec = get_specialized_coupling(
                            mol=molH,
                            h1_idx=h_rep,
                            h2_idx=h2,
                            distance=dist + 1
                        )

                        if J_spec is not None and J_spec >= min3:
                            J_vals.append(J_spec)
                            J_sources.append(other_heavy)
                            J_orders.append(dist + 1)
                            partners.append(h2)
                            homo_J_vals_for_second_order.append(J_spec)
                            homo_pairs.append((h2, J_spec, "structural"))
                            processed_partners.add(h2)

                            partner_atom = base.GetAtomWithIdx(other_heavy)
                            partner_heavy_idx_kept.append(other_heavy)
                            partner_is_sp2_kept.append(
                                partner_atom.GetIsAromatic() or _is_sp2_like(partner_atom)
                            )

        except Exception as e:

            pass

    if PEPTIDE_MODULE_AVAILABLE:

        try:

            is_peptide_env = False

            for nb in heavy_atom.GetNeighbors():
                if nb.GetSymbol() == 'N':

                    for nn in nb.GetNeighbors():
                        if nn.GetSymbol() == 'C':
                            for nnn in nn.GetNeighbors():
                                if nnn.GetSymbol() == 'O':
                                    bond = base.GetBondBetweenAtoms(nn.GetIdx(), nnn.GetIdx())
                                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                                        is_peptide_env = True
                                        break

            if not is_peptide_env and heavy_atom.GetSymbol() == 'C':
                has_amino = False
                has_carboxyl = False

                for nb in heavy_atom.GetNeighbors():

                    if nb.GetSymbol() == 'N':

                        if not nb.GetIsAromatic():
                            has_amino = True

                    if nb.GetSymbol() == 'C':
                        has_carbonyl_O = False
                        has_hydroxyl_O = False
                        for nnn in nb.GetNeighbors():
                            if nnn.GetSymbol() == 'O':
                                bond = base.GetBondBetweenAtoms(nb.GetIdx(), nnn.GetIdx())
                                if bond:
                                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                                        has_carbonyl_O = True
                                    elif bond.GetBondType() == Chem.BondType.SINGLE:
                                        has_hydroxyl_O = True
                        if has_carbonyl_O and has_hydroxyl_O:
                            has_carboxyl = True

                if has_amino and has_carboxyl:
                    is_peptide_env = True

            if not is_peptide_env and heavy_atom.GetSymbol() == 'C':

                for nb in heavy_atom.GetNeighbors():
                    if nb.GetSymbol() == 'N' and not nb.GetIsAromatic():

                        is_peptide_env = True
                        break
                    elif nb.GetSymbol() == 'C':

                        visited = {heavy_idx, nb.GetIdx()}
                        queue = [(nb, 1)]
                        while queue and not is_peptide_env:
                            current, depth = queue.pop(0)
                            if depth > 4:
                                continue
                            for nn in current.GetNeighbors():
                                if nn.GetIdx() in visited:
                                    continue
                                visited.add(nn.GetIdx())
                                if nn.GetSymbol() == 'N' and not nn.GetIsAromatic():

                                    is_peptide_env = True
                                    break
                                elif nn.GetSymbol() == 'C':
                                    queue.append((nn, depth + 1))

            use_peptide_tables = is_peptide_env and (conformer_cache is None)

            if use_peptide_tables:
                n_heavy = base.GetNumAtoms()
                candidate_heavies_pep = [j for j in range(n_heavy) if 1 <= int(D[heavy_idx, j]) <= 3]

                for other_heavy in candidate_heavies_pep:
                    other_hs = heavy_to_H.get(other_heavy, [])
                    if not other_hs:
                        continue

                    exch_pep = _is_exchangeable_H_on_X(base, other_heavy)
                    if exch_pep and solvent.get("name") != "DMSO-d6":
                        continue

                    for h2 in other_hs:
                        if h2 in h_atoms:
                            continue

                        if h2 in processed_partners:
                            continue

                        try:
                            path = Chem.GetShortestPath(molH, h_rep, h2)
                            n_bonds = len(path) - 1 if path else 99
                        except Exception:
                            n_bonds = 99

                        if n_bonds > 3:

                            continue

                        J_pep = None
                        j_source_pep = "peptide"

                        if ROTAMER_AVERAGING_AVAILABLE and n_bonds == 3:

                            current_h_count = len(h_atoms)
                            partner_h_count = len(other_hs)

                            if current_h_count == 1 and partner_h_count == 2:

                                is_hb2 = (h2 == other_hs[0])
                                rotamer_result = _get_rotamer_averaged_J_for_peptide(
                                    molH, h_rep, h2, other_heavy, is_hb2=is_hb2
                                )
                                if rotamer_result:
                                    J_pep, j_source_pep = rotamer_result

                            elif current_h_count == 2 and partner_h_count == 1:

                                is_hb2 = (h_rep == h_atoms[0])
                                rotamer_result = _get_rotamer_averaged_J_for_peptide(
                                    molH, h2, h_rep, heavy_idx, is_hb2=is_hb2
                                )
                                if rotamer_result:
                                    J_pep, j_source_pep = rotamer_result

                        if J_pep is None:
                            J_pep = get_peptide_coupling_robust(
                                mol=molH,
                                h1_idx=h_rep,
                                h2_idx=h2,
                                phi=None,
                                secondary_structure=None
                            )
                            j_source_pep = "peptide"

                        if J_pep is not None and J_pep >= min3:
                            J_vals.append(J_pep)
                            J_sources.append(other_heavy)
                            J_orders.append(3)
                            partners.append(h2)
                            homo_J_vals_for_second_order.append(J_pep)
                            homo_pairs.append((h2, J_pep, j_source_pep))
                            processed_partners.add(h2)

                            partner_atom = base.GetAtomWithIdx(other_heavy)
                            partner_heavy_idx_kept.append(other_heavy)
                            partner_is_sp2_kept.append(
                                partner_atom.GetIsAromatic() or _is_sp2_like(partner_atom)
                            )

        except Exception as e:
            pass

    if COUPLING_EXTENSIONS_AVAILABLE:
        try:
            n_heavy = base.GetNumAtoms()
            candidate_heavies_ext = [j for j in range(n_heavy) if 1 <= int(D[heavy_idx, j]) <= 3]

            for other_heavy in candidate_heavies_ext:

                if config.collapse_equivalent_by_rank and ranks is not None:
                    if _are_protons_on_equivalent_carbons(base, ranks, heavy_idx, other_heavy):
                        continue

                other_hs = heavy_to_H.get(other_heavy, [])
                if not other_hs:
                    continue

                exch_ext = _is_exchangeable_H_on_X(base, other_heavy)
                if exch_ext and solvent.get("name") != "DMSO-d6":
                    continue

                dist = int(D[heavy_idx, other_heavy])

                for h2 in other_hs:
                    if h2 in h_atoms:
                        continue

                    if h2 in processed_partners:
                        continue

                    J_ext = get_extended_coupling(
                        mol=molH,
                        h1_idx=h_rep,
                        h2_idx=h2,
                        distance=dist + 1
                    )

                    if J_ext is not None and abs(J_ext) >= min3:
                        J_vals.append(abs(J_ext))
                        J_sources.append(other_heavy)
                        J_orders.append(dist + 1)
                        partners.append(h2)
                        homo_J_vals_for_second_order.append(abs(J_ext))
                        homo_pairs.append((h2, abs(J_ext), "coupling_extensions"))
                        processed_partners.add(h2)

                        partner_atom = base.GetAtomWithIdx(other_heavy)
                        partner_heavy_idx_kept.append(other_heavy)
                        partner_is_sp2_kept.append(
                            partner_atom.GetIsAromatic() or _is_sp2_like(partner_atom)
                        )

        except Exception as e:
            pass

    if geminal_partner is not None:

        Jg = _get_geminal_2J(heavy_atom, base)
        if Jg >= min2:

            geminal_shifts_differ = False
            if proton_shifts_ppm is not None:
                try:
                    d_rep = proton_shifts_ppm.get(h_rep, None)
                    d_gem = proton_shifts_ppm.get(geminal_partner, None)
                    if d_rep is not None and d_gem is not None:
                        if math.isfinite(d_rep) and math.isfinite(d_gem):
                            geminal_shifts_differ = abs(d_rep - d_gem) > 0.005
                except Exception:
                    pass

            partners.append(geminal_partner)
            homo_pairs.append((geminal_partner, Jg, "geminal_calc"))
            partner_heavy_idx_kept.append(heavy_idx)
            partner_is_sp2_kept.append(is_sp2)

            if geminal_shifts_differ:
                J_vals.append(Jg)
                J_sources.append(None)
                J_orders.append(2)
                homo_J_vals_for_second_order.append(Jg)

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
                    non_equiv_by_shift = (abs(d0 - d1) > 0.005)
            except Exception:
                non_equiv_by_shift = False

        prochiral_ch2 = False
        try:
            prochiral_ch2, _ = _is_prochiral_CH2(molH, heavy_idx)
        except Exception:
            pass
        if chiral_neigh or non_equiv_by_shift or prochiral_ch2:

            Jg = _get_geminal_2J(heavy_atom, base)
            if Jg >= min2:

                partners.append(h_atoms[1])
                homo_pairs.append((h_atoms[1], Jg, "geminal_calc"))
                partner_heavy_idx_kept.append(heavy_idx)
                partner_is_sp2_kept.append(is_sp2)

                if non_equiv_by_shift:
                    J_vals.append(Jg)
                    J_sources.append(None)
                    J_orders.append(2)
                    homo_J_vals_for_second_order.append(Jg)

    if include_heteronuclear:
        for X in base.GetAtoms():
            xsym = X.GetSymbol()
            if xsym not in config.hetero_j_base: continue

            if xsym == "P":
                dist_P = int(D[heavy_idx, X.GetIdx()]) if D[heavy_idx, X.GetIdx()] > 0 else 99
                if dist_P <= 4:
                    J_31P = _estimate_31P_1H_coupling(base, h_rep, X.GetIdx(), dist_P)
                    if J_31P is not None and J_31P >= 0.5:
                        J_vals.append(J_31P); J_sources.append(X.GetIdx()); J_orders.append(dist_P + 1)
                        continue

            for j in _hetero_J_contextual(base, D, aromatic_dmap, heavy_idx, X, config):
                if j >= 0.5:
                    J_vals.append(j); J_sources.append(None); J_orders.append(3)

    if config.peptide_HNCH_coupling:

        for nb in heavy_atom.GetNeighbors():
            if nb.GetSymbol() == 'N':
                n_idx = nb.GetIdx()

                for n_neighbor in nb.GetNeighbors():
                    if n_neighbor.GetSymbol() == 'H':
                        hn_idx = n_neighbor.GetIdx()

                        for h in h_atoms:
                            try:
                                J_hnch = _estimate_HNCH_J(molH, hn_idx, n_idx, heavy_idx, h)
                                if J_hnch is not None and J_hnch >= min3:
                                    J_vals.append(J_hnch)
                                    J_sources.append(n_idx)
                                    J_orders.append(3)
                                    partners.append(hn_idx)
                                    homo_J_vals_for_second_order.append(J_hnch)
                                    homo_pairs.append((hn_idx, J_hnch, "peptide_HNCH"))
                            except Exception:
                                pass

    n_heavy = base.GetNumAtoms()

    is_saturated_heterocycle = False
    if heavy_atom.GetSymbol() in ['C', 'N'] and not heavy_atom.GetIsAromatic():

        for ring in base.GetRingInfo().AtomRings():
            if heavy_idx in ring:

                has_hetero = any(base.GetAtomWithIdx(i).GetSymbol() in ['N', 'O', 'S'] for i in ring)
                if has_hetero and len(ring) in [5, 6]:

                    all_saturated = all(not base.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
                    if all_saturated:
                        is_saturated_heterocycle = True
                        break

    if is_saturated_heterocycle:

        candidate_heavies = [j for j in range(n_heavy) if 1 <= int(D[heavy_idx, j]) <= 2]
    else:

        candidate_heavies = [j for j in range(n_heavy) if 1 <= int(D[heavy_idx, j]) <= 3]

    anomeric_rule_applied = False

    for other_heavy in candidate_heavies:
        if suppressed_heavy and other_heavy in suppressed_heavy: continue
        dist = int(D[heavy_idx, other_heavy])
        other_hs = heavy_to_H.get(other_heavy, [])
        if not other_hs: continue

        exch_partner_global = _is_exchangeable_H_on_X(base, other_heavy)
        if exch_partner_global and solvent.get("name") != "DMSO-d6":
            continue

        if config.collapse_equivalent_by_rank and ranks is not None:

            if _are_protons_on_equivalent_carbons(base, ranks, heavy_idx, other_heavy):
                continue

            heavy_atom_check = base.GetAtomWithIdx(heavy_idx)
            other_heavy_atom = base.GetAtomWithIdx(other_heavy)
            both_aromatic = heavy_atom_check.GetIsAromatic() and other_heavy_atom.GetIsAromatic()

            if both_aromatic:
                if _are_heavy_equivalent(ranks, heavy_idx, other_heavy, heavy_to_H, proton_shifts_ppm, config.equiv_delta_ppm_tol, mol=base):
                    continue

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
            Js_tmp_sources: List[str] = []
            partners_tmp: List[int] = []
            partners_tmp_heavy: List[int] = []
            partners_tmp_is_sp2: List[bool] = []

            Js_all_for_CH3: List[float] = []
            Js_all_for_CH3_sources: List[str] = []
            partners_all_for_CH3: List[int] = []

            forbid_avg = _neighbor_is_diastereotopic_CH2(base, molH, heavy_to_H, other_heavy)
            if geminal_partner is not None:

                forbid_avg = True

            is_CH3_partner = (partner_atom.GetSymbol() == "C" and len(other_hs) == 3 and not partner_is_sp2)

            if (is_sp2 and partner_is_sp2
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
                        partners.append(ptmp); homo_pairs.append((ptmp, jv, "aromatic_exact"))
                        partner_heavy_idx_kept.append(phv); partner_is_sp2_kept.append(pis2)
                        processed_partners.add(ptmp)
                    continue

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
                        partners.append(ptmp); homo_pairs.append((ptmp, jv, "vinylic"))
                        partner_heavy_idx_kept.append(phv); partner_is_sp2_kept.append(pis2)
                    continue

            for h in other_hs:
                if h in h_atoms:
                    continue

                if h in processed_partners:
                    continue

                J = None
                j_source = "unknown"
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

                    is_allylic_coupling = (
                        flexible_cc
                        and partner_is_sp2
                        and len(heavy_to_H.get(other_heavy, [])) == 1
                    )

                    if is_allylic_coupling:

                        J = 7.0
                        j_source = "fallback_average"

                    is_diastereotopic_CH2 = False
                    if ultra_fast_mode:

                        is_diastereotopic_CH2 = False
                    elif flexible_cc:

                        current_h_count = len(heavy_to_H.get(heavy_idx, []))
                        partner_h_count = len(heavy_to_H.get(other_heavy, []))

                        if current_h_count == 2:
                            if diastereotopic_cache and heavy_idx in diastereotopic_cache:
                                is_diast, _ = diastereotopic_cache[heavy_idx]
                                is_diastereotopic_CH2 = is_diast
                            elif not ultra_fast_mode:

                                is_diast, _ = _detect_diastereotopic_protons(molH, heavy_idx)
                                if is_diast:
                                    is_diastereotopic_CH2 = True
                                else:
                                    is_prochiral, _ = _is_prochiral_CH2(molH, heavy_idx)
                                    if is_prochiral:
                                        is_diastereotopic_CH2 = True

                        if not is_diastereotopic_CH2 and partner_h_count == 2:
                            if diastereotopic_cache and other_heavy in diastereotopic_cache:
                                is_diast_partner, _ = diastereotopic_cache[other_heavy]
                                is_diastereotopic_CH2 = is_diast_partner
                            elif not ultra_fast_mode:

                                is_diast_partner, _ = _detect_diastereotopic_protons(molH, other_heavy)
                                if is_diast_partner:
                                    is_diastereotopic_CH2 = True
                                else:
                                    is_prochiral_partner, _ = _is_prochiral_CH2(molH, other_heavy)
                                    if is_prochiral_partner:
                                        is_diastereotopic_CH2 = True

                    if is_allylic_coupling:
                        J = 7.0
                        j_source = "fallback_average"

                    elif ultra_fast_mode and flexible_cc:
                        J = config.ultra_fast_default_vicinal_J
                        j_source = "fallback_average"

                    elif flexible_cc:

                        amino_acid_J = None
                        if is_diastereotopic_CH2 and ROTAMER_AVERAGING_AVAILABLE:

                            current_h_count = len(heavy_to_H.get(heavy_idx, []))
                            partner_h_count = len(heavy_to_H.get(other_heavy, []))

                            if current_h_count == 1 and partner_h_count == 2:
                                ch2_idx = other_heavy
                                ca_idx = heavy_idx
                                is_current_h_hb2 = (h == heavy_to_H.get(other_heavy, [])[0])

                            elif current_h_count == 2 and partner_h_count == 1:
                                ch2_idx = heavy_idx
                                ca_idx = other_heavy
                                is_current_h_hb2 = (h_rep == heavy_to_H.get(heavy_idx, [])[0])
                            else:
                                ch2_idx = None
                                ca_idx = None
                                is_current_h_hb2 = True

                            if ch2_idx is not None:

                                aa_type = _detect_amino_acid_type(molH, ch2_idx)
                                if aa_type:
                                    try:

                                        if ca_idx is not None and _is_free_amino_acid(molH, ca_idx):
                                            secondary_structure = 'free_amino_acid'
                                        else:
                                            secondary_structure = 'random_coil'

                                        J_HB2, J_HB3, aa_source = get_diastereotopic_J_pair(aa_type, secondary_structure)
                                        if is_current_h_hb2:
                                            amino_acid_J = J_HB2
                                        else:
                                            amino_acid_J = J_HB3

                                        j_source = f"rotamer_pair:{aa_type}:{secondary_structure}"
                                    except Exception:
                                        pass

                        if amino_acid_J is not None:
                            J = amino_acid_J
                        else:

                            J, j_source = _get_flexible_chain_fallback_J(base, heavy_idx, other_heavy)
                    else:
                        EN_corr = 0.4 if any(n.GetSymbol() in {"O","N"} for n in base.GetAtomWithIdx(heavy_idx).GetNeighbors()) else 0.0

                        Jk = None

                        macrocycle_J = None
                        if base.GetAtomWithIdx(heavy_idx).IsInRing() and base.GetAtomWithIdx(other_heavy).IsInRing():
                            ring_info = base.GetRingInfo()
                            for ring in ring_info.AtomRings():
                                ring_size = len(ring)
                                if ring_size >= 8 and heavy_idx in ring and other_heavy in ring:

                                    try:
                                        macrocycle_J = get_macrocycle_coupling(
                                            base, h_rep, h, ring_size,
                                            macrocycle_type='macrolactam' if _has_amide_in_ring(base, ring) else 'macrolactone'
                                        )
                                    except Exception:
                                        macrocycle_J = 7.0
                                    break

                        if macrocycle_J is not None:
                            Jk = macrocycle_J
                            j_source = "macrocycle_table"

                        steroid_J = None
                        if Jk is None and COUPLING_TABLES_MODULE_AVAILABLE:
                            if base.GetAtomWithIdx(heavy_idx).IsInRing() and base.GetAtomWithIdx(other_heavy).IsInRing():
                                ring_info = base.GetRingInfo()
                                n_fused_5_6 = sum(1 for r in ring_info.AtomRings() if len(r) in (5, 6))
                                n_stereocenters = sum(1 for a in base.GetAtoms()
                                                     if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED)

                                if n_fused_5_6 >= 4 or (n_fused_5_6 >= 3 and n_stereocenters >= 4):
                                    try:
                                        steroid_J = get_steroid_coupling(base, h_rep, h)
                                    except Exception:
                                        pass

                        if steroid_J is not None:
                            Jk = steroid_J
                            j_source = "steroid_table"

                        terpene_J = None
                        if Jk is None and COUPLING_TABLES_MODULE_AVAILABLE:
                            if base.GetAtomWithIdx(heavy_idx).IsInRing() and base.GetAtomWithIdx(other_heavy).IsInRing():

                                ring_info = base.GetRingInfo()
                                has_gem_dimethyl = False
                                for atom in base.GetAtoms():
                                    if atom.GetSymbol() == 'C' and not atom.IsInRing():

                                        heavy_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != 'H']
                                        if len(heavy_neighbors) == 1 and atom.GetTotalNumHs() == 3:
                                            parent = heavy_neighbors[0]

                                            ch3_count = sum(1 for n in parent.GetNeighbors()
                                                          if n.GetSymbol() == 'C' and not n.IsInRing()
                                                          and n.GetTotalNumHs() == 3)
                                            if ch3_count >= 2 and parent.IsInRing():
                                                has_gem_dimethyl = True
                                                break

                                if has_gem_dimethyl:

                                    ring_size = None
                                    for ring in ring_info.AtomRings():
                                        if heavy_idx in ring and other_heavy in ring:
                                            ring_size = len(ring)
                                            break
                                    if ring_size == 6:
                                        try:
                                            terpene_J = get_terpene_coupling(base, h_rep, h, 'menthane')
                                        except Exception:
                                            pass
                                    elif ring_size in (4, 5):
                                        try:
                                            terpene_J = get_terpene_coupling(base, h_rep, h, 'pinane')
                                        except Exception:
                                            pass

                        if terpene_J is not None:
                            Jk = terpene_J
                            j_source = "terpene_table"

                        elif conformer_cache is not None and config.use_conformational_averaging:
                            Jk = conformer_cache.get_J_vicinal_averaged(
                                h_rep, h,
                                use_haasnoot=(config.karplus_haasnoot and base.GetAtomWithIdx(heavy_idx).IsInRingSize(6)),
                                use_haasnoot_7_params=config.use_haasnoot_7_params,
                                electroneg_corr=EN_corr
                            )
                            j_source = "karplus_3d"

                            n_conf_used = conformer_cache.get_num_conformers() if hasattr(conformer_cache, 'get_num_conformers') else 1
                            if flexible_cc and is_diastereotopic_CH2 and Jk is not None and Jk < 1.5 and n_conf_used < 10:

                                Jk_direct = _estimate_vicinal_sp3_J(
                                    molH, h_rep, h,
                                    use_haasnoot=False,
                                    electroneg_corr=EN_corr,
                                    use_full_haasnoot=False
                                )
                                if Jk_direct is not None:
                                    Jk = Jk_direct

                            current_h_count = len(heavy_to_H.get(heavy_idx, []))
                            partner_h_count = len(heavy_to_H.get(other_heavy, []))
                            is_true_CH2_CH2 = (current_h_count == 2 and partner_h_count == 2)

                            if flexible_cc and is_diastereotopic_CH2 and is_true_CH2_CH2 and Jk is not None:

                                if Jk < 3.0:

                                    Jk = max(Jk, 3.0)
                                elif Jk > 9.0:

                                    Jk = 6.5
                        else:
                            Jk = _estimate_vicinal_sp3_J(
                                molH, h_rep, h,
                                use_haasnoot=(config.karplus_haasnoot and base.GetAtomWithIdx(heavy_idx).IsInRingSize(6)),
                                electroneg_corr=EN_corr,
                                use_full_haasnoot=config.use_haasnoot_7_params
                            )
                            j_source = "karplus_3d"
                        if Jk is None:
                            J = 7.0 if is_carbon_partner else None
                            j_source = "fallback_average"
                        else:
                            J = Jk
                else:
                    rel = _aromatic_rel_from_map(aromatic_dmap, heavy_idx, other_heavy)
                    if rel == "ortho":
                        J = config.aro_j_forced["ortho"]
                        j_source = "aromatic_exact"
                    elif rel == "meta":
                        J = config.aro_j_forced["meta"]
                        j_source = "aromatic_exact"
                    elif rel == "para":
                        J = config.aro_j_forced["para"]
                        j_source = "aromatic_exact"
                    else:
                        Jv = _estimate_vinylic_J(molH, h, h_rep, config=config)
                        if Jv is None:
                            J = _estimate_vinylic_J(molH, h_rep, h, config=config)
                        else:
                            J = Jv
                        if J is not None:
                            j_source = "vinylic"
                        if J is None:

                            J_alkyne = _estimate_alkyne_J(molH, h_rep, h)
                            if J_alkyne is not None:
                                J = J_alkyne
                                j_source = "alkyne"
                            else:

                                if not partner_is_sp2:

                                    bond = base.GetBondBetweenAtoms(heavy_idx, other_heavy)
                                    if (bond is not None
                                        and bond.GetBondTypeAsDouble() == 1.0
                                        and not base.GetAtomWithIdx(other_heavy).IsInRing()):

                                        J = 7.0
                                        j_source = "fallback_average"
                                    else:
                                        J = config.vinylic_J_defaults.get('cis', 10.0)
                                        j_source = "vinylic_default"
                                else:
                                    J = config.vinylic_J_defaults.get('cis', 10.0)
                                    j_source = "vinylic_default"

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

                    A = base.GetAtomWithIdx(heavy_idx)
                    is_formyl_this = (A.GetSymbol()=="C"
                        and any(b.GetBondTypeAsDouble()==2 and b.GetOtherAtom(A).GetSymbol()=="O" for b in A.GetBonds())
                        and len(heavy_to_H.get(heavy_idx, [])) == 1)

                    B = base.GetAtomWithIdx(other_heavy)
                    is_formyl_partner = (B.GetSymbol()=="C"
                        and any(b.GetBondTypeAsDouble()==2 and b.GetOtherAtom(B).GetSymbol()=="O" for b in B.GetBonds())
                        and len(heavy_to_H.get(other_heavy, [])) == 1)

                    if is_formyl_this or is_formyl_partner:

                        if J is None: J = 2.5
                        J = max(1.5, min(3.0, J))

                if config.amide_mode:
                    def _is_amide_like(idx_c):
                        C = base.GetAtomWithIdx(idx_c)
                        return (C.GetSymbol()=="C" and any(b.GetBondTypeAsDouble()==2 and b.GetOtherAtom(C).GetSymbol()=="O" for b in C.GetBonds()))

                    def _is_formyl(idx_c):

                        C = base.GetAtomWithIdx(idx_c)
                        return (C.GetSymbol()=="C"
                            and any(b.GetBondTypeAsDouble()==2 and b.GetOtherAtom(C).GetSymbol()=="O" for b in C.GetBonds())
                            and len(heavy_to_H.get(idx_c, [])) == 1)

                    if (_is_amide_like(heavy_idx) or _is_amide_like(other_heavy)) and not (_is_formyl(heavy_idx) or _is_formyl(other_heavy)):
                        if J is None: J = 6.0
                        J = max(3.5, min(9.0, J))

                if is_CH3_partner and J is not None:
                    Js_all_for_CH3.append(J)
                    Js_all_for_CH3_sources.append(j_source)
                    partners_all_for_CH3.append(h)

                if J is not None and J >= min3:
                    Js_tmp.append(J); Js_tmp_sources.append(j_source)
                    partners_tmp.append(h)
                    partners_tmp_heavy.append(other_heavy); partners_tmp_is_sp2.append(partner_is_sp2)

            if is_CH3_partner and not forbid_avg and len(Js_all_for_CH3) == 3:
                Jmean = float(np.mean(Js_all_for_CH3))

                ch3_source = Js_all_for_CH3_sources[0] if Js_all_for_CH3_sources else "karplus_3d"
                if Jmean >= min3:

                    for h_ch3 in other_hs:
                        if h_ch3 in h_atoms:
                            continue
                        if h_ch3 in processed_partners:
                            continue
                        J_vals.append(Jmean); J_sources.append(other_heavy); J_orders.append(3)
                        homo_J_vals_for_second_order.append(Jmean)
                        partners.append(h_ch3); homo_pairs.append((h_ch3, Jmean, ch3_source))
                        partner_heavy_idx_kept.append(other_heavy); partner_is_sp2_kept.append(partner_is_sp2)
                        processed_partners.add(h_ch3)
                continue

            if Js_tmp:
                dyn_tol = max(TOLERANCE_J_HZ, config.grouping_tol_scale * config.line_width_hz)

                close_enough = (max(Js_tmp) - min(Js_tmp) <= 0.5 * dyn_tol)

                partners_are_equivalent = False
                if len(partners_tmp) >= 2 and proton_shifts_ppm is not None:
                    partner_shifts = [proton_shifts_ppm.get(p, None) for p in partners_tmp]
                    if all(s is not None for s in partner_shifts):
                        shift_spread = max(partner_shifts) - min(partner_shifts)

                        partners_are_equivalent = (shift_spread <= 0.02)

                should_average = (len(Js_tmp) >= 2) and (not is_sp2) and (not partner_is_sp2)
                if partners_are_equivalent:

                    should_average = should_average and True
                else:

                    should_average = should_average and (not forbid_avg) and close_enough

                if should_average:
                    Jmean = float(np.mean(Js_tmp))

                    avg_source = Js_tmp_sources[0] if Js_tmp_sources else "karplus_3d"
                    for ptmp, phv, pis2 in zip(partners_tmp, partners_tmp_heavy, partners_tmp_is_sp2):
                        if ptmp in processed_partners:
                            continue
                        J_vals.append(Jmean); J_sources.append(other_heavy); J_orders.append(3)
                        homo_J_vals_for_second_order.append(Jmean)
                        partners.append(ptmp); homo_pairs.append((ptmp, Jmean, avg_source))
                        partner_heavy_idx_kept.append(phv); partner_is_sp2_kept.append(pis2)
                        processed_partners.add(ptmp)
                else:
                    for jv, jsrc, ptmp, phv, pis2 in zip(Js_tmp, Js_tmp_sources, partners_tmp, partners_tmp_heavy, partners_tmp_is_sp2):
                        if ptmp in processed_partners:
                            continue
                        J_vals.append(jv); J_sources.append(other_heavy); J_orders.append(3)
                        homo_J_vals_for_second_order.append(jv)
                        partners.append(ptmp); homo_pairs.append((ptmp, jv, jsrc))
                        partner_heavy_idx_kept.append(phv); partner_is_sp2_kept.append(pis2)
                        processed_partners.add(ptmp)
            continue

        if dist == 2:

            if config.collapse_equivalent_by_rank and ranks is not None:
                if _are_protons_on_equivalent_carbons(base, ranks, heavy_idx, other_heavy):
                    continue

            if _is_isolated_methyl(base, heavy_idx, heavy_to_H):
                continue

            if _is_NCH3_coupling_to_formamide_CHO(base, heavy_idx, other_heavy, heavy_to_H):
                continue

            if _is_aldehyde_4J_coupling(base, heavy_idx, other_heavy, heavy_to_H, D):
                continue

            if _is_benzylic_CH2_4J_coupling(base, heavy_idx, other_heavy, heavy_to_H):
                continue

            Jctx, rel = _aromatic_J_contextual(base, aromatic_dmap, heavy_idx, other_heavy,
                                               J0_meta=config.aro_j_forced["meta"], J0_para=config.aro_j_forced["para"])
            if rel == "meta":

                try:
                    symmetry_classes = list(Chem.CanonicalRankAtoms(base, breakTies=False))
                    carbons_equivalent = (symmetry_classes[heavy_idx] == symmetry_classes[other_heavy])
                except Exception:
                    carbons_equivalent = False

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

                if is_pure_benzene and carbons_equivalent and Jctx is not None and Jctx < 2.0:
                    continue
                J = Jctx
            else:

                is_homoallyl, J_homoallyl = _is_homoallylic_4J(base, heavy_idx, other_heavy, mol3d=molH)
                if is_homoallyl:

                    J = J_homoallyl

                else:

                    is_allyl, J_allyl = _is_allylic_4J(base, heavy_idx, other_heavy, mol3d=molH)
                    if is_allyl:

                        if _is_geminal_methyl(base, heavy_idx, heavy_to_H):
                            continue
                        if _is_geminal_methyl(base, other_heavy, heavy_to_H):
                            continue
                        J = J_allyl
                    else:

                        isW, isRigid = _is_W_allylic(base, heavy_idx, other_heavy)
                        if isW:
                            J = 2.8 if isRigid else 1.6
                        else:
                            J = 0.3
                        if config.rigid_W_boost and isW:
                            J = min(4.0, max(J, 2.0))

            adaptive_min4 = _get_adaptive_4J_threshold(base, heavy_idx, other_heavy, min4)
            if J is not None and J >= adaptive_min4:

                if other_hs[0] in processed_partners:
                    continue
                J_vals.append(J); J_sources.append(None); J_orders.append(4)
                partners.append(other_hs[0]); homo_J_vals_for_second_order.append(J)
                homo_pairs.append((other_hs[0], J, "long_range_4J"))
                partner_heavy_idx_kept.append(other_heavy); partner_is_sp2_kept.append(_is_sp2_like(base.GetAtomWithIdx(other_heavy)))
                processed_partners.add(other_hs[0])
            continue

        if dist == 3:

            if config.collapse_equivalent_by_rank and ranks is not None:
                if _are_protons_on_equivalent_carbons(base, ranks, heavy_idx, other_heavy):
                    continue

            if _is_isolated_methyl(base, heavy_idx, heavy_to_H):
                continue

            Jctx, rel = _aromatic_J_contextual(base, aromatic_dmap, heavy_idx, other_heavy,
                                               J0_meta=config.aro_j_forced["meta"], J0_para=config.aro_j_forced["para"])
            if rel == "para":

                try:
                    symmetry_classes = list(Chem.CanonicalRankAtoms(base, breakTies=False))
                    carbons_equivalent = (symmetry_classes[heavy_idx] == symmetry_classes[other_heavy])
                except Exception:
                    carbons_equivalent = False

                if carbons_equivalent:
                    J = Jctx
                    para_buf_J.append(J); para_buf_partners.append(other_hs[0])
            else:

                is_W_5J, J_W_5J = _is_W_coupling_5J(base, heavy_idx, other_heavy)
                if is_W_5J:
                    J = J_W_5J
                else:
                    J = 0.2
                    if config.rigid_W_boost:
                        if (base.GetAtomWithIdx(heavy_idx).IsInRing() and base.GetAtomWithIdx(other_heavy).IsInRing()):
                            J = 0.7
                if J >= min5:

                    if other_hs[0] in processed_partners:
                        continue
                    J_vals.append(J); J_sources.append(None); J_orders.append(5)
                    partners.append(other_hs[0]); homo_J_vals_for_second_order.append(J)
                    homo_pairs.append((other_hs[0], J, "long_range_5J"))
                    partner_heavy_idx_kept.append(other_heavy); partner_is_sp2_kept.append(_is_sp2_like(base.GetAtomWithIdx(other_heavy)))
                    processed_partners.add(other_hs[0])
            continue

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
                homo_pairs.append((p_sel, j_sel, "aromatic_para"))
                try:
                    hh = [n.GetIdx() for n in molH.GetAtomWithIdx(p_sel).GetNeighbors() if n.GetSymbol()!='H'][0]
                    partner_heavy_idx_kept.append(hh)
                    partner_is_sp2_kept.append(_is_sp2_like(base.GetAtomWithIdx(hh)))
                except Exception:
                    pass

    visco = solvent.get("viscosity", "low")
    low_cut = config.viscosity_low_cuts.get(visco, config.viscosity_low_cuts["low"])
    J_vals = [j for j in J_vals if j >= low_cut]
    if _has_quadrupolar_neighbor(base, heavy_idx):
        J_vals = [j for j in J_vals if j >= max(low_cut, 0.5)]

    if homo_J_vals_for_second_order:
        real_ppm_vals: List[float] = []
        if proton_shifts_ppm is not None:
            try:
                d0 = proton_shifts_ppm.get(h_rep, None)
                if d0 is not None and math.isfinite(d0):
                    real_ppm_vals.append(d0)

                for p in partners:

                    if p in h_atoms:
                        continue
                    val = proton_shifts_ppm.get(p, None)
                    if val is not None and math.isfinite(val):
                        real_ppm_vals.append(val)
            except Exception:
                real_ppm_vals = []

        use_real_vals = False
        if len(real_ppm_vals) >= 2:
            unique_vals = set(round(v, 3) for v in real_ppm_vals)
            if len(unique_vals) > 1:

                ppm_vals = real_ppm_vals
                use_real_vals = True

        if not use_real_vals:

            if len(h_atoms) == 2:

                delta_ppm_default = DEFAULT_DELTA.get("ch2_near_chiral_1", 0.12) if DIASTEREOTOPIC_TABLE_AVAILABLE else 0.12
            elif is_sp2:

                delta_ppm_default = DEFAULT_DELTA.get("vinyl_geminal", 0.30) if DIASTEREOTOPIC_TABLE_AVAILABLE else 0.30
            else:

                delta_ppm_default = DEFAULT_DELTA.get("unknown", 0.10) if DIASTEREOTOPIC_TABLE_AVAILABLE else 0.10
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

    if system_type in ["AA'BB'", "AA'XX'"]:
        if not second_order:
            second_order = True
            if second_order_type == "first_order":
                second_order_type = "strong_coupling"

    if config.use_solvent_J_corrections and J_vals:
        solvent_name = solvent.get("name", "CDCl3")
        corrected_J_vals = []
        for i, J in enumerate(J_vals):
            order = J_orders[i] if i < len(J_orders) else 3

            if order == 2:
                j_type = "2J_geminal"
            elif order == 4:
                j_type = "4J_allylic"
            elif order >= 5:
                j_type = "3J_aromatic"
            elif is_sp2:
                j_type = "3J_vinyl"
            else:
                j_type = "3J_vicinal"

            J_corrected = _apply_solvent_J_correction(J, solvent_name, j_type)
            corrected_J_vals.append(J_corrected)

        J_vals = corrected_J_vals

        homo_pairs_corrected: List[Tuple[int, float, str]] = []
        for i, hp_tuple in enumerate(homo_pairs):
            partner = hp_tuple[0]
            hp_source = hp_tuple[2] if len(hp_tuple) > 2 else "unknown"
            if i < len(J_vals):
                homo_pairs_corrected.append((partner, J_vals[i], hp_source))
            else:
                homo_pairs_corrected.append((partner, 0.0, hp_source))
        homo_pairs = homo_pairs_corrected

    return (J_vals, second_order, partners, satellites, homo_pairs,
            partner_heavy_idx_kept, partner_is_sp2_kept,
            second_order_type, system_type, delta_nu_min_hz, J_orders)

def _effective_exchange_regime(solvent: dict, temperature_k: float) -> str:
    base_reg = solvent.get("exchange_regime", "moderate")
    rate = solvent.get("exchange_rate_factor", 1.0) * (temperature_k / 298.15)
    if solvent.get("name") in {"D2O"}: return "fast"
    if rate >= 1.5: return "fast"
    if rate <= 0.75: return "slow" if base_reg != "fast" else "moderate"
    return base_reg

def _infer_j_source_type(J_hz: float, is_flexible_system: bool = False) -> Tuple[str, int]:

    J_str = f"{J_hz:.6f}".rstrip('0')
    if '.' in J_str:
        n_decimals = len(J_str.split('.')[1])
    else:
        n_decimals = 0

    J_frac = abs(J_hz - round(J_hz))

    precision_score = min(1.0, n_decimals / 3.0 + J_frac)

    is_very_precise = precision_score > 0.5
    is_moderately_precise = precision_score > 0.2
    is_round = precision_score < 0.1

    is_geminal_range = 10.0 <= J_hz <= 16.0
    is_vicinal_mean_range = 5.5 <= J_hz <= 8.5
    is_gauche_range = 2.0 <= J_hz <= 5.0
    is_anti_range = 8.5 <= J_hz <= 12.0
    is_extreme_range = is_gauche_range or is_anti_range

    if is_geminal_range:
        if is_very_precise:
            return "geminal_calculated", 9
        elif is_moderately_precise:
            return "geminal_structural", 7
        else:
            return "geminal_default", 5

    if is_flexible_system:

        if is_extreme_range:

            if is_very_precise:

                return "karplus_biased", 3
            else:
                return "conformer_biased", 2

        elif is_vicinal_mean_range:

            if is_very_precise:
                return "karplus_averaged", 8
            elif is_moderately_precise:
                return "rotational_estimate", 7
            else:

                return "rotational_default", 6

    else:
        if is_very_precise:

            return "karplus_3d", 10

        elif is_moderately_precise:

            if is_extreme_range:
                return "karplus_likely", 8
            else:
                return "structural", 6

        else:

            if is_geminal_range:
                return "geminal_default", 5
            elif is_vicinal_mean_range:
                return "fallback_vicinal", 3
            else:
                return "fallback_generic", 2

    return "unknown", 1

def _post_process_reciprocity_and_partners(
    result: Dict[int, Dict[str, Any]],
    proton_shifts_ppm: Dict[int, float],
    heavy_to_H: Dict[int, List[int]],
    is_large_molecule: bool = False,
) -> Dict[int, Dict[str, Any]]:

    proton_to_heavy: Dict[int, int] = {}
    for heavy_idx, h_list in heavy_to_H.items():
        for h in h_list:
            proton_to_heavy[h] = heavy_idx

    coupling_registry: Dict[Tuple[int, int], List[Tuple[float, str]]] = {}

    source_priority = {

        "vinylic": 20,
        "alkyne": 20,
        "aromatic_exact": 19,
        "heterocycle": 19,
        "heterocycle_table": 19,
        "heterocycle_expected_mult": 19,
        "heterocycle_expected_singlet": 19,
        "heterocycle_singlet": 19,
        "saturated_heterocycle": 18,
        "dkp_pyrrolidine_fused_table": 18,
        "pyrrolidine_fused_table": 18,
        "diketopiperazine_experimental": 18,
        "steroid_table": 18,
        "terpene_table": 18,
        "macrocycle_table": 18,
        "nucleotide_table": 18,
        "coupling_tables": 18,
        "peptide_HNCH": 17,
        "sugar": 16,
        "structural": 15,
        "coupling_extensions": 18,
        "geminal_calc": 14,
        "long_range_4J": 12,
        "long_range_5J": 12,
        "fallback_ortho": 10,
        "fallback_meta": 10,

        "karplus_3d": 5,
        "aromatic_para": 4,
        "gauche_effect_double": 4,
        "gauche_effect_O": 4,
        "gauche_effect_F": 4,
        "extended_chain_cation": 4,
        "fallback_average": 3,
        "diastereotopic_estimated": 2,
        "vinylic_default": 2,
        "peptide": 2,
        "unknown": 1
    }

    equivalent_groups: Dict[int, List[int]] = {}
    for heavy_idx, h_list in heavy_to_H.items():

        if len(h_list) == 3:
            equivalent_groups[heavy_idx] = h_list

    for heavy_idx, equiv_h_list in equivalent_groups.items():
        if len(equiv_h_list) < 2:
            continue

        partner_to_J_list: Dict[int, List[Tuple[float, str]]] = {}

        for h_idx in equiv_h_list:
            if h_idx not in result:
                continue
            data = result[h_idx]
            for p in data.get("partners", []):
                partner_h = p.get("partner_h_index")
                J_hz = p.get("J_hz", 0)
                J_source = p.get("J_source", "unknown")
                if partner_h is not None and J_hz > 0:

                    if partner_h in equiv_h_list:
                        continue
                    if partner_h not in partner_to_J_list:
                        partner_to_J_list[partner_h] = []
                    partner_to_J_list[partner_h].append((J_hz, J_source))

        for partner_h, J_source_list in partner_to_J_list.items():
            if len(J_source_list) == 0:
                continue

            J_values = [js[0] for js in J_source_list]
            J_mean = sum(J_values) / len(J_values)

            best_source = max(J_source_list, key=lambda x: source_priority.get(x[1], 1))[1]

            for h_idx in equiv_h_list:
                if h_idx not in result:
                    continue
                data = result[h_idx]
                partners = data.get("partners", [])

                found = False
                for p in partners:
                    if p.get("partner_h_index") == partner_h:
                        p["J_hz"] = round(J_mean, 2)
                        p["J_source"] = best_source
                        found = True
                        break

                if not found and partner_h in proton_shifts_ppm:
                    partners.append({
                        "partner_h_index": int(partner_h),
                        "partner_shift_ppm": round(float(proton_shifts_ppm.get(partner_h, 0.0)), 3),
                        "J_hz": round(J_mean, 2),
                        "J_source": best_source,
                    })

    for h_idx, data in result.items():
        partners = data.get("partners", [])
        for p in partners:
            partner_h = p.get("partner_h_index")
            J_hz = p.get("J_hz", 0)
            J_source = p.get("J_source", "unknown")
            if partner_h is not None and J_hz > 0:

                key = (min(h_idx, partner_h), max(h_idx, partner_h))
                if key not in coupling_registry:
                    coupling_registry[key] = []
                coupling_registry[key].append((J_hz, J_source))

    resolved_couplings: Dict[Tuple[int, int], float] = {}

    for key, J_source_pairs in coupling_registry.items():
        if len(J_source_pairs) == 1:
            resolved_couplings[key] = J_source_pairs[0][0]
        elif len(J_source_pairs) >= 2:

            J_values_only = [js[0] for js in J_source_pairs]
            J_min, J_max = min(J_values_only), max(J_values_only)

            h_a, h_b = key
            heavy_a = proton_to_heavy.get(h_a)
            heavy_b = proton_to_heavy.get(h_b)
            is_diastereotopic_pair = False

            if heavy_a is not None:
                siblings_a = heavy_to_H.get(heavy_a, [])
                if len(siblings_a) == 2:
                    is_diastereotopic_pair = True
            if heavy_b is not None:
                siblings_b = heavy_to_H.get(heavy_b, [])
                if len(siblings_b) == 2:
                    is_diastereotopic_pair = True

            if is_diastereotopic_pair and abs(J_max - J_min) > 2.0:

                priorities = [(js[0], source_priority.get(js[1], 1), js[1]) for js in J_source_pairs]
                max_priority = max(p[1] for p in priorities)
                min_priority = min(p[1] for p in priorities)

                TABLE_THRESHOLD = 7
                ESTIMATED_THRESHOLD = 3

                has_table_source = any(p[1] >= TABLE_THRESHOLD for p in priorities)
                has_estimated_source = any(p[1] <= ESTIMATED_THRESHOLD for p in priorities)

                if has_table_source and has_estimated_source and (max_priority - min_priority) >= 4:

                    best_entry = max(priorities, key=lambda p: p[1])
                    resolved_couplings[key] = best_entry[0]
                else:

                    if 2.0 <= J_min <= 5.5 and 5.5 < J_max <= 8.5:

                        resolved_couplings[key] = J_min
                    elif 8.0 <= J_max <= 12.0 and J_min < 8.0:

                        resolved_couplings[key] = J_max
                    elif 2.0 <= J_min <= 5.5:

                        resolved_couplings[key] = J_min
                    elif 8.0 <= J_max <= 12.0:

                        resolved_couplings[key] = J_max
                    else:

                        non_standard = [j for j in J_values_only if j < 5.5 or j > 8.5]
                        if non_standard:
                            resolved_couplings[key] = non_standard[0]
                        else:
                            resolved_couplings[key] = J_min

            elif abs(J_max - J_min) < 1.0:
                resolved_couplings[key] = round(sum(J_values_only) / len(J_values_only), 3)

            else:

                J_source_sorted = sorted(
                    J_source_pairs,
                    key=lambda x: source_priority.get(x[1], 1),
                    reverse=True
                )

                best_J, best_source = J_source_sorted[0]
                best_priority = source_priority.get(best_source, 1)
                worst_priority = source_priority.get(J_source_sorted[-1][1], 1)

                if best_priority > worst_priority:

                    resolved_couplings[key] = best_J
                else:

                    resolved_couplings[key] = round(sum(J_values_only) / len(J_values_only), 2)

    for h_idx, data in result.items():
        partners = data.get("partners", [])
        J_couplings = data.get("J_couplings", [])

        if not partners and J_couplings:
            new_partners = []
            for key, J_resolved in resolved_couplings.items():
                h_a, h_b = key
                if h_a == h_idx:
                    partner_h = h_b
                elif h_b == h_idx:
                    partner_h = h_a
                else:
                    continue

                for jc in J_couplings:
                    if abs(jc.get("J_hz", 0) - J_resolved) < 1.0:

                        J_source_list = coupling_registry.get(key, [])
                        J_source = "unknown"
                        if J_source_list:

                            best_src = max(J_source_list, key=lambda x: source_priority.get(x[1], 1))
                            J_source = best_src[1]
                        new_partners.append({
                            "partner_h_index": int(partner_h),
                            "partner_shift_ppm": round(float(proton_shifts_ppm.get(partner_h, 0.0)), 3),
                            "J_hz": round(float(J_resolved), 2),
                            "J_source": J_source,
                        })
                        break

            if new_partners:
                data["partners"] = new_partners

        else:
            updated_partners = []
            for p in partners:
                partner_h = p.get("partner_h_index")
                if partner_h is None:
                    updated_partners.append(p)
                    continue

                key = (min(h_idx, partner_h), max(h_idx, partner_h))
                J_resolved = resolved_couplings.get(key, p.get("J_hz", 0))

                existing_source = p.get("J_source", "unknown")
                if existing_source == "unknown":
                    J_source_list = coupling_registry.get(key, [])
                    if J_source_list:
                        best_src = max(J_source_list, key=lambda x: source_priority.get(x[1], 1))
                        existing_source = best_src[1]

                updated_partners.append({
                    "partner_h_index": int(partner_h),
                    "partner_shift_ppm": round(float(proton_shifts_ppm.get(partner_h, p.get("partner_shift_ppm", 0.0))), 3),
                    "J_hz": round(float(J_resolved), 2),
                    "J_source": existing_source,
                })

            data["partners"] = updated_partners

    for heavy_idx, h_list in heavy_to_H.items():
        if len(h_list) != 2:
            continue

        h_a, h_b = h_list[0], h_list[1]

        data_a = result.get(h_a, {})
        data_b = result.get(h_b, {})

        if not data_a or not data_b:
            continue

        partners_a = data_a.get("partners", [])
        partners_b = data_b.get("partners", [])

        is_diastereotopic = False

        for p in partners_a:
            if p.get("partner_h_index") == h_b:
                J_hz = p.get("J_hz", 0)
                if 10.0 <= J_hz <= 16.0:
                    is_diastereotopic = True
                    break

        if not is_diastereotopic:
            for p in partners_b:
                if p.get("partner_h_index") == h_a:
                    J_hz = p.get("J_hz", 0)
                    if 10.0 <= J_hz <= 16.0:
                        is_diastereotopic = True
                        break

        if not is_diastereotopic:
            jc_a = data_a.get("J_couplings", [])
            jc_b = data_b.get("J_couplings", [])
            if len(jc_a) != len(jc_b):
                is_diastereotopic = True

        if not is_diastereotopic:
            continue

        partners_a = data_a.get("partners", [])
        partners_b = data_b.get("partners", [])

        vicinal_partners_a = set()
        vicinal_partners_b = set()

        for p in partners_a:
            partner_h = p.get("partner_h_index")
            J_hz = p.get("J_hz", 0)

            if partner_h and partner_h != h_b and J_hz < 11.0:
                vicinal_partners_a.add(partner_h)

        for p in partners_b:
            partner_h = p.get("partner_h_index")
            J_hz = p.get("J_hz", 0)
            if partner_h and partner_h != h_a and J_hz < 11.0:
                vicinal_partners_b.add(partner_h)

        missing_for_a = vicinal_partners_b - vicinal_partners_a
        missing_for_b = vicinal_partners_a - vicinal_partners_b

        for missing_h in missing_for_a:

            J_reciprocal = None
            J_reciprocal_source = None
            if missing_h in result:
                for p in result[missing_h].get("partners", []):
                    if p.get("partner_h_index") == h_a:
                        J_reciprocal = p.get("J_hz", 0)
                        J_reciprocal_source = p.get("J_source", "unknown")
                        break

            J_ref = None
            for p in partners_b:
                if p.get("partner_h_index") == missing_h:
                    J_ref = p.get("J_hz", 0)
                    break

            if J_reciprocal and J_reciprocal > 0:
                reciprocal_priority = source_priority.get(J_reciprocal_source, 1)

                if reciprocal_priority > 3:
                    J_final = J_reciprocal
                    J_source = J_reciprocal_source
                else:

                    if J_ref:
                        if J_ref > 8.0:
                            J_final = 3.0
                        elif J_ref < 5.0:
                            J_final = 10.0
                        else:
                            J_final = J_ref
                    else:
                        J_final = J_reciprocal
                    J_source = "diastereotopic_estimated"
            elif J_ref:

                if J_ref > 8.0:
                    J_final = 3.0
                elif J_ref < 5.0:
                    J_final = 10.0
                else:
                    J_final = J_ref
                J_source = "diastereotopic_estimated"
            else:
                continue

            if h_a in result:
                result[h_a]["partners"].append({
                    "partner_h_index": int(missing_h),
                    "partner_shift_ppm": round(float(proton_shifts_ppm.get(missing_h, 0.0)), 3),
                    "J_hz": round(float(J_final), 2),
                    "J_source": J_source,
                })

                j_couplings_a = result[h_a].get("J_couplings", [])
                j_couplings_a.append({"J_hz": round(float(J_final), 3), "n_equiv": 1})
                result[h_a]["J_couplings"] = j_couplings_a

        for missing_h in missing_for_b:

            J_reciprocal = None
            J_reciprocal_source = None
            if missing_h in result:
                for p in result[missing_h].get("partners", []):
                    if p.get("partner_h_index") == h_b:
                        J_reciprocal = p.get("J_hz", 0)
                        J_reciprocal_source = p.get("J_source", "unknown")
                        break

            J_ref = None
            for p in partners_a:
                if p.get("partner_h_index") == missing_h:
                    J_ref = p.get("J_hz", 0)
                    break

            if J_reciprocal and J_reciprocal > 0:
                reciprocal_priority = source_priority.get(J_reciprocal_source, 1)
                if reciprocal_priority > 3:
                    J_final = J_reciprocal
                    J_source = J_reciprocal_source
                else:
                    if J_ref:
                        if J_ref > 8.0:
                            J_final = 3.0
                        elif J_ref < 5.0:
                            J_final = 10.0
                        else:
                            J_final = J_ref
                    else:
                        J_final = J_reciprocal
                    J_source = "diastereotopic_estimated"
            elif J_ref:
                if J_ref > 8.0:
                    J_final = 3.0
                elif J_ref < 5.0:
                    J_final = 10.0
                else:
                    J_final = J_ref
                J_source = "diastereotopic_estimated"
            else:
                continue

            if h_b in result:
                result[h_b]["partners"].append({
                    "partner_h_index": int(missing_h),
                    "partner_shift_ppm": round(float(proton_shifts_ppm.get(missing_h, 0.0)), 3),
                    "J_hz": round(float(J_final), 2),
                    "J_source": J_source,
                })
                j_couplings_b = result[h_b].get("J_couplings", [])
                j_couplings_b.append({"J_hz": round(float(J_final), 3), "n_equiv": 1})
                result[h_b]["J_couplings"] = j_couplings_b

    for h_idx, data in result.items():
        partners = data.get("partners", [])
        if not partners:
            continue

        current_shift = data.get("corrected_shift_ppm")

        def is_invisible_geminal(partner):

            p_shift = partner.get("partner_shift_ppm")
            j_source = partner.get("J_source", "")

            if "geminal" not in j_source:
                return False

            if p_shift is not None and current_shift is not None:
                if abs(p_shift - current_shift) <= 0.005:
                    return True
            return False

        j_values_from_partners = [
            p.get("J_hz", 0) for p in partners
            if p.get("J_hz", 0) > 0 and not is_invisible_geminal(p)
        ]

        if j_values_from_partners:

            grouped = _group_J_values(j_values_from_partners, tol=0.25)
            new_J_couplings = []
            for g in grouped:
                if g:
                    J_mean = sum(g) / len(g)
                    new_J_couplings.append({
                        "J_hz": round(float(J_mean), 3),
                        "n_equiv": len(g),
                    })

            current_J = data.get("J_couplings", [])
            current_J_set = set(round(jc.get("J_hz", 0), 1) for jc in current_J)
            new_J_set = set(round(jc.get("J_hz", 0), 1) for jc in new_J_couplings)

            if len(new_J_couplings) >= len(current_J) or not current_J_set.issubset(new_J_set):

                if len(new_J_couplings) > 0:
                    data["J_couplings"] = new_J_couplings

    for h_idx, data in result.items():
        J_couplings = data.get("J_couplings", [])
        if not J_couplings:
            continue

        groups = []
        for jc in J_couplings:
            J_hz = jc.get("J_hz", 0)
            n_equiv = jc.get("n_equiv", 1)
            if J_hz > 0:
                groups.append([J_hz] * n_equiv)

        if not groups:
            continue

        new_label = _first_order_label_from_groups(groups, force_singlet_if_empty=False)
        old_label = data.get("label", "")

        if new_label != old_label:

            total_J = sum(jc.get("n_equiv", 1) for jc in J_couplings if jc.get("J_hz", 0) > 0)

            n_distinct_groups = len(groups)

            if total_J == 1 and old_label == "d":
                pass
            elif total_J == 2 and old_label == "d":

                data["label"] = new_label
                data["label_idx"] = LABEL2IDX.get(new_label, LABEL2IDX["m"])
            elif total_J >= 2 and old_label in ("d", "s"):

                data["label"] = new_label
                data["label_idx"] = LABEL2IDX.get(new_label, LABEL2IDX["m"])
            elif total_J >= 3 and old_label in ("d", "dd", "t"):

                data["label"] = new_label
                data["label_idx"] = LABEL2IDX.get(new_label, LABEL2IDX["m"])

            elif old_label == "t" and n_distinct_groups >= 2:

                J_vals_unique = [g[0] for g in groups if g]
                if len(J_vals_unique) >= 2:
                    J_diff = abs(J_vals_unique[0] - J_vals_unique[1])

                    if J_diff > 0.3:
                        data["label"] = "dd"
                        data["label_idx"] = LABEL2IDX.get("dd", LABEL2IDX["m"])

    ch3_heavy_atoms = set()
    for heavy_idx, h_list in heavy_to_H.items():
        if len(h_list) == 3:
            ch3_heavy_atoms.add(heavy_idx)

    if ch3_heavy_atoms:

        for h_idx, data in result.items():
            partners = data.get("partners", [])
            if not partners:
                continue

            partner_by_heavy: Dict[int, List[Dict]] = {}
            for p in partners:
                partner_h = p.get("partner_h_index")
                if partner_h is None:
                    continue

                for hv, h_list in heavy_to_H.items():
                    if partner_h in h_list:
                        if hv not in partner_by_heavy:
                            partner_by_heavy[hv] = []
                        partner_by_heavy[hv].append(p)
                        break

            for hv, partner_group in partner_by_heavy.items():
                if hv not in ch3_heavy_atoms:
                    continue
                if len(partner_group) < 2:
                    continue

                J_values = [p.get("J_hz", 0) for p in partner_group if p.get("J_hz", 0) > 0]
                if not J_values:
                    continue

                from collections import Counter, defaultdict

                J_groups: Dict[float, List[float]] = defaultdict(list)
                tolerance = 1.0

                for j_val in J_values:

                    found_group = None
                    for group_key in J_groups:
                        if abs(j_val - group_key) <= tolerance:
                            found_group = group_key
                            break

                    if found_group is not None:
                        J_groups[found_group].append(j_val)
                    else:

                        J_groups[j_val].append(j_val)

                best_group_key = max(J_groups.keys(), key=lambda k: len(J_groups[k]))
                best_group_values = J_groups[best_group_key]
                count = len(best_group_values)

                if count > len(J_values) / 2:

                    J_harmonized = sum(best_group_values) / len(best_group_values)
                else:

                    typical_groups = [(k, v) for k, v in J_groups.items() if 5.5 <= k <= 8.5]
                    if typical_groups:

                        best_typical = max(typical_groups, key=lambda x: len(x[1]))
                        J_harmonized = sum(best_typical[1]) / len(best_typical[1])
                    else:

                        J_harmonized = max(J_values)

                for p in partner_group:
                    old_J = p.get("J_hz", 0)
                    if abs(old_J - J_harmonized) > 0.01:
                        p["J_hz"] = round(J_harmonized, 2)

                for ch3_h in heavy_to_H.get(hv, []):
                    if ch3_h not in result:
                        continue
                    ch3_partners = result[ch3_h].get("partners", [])
                    for ch3_p in ch3_partners:
                        if ch3_p.get("partner_h_index") == h_idx:
                            old_ch3_J = ch3_p.get("J_hz", 0)
                            if abs(old_ch3_J - J_harmonized) > 0.01:
                                ch3_p["J_hz"] = round(J_harmonized, 2)

    for h_idx, data in result.items():
        partners = data.get("partners", [])
        if not partners:
            continue

        current_shift = data.get("corrected_shift_ppm")

        def is_invisible_geminal_p8(partner):

            p_shift = partner.get("partner_shift_ppm")
            j_source = partner.get("J_source", "")

            if "geminal" not in j_source:
                return False

            if p_shift is not None and current_shift is not None:
                if abs(p_shift - current_shift) <= 0.005:
                    return True
            return False

        j_values_from_partners = [
            p.get("J_hz", 0) for p in partners
            if p.get("J_hz", 0) > 0 and not is_invisible_geminal_p8(p)
        ]

        if j_values_from_partners:

            grouped = _group_J_values(j_values_from_partners, tol=0.25)
            new_J_couplings = []
            for g in grouped:
                if g:

                    J_mean = sum(g) / len(g)
                    new_J_couplings.append({
                        "J_hz": round(float(J_mean), 2),
                        "n_equiv": len(g),
                    })

            if new_J_couplings:
                data["J_couplings"] = new_J_couplings

    return result

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

    enable_dynamics: bool = False,
    ph: float = 7.4,
) -> Dict[int, Dict[str, Any]]:

    if config is None:
        config = SimConfig()

    if config.infer_double_bond_stereo:
        defined, undefined = _count_undefined_double_bonds(mol)
        if undefined > 0:
            if config.infer_stereo_log:
                import logging
                logging.info(f" Inférence stéréochimie: {undefined} double(s) liaison(s) non définies → {config.infer_stereo_default}")
            mol = _infer_double_bond_stereo(mol, config.infer_stereo_default)

    has_explicit_h = any(a.GetSymbol() == "H" for a in mol.GetAtoms())
    if has_explicit_h:
        molH = mol
    else:
        molH = Chem.AddHs(mol, addCoords=True)

    _validate_connectivity(molH, x_h_connectivity)

    H_ALL = [a.GetIdx() for a in molH.GetAtoms() if a.GetSymbol() == "H"]
    heavy_to_H: Dict[int, List[int]] = defaultdict(list)
    for h in H_ALL:
        nb = [n.GetIdx() for n in molH.GetAtomWithIdx(h).GetNeighbors()]
        if nb:
            heavy_to_H[nb[0]].append(h)

    expanded_connectivity: Dict[int, List[int]] = {}
    proton_to_heavy: Dict[int, int] = {}
    geminal_partners: Dict[int, Optional[int]] = {}

    diastereotopic_shift_deltas: Dict[int, float] = {}

    for heavy_idx, h_list in sorted(x_h_connectivity.items(), key=lambda kv: kv[0]):
        if not h_list:
            continue

        heavy_atom = molH.GetAtomWithIdx(heavy_idx)
        actual_h_neighbors = [n for n in heavy_atom.GetNeighbors() if n.GetSymbol() == "H"]
        n_h_attached = len(actual_h_neighbors)

        is_terminal_vinyl_CH2 = False
        if n_h_attached == 2 and len(h_list) == 2:
            hyb = heavy_atom.GetHybridization()
            if hyb == rdchem.HybridizationType.SP2:

                c_neighbors = [n for n in heavy_atom.GetNeighbors() if n.GetSymbol() != "H"]
                if len(c_neighbors) == 1:

                    for bond in heavy_atom.GetBonds():
                        if bond.GetBondType() == rdchem.BondType.DOUBLE:
                            is_terminal_vinyl_CH2 = True
                            break

        if is_terminal_vinyl_CH2:

            h_a, h_b = sorted(h_list)
            expanded_connectivity[heavy_idx] = [h_a]
            pseudo_key = -(h_b + 1)
            expanded_connectivity[pseudo_key] = [h_b]
            proton_to_heavy[h_a] = heavy_idx
            proton_to_heavy[h_b] = heavy_idx
            geminal_partners[h_a] = h_b
            geminal_partners[h_b] = h_a

            delta_ppm = 0.10

            if DIASTEREOTOPIC_TABLE_AVAILABLE:
                try:

                    delta_ppm, delta_reason = get_diastereotopic_delta(
                        mol=molH,
                        h1_idx=h_a,
                        h2_idx=h_b,
                        heavy_idx=heavy_idx,
                        spectrometer_mhz=spectrometer_mhz,
                        solvent=solvent
                    )
                except Exception:

                    delta_ppm = 0.10

            if delta_ppm > 0 and proton_shifts_ppm:
                shift_a = proton_shifts_ppm.get(h_a, 0.0)
                shift_b = proton_shifts_ppm.get(h_b, 0.0)

                if abs(shift_a - shift_b) < 0.01:
                    diastereotopic_shift_deltas[h_a] = -delta_ppm / 2.0
                    diastereotopic_shift_deltas[h_b] = +delta_ppm / 2.0
                    if config.verbose:
                        print(f"[CH2_vinyl] Applied Δδ={delta_ppm:.3f} ppm to H{h_a}/H{h_b} (shifts identical: {shift_a:.3f} ppm)")
                else:
                    if config.verbose:
                        print(f"[CH2_vinyl] Skipped Δδ correction for H{h_a}/H{h_b}: ML already predicts different shifts ({shift_a:.3f} vs {shift_b:.3f} ppm, Δ={abs(shift_a-shift_b):.3f} ppm)")

            continue

        if n_h_attached != 2:
            expanded_connectivity[heavy_idx] = list(h_list)
            for h in h_list:
                proton_to_heavy[h] = heavy_idx
                geminal_partners[h] = None
            continue

        is_diastereotopic = False
        diastereotopic_reason = ""
        should_separate = True

        if config.separate_diastereotopic_protons and len(h_list) == 2:

            neighbors = [n for n in heavy_atom.GetNeighbors() if n.GetSymbol() != "H"]
            for nb in neighbors:
                if nb.GetSymbol() == "C":
                    tag = nb.GetChiralTag()
                    if tag in (rdchem.ChiralType.CHI_TETRAHEDRAL_CW, rdchem.ChiralType.CHI_TETRAHEDRAL_CCW):
                        is_diastereotopic = True
                        diastereotopic_reason = "adjacent_chiral"
                        break

            if not is_diastereotopic:
                is_diast_bfs, delta_ppm_chiral = _detect_diastereotopic_protons(molH, heavy_idx)
                if is_diast_bfs:
                    is_diastereotopic = True
                    diastereotopic_reason = "chiral_bfs"

            if not is_diastereotopic:
                is_prochiral, prochiral_reason = _is_prochiral_CH2(molH, heavy_idx)
                if is_prochiral:
                    is_diastereotopic = True
                    diastereotopic_reason = prochiral_reason

            if not is_diastereotopic and len(h_list) == 2:
                h_a, h_b = h_list
                if _is_diastereotopic_by_prochirality(molH, heavy_idx, h_a, h_b):
                    is_diastereotopic = True
                    diastereotopic_reason = "prochirality_analysis"

            if is_diastereotopic and len(h_list) == 2:
                h_a, h_b = h_list
                try:
                    should_merge, merge_reason, merge_details = should_merge_diastereotopic_CH2(
                        molH, heavy_idx, h_a, h_b,
                        temperature_k=temperature_k,
                        spectrometer_mhz=spectrometer_mhz
                    )
                    if should_merge:

                        should_separate = False
                        if config.verbose:
                            print(f"[dynamics] CH₂ {heavy_idx}: {merge_reason}")
                except Exception as e:
                    if config.verbose:
                        print(f"[dynamics] Erreur pour CH₂ {heavy_idx}: {e}")

                    should_separate = True

        if is_diastereotopic and should_separate:

            h_a, h_b = h_list
            expanded_connectivity[heavy_idx] = [h_a]
            pseudo_key = -(h_b + 1)
            expanded_connectivity[pseudo_key] = [h_b]
            proton_to_heavy[h_a] = heavy_idx
            proton_to_heavy[h_b] = heavy_idx
            geminal_partners[h_a] = h_b
            geminal_partners[h_b] = h_a

            delta_ppm = 0.10

            if DIASTEREOTOPIC_TABLE_AVAILABLE:
                try:

                    delta_ppm, delta_reason = get_diastereotopic_delta(
                        mol=molH,
                        h1_idx=h_a,
                        h2_idx=h_b,
                        heavy_idx=heavy_idx,
                        spectrometer_mhz=spectrometer_mhz,
                        solvent=solvent
                    )
                except Exception:

                    is_diast, delta_from_dist = _detect_diastereotopic_protons(molH, heavy_idx)
                    delta_ppm = delta_from_dist if is_diast else 0.10
            else:

                is_diast, delta_from_dist = _detect_diastereotopic_protons(molH, heavy_idx)
                delta_ppm = delta_from_dist if is_diast else 0.10

            delta_ppm = max(delta_ppm, 0.005)

            if delta_ppm > 0 and proton_shifts_ppm:
                shift_a = proton_shifts_ppm.get(h_a, 0.0)
                shift_b = proton_shifts_ppm.get(h_b, 0.0)

                if abs(shift_a - shift_b) < 0.01:

                    diastereotopic_shift_deltas[h_a] = -delta_ppm / 2.0
                    diastereotopic_shift_deltas[h_b] = +delta_ppm / 2.0
                    if config.verbose:
                        print(f"[CH2_diast] Applied Δδ={delta_ppm:.3f} ppm to H{h_a}/H{h_b} (reason: {diastereotopic_reason}, shifts identical: {shift_a:.3f} ppm)")
                else:
                    if config.verbose:
                        print(f"[CH2_diast] Skipped Δδ correction for H{h_a}/H{h_b} (reason: {diastereotopic_reason}): ML already predicts different shifts ({shift_a:.3f} vs {shift_b:.3f} ppm, Δ={abs(shift_a-shift_b):.3f} ppm)")
        elif is_diastereotopic and not should_separate:

            expanded_connectivity[heavy_idx] = list(h_list)
            for h in h_list:
                proton_to_heavy[h] = heavy_idx
                geminal_partners[h] = None

        else:

            expanded_connectivity[heavy_idx] = list(h_list)
            for h in h_list:
                proton_to_heavy[h] = heavy_idx
                geminal_partners[h] = None

    solv = SOLVENT_DB.get(solvent, SOLVENT_DB[DEFAULT_SOLVENT])
    eff_regime = _effective_exchange_regime(solv, temperature_k)

    corrected_shifts = proton_shifts_ppm
    if config.use_solvent_shift_corrections and proton_shifts_ppm:
        corrected_shifts = _apply_solvent_corrections(
            proton_shifts_ppm, molH, heavy_to_H, solv.get("name", solvent)
        )

    if corrected_shifts and diastereotopic_shift_deltas:
        corrected_shifts = dict(corrected_shifts)
        for h_idx, delta in diastereotopic_shift_deltas.items():
            if h_idx in corrected_shifts:
                corrected_shifts[h_idx] = corrected_shifts[h_idx] + delta

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

    n_heavy = base.GetNumHeavyAtoms()
    n_rings = len(base.GetRingInfo().AtomRings())
    skip_3d_reason = None

    if config.skip_3d_for_large_molecules:

        if n_heavy > config.skip_3d_above_heavy_atoms:
            skip_3d_reason = f"skip_3d: {n_heavy} atomes > {config.skip_3d_above_heavy_atoms}"

        elif n_rings > config.skip_3d_above_fused_rings:
            skip_3d_reason = f"skip_3d: {n_rings} cycles > {config.skip_3d_above_fused_rings}"

        if skip_3d_reason is None:
            coverage = _check_tables_coverage(base)
            if coverage.get('all_covered', False):
                skip_3d_reason = f"skip_3d: tables couvrent tout"

    if skip_3d_reason:
        use_3d = False
        if DEBUG_LOGGING:
            print(f"[FAST MODE] {skip_3d_reason}")
            print(f"[FAST MODE] Utilisation des tables de J uniquement (pas de génération 3D)")

    if use_3d and molH.GetNumConformers() == 0:
        try:
            try:
                ps = AllChem.ETKDGv3()
            except Exception:
                ps = AllChem.ETKDG()
            ps.randomSeed = int(config.random_seed)

            if n_heavy > 50:
                ps.maxAttempts = 5
                ps.useRandomCoords = True
            res = AllChem.EmbedMolecule(molH, ps)
            if res == 0:
                _optimize_molecule_geometry(molH, config.force_field, config.force_field_max_iters, 0)
        except Exception:
            pass

    D = _topological_distances(base)
    aromatic_dmap = _precompute_aromatic_distance_map(base)
    heavy_is_sp2_cache = [_is_sp2_like(base.GetAtomWithIdx(i)) for i in range(base.GetNumAtoms())]
    ranks = _canonical_ranks(base)

    n_heavy_atoms = base.GetNumHeavyAtoms()
    use_averaging = config.use_conformational_averaging

    ultra_fast_mode = (config.ultra_fast_mode and n_heavy_atoms > config.ultra_fast_above_heavy_atoms)

    if ultra_fast_mode:
        use_averaging = False
        if DEBUG_LOGGING:
            print(f"[ULTRA-FAST MODE] Disabling conformational averaging (ultra_fast_mode=True)")

    if config.fast_mode_for_large_molecules and n_heavy_atoms > config.auto_disable_averaging_above_atoms:
        use_averaging = False
        if DEBUG_LOGGING:
            print(f"[FAST MODE] Molecule has {n_heavy_atoms} heavy atoms > {config.auto_disable_averaging_above_atoms} threshold")
            print(f"[FAST MODE] Disabling conformational averaging for speed (tables-only mode)")

    smart_skip_reason = None
    if use_averaging and config.smart_conformer_mode:
        needs_conf, smart_reason = _molecule_needs_conformers(base)
        if not needs_conf:
            use_averaging = False
            smart_skip_reason = smart_reason
            if DEBUG_LOGGING:
                print(f"[SMART MODE] Conformers NOT needed: {smart_reason}")
                print(f"[SMART MODE] Using tables + fixed J values for flexible bonds")

    conformer_cache = None
    if use_averaging:

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

    if ultra_fast_mode and DEBUG_LOGGING:
        print(f"[ULTRA-FAST MODE] Molécule {n_heavy_atoms} atomes > {config.ultra_fast_above_heavy_atoms} seuil")
        print(f"[ULTRA-FAST MODE] Skip détections diastéréotopiques coûteuses, J vicinal = {config.ultra_fast_default_vicinal_J} Hz")

    diastereotopic_cache: Dict[int, Tuple[bool, float]] = {}
    chiral_indices_cache: set = set()

    if not ultra_fast_mode:

        try:
            chiral_centers_list = Chem.FindMolChiralCenters(base, includeUnassigned=True)
            chiral_indices_cache = {idx for idx, _ in chiral_centers_list}
        except Exception:
            chiral_indices_cache = set()

        for c_idx in heavy_to_H.keys():
            if len(heavy_to_H[c_idx]) == 2:

                is_diast, delta = _detect_diastereotopic_protons(
                    molH, c_idx, chiral_indices=chiral_indices_cache
                )

                if not is_diast:
                    is_prochiral, prochiral_reason = _is_prochiral_CH2(molH, c_idx)
                    if is_prochiral:
                        is_diast = True
                        delta = 0.10

                diastereotopic_cache[c_idx] = (is_diast, delta)

    result: Dict[int, Dict[str, Any]] = {}

    for key, h_list in sorted(expanded_connectivity.items(), key=lambda kv: kv[0]):
        if not h_list:
            continue

        h_idx = h_list[0]
        true_heavy_idx = proton_to_heavy[h_idx]
        geminal_partner = geminal_partners.get(h_idx, None)

        broadening = None
        if config.mark_labile_protons_broad:
            broadening = _estimate_labile_proton_broadening(
                molH, true_heavy_idx, solv.get("name", solvent), temperature_k
            )

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
                "heterocycle_position": 2,
                "sugar_position": None,
            }
            continue

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
                "broadening": broadening,
                "heterocycle_position": None,
                "sugar_position": None,
            }
            continue

        (J_vals, second_order, partners, sat, homo_pairs,
         partner_heavy_idx_kept, partner_is_sp2_kept,
         second_order_type, system_type, delta_nu_min_hz, J_orders) = _estimate_Js_for_group(
            molH, base, D, aromatic_dmap, true_heavy_idx, h_list, heavy_to_H,
            spectrometer_mhz, solv, corrected_shifts, include_heteronuclear,
            suppressed_heavy, heavy_is_sp2_cache, config, ranks=ranks,
            delta_cluster_id={},
            geminal_partner=geminal_partner,
            conformer_cache=conformer_cache,
            ultra_fast_mode=ultra_fast_mode,
            diastereotopic_cache=diastereotopic_cache
        )

        if enable_dynamics and DYNAMICS_MODULE_AVAILABLE and J_vals and partners:
            J_vals_corrected = []
            partners_corrected = []
            for J_static, partner_h_idx in zip(J_vals, partners):
                J_effective = apply_dynamics_to_coupling(
                    mol=molH,
                    h1_idx=h_idx,
                    h2_idx=partner_h_idx,
                    J_static=J_static,
                    ph=ph,
                    solvent=solv.get("name", solvent),
                    temperature_k=temperature_k,
                    use_cache=True
                )

                if J_effective is not None and abs(J_effective) > 0.001:
                    J_vals_corrected.append(J_effective)
                    partners_corrected.append(partner_h_idx)

            J_vals = J_vals_corrected
            partners = partners_corrected

        heterocycle_position = None

        if config.advanced_spin_system_detection and true_heavy_idx in spin_systems:
            detected_sys, spin_info = spin_systems[true_heavy_idx]
            if detected_sys and system_type is None:
                system_type = detected_sys

            if detected_sys in ("AB", "AB₂", "AA'BB'", "A₂B₃", "ABX", "complex_second_order"):
                if not second_order:
                    second_order = True
                    second_order_type = "strong_coupling"

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

        heavy = base.GetAtomWithIdx(true_heavy_idx)
        if heavy.GetSymbol() == "C" and heavy.GetIsAromatic():
            if _is_symmetric_para_disubstituted_benzene(base, true_heavy_idx, aromatic_dmap, ranks):

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

        if heavy_is_sp2_cache[true_heavy_idx] and not J_vals and heavy.GetSymbol() == "C" and heavy.GetIsAromatic() and not is_heterocycle_managed:
            J_fallback = _fallback_aromatic_J_when_empty(base, aromatic_dmap, true_heavy_idx, heavy_to_H, config, ranks=ranks)
            if J_fallback:
                J_vals = J_fallback
                second_order = False
                second_order_type = "first_order"

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

        if _is_isolated_methyl(base, true_heavy_idx, heavy_to_H):
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

        J_vals_original = list(J_vals)

        if config.use_frequency_dependent_resolution and J_vals:

            J_vals_filtered, kept_indices = filter_J_by_frequency(
                J_vals, spectrometer_mhz, coupling_orders=None, config=config
            )

            if len(J_vals_filtered) < len(J_vals):
                if config.verbose:
                    removed = [j for j in J_vals if j not in J_vals_filtered]
                    print(f"[FREQ_FILTER] H{h_idx} @ {spectrometer_mhz}MHz: filtré {removed} (sous seuil)")

                homo_pairs_filtered: List[Tuple[int, float, str]] = []
                for hp_tuple in homo_pairs:
                    partner_h_idx, J_val = hp_tuple[0], hp_tuple[1]
                    hp_source = hp_tuple[2] if len(hp_tuple) > 2 else "unknown"
                    if J_val in J_vals_filtered or abs(J_val) >= get_resolution_thresholds(spectrometer_mhz)["min_resolvable"]:
                        homo_pairs_filtered.append((partner_h_idx, J_val, hp_source))
                homo_pairs = homo_pairs_filtered

                J_vals = J_vals_filtered

        if not J_vals:
            result[h_idx] = {
                "label": "s",
                "label_idx": LABEL2IDX["s"],
                "J_couplings": [],
                "peak_positions_hz": [0.0],
                "peak_intensities": [1.0],
                "second_order_type": "first_order",
                "system_type": system_type,
                "virtual_coupling": False,
                "partners": [],
                "satellites_13C_hz": sat or [],
                "broadening": broadening,
                "heterocycle_position": None,
                "sugar_position": None,
                "unresolved_J": J_vals_original if config.verbose else None,
            }
            continue

        anis_high = (solv.get("anisotropy", "normal") == "high")

        heavy = base.GetAtomWithIdx(true_heavy_idx)
        if heavy.GetSymbol() == "C" and heavy.GetIsAromatic() and homo_pairs and J_vals:

            suspect_J = 12.0
            tol_suspect = 0.5
            has_suspect = any(abs(j - suspect_J) < tol_suspect for j in J_vals)
            correct_aromatic_J = [7.8, 7.5, 8.0, 8.2, 2.3, 2.0, 1.5, 1.8, 0.6, 0.5]
            has_correct = any(any(abs(j - av) < 0.5 for av in correct_aromatic_J) for j in J_vals)

            if has_suspect and has_correct:

                filtered_J_vals = []
                filtered_homo_pairs: List[Tuple[int, float, str]] = []
                filtered_J_orders = []
                filtered_partners = []
                seen_partners_filter: set = set()

                for i, hp_tuple in enumerate(homo_pairs):
                    partner_h, j_val = hp_tuple[0], hp_tuple[1]
                    j_source_hp = hp_tuple[2] if len(hp_tuple) > 2 else "unknown"
                    is_suspect = abs(j_val - suspect_J) < tol_suspect
                    is_correct = any(abs(j_val - av) < 0.5 for av in correct_aromatic_J)

                    if is_correct or partner_h not in seen_partners_filter:
                        if not is_suspect or partner_h not in seen_partners_filter:
                            filtered_homo_pairs.append((partner_h, j_val, j_source_hp))
                            if i < len(J_vals):
                                filtered_J_vals.append(J_vals[i])
                            if J_orders and i < len(J_orders):
                                filtered_J_orders.append(J_orders[i])
                            if i < len(partners):
                                filtered_partners.append(partners[i])
                            seen_partners_filter.add(partner_h)

                if len(filtered_J_vals) > 0:
                    homo_pairs = filtered_homo_pairs
                    J_vals = filtered_J_vals
                    if filtered_J_orders:
                        J_orders = filtered_J_orders
                    if filtered_partners:
                        partners = filtered_partners

        if geminal_partner is not None and config.separate_geminal_vicinal_J and J_orders:
            groups = _group_J_values_with_orders(
                J_vals, J_orders,
                tol=_field_based_tol,
                anisotropy_high=anis_high,
                config=config,
                separate_geminal_vicinal=True
            )
        else:
            groups = _group_J_values(J_vals, tol=_field_based_tol, anisotropy_high=anis_high, config=config)

        if not groups:
            groups = [[]]

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

        heterocycle_expected_mult = None
        heterocycle_position = None
        sugar_expected_mult = None
        sugar_position = None

        if is_heterocycle_managed and system_type and ':' in system_type:

            parts = system_type.split(':')
            if len(parts) >= 2:
                heterocycle_expected_mult = parts[1]
            if len(parts) >= 3:
                try:
                    heterocycle_position = int(parts[2])
                except (ValueError, TypeError):
                    heterocycle_position = parts[2] if parts[2] else None

            system_type = parts[0]

        if system_type and system_type.startswith('sugar_') and ':' in system_type:

            parts = system_type.split(':')
            if len(parts) >= 2:
                sugar_expected_mult = parts[1]
            if len(parts) >= 3:
                try:
                    sugar_position = int(parts[2])
                except (ValueError, TypeError):
                    sugar_position = parts[2] if parts[2] else None

            system_type = parts[0]

        if heterocycle_expected_mult and heterocycle_expected_mult in LABEL2IDX:

            label = heterocycle_expected_mult
        elif sugar_expected_mult and sugar_expected_mult in LABEL2IDX:

            label = sugar_expected_mult
        elif second_order:

            if second_order_type == "strong_coupling":
                label = "m"
            else:

                label2 = _first_order_label_from_groups(groups) if groups else "m"
                if heavy_is_sp2_cache[true_heavy_idx] and heavy.GetSymbol() == "C" and J_vals:
                    label2 = _apparent_sp2_label_from_J(J_vals, minJ=0.4)
                label = label2 or "m"
        else:
            label = _first_order_label_from_groups(groups)
            if heavy_is_sp2_cache[true_heavy_idx] and heavy.GetSymbol() == "C" and J_vals:
                label = _apparent_sp2_label_from_J(J_vals, minJ=0.4)

        label_idx = LABEL2IDX.get(label, LABEL2IDX["m"])

        if geminal_partner is not None and len(h_list) == 1:

            n_j_groups = len(groups) if groups else 0
            n_j_total = len(J_vals) if J_vals else 0

            has_geminal = any(abs(J - _get_geminal_2J(heavy, base)) < 2.0 for J in J_vals) if J_vals else False

            if n_j_total >= 2:

                if label in ("t", "q") and n_j_groups >= 1:

                    if n_j_total == 2 and len(groups) == 1:

                        label = "dd"
                        label_idx = LABEL2IDX["dd"]
                    elif n_j_total >= 3:

                        if len(groups) == 1:

                            label = "ddd"
                            label_idx = LABEL2IDX["ddd"]
                        elif len(groups) == 2:

                            g_sizes = sorted([len(g) for g in groups])
                            if g_sizes == [1, 2]:
                                label = "dt"
                                label_idx = LABEL2IDX["dt"]
                            else:
                                label = "ddd"
                                label_idx = LABEL2IDX["ddd"]
                        else:
                            label = "ddd"
                            label_idx = LABEL2IDX["ddd"]

        if label == "t" and len(homo_pairs) >= 2 and geminal_partner is None:

            partner_h_indices = [hp[0] for hp in homo_pairs if hp[1] > 2.5]
            if len(partner_h_indices) >= 2:

                partner_heavies = []

                for partner_h_idx in partner_h_indices:
                    for hv, hs in heavy_to_H.items():
                        if partner_h_idx in hs:
                            partner_heavies.append(hv)
                            break

                if len(partner_heavies) >= 2 and len(set(partner_heavies)) == 1:

                    shared_heavy = partner_heavies[0]
                    is_diast, _ = _detect_diastereotopic_protons(molH, shared_heavy)
                    if is_diast:

                        diast_J_vals = [hp[1] for hp in homo_pairs if hp[0] in partner_h_indices[:2]]
                        if len(diast_J_vals) >= 2:

                            J_diff = abs(diast_J_vals[0] - diast_J_vals[1])
                            if J_diff > 2.0:

                                label = "dd"
                                label_idx = LABEL2IDX["dd"]
                            else:

                                pass
                        else:

                            label = "dd"
                            label_idx = LABEL2IDX["dd"]

        complicated = label in {"ddt", "dtt", "ddd", "dddd", "ttq", "qqt", "tqq", "qtq"}
        if solv.get("viscosity", "low") != "low" and complicated:
            if (J_vals and max(J_vals) <= config.viscosity_simplify_maxJ_hz) or (len(groups) > config.viscosity_simplify_max_groups):
                label = "m"
                label_idx = LABEL2IDX["m"]

        if homo_pairs:

            partner_heavy_map: Dict[int, List[int]] = {}
            for hp_tuple in homo_pairs:
                partner_h_idx = hp_tuple[0]
                for hv_idx, hv_list in heavy_to_H.items():
                    if partner_h_idx in hv_list:
                        partner_heavy_map[partner_h_idx] = hv_idx
                        break

            ch3_groups: Dict[int, List[Tuple[int, float, str]]] = {}
            for i, hp_tuple in enumerate(homo_pairs):
                partner_h_idx, J_val = hp_tuple[0], hp_tuple[1]
                hp_source = hp_tuple[2] if len(hp_tuple) > 2 else "unknown"

                if partner_h_idx in partner_heavy_map:
                    hv = partner_heavy_map[partner_h_idx]

                    if hv in heavy_to_H and len(heavy_to_H[hv]) == 3:
                        if hv not in ch3_groups:
                            ch3_groups[hv] = []
                        ch3_groups[hv].append((i, J_val, hp_source))

            homo_pairs_list = list(homo_pairs)
            for hv_idx, entries in ch3_groups.items():
                if len(entries) < 2:
                    continue

                J_vals_ch3 = [e[1] for e in entries if e[1] > 0]
                if not J_vals_ch3:
                    continue

                J_mean = sum(J_vals_ch3) / len(J_vals_ch3)

                J_filtered = []
                for j in J_vals_ch3:
                    if abs(j - J_mean) < 3.0:
                        J_filtered.append(j)

                if J_filtered:
                    J_harmonized = sum(J_filtered) / len(J_filtered)
                else:
                    J_harmonized = J_mean

                source_priority = {
                    "vinylic": 20, "alkyne": 20, "aromatic_exact": 19, "heterocycle": 19,
                    "heterocycle_table": 19, "heterocycle_expected_mult": 19, "heterocycle_singlet": 19,
                    "saturated_heterocycle": 18, "dkp_pyrrolidine_fused_table": 18, "pyrrolidine_fused_table": 18,
                    "diketopiperazine_experimental": 18, "steroid_table": 18, "terpene_table": 18,
                    "macrocycle_table": 18, "nucleotide_table": 18, "coupling_tables": 18, "coupling_extensions": 18,
                    "peptide_HNCH": 17, "sugar": 16, "structural": 15,
                    "geminal_calc": 14, "long_range_4J": 12, "long_range_5J": 12,
                    "fallback_ortho": 10, "fallback_meta": 10,
                    "karplus_3d": 5,
                    "aromatic_para": 4, "gauche_effect_double": 4, "gauche_effect_O": 4, "gauche_effect_F": 4,
                    "extended_chain_cation": 4,
                    "fallback_average": 3, "diastereotopic_estimated": 2, "vinylic_default": 2,
                    "peptide": 2, "unknown": 1
                }
                best_source = max(entries, key=lambda e: source_priority.get(e[2], 1))[2]

                for idx, _, _ in entries:
                    old_tuple = homo_pairs_list[idx]
                    homo_pairs_list[idx] = (old_tuple[0], round(J_harmonized, 2), best_source)

            homo_pairs = homo_pairs_list

        partner_info: List[Dict[str, Any]] = []
        seen_partners: Dict[int, Tuple[int, float]] = {}

        for hp_tuple in homo_pairs:
            partner_h_idx, J_val = hp_tuple[0], hp_tuple[1]
            hp_source = hp_tuple[2] if len(hp_tuple) > 2 else "unknown"
            if partner_h_idx not in proton_shifts_ppm:
                continue

            if partner_h_idx not in seen_partners:

                idx = len(partner_info)

                partner_shift = corrected_shifts.get(partner_h_idx, proton_shifts_ppm.get(partner_h_idx, 0.0))
                partner_info.append({
                    "partner_h_index": int(partner_h_idx),
                    "partner_shift_ppm": round(float(partner_shift), 3),
                    "J_hz": round(float(J_val), 2),
                    "J_source": hp_source,
                })
                seen_partners[partner_h_idx] = (idx, J_val, hp_source)
            else:

                prev_idx, prev_J, prev_source = seen_partners[partner_h_idx]

                source_priority = {
                    "vinylic": 20, "alkyne": 20, "aromatic_exact": 19, "heterocycle": 19,
                    "heterocycle_table": 19, "heterocycle_expected_mult": 19, "heterocycle_singlet": 19,
                    "saturated_heterocycle": 18, "dkp_pyrrolidine_fused_table": 18, "pyrrolidine_fused_table": 18,
                    "diketopiperazine_experimental": 18, "steroid_table": 18, "terpene_table": 18,
                    "macrocycle_table": 18, "nucleotide_table": 18, "coupling_tables": 18, "coupling_extensions": 18,
                    "peptide_HNCH": 17, "sugar": 16, "structural": 15,
                    "geminal_calc": 14, "long_range_4J": 12, "long_range_5J": 12,
                    "fallback_ortho": 10, "fallback_meta": 10,
                    "karplus_3d": 5,
                    "aromatic_para": 4, "gauche_effect_double": 4, "gauche_effect_O": 4, "gauche_effect_F": 4,
                    "extended_chain_cation": 4,
                    "fallback_average": 3, "diastereotopic_estimated": 2, "vinylic_default": 2,
                    "peptide": 2, "unknown": 1
                }
                new_priority = source_priority.get(hp_source, 1)
                old_priority = source_priority.get(prev_source, 1)

                if new_priority > old_priority:
                    partner_info[prev_idx]["J_hz"] = round(float(J_val), 2)
                    partner_info[prev_idx]["J_source"] = hp_source
                    seen_partners[partner_h_idx] = (prev_idx, J_val, hp_source)

                elif new_priority == old_priority and abs(J_val - prev_J) > 3.0:

                    if 5.5 < prev_J < 8.5 and not (5.5 < J_val < 8.5):
                        partner_info[prev_idx]["J_hz"] = round(float(J_val), 2)
                        partner_info[prev_idx]["J_source"] = hp_source
                        seen_partners[partner_h_idx] = (prev_idx, J_val, hp_source)

        exact_transitions = None
        if (config.use_exact_spin_simulation and system_type is not None
            and system_type in ("AB", "AB_weak", "AB₂", "AA'BB'", "A₂X₂", "ABX", "A₂X₃", "A₂B₃")):

            h_indices = list(h_list)

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

        if exact_transitions is not None and len(exact_transitions) > 0:

            peaks = [(t[0], t[1]) for t in exact_transitions]
        elif second_order and second_order_type in {
            "strong_coupling", "moderate_second_order", "weak_second_order"
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

        corrected_shift = corrected_shifts.get(h_idx, proton_shifts_ppm.get(h_idx, 0.0)) if corrected_shifts else proton_shifts_ppm.get(h_idx, 0.0)

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
            "broadening": broadening,
            "heterocycle_position": heterocycle_position,
            "sugar_position": sugar_position,
            "corrected_shift_ppm": round(float(corrected_shift), 6),
        }

    n_heavy_atoms = mol.GetNumHeavyAtoms() if mol else 0
    is_large_mol = n_heavy_atoms > 40

    result = _post_process_reciprocity_and_partners(
        result, corrected_shifts, heavy_to_H,
        is_large_molecule=is_large_mol
    )

    return result

def prepare_molecule_for_nmr(
    smiles: str,
    generate_3d: bool = True,
    random_seed: int = 42,
    skip_validation: bool = False,
    force_field: str = "MMFF94s",
) -> Tuple[Chem.Mol, Dict[int, List[int]], Dict[int, int]]:

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES invalide: {smiles}")

    if not skip_validation:
        _validate_molecular_structure(mol, smiles)

    molH = Chem.AddHs(mol, addCoords=True)

    if generate_3d and molH.GetNumConformers() == 0:

        n_heavy = mol.GetNumHeavyAtoms()
        n_rings = len(mol.GetRingInfo().AtomRings())
        should_skip_3d = n_heavy > 60 or n_rings > 6

        if should_skip_3d:

            pass
        else:
            try:
                try:
                    ps = AllChem.ETKDGv3()
                except AttributeError:
                    ps = AllChem.ETKDG()
                ps.randomSeed = int(random_seed)

                ps.maxAttempts = 10 if n_heavy > 50 else 50
                if n_heavy > 50:
                    ps.useRandomCoords = True
                res = AllChem.EmbedMolecule(molH, ps)
                if res == 0:
                    _optimize_molecule_geometry(molH, force_field, 200, 0)
            except Exception as e:
                if SimConfig().verbose:
                    print(f"[prepare_molecule_for_nmr] 3D generation failed: {e}")

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

    return molH, dict(heavy_to_H), proton_to_heavy

def _estimate_shifts_heuristic(molH: Chem.Mol, proton_to_heavy: Dict[int, int]) -> Dict[int, float]:

    shifts = {}

    for h_idx, heavy_idx in proton_to_heavy.items():
        heavy = molH.GetAtomWithIdx(heavy_idx)
        symbol = heavy.GetSymbol()
        hyb = heavy.GetHybridization()

        if symbol == 'C':
            if hyb == Chem.rdchem.HybridizationType.SP3:
                shift = 1.5

                for n in heavy.GetNeighbors():
                    if n.GetAtomicNum() in [7, 8, 9, 17, 35]:
                        shift += 1.5
                    elif n.GetSymbol() == 'C' and n.GetIsAromatic():
                        shift += 0.8
                shift = min(shift, 5.5)
            elif hyb == Chem.rdchem.HybridizationType.SP2:
                if heavy.GetIsAromatic():
                    shift = 7.2
                else:
                    shift = 5.5
            elif hyb == Chem.rdchem.HybridizationType.SP:
                shift = 2.5
            else:
                shift = 2.0
        elif symbol == 'N':
            shift = 6.5
        elif symbol == 'O':
            shift = 4.0
        else:
            shift = 2.0

        shifts[h_idx] = shift

    return shifts

def simulate_molecule(
    smiles: str,
    solvent: str = DEFAULT_SOLVENT,
    spectrometer_mhz: float = DEFAULT_SPECTROMETER_FREQ_MHZ,
    temperature_k: float = DEFAULT_TEMPERATURE_K,
    generate_3d: bool = True,
    config: Optional[SimConfig] = None,

    enable_dynamics: bool = False,
    ph: float = 7.4,
) -> Dict[int, Dict[str, Any]]:

    if config is None:
        config = SimConfig()

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES invalide ou non supporté: '{smiles}'")

    n_heavy = mol.GetNumHeavyAtoms()
    n_rings = len(mol.GetRingInfo().AtomRings())

    actual_generate_3d = generate_3d
    if config.skip_3d_for_large_molecules:
        if n_heavy > config.skip_3d_above_heavy_atoms:
            actual_generate_3d = False
            if config.verbose:
                print(f"[FAST MODE] Skip 3D: {n_heavy} atomes > {config.skip_3d_above_heavy_atoms} seuil")
        elif n_rings > config.skip_3d_above_fused_rings:
            actual_generate_3d = False
            if config.verbose:
                print(f"[FAST MODE] Skip 3D: {n_rings} cycles > {config.skip_3d_above_fused_rings} seuil")

    molH, connectivity, proton_to_heavy = prepare_molecule_for_nmr(
        smiles,
        generate_3d=actual_generate_3d,
        random_seed=config.random_seed,
    )

    shifts = _estimate_shifts_heuristic(molH, proton_to_heavy)

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
        enable_dynamics=enable_dynamics,
        ph=ph,
    )
