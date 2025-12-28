# =============================================================================
# sugars_complete.py
# Module COMPLET pour la simulation RMN des sucres
# Version finale avec toutes les corrections
# =============================================================================

"""
Module de simulation RMN pour les sucres (carbohydrates).

Fonctionnalités:
- Pyranoses (6 membres): glucose, galactose, mannose, etc.
- Furanoses (5 membres): ribose, fructose, arabinose, etc.
- Anomères α et β avec détection automatique
- Conformations: 4C1, 1C4, bateaux, enveloppes
- Fallback intelligent sans géométrie 3D (utilise stéréochimie)
- Sucres modifiés: déoxy, amino, acétyl, phosphate
- Multiplicités correctes par position
"""

from typing import Dict, List, Tuple, Optional, Set, Any
from dataclasses import dataclass, field
from enum import Enum
import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    raise ImportError("RDKit est requis.")


# =============================================================================
# SECTION 1: Énumérations
# =============================================================================

class RingType(Enum):
    """Type de cycle de sucre."""
    PYRANOSE = "pyranose"
    FURANOSE = "furanose"
    SEPTANOSE = "septanose"
    UNKNOWN = "unknown"


class Conformation(Enum):
    """Conformation du cycle de sucre."""
    CHAIR_4C1 = "4C1"
    CHAIR_1C4 = "1C4"
    BOAT_B14 = "B1,4"
    BOAT_B25 = "B2,5"
    BOAT_B36 = "B3,O"
    SKEW_1S3 = "1S3"
    SKEW_1S5 = "1S5"
    SKEW_3S1 = "3S1"
    SKEW_5S1 = "5S1"
    HALF_CHAIR = "H"
    ENVELOPE_E = "E"
    TWIST_T = "T"
    C2_ENDO = "C2-endo"
    C3_ENDO = "C3-endo"
    C2_EXO = "C2-exo"
    C3_EXO = "C3-exo"
    UNKNOWN = "unknown"


class ProtonOrientation(Enum):
    """Orientation d'un proton par rapport au plan du cycle."""
    AXIAL = "axial"
    EQUATORIAL = "equatorial"
    PSEUDO_AXIAL = "pseudo-axial"
    PSEUDO_EQUATORIAL = "pseudo-equatorial"
    UNKNOWN = "unknown"


class AnomericConfig(Enum):
    """Configuration anomérique du sucre."""
    ALPHA = "alpha"
    BETA = "beta"
    UNKNOWN = "unknown"


# =============================================================================
# SECTION 2: Tables de couplages
# =============================================================================

_PYRANOSE_3J: Dict[Tuple[str, str], float] = {
    ('axial', 'axial'): 9.5,
    ('axial', 'equatorial'): 3.8,
    ('equatorial', 'axial'): 3.8,
    ('equatorial', 'equatorial'): 2.5,
    ('pseudo-axial', 'pseudo-axial'): 7.0,
    ('pseudo-axial', 'pseudo-equatorial'): 3.5,
    ('pseudo-equatorial', 'pseudo-axial'): 3.5,
    ('pseudo-equatorial', 'pseudo-equatorial'): 3.0,
    ('unknown', 'unknown'): 5.0,
    ('axial', 'unknown'): 6.5,
    ('unknown', 'axial'): 6.5,
    ('equatorial', 'unknown'): 3.2,
    ('unknown', 'equatorial'): 3.2,
}

_FURANOSE_3J: Dict[Tuple[str, str], float] = {
    ('pseudo-axial', 'pseudo-axial'): 5.5,
    ('pseudo-axial', 'pseudo-equatorial'): 3.0,
    ('pseudo-equatorial', 'pseudo-axial'): 3.0,
    ('pseudo-equatorial', 'pseudo-equatorial'): 2.0,
    ('unknown', 'unknown'): 4.0,
}

_GEMINAL_2J: Dict[str, float] = {
    'pyranose_CH2': 12.5,
    'furanose_CH2': 12.0,
    'deoxysugar_CH2': 13.0,
}

_ANOMERIC_J: Dict[str, Dict[str, float]] = {
    'pyranose': {
        'alpha_4C1': 3.5,
        'beta_4C1': 7.8,
        'alpha_1C4': 7.8,
        'beta_1C4': 3.5,
        'unknown': 5.5,
    },
    'furanose': {
        'alpha': 4.5,
        'beta': 1.5,
        'unknown': 3.0,
    },
}

_CONFORMER_POPULATIONS: Dict[str, Dict[str, float]] = {
    'glucose': {'4C1': 0.98, '1C4': 0.02},
    'galactose': {'4C1': 0.95, '1C4': 0.05},
    'mannose': {'4C1': 0.97, '1C4': 0.03},
    'allose': {'4C1': 0.90, '1C4': 0.10},
    'altrose': {'4C1': 0.65, '1C4': 0.35},
    'gulose': {'4C1': 0.80, '1C4': 0.20},
    'idose': {'4C1': 0.40, '1C4': 0.60},
    'talose': {'4C1': 0.75, '1C4': 0.25},
    'ribose': {'C2-endo': 0.60, 'C3-endo': 0.40},
    'deoxyribose': {'C2-endo': 0.65, 'C3-endo': 0.35},
    'arabinose': {'C2-endo': 0.55, 'C3-endo': 0.45},
    'fructose': {'C2-exo': 0.55, 'C3-exo': 0.45},
    'default_pyranose': {'4C1': 0.95, '1C4': 0.05},
    'default_furanose': {'C2-endo': 0.50, 'C3-endo': 0.50},
}

# Multiplicités attendues par position
_TYPICAL_SUGAR_MULTIPLICITIES: Dict[str, Dict[int, Any]] = {
    'pyranose': {
        1: {'alpha': 'd', 'beta': 'd', 'unknown': 'd'},
        2: {'glucose': 'dd', 'mannose': 'dd', 'galactose': 'dd', 'default': 'dd'},
        3: {'glucose': 't', 'mannose': 'dd', 'galactose': 'dd', 'default': 't'},
        4: {'glucose': 't', 'mannose': 't', 'galactose': 'd', 'default': 't'},
        5: {'default': 'ddd'},
        6: {'default': 'dd'},
    },
    'furanose': {
        1: {'alpha': 'd', 'beta': 'd', 'unknown': 'd'},
        2: {'ribose': 'dd', 'deoxyribose': 'ddd', 'default': 'dd'},
        3: {'ribose': 'dd', 'deoxyribose': 'dt', 'default': 'dd'},
        4: {'default': 'ddd'},
        5: {'default': 'dd'},
    },
}


# =============================================================================
# SECTION 3: Structures de données
# =============================================================================

@dataclass
class SugarRing:
    """Représente un cycle de sucre détecté."""
    ring_type: RingType
    ring_atoms: List[int]
    ring_oxygen: int
    anomeric_carbon: Optional[int]
    atom_positions: Dict[int, int]
    conformation: Conformation = Conformation.UNKNOWN
    conformer_populations: Dict[str, float] = field(default_factory=dict)
    proton_orientations: Dict[int, ProtonOrientation] = field(default_factory=dict)
    anomeric_config: AnomericConfig = AnomericConfig.UNKNOWN
    sugar_type: str = "unknown"


@dataclass
class SugarModification:
    """Modifications chimiques sur le sucre."""
    deoxy_positions: List[int] = field(default_factory=list)
    amino_positions: List[int] = field(default_factory=list)
    acetyl_positions: List[int] = field(default_factory=list)
    phosphate_positions: List[int] = field(default_factory=list)
    sulfate_positions: List[int] = field(default_factory=list)
    methyl_positions: List[int] = field(default_factory=list)


# =============================================================================
# SECTION 4: Détection des cycles de sucre
# =============================================================================

def _is_sugar_ring_oxygen(mol: Chem.Mol, atom_idx: int, ring_atoms: List[int]) -> bool:
    """Vérifie si un oxygène est l'oxygène d'un cycle de sucre."""
    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.GetSymbol() != 'O':
        return False
    ring_neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() in ring_atoms]
    return len(ring_neighbors) == 2


def _is_pyranose_ring(mol: Chem.Mol, ring_atoms: List[int]) -> bool:
    """Vérifie si un cycle de 6 membres est un pyranose."""
    if len(ring_atoms) != 6:
        return False
    
    n_oxygen = 0
    n_carbon = 0
    n_nitrogen = 0
    n_hydroxyl = 0
    has_anomeric_oxygen = False
    ring_oxygen_idx = None
    
    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        symbol = atom.GetSymbol()
        
        if symbol == 'O':
            n_oxygen += 1
            ring_oxygen_idx = atom_idx
        elif symbol == 'C':
            n_carbon += 1
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                if n_idx not in ring_atoms and neighbor.GetSymbol() == 'O':
                    bond = mol.GetBondBetweenAtoms(atom_idx, n_idx)
                    if bond and bond.GetBondType() == Chem.BondType.SINGLE:
                        n_hydroxyl += 1
        elif symbol == 'N':
            n_nitrogen += 1
        else:
            return False
    
    if ring_oxygen_idx is not None:
        ring_o_atom = mol.GetAtomWithIdx(ring_oxygen_idx)
        for neighbor in ring_o_atom.GetNeighbors():
            if neighbor.GetIdx() in ring_atoms and neighbor.GetSymbol() == 'C':
                for nn in neighbor.GetNeighbors():
                    if nn.GetIdx() not in ring_atoms and nn.GetSymbol() == 'O':
                        has_anomeric_oxygen = True
                        break
    
    valid_carbon = (n_carbon == 5) or (n_carbon == 4 and n_nitrogen == 1)
    valid_oxygen = (n_oxygen == 1)
    valid_hydroxyls = (n_hydroxyl >= 1) or has_anomeric_oxygen
    
    return valid_oxygen and valid_carbon and valid_hydroxyls


def _is_furanose_ring(mol: Chem.Mol, ring_atoms: List[int]) -> bool:
    """Vérifie si un cycle de 5 membres est un furanose."""
    if len(ring_atoms) != 5:
        return False
    
    n_oxygen = 0
    n_carbon = 0
    n_hydroxyl = 0
    has_anomeric_oxygen = False
    ring_oxygen_idx = None
    
    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        symbol = atom.GetSymbol()
        
        if symbol == 'O':
            n_oxygen += 1
            ring_oxygen_idx = atom_idx
        elif symbol == 'C':
            n_carbon += 1
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                if n_idx not in ring_atoms and neighbor.GetSymbol() == 'O':
                    bond = mol.GetBondBetweenAtoms(atom_idx, n_idx)
                    if bond and bond.GetBondType() == Chem.BondType.SINGLE:
                        n_hydroxyl += 1
        else:
            return False
    
    if ring_oxygen_idx is not None:
        ring_o_atom = mol.GetAtomWithIdx(ring_oxygen_idx)
        for neighbor in ring_o_atom.GetNeighbors():
            if neighbor.GetIdx() in ring_atoms and neighbor.GetSymbol() == 'C':
                for nn in neighbor.GetNeighbors():
                    if nn.GetIdx() not in ring_atoms and nn.GetSymbol() == 'O':
                        has_anomeric_oxygen = True
                        break
    
    return (n_oxygen == 1) and (n_carbon == 4) and ((n_hydroxyl >= 1) or has_anomeric_oxygen)


def _identify_anomeric_carbon(mol: Chem.Mol, ring_atoms: List[int], ring_oxygen: int) -> Optional[int]:
    """Identifie le carbone anomérique (C1)."""
    oxygen_atom = mol.GetAtomWithIdx(ring_oxygen)
    adjacent_carbons = []
    
    for neighbor in oxygen_atom.GetNeighbors():
        if neighbor.GetIdx() in ring_atoms and neighbor.GetSymbol() == 'C':
            adjacent_carbons.append(neighbor.GetIdx())
    
    if not adjacent_carbons:
        return None
    
    if len(adjacent_carbons) == 1:
        return adjacent_carbons[0]
    
    for c_idx in adjacent_carbons:
        c_atom = mol.GetAtomWithIdx(c_idx)
        for neighbor in c_atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in ring_atoms and neighbor.GetSymbol() == 'O':
                bond = mol.GetBondBetweenAtoms(c_idx, n_idx)
                if bond and bond.GetBondType() == Chem.BondType.SINGLE:
                    return c_idx
    
    return adjacent_carbons[0]


def _assign_ring_positions(mol: Chem.Mol, ring_atoms: List[int], ring_oxygen: int, anomeric_carbon: Optional[int]) -> Dict[int, int]:
    """Assigne les positions (1, 2, 3, ...) aux atomes du cycle."""
    n = len(ring_atoms)
    
    if anomeric_carbon is None:
        o_idx = ring_atoms.index(ring_oxygen)
        positions = {}
        for i in range(n):
            atom_idx = ring_atoms[(o_idx + 1 + i) % n]
            if atom_idx == ring_oxygen:
                positions[atom_idx] = n
            else:
                positions[atom_idx] = i + 1
        return positions
    
    start_idx = ring_atoms.index(anomeric_carbon)
    ordered = ring_atoms[start_idx:] + ring_atoms[:start_idx]
    
    oxygen_pos = ordered.index(ring_oxygen)
    expected_o_pos = n - 1
    
    if oxygen_pos != expected_o_pos and oxygen_pos != 1:
        ordered = [ordered[0]] + ordered[1:][::-1]
    
    return {atom_idx: i + 1 for i, atom_idx in enumerate(ordered)}


def find_sugar_rings(mol: Chem.Mol) -> List[SugarRing]:
    """Trouve tous les cycles de sucre dans la molécule."""
    ring_info = mol.GetRingInfo()
    sugar_rings = []
    processed_atoms: Set[int] = set()
    
    for ring in ring_info.AtomRings():
        ring_list = list(ring)
        
        if any(a in processed_atoms for a in ring_list):
            continue
        
        ring_type = RingType.UNKNOWN
        
        if len(ring_list) == 6 and _is_pyranose_ring(mol, ring_list):
            ring_type = RingType.PYRANOSE
        elif len(ring_list) == 5 and _is_furanose_ring(mol, ring_list):
            ring_type = RingType.FURANOSE
        
        if ring_type != RingType.UNKNOWN:
            ring_oxygen = None
            for atom_idx in ring_list:
                if mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'O':
                    ring_oxygen = atom_idx
                    break
            
            if ring_oxygen is None:
                continue
            
            anomeric_carbon = _identify_anomeric_carbon(mol, ring_list, ring_oxygen)
            atom_positions = _assign_ring_positions(mol, ring_list, ring_oxygen, anomeric_carbon)
            
            sugar_ring = SugarRing(
                ring_type=ring_type,
                ring_atoms=ring_list,
                ring_oxygen=ring_oxygen,
                anomeric_carbon=anomeric_carbon,
                atom_positions=atom_positions
            )
            
            sugar_rings.append(sugar_ring)
            processed_atoms.update(ring_list)
    
    return sugar_rings


# =============================================================================
# SECTION 5: Identification du type de sucre
# =============================================================================

def _identify_sugar_type(mol: Chem.Mol, sugar_ring: SugarRing) -> str:
    """Identifie le type de sucre par analyse de la stéréochimie."""
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    
    ring_chirality: Dict[int, str] = {}
    for center_idx, chirality in chiral_centers:
        if center_idx in sugar_ring.ring_atoms:
            pos = sugar_ring.atom_positions.get(center_idx)
            if pos:
                ring_chirality[pos] = chirality
    
    aldohexose_patterns = {
        'glucose': {2: 'R', 3: 'S', 4: 'R', 5: 'R'},
        'galactose': {2: 'R', 3: 'S', 4: 'S', 5: 'R'},
        'mannose': {2: 'S', 3: 'S', 4: 'R', 5: 'R'},
        'allose': {2: 'R', 3: 'R', 4: 'R', 5: 'R'},
        'altrose': {2: 'S', 3: 'R', 4: 'R', 5: 'R'},
        'gulose': {2: 'R', 3: 'R', 4: 'S', 5: 'R'},
        'idose': {2: 'S', 3: 'R', 4: 'S', 5: 'R'},
        'talose': {2: 'S', 3: 'S', 4: 'S', 5: 'R'},
    }
    
    furanose_patterns = {
        'ribose': {2: 'R', 3: 'S', 4: 'R'},
        'arabinose': {2: 'S', 3: 'S', 4: 'R'},
        'xylose': {2: 'R', 3: 'S', 4: 'S'},
        'lyxose': {2: 'S', 3: 'S', 4: 'S'},
        'deoxyribose': {3: 'S', 4: 'R'},
    }
    
    patterns = aldohexose_patterns if sugar_ring.ring_type == RingType.PYRANOSE else furanose_patterns
    
    best_match = 'unknown'
    best_score = 0
    
    for sugar_name, pattern in patterns.items():
        score = 0
        total = len(pattern)
        
        for pos, expected_chir in pattern.items():
            if pos in ring_chirality and ring_chirality[pos] == expected_chir:
                score += 1
        
        if score > best_score and score >= total - 1:
            best_score = score
            best_match = sugar_name
    
    return best_match


# =============================================================================
# SECTION 6: Analyse conformationnelle
# =============================================================================

def determine_conformation_from_3d(mol: Chem.Mol, sugar_ring: SugarRing, conf_id: int = 0) -> Conformation:
    """Détermine la conformation à partir de la géométrie 3D."""
    if mol.GetNumConformers() == 0:
        return Conformation.UNKNOWN
    
    conf = mol.GetConformer(conf_id)
    ring_atoms = sugar_ring.ring_atoms
    
    positions = [np.array(conf.GetAtomPosition(i)) for i in ring_atoms]
    centroid = np.mean(positions, axis=0)
    
    v1 = positions[1] - positions[0]
    v2 = positions[2] - positions[0]
    normal = np.cross(v1, v2)
    norm = np.linalg.norm(normal)
    
    if norm < 1e-6:
        return Conformation.UNKNOWN
    
    normal = normal / norm
    deviations = [np.dot(pos - centroid, normal) for pos in positions]
    
    positive = sum(1 for d in deviations if d > 0.2)
    negative = sum(1 for d in deviations if d < -0.2)
    
    if sugar_ring.ring_type == RingType.PYRANOSE:
        if positive >= 2 and negative >= 2:
            max_dev_idx = max(range(len(deviations)), key=lambda i: abs(deviations[i]))
            return Conformation.CHAIR_4C1 if deviations[max_dev_idx] > 0 else Conformation.CHAIR_1C4
        elif (positive >= 4 and negative <= 2) or (negative >= 4 and positive <= 2):
            return Conformation.BOAT_B14
        return Conformation.CHAIR_4C1
    
    elif sugar_ring.ring_type == RingType.FURANOSE:
        max_deviation = 0
        out_of_plane_idx = 0
        
        for i in range(len(positions)):
            other_pos = [p for j, p in enumerate(positions) if j != i]
            if len(other_pos) < 3:
                continue
            
            v1 = other_pos[1] - other_pos[0]
            v2 = other_pos[2] - other_pos[0]
            n = np.cross(v1, v2)
            n_norm = np.linalg.norm(n)
            if n_norm < 1e-6:
                continue
            n = n / n_norm
            
            plane_centroid = np.mean(other_pos, axis=0)
            deviation = abs(np.dot(positions[i] - plane_centroid, n))
            
            if deviation > max_deviation:
                max_deviation = deviation
                out_of_plane_idx = i
        
        pos = sugar_ring.atom_positions.get(ring_atoms[out_of_plane_idx], 0)
        
        if max_deviation > 0.3:
            if pos == 2:
                return Conformation.C2_ENDO
            elif pos == 3:
                return Conformation.C3_ENDO
        
        return Conformation.TWIST_T
    
    return Conformation.UNKNOWN


def determine_conformation_fallback(mol: Chem.Mol, sugar_ring: SugarRing) -> Tuple[Conformation, Dict[int, ProtonOrientation]]:
    """Détermine la conformation SANS géométrie 3D via stéréochimie."""
    proton_orientations: Dict[int, ProtonOrientation] = {}
    
    sugar_type = _identify_sugar_type(mol, sugar_ring)
    sugar_ring.sugar_type = sugar_type
    
    chiral_centers = {c: cfg for c, cfg in Chem.FindMolChiralCenters(mol, includeUnassigned=True)}
    
    if sugar_ring.ring_type == RingType.PYRANOSE:
        c5_idx = None
        for atom_idx, pos in sugar_ring.atom_positions.items():
            if pos == 5:
                c5_idx = atom_idx
                break
        
        is_D_sugar = True
        if c5_idx is not None and c5_idx in chiral_centers:
            is_D_sugar = (chiral_centers[c5_idx] == 'R')
        
        if is_D_sugar:
            conformation = Conformation.CHAIR_1C4 if sugar_type == 'idose' else Conformation.CHAIR_4C1
        else:
            conformation = Conformation.CHAIR_1C4
        
        for atom_idx, pos in sugar_ring.atom_positions.items():
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() != 'C':
                continue
            
            h_indices = [n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == 'H']
            if not h_indices:
                continue
            
            chirality = chiral_centers.get(atom_idx)
            orient = ProtonOrientation.UNKNOWN
            
            if conformation == Conformation.CHAIR_4C1:
                if pos == 1:
                    orient = ProtonOrientation.EQUATORIAL if chirality == 'S' else (ProtonOrientation.AXIAL if chirality == 'R' else ProtonOrientation.UNKNOWN)
                elif pos in [2, 4]:
                    orient = ProtonOrientation.AXIAL if chirality == 'R' else (ProtonOrientation.EQUATORIAL if chirality == 'S' else ProtonOrientation.UNKNOWN)
                elif pos in [3, 5]:
                    orient = ProtonOrientation.EQUATORIAL if chirality == 'R' else (ProtonOrientation.AXIAL if chirality == 'S' else ProtonOrientation.UNKNOWN)
            elif conformation == Conformation.CHAIR_1C4:
                if pos == 1:
                    orient = ProtonOrientation.AXIAL if chirality == 'S' else (ProtonOrientation.EQUATORIAL if chirality == 'R' else ProtonOrientation.UNKNOWN)
                elif pos in [2, 4]:
                    orient = ProtonOrientation.EQUATORIAL if chirality == 'R' else (ProtonOrientation.AXIAL if chirality == 'S' else ProtonOrientation.UNKNOWN)
                elif pos in [3, 5]:
                    orient = ProtonOrientation.AXIAL if chirality == 'R' else (ProtonOrientation.EQUATORIAL if chirality == 'S' else ProtonOrientation.UNKNOWN)
            
            for h_idx in h_indices:
                proton_orientations[h_idx] = orient
    
    elif sugar_ring.ring_type == RingType.FURANOSE:
        if sugar_type in ['ribose', 'deoxyribose']:
            conformation = Conformation.C2_ENDO
        elif sugar_type == 'fructose':
            conformation = Conformation.C2_EXO
        elif sugar_type == 'arabinose':
            conformation = Conformation.C3_ENDO
        else:
            conformation = Conformation.TWIST_T
        
        for atom_idx, pos in sugar_ring.atom_positions.items():
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() != 'C':
                continue
            h_indices = [n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == 'H']
            for h_idx in h_indices:
                proton_orientations[h_idx] = ProtonOrientation.PSEUDO_EQUATORIAL
    else:
        conformation = Conformation.UNKNOWN
    
    return conformation, proton_orientations


def determine_proton_orientation_3d(mol: Chem.Mol, h_idx: int, sugar_ring: SugarRing, conf_id: int = 0) -> ProtonOrientation:
    """Détermine l'orientation d'un proton via géométrie 3D."""
    if mol.GetNumConformers() == 0:
        return ProtonOrientation.UNKNOWN
    
    conf = mol.GetConformer(conf_id)
    
    h_atom = mol.GetAtomWithIdx(h_idx)
    neighbors = list(h_atom.GetNeighbors())
    if not neighbors:
        return ProtonOrientation.UNKNOWN
    
    c_atom = neighbors[0]
    c_idx = c_atom.GetIdx()
    
    if c_idx not in sugar_ring.ring_atoms:
        return ProtonOrientation.UNKNOWN
    
    positions = [np.array(conf.GetAtomPosition(i)) for i in sugar_ring.ring_atoms]
    centroid = np.mean(positions, axis=0)
    
    v1 = positions[1] - positions[0]
    v2 = positions[2] - positions[0]
    normal = np.cross(v1, v2)
    norm = np.linalg.norm(normal)
    
    if norm < 1e-6:
        return ProtonOrientation.UNKNOWN
    
    normal = normal / norm
    
    c_pos = np.array(conf.GetAtomPosition(c_idx))
    h_pos = np.array(conf.GetAtomPosition(h_idx))
    ch_vector = h_pos - c_pos
    ch_norm = np.linalg.norm(ch_vector)
    
    if ch_norm < 1e-6:
        return ProtonOrientation.UNKNOWN
    
    ch_vector = ch_vector / ch_norm
    cos_angle = abs(np.dot(ch_vector, normal))
    
    if sugar_ring.ring_type == RingType.PYRANOSE:
        if cos_angle > 0.75:
            return ProtonOrientation.AXIAL
        elif cos_angle < 0.35:
            return ProtonOrientation.EQUATORIAL
        return ProtonOrientation.UNKNOWN
    
    elif sugar_ring.ring_type == RingType.FURANOSE:
        if cos_angle > 0.6:
            return ProtonOrientation.PSEUDO_AXIAL
        elif cos_angle < 0.4:
            return ProtonOrientation.PSEUDO_EQUATORIAL
        return ProtonOrientation.UNKNOWN
    
    return ProtonOrientation.UNKNOWN


def get_proton_orientation_smart(mol: Chem.Mol, h_idx: int, sugar_ring: SugarRing, conf_id: int = 0) -> ProtonOrientation:
    """Obtient l'orientation avec fallback intelligent."""
    if h_idx in sugar_ring.proton_orientations:
        cached = sugar_ring.proton_orientations[h_idx]
        if cached != ProtonOrientation.UNKNOWN:
            return cached
    
    if mol.GetNumConformers() > 0:
        orient = determine_proton_orientation_3d(mol, h_idx, sugar_ring, conf_id)
        if orient != ProtonOrientation.UNKNOWN:
            sugar_ring.proton_orientations[h_idx] = orient
            return orient
    
    return sugar_ring.proton_orientations.get(h_idx, ProtonOrientation.UNKNOWN)


# =============================================================================
# SECTION 7: Configuration anomérique
# =============================================================================

def determine_anomeric_configuration(mol: Chem.Mol, sugar_ring: SugarRing, conf_id: int = 0) -> AnomericConfig:
    """Détermine si le sucre est α ou β."""
    if sugar_ring.anomeric_carbon is None:
        return AnomericConfig.UNKNOWN
    
    c1 = sugar_ring.anomeric_carbon
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
    for center_idx, chirality in chiral_centers:
        if center_idx == c1:
            if chirality == 'R':
                return AnomericConfig.BETA
            elif chirality == 'S':
                return AnomericConfig.ALPHA
    
    if mol.GetNumConformers() > 0:
        c1_atom = mol.GetAtomWithIdx(c1)
        h_anom = None
        
        for neighbor in c1_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'H':
                h_anom = neighbor.GetIdx()
                break
        
        if h_anom is not None:
            orient = determine_proton_orientation_3d(mol, h_anom, sugar_ring, conf_id)
            
            if sugar_ring.conformation in [Conformation.CHAIR_4C1, Conformation.UNKNOWN]:
                if orient == ProtonOrientation.AXIAL:
                    return AnomericConfig.BETA
                elif orient == ProtonOrientation.EQUATORIAL:
                    return AnomericConfig.ALPHA
            elif sugar_ring.conformation == Conformation.CHAIR_1C4:
                if orient == ProtonOrientation.AXIAL:
                    return AnomericConfig.ALPHA
                elif orient == ProtonOrientation.EQUATORIAL:
                    return AnomericConfig.BETA
    
    return AnomericConfig.UNKNOWN


# =============================================================================
# SECTION 8: Calcul des couplages
# =============================================================================

def _get_J_from_orientations(orient1: ProtonOrientation, orient2: ProtonOrientation, ring_type: RingType) -> float:
    """Obtient J à partir des orientations."""
    table = _PYRANOSE_3J if ring_type == RingType.PYRANOSE else _FURANOSE_3J
    
    key = (orient1.value, orient2.value)
    if key in table:
        return table[key]
    
    key_rev = (orient2.value, orient1.value)
    if key_rev in table:
        return table[key_rev]
    
    if orient1 == ProtonOrientation.AXIAL and orient2 == ProtonOrientation.AXIAL:
        return 9.5
    elif orient1 == ProtonOrientation.EQUATORIAL and orient2 == ProtonOrientation.EQUATORIAL:
        return 2.5
    elif ProtonOrientation.AXIAL in [orient1, orient2] and ProtonOrientation.EQUATORIAL in [orient1, orient2]:
        return 3.8
    
    return 5.0 if ring_type == RingType.PYRANOSE else 4.0


def get_anomeric_coupling(mol: Chem.Mol, sugar_ring: SugarRing) -> float:
    """Retourne le couplage J1,2 anomérique."""
    config = sugar_ring.anomeric_config
    conformation = sugar_ring.conformation
    
    if sugar_ring.ring_type == RingType.PYRANOSE:
        table = _ANOMERIC_J['pyranose']
        if conformation in [Conformation.CHAIR_4C1, Conformation.UNKNOWN]:
            if config == AnomericConfig.ALPHA:
                return table['alpha_4C1']
            elif config == AnomericConfig.BETA:
                return table['beta_4C1']
            return table['unknown']
        elif conformation == Conformation.CHAIR_1C4:
            if config == AnomericConfig.ALPHA:
                return table['alpha_1C4']
            elif config == AnomericConfig.BETA:
                return table['beta_1C4']
    
    elif sugar_ring.ring_type == RingType.FURANOSE:
        table = _ANOMERIC_J['furanose']
        if config == AnomericConfig.ALPHA:
            return table['alpha']
        elif config == AnomericConfig.BETA:
            return table['beta']
        return table['unknown']
    
    return 5.5


def get_geminal_coupling(sugar_ring: SugarRing, position: int) -> float:
    """Retourne le couplage 2J géminal pour un CH2."""
    if sugar_ring.ring_type == RingType.PYRANOSE:
        return _GEMINAL_2J['pyranose_CH2']
    elif sugar_ring.ring_type == RingType.FURANOSE:
        return _GEMINAL_2J['furanose_CH2']
    return 12.0


# =============================================================================
# SECTION 9: Multiplicité
# =============================================================================

def _determine_sugar_multiplicity(J_vals: List[float], coupling_partners: List[int], sugar_ring: SugarRing) -> str:
    """Détermine la multiplicité correcte pour les sucres."""
    if not J_vals or not coupling_partners:
        return 's'
    
    n_couplings = len(J_vals)
    
    if n_couplings == 1:
        return 'd'
    
    J_groups: List[List[float]] = []
    used = [False] * len(J_vals)
    
    for i, J in enumerate(J_vals):
        if used[i]:
            continue
        group = [J]
        used[i] = True
        for j, J2 in enumerate(J_vals):
            if not used[j] and abs(J - J2) < 1.0:
                group.append(J2)
                used[j] = True
        J_groups.append(group)
    
    n_groups = len(J_groups)
    group_sizes = sorted([len(g) for g in J_groups], reverse=True)
    
    if n_couplings == 2:
        return 't' if n_groups == 1 else 'dd'
    elif n_couplings == 3:
        if n_groups == 1:
            return 'q'
        elif n_groups == 2:
            return 'dt' if group_sizes == [2, 1] else 'ddd'
        return 'ddd'
    elif n_couplings == 4:
        if n_groups == 1:
            return 'quint'
        elif n_groups == 2:
            if group_sizes == [2, 2]:
                return 'tt'
            elif group_sizes == [3, 1]:
                return 'dq'
            return 'dddd'
        elif n_groups == 3:
            return 'ddt' if 2 in group_sizes else 'dddd'
        return 'dddd'
    elif n_couplings >= 5:
        return 'm'
    
    return 'm'


def _get_expected_multiplicity(sugar_ring: SugarRing, position: int) -> Optional[str]:
    """Retourne la multiplicité attendue pour validation."""
    ring_key = 'pyranose' if sugar_ring.ring_type == RingType.PYRANOSE else 'furanose'
    
    if ring_key not in _TYPICAL_SUGAR_MULTIPLICITIES:
        return None
    
    pos_table = _TYPICAL_SUGAR_MULTIPLICITIES[ring_key]
    
    if position not in pos_table:
        return None
    
    mult_entry = pos_table[position]
    
    if isinstance(mult_entry, dict):
        if position == 1:
            config_key = sugar_ring.anomeric_config.value if sugar_ring.anomeric_config != AnomericConfig.UNKNOWN else 'unknown'
            return mult_entry.get(config_key, mult_entry.get('default', 'd'))
        return mult_entry.get(sugar_ring.sugar_type, mult_entry.get('default', 'dd'))
    
    return mult_entry


# =============================================================================
# SECTION 10: Détection des modifications
# =============================================================================

def detect_sugar_modifications(mol: Chem.Mol, sugar_ring: SugarRing) -> SugarModification:
    """Détecte les modifications chimiques sur le sucre."""
    modifications = SugarModification()
    
    for atom_idx in sugar_ring.ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != 'C':
            continue
        
        position = sugar_ring.atom_positions.get(atom_idx)
        if position is None:
            continue
        
        has_hydroxyl = False
        has_amino = False
        has_acetyl = False
        has_phosphate = False
        has_methyl = False
        
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx in sugar_ring.ring_atoms:
                continue
            
            symbol = neighbor.GetSymbol()
            
            if symbol == 'O':
                bond = mol.GetBondBetweenAtoms(atom_idx, n_idx)
                if bond and bond.GetBondType() == Chem.BondType.SINGLE:
                    o_neighbors = [n for n in neighbor.GetNeighbors() if n.GetIdx() != atom_idx]
                    
                    if not o_neighbors:
                        has_hydroxyl = True
                    elif len(o_neighbors) == 1:
                        on = o_neighbors[0]
                        if on.GetSymbol() == 'H':
                            has_hydroxyl = True
                        elif on.GetSymbol() == 'P':
                            has_phosphate = True
                        elif on.GetSymbol() == 'C':
                            for cc in on.GetNeighbors():
                                if cc.GetIdx() not in [n_idx, neighbor.GetIdx()]:
                                    if cc.GetSymbol() == 'O':
                                        cb = mol.GetBondBetweenAtoms(on.GetIdx(), cc.GetIdx())
                                        if cb and cb.GetBondType() == Chem.BondType.DOUBLE:
                                            has_acetyl = True
                                            break
                            if not has_acetyl:
                                has_methyl = True
            elif symbol == 'N':
                has_amino = True
        
        if not has_hydroxyl and position not in [1, 5]:
            n_H = atom.GetTotalNumHs()
            if n_H >= 2:
                modifications.deoxy_positions.append(position)
        
        if has_amino:
            modifications.amino_positions.append(position)
        if has_acetyl:
            modifications.acetyl_positions.append(position)
        if has_phosphate:
            modifications.phosphate_positions.append(position)
        if has_methyl:
            modifications.methyl_positions.append(position)
    
    return modifications


# =============================================================================
# SECTION 11: Fonction principale
# =============================================================================

def estimate_sugar_couplings_complete(
    mol: Chem.Mol,
    heavy_idx: int,
    h_list: List[int],
    heavy_to_H: Dict[int, List[int]],
    proton_to_heavy: Dict[int, int],
    D: np.ndarray,
    config: Any,
    use_boltzmann: bool = False
) -> Optional[Dict[str, Any]]:
    """
    Fonction principale pour estimer les couplages dans les sucres.
    
    Returns:
        Dict avec J_values, coupling_partners, multiplicity, sugar_info
        ou None si pas un sucre
    """
    sugar_rings = find_sugar_rings(mol)
    
    target_ring: Optional[SugarRing] = None
    for ring in sugar_rings:
        if heavy_idx in ring.ring_atoms:
            target_ring = ring
            break
    
    if target_ring is None:
        return None
    
    # Conformation avec fallback
    if mol.GetNumConformers() > 0:
        target_ring.conformation = determine_conformation_from_3d(mol, target_ring)
        for h_idx in h_list:
            orient = determine_proton_orientation_3d(mol, h_idx, target_ring)
            target_ring.proton_orientations[h_idx] = orient
    
    fallback_conf, fallback_orientations = determine_conformation_fallback(mol, target_ring)
    
    if target_ring.conformation == Conformation.UNKNOWN:
        target_ring.conformation = fallback_conf
    
    for h_idx, orient in fallback_orientations.items():
        if target_ring.proton_orientations.get(h_idx, ProtonOrientation.UNKNOWN) == ProtonOrientation.UNKNOWN:
            target_ring.proton_orientations[h_idx] = orient
    
    target_ring.anomeric_config = determine_anomeric_configuration(mol, target_ring)
    
    if target_ring.sugar_type == "unknown":
        target_ring.sugar_type = _identify_sugar_type(mol, target_ring)
    
    modifications = detect_sugar_modifications(mol, target_ring)
    
    position = target_ring.atom_positions.get(heavy_idx)
    is_anomeric = (heavy_idx == target_ring.anomeric_carbon)
    
    coupling_data: List[Tuple[int, float]] = []
    seen_carbon_pairs: Set[Tuple[int, int]] = set()
    
    # Couplages vicinaux dans le cycle
    for other_heavy in target_ring.ring_atoms:
        if other_heavy == heavy_idx:
            continue
        
        other_atom = mol.GetAtomWithIdx(other_heavy)
        if other_atom.GetSymbol() not in ['C', 'N']:
            continue
        
        dist = int(D[heavy_idx, other_heavy])
        if dist != 1:
            continue
        
        pair_key = tuple(sorted([heavy_idx, other_heavy]))
        if pair_key in seen_carbon_pairs:
            continue
        seen_carbon_pairs.add(pair_key)
        
        other_hs = heavy_to_H.get(other_heavy, [])
        if not other_hs:
            continue
        
        if is_anomeric:
            J = get_anomeric_coupling(mol, target_ring)
        else:
            h_ref = h_list[0]
            other_h_ref = other_hs[0]
            orient1 = get_proton_orientation_smart(mol, h_ref, target_ring)
            orient2 = get_proton_orientation_smart(mol, other_h_ref, target_ring)
            J = _get_J_from_orientations(orient1, orient2, target_ring.ring_type)
        
        for other_h in other_hs:
            coupling_data.append((other_h, round(J, 1)))
    
    # Couplages avec substituants hors cycle
    atom = mol.GetAtomWithIdx(heavy_idx)
    for neighbor in atom.GetNeighbors():
        n_idx = neighbor.GetIdx()
        if n_idx in target_ring.ring_atoms:
            continue
        if neighbor.GetSymbol() != 'C':
            continue
        
        sub_hs = [n.GetIdx() for n in neighbor.GetNeighbors() if n.GetSymbol() == 'H']
        if sub_hs:
            J_sub = 5.5
            for sub_h in sub_hs:
                coupling_data.append((sub_h, round(J_sub, 1)))
    
    # Filtrer par seuil
    min_j = 0.5
    if hasattr(config, 'min_j_detect'):
        if isinstance(config.min_j_detect, dict):
            min_j = config.min_j_detect.get(3, 0.5)
        else:
            min_j = float(config.min_j_detect)
    
    J_vals: List[float] = []
    coupling_partners: List[int] = []
    
    for partner_h, J in coupling_data:
        if J >= min_j:
            J_vals.append(J)
            coupling_partners.append(partner_h)
    
    multiplicity = _determine_sugar_multiplicity(J_vals, coupling_partners, target_ring)
    expected = _get_expected_multiplicity(target_ring, position) if position else None
    
    if expected and multiplicity == 'm' and expected != 'm':
        multiplicity = expected
    
    return {
        'J_values': J_vals,
        'coupling_partners': coupling_partners,
        'multiplicity': multiplicity,
        'sugar_info': {
            'ring_type': target_ring.ring_type.value,
            'conformation': target_ring.conformation.value,
            'anomeric_config': target_ring.anomeric_config.value,
            'sugar_type': target_ring.sugar_type,
            'position': position,
            'is_anomeric': is_anomeric,
            'expected_multiplicity': expected,
            'modifications': {
                'deoxy': modifications.deoxy_positions,
                'amino': modifications.amino_positions,
                'acetyl': modifications.acetyl_positions,
                'phosphate': modifications.phosphate_positions,
            },
        },
        'source': 'sugar_complete'
    }


# =============================================================================
# SECTION 12: Fonctions utilitaires
# =============================================================================

def is_sugar_carbon(mol: Chem.Mol, atom_idx: int) -> bool:
    """Vérifie rapidement si un atome fait partie d'un cycle de sucre."""
    atom = mol.GetAtomWithIdx(atom_idx)
    
    if atom.GetSymbol() != 'C':
        return False
    
    if not atom.IsInRing():
        return False
    
    ring_info = mol.GetRingInfo()
    
    for ring in ring_info.AtomRings():
        if atom_idx in ring:
            ring_list = list(ring)
            if len(ring_list) == 6 and _is_pyranose_ring(mol, ring_list):
                return True
            if len(ring_list) == 5 and _is_furanose_ring(mol, ring_list):
                return True
    
    return False


def get_sugar_ring_for_atom(mol: Chem.Mol, atom_idx: int) -> Optional[SugarRing]:
    """Retourne le SugarRing contenant un atome donné."""
    sugar_rings = find_sugar_rings(mol)
    
    for ring in sugar_rings:
        if atom_idx in ring.ring_atoms:
            return ring
    
    return None


if __name__ == "__main__":
    print("Module sugars_complete.py chargé avec succès.")
