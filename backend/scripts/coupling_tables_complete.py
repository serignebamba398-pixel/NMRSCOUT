# =============================================================================
# coupling_tables_complete.py
# Tables de couplages complètes pour couverture 100%
# =============================================================================

"""
Module de tables empiriques pour couplages J RMN.

Couverture:
- Alcynes (C≡C-H)
- Aldéhydes (CHO)
- Amides et peptides (CONH)
- Phosphines (P-H)
- Boranes (B-H)
- Silanes (Si-H)
- Organométalliques
- Stéroïdes
- Terpènes
- Macrocycles
- Carbohydrates (déjà couvert, référence)
- Nucléosides/Nucléotides

Auteur: [Votre nom]
Date: 2024
"""

from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
from enum import Enum
import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
except ImportError:
    raise ImportError("RDKit requis")

# Import du module peptides spécialisé (évite duplication)
try:
    from peptide_couplings import (
        get_peptide_coupling_robust,
        calculate_3J_HN_HA,
        calculate_3J_HA_HB,
        calculate_2J_geminal,
        KARPLUS_PARAMETERS,
        AMINO_ACID_SPECIFIC_COUPLINGS,
        AminoAcid,
        SecondaryStructure,
    )
    _PEPTIDE_MODULE_AVAILABLE = True
except ImportError:
    _PEPTIDE_MODULE_AVAILABLE = False


# =============================================================================
# SECTION 1: ALCYNES (C≡C-H)
# =============================================================================

_ALKYNE_COUPLING_TABLE: Dict[str, Dict[str, float]] = {
    # Couplages ²J à travers la triple liaison
    '2J_terminal': {
        'C≡C-H': 2.5,           # Alcyne terminal basique
        'HC≡C-CH3': 2.7,        # Propyne
        'HC≡C-Ph': 2.4,         # Phénylacétylène
        'HC≡C-CH2': 2.6,        # Alcyne avec CH₂
        'HC≡C-C(O)': 2.2,       # Alcyne conjugué carbonyle
    },
    
    # Couplages ⁴J propargiliques (H-C-C≡C-H)
    '4J_propargylic': {
        'default': 2.5,          # H-CH₂-C≡C-H typique
        'allenic': 7.0,          # Système allénique H₂C=C=CH₂
        'conjugated': 3.0,       # Conjugué avec aromatique
    },
    
    # Couplages ⁵J
    '5J_alkyne': {
        'extended_W': 1.5,       # W étendu à travers triple liaison
        'default': 0.5,
    },
    
    # Couplages avec substituants
    '3J_alkyne_substituent': {
        'CH3-C≡C': 0.0,          # Pas de couplage (pas de H sur C≡)
        'CH2-C≡C-H': 2.8,        # À travers CH₂
    },
}

# Allènes (cumulènes)
_ALLENE_COUPLING_TABLE: Dict[str, float] = {
    '4J_allenic': 7.0,           # H₂C=C=CH₂ caractéristique
    '5J_allenic': 3.0,           # Couplage étendu
    '2J_allenic_gem': -6.0,      # Geminal sur allène
}


def get_alkyne_coupling(
    mol: Chem.Mol, 
    h1_idx: int, 
    h2_idx: int,
    bond_path: List[int]
) -> Optional[float]:
    """
    Calcule le couplage pour systèmes alcynes.
    
    Args:
        mol: Molécule RDKit
        h1_idx, h2_idx: Indices des protons
        bond_path: Chemin de liaisons entre les protons
    
    Returns:
        J en Hz ou None
    """
    # Détecter si triple liaison dans le chemin
    has_triple = False
    for i in range(len(bond_path) - 1):
        bond = mol.GetBondBetweenAtoms(bond_path[i], bond_path[i+1])
        if bond and bond.GetBondType() == Chem.BondType.TRIPLE:
            has_triple = True
            break
    
    if not has_triple:
        return None
    
    n_bonds = len(bond_path) - 1
    
    if n_bonds == 2:
        # ²J à travers C≡C
        return _ALKYNE_COUPLING_TABLE['2J_terminal'].get('C≡C-H', 2.5)
    
    elif n_bonds == 4:
        # ⁴J propargilique
        # Vérifier si allénique
        is_allenic = _is_allenic_system(mol, bond_path)
        if is_allenic:
            return _ALLENE_COUPLING_TABLE['4J_allenic']
        return _ALKYNE_COUPLING_TABLE['4J_propargylic'].get('default', 2.5)
    
    elif n_bonds == 5:
        return _ALKYNE_COUPLING_TABLE['5J_alkyne'].get('default', 0.5)
    
    return None


def _is_allenic_system(mol: Chem.Mol, bond_path: List[int]) -> bool:
    """Détecte un système allénique (C=C=C)."""
    for i in range(len(bond_path) - 2):
        a1, a2, a3 = bond_path[i], bond_path[i+1], bond_path[i+2]
        
        bond1 = mol.GetBondBetweenAtoms(a1, a2)
        bond2 = mol.GetBondBetweenAtoms(a2, a3)
        
        if (bond1 and bond2 and 
            bond1.GetBondType() == Chem.BondType.DOUBLE and
            bond2.GetBondType() == Chem.BondType.DOUBLE):
            
            # Vérifier que l'atome central a hybridation sp
            central = mol.GetAtomWithIdx(a2)
            if central.GetHybridization() == Chem.HybridizationType.SP:
                return True
    
    return False


# =============================================================================
# SECTION 2: ALDÉHYDES (CHO)
# =============================================================================

_ALDEHYDE_COUPLING_TABLE: Dict[str, Dict[str, float]] = {
    # ³J aldéhyde-αH
    '3J_CHO_alphaH': {
        'default': 2.8,              # CHO-CH₂ basique
        'branched': 2.5,             # CHO-CH<
        'conjugated': 3.2,           # CHO conjugué
        'aromatic': 0.0,             # Benzaldéhyde (pas de αH)
    },
    
    # ⁴J aldéhyde-βH
    '4J_CHO_betaH': {
        'trans_planar': 1.5,         # Configuration W
        'gauche': 0.5,
        'default': 0.8,
    },
    
    # Couplages dans aldéhydes aromatiques
    'aromatic_aldehyde': {
        '4J_CHO_ortho': 0.8,         # Couplage avec H ortho
        '5J_CHO_meta': 0.3,          # Couplage avec H méta (faible)
    },
    
    # Aldéhydes α,β-insaturés
    'unsaturated_aldehyde': {
        '3J_CHO_vinyl_cis': 7.5,     # CHO-CH=CH cis
        '3J_CHO_vinyl_trans': 15.5,  # CHO-CH=CH trans
        '4J_CHO_vinyl_beta': 1.2,    # Couplage longue portée
    },
}


def get_aldehyde_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    distance: int
) -> Optional[float]:
    """
    Calcule le couplage impliquant un proton aldéhydique.
    
    Note: Pour les aldéhydes aromatiques, les couplages 4J et 5J avec
    les H du cycle sont très faibles (<1 Hz) et généralement non résolus.
    Le CHO apparaît typiquement comme un singulet (éventuellement élargi).
    """
    # Identifier le proton aldéhydique et le carbone du CHO
    aldehyde_h = None
    aldehyde_c = None
    other_h = None
    
    for h_idx in [h1_idx, h2_idx]:
        h_atom = mol.GetAtomWithIdx(h_idx)
        for neighbor in h_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                # Vérifier si C=O (aldéhyde)
                for nn in neighbor.GetNeighbors():
                    if nn.GetSymbol() == 'O':
                        bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nn.GetIdx())
                        if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                            aldehyde_h = h_idx
                            aldehyde_c = neighbor
                            other_h = h2_idx if h_idx == h1_idx else h1_idx
                            break
    
    if aldehyde_h is None or aldehyde_c is None:
        return None
    
    # Vérifier si c'est un aldéhyde aromatique (attaché à cycle aromatique)
    is_aromatic_aldehyde = False
    for neighbor in aldehyde_c.GetNeighbors():
        if neighbor.GetSymbol() != 'O' and neighbor.GetIsAromatic():
            is_aromatic_aldehyde = True
            break
    
    # Pour les aldéhydes aromatiques, ne PAS retourner de couplages
    # car le CHO apparaît comme un singulet (les 4J/5J sont trop faibles)
    if is_aromatic_aldehyde:
        return None
    
    # Déterminer le type de couplage pour les autres aldéhydes
    if distance == 3:
        # ³J avec α-H (aldéhydes aliphatiques)
        return _ALDEHYDE_COUPLING_TABLE['3J_CHO_alphaH'].get('default', 2.8)
    
    elif distance == 4:
        # ⁴J avec β-H - seulement si configuration W favorable
        j_val = _ALDEHYDE_COUPLING_TABLE['4J_CHO_betaH'].get('default', 0.8)
        # Ne retourner que si >= 1.0 Hz (seuil de détection typique)
        return j_val if j_val >= 1.0 else None
    
    return None


# =============================================================================
# SECTION 3: AMIDES ET PEPTIDES
# =============================================================================
# NOTE: La logique complète des peptides est dans peptide_couplings.py
# Ici on garde seulement une table simplifiée et un wrapper
# =============================================================================

# Table simplifiée pour fallback si peptide_couplings.py n'est pas disponible
_AMIDE_COUPLING_TABLE_FALLBACK: Dict[str, Dict[str, float]] = {
    '3J_NH_alphaH': {
        'helix_alpha': 4.0,
        'sheet_beta': 8.5,
        'random_coil': 7.0,
        'default': 6.5,
    },
    '3J_alphaH_betaH': {
        'gauche': 3.5,
        'trans': 12.0,
        'average': 7.0,
        'default': 6.5,
    },
    '2J_geminal': {
        'glycine': -17.5,
        'default': -12.5,
    },
}


def get_peptide_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    dihedral_angle: Optional[float] = None,
    secondary_structure: str = 'random_coil'
) -> Optional[float]:
    """
    Calcule le couplage dans un peptide.
    
    Utilise peptide_couplings.py si disponible, sinon fallback simplifié.
    
    Args:
        mol: Molécule RDKit
        h1_idx, h2_idx: Indices des protons
        dihedral_angle: Angle dièdre en degrés (si connu)
        secondary_structure: 'helix_alpha', 'sheet_beta', 'random_coil'
    
    Returns:
        J en Hz ou None
    """
    # Utiliser le module peptide robuste si disponible
    if _PEPTIDE_MODULE_AVAILABLE:
        try:
            return get_peptide_coupling_robust(
                mol, h1_idx, h2_idx,
                phi_angle=dihedral_angle,
                secondary_structure=secondary_structure
            )
        except Exception:
            pass  # Fallback ci-dessous
    
    # Fallback simplifié si le module n'est pas disponible
    h1 = mol.GetAtomWithIdx(h1_idx)
    h2 = mol.GetAtomWithIdx(h2_idx)
    
    # Trouver les atomes lourds
    c1 = list(h1.GetNeighbors())[0]
    c2 = list(h2.GetNeighbors())[0]
    
    # Vérifier si c'est un couplage NH-CαH
    is_NH = False
    is_CaH = False
    
    for h_idx, c in [(h1_idx, c1), (h2_idx, c2)]:
        # NH amide
        if c.GetSymbol() == 'N':
            for neighbor in c.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    for nn in neighbor.GetNeighbors():
                        if nn.GetSymbol() == 'O':
                            bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nn.GetIdx())
                            if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                                is_NH = True
                                break
        
        # Cα (attaché à N et C=O)
        elif c.GetSymbol() == 'C':
            has_N = False
            has_CO = False
            for neighbor in c.GetNeighbors():
                if neighbor.GetSymbol() == 'N':
                    has_N = True
                elif neighbor.GetSymbol() == 'C':
                    for nn in neighbor.GetNeighbors():
                        if nn.GetSymbol() == 'O':
                            bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nn.GetIdx())
                            if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                                has_CO = True
            if has_N and has_CO:
                is_CaH = True
    
    if is_NH and is_CaH:
        # Couplage NH-CαH - utiliser Karplus si angle disponible
        if dihedral_angle is not None:
            # Bystrov: ³J = 6.4 cos²(φ-60°) - 1.4 cos(φ-60°) + 1.9
            phi_rad = np.radians(dihedral_angle - 60)
            J = 6.4 * np.cos(phi_rad)**2 - 1.4 * np.cos(phi_rad) + 1.9
            return max(0, J)
        else:
            # Valeur par structure secondaire
            return _AMIDE_COUPLING_TABLE_FALLBACK['3J_NH_alphaH'].get(
                secondary_structure, 6.5
            )
    
    # Couplage Hα-Hβ par défaut
    return _AMIDE_COUPLING_TABLE_FALLBACK['3J_alphaH_betaH'].get('default', 6.5)


# =============================================================================
# SECTION 4: PHOSPHINES (P-H)
# =============================================================================

_PHOSPHINE_COUPLING_TABLE: Dict[str, Dict[str, float]] = {
    # ¹J(P-H) direct
    '1J_PH': {
        'PH3': 188.0,                # Phosphine
        'RPH2': 195.0,               # Phosphine primaire
        'R2PH': 205.0,               # Phosphine secondaire
        'HPO3_2-': 680.0,            # Phosphonate
        'H3PO4': 700.0,              # Acide phosphorique
        'phosphite': 690.0,
        'default': 200.0,
    },
    
    # ²J(P-C-H)
    '2J_PCH': {
        'phosphine_CH3': 2.5,
        'phosphine_CH2': 4.0,
        'phosphonate_CH': 20.0,      # O=P-O-CH
        'phosphate_CH': 8.0,
        'default': 5.0,
    },
    
    # ³J(P-C-C-H) - dépend de Karplus
    '3J_PCCH': {
        'trans': 20.0,               # Dièdre ~180°
        'gauche': 5.0,               # Dièdre ~60°
        'average': 10.0,
        'phospholipid': 7.0,
        'nucleotide': 6.5,
        'default': 8.0,
    },
    
    # ⁴J longue portée
    '4J_PCCCH': {
        'W_arrangement': 3.0,
        'default': 1.0,
    },
}


def get_phosphine_coupling(
    mol: Chem.Mol,
    h_idx: int,
    p_idx: int,
    distance: int
) -> Optional[float]:
    """Calcule le couplage H-P."""
    p_atom = mol.GetAtomWithIdx(p_idx)
    
    # Déterminer le type de phosphore
    n_bonds = p_atom.GetDegree()
    has_double_O = False
    
    for neighbor in p_atom.GetNeighbors():
        if neighbor.GetSymbol() == 'O':
            bond = mol.GetBondBetweenAtoms(p_idx, neighbor.GetIdx())
            if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                has_double_O = True
    
    if distance == 1:
        # ¹J(P-H) direct
        if has_double_O:
            return _PHOSPHINE_COUPLING_TABLE['1J_PH'].get('phosphonate', 690.0)
        else:
            return _PHOSPHINE_COUPLING_TABLE['1J_PH'].get('default', 200.0)
    
    elif distance == 2:
        return _PHOSPHINE_COUPLING_TABLE['2J_PCH'].get('default', 5.0)
    
    elif distance == 3:
        return _PHOSPHINE_COUPLING_TABLE['3J_PCCH'].get('default', 8.0)
    
    elif distance == 4:
        return _PHOSPHINE_COUPLING_TABLE['4J_PCCCH'].get('default', 1.0)
    
    return None


# =============================================================================
# SECTION 5: BORANES (B-H)
# =============================================================================

_BORANE_COUPLING_TABLE: Dict[str, Dict[str, float]] = {
    # ¹J(B-H) - Note: ¹¹B a I=3/2
    '1J_BH': {
        'BH3': 103.0,
        'BH4-': 81.0,                # Borohydride
        'R-BH2': 110.0,
        'R2-BH': 120.0,
        'catechol_borane': 140.0,
        'pinacol_borane': 145.0,
        'default': 110.0,
    },
    
    # ²J(B-C-H)
    '2J_BCH': {
        'default': 5.0,
        'vinyl': 7.0,
    },
    
    # ³J(B-C-C-H)
    '3J_BCCH': {
        'default': 3.0,
    },
}


def get_borane_coupling(
    mol: Chem.Mol,
    h_idx: int,
    b_idx: int,
    distance: int
) -> Optional[float]:
    """Calcule le couplage H-B."""
    if distance == 1:
        return _BORANE_COUPLING_TABLE['1J_BH'].get('default', 110.0)
    elif distance == 2:
        return _BORANE_COUPLING_TABLE['2J_BCH'].get('default', 5.0)
    elif distance == 3:
        return _BORANE_COUPLING_TABLE['3J_BCCH'].get('default', 3.0)
    
    return None


# =============================================================================
# SECTION 6: SILANES (Si-H)
# =============================================================================

_SILANE_COUPLING_TABLE: Dict[str, Dict[str, float]] = {
    # ¹J(Si-H)
    '1J_SiH': {
        'SiH4': 202.5,
        'R-SiH3': 185.0,
        'R2-SiH2': 195.0,
        'R3-SiH': 210.0,
        'chlorosilane': 235.0,
        'alkoxysilane': 220.0,
        'default': 200.0,
    },
    
    # ²J(Si-C-H)
    '2J_SiCH': {
        'SiCH3': 6.5,
        'SiCH2': 5.0,
        'default': 6.0,
    },
    
    # ³J(Si-C-C-H)
    '3J_SiCCH': {
        'trans': 12.0,
        'gauche': 3.0,
        'average': 6.0,
        'default': 6.0,
    },
    
    # ⁴J
    '4J_SiCCCH': {
        'default': 1.0,
    },
    
    # ²⁹Si satellites (4.7% abondance)
    'Si29_satellites': {
        '1J_29Si_H': 200.0,          # Visible comme satellites
    },
}


def get_silane_coupling(
    mol: Chem.Mol,
    h_idx: int,
    si_idx: int,
    distance: int
) -> Optional[float]:
    """Calcule le couplage H-Si."""
    si_atom = mol.GetAtomWithIdx(si_idx)
    n_bonds = si_atom.GetDegree()
    
    if distance == 1:
        # ¹J(Si-H)
        # Détecter le type de silane
        n_H = sum(1 for n in si_atom.GetNeighbors() if n.GetSymbol() == 'H')
        n_Cl = sum(1 for n in si_atom.GetNeighbors() if n.GetSymbol() == 'Cl')
        
        if n_Cl > 0:
            return _SILANE_COUPLING_TABLE['1J_SiH'].get('chlorosilane', 235.0)
        elif n_H == 3:
            return _SILANE_COUPLING_TABLE['1J_SiH'].get('R-SiH3', 185.0)
        elif n_H == 2:
            return _SILANE_COUPLING_TABLE['1J_SiH'].get('R2-SiH2', 195.0)
        else:
            return _SILANE_COUPLING_TABLE['1J_SiH'].get('R3-SiH', 210.0)
    
    elif distance == 2:
        return _SILANE_COUPLING_TABLE['2J_SiCH'].get('default', 6.0)
    
    elif distance == 3:
        return _SILANE_COUPLING_TABLE['3J_SiCCH'].get('default', 6.0)
    
    elif distance == 4:
        return _SILANE_COUPLING_TABLE['4J_SiCCCH'].get('default', 1.0)
    
    return None


# =============================================================================
# SECTION 7: ORGANOMÉTALLIQUES
# =============================================================================

_ORGANOMETALLIC_COUPLING_TABLE: Dict[str, Dict[str, Any]] = {
    # Ferrocène et métallocènes
    'ferrocene': {
        '3J_HH_ring': 2.5,           # Couplage sur cycle Cp
        '4J_HH_across': 1.5,         # Entre cycles
    },
    
    # Complexes carbonyles
    'metal_carbonyl': {
        '2J_M_CO_H': 3.0,            # M-CO couplage (faible)
    },
    
    # Complexes phosphine
    'metal_phosphine': {
        '2J_P_M_H': 5.0,
        '3J_P_M_C_H': 8.0,
    },
    
    # Alkyls métalliques
    'metal_alkyl': {
        '2J_M_CH': 8.0,              # Ex: Ti-CH₃
        '3J_M_C_CH': 6.0,
    },
    
    # Hydrures métalliques
    'metal_hydride': {
        '1J_M_H': {
            'default': 30.0,
            'Pt': 850.0,             # ¹⁹⁵Pt-H très grand
            'Rh': 25.0,
            'Ir': 35.0,
            'Ru': 30.0,
        },
        '2J_H_M_H': -5.0,            # Couplage geminal cis
        '2J_H_M_H_trans': -15.0,     # Trans plus grand
    },
    
    # Couplages η²-complexes (oléfines)
    'olefin_complex': {
        '3J_HH_coordinated': 8.0,    # Oléfine coordonnée
    },
}


def get_organometallic_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    metal_symbol: str
) -> Optional[float]:
    """Calcule le couplage dans complexes organométalliques."""
    # Pour l'instant, retourner une valeur par défaut
    # Une implémentation complète nécessiterait la détection du type de complexe
    
    return 5.0  # Valeur générique


# =============================================================================
# SECTION 8: STÉROÏDES
# =============================================================================

_STEROID_COUPLING_TABLE: Dict[str, Dict[str, float]] = {
    # Jonctions de cycles A/B (cis ou trans)
    'AB_junction': {
        '3J_H5_H6ax_trans': 12.0,    # Trans diaxial
        '3J_H5_H6eq_trans': 4.5,
        '3J_H5_H6ax_cis': 5.0,       # Jonction cis
        '3J_H5_H6eq_cis': 3.5,
    },
    
    # Jonction B/C (toujours trans)
    'BC_junction': {
        '3J_H8_H9': 11.0,            # Trans diaxial
        '3J_H8_H14': 10.5,
    },
    
    # Jonction C/D
    'CD_junction': {
        '3J_H14_H15ax': 12.0,
        '3J_H14_H15eq': 4.0,
    },
    
    # W-couplings caractéristiques
    'W_couplings': {
        '4J_H1ax_H3ax': 2.5,         # W dans cycle A
        '4J_H2ax_H4ax': 2.5,
        '4J_H11ax_H12ax': 2.0,       # W dans cycle C
        '4J_H15ax_H17': 1.5,
    },
    
    # Méthyles angulaires
    'angular_methyl': {
        '4J_CH3_18_H12': 0.8,        # CH₃ C-18 avec H-12
        '4J_CH3_19_H1ax': 0.5,       # CH₃ C-19 avec H-1ax
        '4J_CH3_19_H2ax': 0.5,
    },
    
    # Cycle A aromatique (estrogènes)
    'aromatic_A_ring': {
        '3J_H1_H2': 8.5,             # Ortho aromatique
        '4J_H1_H4': 2.5,             # Meta
        '4J_H2_H4': 2.5,
    },
    
    # Chaîne latérale (cholestérol type)
    'sidechain': {
        '3J_H17_H20': 9.0,           # Trans
        '3J_H20_H22': 6.5,
        '3J_CH3_21_H20': 6.5,        # Méthyle-21 avec H-20
    },
}


def get_steroid_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    ring_system: str,
    coupling_type: str
) -> Optional[float]:
    """
    Calcule les couplages spécifiques aux stéroïdes.
    
    IMPORTANT: Cette fonction ne retourne des valeurs que pour des protons
    confirmés comme étant dans le système de cycles fusionnés du stéroïde.
    
    Args:
        ring_system: 'AB_junction', 'BC_junction', 'CD_junction', etc.
        coupling_type: Clé spécifique dans la table
    """
    # Vérifier que les deux protons sont attachés à des carbones dans des CYCLES
    h1 = mol.GetAtomWithIdx(h1_idx)
    h2 = mol.GetAtomWithIdx(h2_idx)
    
    c1_neighbors = [n for n in h1.GetNeighbors() if n.GetSymbol() == 'C']
    c2_neighbors = [n for n in h2.GetNeighbors() if n.GetSymbol() == 'C']
    
    if not c1_neighbors or not c2_neighbors:
        return None
    
    c1 = c1_neighbors[0]
    c2 = c2_neighbors[0]
    
    # Les deux carbones doivent être dans un cycle pour appliquer les couplages stéroïde
    if not c1.IsInRing() or not c2.IsInRing():
        return None
    
    if ring_system in _STEROID_COUPLING_TABLE:
        return _STEROID_COUPLING_TABLE[ring_system].get(coupling_type)
    
    return None


# =============================================================================
# SECTION 9: TERPÈNES
# =============================================================================

_TERPENE_COUPLING_TABLE: Dict[str, Dict[str, float]] = {
    # Monoterpènes cycliques
    'cyclohexene': {
        '3J_vinyl': 10.0,            # =CH-CH
        '3J_allylic_ax': 5.0,        # =CH-CH₂ axial
        '3J_allylic_eq': 2.5,        # =CH-CH₂ équatorial
        '4J_homoallylic': 1.5,
    },
    
    'pinane': {
        # Système bicyclique [3.1.1]
        '3J_bridge_endo': 9.0,
        '3J_bridge_exo': 6.0,
        '4J_W_gem_dimethyl': 0.5,
    },
    
    'menthane': {
        # Cyclohexane substitué
        '3J_ax_ax': 11.0,
        '3J_ax_eq': 4.0,
        '3J_eq_eq': 3.5,
        '4J_1_3_diaxial': 1.0,
    },
    
    # Sesquiterpènes
    'decalin': {
        '3J_trans_diaxial': 11.5,
        '3J_cis_junction': 4.0,
    },
    
    'cadinane': {
        '3J_H1_H2ax': 12.0,
        '3J_H1_H2eq': 4.5,
    },
    
    # Isoprène unit
    'isoprene_unit': {
        '3J_CH3_vinyl': 1.3,         # CH₃ avec =CH
        '3J_vinyl_vinyl_cis': 10.0,
        '3J_vinyl_vinyl_trans': 16.0,
        '4J_CH3_methylene': 1.0,
    },
    
    # Gem-diméthyle
    'gem_dimethyl': {
        '4J_CH3_CH3': 0.0,           # Pas de couplage
        '4J_CH3_CH2': 0.5,           # Faible
    },
}


def get_terpene_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    skeleton_type: str
) -> Optional[float]:
    """
    Calcule les couplages dans les terpènes.
    
    IMPORTANT: Cette fonction ne retourne des valeurs que pour des protons
    confirmés comme étant dans un cycle terpénique (cyclohexane/cyclopentane).
    Pour les chaînes alkyles simples, retourne None pour laisser le simulateur
    principal calculer les couplages normaux.
    """
    # Vérifier que les deux protons sont attachés à des carbones dans un CYCLE
    # pour éviter d'appliquer les couplages terpène à des chaînes alkyles
    h1 = mol.GetAtomWithIdx(h1_idx)
    h2 = mol.GetAtomWithIdx(h2_idx)
    
    # Obtenir les carbones parents
    c1_neighbors = [n for n in h1.GetNeighbors() if n.GetSymbol() == 'C']
    c2_neighbors = [n for n in h2.GetNeighbors() if n.GetSymbol() == 'C']
    
    if not c1_neighbors or not c2_neighbors:
        return None
    
    c1 = c1_neighbors[0]
    c2 = c2_neighbors[0]
    
    # Les deux carbones doivent être dans un cycle pour appliquer les couplages terpène
    if not c1.IsInRing() or not c2.IsInRing():
        return None
    
    # Les deux carbones doivent être dans le MÊME cycle
    ring_info = mol.GetRingInfo()
    c1_rings = set()
    c2_rings = set()
    
    for i, ring in enumerate(ring_info.AtomRings()):
        if c1.GetIdx() in ring:
            c1_rings.add(i)
        if c2.GetIdx() in ring:
            c2_rings.add(i)
    
    # Pas de cycle commun = pas de couplage terpène
    if not c1_rings.intersection(c2_rings):
        return None
    
    # OK, les deux protons sont dans le même cycle - appliquer la table
    if skeleton_type in _TERPENE_COUPLING_TABLE:
        table = _TERPENE_COUPLING_TABLE[skeleton_type]
        return table.get('3J_ax_ax', 7.0)
    
    return None


# =============================================================================
# SECTION 10: MACROCYCLES
# =============================================================================

_MACROCYCLE_COUPLING_TABLE: Dict[str, Dict[str, float]] = {
    # Macrolactones (lactones >8 membres)
    'macrolactone': {
        '3J_average': 7.0,           # Rotation partielle
        '3J_rigid': 10.0,            # Si rigide
        '3J_flexible': 5.0,          # Très flexible
    },
    
    # Macrolides (érythromycine type)
    'macrolide': {
        '3J_backbone': 6.5,
        '3J_sugar_aglycone': 7.5,
        '2J_CH2': -12.0,
    },
    
    # Lactames macrocycliques
    'macrolactam': {
        '3J_NH_CH': 6.0,
        '3J_backbone': 7.0,
    },
    
    # Éthers couronne
    'crown_ether': {
        '3J_OCH2CH2O_gauche': 4.0,
        '3J_OCH2CH2O_anti': 10.0,
        '3J_average': 5.5,           # Rotation
    },
    
    # Cyclodextrines
    'cyclodextrin': {
        # Similaire aux sucres
        '3J_H1_H2': 3.5,
        '3J_H2_H3': 9.5,
        '3J_H3_H4': 9.0,
        '3J_H4_H5': 9.5,
    },
    
    # Peptides cycliques
    'cyclic_peptide': {
        '3J_NH_alphaH': 7.5,         # Plus rigide
        '3J_alpha_beta': 6.5,
    },
    
    # Depsipeptides
    'depsipeptide': {
        '3J_ester_OCH': 3.5,
        '3J_amide_NCH': 7.0,
    },
    
    # Polyènes macrocycliques
    'macrocyclic_polyene': {
        '3J_vinyl_cis': 11.0,
        '3J_vinyl_trans': 15.0,
        '4J_across_double': 1.0,
    },
}


def get_macrocycle_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    ring_size: int,
    macrocycle_type: str = 'default'
) -> Optional[float]:
    """
    Calcule les couplages dans les macrocycles.
    
    IMPORTANT: Cette fonction ne retourne des valeurs que pour des protons
    confirmés comme étant dans le macrocycle.
    
    Args:
        ring_size: Taille du macrocycle (>8)
        macrocycle_type: Type de macrocycle
    """
    # Vérifier que les deux protons sont attachés à des carbones dans le MACROCYCLE
    h1 = mol.GetAtomWithIdx(h1_idx)
    h2 = mol.GetAtomWithIdx(h2_idx)
    
    c1_neighbors = [n for n in h1.GetNeighbors() if n.GetSymbol() == 'C']
    c2_neighbors = [n for n in h2.GetNeighbors() if n.GetSymbol() == 'C']
    
    if not c1_neighbors or not c2_neighbors:
        return None
    
    c1 = c1_neighbors[0]
    c2 = c2_neighbors[0]
    
    # Les deux carbones doivent être dans un cycle (et dans le même grand cycle)
    if not c1.IsInRing() or not c2.IsInRing():
        return None
    
    # Vérifier qu'ils sont dans un cycle de taille >= 8 (macrocycle)
    ring_info = mol.GetRingInfo()
    in_macrocycle = False
    
    for ring in ring_info.AtomRings():
        if len(ring) >= 8 and c1.GetIdx() in ring and c2.GetIdx() in ring:
            in_macrocycle = True
            break
    
    if not in_macrocycle:
        return None
    
    # Ajustement selon la taille
    flexibility_factor = 1.0
    
    if ring_size <= 10:
        flexibility_factor = 0.9  # Plus rigide
    elif ring_size <= 14:
        flexibility_factor = 1.0
    elif ring_size <= 20:
        flexibility_factor = 1.1  # Plus flexible
    else:
        flexibility_factor = 1.2  # Très flexible
    
    if macrocycle_type in _MACROCYCLE_COUPLING_TABLE:
        base_J = _MACROCYCLE_COUPLING_TABLE[macrocycle_type].get('3J_average', 7.0)
        return base_J * flexibility_factor
    
    # Valeur par défaut pour macrocycles
    return 7.0 * flexibility_factor


# =============================================================================
# SECTION 11: NUCLÉOSIDES ET NUCLÉOTIDES
# =============================================================================

_NUCLEOTIDE_COUPLING_TABLE: Dict[str, Dict[str, float]] = {
    # Ribose (furanose)
    'ribose': {
        # Dépend de la pseudorotation
        '3J_H1_H2_C2endo': 6.5,      # C2'-endo (S-type)
        '3J_H1_H2_C3endo': 1.5,      # C3'-endo (N-type)
        '3J_H2_H3_C2endo': 5.0,
        '3J_H2_H3_C3endo': 5.5,
        '3J_H3_H4_C2endo': 3.0,
        '3J_H3_H4_C3endo': 9.0,
        '3J_H4_H5a': 3.0,
        '3J_H4_H5b': 4.0,
        '2J_H5a_H5b': -12.0,
    },
    
    # Désoxyribose
    'deoxyribose': {
        '3J_H1_H2a': 6.5,
        '3J_H1_H2b': 7.0,
        '2J_H2a_H2b': -14.0,
        '3J_H2a_H3': 6.0,
        '3J_H2b_H3': 3.5,
        '3J_H3_H4': 3.5,
    },
    
    # Couplages avec phosphate
    'phosphate': {
        '3J_P_H3': 8.0,              # ³¹P-H3'
        '3J_P_H5a': 5.0,
        '3J_P_H5b': 3.0,
        '4J_P_H4': 2.0,
    },
    
    # Bases azotées
    'adenine': {
        '3J_H2_H8': 0.0,             # Pas de couplage (séparés)
    },
    
    'guanine': {
        '1J_N1H': 90.0,              # Si détectable
    },
    
    'cytosine': {
        '3J_H5_H6': 7.5,             # Couplage vinylique
    },
    
    'thymine': {
        '4J_CH3_H6': 1.2,            # Couplage allylique
    },
    
    'uracil': {
        '3J_H5_H6': 8.0,
    },
    
    # Couplage glycosidique H1'-base
    'glycosidic': {
        # Dépend de l'angle χ
        '4J_H1_H8_syn': 1.0,         # Conformation syn
        '4J_H1_H8_anti': 0.3,        # Conformation anti
    },
}

# Équations de pseudorotation
def _calculate_ribose_coupling(P: float, coupling_type: str) -> float:
    """
    Calcule le couplage ribose selon la phase de pseudorotation P.
    
    P = 0° → C3'-endo (N-type, ARN)
    P = 180° → C2'-endo (S-type, ADN)
    
    Équation: J = A + B*cos(P + phase) + C*cos(2*(P + phase))
    """
    # Paramètres empiriques (Altona-Sundaralingam)
    params = {
        '3J_H1_H2': {'A': 4.0, 'B': 2.5, 'C': 0.5, 'phase': 0},
        '3J_H2_H3': {'A': 5.2, 'B': 0.0, 'C': 0.5, 'phase': 0},
        '3J_H3_H4': {'A': 6.0, 'B': -3.0, 'C': 0.5, 'phase': 0},
    }
    
    if coupling_type not in params:
        return 5.0
    
    p = params[coupling_type]
    P_rad = np.radians(P + p['phase'])
    
    J = p['A'] + p['B'] * np.cos(P_rad) + p['C'] * np.cos(2 * P_rad)
    
    return max(0, J)


def get_nucleotide_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    sugar_type: str = 'ribose',
    pseudorotation_phase: Optional[float] = None
) -> Optional[float]:
    """Calcule les couplages dans nucléosides/nucléotides."""
    if sugar_type in _NUCLEOTIDE_COUPLING_TABLE:
        table = _NUCLEOTIDE_COUPLING_TABLE[sugar_type]
        
        # Si phase de pseudorotation connue, calculer précisément
        if pseudorotation_phase is not None and sugar_type == 'ribose':
            # Déterminer quel couplage
            return _calculate_ribose_coupling(pseudorotation_phase, '3J_H1_H2')
        
        # Sinon utiliser valeur moyenne
        return table.get('3J_H1_H2_C2endo', 5.0)
    
    return None


# =============================================================================
# SECTION 12: DÉTECTION AUTOMATIQUE ET DISPATCH
# =============================================================================

class MoleculeClassifier:
    """Classifie les molécules pour choisir la bonne table de couplages."""
    
    @staticmethod
    def classify(mol: Chem.Mol) -> List[str]:
        """
        Retourne les types de structures présentes dans la molécule.
        
        Returns:
            Liste de types: ['alkyne', 'aldehyde', 'peptide', ...]
        """
        types = []
        
        # Parcourir les atomes
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            
            # Phosphore
            if symbol == 'P':
                types.append('phosphorus')
            
            # Bore
            elif symbol == 'B':
                types.append('borane')
            
            # Silicium
            elif symbol == 'Si':
                types.append('silane')
            
            # Métaux
            elif symbol in ['Fe', 'Ru', 'Os', 'Co', 'Rh', 'Ir', 'Ni', 'Pd', 'Pt',
                           'Cu', 'Ag', 'Au', 'Zn', 'Ti', 'Zr', 'Hf', 'Mo', 'W']:
                types.append('organometallic')
        
        # Patterns structurels via SMARTS
        smarts_patterns = {
            'alkyne': '[C]#[C]',
            'allene': '[C]=[C]=[C]',
            'aldehyde': '[CH]=O',
            'amide': '[NX3][CX3](=[OX1])',
            'ester': '[OX2][CX3](=[OX1])',
            'macrolactone': '[OR1](=[O])',  # Approximatif
        }
        
        for name, smarts in smarts_patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                types.append(name)
        
        # Détection de cycles
        ring_info = mol.GetRingInfo()
        ring_sizes = [len(r) for r in ring_info.AtomRings()]
        
        if any(s > 8 for s in ring_sizes):
            types.append('macrocycle')
        
        if 'amide' in types and len(types) > 3:
            # Probablement un peptide
            types.append('peptide')
        
        # Stéroïde: 4 cycles fusionnés (6-6-6-5 ou similaire)
        if len(ring_sizes) >= 4:
            fused_count = sum(1 for s in ring_sizes if s in [5, 6])
            if fused_count >= 4:
                types.append('steroid')
        
        # Terpène: patterns isoprène
        isoprene = Chem.MolFromSmarts('[CH3]C(=C)C')
        if isoprene and mol.HasSubstructMatch(isoprene):
            types.append('terpene')
        
        # Nucléoside
        nucleobase = Chem.MolFromSmarts('[nR1]1[cR1][nR1][cR1][cR1]1')  # Purine-like
        if nucleobase and mol.HasSubstructMatch(nucleobase):
            types.append('nucleoside')
        
        return list(set(types))


def get_specialized_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    distance: int,
    bond_path: Optional[List[int]] = None
) -> Optional[float]:
    """
    Fonction principale de dispatch vers les tables spécialisées.
    
    Args:
        mol: Molécule RDKit
        h1_idx, h2_idx: Indices des protons
        distance: Distance topologique (nombre de liaisons)
        bond_path: Chemin d'atomes entre h1 et h2
    
    Returns:
        J en Hz ou None si non applicable
    """
    # Classifier la molécule
    mol_types = MoleculeClassifier.classify(mol)
    
    # Essayer chaque table spécialisée
    
    # 1. Alcynes
    if 'alkyne' in mol_types and bond_path:
        J = get_alkyne_coupling(mol, h1_idx, h2_idx, bond_path)
        if J is not None:
            return J
    
    # 2. Aldéhydes
    if 'aldehyde' in mol_types:
        J = get_aldehyde_coupling(mol, h1_idx, h2_idx, distance)
        if J is not None:
            return J
    
    # 3. Peptides/Amides
    if 'peptide' in mol_types or 'amide' in mol_types:
        J = get_peptide_coupling(mol, h1_idx, h2_idx)
        if J is not None:
            return J
    
    # 4. Phosphore
    if 'phosphorus' in mol_types:
        # Trouver le P
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'P':
                J = get_phosphine_coupling(mol, h1_idx, atom.GetIdx(), distance)
                if J is not None:
                    return J
    
    # 5. Bore
    if 'borane' in mol_types:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'B':
                J = get_borane_coupling(mol, h1_idx, atom.GetIdx(), distance)
                if J is not None:
                    return J
    
    # 6. Silicium
    if 'silane' in mol_types:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'Si':
                J = get_silane_coupling(mol, h1_idx, atom.GetIdx(), distance)
                if J is not None:
                    return J
    
    # 7. Stéroïdes
    if 'steroid' in mol_types:
        J = get_steroid_coupling(mol, h1_idx, h2_idx, 'AB_junction', '3J_H5_H6ax_trans')
        if J is not None:
            return J
    
    # 8. Terpènes
    if 'terpene' in mol_types:
        J = get_terpene_coupling(mol, h1_idx, h2_idx, 'menthane')
        if J is not None:
            return J
    
    # 9. Macrocycles
    if 'macrocycle' in mol_types:
        # Trouver la taille du plus grand cycle
        ring_info = mol.GetRingInfo()
        max_ring = max(len(r) for r in ring_info.AtomRings()) if ring_info.AtomRings() else 0
        
        J = get_macrocycle_coupling(mol, h1_idx, h2_idx, max_ring)
        if J is not None:
            return J
    
    # 10. Nucléosides
    if 'nucleoside' in mol_types:
        J = get_nucleotide_coupling(mol, h1_idx, h2_idx, 'ribose')
        if J is not None:
            return J
    
    # Pas de table spécialisée applicable
    return None


# =============================================================================
# SECTION 13: TESTS
# =============================================================================

def _test_coupling_tables():
    """Tests des tables de couplages."""
    print("=" * 60)
    print("TEST DES TABLES DE COUPLAGES COMPLÈTES")
    print("=" * 60)
    
    test_cases = [
        ("Propyne (alcyne)", "CC#C"),
        ("Allène", "C=C=C"),
        ("Benzaldéhyde", "O=Cc1ccccc1"),
        ("Propanal", "CCC=O"),
        ("N-méthylacétamide", "CC(=O)NC"),
        ("Alanine", "C[C@H](N)C(=O)O"),
        ("Triméthylphosphine", "CP(C)C"),
        ("Triméthylsilane", "C[SiH](C)C"),
        ("Pinène", "CC1=CCC2CC1C2(C)C"),
        ("Macrolide 12", "C1CCCCCCCCCCC(=O)O1"),  # 12-membered lactone
        ("Adénosine", "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"),
    ]
    
    for name, smiles in test_cases:
        print(f"\n{name}")
        print("-" * 40)
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print("  ERREUR: SMILES invalide")
            continue
        
        # Classifier
        types = MoleculeClassifier.classify(mol)
        print(f"  Types détectés: {types}")
        
        # Tester le dispatch
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        
        # Trouver 2 protons
        h_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'H']
        
        if len(h_indices) >= 2:
            J = get_specialized_coupling(mol, h_indices[0], h_indices[1], 3)
            print(f"  J test: {J} Hz")
    
    print("\n" + "=" * 60)
    print("TESTS TERMINÉS")
    print("=" * 60)


if __name__ == "__main__":
    _test_coupling_tables()