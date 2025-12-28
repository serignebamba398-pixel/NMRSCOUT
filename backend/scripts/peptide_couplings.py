# =============================================================================
# peptide_couplings.py
# Module robuste pour les couplages dans les peptides et protéines
# =============================================================================

"""
Module de couplages RMN pour peptides et protéines.

Fonctionnalités:
- Équations de Karplus multiples (Bystrov, Pardi, Vuister, Habeck)
- Détection automatique backbone vs sidechain
- Tous les acides aminés standards
- Structures secondaires (α-helix, β-sheet, turns)
- Couplages ³J(HN-Hα), ³J(Hα-Hβ), ³J(N-Cγ), etc.
- Proline (cis/trans), Glycine (²J geminal)

Références:
- Bystrov VF (1976) Prog NMR Spectrosc 10:41-82
- Pardi A et al. (1984) J Mol Biol 180:741-751
- Vuister GW & Bax A (1993) JACS 115:7772-7777
- Habeck M et al. (2005) J Magn Reson 177:160-165

Auteur: [Votre nom]
Date: 2024
"""

from typing import Dict, List, Tuple, Optional, Any, Union
from dataclasses import dataclass, field
from enum import Enum
import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
except ImportError:
    raise ImportError("RDKit requis")


# =============================================================================
# SECTION 1: Énumérations et constantes
# =============================================================================

class AminoAcid(Enum):
    """Les 20 acides aminés standards."""
    ALA = "ALA"  # Alanine
    ARG = "ARG"  # Arginine
    ASN = "ASN"  # Asparagine
    ASP = "ASP"  # Aspartate
    CYS = "CYS"  # Cystéine
    GLN = "GLN"  # Glutamine
    GLU = "GLU"  # Glutamate
    GLY = "GLY"  # Glycine
    HIS = "HIS"  # Histidine
    ILE = "ILE"  # Isoleucine
    LEU = "LEU"  # Leucine
    LYS = "LYS"  # Lysine
    MET = "MET"  # Méthionine
    PHE = "PHE"  # Phénylalanine
    PRO = "PRO"  # Proline
    SER = "SER"  # Sérine
    THR = "THR"  # Thréonine
    TRP = "TRP"  # Tryptophane
    TYR = "TYR"  # Tyrosine
    VAL = "VAL"  # Valine


class SecondaryStructure(Enum):
    """Types de structures secondaires."""
    ALPHA_HELIX = "alpha_helix"
    BETA_SHEET_PARALLEL = "beta_parallel"
    BETA_SHEET_ANTIPARALLEL = "beta_antiparallel"
    TURN_TYPE_I = "turn_I"
    TURN_TYPE_II = "turn_II"
    TURN_TYPE_I_PRIME = "turn_I'"
    TURN_TYPE_II_PRIME = "turn_II'"
    THREE_TEN_HELIX = "310_helix"
    PI_HELIX = "pi_helix"
    POLYPROLINE_II = "PPII"
    RANDOM_COIL = "random_coil"


class CouplingType(Enum):
    """Types de couplages dans les peptides."""
    HN_HA = "3J_HN_HA"           # ³J(HN-Hα)
    HA_HB = "3J_HA_HB"           # ³J(Hα-Hβ)
    HB_HG = "3J_HB_HG"           # ³J(Hβ-Hγ)
    HG_HD = "3J_HG_HD"           # ³J(Hγ-Hδ)
    N_CG = "3J_N_CG"             # ³J(N-Cγ) pour Asn/Asp
    CO_CG = "3J_CO_CG"           # ³J(CO-Cγ)
    HA_HA_GLY = "2J_HA_HA"       # ²J geminal Gly
    HB_HB = "2J_HB_HB"           # ²J geminal Hβ
    NH2_SIDECHAIN = "3J_NH2"    # Couplages NH₂ (Asn, Gln)


# =============================================================================
# SECTION 2: Paramètres de Karplus
# =============================================================================

@dataclass
class KarplusParameters:
    """Paramètres pour l'équation de Karplus généralisée.
    
    J(θ) = A·cos²(θ + phase) + B·cos(θ + phase) + C
    
    Ou forme étendue (6 paramètres):
    J(θ) = C0 + C1·cos(θ) + C2·cos(2θ) + C3·sin(θ) + C4·sin(2θ)
    """
    A: float = 0.0
    B: float = 0.0
    C: float = 0.0
    phase: float = 0.0  # En degrés
    
    # Forme étendue (optionnelle)
    C0: Optional[float] = None
    C1: Optional[float] = None
    C2: Optional[float] = None
    C3: Optional[float] = None
    C4: Optional[float] = None
    
    # Métadonnées
    name: str = ""
    reference: str = ""
    
    def calculate(self, theta_deg: float) -> float:
        """Calcule J pour un angle θ donné en degrés."""
        theta = np.radians(theta_deg + self.phase)
        
        if self.C0 is not None:
            # Forme étendue à 5 paramètres
            J = (self.C0 + 
                 self.C1 * np.cos(theta) + 
                 self.C2 * np.cos(2 * theta))
            
            if self.C3 is not None:
                J += self.C3 * np.sin(theta)
            if self.C4 is not None:
                J += self.C4 * np.sin(2 * theta)
            
            return max(0, J)
        else:
            # Forme standard à 3 paramètres
            J = self.A * np.cos(theta)**2 + self.B * np.cos(theta) + self.C
            return max(0, J)


# Dictionnaire des paramètres de Karplus par type de couplage
KARPLUS_PARAMETERS: Dict[str, Dict[str, KarplusParameters]] = {
    
    # ═══════════════════════════════════════════════════════════════════════
    # ³J(HN-Hα) - Couplage backbone principal
    # ═══════════════════════════════════════════════════════════════════════
    '3J_HN_HA': {
        # Bystrov (1976) - Classique
        'bystrov': KarplusParameters(
            A=6.4, B=-1.4, C=1.9, phase=-60,
            name="Bystrov 1976",
            reference="Prog NMR Spectrosc 10:41-82"
        ),
        
        # Pardi (1984) - Affiné pour protéines
        'pardi': KarplusParameters(
            A=6.98, B=-1.38, C=1.72, phase=-60,
            name="Pardi 1984",
            reference="J Mol Biol 180:741-751"
        ),
        
        # Vuister & Bax (1993) - Haute précision
        'vuister': KarplusParameters(
            A=6.51, B=-1.76, C=1.60, phase=-60,
            name="Vuister-Bax 1993",
            reference="JACS 115:7772-7777"
        ),
        
        # Wang & Bhm (1996)
        'wang': KarplusParameters(
            A=7.09, B=-1.42, C=1.55, phase=-60,
            name="Wang-Bax 1996",
            reference="JACS 118:2483-2494"
        ),
        
        # Habeck (2005) - Bayésien
        'habeck': KarplusParameters(
            A=6.7, B=-1.3, C=1.5, phase=-60,
            name="Habeck 2005",
            reference="J Magn Reson 177:160-165"
        ),
        
        # Schmidt (1999) - Pour protéines non-structurées
        'schmidt': KarplusParameters(
            A=6.0, B=-1.0, C=2.0, phase=-60,
            name="Schmidt 1999",
            reference="J Biomol NMR 14:1-12"
        ),
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # ³J(Hα-Hβ) - Couplage Cα-Cβ (angle χ1)
    # ═══════════════════════════════════════════════════════════════════════
    '3J_HA_HB': {
        # DeMarco (1978)
        'demarco': KarplusParameters(
            A=9.5, B=-1.6, C=1.8, phase=0,
            name="DeMarco 1978",
            reference="Biopolymers 17:617"
        ),
        
        # Pachler (1964) - Standard
        'pachler': KarplusParameters(
            A=9.4, B=-1.1, C=0.4, phase=0,
            name="Pachler 1964",
            reference="Spectrochim Acta 20:581"
        ),
        
        # Haasnoot (1980) - Avec électronégativité
        'haasnoot': KarplusParameters(
            A=13.22, B=-0.99, C=0.0, phase=0,
            name="Haasnoot 1980",
            reference="Tetrahedron 36:2783"
        ),
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # ³J(N-Cγ) - Pour Asn, Asp, His
    # ═══════════════════════════════════════════════════════════════════════
    '3J_N_CG': {
        'perez': KarplusParameters(
            A=4.5, B=-1.2, C=0.5, phase=0,
            name="Perez 2001",
            reference="Magn Reson Chem 39:S24"
        ),
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # ³J(CO-Cγ) - Pour déterminer χ1
    # ═══════════════════════════════════════════════════════════════════════
    '3J_CO_CG': {
        'hu': KarplusParameters(
            A=4.4, B=-1.2, C=0.1, phase=0,
            name="Hu-Bax 1997",
            reference="JACS 119:6360"
        ),
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # ³J(Hβ-Hγ) - Chaîne latérale
    # ═══════════════════════════════════════════════════════════════════════
    '3J_HB_HG': {
        'standard': KarplusParameters(
            A=9.4, B=-1.1, C=0.4, phase=0,
            name="Standard aliphatique",
            reference="Pachler 1964"
        ),
    },
}


# =============================================================================
# SECTION 3: Angles dièdres par structure secondaire
# =============================================================================

# Angles φ (phi) et ψ (psi) typiques en degrés
SECONDARY_STRUCTURE_ANGLES: Dict[SecondaryStructure, Dict[str, float]] = {
    SecondaryStructure.ALPHA_HELIX: {
        'phi': -57.0,
        'psi': -47.0,
        'omega': 180.0,
    },
    SecondaryStructure.BETA_SHEET_PARALLEL: {
        'phi': -119.0,
        'psi': 113.0,
        'omega': 180.0,
    },
    SecondaryStructure.BETA_SHEET_ANTIPARALLEL: {
        'phi': -139.0,
        'psi': 135.0,
        'omega': 180.0,
    },
    SecondaryStructure.THREE_TEN_HELIX: {
        'phi': -49.0,
        'psi': -26.0,
        'omega': 180.0,
    },
    SecondaryStructure.PI_HELIX: {
        'phi': -57.0,
        'psi': -70.0,
        'omega': 180.0,
    },
    SecondaryStructure.TURN_TYPE_I: {
        'phi': -60.0,  # Position i+1
        'psi': -30.0,
        'omega': 180.0,
    },
    SecondaryStructure.TURN_TYPE_II: {
        'phi': -60.0,
        'psi': 120.0,
        'omega': 180.0,
    },
    SecondaryStructure.POLYPROLINE_II: {
        'phi': -75.0,
        'psi': 145.0,
        'omega': 180.0,
    },
    SecondaryStructure.RANDOM_COIL: {
        'phi': -120.0,  # Moyenne approximative
        'psi': 120.0,
        'omega': 180.0,
    },
}

# ³J(HN-Hα) attendu par structure secondaire
EXPECTED_3J_HN_HA: Dict[SecondaryStructure, float] = {
    SecondaryStructure.ALPHA_HELIX: 3.9,
    SecondaryStructure.BETA_SHEET_PARALLEL: 8.9,
    SecondaryStructure.BETA_SHEET_ANTIPARALLEL: 8.5,
    SecondaryStructure.THREE_TEN_HELIX: 4.2,
    SecondaryStructure.PI_HELIX: 4.5,
    SecondaryStructure.TURN_TYPE_I: 5.5,
    SecondaryStructure.TURN_TYPE_II: 6.5,
    SecondaryStructure.POLYPROLINE_II: 6.0,
    SecondaryStructure.RANDOM_COIL: 7.0,
}


# =============================================================================
# SECTION 4: Couplages spécifiques par acide aminé
# =============================================================================

AMINO_ACID_SPECIFIC_COUPLINGS: Dict[AminoAcid, Dict[str, Any]] = {
    
    # ═══════════════════════════════════════════════════════════════════════
    # GLYCINE - Cas spécial avec ²J geminal
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.GLY: {
        '2J_HA_HA': -17.5,  # Couplage geminal fort
        'notes': "Pas de Hβ, deux Hα équivalents",
        'equivalent_Ha': True,
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # PROLINE - Pas de NH, cycle pyrrolidine
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.PRO: {
        'has_NH': False,
        'ring_couplings': {
            '3J_HA_HB2': {'endo': 6.0, 'exo': 8.5},
            '3J_HA_HB3': {'endo': 8.5, 'exo': 6.0},
            '3J_HB_HG': 7.0,
            '3J_HG_HD': 7.0,
            '2J_HB2_HB3': -12.5,
            '2J_HG2_HG3': -12.5,
            '2J_HD2_HD3': -12.5,
        },
        'cis_trans': {
            'omega_cis': 0.0,      # Rare (~5%)
            'omega_trans': 180.0,  # Majoritaire (~95%)
            '3J_HA_HD_cis': 2.0,
            '3J_HA_HD_trans': 0.5,
        },
        'notes': "Cycle 5 membres, pas de HN amide",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # ALANINE - Simple, référence
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.ALA: {
        '3J_HA_HB': 7.0,  # Rotation libre du CH₃
        'notes': "Méthyle, rotation libre",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # VALINE - Branché au Cβ
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.VAL: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB': 11.0, 'chi1': -60},
            'gauche_plus': {'3J_HA_HB': 3.5, 'chi1': 60},
            'trans': {'3J_HA_HB': 3.5, 'chi1': 180},
        },
        '3J_HB_HG': 6.5,  # Méthyles
        'notes': "Isopropyle, 3 rotamères χ1",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # LEUCINE
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.LEU: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},
            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},
            'trans': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 3.5},
        },
        '2J_HB2_HB3': -14.0,
        '3J_HB_HG': 6.5,
        '3J_HG_HD': 6.5,
        'notes': "Isobutyle",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # ISOLEUCINE
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.ILE: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB': 11.0},
            'gauche_plus': {'3J_HA_HB': 3.5},
            'trans': {'3J_HA_HB': 3.5},
        },
        '3J_HB_HG12': 7.0,
        '3J_HB_HG13': 7.0,
        '2J_HG12_HG13': -14.0,
        '3J_HG_HD': 6.5,
        'notes': "sec-butyle",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # THRÉONINE
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.THR: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB': 11.0},
            'gauche_plus': {'3J_HA_HB': 3.5},
            'trans': {'3J_HA_HB': 3.5},
        },
        '3J_HB_HG': 6.5,  # Méthyle
        'notes': "Hydroxyle sur Cβ",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # SÉRINE
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.SER: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},
            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},
            'trans': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 3.5},
        },
        '2J_HB2_HB3': -12.0,
        'notes': "Hydroxyle sur Cβ, CH₂",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # CYSTÉINE
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.CYS: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},
            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},
            'trans': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 3.5},
        },
        '2J_HB2_HB3': -14.5,  # Plus grand à cause du S
        'notes': "Thiol, peut former ponts disulfure",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # MÉTHIONINE
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.MET: {
        '3J_HA_HB': 7.0,  # Moyenne
        '3J_HB_HG': 7.0,
        '2J_HB2_HB3': -14.5,
        '2J_HG2_HG3': -14.5,
        '3J_HG_CH3': 0.0,  # S-CH₃ pas de couplage visible
        'notes': "Thioéther, S-CH₃",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # PHÉNYLALANINE
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.PHE: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},
            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},
            'trans': {'3J_HA_HB2': 5.5, '3J_HA_HB3': 5.5},
        },
        '2J_HB2_HB3': -14.0,
        'aromatic': {
            '3J_ortho': 7.5,
            '4J_meta': 1.5,
            '5J_para': 0.5,
        },
        'notes': "Benzyle",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # TYROSINE
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.TYR: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},
            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},
            'trans': {'3J_HA_HB2': 5.5, '3J_HA_HB3': 5.5},
        },
        '2J_HB2_HB3': -14.0,
        'aromatic': {
            '3J_ortho': 8.5,  # Système AA'BB'
            '4J_meta': 2.5,
        },
        'notes': "p-hydroxybenzyle, AA'BB'",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # TRYPTOPHANE
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.TRP: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},
            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},
            'trans': {'3J_HA_HB2': 5.5, '3J_HA_HB3': 5.5},
        },
        '2J_HB2_HB3': -14.5,
        'indole': {
            '3J_H2_H7a': 0.0,  # Pas de couplage direct
            '3J_H4_H5': 8.0,
            '3J_H5_H6': 7.5,
            '3J_H6_H7': 8.0,
            '4J_H4_H6': 1.0,
            '4J_H5_H7': 1.2,
        },
        'notes': "Indole fusionné",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # HISTIDINE
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.HIS: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},
            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},
            'trans': {'3J_HA_HB2': 5.5, '3J_HA_HB3': 5.5},
        },
        '2J_HB2_HB3': -14.5,
        'imidazole': {
            '4J_H2_H4': 1.0,  # Couplage à travers le cycle
            '4J_H2_H5': 1.5,
            '4J_H4_H5': 1.0,
        },
        'tautomers': ['Nε-H', 'Nδ-H'],
        'notes': "Imidazole, pKa ~6",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # ASPARAGINE
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.ASN: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},
            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},
            'trans': {'3J_HA_HB2': 5.5, '3J_HA_HB3': 5.5},
        },
        '2J_HB2_HB3': -15.5,  # α-carbonyle
        'amide_sidechain': {
            '2J_NH2': -14.0,  # Geminal NH₂ (inéquivalents)
            '3J_HB_NH': 6.0,
        },
        'notes': "Amide primaire en chaîne latérale",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # GLUTAMINE
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.GLN: {
        '3J_HA_HB': 7.0,
        '3J_HB_HG': 7.0,
        '2J_HB2_HB3': -14.0,
        '2J_HG2_HG3': -15.5,  # α-carbonyle
        'amide_sidechain': {
            '2J_NH2': -14.0,
            '3J_HG_NH': 6.0,
        },
        'notes': "Amide primaire, chaîne plus longue",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # ASPARTATE
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.ASP: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},
            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},
            'trans': {'3J_HA_HB2': 5.5, '3J_HA_HB3': 5.5},
        },
        '2J_HB2_HB3': -16.5,  # α-carboxyle (effet électroattracteur)
        'notes': "Carboxyle, chargé négativement à pH >4",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # GLUTAMATE
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.GLU: {
        '3J_HA_HB': 7.0,
        '3J_HB_HG': 7.0,
        '2J_HB2_HB3': -14.0,
        '2J_HG2_HG3': -16.5,  # α-carboxyle
        'notes': "Carboxyle, chaîne plus longue",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # LYSINE
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.LYS: {
        '3J_HA_HB': 7.0,
        '3J_HB_HG': 7.0,
        '3J_HG_HD': 7.0,
        '3J_HD_HE': 7.0,
        '2J_CH2': -12.5,  # Tous les CH₂
        '3J_HE_NH3': 6.0,  # NH₃⁺ en échange
        'notes': "Longue chaîne, NH₃⁺ terminal",
    },
    
    # ═══════════════════════════════════════════════════════════════════════
    # ARGININE
    # ═══════════════════════════════════════════════════════════════════════
    AminoAcid.ARG: {
        '3J_HA_HB': 7.0,
        '3J_HB_HG': 7.0,
        '3J_HG_HD': 7.0,
        '2J_CH2': -12.5,
        'guanidinium': {
            '3J_HD_NE': 6.5,
            'NH_equivalent': True,  # En échange rapide
        },
        'notes': "Guanidinium, toujours chargé +",
    },
}


# =============================================================================
# SECTION 5: Fonctions de calcul principales
# =============================================================================

def calculate_3J_HN_HA(
    phi: Optional[float] = None,
    secondary_structure: Optional[SecondaryStructure] = None,
    karplus_set: str = 'vuister',
    amino_acid: Optional[AminoAcid] = None,
) -> float:
    """
    Calcule ³J(HN-Hα) à partir de l'angle φ ou de la structure secondaire.
    
    Args:
        phi: Angle φ en degrés (prioritaire si fourni)
        secondary_structure: Type de structure secondaire
        karplus_set: Jeu de paramètres ('bystrov', 'pardi', 'vuister', 'habeck')
        amino_acid: Acide aminé (pour cas spéciaux)
    
    Returns:
        ³J(HN-Hα) en Hz
    
    Example:
        >>> calculate_3J_HN_HA(phi=-57)  # Hélice α
        3.9
        >>> calculate_3J_HN_HA(secondary_structure=SecondaryStructure.ALPHA_HELIX)
        3.9
    """
    # Cas spéciaux
    if amino_acid == AminoAcid.PRO:
        return 0.0  # Proline n'a pas de HN
    
    if amino_acid == AminoAcid.GLY:
        # Glycine a des valeurs légèrement différentes
        if phi is not None:
            params = KARPLUS_PARAMETERS['3J_HN_HA'].get(karplus_set)
            if params:
                return params.calculate(phi)
        elif secondary_structure:
            return EXPECTED_3J_HN_HA.get(secondary_structure, 7.0)
        return 6.5  # Valeur moyenne pour Gly
    
    # Calcul avec angle φ
    if phi is not None:
        params = KARPLUS_PARAMETERS['3J_HN_HA'].get(karplus_set)
        if params:
            return params.calculate(phi)
        # Fallback
        theta = np.radians(phi - 60)
        return max(0, 6.4 * np.cos(theta)**2 - 1.4 * np.cos(theta) + 1.9)
    
    # Utiliser la structure secondaire
    if secondary_structure:
        return EXPECTED_3J_HN_HA.get(secondary_structure, 7.0)
    
    # Valeur par défaut (random coil)
    return 7.0


def calculate_3J_HA_HB(
    chi1: Optional[float] = None,
    amino_acid: Optional[AminoAcid] = None,
    hb_type: str = 'HB2',  # 'HB2', 'HB3', ou 'HB' pour méthyle
    karplus_set: str = 'demarco',
) -> float:
    """
    Calcule ³J(Hα-Hβ) à partir de l'angle χ1.
    
    Args:
        chi1: Angle χ1 en degrés
        amino_acid: Acide aminé
        hb_type: Type de Hβ ('HB2', 'HB3', 'HB' pour CH unique)
        karplus_set: Jeu de paramètres
    
    Returns:
        ³J(Hα-Hβ) en Hz
    """
    # Cas spéciaux
    if amino_acid == AminoAcid.GLY:
        return 0.0  # Pas de Hβ
    
    if amino_acid == AminoAcid.ALA:
        return 7.0  # Méthyle, rotation libre
    
    # Utiliser les données spécifiques à l'acide aminé
    if amino_acid and amino_acid in AMINO_ACID_SPECIFIC_COUPLINGS:
        aa_data = AMINO_ACID_SPECIFIC_COUPLINGS[amino_acid]
        
        # Vérifier si χ1 donné et rotamères disponibles
        if chi1 is not None and 'chi1_rotamers' in aa_data:
            # Déterminer le rotamère
            if -90 <= chi1 < 30:
                rotamer = 'gauche_minus'
            elif 30 <= chi1 < 150:
                rotamer = 'trans'
            else:
                rotamer = 'gauche_plus'
            
            rot_data = aa_data['chi1_rotamers'].get(rotamer, {})
            key = f'3J_HA_{hb_type}'
            if key in rot_data:
                return rot_data[key]
    
    # Calcul avec Karplus si χ1 fourni
    if chi1 is not None:
        params = KARPLUS_PARAMETERS['3J_HA_HB'].get(karplus_set)
        if params:
            # Ajuster la phase selon Hβ2 ou Hβ3
            phase_adjustment = 0 if hb_type == 'HB2' else 120
            return params.calculate(chi1 + phase_adjustment)
    
    # Valeur par défaut (rotation libre)
    return 7.0


def calculate_2J_geminal(
    amino_acid: Optional[AminoAcid] = None,
    position: str = 'HB',  # 'HA' pour Gly, 'HB', 'HG', 'HD'
) -> float:
    """
    Calcule le couplage ²J geminal.
    
    Args:
        amino_acid: Acide aminé
        position: Position des protons
    
    Returns:
        ²J en Hz (négatif)
    """
    # Glycine Hα-Hα
    if amino_acid == AminoAcid.GLY and position == 'HA':
        return AMINO_ACID_SPECIFIC_COUPLINGS[AminoAcid.GLY]['2J_HA_HA']
    
    # Couplages spécifiques par acide aminé
    if amino_acid and amino_acid in AMINO_ACID_SPECIFIC_COUPLINGS:
        aa_data = AMINO_ACID_SPECIFIC_COUPLINGS[amino_acid]
        key = f'2J_{position}2_{position}3'
        if key in aa_data:
            return aa_data[key]
    
    # Valeurs par défaut selon la position
    defaults = {
        'HA': -17.5,  # Gly
        'HB': -14.0,
        'HG': -14.0,
        'HD': -12.5,
        'HE': -12.5,
    }
    
    return defaults.get(position, -12.5)


def get_peptide_coupling_robust(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    distance: Optional[int] = None,
    phi: Optional[float] = None,
    psi: Optional[float] = None,
    chi1: Optional[float] = None,
    amino_acid: Optional[AminoAcid] = None,
    secondary_structure: Optional[SecondaryStructure] = None,
) -> Optional[float]:
    """
    Fonction principale robuste pour calculer les couplages peptidiques.
    
    Args:
        mol: Molécule RDKit
        h1_idx, h2_idx: Indices des protons
        distance: Distance topologique (auto-détectée si None)
        phi, psi, chi1: Angles dièdres (extraits de la géométrie si None)
        amino_acid: Type d'acide aminé (auto-détecté si None)
        secondary_structure: Structure secondaire
    
    Returns:
        J en Hz, ou None si non applicable
    """
    # Obtenir les atomes parents
    h1 = mol.GetAtomWithIdx(h1_idx)
    h2 = mol.GetAtomWithIdx(h2_idx)
    
    heavy1_list = [n for n in h1.GetNeighbors() if n.GetSymbol() != 'H']
    heavy2_list = [n for n in h2.GetNeighbors() if n.GetSymbol() != 'H']
    
    if not heavy1_list or not heavy2_list:
        return None
    
    heavy1 = heavy1_list[0]
    heavy2 = heavy2_list[0]
    
    # Auto-détection de la distance
    if distance is None:
        try:
            path = Chem.GetShortestPath(mol, h1_idx, h2_idx)
            distance = len(path) - 1
        except:
            return None
    
    # Détecter le type de couplage
    coupling_type = _detect_peptide_coupling_type(mol, h1_idx, h2_idx, heavy1, heavy2)
    
    if coupling_type is None:
        return None
    
    # Extraire les angles dièdres si disponibles
    if mol.GetNumConformers() > 0 and phi is None:
        try:
            phi = _extract_phi_angle(mol, heavy1.GetIdx())
        except:
            pass
    
    if mol.GetNumConformers() > 0 and chi1 is None:
        try:
            chi1 = _extract_chi1_angle(mol, heavy1.GetIdx())
        except:
            pass
    
    # Calculer selon le type
    if coupling_type == CouplingType.HN_HA:
        return calculate_3J_HN_HA(
            phi=phi,
            secondary_structure=secondary_structure,
            amino_acid=amino_acid,
        )
    
    elif coupling_type == CouplingType.HA_HB:
        return calculate_3J_HA_HB(
            chi1=chi1,
            amino_acid=amino_acid,
        )
    
    elif coupling_type == CouplingType.HA_HA_GLY:
        return calculate_2J_geminal(AminoAcid.GLY, 'HA')
    
    elif coupling_type == CouplingType.HB_HB:
        return calculate_2J_geminal(amino_acid, 'HB')
    
    elif coupling_type == CouplingType.HB_HG:
        return 7.0  # Standard aliphatique
    
    elif coupling_type == CouplingType.HG_HD:
        return 7.0
    
    return None


def _detect_peptide_coupling_type(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    heavy1: Chem.Atom,
    heavy2: Chem.Atom
) -> Optional[CouplingType]:
    """Détecte automatiquement le type de couplage peptidique."""
    
    # Vérifier si geminal (même atome lourd)
    if heavy1.GetIdx() == heavy2.GetIdx():
        carbon = heavy1
        
        # Glycine Hα-Hα ?
        if _is_alpha_carbon(mol, carbon.GetIdx()):
            return CouplingType.HA_HA_GLY
        
        # Autre CH₂
        return CouplingType.HB_HB
    
    # Couplage vicinal
    # Vérifier HN-Hα
    is_h1_on_N = heavy1.GetSymbol() == 'N'
    is_h2_on_N = heavy2.GetSymbol() == 'N'
    is_h1_on_Ca = _is_alpha_carbon(mol, heavy1.GetIdx())
    is_h2_on_Ca = _is_alpha_carbon(mol, heavy2.GetIdx())
    
    if (is_h1_on_N and is_h2_on_Ca) or (is_h2_on_N and is_h1_on_Ca):
        return CouplingType.HN_HA
    
    # Hα-Hβ
    is_h1_on_Cb = _is_beta_carbon(mol, heavy1.GetIdx())
    is_h2_on_Cb = _is_beta_carbon(mol, heavy2.GetIdx())
    
    if (is_h1_on_Ca and is_h2_on_Cb) or (is_h2_on_Ca and is_h1_on_Cb):
        return CouplingType.HA_HB
    
    # Hβ-Hγ
    is_h1_on_Cg = _is_gamma_carbon(mol, heavy1.GetIdx())
    is_h2_on_Cg = _is_gamma_carbon(mol, heavy2.GetIdx())
    
    if (is_h1_on_Cb and is_h2_on_Cg) or (is_h2_on_Cb and is_h1_on_Cg):
        return CouplingType.HB_HG
    
    return None


def _is_alpha_carbon(mol: Chem.Mol, c_idx: int) -> bool:
    """Vérifie si un carbone est un Cα (attaché à N et C=O)."""
    atom = mol.GetAtomWithIdx(c_idx)
    
    if atom.GetSymbol() != 'C':
        return False
    
    has_N = False
    has_CO = False
    
    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() == 'N':
            has_N = True
        elif neighbor.GetSymbol() == 'C':
            for nn in neighbor.GetNeighbors():
                if nn.GetSymbol() == 'O':
                    bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nn.GetIdx())
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                        has_CO = True
    
    return has_N and has_CO


def _is_beta_carbon(mol: Chem.Mol, c_idx: int) -> bool:
    """Vérifie si un carbone est un Cβ (attaché à Cα)."""
    atom = mol.GetAtomWithIdx(c_idx)
    
    if atom.GetSymbol() != 'C':
        return False
    
    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() == 'C' and _is_alpha_carbon(mol, neighbor.GetIdx()):
            return True
    
    return False


def _is_gamma_carbon(mol: Chem.Mol, c_idx: int) -> bool:
    """Vérifie si un carbone est un Cγ (attaché à Cβ)."""
    atom = mol.GetAtomWithIdx(c_idx)
    
    if atom.GetSymbol() != 'C':
        return False
    
    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() == 'C' and _is_beta_carbon(mol, neighbor.GetIdx()):
            return True
    
    return False


def _extract_phi_angle(mol: Chem.Mol, ca_idx: int) -> Optional[float]:
    """Extrait l'angle φ (C-N-Cα-C) de la géométrie 3D."""
    if mol.GetNumConformers() == 0:
        return None
    
    ca = mol.GetAtomWithIdx(ca_idx)
    
    # Trouver N, C(O) du backbone
    n_idx = None
    co_idx = None
    prev_c_idx = None
    
    for neighbor in ca.GetNeighbors():
        if neighbor.GetSymbol() == 'N':
            n_idx = neighbor.GetIdx()
            # Trouver le C précédent
            for nn in neighbor.GetNeighbors():
                if nn.GetSymbol() == 'C' and nn.GetIdx() != ca_idx:
                    for nnn in nn.GetNeighbors():
                        if nnn.GetSymbol() == 'O':
                            bond = mol.GetBondBetweenAtoms(nn.GetIdx(), nnn.GetIdx())
                            if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                                prev_c_idx = nn.GetIdx()
        elif neighbor.GetSymbol() == 'C':
            for nn in neighbor.GetNeighbors():
                if nn.GetSymbol() == 'O':
                    bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nn.GetIdx())
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                        co_idx = neighbor.GetIdx()
    
    if all(x is not None for x in [prev_c_idx, n_idx, co_idx]):
        try:
            conf = mol.GetConformer()
            phi = AllChem.GetDihedralDeg(conf, prev_c_idx, n_idx, ca_idx, co_idx)
            return phi
        except:
            pass
    
    return None


def _extract_chi1_angle(mol: Chem.Mol, ca_idx: int) -> Optional[float]:
    """Extrait l'angle χ1 (N-Cα-Cβ-Xγ) de la géométrie 3D."""
    if mol.GetNumConformers() == 0:
        return None
    
    ca = mol.GetAtomWithIdx(ca_idx)
    
    n_idx = None
    cb_idx = None
    xg_idx = None
    
    for neighbor in ca.GetNeighbors():
        if neighbor.GetSymbol() == 'N':
            n_idx = neighbor.GetIdx()
        elif neighbor.GetSymbol() == 'C':
            # Vérifier si c'est Cβ (pas le carbonyle)
            is_carbonyl = False
            for nn in neighbor.GetNeighbors():
                if nn.GetSymbol() == 'O':
                    bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nn.GetIdx())
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                        is_carbonyl = True
            
            if not is_carbonyl:
                cb_idx = neighbor.GetIdx()
                # Trouver Xγ (premier atome lourd après Cβ)
                for nn in neighbor.GetNeighbors():
                    if nn.GetIdx() != ca_idx and nn.GetSymbol() != 'H':
                        xg_idx = nn.GetIdx()
                        break
    
    if all(x is not None for x in [n_idx, cb_idx, xg_idx]):
        try:
            conf = mol.GetConformer()
            chi1 = AllChem.GetDihedralDeg(conf, n_idx, ca_idx, cb_idx, xg_idx)
            return chi1
        except:
            pass
    
    return None


# =============================================================================
# SECTION 6: Utilitaires et export
# =============================================================================

def predict_secondary_structure_from_3J(j_hn_ha: float) -> SecondaryStructure:
    """
    Prédit la structure secondaire à partir de ³J(HN-Hα).
    
    Args:
        j_hn_ha: Valeur de ³J(HN-Hα) en Hz
    
    Returns:
        Structure secondaire la plus probable
    """
    if j_hn_ha < 5.0:
        return SecondaryStructure.ALPHA_HELIX
    elif j_hn_ha < 6.0:
        return SecondaryStructure.THREE_TEN_HELIX
    elif j_hn_ha < 7.5:
        return SecondaryStructure.RANDOM_COIL
    else:
        return SecondaryStructure.BETA_SHEET_ANTIPARALLEL


def get_all_backbone_couplings(
    phi: float,
    psi: float,
    amino_acid: AminoAcid = AminoAcid.ALA,
) -> Dict[str, float]:
    """
    Calcule tous les couplages backbone pour un résidu.
    
    Args:
        phi: Angle φ en degrés
        psi: Angle ψ en degrés
        amino_acid: Type d'acide aminé
    
    Returns:
        Dict de tous les couplages
    """
    result = {}
    
    # ³J(HN-Hα)
    if amino_acid != AminoAcid.PRO:
        result['3J_HN_HA'] = calculate_3J_HN_HA(phi=phi, amino_acid=amino_acid)
    
    # ²J(Hα-Hα) pour Glycine
    if amino_acid == AminoAcid.GLY:
        result['2J_HA_HA'] = calculate_2J_geminal(AminoAcid.GLY, 'HA')
    
    return result


def get_sidechain_couplings(
    amino_acid: AminoAcid,
    chi1: Optional[float] = None,
    chi2: Optional[float] = None,
) -> Dict[str, float]:
    """
    Calcule les couplages de la chaîne latérale.
    
    Args:
        amino_acid: Type d'acide aminé
        chi1, chi2: Angles χ1 et χ2
    
    Returns:
        Dict de tous les couplages sidechain
    """
    result = {}
    
    if amino_acid not in AMINO_ACID_SPECIFIC_COUPLINGS:
        return result
    
    aa_data = AMINO_ACID_SPECIFIC_COUPLINGS[amino_acid]
    
    # ³J(Hα-Hβ)
    if 'chi1_rotamers' in aa_data and chi1 is not None:
        result['3J_HA_HB2'] = calculate_3J_HA_HB(chi1, amino_acid, 'HB2')
        result['3J_HA_HB3'] = calculate_3J_HA_HB(chi1, amino_acid, 'HB3')
    elif '3J_HA_HB' in aa_data:
        result['3J_HA_HB'] = aa_data['3J_HA_HB']
    
    # ²J geminaux
    for key, value in aa_data.items():
        if key.startswith('2J_'):
            result[key] = value
    
    # Couplages aromatiques
    if 'aromatic' in aa_data:
        for key, value in aa_data['aromatic'].items():
            result[f'aromatic_{key}'] = value
    
    # Couplages hétérocycliques
    for ring_type in ['imidazole', 'indole']:
        if ring_type in aa_data:
            for key, value in aa_data[ring_type].items():
                result[f'{ring_type}_{key}'] = value
    
    return result


# =============================================================================
# SECTION 7: Tests
# =============================================================================

def _test_peptide_module():
    """Tests du module peptide."""
    print("=" * 70)
    print("TEST MODULE PEPTIDES ROBUSTE")
    print("=" * 70)
    
    # Test 1: Karplus pour différentes structures secondaires
    print("\n1. ³J(HN-Hα) par structure secondaire:")
    print("-" * 50)
    
    for ss in SecondaryStructure:
        angles = SECONDARY_STRUCTURE_ANGLES.get(ss, {})
        phi = angles.get('phi')
        
        if phi:
            j_calc = calculate_3J_HN_HA(phi=phi)
            j_expected = EXPECTED_3J_HN_HA.get(ss, 7.0)
            
            print(f"  {ss.value:25s} φ={phi:6.1f}°  "
                  f"J_calc={j_calc:4.1f} Hz  J_exp={j_expected:4.1f} Hz")
    
    # Test 2: Couplages par acide aminé
    print("\n2. Couplages spécifiques par acide aminé:")
    print("-" * 50)
    
    for aa in [AminoAcid.GLY, AminoAcid.ALA, AminoAcid.PRO, 
               AminoAcid.VAL, AminoAcid.PHE]:
        print(f"\n  {aa.value}:")
        
        # Backbone
        j_hn_ha = calculate_3J_HN_HA(phi=-120, amino_acid=aa)
        print(f"    ³J(HN-Hα) = {j_hn_ha:.1f} Hz")
        
        # Sidechain
        sidechain = get_sidechain_couplings(aa, chi1=-60)
        for key, value in sidechain.items():
            print(f"    {key} = {value:.1f} Hz")
    
    # Test 3: Comparaison des paramètres Karplus
    print("\n3. Comparaison des jeux de paramètres Karplus:")
    print("-" * 50)
    
    test_angles = [-60, -120, 180]
    
    print(f"  {'φ (°)':<10}", end="")
    for name in KARPLUS_PARAMETERS['3J_HN_HA'].keys():
        print(f"{name:<12}", end="")
    print()
    
    for phi in test_angles:
        print(f"  {phi:<10}", end="")
        for name, params in KARPLUS_PARAMETERS['3J_HN_HA'].items():
            j = params.calculate(phi)
            print(f"{j:<12.2f}", end="")
        print()
    
    # Test 4: Prédiction structure secondaire
    print("\n4. Prédiction structure secondaire depuis ³J:")
    print("-" * 50)
    
    test_j = [3.5, 5.0, 7.0, 9.0]
    for j in test_j:
        ss = predict_secondary_structure_from_3J(j)
        print(f"  ³J = {j:.1f} Hz  →  {ss.value}")
    
    print("\n" + "=" * 70)
    print("TESTS TERMINÉS")
    print("=" * 70)


if __name__ == "__main__":
    _test_peptide_module()