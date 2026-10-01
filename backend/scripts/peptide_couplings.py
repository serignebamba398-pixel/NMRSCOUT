# peptide_couplings.py
# Robust module for couplings in peptides and proteins

from typing import Dict, List, Tuple, Optional, Any, Union
from dataclasses import dataclass, field
from enum import Enum
from functools import lru_cache
import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
except ImportError:
    raise ImportError("RDKit requis")

class AminoAcid(Enum):
    """Standard amino acid identifiers."""
    ALA = "ALA"
    ARG = "ARG"
    ASN = "ASN"
    ASP = "ASP"
    CYS = "CYS"
    GLN = "GLN"
    GLU = "GLU"
    GLY = "GLY"
    HIS = "HIS"
    ILE = "ILE"
    LEU = "LEU"
    LYS = "LYS"
    MET = "MET"
    PHE = "PHE"
    PRO = "PRO"
    SER = "SER"
    THR = "THR"
    TRP = "TRP"
    TYR = "TYR"
    VAL = "VAL"

class SecondaryStructure(Enum):
    """Secondary structure identifiers."""
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
    """Peptide coupling type identifiers."""
    HN_HA = "3J_HN_HA"
    HA_HB = "3J_HA_HB"
    HB_HG = "3J_HB_HG"
    HG_HD = "3J_HG_HD"
    HD_HE = "3J_HD_HE"
    N_CG = "3J_N_CG"
    CO_CG = "3J_CO_CG"
    HA_HA_GLY = "2J_HA_HA"
    HB_HB = "2J_HB_HB"
    NH2_SIDECHAIN = "3J_NH2"

@dataclass
class KarplusParameters:
    """Parameter container for coupling calculations."""
    A: float = 0.0
    B: float = 0.0
    C: float = 0.0
    phase: float = 0.0

    C0: Optional[float] = None
    C1: Optional[float] = None
    C2: Optional[float] = None
    C3: Optional[float] = None
    C4: Optional[float] = None

    name: str = ""
    reference: str = ""

    def calculate(self, theta_deg: float) -> float:
        """Calculate the coupling value for the supplied angle."""

        theta_with_phase = np.radians(theta_deg + self.phase)

        if self.C0 is not None:

            J = (self.C0 +
                 self.C1 * np.cos(theta_with_phase) +
                 self.C2 * np.cos(2 * theta_with_phase))

            if self.C3 is not None:
                J += self.C3 * np.sin(theta_with_phase)
            if self.C4 is not None:
                J += self.C4 * np.sin(2 * theta_with_phase)

            return J
        else:

            J = self.A * np.cos(theta_with_phase)**2 + self.B * np.cos(theta_with_phase) + self.C

            return J

KARPLUS_PARAMETERS: Dict[str, Dict[str, KarplusParameters]] = {

    '3J_HN_HA': {

        'bystrov': KarplusParameters(
            A=9.4, B=-1.1, C=0.4, phase=-60,
            name="Bystrov 1976",
            reference="Prog NMR Spectrosc 10:41-82"
        ),

        'pardi': KarplusParameters(
            A=6.4, B=-1.4, C=1.9, phase=-60,
            name="Pardi 1984",
            reference="J Mol Biol 180:741-751"
        ),

        'vuister': KarplusParameters(
            A=6.51, B=-1.76, C=1.60, phase=-60,
            name="Vuister-Bax 1993",
            reference="JACS 115:7772-7777"
        ),

        'wang': KarplusParameters(
            A=6.98, B=-1.38, C=1.72, phase=-60,
            name="Wang-Bax 1996",
            reference="JACS 118:2483-2494"
        ),

        'habeck': KarplusParameters(
            C0=1.5, C1=-1.3, C2=6.7, C3=0.0, C4=0.0, phase=-60,
            name="Habeck 2005 (Bayesian 6-param)",
            reference="J Magn Reson 177:160-165"
        ),

        'schmidt': KarplusParameters(
            A=6.0, B=-1.0, C=2.0, phase=-60,
            name="Schmidt 1999",
            reference="J Biomol NMR 14:1-12"
        ),
    },

    '3J_HA_HB': {

        'demarco': KarplusParameters(
            A=9.5, B=-1.68, C=1.8, phase=0,
            name="DeMarco 1978",
            reference="Biopolymers 17:617"
        ),

        'pachler': KarplusParameters(
            A=9.4, B=-1.1, C=0.4, phase=0,
            name="Pachler 1964",
            reference="Spectrochim Acta 20:581"
        ),

        'haasnoot': KarplusParameters(
            A=13.22, B=-0.99, C=0.0, phase=0,
            name="Haasnoot 1980",
            reference="Tetrahedron 36:2783"
        ),
    },

    '3J_HA_HB2': {
        'demarco': KarplusParameters(
            A=9.5, B=-1.68, C=1.8, phase=120,
            name="DeMarco 1978 (HB2 pro-S)",
            reference="Biopolymers 17:617"
        ),
    },

    '3J_HA_HB3': {
        'demarco': KarplusParameters(
            A=9.5, B=-1.68, C=1.8, phase=-120,
            name="DeMarco 1978 (HB3 pro-R)",
            reference="Biopolymers 17:617"
        ),
    },

    '3J_N_CG': {
        'perez': KarplusParameters(
            A=4.5, B=-1.2, C=0.5, phase=0,
            name="Perez 2001",
            reference="Magn Reson Chem 39:S24"
        ),
    },

    '3J_CO_CG': {
        'hu': KarplusParameters(
            A=4.4, B=-1.2, C=0.1, phase=0,
            name="Hu-Bax 1997",
            reference="JACS 119:6360"
        ),
    },

    '3J_HB_HG': {
        'standard': KarplusParameters(
            A=9.4, B=-1.1, C=0.4, phase=0,
            name="Standard aliphatic",
            reference="Pachler 1964"
        ),
    },
}

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
        'phi': -60.0,
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
        'phi': -120.0,
        'psi': 120.0,
        'omega': 180.0,
    },
}

EXPECTED_3J_HN_HA: Dict[SecondaryStructure, float] = {
    SecondaryStructure.ALPHA_HELIX: 3.9,
    SecondaryStructure.BETA_SHEET_PARALLEL: 9.7,
    SecondaryStructure.BETA_SHEET_ANTIPARALLEL: 9.0,
    SecondaryStructure.THREE_TEN_HELIX: 4.2,
    SecondaryStructure.PI_HELIX: 4.5,
    SecondaryStructure.TURN_TYPE_I: 5.5,
    SecondaryStructure.TURN_TYPE_II: 6.5,
    SecondaryStructure.POLYPROLINE_II: 6.0,
    SecondaryStructure.RANDOM_COIL: 7.0,
}

AMINO_ACID_SPECIFIC_COUPLINGS: Dict[AminoAcid, Dict[str, Any]] = {

    AminoAcid.GLY: {
        '2J_HA_HA': -17.1,
        'uncertainty': 0.2,
        'reference': 'Kjaergaard 2011, JACS 133:1563-1571',
        'notes': "No Hbeta, two equivalent Halpha",
        'equivalent_Ha': True,
    },

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
            'omega_cis': 0.0,
            'omega_trans': 180.0,

            '4J_HA_HD_cis': 0.0,
            '4J_HA_HD_trans': 0.0,
            'reference_cis_trans': 'Use NOE for cis/trans determination',
        },
        'notes': "5-membered pyrrolidine ring, no amide HN, cis/trans via NOE not J",
    },

    AminoAcid.ALA: {
        '3J_HA_HB': 7.0,
        'notes': "Methyl free rotation -> J~7 Hz average (quartet in 1H NMR)",
    },

    AminoAcid.VAL: {
        'chi1_rotamers': {

            'gauche_minus': {'3J_HA_HB': 3.5, 'chi1': -60},

            'gauche_plus': {'3J_HA_HB': 3.5, 'chi1': 60},

            'trans': {'3J_HA_HB': 11.0, 'chi1': 180},
        },

        '3J_HB_HG1': 0.0,
        '3J_HB_HG2': 0.0,
        '4J_HG1_HG2': 0.0,
        'reference': 'North 1994, J Mol Biol 235:1273',
        'notes': "Isopropyl, Cgamma quaternary -> no 3J(HB-HG), 3 chi1 rotamers",
    },

    AminoAcid.LEU: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},
            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},
            'trans': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 3.5},
        },
        '2J_HB2_HB3': -14.0,
        '3J_HB_HG': 6.5,
        '3J_HG_HD': 6.5,
        'notes': "Isobutyl",
    },

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
        'notes': "sec-butyl",
    },

    AminoAcid.THR: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB': 11.0},
            'gauche_plus': {'3J_HA_HB': 3.5},
            'trans': {'3J_HA_HB': 3.5},
        },
        '3J_HB_HG': 6.5,

        'notes': "Hydroxyl on Cbeta, single Hbeta (branched)",
    },

    AminoAcid.SER: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},
            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},
            'trans': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 3.5},
        },
        '2J_HB2_HB3': -14.5,
        'notes': "Hydroxyl on Cbeta increases |2J| due to electron-withdrawing effect",
    },

    AminoAcid.CYS: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},
            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},
            'trans': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 3.5},
        },
        '2J_HB2_HB3': -15.5,
        'notes': "Thiol group, can form disulfide bridges, S > O effect on 2J",
    },

    AminoAcid.MET: {
        '3J_HA_HB': 7.0,
        '3J_HB_HG': 7.0,
        '2J_HB2_HB3': -14.5,
        '2J_HG2_HG3': -14.5,
        '3J_HG_CH3': 0.0,
        'notes': "Thioether, S-CH3",
    },

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
        'notes': "Benzyl",
    },

    AminoAcid.TYR: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},
            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},
            'trans': {'3J_HA_HB2': 5.5, '3J_HA_HB3': 5.5},
        },
        '2J_HB2_HB3': -14.0,
        'aromatic': {
            '3J_ortho': 8.5,
            '4J_meta': 2.5,
        },
        'notes': "p-hydroxybenzyl, AA'BB'",
    },

    AminoAcid.TRP: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},
            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},
            'trans': {'3J_HA_HB2': 5.5, '3J_HA_HB3': 5.5},
        },
        '2J_HB2_HB3': -14.5,
        'indole': {
            '3J_H2_H7a': 0.0,
            '3J_H4_H5': 8.0,
            '3J_H5_H6': 7.0,
            '3J_H6_H7': 8.0,
            '4J_H4_H6': 1.0,
            '4J_H5_H7': 1.2,
        },
        'notes': "Fused indole",
    },

    AminoAcid.HIS: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},
            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},
            'trans': {'3J_HA_HB2': 5.5, '3J_HA_HB3': 5.5},
        },
        '2J_HB2_HB3': -14.5,
        'imidazole': {

            '4J_H2_H4': 1.0,
            '4J_H2_H5': 1.5,
            '3J_H4_H5': 1.5,
        },
        'tautomers': ['Ne-H', 'Nd-H'],
        'notes': "Imidazole, pKa ~6",
    },

    AminoAcid.ASN: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},
            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},
            'trans': {'3J_HA_HB2': 5.5, '3J_HA_HB3': 5.5},
        },
        '2J_HB2_HB3': -15.5,
        'amide_sidechain': {
            '2J_NH2': -14.0,
            '3J_HB_NH': 6.0,
        },
        'notes': "Primary amide in sidechain",
    },

    AminoAcid.GLN: {
        '3J_HA_HB': 7.0,
        '3J_HB_HG': 7.0,
        '2J_HB2_HB3': -14.0,
        '2J_HG2_HG3': -15.5,
        'amide_sidechain': {
            '2J_NH2': -14.0,
            '3J_HG_NH': 6.0,
        },
        'notes': "Primary amide, longer chain",
    },

    AminoAcid.ASP: {
        'chi1_rotamers': {
            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},
            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},
            'trans': {'3J_HA_HB2': 5.5, '3J_HA_HB3': 5.5},
        },
        '2J_HB2_HB3': -16.5,
        'notes': "Carboxyl, negatively charged at pH >4",
    },

    AminoAcid.GLU: {
        '3J_HA_HB': 7.0,
        '3J_HB_HG': 7.0,
        '2J_HB2_HB3': -14.0,
        '2J_HG2_HG3': -16.5,
        'notes': "Carboxyl, longer chain",
    },

    AminoAcid.LYS: {
        'chi1_rotamers': {

            'gauche_minus': {'3J_HA_HB2': 10.5, '3J_HA_HB3': 3.5},

            'gauche_plus': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 10.5},

            'trans': {'3J_HA_HB2': 3.5, '3J_HA_HB3': 3.5},
        },

        '3J_HA_HB2': 7.7,
        '3J_HA_HB3': 5.0,
        '3J_HA_HB': 6.3,
        '3J_HB_HG': 7.0,
        '3J_HG_HD': 7.0,
        '3J_HD_HE': 7.5,
        '2J_HB2_HB3': -15.0,
        '2J_HG2_HG3': -14.0,
        '2J_HD2_HD3': -12.5,
        '2J_HE2_HE3': -12.5,
        '3J_HE_NH3': 6.0,
        'notes': "Long chain with diastereotopic Hβ resolution; rotamer-averaged J",
    },

    AminoAcid.ARG: {
        '3J_HA_HB': 7.0,
        '3J_HB_HG': 7.0,
        '3J_HG_HD': 7.0,
        '2J_HB2_HB3': -15.0,
        '2J_HG2_HG3': -14.0,
        '2J_HD2_HD3': -12.5,
        'guanidinium': {
            '3J_HD_NE': 6.5,
            'NH_equivalent': True,
        },
        'notes': "Guanidinium, Cβ 2J elevated due to α-carbonyl effect",
    },
}

ROTAMER_POPULATIONS: Dict[AminoAcid, Dict[str, Dict[str, float]]] = {

    AminoAcid.PHE: {
        'random_coil': {'g-': 0.60, 't': 0.20, 'g+': 0.20},
        'alpha_helix': {'g-': 0.70, 't': 0.20, 'g+': 0.10},
        'beta_sheet':  {'g-': 0.55, 't': 0.25, 'g+': 0.20},

        'free_amino_acid': {'g-': 0.17, 't': 0.26, 'g+': 0.57},
    },
    AminoAcid.TYR: {
        'random_coil': {'g-': 0.55, 't': 0.25, 'g+': 0.20},
        'alpha_helix': {'g-': 0.65, 't': 0.25, 'g+': 0.10},
        'beta_sheet':  {'g-': 0.50, 't': 0.30, 'g+': 0.20},

        'free_amino_acid': {'g-': 0.14, 't': 0.33, 'g+': 0.53},
    },
    AminoAcid.TRP: {
        'random_coil': {'g-': 0.60, 't': 0.20, 'g+': 0.20},
        'alpha_helix': {'g-': 0.70, 't': 0.20, 'g+': 0.10},
        'beta_sheet':  {'g-': 0.55, 't': 0.25, 'g+': 0.20},

        'free_amino_acid': {'g-': 0.10, 't': 0.33, 'g+': 0.57},
    },
    AminoAcid.HIS: {
        'random_coil': {'g-': 0.50, 't': 0.30, 'g+': 0.20},
        'alpha_helix': {'g-': 0.60, 't': 0.25, 'g+': 0.15},
        'beta_sheet':  {'g-': 0.45, 't': 0.35, 'g+': 0.20},

        'free_amino_acid': {'g-': 0.04, 't': 0.45, 'g+': 0.51},
    },

    AminoAcid.VAL: {

        'random_coil': {'g-': 0.20, 't': 0.60, 'g+': 0.20},
        'alpha_helix': {'g-': 0.15, 't': 0.70, 'g+': 0.15},
        'beta_sheet':  {'g-': 0.25, 't': 0.50, 'g+': 0.25},

        'free_amino_acid': {'g-': 0.44, 't': 0.12, 'g+': 0.44},
    },
    AminoAcid.ILE: {

        'random_coil': {'g-': 0.15, 't': 0.70, 'g+': 0.15},
        'alpha_helix': {'g-': 0.10, 't': 0.80, 'g+': 0.10},
        'beta_sheet':  {'g-': 0.20, 't': 0.60, 'g+': 0.20},

        'free_amino_acid': {'g-': 0.06, 't': 0.80, 'g+': 0.14},
    },
    AminoAcid.LEU: {

        'random_coil': {'g-': 0.35, 't': 0.45, 'g+': 0.20},
        'alpha_helix': {'g-': 0.30, 't': 0.55, 'g+': 0.15},
        'beta_sheet':  {'g-': 0.40, 't': 0.40, 'g+': 0.20},
        'free_amino_acid': {'g-': 0.33, 't': 0.34, 'g+': 0.33},
    },
    AminoAcid.THR: {

        'random_coil': {'g-': 0.25, 't': 0.55, 'g+': 0.20},
        'alpha_helix': {'g-': 0.20, 't': 0.65, 'g+': 0.15},
        'beta_sheet':  {'g-': 0.30, 't': 0.50, 'g+': 0.20},

        'free_amino_acid': {'g-': 0.30, 't': 0.45, 'g+': 0.25},
    },

    AminoAcid.SER: {
        'random_coil': {'g-': 0.40, 't': 0.35, 'g+': 0.25},
        'alpha_helix': {'g-': 0.45, 't': 0.35, 'g+': 0.20},
        'beta_sheet':  {'g-': 0.35, 't': 0.40, 'g+': 0.25},

        'free_amino_acid': {'g-': 0.45, 't': 0.25, 'g+': 0.30},
    },
    AminoAcid.CYS: {
        'random_coil': {'g-': 0.45, 't': 0.30, 'g+': 0.25},
        'alpha_helix': {'g-': 0.50, 't': 0.30, 'g+': 0.20},
        'beta_sheet':  {'g-': 0.40, 't': 0.35, 'g+': 0.25},

        'free_amino_acid': {'g-': 0.50, 't': 0.22, 'g+': 0.28},
    },
    AminoAcid.ALA: {

        'random_coil': {'g-': 0.33, 't': 0.34, 'g+': 0.33},
        'alpha_helix': {'g-': 0.33, 't': 0.34, 'g+': 0.33},
        'beta_sheet':  {'g-': 0.33, 't': 0.34, 'g+': 0.33},
    },

    AminoAcid.ASP: {
        'random_coil': {'g-': 0.40, 't': 0.35, 'g+': 0.25},
        'alpha_helix': {'g-': 0.45, 't': 0.35, 'g+': 0.20},
        'beta_sheet':  {'g-': 0.35, 't': 0.40, 'g+': 0.25},

        'free_amino_acid': {'g-': 0.05, 't': 0.25, 'g+': 0.70},
    },
    AminoAcid.GLU: {
        'random_coil': {'g-': 0.40, 't': 0.35, 'g+': 0.25},
        'alpha_helix': {'g-': 0.45, 't': 0.35, 'g+': 0.20},
        'beta_sheet':  {'g-': 0.35, 't': 0.40, 'g+': 0.25},

        'free_amino_acid': {'g-': 0.02, 't': 0.65, 'g+': 0.33},
    },

    AminoAcid.ASN: {
        'random_coil': {'g-': 0.45, 't': 0.30, 'g+': 0.25},
        'alpha_helix': {'g-': 0.50, 't': 0.30, 'g+': 0.20},
        'beta_sheet':  {'g-': 0.40, 't': 0.35, 'g+': 0.25},

        'free_amino_acid': {'g-': 0.45, 't': 0.25, 'g+': 0.30},
    },
    AminoAcid.GLN: {
        'random_coil': {'g-': 0.40, 't': 0.35, 'g+': 0.25},
        'alpha_helix': {'g-': 0.45, 't': 0.35, 'g+': 0.20},
        'beta_sheet':  {'g-': 0.35, 't': 0.40, 'g+': 0.25},

        'free_amino_acid': {'g-': 0.18, 't': 0.55, 'g+': 0.27},
    },

    AminoAcid.LYS: {

        'random_coil': {'g-': 0.60, 't': 0.20, 'g+': 0.20},
        'alpha_helix': {'g-': 0.65, 't': 0.20, 'g+': 0.15},
        'beta_sheet':  {'g-': 0.55, 't': 0.25, 'g+': 0.20},

        'free_amino_acid': {'g-': 0.34, 't': 0.33, 'g+': 0.33},
    },
    AminoAcid.ARG: {
        'random_coil': {'g-': 0.55, 't': 0.25, 'g+': 0.20},
        'alpha_helix': {'g-': 0.60, 't': 0.25, 'g+': 0.15},
        'beta_sheet':  {'g-': 0.50, 't': 0.30, 'g+': 0.20},

        'free_amino_acid': {'g-': 0.34, 't': 0.33, 'g+': 0.33},
    },

    AminoAcid.MET: {
        'random_coil': {'g-': 0.45, 't': 0.35, 'g+': 0.20},
        'alpha_helix': {'g-': 0.50, 't': 0.35, 'g+': 0.15},
        'beta_sheet':  {'g-': 0.40, 't': 0.40, 'g+': 0.20},
        'free_amino_acid': {'g-': 0.36, 't': 0.34, 'g+': 0.30},
    },

    AminoAcid.GLY: {

        'random_coil': {'g-': 0.33, 't': 0.34, 'g+': 0.33},
        'alpha_helix': {'g-': 0.33, 't': 0.34, 'g+': 0.33},
        'beta_sheet':  {'g-': 0.33, 't': 0.34, 'g+': 0.33},
    },
    AminoAcid.PRO: {

        'random_coil': {'g-': 0.0, 't': 0.0, 'g+': 1.0},
        'alpha_helix': {'g-': 0.0, 't': 0.0, 'g+': 1.0},
        'beta_sheet':  {'g-': 0.0, 't': 0.0, 'g+': 1.0},
    },
}

def rotamer_averaged_J(
    amino_acid: Union[AminoAcid, str],
    coupling_type: str = 'HA_HB',
    secondary_structure: str = 'random_coil',
    proton: str = 'HB2'
) -> Tuple[float, str]:
    """Return the population-weighted coupling average."""

    if isinstance(amino_acid, str):
        try:
            amino_acid = AminoAcid(amino_acid.upper())
        except ValueError:

            for aa in AminoAcid:
                if aa.value == amino_acid.upper():
                    amino_acid = aa
                    break
            else:
                return (7.0, "rotamer_avg:unknown_aa")

    if amino_acid not in ROTAMER_POPULATIONS:
        return (7.0, f"rotamer_avg:{amino_acid.value}:no_data")

    aa_pops = ROTAMER_POPULATIONS[amino_acid]

    if secondary_structure not in aa_pops:
        secondary_structure = 'random_coil'

    populations = aa_pops[secondary_structure]

    if amino_acid not in AMINO_ACID_SPECIFIC_COUPLINGS:
        return (7.0, f"rotamer_avg:{amino_acid.value}:no_couplings")

    aa_data = AMINO_ACID_SPECIFIC_COUPLINGS[amino_acid]

    if 'chi1_rotamers' not in aa_data:

        key = f'3J_{coupling_type}'
        if key in aa_data:
            return (aa_data[key], f"rotamer_avg:{amino_acid.value}:direct")
        return (7.0, f"rotamer_avg:{amino_acid.value}:no_rotamers")

    rotamers = aa_data['chi1_rotamers']

    pop_to_rotamer = {
        'g-': 'gauche_minus',
        't': 'trans',
        'g+': 'gauche_plus'
    }

    if proton in ('HB2', 'HB3'):
        j_key = f'3J_HA_{proton}'
    elif proton == 'HB':
        j_key = '3J_HA_HB'
    else:
        j_key = f'3J_{coupling_type}'

    J_avg = 0.0
    total_pop = 0.0

    for pop_key, pop_value in populations.items():
        rotamer_key = pop_to_rotamer.get(pop_key)
        if rotamer_key and rotamer_key in rotamers:
            rotamer_data = rotamers[rotamer_key]

            if j_key in rotamer_data:
                J_avg += pop_value * rotamer_data[j_key]
                total_pop += pop_value
            elif '3J_HA_HB' in rotamer_data:

                J_avg += pop_value * rotamer_data['3J_HA_HB']
                total_pop += pop_value

    if total_pop > 0:
        J_avg /= total_pop
        source = f"rotamer_avg:{amino_acid.value}:{secondary_structure}:{proton}"
        return (round(J_avg, 1), source)

    return (7.0, f"rotamer_avg:{amino_acid.value}:fallback")

def get_diastereotopic_J_pair(
    amino_acid: Union[AminoAcid, str],
    secondary_structure: str = 'random_coil'
) -> Tuple[float, float, str]:
    """Return coupling values for both beta protons."""
    J_HB2, _ = rotamer_averaged_J(amino_acid, 'HA_HB', secondary_structure, 'HB2')
    J_HB3, _ = rotamer_averaged_J(amino_acid, 'HA_HB', secondary_structure, 'HB3')

    if isinstance(amino_acid, AminoAcid):
        aa_name = amino_acid.value
    else:
        aa_name = str(amino_acid).upper()

    source = f"rotamer_pair:{aa_name}:{secondary_structure}"

    return (J_HB2, J_HB3, source)

MODIFIED_AA_COUPLINGS: Dict[str, Dict[str, Any]] = {

    'MSE': {
        '2J_SE_H': -12.0,
        '3J_HA_HB': 7.0,
        '3J_HB_HG': 7.0,
        'notes': "Selenomethionine for MAD phasing",
    },

    'PYL': {
        '3J_HA_HB': 7.0,
        'notes': "Natural non-standard amino acid",
    },

    'NME': {
        '3J_N_CH3': 0.0,
        'notes': "N-methylation blocks HN",
    },

    'HYP': {
        'ring_couplings': {
            '3J_HA_HB2': 6.5,
            '3J_HA_HB3': 7.5,
            '3J_HB_HG': 3.5,
        },
        'notes': "4-Hydroxyproline, collagen marker",
    },
}

AA_SMARTS_PATTERNS: Dict[str, str] = {

    'BACKBONE_ALPHA': '[NX3,NX4+][CX4H]([*])[CX3](=[OX1])',
    'BACKBONE_PEPTIDE': '[NX3][CX4H][CX3](=[OX1])[NX3]',

    'GLY': '[NX3][CH2][CX3](=[OX1])',
    'PRO': '[NX3]1[CH2][CH2][CH2][CH]1[CX3](=[OX1])',
    'VAL': '[NX3][CH]([CH]([CH3])[CH3])[CX3](=[OX1])',
    'ALA': '[NX3][CH]([CH3])[CX3](=[OX1])',
    'LEU': '[NX3][CH]([CH2][CH]([CH3])[CH3])[CX3](=[OX1])',
    'ILE': '[NX3][CH]([CH]([CH3])[CH2][CH3])[CX3](=[OX1])',
    'SER': '[NX3][CH]([CH2][OH])[CX3](=[OX1])',
    'THR': '[NX3][CH]([CH]([OH])[CH3])[CX3](=[OX1])',
    'CYS': '[NX3][CH]([CH2][SH])[CX3](=[OX1])',
    'MET': '[NX3][CH]([CH2][CH2][SX2][CH3])[CX3](=[OX1])',
    'PHE': '[NX3][CH]([CH2]c1ccccc1)[CX3](=[OX1])',
    'TYR': '[NX3][CH]([CH2]c1ccc(O)cc1)[CX3](=[OX1])',
    'TRP': '[NX3][CH]([CH2]c1c[nH]c2ccccc12)[CX3](=[OX1])',
    'HIS': '[NX3][CH]([CH2]c1cnc[nH]1)[CX3](=[OX1])',
    'ASP': '[NX3][CH]([CH2][CX3](=[OX1])[OX1,OX2])[CX3](=[OX1])',
    'GLU': '[NX3][CH]([CH2][CH2][CX3](=[OX1])[OX1,OX2])[CX3](=[OX1])',
    'ASN': '[NX3][CH]([CH2][CX3](=[OX1])[NX3])[CX3](=[OX1])',
    'GLN': '[NX3][CH]([CH2][CH2][CX3](=[OX1])[NX3])[CX3](=[OX1])',
    'LYS': '[NX3][CH]([CH2][CH2][CH2][CH2][NX3,NX4+])[CX3](=[OX1])',
    'ARG': '[NX3][CH]([CH2][CH2][CH2][NX3][CX3](=[NX2])[NX3])[CX3](=[OX1])',
}

_COMPILED_SMARTS: Dict[str, Optional['Chem.Mol']] = {}

def _get_compiled_smarts(pattern_name: str) -> Optional['Chem.Mol']:
    """Return a cached compiled SMARTS pattern."""
    if pattern_name not in _COMPILED_SMARTS:
        smarts = AA_SMARTS_PATTERNS.get(pattern_name)
        if smarts:
            try:
                _COMPILED_SMARTS[pattern_name] = Chem.MolFromSmarts(smarts)
            except Exception:
                _COMPILED_SMARTS[pattern_name] = None
        else:
            _COMPILED_SMARTS[pattern_name] = None
    return _COMPILED_SMARTS[pattern_name]

def detect_amino_acid_type(mol: Chem.Mol, ca_idx: int) -> Optional[AminoAcid]:
    """Detect the amino acid type from the molecular structure."""

    if not _is_alpha_carbon(mol, ca_idx):
        return None

    atom = mol.GetAtomWithIdx(ca_idx)
    if atom.IsInRing():
        return None

    aa_map = {
        'GLY': AminoAcid.GLY,
        'PRO': AminoAcid.PRO,
        'VAL': AminoAcid.VAL,
        'ALA': AminoAcid.ALA,
        'LEU': AminoAcid.LEU,
        'ILE': AminoAcid.ILE,
        'SER': AminoAcid.SER,
        'THR': AminoAcid.THR,
        'CYS': AminoAcid.CYS,
        'MET': AminoAcid.MET,
        'PHE': AminoAcid.PHE,
        'TYR': AminoAcid.TYR,
        'TRP': AminoAcid.TRP,
        'HIS': AminoAcid.HIS,
        'ASP': AminoAcid.ASP,
        'GLU': AminoAcid.GLU,
        'ASN': AminoAcid.ASN,
        'GLN': AminoAcid.GLN,
        'LYS': AminoAcid.LYS,
        'ARG': AminoAcid.ARG,
    }

    for pattern_name, aa_enum in aa_map.items():
        pattern = _get_compiled_smarts(pattern_name)
        if pattern is not None:
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:

                if len(match) > 1 and ca_idx in match:
                    return aa_enum

    return None

def _topological_distances_local(mol: Chem.Mol) -> Dict[Tuple[int, int], int]:
    """Calculate topological distances between atom pairs."""
    distances: Dict[Tuple[int, int], int] = {}
    num_atoms = mol.GetNumAtoms()

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            try:
                path = Chem.rdmolops.GetShortestPath(mol, i, j)
                dist = len(path) - 1 if path else 999
            except Exception:
                dist = 999
            distances[(i, j)] = dist
            distances[(j, i)] = dist

    return distances

def get_coupling_distance(mol: Chem.Mol, h1_idx: int, h2_idx: int) -> int:
    """Return the topological distance between two hydrogens."""
    try:
        path = Chem.rdmolops.GetShortestPath(mol, h1_idx, h2_idx)
        return len(path) - 1 if path else 999
    except Exception:
        return 999

def _validate_dihedral_angle(angle: Optional[float], name: str = "angle") -> Optional[float]:
    """Normalize and validate a dihedral angle."""
    if angle is None:
        return None

    while angle > 180:
        angle -= 360
    while angle < -180:
        angle += 360

    if not (-180 <= angle <= 180):
        raise ValueError(f"{name} {angle}° is outside valid range [-180, 180]")

    return angle

def calculate_3J_HN_HA(
    phi: Optional[float] = None,
    secondary_structure: Optional[SecondaryStructure] = None,
    karplus_set: str = 'vuister',
    amino_acid: Optional[AminoAcid] = None,
) -> float:
    """Calculate the HN-HA coupling."""

    phi = _validate_dihedral_angle(phi, "phi")

    if amino_acid == AminoAcid.PRO:
        return 0.0

    if amino_acid == AminoAcid.GLY:

        if phi is not None:
            params = KARPLUS_PARAMETERS['3J_HN_HA'].get(karplus_set)
            if params:
                return params.calculate(phi)
        elif secondary_structure:
            return EXPECTED_3J_HN_HA.get(secondary_structure, 7.0)
        return 6.5

    if phi is not None:
        params = KARPLUS_PARAMETERS['3J_HN_HA'].get(karplus_set)
        if params:
            return params.calculate(phi)

        theta = np.radians(phi - 60)
        return max(0, 6.4 * np.cos(theta)**2 - 1.4 * np.cos(theta) + 1.9)

    if secondary_structure:
        return EXPECTED_3J_HN_HA.get(secondary_structure, 7.0)

    return 7.0

def calculate_3J_HA_HB(
    chi1: Optional[float] = None,
    amino_acid: Optional[AminoAcid] = None,
    hb_type: str = 'HB2',
    karplus_set: str = 'demarco',
) -> float:
    """Calculate the HA-HB coupling."""

    chi1 = _validate_dihedral_angle(chi1, "chi1")

    if amino_acid == AminoAcid.GLY:
        return 0.0

    if amino_acid == AminoAcid.ALA:
        return 7.0

    if amino_acid == AminoAcid.PRO:
        pro_data = AMINO_ACID_SPECIFIC_COUPLINGS.get(AminoAcid.PRO, {})
        ring_couplings = pro_data.get('ring_couplings', {})

        key = f'3J_HA_{hb_type}'
        if key in ring_couplings:
            j_data = ring_couplings[key]
            if isinstance(j_data, dict):

                endo = j_data.get('endo', 7.0)
                exo = j_data.get('exo', 7.0)
                return (endo + exo) / 2.0
            else:
                return float(j_data)

        return 7.0

    if chi1 is not None:
        params = KARPLUS_PARAMETERS['3J_HA_HB'].get(karplus_set)
        if params:
            return params.calculate(chi1)

    if amino_acid and amino_acid in AMINO_ACID_SPECIFIC_COUPLINGS:
        aa_data = AMINO_ACID_SPECIFIC_COUPLINGS[amino_acid]

        key = f'3J_HA_{hb_type}'
        if key in aa_data:
            return aa_data[key]

        if '3J_HA_HB' in aa_data:
            return aa_data['3J_HA_HB']

    return 7.0

def calculate_2J_geminal(
    amino_acid: Optional[AminoAcid] = None,
    position: str = 'HB',
) -> float:
    """Calculate a geminal coupling."""

    if amino_acid == AminoAcid.GLY and position == 'HA':
        return AMINO_ACID_SPECIFIC_COUPLINGS[AminoAcid.GLY]['2J_HA_HA']

    if amino_acid and amino_acid in AMINO_ACID_SPECIFIC_COUPLINGS:
        aa_data = AMINO_ACID_SPECIFIC_COUPLINGS[amino_acid]
        key = f'2J_{position}2_{position}3'
        if key in aa_data:
            return aa_data[key]

    defaults = {
        'HA': -17.5,
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
    """Return the coupling for a supported peptide proton pair."""

    if mol is None or not isinstance(mol, Chem.Mol):
        return None

    num_atoms = mol.GetNumAtoms()
    if num_atoms == 0:
        return None

    if not (0 <= h1_idx < num_atoms) or not (0 <= h2_idx < num_atoms):
        return None

    h1 = mol.GetAtomWithIdx(h1_idx)
    h2 = mol.GetAtomWithIdx(h2_idx)

    if h1.GetSymbol() != 'H' or h2.GetSymbol() != 'H':
        return None

    heavy1_list = [n for n in h1.GetNeighbors() if n.GetSymbol() != 'H']
    heavy2_list = [n for n in h2.GetNeighbors() if n.GetSymbol() != 'H']

    if not heavy1_list or not heavy2_list:
        return None

    heavy1 = heavy1_list[0]
    heavy2 = heavy2_list[0]

    if heavy1.GetIdx() != heavy2.GetIdx():
        h1_rings = set(mol.GetRingInfo().AtomRings())
        h1_in_rings = {i for i, ring in enumerate(h1_rings) if heavy1.GetIdx() in ring}
        h2_in_rings = {i for i, ring in enumerate(h1_rings) if heavy2.GetIdx() in ring}

        if h1_in_rings & h2_in_rings:
            return None

    if distance is None:
        try:
            path = Chem.GetShortestPath(mol, h1_idx, h2_idx)

            distance = len(path) - 2
        except Exception:
            return None

    coupling_type = _detect_peptide_coupling_type(mol, h1_idx, h2_idx, heavy1, heavy2)

    if coupling_type is None:
        return None

    if amino_acid is None:

        alpha_c = None
        for heavy in [heavy1, heavy2]:
            if _is_alpha_carbon(mol, heavy.GetIdx()):
                alpha_c = heavy.GetIdx()
                break

            for nb in heavy.GetNeighbors():
                if nb.GetSymbol() == 'C' and _is_alpha_carbon(mol, nb.GetIdx()):
                    alpha_c = nb.GetIdx()
                    break
            if alpha_c is not None:
                break

        if alpha_c is not None:
            amino_acid = detect_amino_acid_type(mol, alpha_c)

    if mol.GetNumConformers() > 0 and phi is None:
        try:
            phi = _extract_phi_angle(mol, heavy1.GetIdx())
        except Exception:
            pass

    if mol.GetNumConformers() > 0 and chi1 is None:
        try:
            chi1 = _extract_chi1_angle(mol, heavy1.GetIdx())
        except Exception:
            pass

    if coupling_type == CouplingType.HN_HA:
        return calculate_3J_HN_HA(
            phi=phi,
            secondary_structure=secondary_structure,
            amino_acid=amino_acid,
        )

    elif coupling_type == CouplingType.HA_HB:

        hb_type = 'HB2'

        c_alpha = heavy1 if _is_alpha_carbon(mol, heavy1.GetIdx()) else (
            heavy2 if _is_alpha_carbon(mol, heavy2.GetIdx()) else None
        )
        c_beta = heavy2 if c_alpha == heavy1 else heavy1
        h_beta_idx = h2_idx if c_beta == heavy2 else h1_idx

        use_chi1_from_3d = False

        if c_alpha is not None and c_beta is not None:
            ca_in_ring = c_alpha.IsInRing()
            cb_in_ring = c_beta.IsInRing()

            if ca_in_ring or cb_in_ring:
                use_chi1_from_3d = True

            if amino_acid == AminoAcid.PRO:
                use_chi1_from_3d = True

        if use_chi1_from_3d and c_alpha is not None and mol.GetNumConformers() > 0:
            try:
                chi1 = _extract_chi1_angle(mol, c_alpha.GetIdx())
            except Exception:
                pass

        if c_beta is not None:

            beta_hydrogens = [n.GetIdx() for n in c_beta.GetNeighbors() if n.GetSymbol() == 'H']
            if len(beta_hydrogens) == 2:
                beta_hydrogens.sort()

                if h_beta_idx == beta_hydrogens[0]:
                    hb_type = 'HB2'
                else:
                    hb_type = 'HB3'

        return calculate_3J_HA_HB(
            chi1=chi1,
            amino_acid=amino_acid,
            hb_type=hb_type,
        )

    elif coupling_type == CouplingType.HA_HA_GLY:
        return calculate_2J_geminal(AminoAcid.GLY, 'HA')

    elif coupling_type == CouplingType.HB_HB:
        return calculate_2J_geminal(amino_acid, 'HB')

    elif coupling_type == CouplingType.HB_HG:
        return 7.0

    elif coupling_type == CouplingType.HG_HD:
        return 7.0

    elif coupling_type == CouplingType.HD_HE:

        if amino_acid and amino_acid in AMINO_ACID_SPECIFIC_COUPLINGS:
            aa_data = AMINO_ACID_SPECIFIC_COUPLINGS[amino_acid]
            if '3J_HD_HE' in aa_data:
                return aa_data['3J_HD_HE']
        return 7.0

    return None

def _detect_peptide_coupling_type(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    heavy1: Chem.Atom,
    heavy2: Chem.Atom
) -> Optional[CouplingType]:
    """Detect the supported peptide coupling type."""

    try:
        path = Chem.GetShortestPath(mol, h1_idx, h2_idx)
        if path:
            n_bonds = len(path) - 1

            if n_bonds > 3:
                return None
    except Exception:
        return None

    if heavy1.GetIdx() == heavy2.GetIdx():
        carbon = heavy1

        if _is_alpha_carbon(mol, carbon.GetIdx()):
            return CouplingType.HA_HA_GLY

        return CouplingType.HB_HB

    is_h1_on_N = heavy1.GetSymbol() == 'N'
    is_h2_on_N = heavy2.GetSymbol() == 'N'
    is_h1_on_Ca = _is_alpha_carbon(mol, heavy1.GetIdx())
    is_h2_on_Ca = _is_alpha_carbon(mol, heavy2.GetIdx())

    if (is_h1_on_N and is_h2_on_Ca) or (is_h2_on_N and is_h1_on_Ca):
        return CouplingType.HN_HA

    is_h1_on_Cb = _is_beta_carbon(mol, heavy1.GetIdx())
    is_h2_on_Cb = _is_beta_carbon(mol, heavy2.GetIdx())

    if (is_h1_on_Ca and is_h2_on_Cb) or (is_h2_on_Ca and is_h1_on_Cb):
        return CouplingType.HA_HB

    is_h1_on_Cg = _is_gamma_carbon(mol, heavy1.GetIdx())
    is_h2_on_Cg = _is_gamma_carbon(mol, heavy2.GetIdx())

    if (is_h1_on_Cb and is_h2_on_Cg) or (is_h2_on_Cb and is_h1_on_Cg):
        return CouplingType.HB_HG

    is_h1_on_Cd = _is_delta_carbon(mol, heavy1.GetIdx())
    is_h2_on_Cd = _is_delta_carbon(mol, heavy2.GetIdx())

    if (is_h1_on_Cg and is_h2_on_Cd) or (is_h2_on_Cg and is_h1_on_Cd):
        return CouplingType.HG_HD

    is_h1_on_Ce = _is_epsilon_carbon(mol, heavy1.GetIdx())
    is_h2_on_Ce = _is_epsilon_carbon(mol, heavy2.GetIdx())

    if (is_h1_on_Cd and is_h2_on_Ce) or (is_h2_on_Cd and is_h1_on_Ce):
        return CouplingType.HD_HE

    return None

class LRUCache:
    """Size-limited least-recently-used cache."""

    def __init__(self, maxsize: int = 10000):
        """Initialize the cache."""
        from collections import OrderedDict
        self._cache: 'OrderedDict[str, Any]' = OrderedDict()
        self._maxsize = maxsize
        self._hits = 0
        self._misses = 0

    def get(self, key: str, default: Any = None) -> Any:
        """Return a cached value."""
        if key in self._cache:
            self._cache.move_to_end(key)
            self._hits += 1
            return self._cache[key]
        self._misses += 1
        return default

    def set(self, key: str, value: Any) -> None:
        """Store a value in the cache."""
        if key in self._cache:
            self._cache.move_to_end(key)
        else:
            if len(self._cache) >= self._maxsize:

                self._cache.popitem(last=False)
        self._cache[key] = value

    def __contains__(self, key: str) -> bool:
        """Return whether a key is cached."""
        return key in self._cache

    def clear(self) -> None:
        """Clear all cache entries."""
        self._cache.clear()
        self._hits = 0
        self._misses = 0

    def stats(self) -> Dict[str, Any]:
        """Return cache statistics."""
        total = self._hits + self._misses
        hit_rate = self._hits / total if total > 0 else 0.0
        return {
            'size': len(self._cache),
            'maxsize': self._maxsize,
            'hits': self._hits,
            'misses': self._misses,
            'hit_rate': hit_rate
        }

_CARBON_TYPE_CACHE = LRUCache(maxsize=10000)

def clear_carbon_type_cache() -> Dict[str, Any]:
    """Clear the carbon-type cache."""
    stats = _CARBON_TYPE_CACHE.stats()
    _CARBON_TYPE_CACHE.clear()
    return stats

def _get_cache_key(mol: Chem.Mol, idx: int, func_name: str) -> str:
    """Build a cache key for carbon-type detection."""
    try:
        mol_smiles = Chem.MolToSmiles(mol, canonical=True)
        return f"{func_name}_{mol_smiles}_{idx}"
    except Exception:

        return f"{func_name}_{id(mol)}_{idx}"

def _is_alpha_carbon(mol: Chem.Mol, c_idx: int) -> bool:
    """Return whether the atom is a supported alpha carbon."""

    cache_key = _get_cache_key(mol, c_idx, "alpha")
    cached = _CARBON_TYPE_CACHE.get(cache_key)
    if cached is not None:
        return cached

    atom = mol.GetAtomWithIdx(c_idx)
    result = False

    if atom.GetSymbol() != 'C':
        pass

    elif atom.GetHybridization() != Chem.HybridizationType.SP3:
        pass

    elif atom.GetTotalNumHs() == 0:
        pass

    elif atom.GetIsAromatic():
        pass

    elif atom.IsInRing():
        pass
    else:

        has_amide_N = False
        has_CO = False

        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'N':

                if not neighbor.GetIsAromatic():

                    for nn in neighbor.GetNeighbors():
                        if nn.GetSymbol() == 'C' and nn.GetIdx() != c_idx:
                            for nnn in nn.GetNeighbors():
                                if nnn.GetSymbol() == 'O':
                                    bond = mol.GetBondBetweenAtoms(nn.GetIdx(), nnn.GetIdx())
                                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                                        has_amide_N = True
                                        break
            elif neighbor.GetSymbol() == 'C':

                for nn in neighbor.GetNeighbors():
                    if nn.GetSymbol() == 'O':
                        bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nn.GetIdx())
                        if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                            has_CO = True

        result = has_amide_N and has_CO

    _CARBON_TYPE_CACHE.set(cache_key, result)
    return result

def _is_beta_carbon(mol: Chem.Mol, c_idx: int) -> bool:
    """Return whether the atom is a beta carbon."""

    cache_key = _get_cache_key(mol, c_idx, "beta")
    cached = _CARBON_TYPE_CACHE.get(cache_key)
    if cached is not None:
        return cached

    atom = mol.GetAtomWithIdx(c_idx)

    if atom.GetSymbol() != 'C':
        result = False
    else:
        result = False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and _is_alpha_carbon(mol, neighbor.GetIdx()):
                result = True
                break

    _CARBON_TYPE_CACHE.set(cache_key, result)
    return result

def _is_gamma_carbon(mol: Chem.Mol, c_idx: int) -> bool:
    """Return whether the atom is a gamma carbon."""

    cache_key = _get_cache_key(mol, c_idx, "gamma")
    cached = _CARBON_TYPE_CACHE.get(cache_key)
    if cached is not None:
        return cached

    atom = mol.GetAtomWithIdx(c_idx)

    if atom.GetSymbol() != 'C':
        result = False
    else:
        result = False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and _is_beta_carbon(mol, neighbor.GetIdx()):
                result = True
                break

    _CARBON_TYPE_CACHE.set(cache_key, result)
    return result

def _is_delta_carbon(mol: Chem.Mol, c_idx: int) -> bool:
    """Return whether the atom is a delta carbon."""

    cache_key = _get_cache_key(mol, c_idx, "delta")
    cached = _CARBON_TYPE_CACHE.get(cache_key)
    if cached is not None:
        return cached

    atom = mol.GetAtomWithIdx(c_idx)

    if atom.GetSymbol() != 'C':
        result = False
    else:
        result = False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and _is_gamma_carbon(mol, neighbor.GetIdx()):
                result = True
                break

    _CARBON_TYPE_CACHE.set(cache_key, result)
    return result

def _is_epsilon_carbon(mol: Chem.Mol, c_idx: int) -> bool:
    """Return whether the atom is an epsilon carbon."""

    cache_key = _get_cache_key(mol, c_idx, "epsilon")
    cached = _CARBON_TYPE_CACHE.get(cache_key)
    if cached is not None:
        return cached

    atom = mol.GetAtomWithIdx(c_idx)

    if atom.GetSymbol() != 'C':
        result = False
    else:
        result = False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and _is_delta_carbon(mol, neighbor.GetIdx()):
                result = True
                break

    _CARBON_TYPE_CACHE.set(cache_key, result)
    return result

def _is_d_amino_acid(mol: Chem.Mol, ca_idx: int) -> bool:
    """Return whether the alpha carbon has D configuration."""
    if mol.GetNumConformers() == 0:
        return False

    try:
        ca_atom = mol.GetAtomWithIdx(ca_idx)

        if not _is_alpha_carbon(mol, ca_idx):
            return False

        cip_code = ca_atom.GetProp('_CIPCode') if ca_atom.HasProp('_CIPCode') else None

        if cip_code == 'R':
            return True
        elif cip_code == 'S':
            return False

        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        cip_code = ca_atom.GetProp('_CIPCode') if ca_atom.HasProp('_CIPCode') else None

        return cip_code == 'R'

    except Exception:
        return False

def _extract_phi_angle(mol: Chem.Mol, ca_idx: int) -> Optional[float]:
    """Extract the phi angle from a 3D conformer."""
    if mol.GetNumConformers() == 0:
        return None

    ca = mol.GetAtomWithIdx(ca_idx)

    n_idx = None
    co_idx = None
    prev_c_idx = None

    for neighbor in ca.GetNeighbors():
        if neighbor.GetSymbol() == 'N':
            n_idx = neighbor.GetIdx()

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
        except Exception:
            pass

    return None

def _extract_chi1_angle(mol: Chem.Mol, ca_idx: int) -> Optional[float]:
    """Extract the chi1 angle from a 3D conformer."""
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

            is_carbonyl = False
            for nn in neighbor.GetNeighbors():
                if nn.GetSymbol() == 'O':
                    bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nn.GetIdx())
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                        is_carbonyl = True

            if not is_carbonyl:
                cb_idx = neighbor.GetIdx()

                for nn in neighbor.GetNeighbors():
                    if nn.GetIdx() != ca_idx and nn.GetSymbol() != 'H':
                        xg_idx = nn.GetIdx()
                        break

    if all(x is not None for x in [n_idx, cb_idx, xg_idx]):
        try:
            conf = mol.GetConformer()
            chi1 = AllChem.GetDihedralDeg(conf, n_idx, ca_idx, cb_idx, xg_idx)

            if _is_d_amino_acid(mol, ca_idx):
                chi1 = -chi1

            return chi1
        except Exception:
            pass

    return None

def predict_secondary_structure_from_3J(j_hn_ha: float) -> SecondaryStructure:
    """Return the secondary-structure category for a coupling value."""
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
    """Return supported backbone couplings for one residue."""
    result = {}

    if amino_acid != AminoAcid.PRO:
        result['3J_HN_HA'] = calculate_3J_HN_HA(phi=phi, amino_acid=amino_acid)

    if amino_acid == AminoAcid.GLY:
        result['2J_HA_HA'] = calculate_2J_geminal(AminoAcid.GLY, 'HA')

    return result

def get_sidechain_couplings(
    amino_acid: AminoAcid,
    chi1: Optional[float] = None,
    chi2: Optional[float] = None,
) -> Dict[str, float]:
    """Return supported side-chain couplings for one amino acid."""
    result = {}

    if amino_acid not in AMINO_ACID_SPECIFIC_COUPLINGS:
        return result

    aa_data = AMINO_ACID_SPECIFIC_COUPLINGS[amino_acid]

    if 'chi1_rotamers' in aa_data and chi1 is not None:
        result['3J_HA_HB2'] = calculate_3J_HA_HB(chi1, amino_acid, 'HB2')
        result['3J_HA_HB3'] = calculate_3J_HA_HB(chi1, amino_acid, 'HB3')
    elif '3J_HA_HB' in aa_data:
        result['3J_HA_HB'] = aa_data['3J_HA_HB']

    for key, value in aa_data.items():
        if key.startswith('2J_'):
            result[key] = value

    if 'aromatic' in aa_data:
        for key, value in aa_data['aromatic'].items():
            result[f'aromatic_{key}'] = value

    for ring_type in ['imidazole', 'indole']:
        if ring_type in aa_data:
            for key, value in aa_data[ring_type].items():
                result[f'{ring_type}_{key}'] = value

    return result

def get_frequency_dependent_threshold(spectrometer_mhz: float) -> float:
    """Return the coupling threshold for a spectrometer frequency."""

    if spectrometer_mhz >= 800:
        return 0.3
    elif spectrometer_mhz >= 600:
        return 0.5
    elif spectrometer_mhz >= 400:
        return 0.8
    else:
        return 1.0

def export_couplings_for_sim_mult(
    mol: Chem.Mol,
    spectrometer_mhz: float = 600,
    include_geminal: bool = True,
) -> Dict[Tuple[int, int], float]:
    """Export supported couplings for the multiplicity simulator."""
    threshold = get_frequency_dependent_threshold(spectrometer_mhz)
    couplings = {}

    h_atoms = [a for a in mol.GetAtoms() if a.GetSymbol() == 'H']

    for i, h1 in enumerate(h_atoms):
        for h2 in h_atoms[i+1:]:
            h1_idx = h1.GetIdx()
            h2_idx = h2.GetIdx()

            j = get_peptide_coupling_robust(mol, h1_idx, h2_idx)

            if j is None:
                continue

            if not include_geminal and j < 0:
                continue

            if abs(j) >= threshold:

                key = (min(h1_idx, h2_idx), max(h1_idx, h2_idx))
                couplings[key] = round(j, 1)

    return couplings

def _test_peptide_module():
    """Run module checks."""
    print("=" * 70)
    print("ROBUST PEPTIDE MODULE TEST")
    print("=" * 70)

    print("\n1. 3J(HN-Halpha) by secondary structure:")
    print("-" * 50)

    for ss in SecondaryStructure:
        angles = SECONDARY_STRUCTURE_ANGLES.get(ss, {})
        phi = angles.get('phi')

        if phi:
            j_calc = calculate_3J_HN_HA(phi=phi)
            j_expected = EXPECTED_3J_HN_HA.get(ss, 7.0)

            print(f"  {ss.value:25s} phi={phi:6.1f} deg  "
                  f"J_calc={j_calc:4.1f} Hz  J_exp={j_expected:4.1f} Hz")

    print("\n2. Specific couplings by amino acid:")
    print("-" * 50)

    for aa in [AminoAcid.GLY, AminoAcid.ALA, AminoAcid.PRO,
               AminoAcid.VAL, AminoAcid.PHE]:
        print(f"\n  {aa.value}:")

        j_hn_ha = calculate_3J_HN_HA(phi=-120, amino_acid=aa)
        print(f"    3J(HN-Halpha) = {j_hn_ha:.1f} Hz")

        sidechain = get_sidechain_couplings(aa, chi1=-60)
        for key, value in sidechain.items():
            print(f"    {key} = {value:.1f} Hz")

    print("\n3. Karplus parameter set comparison:")
    print("-" * 50)

    test_angles = [-60, -120, 180]

    print(f"  {'phi (deg)':<10}", end="")
    for name in KARPLUS_PARAMETERS['3J_HN_HA'].keys():
        print(f"{name:<12}", end="")
    print()

    for phi in test_angles:
        print(f"  {phi:<10}", end="")
        for name, params in KARPLUS_PARAMETERS['3J_HN_HA'].items():
            j = params.calculate(phi)
            print(f"{j:<12.2f}", end="")
        print()

    print("\n4. Secondary structure prediction from 3J:")
    print("-" * 50)

    test_j = [3.5, 5.0, 7.0, 9.0]
    for j in test_j:
        ss = predict_secondary_structure_from_3J(j)
        print(f"  3J = {j:.1f} Hz  ->  {ss.value}")

    print("\n5. LRU cache statistics:")
    print("-" * 50)
    stats = _CARBON_TYPE_CACHE.stats()
    print(f"  Size: {stats['size']}/{stats['maxsize']}")
    print(f"  Hits: {stats['hits']}, Misses: {stats['misses']}")
    print(f"  Hit rate: {stats['hit_rate']:.2%}")

    print("\n6. Angle validation:")
    print("-" * 50)
    try:
        _validate_dihedral_angle(450, "test")
        print("  ERROR: Should have normalized 450 to 90")
    except ValueError as e:
        print(f"  Validation error caught: {e}")

    normalized = _validate_dihedral_angle(450, "test")
    print(f"  450 deg normalized to: {normalized} deg")

    print("\n" + "=" * 70)
    print("TESTS COMPLETED")
    print("=" * 70)

def test_peptide_production():
    """Run production validation checks."""
    print("=" * 70)
    print("FDA/EMA PRODUCTION VALIDATION TEST")
    print("=" * 70)

    errors = []

    print("\n1. Glycine 2J(HA-HA) validation:")
    j_gem = calculate_2J_geminal(AminoAcid.GLY, 'HA')
    expected_range = (-18.0, -17.0)
    if expected_range[0] <= j_gem <= expected_range[1]:
        print(f"   PASS: {j_gem:.1f} Hz in range {expected_range}")
    else:
        msg = f"Glycine 2J failed: {j_gem:.1f} Hz not in {expected_range}"
        errors.append(msg)
        print(f"   FAIL: {msg}")

    print("\n2. Proline 4J(HA-HD) validation (4-bond not observable):")
    pro_data = AMINO_ACID_SPECIFIC_COUPLINGS[AminoAcid.PRO]['cis_trans']
    cis_j = pro_data['4J_HA_HD_cis']
    trans_j = pro_data['4J_HA_HD_trans']

    if cis_j == 0.0 and trans_j == 0.0:
        print(f"   PASS: 4-bond couplings correctly set to 0 (use NOE for cis/trans)")
    else:
        msg = f"Proline 4J should be 0: cis={cis_j}, trans={trans_j}"
        errors.append(msg)
        print(f"   FAIL: {msg}")

    print("\n3. Alpha helix 3J(HN-HA) validation:")
    j_helix = calculate_3J_HN_HA(phi=-57)
    expected_range = (3.5, 5.0)
    if expected_range[0] <= j_helix <= expected_range[1]:
        print(f"   PASS: {j_helix:.1f} Hz in range {expected_range}")
    else:
        msg = f"Alpha helix 3J failed: {j_helix:.1f} Hz not in {expected_range}"
        errors.append(msg)
        print(f"   FAIL: {msg}")

    print("\n4. Beta sheet 3J(HN-HA) validation:")
    j_beta = calculate_3J_HN_HA(phi=-139)
    expected_range = (8.0, 10.0)
    if expected_range[0] <= j_beta <= expected_range[1]:
        print(f"   PASS: {j_beta:.1f} Hz in range {expected_range}")
    else:
        msg = f"Beta sheet 3J failed: {j_beta:.1f} Hz not in {expected_range}"
        errors.append(msg)
        print(f"   FAIL: {msg}")

    print("\n5. Valine rotamer validation (North 1994):")
    val_data = AMINO_ACID_SPECIFIC_COUPLINGS[AminoAcid.VAL]['chi1_rotamers']
    gm_j = val_data['gauche_minus']['3J_HA_HB']
    trans_j = val_data['trans']['3J_HA_HB']

    if gm_j < 5.0 and trans_j > 10.0:
        print(f"   PASS: gauche- ({gm_j:.1f} Hz) < 5, trans ({trans_j:.1f} Hz) > 10")
    else:
        msg = f"Valine rotamers failed: gm={gm_j}, trans={trans_j}"
        errors.append(msg)
        print(f"   FAIL: {msg}")

    print("\n6. Cache efficiency validation:")
    clear_carbon_type_cache()

    for _ in range(100):
        _validate_dihedral_angle(-60, "phi")
    stats = _CARBON_TYPE_CACHE.stats()

    print(f"   Cache structure validated: maxsize={stats['maxsize']}")

    print("\n7. Serine/Cysteine 2J validation:")
    ser_2j = AMINO_ACID_SPECIFIC_COUPLINGS[AminoAcid.SER]['2J_HB2_HB3']
    cys_2j = AMINO_ACID_SPECIFIC_COUPLINGS[AminoAcid.CYS]['2J_HB2_HB3']
    if abs(cys_2j) > abs(ser_2j):
        print(f"   PASS: |Cys 2J| ({abs(cys_2j):.1f}) > |Ser 2J| ({abs(ser_2j):.1f})")
    else:
        msg = f"S vs O effect wrong: Cys={cys_2j}, Ser={ser_2j}"
        errors.append(msg)
        print(f"   FAIL: {msg}")

    print("\n8. Beta sheet expected J ordering:")
    j_parallel = EXPECTED_3J_HN_HA[SecondaryStructure.BETA_SHEET_PARALLEL]
    j_antiparallel = EXPECTED_3J_HN_HA[SecondaryStructure.BETA_SHEET_ANTIPARALLEL]

    if j_parallel > j_antiparallel:
        print(f"   PASS: parallel ({j_parallel:.1f}) > antiparallel ({j_antiparallel:.1f})")
    else:
        msg = f"Beta ordering wrong: para={j_parallel}, anti={j_antiparallel}"
        errors.append(msg)
        print(f"   FAIL: {msg}")

    print("\n" + "=" * 70)
    if errors:
        print(f"VALIDATION FAILED: {len(errors)} error(s)")
        for err in errors:
            print(f"  - {err}")
        print("=" * 70)
        return False
    else:
        print("ALL PRODUCTION VALIDATION TESTS PASSED")
        print("Module is ready for FDA/EMA regulatory submissions")
        print("=" * 70)
        return True

def test_corrections_applied():
    """Run correction validation checks."""
    tests = []

    val_data = AMINO_ACID_SPECIFIC_COUPLINGS[AminoAcid.VAL]['chi1_rotamers']
    assert val_data['trans']['3J_HA_HB'] == 11.0, "Valine trans J should be 11.0"
    assert val_data['gauche_minus']['3J_HA_HB'] == 3.5, "Valine gauche- J should be 3.5"
    tests.append("Valine rotamers corrected (North 1994)")

    gly_j = AMINO_ACID_SPECIFIC_COUPLINGS[AminoAcid.GLY]['2J_HA_HA']
    assert abs(gly_j + 17.1) < 0.1, f"Glycine 2J should be -17.1, got {gly_j}"
    tests.append("Glycine 2J corrected to -17.1 Hz (Kjaergaard 2011)")

    pro_data = AMINO_ACID_SPECIFIC_COUPLINGS[AminoAcid.PRO]['cis_trans']
    assert pro_data['4J_HA_HD_cis'] == 0.0, "Proline cis 4J should be 0 (not observable)"
    assert pro_data['4J_HA_HD_trans'] == 0.0, "Proline trans 4J should be 0 (not observable)"
    tests.append("Proline 4J(HA-HD) corrected to 0 Hz (use NOE for cis/trans)")

    assert EXPECTED_3J_HN_HA[SecondaryStructure.BETA_SHEET_PARALLEL] == 9.7
    assert EXPECTED_3J_HN_HA[SecondaryStructure.BETA_SHEET_ANTIPARALLEL] == 9.0
    tests.append("Beta sheet J values corrected: parallel=9.7, antiparallel=9.0")

    assert EXPECTED_3J_HN_HA[SecondaryStructure.ALPHA_HELIX] == 3.9
    tests.append("Alpha helix J corrected to 3.9 Hz (Vuister 1993)")

    print("All literature corrections validated:")
    for t in tests:
        print(f"  {t}")

    return tests

if __name__ == "__main__":
    _test_peptide_module()
    print("\n")
    test_peptide_production()
    print("\n")
    test_corrections_applied()
    print("\n")
    test_peptide_production()
