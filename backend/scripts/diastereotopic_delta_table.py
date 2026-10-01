# -*- coding: utf-8 -*-
"""Diastereotopic proton chemical-shift difference lookup and estimation utilities."""

from typing import Optional, Tuple, Dict, Any, List
from rdkit import Chem
from rdkit.Chem import rdchem

# Solvent configuration
SOLVENT_ADDITIVE_SHIFT = {

    "CDCl3": 0.0,
    "chloroform": 0.0,
    "chloroform-d": 0.0,
    "CD2Cl2": 0.0,
    "dichloromethane": 0.0,
    "CCl4": 0.0,

    "D2O": -0.02,
    "water": -0.02,
    "MeOD": -0.01,
    "methanol-d4": -0.01,
    "MeOH-d4": -0.01,

    "DMSO-d6": 0.0,
    "DMSO": 0.0,
    "acetone-d6": 0.0,
    "CD3CN": 0.0,
    "DMF-d7": 0.0,
    "THF-d8": 0.0,

    "C6D6": 0.6,
    "benzene-d6": 0.6,
    "benzene": 0.6,
    "toluene-d8": 0.5,
    "pyridine-d5": 0.7,
    "pyridine": 0.7,

    "default": 0.0,
}

SOLVENT_MULTIPLIER = {
    "CDCl3": 1.0, "C6D6": 1.0, "pyridine-d5": 1.0,
    "D2O": 1.0, "MeOD": 1.0, "DMSO-d6": 1.0,
    "default": 1.0,
}

SOLVENT_PROPERTIES = {

    "CDCl3": {"type": "halogenated", "polarity": "low", "protic": False, "aromatic": False},
    "DMSO-d6": {"type": "polar_aprotic", "polarity": "high", "protic": False, "aromatic": False},
    "D2O": {"type": "polar_protic", "polarity": "very_high", "protic": True, "aromatic": False},
    "MeOD": {"type": "polar_protic", "polarity": "high", "protic": True, "aromatic": False},
    "C6D6": {"type": "aromatic", "polarity": "low", "protic": False, "aromatic": True},
    "pyridine-d5": {"type": "aromatic", "polarity": "medium", "protic": False, "aromatic": True},
}

# Diastereotopic delta lookup tables
VINYL_DELTA_BY_SUBSTITUENT = {

    "H": (0.0, 0.0, "éthylène référence"),
    "alkyl_quasi_symmetric": (0.03, 0.05, "C(CH3)2-CH=CH2 achiral, Neoechinulin A"),
    "alkyl": (0.08, 0.12, "alcènes simples, propène ~0.1"),
    "methyl": (0.10, 0.12, "propène"),
    "phenyl": (0.50, 0.55, "styrène: 5.24 vs 5.74 ppm"),
    "aryl": (0.50, 0.55, "aryles en général"),
    "naphthyl": (0.55, 0.60, "naphtalène"),

    "ester": (0.45, 0.55, "acrylate: 0.4-0.5"),
    "COOR": (0.45, 0.55, "esters"),
    "COOR_conjugated": (0.52, 0.55, "MMA: 0.52 ppm confirmé"),
    "acid": (0.45, 0.55, "acide acrylique"),
    "COOH": (0.45, 0.55, "acides"),
    "amide": (0.40, 0.50, "acrylamide"),
    "CONR2": (0.40, 0.50, "amides"),
    "ketone": (0.45, 0.55, "MVK"),
    "COR": (0.45, 0.55, "cétones"),
    "aldehyde": (0.50, 0.60, "acroléine"),
    "CHO": (0.50, 0.60, "aldéhydes"),

    "nitrile": (0.40, 0.45, "acrylonitrile"),
    "CN": (0.40, 0.45, "nitriles"),
    "nitro": (0.60, 0.70, "nitroéthylène - fort effet"),
    "NO2": (0.60, 0.70, "nitro"),

    "F": (0.08, 0.12, "fluor - faible"),
    "Cl": (0.10, 0.15, "chlorure de vinyle"),
    "Br": (0.30, 0.35, "bromure de vinyle"),
    "I": (0.40, 0.45, "iodure de vinyle"),
    "halogen": (0.20, 0.35, "halogènes en général"),

    "ether": (0.15, 0.20, "éther vinylique"),
    "OR": (0.15, 0.20, "alkoxy"),
    "OAr": (0.20, 0.25, "aryloxy"),
    "amine": (0.15, 0.20, "énamine"),
    "NR2": (0.15, 0.20, "amines"),

    "sulfone": (0.50, 0.60, "vinyl sulfone"),
    "SO2R": (0.50, 0.60, "sulfones"),
    "phosphonate": (0.40, 0.50, "vinyl phosphonate"),
    "silyl": (0.15, 0.20, "vinyl silane"),
    "SiR3": (0.15, 0.20, "silyles"),

    "default": (0.30, 0.40, "valeur moyenne vinylique"),
}

BENZYLIC_DELTA = {

    "alkyl": (0.10, 0.15, "Ph-CH₂-CH₃, toluène dérivés"),
    "aryl": (0.06, 0.10, "Ph-CH₂-Ph, diphénylméthane"),
    "alcohol": (0.10, 0.18, "Ph-CH₂-OH benzyl alcool"),
    "amine": (0.08, 0.15, "Ph-CH₂-NH₂ benzylamine"),
    "ether": (0.08, 0.12, "Ph-CH₂-OR benzyl éther"),
    "ester": (0.15, 0.22, "Ph-CH₂-OCOR benzyl ester"),
    "carboxyl": (0.12, 0.18, "Ph-CH₂-COOH acide phénylacétique"),
    "nitrile": (0.18, 0.25, "Ph-CH₂-CN benzyl nitrile - EWG fort"),
    "nitro": (0.20, 0.28, "Ph-CH₂-NO₂ nitro - effet fort"),
    "halogen": (0.12, 0.18, "Ph-CH₂-X halogénure de benzyle"),
    "ionized": (0.03, 0.06, "Ph-CH₂-COO⁻ en D₂O - très faible"),
    "amino_acid": (0.057, 0.08, "Phénylalanine Hβ (BMRB)"),
    "default": (0.10, 0.18, "valeur moyenne benzylique"),
}

ALLYLIC_DELTA = {

    "simple": (0.08, 0.12, "CH₂=CH-CH₂-R simple"),
    "alcohol": (0.10, 0.15, "allyl alcool"),
    "ester": (0.12, 0.18, "allyl ester"),
    "ether": (0.08, 0.12, "allyl éther"),
    "aryl": (0.12, 0.18, "allylbenzène"),
    "terpene": (0.10, 0.15, "géraniol, nérol, terpènes"),
    "steroid": (0.08, 0.12, "stéroïdes insaturés"),
    "cyclic": (0.15, 0.20, "cyclohexène, allylique en cycle"),
    "conjugated": (0.15, 0.22, "butadiène, conjugué"),
    "default": (0.10, 0.15, "valeur moyenne allylique"),
}

EPOXIDE_DELTA = {

    "terminal": (0.30, 0.40, "époxyde terminal R-CH-O-CH₂"),
    "internal": (0.25, 0.35, "époxyde interne R-CH-O-CH-R"),
    "styrene_oxide": (0.28, 0.35, "styrène oxyde"),
    "glycidyl": (0.35, 0.45, "glycidyl éther/ester"),
    "spiro": (0.40, 0.55, "spiro-époxyde - rigide"),
    "bis_epoxy": (0.50, 0.65, "bis-époxyde, ex: 1,3-di(époxypropoxy)benzène"),
    "default": (0.30, 0.40, "valeur moyenne époxyde"),
}

AZIRIDINE_DELTA = {

    "N_H": (0.18, 0.25, "aziridine NH"),
    "N_alkyl": (0.15, 0.22, "N-méthyl aziridine"),
    "N_aryl": (0.20, 0.28, "N-phényl aziridine"),
    "N_acyl": (0.22, 0.30, "N-acyl aziridine"),
    "default": (0.18, 0.25, "valeur moyenne aziridine"),
}

THIIRANE_DELTA = {
    "simple": (0.20, 0.28, "thiirane simple"),
    "substituted": (0.22, 0.30, "thiirane substitué"),
    "default": (0.20, 0.28, "valeur moyenne thiirane"),
}

SULFOXIDE_SULFONE_DELTA = {

    "sulfoxide_alpha": (0.12, 0.18, "CH₂ alpha au S=O"),
    "sulfoxide_beta": (0.06, 0.10, "CH₂ beta au S=O"),
    "methionine": (0.029, 0.05, "Méthionine Hβ (BMRB) - thioéther"),

    "sulfone_alpha": (0.18, 0.25, "CH₂ alpha au SO₂"),
    "sulfone_beta": (0.08, 0.12, "CH₂ beta au SO₂"),

    "sulfonate_ester": (0.15, 0.22, "R-O-SO₂-CH₂"),
    "sulfonamide": (0.12, 0.18, "R-NH-SO₂-CH₂"),

    "default_sulfoxide": (0.12, 0.18, "sulfoxyde défaut"),
    "default_sulfone": (0.18, 0.25, "sulfone défaut"),
}

PHOSPHORUS_DELTA = {

    "phosphine_alkyl": (0.08, 0.12, "R₂P-CH₂ phosphine alkyle"),
    "phosphine_aryl": (0.10, 0.15, "Ar₂P-CH₂ phosphine aryle"),

    "phosphonate": (0.15, 0.20, "(RO)₂P(=O)-CH₂"),
    "phosphonate_dimethyl": (0.15, 0.18, "diméthylphosphonate"),
    "phosphonate_diethyl": (0.14, 0.18, "diéthylphosphonate"),

    "phosphate_ester": (0.12, 0.18, "RO-P(=O)(OR)₂"),

    "phosphine_oxide": (0.18, 0.25, "R₃P(=O)-CH₂"),

    "default": (0.12, 0.18, "valeur moyenne P"),
}

N_OXIDE_DELTA = {

    "pyridine_N_oxide": (0.18, 0.25, "CH₂ adjacent à pyridine N-oxyde"),
    "amine_N_oxide": (0.20, 0.28, "R₃N⁺-O⁻ amine N-oxyde"),
    "azine_N_oxide": (0.15, 0.22, "quinoline N-oxyde, etc."),
    "default": (0.18, 0.25, "valeur moyenne N-oxyde"),
}

SILANE_DELTA = {
    "trialkylsilyl": (0.05, 0.08, "R₃Si-CH₂ TMS, TES, etc."),
    "triarylsilyl": (0.08, 0.12, "Ar₃Si-CH₂"),
    "siloxane": (0.06, 0.10, "(R₃Si)₂O polymères"),
    "chlorosilane": (0.10, 0.15, "Cl₃Si-CH₂ plus électrophile"),
    "default": (0.06, 0.10, "valeur moyenne Si"),
}

BORON_DELTA = {
    "boronate_ester": (0.08, 0.12, "RB(OR')₂ pinacol boronate"),
    "boronic_acid": (0.10, 0.15, "RB(OH)₂"),
    "borane": (0.12, 0.18, "R₃B"),
    "borate": (0.06, 0.10, "B(OR)₃ borate"),
    "default": (0.08, 0.12, "valeur moyenne B"),
}

NUCLEOSIDE_DELTA = {

    "H5_prime": (0.15, 0.20, "H5'/H5\" ribose/désoxyribose"),
    "H5_prime_DNA": (0.12, 0.18, "ADN désoxyribose"),
    "H5_prime_RNA": (0.15, 0.22, "ARN ribose - plus rigide"),
    "C2_endo": (0.12, 0.16, "conformation C2'-endo (ADN-B)"),
    "C3_endo": (0.16, 0.22, "conformation C3'-endo (ARN-A)"),
    "default": (0.14, 0.20, "valeur moyenne nucléoside"),
}

CROWN_ETHER_DELTA = {
    "18_crown_6": (0.05, 0.08, "18-crown-6 symétrique"),
    "15_crown_5": (0.06, 0.10, "15-crown-5"),
    "12_crown_4": (0.08, 0.12, "12-crown-4 plus tendu"),
    "substituted": (0.12, 0.18, "crown substitué asymétrique"),
    "complexed": (0.15, 0.25, "crown complexé avec métal"),
    "cryptand": (0.10, 0.18, "cryptands [2.2.1], [2.2.2]"),
    "default": (0.08, 0.15, "valeur moyenne crown"),
}

CYCLODEXTRIN_DELTA = {
    "alpha_CD": (0.06, 0.10, "α-cyclodextrine (6 unités)"),
    "beta_CD": (0.08, 0.12, "β-cyclodextrine (7 unités)"),
    "gamma_CD": (0.10, 0.15, "γ-cyclodextrine (8 unités)"),
    "inclusion_weak": (0.15, 0.25, "inclusion faible"),
    "inclusion_strong": (0.30, 0.50, "inclusion forte - aromatiques"),
    "default": (0.08, 0.15, "valeur moyenne CD"),
}

ROTAMER_DELTA = {

    "amide_tertiary": (0.08, 0.15, "R₂N-C(=O)-R rotation lente à RT"),
    "DMF_type": (0.05, 0.10, "formamides N,N-disubstitués"),
    "carbamate": (0.10, 0.18, "R₂N-C(=O)-OR"),
    "urea": (0.08, 0.15, "R₂N-C(=O)-NR₂"),
    "thioamide": (0.12, 0.20, "R₂N-C(=S)-R barrière plus haute"),
    "atropisomer": (0.15, 0.30, "biphényles substitués, BINOL"),
    "default": (0.10, 0.18, "valeur moyenne rotamère"),
}

METAL_COMPLEX_DELTA = {

    "ferrocene": (1.5, 2.5, "ferrocène et dérivés"),
    "ruthenocene": (1.2, 2.0, "ruthénocène"),
    "chromocene": (1.0, 1.8, "chromocène"),
    "Pt_II_cis_trans": (3.0, 5.0, "complexes Pt(II) cis/trans"),
    "Pd_II": (2.0, 4.0, "complexes Pd(II)"),
    "Rh_I": (1.5, 3.0, "complexes Rh(I)"),
    "Ti_IV": (0.8, 1.5, "complexes Ti(IV)"),
    "lanthanide": (5.0, 15.0, "lanthanides paramagnétiques - extrême"),
    "default": (1.5, 3.0, "valeur moyenne complexe"),
}

CH3_DIASTEREOTOPIC_DELTA = {

    "isopropyl_chiral_dist1": (0.08, 0.15, "iPr à 1 liaison du C* (ex: 2-butanol 0.12 ppm)"),
    "isopropyl_chiral_dist2": (0.04, 0.08, "iPr à 2 liaisons du C* (ex: menthol)"),
    "isopropyl_chiral_dist3": (0.02, 0.05, "iPr à 3+ liaisons du C*"),

    "gem_dimethyl_chiral": (0.05, 0.10, "C(CH₃)₂ près d'un chiral (pinane, camphre)"),

    "valine_CH3": (0.01, 0.03, "Valine γ-CH₃ BMRB bmse000860: 0.80 vs 0.79 ppm"),
    "leucine_CH3": (0.04, 0.07, "Leucine δ-CH₃ BMRB bmse000052: médiane 0.04 ppm"),
    "isoleucine_CH3": (0.07, 0.10, "Isoleucine γ-CH₃ vs δ-CH₃ (séc-butyle, pas iPr)"),

    "terpene_gem": (0.04, 0.08, "géraniol, menthol, etc."),
    "steroid_angular_CH3": (0.02, 0.05, "méthyles angulaires stéroïdes"),

    "N_methyl_chiral": (0.03, 0.08, "N-CH₃ sur amine chirale stable (N-oxydes)"),
    "S_methyl_chiral": (0.05, 0.10, "S-CH₃ sur sulfonium chiral"),

    "default": (0.05, 0.10, "valeur moyenne CH₃ diastéréotopique"),
}

CH3_CONFORMATIONAL_DELTA = {
    "cyclohexyl_CH3_ax_eq": (0.10, 0.20, "CH₃ axial vs équatorial sur cyclohexane"),
    "menthol_type": (0.03, 0.06, "Menthol: rotation rapide atténue"),
    "steroid_ring_CH3": (0.05, 0.12, "CH₃ sur cycle stéroïde rigide"),
    "default": (0.08, 0.15, "valeur moyenne conformationnelle"),
}

CH3_GEOMETRIC_EZ_DELTA = {

    "vinyl_cis_trans": (0.15, 0.35, "CH₃ cis vs trans sur C=C"),
    "isoprene_type": (0.10, 0.15, "Isoprène: 1.78 vs 1.68 ppm"),
    "beta_ionone_type": (0.12, 0.18, "β-ionone: 1.90 vs 2.02 ppm"),
    "myrcene_type": (0.08, 0.12, "Myrcène et terpènes linéaires"),
    "default": (0.12, 0.25, "valeur moyenne géométrie E/Z"),
}

CH3_STERIC_DELTA = {
    "hindered_rotation": (0.03, 0.08, "rotation restreinte par encombrement"),
    "tert_butyl_compressed": (0.02, 0.05, "t-Bu comprimé (rare)"),
    "default": (0.03, 0.06, "valeur moyenne stérique"),
}

VINYL_CH3_GEOMETRY = {

    "cis_to_H": (1.60, "CH₃ cis au H vinylique"),
    "trans_to_H": (1.75, "CH₃ trans au H vinylique"),
    "delta_cis_trans": (0.15, "Différence typique entre cis et trans"),

    "isoprene_CH3": {"cis": 1.78, "trans": 1.68, "delta": 0.10},
    "beta_ionone_CH3": {"cis": 1.90, "trans": 2.02, "delta": 0.12},
    "myrcene_CH3": {"cis": 1.68, "trans": 1.60, "delta": 0.08},
    "geraniol_CH3": {"cis": 1.68, "trans": 1.60, "delta": 0.08},
    "farnesol_CH3": {"cis": 1.68, "trans": 1.60, "delta": 0.08},
}

CH2_CHIRAL_DELTA_BY_DISTANCE = {

    1: (0.08, 0.25, "1 liaison: médiane 0.08 ppm (BMRB); max 0.25 sauf cas spéciaux"),
    2: (0.06, 0.15, "2 liaisons: effet diminue rapidement"),
    3: (0.04, 0.10, "3 liaisons: effet quasi-nul pour beaucoup de systèmes"),
    4: (0.03, 0.08, "4 liaisons: très faible"),
    5: (0.02, 0.05, "5 liaisons: négligeable"),
    6: (0.01, 0.03, "6+ liaisons: souvent non résolu"),
}

CH2_CHIRAL_SPECIAL_CASES = {
    "citric_acid": 0.14,
    "malic_acid": 0.15,
    "2_butanol": 0.15,
    "lactic_acid": 0.12,
    "mandelic_acid": 0.18,
    "epoxypropoxy": 0.60,
}

CYCLIC_DELTA = {

    "cyclohexane": (0.48, 0.50, "CS₂ -103°C, δax=1.12, δeq=1.60"),
    "cyclohexane_RT": (0.0, 0.0, "moyennage rapide à RT, signal à 1.44 ppm"),
    "cyclohexanone": (0.24, 0.30, "C4-H, <-185°C"),
    "cyclohexanol": (0.45, 0.50, "similaire cyclohexane"),
    "methylcyclohexane": (0.48, 0.52, "légère variation"),

    "tetrahydropyran": (0.40, 0.56, "moyenne THP"),
    "THP": (0.40, 0.56, "alias"),
    "THP_H2": (0.558, 0.60, "H2: 3.439 vs 3.997 ppm"),
    "THP_H3": (0.116, 0.15, "EXCEPTION: δax > δeq, effet γ O"),
    "THP_H4": (0.367, 0.40, "H4: 1.493 vs 1.860 ppm"),

    "dioxane_1_3": (0.50, 0.72, "1,3-dioxane"),
    "dioxane_H4": (0.29, 0.35, "H4: 3.61 vs 3.90"),
    "dioxane_H5": (0.63, 0.72, "H5: 1.15 vs 1.78 ppm - fort!"),
    "2_methyl_dioxane_H5": (0.72, 0.75, "2-méthyl-1,3-dioxane: 1.23 vs 1.95 ppm"),
    "4_methyl_dioxane_H5": (0.75, 0.80, "4-méthyl augmente encore Δδ"),
    "dioxane_1_4": (0.0, 0.0, "1,4-dioxane symétrique - équivalent"),

    "dioxolane": (0.30, 0.40, "5 membres"),

    "pyranose": (0.18, 0.25, "H6a/H6b typique"),
    "glucopyranose": (0.15, 0.21, "β-D-Glucopyranose"),
    "galactopyranose": (0.18, 0.25, "galactose"),
    "mannopyranose": (0.16, 0.22, "mannose"),
    "methyl_glucoside": (0.21, 0.25, "Méthyl 6-O-trityl-α-D-glucopyranoside"),

    "furanose": (0.12, 0.20, "5 membres"),
    "ribofuranose": (0.10, 0.18, "ribose"),

    "piperidine": (0.45, 0.55, "N-alkyl piperidine"),
    "piperazine": (0.40, 0.50, "piperazine"),
    "morpholine": (0.40, 0.50, "morpholine"),
    "pyrrolidine": (0.25, 0.35, "5 membres"),
    "proline": (0.35, 0.45, "proline cyclique"),

    "lactone_6": (0.35, 0.45, "δ-lactone"),
    "lactone_5": (0.25, 0.35, "γ-lactone"),
    "lactam_6": (0.30, 0.40, "δ-lactame"),
    "lactam_5": (0.25, 0.35, "γ-lactame"),
    "caprolactam": (0.35, 0.45, "7 membres"),

    "norbornane": (0.25, 0.35, "exo/endo"),
    "7_oxanorbornane": (0.25, 0.30, "7-Oxanorbornane"),
    "azabicycle_3_3_1": (1.30, 1.40, "[3.3.1]Azabicycle - H endo compressé par lone pair N ~1.3 ppm"),
    "azabicycle_2_2_1": (0.30, 0.40, "Azanorbornane exo/endo"),
    "decalin_trans": (0.45, 0.55, "trans-décaline rigide"),
    "decalin_cis": (0.40, 0.50, "cis-décaline - légèrement moins que trans"),

    "steroid_A_ring": (0.50, 0.60, "cycle A stéroïdes - Δδ ax/eq ~0.5-0.6 ppm"),
    "steroid_B_ring": (0.45, 0.55, "cycle B cis-décaline fusionnée"),
    "steroid_C_ring": (0.40, 0.50, "cycle C"),
    "steroid_D_ring": (0.35, 0.45, "cycle D cyclopentane ~0.3-0.4 ppm"),

    "default_3ring": (0.15, 0.25, "cyclopropane"),
    "default_4ring": (0.20, 0.30, "cyclobutane"),
    "default_5ring": (0.25, 0.35, "5 membres général"),
    "default_6ring": (0.45, 0.55, "6 membres général"),
    "default_7ring": (0.40, 0.50, "7 membres"),
    "default_8ring": (0.35, 0.45, "8+ membres"),
}

AMINO_ACID_BETA_DELTA = {

    "ILE": (0.209, 0.02, "Isoleucine Hγ1/Hγ2: 1.458 vs 1.249 ppm (BMRB bmse000866, D₂O 298K)"),
    "CYS": (0.075, 0.06, "Cystéine - le plus grand Hβ"),
    "PRO": (0.072, 0.004, "Proline"),
    "GLY": (0.068, 0.004, "Glycine Hα1/Hα2"),
    "TRP": (0.060, 0.004, "Tryptophane Hβ"),
    "LEU": (0.059, 0.004, "Leucine"),
    "HIS": (0.058, 0.010, "Histidine"),
    "PHE": (0.057, 0.004, "Phénylalanine"),
    "ASN": (0.056, 0.004, "Asparagine"),
    "TYR": (0.055, 0.004, "Tyrosine"),
    "ASP": (0.049, 0.005, "Aspartate"),
    "PRO_Hdelta": (0.034, 0.003, "Proline Hδ"),
    "ARG": (0.029, 0.003, "Arginine"),
    "MET": (0.029, 0.005, "Méthionine"),
    "GLN": (0.028, 0.003, "Glutamine"),
    "LYS": (0.027, 0.003, "Lysine"),
    "SER": (0.024, 0.003, "Sérine"),
    "GLU": (0.023, 0.002, "Glutamate"),
    "THR": (0.020, 0.002, "Thréonine"),
    "VAL": (0.015, 0.002, "Valine"),
    "ALA": (0.0, 0.0, "Alanine - pas de CH₂"),
    "default": (0.04, 0.01, "valeur moyenne AA"),
}

SECONDARY_STRUCTURE_EFFECT = {

    "helix_alpha": {
        "Ha_shift": -0.39,
        "Hb_shift": -0.10,
        "delta_increase": 0.20,
        "note": "Hα décalé de -0.39 ppm en hélice α (Wishart et al.)",
    },
    "sheet_beta": {
        "Ha_shift": +0.37,
        "Hb_shift": +0.10,
        "delta_increase": 0.25,
        "note": "Hα décalé de +0.37 ppm en feuillet β (Wishart et al.)",
    },
    "turn": {
        "Ha_shift": -0.20,
        "Hb_shift": 0.0,
        "delta_increase": 0.15,
        "note": "Coudes et boucles",
    },
    "random_coil": {
        "Ha_shift": 0.0,
        "Hb_shift": 0.0,
        "delta_increase": 0.0,
        "note": "Référence - peptides non structurés",
    },

    "trp_cage_gly": {
        "Ha_shift": 0.0,
        "Hb_shift": 0.0,
        "delta_increase": 2.30,
        "note": "Trp-cage: Gly rigidifié + ring current aromatique = Δδ extrême",
    },
}

SECONDARY_STRUCTURE_SHIFT = {
    "helix_alpha": (-0.39, -0.10),
    "sheet_beta": (+0.37, +0.10),
    "random_coil": (0.0, 0.0),
}

ANISOTROPIC_SYSTEMS = {

    "calixarene_cone": (1.43, 1.70, "ArCH₂Ar cone: Ha 4.50, He 3.07 ppm (J.Org.Chem.1985)"),
    "calixarene_partial_cone": (0.80, 1.20, "partial cone"),
    "calixarene_1_2_alternate": (0.50, 0.80, "1,2-alternate"),
    "calixarene_1_3_alternate": (0.30, 0.50, "1,3-alternate - plus petit"),

    "paracyclophane_6": (3.00, 3.50, "[6]paracyclophane: endo -0.5, exo 2.5 ppm - Δδ~3.0"),
    "paracyclophane_8": (0.80, 1.00, "[8]paracyclophane - plus flexible"),
    "metacyclophane": (0.60, 0.80, "métacyclophane ~0.6-0.7 ppm"),

    "porphyrin_NH": (11.2, 13.0, "NH internes: -2.8 ppm vs méso 8.4 ppm (JACS 1973)"),
    "porphyrin_meso": (2.0, 3.0, "protons méso"),

    "annulene_18_neutral": (12.1, 12.5, "[18]annulène neutre: int +9.2, ext -2.9 ppm"),

    "annulene_18_dianion": (31.0, 33.0, "[18]annulène²⁻: int +29.5-32, ext -1.13 (JACS 1972) - RECORD"),
    "annulene_18_tetraanion": (16.0, 17.0, "[18]annulène⁴⁻ aromatique"),
    "annulene_14": (8.0, 10.0, "[14]annulène"),

    "fullerene_adduct": (1.0, 2.0, "adduits fullerène"),
    "cryptand": (0.50, 0.80, "cryptands"),

    "carbonyl_alpha": (0.30, 0.50, "α au C=O, dans le plan"),
    "carbonyl_beta": (0.15, 0.25, "β au C=O"),

    "alkyne_propargyl": (0.20, 0.30, "CH₂ propargylique"),
}

ANISOTROPIC_AMPLIFICATION = {

    "benzene_face": 2.5,
    "benzene_edge": 1.3,
    "carbonyl_plane": 1.8,
    "carbonyl_perpendicular": 0.9,
    "triple_bond_axis": 0.7,
    "triple_bond_perp": 1.2,
    "nitro_group": 1.5,
    "calixarene_cavity": 5.0,
    "porphyrin_ring": 10.0,
    "intramol_hbond": 1.4,
}

ADDITIONAL_DELTA_FACTORS = {

    "intramol_hydrogen_bond": (0.25, "Liaison H intramoléculaire figeant une conformation"),
    "steric_compression": (0.20, "Compression stérique (ex: norbornane endo)"),
    "ring_current_face": (0.50, "Proton au-dessus d'un cycle aromatique"),
    "lone_pair_compression": (0.40, "Compression par paire libre N/O (ex: azabicycles)"),
}

DEFAULT_DELTA = {

    "vinyl_geminal": 0.30,
    "vinyl_phenyl": 0.50,
    "vinyl_ester": 0.45,
    "vinyl_symmetric": 0.03,

    "ch2_near_chiral_1": 0.08,
    "ch2_near_chiral_2": 0.06,
    "ch2_near_chiral_3": 0.04,
    "ch2_chiral_median": 0.07,

    "cyclic_ax_eq": 0.45,
    "cyclic_ax_eq_rigid": 0.50,
    "cyclic_ax_eq_RT": 0.0,
    "cyclic_THP": 0.40,
    "cyclic_pyranose": 0.18,

    "amino_acid": 0.04,
    "amino_acid_helix": 0.15,
    "amino_acid_sheet": 0.18,

    "flexible_chain": 0.08,
    "rigid_system": 0.30,

    "calixarene": 1.43,

    "significant_threshold": 0.10,
    "unknown": 0.07,
}

# Detection helpers
def _is_vinyl_carbon(mol: Chem.Mol, c_idx: int) -> bool:
    """Check whether the target carbon is a terminal vinyl carbon."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetSymbol() != 'C':
        return False

    if atom.GetHybridization() != rdchem.HybridizationType.SP2:
        return False

    has_double = False
    for bond in atom.GetBonds():
        if bond.GetBondType() == rdchem.BondType.DOUBLE:
            other = bond.GetOtherAtom(atom)
            if other.GetSymbol() == 'C':
                has_double = True
                break

    if not has_double:
        return False

    h_count = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == 'H')
    return h_count == 2

def _is_gem_dimethyl_achiral(mol: Chem.Mol, c_idx: int) -> bool:
    """Check whether the target carbon belongs to an achiral gem-dimethyl group."""
    atom = mol.GetAtomWithIdx(c_idx)

    if atom.GetSymbol() != 'C':
        return False

    carbon_neighbors = []
    methyl_neighbors = []

    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() == 'C':
            carbon_neighbors.append(neighbor)

            h_count = sum(1 for n in neighbor.GetNeighbors() if n.GetSymbol() == 'H')
            if h_count == 3:
                methyl_neighbors.append(neighbor.GetIdx())

    if len(methyl_neighbors) != 2:
        return False

    if len(carbon_neighbors) in (3, 4) and len(methyl_neighbors) == 2:
        return True

    return False

def _get_vinyl_substituent_type(mol: Chem.Mol, vinyl_c_idx: int) -> str:
    """Identify the substituent type adjacent to a vinyl carbon."""
    atom = mol.GetAtomWithIdx(vinyl_c_idx)

    adjacent_c = None
    for bond in atom.GetBonds():
        if bond.GetBondType() == rdchem.BondType.DOUBLE:
            other = bond.GetOtherAtom(atom)
            if other.GetSymbol() == 'C':
                adjacent_c = other
                break

    if adjacent_c is None:
        return "default"

    for neighbor in adjacent_c.GetNeighbors():
        if neighbor.GetIdx() == vinyl_c_idx:
            continue

        sym = neighbor.GetSymbol()

        if neighbor.GetIsAromatic():
            return "phenyl"

        if sym == 'C':
            for n2 in neighbor.GetNeighbors():
                if n2.GetSymbol() == 'O':
                    for b in mol.GetBondBetweenAtoms(neighbor.GetIdx(), n2.GetIdx()):
                        pass

                    bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), n2.GetIdx())
                    if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                        return "ester"

        if sym == 'C':
            for n2 in neighbor.GetNeighbors():
                if n2.GetSymbol() == 'N':
                    bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), n2.GetIdx())
                    if bond and bond.GetBondType() == rdchem.BondType.TRIPLE:
                        return "nitrile"

        if sym in ('Cl', 'Br', 'I', 'F'):
            return sym if sym in ('Cl', 'Br') else "halogen"

        if sym == 'N':
            o_count = sum(1 for n2 in neighbor.GetNeighbors() if n2.GetSymbol() == 'O')
            if o_count >= 2:
                return "nitro"

        if sym == 'O':
            return "ether"

        if sym == 'C' and not neighbor.GetIsAromatic():

            if _is_gem_dimethyl_achiral(mol, neighbor.GetIdx()):
                return "alkyl_quasi_symmetric"
            return "alkyl"

    return "default"

def _is_benzylic_carbon(mol: Chem.Mol, c_idx: int) -> bool:
    """Check whether the target carbon is benzylic."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetSymbol() != 'C':
        return False

    if atom.GetHybridization() != rdchem.HybridizationType.SP3:
        return False

    h_count = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == 'H')
    if h_count != 2:
        return False

    for neighbor in atom.GetNeighbors():
        if neighbor.GetIsAromatic():
            return True

    return False

def _get_benzylic_substituent_type(mol: Chem.Mol, c_idx: int) -> str:
    """Identify the substituent type on a benzylic carbon."""
    atom = mol.GetAtomWithIdx(c_idx)

    for neighbor in atom.GetNeighbors():
        if neighbor.GetIsAromatic():
            continue
        if neighbor.GetSymbol() == 'H':
            continue

        sym = neighbor.GetSymbol()

        if sym == 'O':

            for n2 in neighbor.GetNeighbors():
                if n2.GetSymbol() == 'C' and n2.GetIdx() != c_idx:

                    for n3 in n2.GetNeighbors():
                        if n3.GetSymbol() == 'O':
                            bond = mol.GetBondBetweenAtoms(n2.GetIdx(), n3.GetIdx())
                            if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                                return "ester"
                    return "ether"
            return "alcohol"

        if sym == 'N':
            return "amine"

        if sym == 'C':
            for n2 in neighbor.GetNeighbors():
                if n2.GetSymbol() == 'O':
                    bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), n2.GetIdx())
                    if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                        return "carboxyl"
                if n2.GetSymbol() == 'N':
                    bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), n2.GetIdx())
                    if bond and bond.GetBondType() == rdchem.BondType.TRIPLE:
                        return "nitrile"
            if neighbor.GetIsAromatic():
                return "aryl"
            return "alkyl"

        if sym in ('Cl', 'Br', 'I', 'F'):
            return "halogen"

        if sym == 'S':
            return "thioether"

    return "default"

def _is_allylic_carbon(mol: Chem.Mol, c_idx: int) -> bool:
    """Check whether the target carbon is allylic."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetSymbol() != 'C':
        return False

    if atom.GetHybridization() != rdchem.HybridizationType.SP3:
        return False

    h_count = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == 'H')
    if h_count != 2:
        return False

    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() == 'C' and neighbor.GetHybridization() == rdchem.HybridizationType.SP2:
            for bond in neighbor.GetBonds():
                if bond.GetBondType() == rdchem.BondType.DOUBLE:
                    other = bond.GetOtherAtom(neighbor)
                    if other.GetSymbol() == 'C':
                        return True

    return False

def _is_epoxide_carbon(mol: Chem.Mol, c_idx: int) -> bool:
    """Check whether the target carbon belongs to an epoxide."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetSymbol() != 'C':
        return False

    if not atom.IsInRing():
        return False

    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if c_idx in ring and len(ring) == 3:

            for idx in ring:
                if mol.GetAtomWithIdx(idx).GetSymbol() == 'O':
                    return True

    return False

def _is_aziridine_carbon(mol: Chem.Mol, c_idx: int) -> bool:
    """Check whether the target carbon belongs to an aziridine."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetSymbol() != 'C':
        return False

    if not atom.IsInRing():
        return False

    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if c_idx in ring and len(ring) == 3:
            for idx in ring:
                if mol.GetAtomWithIdx(idx).GetSymbol() == 'N':
                    return True

    return False

def _is_thiirane_carbon(mol: Chem.Mol, c_idx: int) -> bool:
    """Check whether the target carbon belongs to a thiirane."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetSymbol() != 'C':
        return False

    if not atom.IsInRing():
        return False

    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if c_idx in ring and len(ring) == 3:
            for idx in ring:
                if mol.GetAtomWithIdx(idx).GetSymbol() == 'S':
                    return True

    return False

def _is_sulfoxide_sulfone_adjacent(mol: Chem.Mol, c_idx: int) -> Tuple[bool, str]:
    """Detect adjacency to a sulfoxide or sulfone group."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetSymbol() != 'C':
        return False, ""

    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() == 'S':

            o_count = 0
            for n2 in neighbor.GetNeighbors():
                if n2.GetSymbol() == 'O':
                    o_count += 1

            if o_count == 2:
                return True, "sulfone_alpha"
            elif o_count == 1:
                return True, "sulfoxide_alpha"
            else:
                return True, "thioether"

    return False, ""

def _is_phosphorus_adjacent(mol: Chem.Mol, c_idx: int) -> Tuple[bool, str]:
    """Detect adjacency to a phosphorus-containing group."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetSymbol() != 'C':
        return False, ""

    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() == 'P':

            o_count = 0
            o_double = False
            for n2 in neighbor.GetNeighbors():
                if n2.GetSymbol() == 'O':
                    o_count += 1
                    bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), n2.GetIdx())
                    if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                        o_double = True

            if o_double and o_count >= 2:
                return True, "phosphonate"
            elif o_double:
                return True, "phosphine_oxide"
            else:
                return True, "phosphine_alkyl"

    return False, ""

def _is_n_oxide_adjacent(mol: Chem.Mol, c_idx: int) -> bool:
    """Check whether the target carbon is adjacent to an N-oxide."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetSymbol() != 'C':
        return False

    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() == 'N':

            if neighbor.GetFormalCharge() == 1:
                for n2 in neighbor.GetNeighbors():
                    if n2.GetSymbol() == 'O' and n2.GetFormalCharge() == -1:
                        return True

            for n2 in neighbor.GetNeighbors():
                if n2.GetSymbol() == 'O':

                    if len([x for x in n2.GetNeighbors() if x.GetSymbol() != 'N']) == 0:
                        return True

    return False

def _is_silane_adjacent(mol: Chem.Mol, c_idx: int) -> bool:
    """Check whether the target carbon is adjacent to silicon."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetSymbol() != 'C':
        return False

    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() == 'Si':
            return True

    return False

def _is_boron_adjacent(mol: Chem.Mol, c_idx: int) -> bool:
    """Check whether the target carbon is adjacent to boron."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetSymbol() != 'C':
        return False

    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() == 'B':
            return True

    return False

def _is_nucleoside(mol: Chem.Mol, c_idx: int) -> bool:
    """Check whether the target carbon belongs to a nucleoside motif."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetSymbol() != 'C':
        return False

    o_neighbor = None
    c_ring_neighbor = None

    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() == 'O' and not neighbor.IsInRing():
            o_neighbor = neighbor
        if neighbor.GetSymbol() == 'C' and neighbor.IsInRing():
            c_ring_neighbor = neighbor

    if o_neighbor and c_ring_neighbor:

        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if c_ring_neighbor.GetIdx() in ring and len(ring) == 5:
                has_ring_o = any(mol.GetAtomWithIdx(i).GetSymbol() == 'O' for i in ring)
                if has_ring_o:
                    return True

    return False

def _is_crown_ether_carbon(mol: Chem.Mol, c_idx: int) -> bool:
    """Check whether the target carbon belongs to a crown ether."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetSymbol() != 'C':
        return False

    if not atom.IsInRing():
        return False

    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if c_idx in ring and len(ring) >= 12:
            o_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == 'O')
            if o_count >= 3:
                return True

    return False

def _is_rotamer_blocked(mol: Chem.Mol, c_idx: int) -> Tuple[bool, str]:
    """Detect restricted rotation around the target carbon environment."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetSymbol() != 'C':
        return False, ""

    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() == 'N':

            for n2 in neighbor.GetNeighbors():
                if n2.GetSymbol() == 'C':

                    for n3 in n2.GetNeighbors():
                        if n3.GetSymbol() == 'O':
                            bond = mol.GetBondBetweenAtoms(n2.GetIdx(), n3.GetIdx())
                            if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:

                                n_h_count = sum(1 for x in neighbor.GetNeighbors() if x.GetSymbol() == 'H')
                                if n_h_count == 0:

                                    for n4 in n2.GetNeighbors():
                                        if n4.GetSymbol() == 'O' and n4.GetIdx() != n3.GetIdx():
                                            return True, "carbamate"
                                    return True, "amide_tertiary"

    return False, ""

def _is_geminal_dimethyl(mol: Chem.Mol, c1_idx: int, c2_idx: int) -> bool:
    """Check whether two carbons are geminal methyl groups."""
    atom1 = mol.GetAtomWithIdx(c1_idx)
    atom2 = mol.GetAtomWithIdx(c2_idx)

    if atom1.GetSymbol() != 'C' or atom2.GetSymbol() != 'C':
        return False

    h1_count = sum(1 for n in atom1.GetNeighbors() if n.GetSymbol() == 'H')
    h2_count = sum(1 for n in atom2.GetNeighbors() if n.GetSymbol() == 'H')

    if h1_count != 3 or h2_count != 3:
        return False

    c1_neighbors = [n.GetIdx() for n in atom1.GetNeighbors() if n.GetSymbol() == 'C']
    c2_neighbors = [n.GetIdx() for n in atom2.GetNeighbors() if n.GetSymbol() == 'C']

    if len(c1_neighbors) != 1 or len(c2_neighbors) != 1:
        return False

    return c1_neighbors[0] == c2_neighbors[0]

def _find_chiral_distance_for_ch3(mol: Chem.Mol, ch3_carbon_idx: int) -> int:
    """Return the bond distance from a methyl carbon to the nearest chiral center."""
    try:
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        chiral_indices = {idx for idx, _ in chiral_centers}

        if not chiral_indices:
            return -1

        visited = {ch3_carbon_idx}
        queue = [(ch3_carbon_idx, 0)]

        while queue:
            current_idx, dist = queue.pop(0)

            if dist > 6:
                break

            if current_idx in chiral_indices:
                return dist

            current = mol.GetAtomWithIdx(current_idx)
            for nb in current.GetNeighbors():
                if nb.GetIdx() not in visited and nb.GetSymbol() in ('C', 'N', 'O', 'S'):
                    visited.add(nb.GetIdx())
                    queue.append((nb.GetIdx(), dist + 1))

        return -1

    except Exception:
        return -1

def get_ch3_diastereotopic_delta(
    mol: Chem.Mol,
    ch3_c1_idx: int,
    ch3_c2_idx: int,
    solvent: str = "CDCl3"
) -> Tuple[float, str]:
    """Estimate the chemical-shift difference between geminal methyl groups."""
    try:

        if not _is_geminal_dimethyl(mol, ch3_c1_idx, ch3_c2_idx):
            return 0.0, "not_geminal_dimethyl"

        solvent_shift = get_solvent_additive_shift(solvent)
        solvent_note = f"_solvent_{solvent_shift:+.2f}ppm" if solvent_shift != 0.0 else ""

        dist1 = _find_chiral_distance_for_ch3(mol, ch3_c1_idx)
        dist2 = _find_chiral_distance_for_ch3(mol, ch3_c2_idx)
        chiral_dist = min(dist1, dist2) if dist1 > 0 and dist2 > 0 else max(dist1, dist2)

        if chiral_dist <= 0:

            return 0.0, "no_chiral_center"

        is_aa, aa_code = _is_amino_acid(mol, ch3_c1_idx)
        if is_aa:
            if aa_code == "VAL":
                base_delta = CH3_DIASTEREOTOPIC_DELTA["valine_CH3"][0]
                return base_delta + solvent_shift, f"amino_acid_VAL_CH3{solvent_note}"
            elif aa_code == "LEU":
                base_delta = CH3_DIASTEREOTOPIC_DELTA["leucine_CH3"][0]
                return base_delta + solvent_shift, f"amino_acid_LEU_CH3{solvent_note}"
            elif aa_code == "ILE":
                base_delta = CH3_DIASTEREOTOPIC_DELTA["isoleucine_CH3"][0]
                return base_delta + solvent_shift, f"amino_acid_ILE_CH3{solvent_note}"

        if chiral_dist == 1:
            base_delta = CH3_DIASTEREOTOPIC_DELTA["isopropyl_chiral_dist1"][0]
        elif chiral_dist == 2:
            base_delta = CH3_DIASTEREOTOPIC_DELTA["isopropyl_chiral_dist2"][0]
        else:
            base_delta = CH3_DIASTEREOTOPIC_DELTA["isopropyl_chiral_dist3"][0]

        final_delta = base_delta + solvent_shift
        return final_delta, f"isopropyl_chiral_dist{chiral_dist}{solvent_note}"

    except Exception as e:
        return CH3_DIASTEREOTOPIC_DELTA["default"][0], f"error_{str(e)[:20]}"

def _is_in_ring(mol: Chem.Mol, atom_idx: int) -> bool:
    """Check whether an atom belongs to a ring."""
    return mol.GetAtomWithIdx(atom_idx).IsInRing()

def _get_ring_type(mol: Chem.Mol, c_idx: int) -> str:
    """Identify the ring type containing the target carbon."""
    atom = mol.GetAtomWithIdx(c_idx)

    if not atom.IsInRing():
        return "none"

    ring_info = mol.GetRingInfo()

    for ring in ring_info.AtomRings():
        if c_idx in ring:
            ring_size = len(ring)

            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            o_count = sum(1 for a in ring_atoms if a.GetSymbol() == 'O')
            n_count = sum(1 for a in ring_atoms if a.GetSymbol() == 'N')
            s_count = sum(1 for a in ring_atoms if a.GetSymbol() == 'S')

            if ring_size == 6:
                if o_count == 2:
                    return "dioxane"
                elif o_count == 1 and n_count == 0:
                    return "tetrahydropyran"
                elif n_count == 1 and o_count == 1:
                    return "morpholine"
                elif n_count == 1:
                    return "piperidine"
                elif n_count == 2:
                    return "piperazine"
                elif o_count == 0 and n_count == 0:

                    oh_count = 0
                    for a in ring_atoms:
                        for nb in a.GetNeighbors():
                            if nb.GetSymbol() == 'O' and nb.GetIdx() not in ring:
                                oh_count += 1
                    if oh_count >= 2:
                        return "pyranose"
                    return "cyclohexane"
                else:
                    return "default_6ring"

            elif ring_size == 5:
                if o_count == 2:
                    return "dioxolane"
                elif o_count == 1:

                    return "furanose"
                elif n_count >= 1:
                    return "pyrrolidine"
                else:
                    return "default_5ring"

            else:
                return f"ring_{ring_size}"

    return "unknown_ring"

def _find_chiral_distance(mol: Chem.Mol, ch2_carbon_idx: int) -> int:
    """Return the bond distance to the nearest chiral center."""
    try:

        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        chiral_indices = {idx for idx, _ in chiral_centers}

        if not chiral_indices:
            return -1

        visited = {ch2_carbon_idx}
        queue = [(ch2_carbon_idx, 0)]

        while queue:
            current_idx, dist = queue.pop(0)

            if dist > 6:
                break

            if current_idx in chiral_indices and current_idx != ch2_carbon_idx:
                return dist

            current = mol.GetAtomWithIdx(current_idx)
            for nb in current.GetNeighbors():
                if nb.GetIdx() not in visited and nb.GetSymbol() in ('C', 'N', 'O', 'S'):
                    visited.add(nb.GetIdx())
                    queue.append((nb.GetIdx(), dist + 1))

        return -1

    except Exception:
        return -1

def _has_adjacent_anisotropic_group(mol: Chem.Mol, c_idx: int) -> Tuple[bool, float]:
    """Detect an adjacent anisotropic group and its amplification factor."""
    atom = mol.GetAtomWithIdx(c_idx)

    for neighbor in atom.GetNeighbors():

        if neighbor.GetIsAromatic():
            return True, ANISOTROPIC_AMPLIFICATION.get("benzene_face", 2.0)

        if neighbor.GetSymbol() == 'C':
            for n2 in neighbor.GetNeighbors():
                if n2.GetSymbol() == 'O':
                    bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), n2.GetIdx())
                    if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                        return True, ANISOTROPIC_AMPLIFICATION.get("carbonyl_plane", 1.5)

    return False, 1.0

def _is_amino_acid(mol: Chem.Mol, c_idx: int) -> Tuple[bool, str]:
    """Detect supported amino-acid environments."""
    try:
        atom = mol.GetAtomWithIdx(c_idx)
        if atom.GetSymbol() != 'C':
            return False, ""

        ca_candidate = None
        nitrogen_atom = None
        carbonyl_carbon = None

        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                potential_ca = neighbor
                n_found = None
                co_found = None

                for n2 in potential_ca.GetNeighbors():
                    if n2.GetSymbol() == 'N':
                        n_found = n2
                    if n2.GetSymbol() == 'C' and n2.GetIdx() != c_idx:

                        for n3 in n2.GetNeighbors():
                            if n3.GetSymbol() == 'O':
                                bond = mol.GetBondBetweenAtoms(n2.GetIdx(), n3.GetIdx())
                                if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                                    co_found = n2
                                    break

                if n_found and co_found:
                    ca_candidate = potential_ca
                    nitrogen_atom = n_found
                    carbonyl_carbon = co_found
                    break

        if not ca_candidate:

            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    cb_candidate = neighbor
                    for n2 in cb_candidate.GetNeighbors():
                        if n2.GetSymbol() == 'C' and n2.GetIdx() != c_idx:
                            potential_ca = n2
                            n_found = None
                            co_found = None

                            for n3 in potential_ca.GetNeighbors():
                                if n3.GetSymbol() == 'N':
                                    n_found = n3
                                if n3.GetSymbol() == 'C' and n3.GetIdx() != cb_candidate.GetIdx():
                                    for n4 in n3.GetNeighbors():
                                        if n4.GetSymbol() == 'O':
                                            bond = mol.GetBondBetweenAtoms(n3.GetIdx(), n4.GetIdx())
                                            if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                                                co_found = n3
                                                break

                            if n_found and co_found:
                                ca_candidate = potential_ca
                                nitrogen_atom = n_found
                                carbonyl_carbon = co_found
                                break
                    if ca_candidate:
                        break

        if not ca_candidate:
            return False, ""

        if ca_candidate.IsInRing():

            ring_info = mol.GetRingInfo()
            ca_rings = [r for r in ring_info.AtomRings() if ca_candidate.GetIdx() in r]

            for ring in ca_rings:
                ring_size = len(ring)

                if ring_size == 6:
                    n_count = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == 'N')
                    co_count = 0
                    for i in ring:
                        at = mol.GetAtomWithIdx(i)
                        if at.GetSymbol() == 'C':
                            for nb in at.GetNeighbors():
                                if nb.GetSymbol() == 'O' and nb.GetIdx() not in ring:
                                    bond = mol.GetBondBetweenAtoms(i, nb.GetIdx())
                                    if bond and bond.GetBondTypeAsDouble() >= 2:
                                        co_count += 1
                    if n_count >= 2 and co_count >= 2:
                        return False, ""

                if nitrogen_atom.GetIdx() in ring or carbonyl_carbon.GetIdx() in ring:
                    return False, ""

        if nitrogen_atom.IsInRing():
            ring_info = mol.GetRingInfo()
            n_rings = [r for r in ring_info.AtomRings() if nitrogen_atom.GetIdx() in r]
            for ring in n_rings:

                if carbonyl_carbon.GetIdx() in ring:
                    return False, ""

        for nb in ca_candidate.GetNeighbors():
            if nb.GetIsAromatic() and nb.IsInRing():

                if ca_candidate.GetIsAromatic():
                    return False, ""

        if nitrogen_atom.GetIsAromatic():
            return False, ""

        n_h_count = sum(1 for nb in nitrogen_atom.GetNeighbors() if nb.GetSymbol() == 'H')
        n_c_count = sum(1 for nb in nitrogen_atom.GetNeighbors() if nb.GetSymbol() == 'C')

        is_free_amine = n_h_count >= 2
        is_peptide_n = (n_h_count == 1 and n_c_count >= 2)

        if not (is_free_amine or is_peptide_n):
            return False, ""

        aa_code = _identify_amino_acid_by_sidechain(mol, ca_candidate.GetIdx())
        return True, aa_code

    except Exception:
        return False, ""

def _identify_amino_acid_by_sidechain(mol: Chem.Mol, ca_idx: int) -> str:
    """Identify an amino acid from its side-chain pattern."""
    try:
        ca = mol.GetAtomWithIdx(ca_idx)

        cb = None
        for neighbor in ca.GetNeighbors():
            if neighbor.GetSymbol() == 'C':

                is_co = False
                for n2 in neighbor.GetNeighbors():
                    if n2.GetSymbol() == 'O':
                        bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), n2.GetIdx())
                        if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                            is_co = True
                if not is_co:
                    cb = neighbor
                    break

        if cb is None:
            return "GLY"

        has_aromatic = False
        has_indole = False
        has_phenol = False
        has_sulfur = False
        has_thiol = False
        has_carboxyl = False
        has_amino = False
        has_guanidine = False
        has_imidazole = False
        sidechain_length = 0
        sidechain_carbons = 0

        visited = {ca_idx}
        queue = [(cb.GetIdx(), 1)]

        while queue:
            idx, depth = queue.pop(0)
            if idx in visited:
                continue
            visited.add(idx)
            sidechain_length = max(sidechain_length, depth)

            if mol.GetAtomWithIdx(idx).GetSymbol() == 'C':
                sidechain_carbons += 1

            atom = mol.GetAtomWithIdx(idx)
            sym = atom.GetSymbol()

            if atom.GetIsAromatic():
                has_aromatic = True

                if sym == 'N' and atom.IsInRing():
                    has_indole = True

                for n in atom.GetNeighbors():
                    if n.GetSymbol() == 'O' and not n.IsInRing():
                        has_phenol = True

            if sym == 'S':
                has_sulfur = True

                h_count = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == 'H')
                if h_count > 0:
                    has_thiol = True

            if sym == 'O' and not atom.GetIsAromatic():

                c_neighbor = None
                for n in atom.GetNeighbors():
                    if n.GetSymbol() == 'C':
                        c_neighbor = n
                if c_neighbor:
                    o_count = sum(1 for n in c_neighbor.GetNeighbors() if n.GetSymbol() == 'O')
                    if o_count >= 2:
                        has_carboxyl = True

            if sym == 'N' and not atom.GetIsAromatic():
                has_amino = True

                c_neighbor = None
                for n in atom.GetNeighbors():
                    if n.GetSymbol() == 'C':
                        n_count = sum(1 for n2 in n.GetNeighbors() if n2.GetSymbol() == 'N')
                        if n_count >= 3:
                            has_guanidine = True

            for n in atom.GetNeighbors():
                if n.GetIdx() not in visited and n.GetIdx() != ca_idx and n.GetSymbol() != 'H':
                    queue.append((n.GetIdx(), depth + 1))

        if has_indole:
            return "TRP"
        elif has_phenol:
            return "TYR"
        elif has_aromatic and not has_indole and not has_phenol:
            if has_imidazole:
                return "HIS"
            return "PHE"
        elif has_thiol:
            return "CYS"
        elif has_sulfur and not has_thiol:
            return "MET"
        elif has_guanidine:
            return "ARG"
        elif has_carboxyl and sidechain_carbons <= 2:

            return "ASP"
        elif has_carboxyl and sidechain_carbons > 2:

            return "GLU"
        elif has_amino and has_carboxyl:
            if sidechain_length <= 2:
                return "ASN"
            else:
                return "GLN"
        elif has_amino and sidechain_length >= 4:
            return "LYS"
        elif sidechain_length == 1:

            cb_neighbors = [n for n in cb.GetNeighbors() if n.GetIdx() != ca_idx]
            if any(n.GetSymbol() == 'O' for n in cb_neighbors):

                if any(n.GetSymbol() == 'C' for n in cb_neighbors):
                    return "THR"
                return "SER"
            return "ALA"
        elif sidechain_length == 2:

            cb_c_neighbors = [n for n in cb.GetNeighbors() if n.GetSymbol() == 'C' and n.GetIdx() != ca_idx]
            if len(cb_c_neighbors) == 2:

                both_terminal = True
                for c_neighbor in cb_c_neighbors:
                    c_c_count = sum(1 for n in c_neighbor.GetNeighbors()
                                   if n.GetSymbol() == 'C' and n.GetIdx() != cb.GetIdx())
                    if c_c_count >= 1:
                        both_terminal = False
                        break
                if both_terminal:
                    return "VAL"
            return "VAL"
        elif sidechain_length >= 3:

            cb_c_neighbors = [n for n in cb.GetNeighbors() if n.GetSymbol() == 'C' and n.GetIdx() != ca_idx]
            if len(cb_c_neighbors) == 2:

                return "ILE"
            elif len(cb_c_neighbors) == 1:

                return "LEU"
            return "LEU"

        return "default"

    except Exception:
        return "default"

def _is_anisotropic_system(mol: Chem.Mol, c_idx: int) -> Tuple[bool, str, float]:
    """Detect supported anisotropic structural environments."""
    try:
        atom = mol.GetAtomWithIdx(c_idx)
        if atom.GetSymbol() != 'C':
            return False, "", 0.0

        aromatic_neighbors = [n for n in atom.GetNeighbors() if n.GetIsAromatic()]
        if len(aromatic_neighbors) >= 2:

            delta = ANISOTROPIC_SYSTEMS.get("calixarene_cone", (1.43, 1.70, ""))[0]
            return True, "calixarene_cone", delta

        if atom.IsInRing():
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if c_idx in ring:

                    if len(ring) >= 8:
                        has_aromatic_in_ring = any(
                            mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring
                        )
                        if has_aromatic_in_ring:
                            delta = ANISOTROPIC_SYSTEMS.get("paracyclophane_8", (0.80, 1.00, ""))[0]
                            return True, "paracyclophane", delta

        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                for n2 in neighbor.GetNeighbors():
                    if n2.GetSymbol() == 'O':
                        bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), n2.GetIdx())
                        if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                            delta = ANISOTROPIC_SYSTEMS.get("carbonyl_alpha", (0.30, 0.50, ""))[0]
                            return True, "carbonyl_alpha", delta

        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                for bond in neighbor.GetBonds():
                    if bond.GetBondType() == rdchem.BondType.TRIPLE:
                        delta = ANISOTROPIC_SYSTEMS.get("alkyne_propargyl", (0.20, 0.30, ""))[0]
                        return True, "alkyne_propargyl", delta

        return False, "", 0.0

    except Exception:
        return False, "", 0.0

# Solvent helpers
def _normalize_solvent_name(solvent: str) -> str:
    """Normalize solvent names used by the lookup tables."""
    if not solvent:
        return "CDCl3"

    s = solvent.strip().lower()

    mappings = {
        "cdcl3": "CDCl3",
        "chloroform": "CDCl3",
        "chloroform-d": "CDCl3",
        "deuterated chloroform": "CDCl3",
        "dmso": "DMSO-d6",
        "dmso-d6": "DMSO-d6",
        "d6-dmso": "DMSO-d6",
        "dimethylsulfoxide": "DMSO-d6",
        "d2o": "D2O",
        "water": "D2O",
        "deuterium oxide": "D2O",
        "meod": "MeOD",
        "meod-d4": "MeOD",
        "methanol": "MeOD",
        "methanol-d4": "MeOD",
        "cd3od": "MeOD",
        "c6d6": "C6D6",
        "benzene": "C6D6",
        "benzene-d6": "C6D6",
        "deuterated benzene": "C6D6",
        "pyridine": "pyridine-d5",
        "pyridine-d5": "pyridine-d5",
        "c5d5n": "pyridine-d5",
        "acetone": "acetone-d6",
        "acetone-d6": "acetone-d6",
        "cd3cn": "CD3CN",
        "acetonitrile": "CD3CN",
        "acetonitrile-d3": "CD3CN",
        "thf": "THF-d8",
        "thf-d8": "THF-d8",
        "toluene": "toluene-d8",
        "toluene-d8": "toluene-d8",
        "cd2cl2": "CD2Cl2",
        "dcm": "CD2Cl2",
        "dichloromethane": "CD2Cl2",
        "dmf": "DMF-d7",
        "dmf-d7": "DMF-d7",
    }

    return mappings.get(s, solvent)

def get_solvent_additive_shift(solvent: str) -> float:
    """Return the additive solvent correction in ppm."""
    normalized = _normalize_solvent_name(solvent)
    return SOLVENT_ADDITIVE_SHIFT.get(normalized, SOLVENT_ADDITIVE_SHIFT.get("default", 0.0))

def get_solvent_multiplier(solvent: str) -> float:
    """Return the legacy solvent multiplier."""
    return 1.0

def is_aromatic_solvent(solvent: str) -> bool:
    """Check whether a solvent is classified as aromatic."""
    normalized = _normalize_solvent_name(solvent)
    props = SOLVENT_PROPERTIES.get(normalized, {})
    return props.get("aromatic", False)

# Public API
def get_diastereotopic_delta(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    heavy_idx: int,
    spectrometer_mhz: float = 700.0,
    solvent: str = "CDCl3",
    secondary_structure: str = "random_coil"
) -> Tuple[float, str]:
    """Estimate the chemical-shift difference between two diastereotopic protons."""
    try:
        heavy_atom = mol.GetAtomWithIdx(heavy_idx)
        heavy_sym = heavy_atom.GetSymbol()

        solvent_shift = get_solvent_additive_shift(solvent)
        solvent_note = ""
        if solvent_shift != 0.0:
            solvent_note = f"_solvent_{solvent_shift:+.2f}ppm"

        is_aniso, aniso_type, aniso_delta = _is_anisotropic_system(mol, heavy_idx)
        if is_aniso and aniso_delta >= 0.80:
            final_delta = aniso_delta + solvent_shift
            return final_delta, f"anisotropic_{aniso_type}{solvent_note}"

        if heavy_sym == 'C' and _is_vinyl_carbon(mol, heavy_idx):
            substituent = _get_vinyl_substituent_type(mol, heavy_idx)
            delta_data = VINYL_DELTA_BY_SUBSTITUENT.get(
                substituent,
                VINYL_DELTA_BY_SUBSTITUENT["default"]
            )
            base_delta = delta_data[0]
            final_delta = base_delta + solvent_shift
            return final_delta, f"vinyl_{substituent}{solvent_note}"

        if heavy_sym == 'C' and _is_epoxide_carbon(mol, heavy_idx):
            delta_data = EPOXIDE_DELTA.get("default", (0.30, 0.40, ""))
            base_delta = delta_data[0]
            final_delta = base_delta + solvent_shift
            return final_delta, f"epoxide{solvent_note}"

        if heavy_sym == 'C' and _is_aziridine_carbon(mol, heavy_idx):
            delta_data = AZIRIDINE_DELTA.get("default", (0.18, 0.25, ""))
            base_delta = delta_data[0]
            final_delta = base_delta + solvent_shift
            return final_delta, f"aziridine{solvent_note}"

        if heavy_sym == 'C' and _is_thiirane_carbon(mol, heavy_idx):
            delta_data = THIIRANE_DELTA.get("default", (0.20, 0.28, ""))
            base_delta = delta_data[0]
            final_delta = base_delta + solvent_shift
            return final_delta, f"thiirane{solvent_note}"

        if _is_nucleoside(mol, heavy_idx):
            delta_data = NUCLEOSIDE_DELTA.get("default", (0.14, 0.20, ""))
            base_delta = delta_data[0]
            final_delta = base_delta + solvent_shift
            return final_delta, f"nucleoside_H5{solvent_note}"

        if _is_crown_ether_carbon(mol, heavy_idx):
            delta_data = CROWN_ETHER_DELTA.get("default", (0.08, 0.15, ""))
            base_delta = delta_data[0]
            final_delta = base_delta + solvent_shift
            return final_delta, f"crown_ether{solvent_note}"

        contributions: List[Tuple[float, str]] = []

        chiral_dist = _find_chiral_distance(mol, heavy_idx)
        if chiral_dist > 0 and chiral_dist <= 6:
            delta_chiral = CH2_CHIRAL_DELTA_BY_DISTANCE.get(
                chiral_dist, (0.02, 0.05, "")
            )[0]
            contributions.append((delta_chiral, f"chiral_d{chiral_dist}"))

        alpha_co_found = False
        beta_co_found = False

        for nb in heavy_atom.GetNeighbors():
            if nb.GetSymbol() == 'C':

                for n2 in nb.GetNeighbors():
                    if n2.GetSymbol() == 'O':
                        bond = mol.GetBondBetweenAtoms(nb.GetIdx(), n2.GetIdx())
                        if bond and bond.GetBondTypeAsDouble() >= 2.0:
                            alpha_co_found = True
                            break

                if not alpha_co_found:
                    for n2 in nb.GetNeighbors():
                        if n2.GetSymbol() == 'C' and n2.GetIdx() != heavy_idx:
                            for n3 in n2.GetNeighbors():
                                if n3.GetSymbol() == 'O':
                                    bond = mol.GetBondBetweenAtoms(n2.GetIdx(), n3.GetIdx())
                                    if bond and bond.GetBondTypeAsDouble() >= 2.0:
                                        beta_co_found = True
                                        break
                            if beta_co_found:
                                break

        if alpha_co_found:
            contributions.append((0.15, "alpha_CO"))
        elif beta_co_found:
            contributions.append((0.12, "beta_CO"))

        for nb in heavy_atom.GetNeighbors():
            if nb.GetIsAromatic():

                ring_info = mol.GetRingInfo()
                is_heteroaromatic = False
                for ring in ring_info.AtomRings():
                    if nb.GetIdx() in ring:
                        if any(mol.GetAtomWithIdx(i).GetSymbol() in ('N', 'O', 'S')
                               for i in ring):
                            is_heteroaromatic = True
                            break
                if is_heteroaromatic:
                    contributions.append((0.12, "heteroAr"))
                else:
                    contributions.append((0.08, "Ar"))
                break

        for nb in heavy_atom.GetNeighbors():
            if nb.GetSymbol() == 'C' and not nb.GetIsAromatic():
                for bond in nb.GetBonds():
                    if bond.GetBondType() == rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nb)
                        if other.GetSymbol() == 'C' and not other.GetIsAromatic():
                            contributions.append((0.06, "allyl"))
                            break

        for nb in heavy_atom.GetNeighbors():
            if nb.GetSymbol() == 'C':
                for bond in nb.GetBonds():
                    if bond.GetBondType() == rdchem.BondType.TRIPLE:
                        contributions.append((0.10, "propargyl"))
                        break

        for nb in heavy_atom.GetNeighbors():
            if nb.GetSymbol() == 'C':
                nb_heavy = [n for n in nb.GetNeighbors() if n.GetSymbol() != 'H']
                if len(nb_heavy) == 4:
                    has_double = any(b.GetBondTypeAsDouble() >= 2 for b in nb.GetBonds())
                    if not has_double:
                        try:
                            ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                            sub_ranks = [ranks[n.GetIdx()] for n in nb_heavy]
                            if len(set(sub_ranks)) >= 3:
                                contributions.append((0.08, "prochiral_quat"))
                        except Exception:
                            pass

        is_sulfur, sulfur_type = _is_sulfoxide_sulfone_adjacent(mol, heavy_idx)
        if is_sulfur:
            delta_s = SULFOXIDE_SULFONE_DELTA.get(
                sulfur_type,
                SULFOXIDE_SULFONE_DELTA.get("default_sulfoxide", (0.12, 0.18, ""))
            )[0]
            contributions.append((delta_s, f"S_{sulfur_type[:4]}"))

        is_phosph, phosph_type = _is_phosphorus_adjacent(mol, heavy_idx)
        if is_phosph:
            delta_p = PHOSPHORUS_DELTA.get(phosph_type, PHOSPHORUS_DELTA["default"])[0]
            contributions.append((delta_p, f"P_{phosph_type[:4]}"))

        if _is_n_oxide_adjacent(mol, heavy_idx):
            contributions.append((0.18, "N_oxide"))

        if _is_silane_adjacent(mol, heavy_idx):
            contributions.append((0.05, "Si"))

        if _is_boron_adjacent(mol, heavy_idx):
            contributions.append((0.08, "B"))

        if heavy_atom.IsInRing():
            ring_type = _get_ring_type(mol, heavy_idx)

            if ring_type in CYCLIC_DELTA:
                delta_ring = CYCLIC_DELTA[ring_type][0]
            elif ring_type.startswith("ring_"):
                delta_ring = DEFAULT_DELTA["cyclic_ax_eq"]
            else:
                ring_size = 6 if "6" in ring_type else 5
                key = f"default_{ring_size}ring"
                delta_ring = CYCLIC_DELTA.get(key, (0.35, 0.45, ""))[0]

            contributions = [(d, r) for d, r in contributions if "chiral" not in r]
            contributions.append((delta_ring, f"ring_{ring_type[:6]}"))

        is_rotamer, rotamer_type = _is_rotamer_blocked(mol, heavy_idx)
        if is_rotamer:
            contributions.append((0.08, f"rotamer"))

        is_aa, aa_code = _is_amino_acid(mol, heavy_idx)
        if is_aa:

            if aa_code in AMINO_ACID_BETA_DELTA:
                delta_aa = AMINO_ACID_BETA_DELTA[aa_code][0]
            else:
                delta_aa = AMINO_ACID_BETA_DELTA.get("default", (0.04, 0.01, ""))[0]

            if secondary_structure in SECONDARY_STRUCTURE_EFFECT:
                ss_effect = SECONDARY_STRUCTURE_EFFECT[secondary_structure]
                delta_aa += ss_effect.get("delta_increase", 0.0)

            contributions = [(d, r) for d, r in contributions
                           if r not in ("chiral_d1", "chiral_d2", "alpha_CO", "beta_CO", "heteroAr", "Ar")]
            contributions.append((delta_aa, f"AA_{aa_code}"))

        if not is_aa:
            for nb in heavy_atom.GetNeighbors():
                if nb.GetSymbol() == 'C' and nb.IsInRing():
                    ring_info = mol.GetRingInfo()
                    for ring in ring_info.AtomRings():
                        if nb.GetIdx() in ring and len(ring) == 6:

                            n_count = sum(1 for i in ring
                                         if mol.GetAtomWithIdx(i).GetSymbol() == 'N')
                            co_count = 0
                            for i in ring:
                                atom_i = mol.GetAtomWithIdx(i)
                                if atom_i.GetSymbol() == 'C':
                                    for n2 in atom_i.GetNeighbors():
                                        if n2.GetSymbol() == 'O' and n2.GetIdx() not in ring:
                                            bond = mol.GetBondBetweenAtoms(i, n2.GetIdx())
                                            if bond and bond.GetBondTypeAsDouble() >= 2:
                                                co_count += 1
                            if n_count >= 2 and co_count >= 2:

                                contributions.append((0.10, "DKP_rigid"))
                                break
                    else:
                        continue
                    break

        if not contributions:

            base_delta = DEFAULT_DELTA.get("flexible_chain", 0.08)
            final_delta = base_delta + solvent_shift
            return final_delta, f"flexible_chain{solvent_note}"

        total_delta = sum(d for d, r in contributions)

        total_delta += solvent_shift

        if heavy_atom.IsInRing():

            total_delta = min(total_delta, 0.70)
        else:

            total_delta = min(total_delta, 0.60)

        reasons = "+".join(r for d, r in contributions)
        if solvent_note:
            reasons += solvent_note

        return round(total_delta, 3), reasons

    except Exception as e:
        return DEFAULT_DELTA["unknown"], f"error_{str(e)[:20]}"

def get_diastereotopic_delta_simple(
    mol: Chem.Mol,
    heavy_idx: int,
    context: str = "auto",
    solvent: str = "CDCl3"
) -> float:
    """Return a simplified diastereotopic chemical-shift difference estimate."""
    solvent_factor = get_solvent_multiplier(solvent)

    if context == "vinyl":
        return DEFAULT_DELTA["vinyl_geminal"] * solvent_factor
    elif context == "cyclic":
        return DEFAULT_DELTA["cyclic_ax_eq"] * solvent_factor
    elif context == "chiral":
        return DEFAULT_DELTA["ch2_near_chiral_1"] * solvent_factor
    elif context == "amino_acid":
        return DEFAULT_DELTA["amino_acid"] * solvent_factor
    elif context == "auto":

        try:
            heavy_atom = mol.GetAtomWithIdx(heavy_idx)

            if _is_vinyl_carbon(mol, heavy_idx):
                return DEFAULT_DELTA["vinyl_geminal"] * solvent_factor

            if _is_in_ring(mol, heavy_idx):
                return DEFAULT_DELTA["cyclic_ax_eq"] * solvent_factor

            chiral_dist = _find_chiral_distance(mol, heavy_idx)
            if chiral_dist > 0:
                key = f"ch2_near_chiral_{min(chiral_dist, 3)}"
                return DEFAULT_DELTA.get(key, 0.12) * solvent_factor

            return DEFAULT_DELTA["flexible_chain"] * solvent_factor
        except Exception:
            pass

    return DEFAULT_DELTA["unknown"] * solvent_factor

# Exported constants
DELTA_CYCLOHEXANE = 0.48

DELTA_CHIRAL_ADJ = 0.08

DELTA_DEFAULT = 0.07
