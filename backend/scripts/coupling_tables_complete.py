# =============================================================================
# coupling_tables_complete.py
# Complete coupling tables for 100% coverage
# =============================================================================

"""
Module of empirical tables for NMR J couplings.

Coverage:
- Alkynes (C≡C-H)
- Aldehydes (CHO)
- Amides and peptides (CONH)
- Phosphines (P-H)
- Boranes (B-H)
- Silanes (Si-H)
- Organometallics
- Steroids
- Terpenes
- Macrocycles
- Carbohydrates (already covered, reference)
- Nucleosides/Nucleotides

Author: [Your name]
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

# Import specialized peptide module (avoids duplication)
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
# SECTION 1: ALKYNES (C≡C-H)
# =============================================================================

_ALKYNE_COUPLING_TABLE: Dict[str, Dict[str, float]] = {
    # 2J couplings through the triple bond
    '2J_terminal': {
        'C≡C-H': 2.5,           # Basic terminal alkyne
        'HC≡C-CH3': 2.7,        # Propyne
        'HC≡C-Ph': 2.4,         # Phenylacetylene
        'HC≡C-CH2': 2.6,        # Alkyne with CH2
        'HC≡C-C(O)': 2.2,       # Conjugated carbonyl alkyne
    },
    
    # 4J propargylic couplings (H-C-C≡C-H)
    '4J_propargylic': {
        'default': 2.5,          # Typical H-CH2-C≡C-H
        'allenic': 7.0,          # Allenic system H2C=C=CH2
        'conjugated': 3.0,       # Conjugated with aromatic
    },
    
    # 5J couplings
    '5J_alkyne': {
        'extended_W': 1.5,       # Extended W through triple bond
        'default': 0.5,
    },
    
    # Couplings with substituents
    '3J_alkyne_substituent': {
        'CH3-C≡C': 0.0,          # No coupling (no H on C triple)
        'CH2-C≡C-H': 2.8,        # Through CH2
    },
    
    # 3J across triple bond (R-C≡C-H type)
    # Note: Triple bonds are linear, no cis/trans isomerism exists
    # Literature: 3J(H-C≡C-H) ≈ 9-10 Hz (single value)
    '3J_across_triple': {
        'HC≡CH': 9.6,            # Acetylene (experimental reference)
        'R-C≡C-H': 9.0,          # Substituted terminal alkynes
        'default': 9.5,          # Average value
    },
}

# Allènes (cumulènes)
_ALLENE_COUPLING_TABLE: Dict[str, float] = {
    '4J_allenic': 7.0,           # Characteristic H2C=C=CH2
    '5J_allenic': 3.0,           # Extended coupling
    '2J_allenic_gem': -6.0,      # Geminal on allene
}


def get_alkyne_coupling(
    mol: Chem.Mol, 
    h1_idx: int, 
    h2_idx: int,
    bond_path: List[int]
) -> Optional[float]:
    """
    Calculate coupling for alkyne systems.
    
    Args:
        mol: RDKit molecule
        h1_idx, h2_idx: Proton indices
        bond_path: Bond path between protons
    
    Returns:
        J in Hz or None
    """
    # Detect if triple bond in path
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
        # Check if allenic
        is_allenic = _is_allenic_system(mol, bond_path)
        if is_allenic:
            return _ALLENE_COUPLING_TABLE['4J_allenic']
        return _ALKYNE_COUPLING_TABLE['4J_propargylic'].get('default', 2.5)
    
    elif n_bonds == 5:
        return _ALKYNE_COUPLING_TABLE['5J_alkyne'].get('default', 0.5)
    
    return None


def _is_allenic_system(mol: Chem.Mol, bond_path: List[int]) -> bool:
    """Detect an allenic system (C=C=C)."""
    for i in range(len(bond_path) - 2):
        a1, a2, a3 = bond_path[i], bond_path[i+1], bond_path[i+2]
        
        bond1 = mol.GetBondBetweenAtoms(a1, a2)
        bond2 = mol.GetBondBetweenAtoms(a2, a3)
        
        if (bond1 and bond2 and 
            bond1.GetBondType() == Chem.BondType.DOUBLE and
            bond2.GetBondType() == Chem.BondType.DOUBLE):
            
            # Check that central atom has sp hybridization
            central = mol.GetAtomWithIdx(a2)
            if central.GetHybridization() == Chem.HybridizationType.SP:
                return True
    
    return False


# =============================================================================
# SECTION 2: ALDEHYDES (CHO)
# =============================================================================

_ALDEHYDE_COUPLING_TABLE: Dict[str, Dict[str, float]] = {
    # 3J aldehyde-alphaH (per JACS 1963 reference data)
    # Literature: acetaldehyde ~2.8 Hz
    '3J_CHO_alphaH': {
        'default': 2.8,              # Acetaldehyde reference: 2.8 Hz
        'branched': 2.0,             # CHO-CH< (pivaldehyde-like)
        'conjugated': 3.2,           # CHO conjugated
        'aromatic': 0.0,             # Benzaldehyde (no alphaH)
    },
    
    # 4J aldehyde-betaH
    '4J_CHO_betaH': {
        'trans_planar': 1.2,         # W configuration
        'gauche': 0.3,
        'default': 0.6,
    },
    
    # Couplings in aromatic aldehydes
    'aromatic_aldehyde': {
        '4J_CHO_ortho': 0.8,         # Coupling with ortho H
        '5J_CHO_meta': 0.3,          # Coupling with meta H (weak)
    },
    
    # Alpha,beta-unsaturated aldehydes
    # Note: 15.5 Hz is the vinyl H-C=C-H coupling, NOT CHO->vinyl coupling
    # Literature: 3J(CHO, =CH) ≈ 7-8 Hz (Pardi et al.)
    'unsaturated_aldehyde': {
        '3J_CHO_vinyl_s_trans': 7.8,  # CHO-CH=CH, s-trans conformer (major)
        '3J_CHO_vinyl_s_cis': 4.5,    # CHO-CH=CH, s-cis conformer
        '4J_CHO_vinyl_beta': 1.2,     # Long-range coupling (unchanged)
    },
}


def get_aldehyde_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    distance: int,
    high_resolution: bool = False
) -> Optional[float]:
    """
    Calculate coupling involving an aldehyde proton.
    
    Note: For aromatic aldehydes, 4J and 5J couplings with
    ring H are very weak (<1 Hz) and generally unresolved.
    The CHO typically appears as a singlet (possibly broadened).
    
    IMPORTANT: Excludes formamides (N-CHO) because CHO has no true alpha-H.
    The path H(CHO)-C=O-N-C(CH3)-H is a 4J through N, not a 3J.
    
    IMPORTANT: For aliphatic aldehydes, we verify that other H
    is indeed an alpha-H (on the carbon directly bonded to CHO), not a beta-H
    (on the next carbon). Only alpha-H gives an observable 3J.
    """
    # Identify the aldehyde proton and the CHO carbon
    aldehyde_h = None
    aldehyde_c = None
    other_h = None
    
    for h_idx in [h1_idx, h2_idx]:
        h_atom = mol.GetAtomWithIdx(h_idx)
        for neighbor in h_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                # Check if C=O (aldehyde)
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
    
    # NEW: Exclude formamides (N-CHO)
    # A formamide has an N attached to carbonyl: N-C(=O)-H
    # The CHO of a formamide does not couple with N-CH3 (4J through N)
    if _is_formamide(mol, aldehyde_c.GetIdx()):
        # The CHO of a formamide is a singlet
        # Couplings with N-alkyl are 4J (<0.5 Hz) not observable
        return None
    
    # Check if it's an aromatic aldehyde (attached to aromatic ring)
    is_aromatic_aldehyde = False
    for neighbor in aldehyde_c.GetNeighbors():
        if neighbor.GetSymbol() != 'O' and neighbor.GetIsAromatic():
            is_aromatic_aldehyde = True
            break
    
    # For aromatic aldehydes, only return couplings in high resolution
    # because CHO appears as singlet (4J/5J are very weak)
    if is_aromatic_aldehyde:
        if high_resolution and distance == 4:
            # Very weak 4J coupling with ortho H of ring
            return _ALDEHYDE_COUPLING_TABLE['aromatic_aldehyde'].get('4J_CHO_ortho', 0.8)
        return None
    
    # =========================================================================
    # NEW: Verify that other_h is indeed an alpha-H (on the carbon
    # directly bonded to carbonyl CHO), not a beta-H or further
    # =========================================================================
    # Find the alpha carbon (directly bonded to carbonyl CHO)
    alpha_carbon = None
    for neighbor in aldehyde_c.GetNeighbors():
        if neighbor.GetSymbol() == 'C':  # The alpha carbon
            alpha_carbon = neighbor
            break
    
    if alpha_carbon is None:
        # No alpha carbon (e.g., formaldehyde H-CHO or CHO attached to O/N)
        return None
    
    # Verify that other_h is on alpha_carbon
    other_h_atom = mol.GetAtomWithIdx(other_h)
    other_h_parent = None
    for neighbor in other_h_atom.GetNeighbors():
        if neighbor.GetSymbol() == 'C':
            other_h_parent = neighbor
            break
    
    if other_h_parent is None:
        return None
    
    # CRITICAL: Only the CHO <-> alpha-H coupling is allowed
    # If other_h is not on alpha_carbon, it's a beta-H or further -> ignore
    if other_h_parent.GetIdx() != alpha_carbon.GetIdx():
        # other_h is on a different carbon than the alpha-carbon
        # It's a beta-H, gamma-H, etc. -> 4J or more, very weak, ignore
        return None
    
    # At this point, other_h is confirmed as alpha-H
    # The 3J(CHO-alphaH) coupling is typically 2.5-3.0 Hz
    if distance == 3:
        # 3J with alpha-H (aliphatic aldehydes) - confirmed
        return _ALDEHYDE_COUPLING_TABLE['3J_CHO_alphaH'].get('default', 2.8)
    
    # For distance == 4 (4J beta-H), we already filtered above
    # But in case we arrive here with distance != 3, ignore
    return None


# =============================================================================
# SECTION 3: AMIDES AND PEPTIDES
# =============================================================================
# NOTE: Complete peptide logic is in peptide_couplings.py
# Here we only keep a simplified table and a wrapper
# =============================================================================

# Simplified table for fallback if peptide_couplings.py is not available
def _is_formamide(mol: Chem.Mol, c_idx: int) -> bool:
    """
    Check if a carbon is in a formamide (N-C(=O)H).
    Includes N,N-disubstituted formamides like DMF.

    Args:
        mol: RDKit molecule
        c_idx: Carbon index to test

    Returns:
        True if it's a formamide CHO
    """
    c = mol.GetAtomWithIdx(c_idx)

    # The carbon must have at least one H (CHO group)
    has_H_on_C = (c.GetTotalNumHs() >= 1)

    # Check if C=O and has N neighbor
    has_double_O = False
    has_N = False

    for neighbor in c.GetNeighbors():
        if neighbor.GetSymbol() == 'O':
            bond = mol.GetBondBetweenAtoms(c_idx, neighbor.GetIdx())
            if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                has_double_O = True
        elif neighbor.GetSymbol() == 'N':
            # Any N linked to C=O (with or without H on N)
            # This includes DMF (N,N-dimethylformamide)
            has_N = True

    return has_double_O and has_N and has_H_on_C


def _is_true_peptide_bond(mol: Chem.Mol, atom_N: Chem.Atom) -> bool:
    """
    Check if N is in a peptide bond N-C(=O)-C
    (and not just a tertiary amide like DMF)
    """
    if atom_N.GetSymbol() != 'N':
        return False
    
    # Check that N is attached to a C=O
    carbonyl_C = None
    for n in atom_N.GetNeighbors():
        if n.GetSymbol() == 'C':
            for nn in n.GetNeighbors():
                if nn.GetSymbol() == 'O':
                    bond = mol.GetBondBetweenAtoms(n.GetIdx(), nn.GetIdx())
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                        carbonyl_C = n
                        break
    
    if carbonyl_C is None:
        return False
    
    # Check that this C=O also has a C-alpha (attached to a C with H)
    for neighbor in carbonyl_C.GetNeighbors():
        if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() != atom_N.GetIdx():
            # Ce C doit avoir au moins un H (Cα)
            if neighbor.GetTotalNumHs() >= 1:
                return True
    
    return False


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
    Calculate coupling in a peptide.

    REFACTORED: Module peptide restrictif - s'active UNIQUEMENT pour les vrais
    résidus peptidiques acycliques. Les cycles (DKP, lactames, proline) et les
    vinyles doivent être traités par leurs modules spécifiques.

    Conditions d'activation (TOUTES doivent être vraies):
    1. Distance exactement 3 liaisons (³J)
    2. Aucun des deux carbones n'est dans un cycle
    3. C'est un vrai Cα peptidique (sp3 portant NH et C=O)
    4. Pas de carbone sp2 (vinyl/aromatique)

    Returns:
        J in Hz or None if not applicable
    """
    h1 = mol.GetAtomWithIdx(h1_idx)
    h2 = mol.GetAtomWithIdx(h2_idx)

    # Get heavy atoms attached to protons
    heavy1_list = [n for n in h1.GetNeighbors() if n.GetSymbol() != 'H']
    heavy2_list = [n for n in h2.GetNeighbors() if n.GetSymbol() != 'H']

    if not heavy1_list or not heavy2_list:
        return None

    heavy1 = heavy1_list[0]
    heavy2 = heavy2_list[0]

    # =========================================================================
    # CONDITION 1: Distance exactement 3 liaisons (³J H-C-C-H ou H-N-C-H)
    # Le module peptide ne gère PAS les ⁴J ou ⁵J
    # =========================================================================
    try:
        path = Chem.GetShortestPath(mol, h1_idx, h2_idx)
        n_bonds = len(path) - 1
        if n_bonds != 3:
            return None  # Pas un ³J → autre module
    except Exception:
        return None

    # =========================================================================
    # CONDITION 2: Pas dans un cycle (DKP, lactames, proline → saturated_heterocycle)
    # Les cycles ont des J fixes par tables, pas des rotamères
    # =========================================================================
    if heavy1.GetSymbol() == 'C' and heavy2.GetSymbol() == 'C':
        if heavy1.IsInRing() and heavy2.IsInRing():
            return None  # Les deux dans un cycle → module saturated_heterocycle

    # =========================================================================
    # CONDITION 3: Pas de carbone sp2 (vinyl, aromatique → vinylic ou aromatic)
    # =========================================================================
    def _is_sp2_or_vinyl(atom) -> bool:
        """Check if atom is sp2 or attached to C=C."""
        if atom.GetSymbol() != 'C':
            return False
        if atom.GetIsAromatic():
            return True
        hyb = str(atom.GetHybridization())
        if 'SP2' in hyb:
            # Check if it's a C=C (vinyl) not C=O (carbonyl)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                        return True  # Vinyl carbon
        return False

    if _is_sp2_or_vinyl(heavy1) or _is_sp2_or_vinyl(heavy2):
        return None  # Vinyl → module vinylic

    # =========================================================================
    # CONDITION 4: C'est un VRAI résidu peptidique acyclique
    # Un Cα peptidique est sp3, attaché à un N-amide ET à un C=O
    # =========================================================================
    def _is_real_alpha_carbon(atom) -> bool:
        """Check if atom is a true peptidic alpha carbon."""
        if atom.GetSymbol() != 'C':
            return False
        if atom.IsInRing():
            return False  # Cα cyclique → pas un peptide linéaire

        has_amide_N = False
        has_carbonyl = False

        for neighbor in atom.GetNeighbors():
            # Check for amide N (N attached to C=O)
            if neighbor.GetSymbol() == 'N':
                for nn in neighbor.GetNeighbors():
                    if nn.GetSymbol() == 'C' and nn.GetIdx() != atom.GetIdx():
                        for nnn in nn.GetNeighbors():
                            if nnn.GetSymbol() == 'O':
                                bond = mol.GetBondBetweenAtoms(nn.GetIdx(), nnn.GetIdx())
                                if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                                    has_amide_N = True
                                    break

            # Check for attached carbonyl (C=O)
            if neighbor.GetSymbol() == 'C':
                for nn in neighbor.GetNeighbors():
                    if nn.GetSymbol() == 'O':
                        bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nn.GetIdx())
                        if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                            has_carbonyl = True
                            break

        return has_amide_N and has_carbonyl

    # At least one of the carbons must be a real Cα
    is_h1_on_alpha = _is_real_alpha_carbon(heavy1) if heavy1.GetSymbol() == 'C' else False
    is_h2_on_alpha = _is_real_alpha_carbon(heavy2) if heavy2.GetSymbol() == 'C' else False

    # For NH-Cα coupling, check if one is N (amide) and other is Cα
    is_hn_ca_coupling = False
    if heavy1.GetSymbol() == 'N' and is_h2_on_alpha:
        is_hn_ca_coupling = True
    elif heavy2.GetSymbol() == 'N' and is_h1_on_alpha:
        is_hn_ca_coupling = True

    # For Cα-Cβ coupling, one must be Cα and other must be adjacent to Cα
    is_ca_cb_coupling = is_h1_on_alpha or is_h2_on_alpha

    if not (is_hn_ca_coupling or is_ca_cb_coupling):
        return None  # Not a real peptide coupling

    # =========================================================================
    # Passed all filters - this is a true peptide coupling
    # Continue with peptide module calculation
    # =========================================================================

    def _is_peptide_proton(h_atom) -> bool:
        """
        Check if a proton is on the peptide backbone.
        """
        neighbors = list(h_atom.GetNeighbors())
        if not neighbors:
            return False
        c = neighbors[0]

        # Proton on amide NH
        if c.GetSymbol() == 'N':
            for neighbor in c.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    for nn in neighbor.GetNeighbors():
                        if nn.GetSymbol() == 'O':
                            bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nn.GetIdx())
                            if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                                return True
            return False

        # Proton on a carbon - already validated above
        if c.GetSymbol() == 'C':
            # Check if carbon is on the amide backbone:
            # - C-alpha: attached to N AND to a C=O
            # - C-beta: attached to C-alpha
            for neighbor in c.GetNeighbors():
                if neighbor.GetSymbol() == 'N':
                    # C-alpha: attached to N - check if this N is an amide
                    for nn in neighbor.GetNeighbors():
                        if nn.GetSymbol() == 'C' and nn.GetIdx() != c.GetIdx():
                            for nnn in nn.GetNeighbors():
                                if nnn.GetSymbol() == 'O':
                                    bond = mol.GetBondBetweenAtoms(nn.GetIdx(), nnn.GetIdx())
                                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                                        return True
                elif neighbor.GetSymbol() == 'C':
                    # Check if neighboring C is a C-alpha (attached to amide N)
                    for nn in neighbor.GetNeighbors():
                        if nn.GetSymbol() == 'N':
                            for nnn in nn.GetNeighbors():
                                if nnn.GetSymbol() == 'C' and nnn.GetIdx() != neighbor.GetIdx():
                                    for nnnn in nnn.GetNeighbors():
                                        if nnnn.GetSymbol() == 'O':
                                            bond = mol.GetBondBetweenAtoms(nnn.GetIdx(), nnnn.GetIdx())
                                            if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                                                return True  # C-beta adjacent to C-alpha
        
        return False
    
    # Verify that BOTH protons are peptidic
    if not (_is_peptide_proton(h1) and _is_peptide_proton(h2)):
        return None  # Not a peptidic coupling
    
    # Exclude simple formamides (N-CHO)
    # A simple formamide has pattern N-C(=O)-H without alpha-carbon
    for h_idx, h_atom in [(h1_idx, h1), (h2_idx, h2)]:
        neighbors = list(h_atom.GetNeighbors())
        if not neighbors:
            continue
        c = neighbors[0]
        
        # Check if it's a formamide CHO (C=O attached to N with an H)
        if c.GetSymbol() == 'C':
            has_double_O = False
            has_N = False
            for neighbor in c.GetNeighbors():
                if neighbor.GetSymbol() == 'O':
                    bond = mol.GetBondBetweenAtoms(c.GetIdx(), neighbor.GetIdx())
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                        has_double_O = True
                elif neighbor.GetSymbol() == 'N':
                    has_N = True
            
            # It's a formamide (N-CHO) if we have C=O and N attached
            if has_double_O and has_N:
                # Check if the H is directly on this C (CHO)
                h_on_c = any(n.GetIdx() == h_idx for n in c.GetNeighbors())
                if h_on_c:
                    # It's a formamide CHO - no peptide coupling applicable
                    return None
    
    # NEW: Verify it's a true peptide bond (not DMF/formamide)
    h1_neighbor = list(h1.GetNeighbors())[0] if list(h1.GetNeighbors()) else None
    h2_neighbor = list(h2.GetNeighbors())[0] if list(h2.GetNeighbors()) else None
    
    has_peptide_bond = False
    if h1_neighbor and h1_neighbor.GetSymbol() == 'N':
        has_peptide_bond = _is_true_peptide_bond(mol, h1_neighbor)
    elif h2_neighbor and h2_neighbor.GetSymbol() == 'N':
        has_peptide_bond = _is_true_peptide_bond(mol, h2_neighbor)
    
    if not has_peptide_bond:
        return None  # Pas une vraie liaison peptidique
    
    # Utiliser le module peptide robuste si disponible
    if _PEPTIDE_MODULE_AVAILABLE:
        try:
            return get_peptide_coupling_robust(
                mol, h1_idx, h2_idx,
                phi_angle=dihedral_angle,
                secondary_structure=secondary_structure
            )
        except Exception as e:
            print(f"Warning: peptide module failed ({e}), using fallback")
            pass  # Fallback below
    
    # Simplified fallback if module is not available
    
    # Find heavy atoms
    c1 = list(h1.GetNeighbors())[0]
    c2 = list(h2.GetNeighbors())[0]
    
    # Check if it's an NH-C-alpha-H coupling
    is_NH = False
    is_CaH = False
    
    for h_idx, c in [(h1_idx, c1), (h2_idx, c2)]:
        # Amide NH
        if c.GetSymbol() == 'N':
            for neighbor in c.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    for nn in neighbor.GetNeighbors():
                        if nn.GetSymbol() == 'O':
                            bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nn.GetIdx())
                            if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                                is_NH = True
                                break
        
        # C-alpha (attached to N and C=O)
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
        # NH-C-alpha-H coupling - use Karplus if angle available
        if dihedral_angle is not None:
            # Pardi, Billeter & Wüthrich (J. Mol. Biol. 1984, 180, 741):
            # 3J(NH,Hα) = 6.4·cos²(φ-60°) - 1.4·cos(φ-60°) + 1.9
            # Note: Bystrov (1976) uses different coefficients (A=9.4, B=-1.1, C=0.4)
            phi_rad = np.radians(dihedral_angle - 60)
            J = 6.4 * np.cos(phi_rad)**2 - 1.4 * np.cos(phi_rad) + 1.9
            return max(0, J)
        else:
            # Value by secondary structure
            return _AMIDE_COUPLING_TABLE_FALLBACK['3J_NH_alphaH'].get(
                secondary_structure, 6.5
            )
    
    # Default H-alpha - H-beta coupling
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
    # ¹J(B-H) - Note: ¹¹B has I=3/2
    '1J_BH': {
        'BH3_adduct': 103.0,         # BH3·L (adducts: BH3·THF, BH3·SMe2)
        'B2H6_terminal': 133.5,      # Diborane, terminal B-H
        'B2H6_bridge': 46.3,         # Diborane, B-H-B bridge
        'BH4-': 81.0,                # Borohydride
        'R-BH2': 110.0,
        'R2-BH': 120.0,
        'catechol_borane': 140.0,
        'pinacol_borane': 145.0,
        'default': 110.0,
    },
    
    # 2J(B-C-H)
    '2J_BCH': {
        'default': 3.0,
        'vinyl': 5.0,
    },
    
    # 3J(B-C-C-H)
    '3J_BCCH': {
        'default': 3.0,
    },
    
    # Polyhedral boranes (carboranes, etc.)
    'polyhedral': {
        '2J_BB': 15.0,           # B-B coupling in clusters
        '3J_BH_vertex': 5.0,     # Vertex B-H coupling
        '2J_CH_carborane': 8.0,  # C-H in carborane
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
    
    # 2J(Si-C-H)
    '2J_SiCH': {
        'SiCH3': 4.5,
        'SiCH2': 3.5,
        'default': 4.0,
    },
    
    # 3J(Si-C-C-H)
    '3J_SiCCH': {
        'trans': 12.0,
        'gauche': 3.0,
        'average': 6.0,
        'default': 6.0,
    },
    
    # 4J
    '4J_SiCCCH': {
        'default': 1.0,
    },
    
    # 29Si satellites (4.7% abundance)
    'Si29_satellites': {
        '1J_29Si_H': 200.0,          # Visible as satellites
    },
    
    # Fluorosilanes
    'fluorosilane': {
        '1J_SiF_H': 275.0,           # Si-F increases 1J(Si-H)
        '2J_FSi_CH': 8.0,            # F-Si-C-H coupling
        '3J_FSi_CCH': 4.0,           # Long-range F-Si coupling
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
        # Detect silane type
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
# SECTION 7: ORGANOMETALLICS
# =============================================================================

_ORGANOMETALLIC_COUPLING_TABLE: Dict[str, Dict[str, Any]] = {
    # Ferrocene and metallocenes
    'ferrocene': {
        '3J_HH_ring': 2.5,           # Coupling on Cp ring
        '4J_HH_across': 1.5,         # Between rings
    },
    
    # Carbonyl complexes
    'metal_carbonyl': {
        '2J_M_CO_H': 3.0,            # M-CO coupling (weak)
    },
    
    # Phosphine complexes
    'metal_phosphine': {
        '2J_P_M_H': 5.0,
        '3J_P_M_C_H': 8.0,
    },
    
    # Metal alkyls
    'metal_alkyl': {
        '2J_M_CH': 8.0,              # E.g.: Ti-CH3
        '3J_M_C_CH': 6.0,
    },
    
    # Metal hydrides
    'metal_hydride': {
        '1J_M_H': {
            'default': 30.0,
            'Pt': 850.0,             # 195Pt-H very large
            'Rh': 25.0,
            'Ir': 35.0,
            'Ru': 30.0,
            'W': 45.0,               # 183W-H coupling
            'Mo': 40.0,              # 95Mo-H coupling
            'Os': 35.0,              # 187Os-H coupling
            'Re': 50.0,              # 187Re-H coupling
        },
        '2J_H_M_H': -5.0,            # Geminal cis coupling
        '2J_H_M_H_trans': -15.0,     # Trans larger
    },
    
    # Eta-2 complexes (olefins)
    'olefin_complex': {
        '3J_HH_coordinated': 8.0,    # Coordinated olefin
    },
}


def get_organometallic_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    metal_symbol: Optional[str] = None
) -> Optional[float]:
    """Calculate coupling in organometallic complexes."""
    
    # Auto-detect metal if not specified
    if metal_symbol is None:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in ['Fe', 'Ru', 'Os', 'Co', 'Rh', 'Ir', 'Ni', 'Pd', 'Pt']:
                metal_symbol = atom.GetSymbol()
                break
    
    # Metal-specific values
    metal_specific = {
        'Pt': 25.0,  # Pt-H complexes show stronger couplings
        'Rh': 15.0,  # Typical Rh complexes
        'Fe': 8.0,   # Ferrocene and derivatives
    }
    
    return metal_specific.get(metal_symbol, 5.0)  # Generic default


# =============================================================================
# SECTION 8: STEROIDS
# =============================================================================

_STEROID_COUPLING_TABLE: Dict[str, Dict[str, float]] = {
    # A/B ring junctions (cis or trans)
    'AB_junction': {
        '3J_H5_H6ax_trans': 12.0,    # Trans diaxial
        '3J_H5_H6eq_trans': 4.5,
        '3J_H5_H6ax_cis': 5.0,       # Cis junction
        '3J_H5_H6eq_cis': 3.5,
    },
    
    # B/C junction (always trans)
    'BC_junction': {
        '3J_H8_H9': 11.0,            # Trans diaxial
        '3J_H8_H14': 10.5,
    },
    
    # C/D junction
    'CD_junction': {
        '3J_H14_H15ax': 12.0,
        '3J_H14_H15eq': 4.0,
    },
    
    # Characteristic W-couplings
    'W_couplings': {
        '4J_H1ax_H3ax': 2.5,         # W in ring A
        '4J_H2ax_H4ax': 2.5,
        '4J_H11ax_H12ax': 2.0,       # W in ring C
        '4J_H15ax_H17': 1.5,
    },
    
    # Angular methyls
    'angular_methyl': {
        '4J_CH3_18_H12': 0.8,        # CH3 C-18 with H-12
        '4J_CH3_19_H1ax': 0.5,       # CH3 C-19 with H-1ax
        '4J_CH3_19_H2ax': 0.5,
    },
    
    # Aromatic A ring (estrogens)
    'aromatic_A_ring': {
        '3J_H1_H2': 8.5,             # Aromatic ortho
        '4J_H1_H4': 2.5,             # Meta
        '4J_H2_H4': 2.5,
    },
    
    # Side chain (cholesterol type)
    'sidechain': {
        '3J_H17_H20': 9.0,           # Trans
        '3J_H20_H22': 6.5,
        '3J_CH3_21_H20': 6.5,        # Methyl-21 with H-20
    },
    
    # ==========================================================================
    # EXTENDED STEROID COUPLINGS - Unusual functional groups
    # ==========================================================================
    
    # Δ4-3-keto steroids (enone system, e.g., testosterone, progesterone)
    'delta4_3keto': {
        '3J_H4_H6ax': 2.5,           # Allylic (H4 is vinylic)
        '3J_H4_H6eq': 2.0,
        '4J_H4_H7ax': 1.2,           # Long-range through enone
        '4J_H2_H4': 1.5,             # W-coupling
        '3J_H1_H2ax': 5.0,           # α to carbonyl
        '3J_H1_H2eq': 2.5,
    },
    
    # Δ5-steroids (double bond C5=C6, e.g., cholesterol)
    'delta5': {
        '3J_H4_H5': 5.0,             # Allylic
        '4J_H4_H7': 2.0,             # Homoallylic
        '4J_H3_H6': 1.5,             # Cross-ring
        '3J_H6_H7ax': 4.5,
        '3J_H6_H7eq': 2.5,
    },
    
    # Δ1-steroids (additional double bond C1=C2)
    'delta1': {
        '3J_H1_H2': 10.0,            # Vinyl coupling
        '4J_H2_H4': 1.8,             # Allylic
    },
    
    # 11β-hydroxy/keto (cortisol-like)
    '11_substituted': {
        '3J_H9_H11ax': 11.0,         # Trans-diaxial
        '3J_H9_H11eq': 4.0,          # Ax-eq
        '3J_H11_H12ax': 12.0,
        '3J_H11_H12eq': 4.5,
    },
    
    # 17α-hydroxy (corticoids)
    '17_hydroxy': {
        '3J_H16ax_H17': 9.0,         # Depends on 17α or 17β
        '3J_H16eq_H17': 5.0,
    },
    
    # Spiroketal (spironolactone-like)
    'spiroketal': {
        '3J_H15_H16': 6.5,           # Across spirocenter
        '3J_ring_to_lactone': 7.0,
    },
    
    # Epoxy steroids (5,6-epoxide)
    'epoxy_56': {
        '3J_H5_H6_alpha': 4.0,       # α-epoxide
        '3J_H5_H6_beta': 3.0,        # β-epoxide
        '3J_H4_H5': 2.0,             # Modified by epoxide
        '3J_H6_H7': 2.5,
    },
    
    # Fluorinated steroids (9α-fluoro like dexamethasone)
    'fluoro_steroid': {
        '3J_H9_F9': 0.0,             # No H9 when F present
        '2J_H11_F9': 48.0,           # Large geminal-like
        '3J_H8_F9': 24.0,            # Vicinal H-F
        '3J_H12_F9': 5.0,            # Through ring
    },
}


def _detect_steroid_junction(mol: Chem.Mol, c1_idx: int, c2_idx: int) -> Optional[str]:
    """
    Automatically detect steroid junction type between two carbons.
    Returns 'AB_junction', 'BC_junction', 'CD_junction' or None.
    """
    c1 = mol.GetAtomWithIdx(c1_idx)
    c2 = mol.GetAtomWithIdx(c2_idx)
    
    if not (c1.IsInRing() and c2.IsInRing()):
        return None

    # Aromatic carbons are not steroid junctions - they have different coupling values
    if c1.GetIsAromatic() or c2.GetIsAromatic():
        return None

    # Detect if both carbons are in the same fused ring system
    ring_info = mol.GetRingInfo()
    c1_rings = set()
    c2_rings = set()
    
    for i, ring in enumerate(ring_info.AtomRings()):
        if c1_idx in ring:
            c1_rings.add(i)
        if c2_idx in ring:
            c2_rings.add(i)
    
    # If they share a ring or are in adjacent rings
    if c1_rings.intersection(c2_rings) or len(c1_rings) + len(c2_rings) >= 3:
        # Simple approximation: return AB_junction by default
        # A complete implementation would analyze exact topology
        return 'AB_junction'
    
    return None


def get_steroid_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    ring_system: Optional[str] = None,
    coupling_type: Optional[str] = None
) -> Optional[float]:
    """
    Calculate steroid-specific couplings.
    
    IMPORTANT: This function only returns values for protons
    confirmed to be in the steroid fused ring system.
    
    Args:
        ring_system: 'AB_junction', 'BC_junction', 'CD_junction', etc.
        coupling_type: Specific key in table
    """
    # Verify that both protons are attached to carbons in RINGS
    h1 = mol.GetAtomWithIdx(h1_idx)
    h2 = mol.GetAtomWithIdx(h2_idx)
    
    c1_neighbors = [n for n in h1.GetNeighbors() if n.GetSymbol() == 'C']
    c2_neighbors = [n for n in h2.GetNeighbors() if n.GetSymbol() == 'C']
    
    if not c1_neighbors or not c2_neighbors:
        return None
    
    c1 = c1_neighbors[0]
    c2 = c2_neighbors[0]

    # Aromatic carbons are not steroid junctions - use aromatic coupling tables instead
    if c1.GetIsAromatic() or c2.GetIsAromatic():
        return None

    # Both carbons must be in a ring to apply steroid couplings
    if not c1.IsInRing() or not c2.IsInRing():
        return None

    # Automatic detection if not specified
    if ring_system is None:
        ring_system = _detect_steroid_junction(mol, c1.GetIdx(), c2.GetIdx())
        if ring_system is None:
            return None
    
    # Default coupling if coupling_type not specified
    if coupling_type is None:
        if ring_system == 'AB_junction':
            coupling_type = '3J_H5_H6ax_trans'
        elif ring_system == 'BC_junction':
            coupling_type = '3J_H8_H9'
        elif ring_system == 'CD_junction':
            coupling_type = '3J_H14_H15ax'
        else:
            coupling_type = '3J_H5_H6ax_trans'  # Default
    
    if ring_system in _STEROID_COUPLING_TABLE:
        return _STEROID_COUPLING_TABLE[ring_system].get(coupling_type)
    
    return None


# =============================================================================
# SECTION 9: TERPENES
# =============================================================================

_TERPENE_COUPLING_TABLE: Dict[str, Dict[str, float]] = {
    # Cyclic monoterpenes
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
    
    # Sesquiterpenes
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
    
    # Gem-dimethyl
    'gem_dimethyl': {
        '4J_CH3_CH3': 0.0,           # No coupling
        '4J_CH3_CH2': 0.5,           # Weak
    },
}


def get_terpene_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    skeleton_type: str
) -> Optional[float]:
    """
    Calculate couplings in terpenes.
    
    IMPORTANT: This function only returns values for protons
    confirmed to be in a terpenic ring (cyclohexane/cyclopentane).
    For simple alkyl chains, returns None to let the main simulator
    calculate normal couplings.
    """
    # Verify that both protons are attached to carbons in a RING
    # to avoid applying terpene couplings to alkyl chains
    h1 = mol.GetAtomWithIdx(h1_idx)
    h2 = mol.GetAtomWithIdx(h2_idx)
    
    # Get parent carbons
    c1_neighbors = [n for n in h1.GetNeighbors() if n.GetSymbol() == 'C']
    c2_neighbors = [n for n in h2.GetNeighbors() if n.GetSymbol() == 'C']
    
    if not c1_neighbors or not c2_neighbors:
        return None
    
    c1 = c1_neighbors[0]
    c2 = c2_neighbors[0]
    
    # Allow vinyl couplings even outside ring (terpenic side chains)
    is_terpenic_relevant = (
        (c1.IsInRing() or c2.IsInRing()) or  # At least one in ring
        (c1.GetIsAromatic() == False and c2.GetIsAromatic() == False and  # Not aromatic
         (str(c1.GetHybridization()) == 'SP2' or str(c2.GetHybridization()) == 'SP2'))  # At least one sp2 (vinyl)
    )
    
    if not is_terpenic_relevant:
        return None
    
    # Both carbons must be in the SAME ring
    ring_info = mol.GetRingInfo()
    c1_rings = set()
    c2_rings = set()
    
    for i, ring in enumerate(ring_info.AtomRings()):
        if c1.GetIdx() in ring:
            c1_rings.add(i)
        if c2.GetIdx() in ring:
            c2_rings.add(i)
    
    # Common rings OR adjacent rings (junctions) for terpene couplings
    common_rings = c1_rings.intersection(c2_rings)
    adjacent_rings = len(c1_rings) > 0 and len(c2_rings) > 0 and not common_rings
    
    if not (common_rings or adjacent_rings):
        return None
    
    # OK, both protons are in the same ring - apply the table
    if skeleton_type in _TERPENE_COUPLING_TABLE:
        table = _TERPENE_COUPLING_TABLE[skeleton_type]
        return table.get('3J_ax_ax', 7.0)
    
    return None


# =============================================================================
# SECTION 10: MACROCYCLES
# =============================================================================

_MACROCYCLE_COUPLING_TABLE: Dict[str, Dict[str, float]] = {
    # Macrolactones (lactones >8 members)
    'macrolactone': {
        '3J_average': 7.0,           # Partial rotation
        '3J_rigid': 10.0,            # If rigid
        '3J_flexible': 5.0,          # Very flexible
    },
    
    # Macrolides (erythromycin type)
    'macrolide': {
        '3J_backbone': 6.5,
        '3J_sugar_aglycone': 7.5,
        '2J_CH2': -12.0,
    },
    
    # Macrocyclic lactams
    'macrolactam': {
        '3J_NH_CH': 6.0,
        '3J_backbone': 7.0,
    },
    
    # Crown ethers
    'crown_ether': {
        '3J_OCH2CH2O_gauche': 4.0,
        '3J_OCH2CH2O_anti': 10.0,
        '3J_average': 5.5,           # Rotation
    },
    
    # Cyclodextrins
    'cyclodextrin': {
        # Similar to sugars
        '3J_H1_H2': 3.5,
        '3J_H2_H3': 9.5,
        '3J_H3_H4': 9.0,
        '3J_H4_H5': 9.5,
    },
    
    # Cyclic peptides
    'cyclic_peptide': {
        '3J_NH_alphaH': 7.5,         # More rigid
        '3J_alpha_beta': 6.5,
    },
    
    # Depsipeptides
    'depsipeptide': {
        '3J_ester_OCH': 3.5,
        '3J_amide_NCH': 7.0,
    },
    
    # Macrocyclic polyenes
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
    Calculate couplings in macrocycles.
    
    IMPORTANT: This function only returns values for protons
    confirmed to be in the macrocycle.
    
    Args:
        ring_size: Size of macrocycle (>8)
        macrocycle_type: Type of macrocycle
    """
    # Verify that both protons are attached to carbons in the MACROCYCLE
    h1 = mol.GetAtomWithIdx(h1_idx)
    h2 = mol.GetAtomWithIdx(h2_idx)
    
    c1_neighbors = [n for n in h1.GetNeighbors() if n.GetSymbol() == 'C']
    c2_neighbors = [n for n in h2.GetNeighbors() if n.GetSymbol() == 'C']
    
    if not c1_neighbors or not c2_neighbors:
        return None
    
    c1 = c1_neighbors[0]
    c2 = c2_neighbors[0]

    # Aromatic carbons use aromatic coupling tables, not macrocycle tables
    if c1.GetIsAromatic() or c2.GetIsAromatic():
        return None

    # Both carbons must be in a ring (and in the same large ring)
    if not c1.IsInRing() or not c2.IsInRing():
        return None

    # Verify they are in a ring of size >= 8 (macrocycle)
    ring_info = mol.GetRingInfo()
    in_macrocycle = False
    
    for ring in ring_info.AtomRings():
        if len(ring) >= 8 and c1.GetIdx() in ring and c2.GetIdx() in ring:
            in_macrocycle = True
            break
    
    if not in_macrocycle:
        return None
    
    # Adjustment based on ring size
    # Larger rings = more flexible = more averaging = smaller J values
    flexibility_factor = 1.0
    
    if ring_size <= 10:
        flexibility_factor = 1.1   # More rigid, larger J
    elif ring_size <= 14:
        flexibility_factor = 1.0   # Reference
    elif ring_size <= 20:
        flexibility_factor = 0.9   # More flexible, smaller J
    else:
        flexibility_factor = 0.8   # Very flexible, significant averaging
    
    if macrocycle_type in _MACROCYCLE_COUPLING_TABLE:
        base_J = _MACROCYCLE_COUPLING_TABLE[macrocycle_type].get('3J_average', 7.0)
        return base_J * flexibility_factor
    
    # Default value for macrocycles
    return 7.0 * flexibility_factor


# =============================================================================
# SECTION 11: NUCLEOSIDES AND NUCLEOTIDES
# =============================================================================

_NUCLEOTIDE_COUPLING_TABLE: Dict[str, Dict[str, float]] = {
    # Ribose (furanose)
    'ribose': {
        # Depends on pseudorotation
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
    
    # Deoxyribose
    'deoxyribose': {
        '3J_H1_H2a': 6.5,
        '3J_H1_H2b': 7.0,
        '2J_H2a_H2b': -14.0,
        '3J_H2a_H3': 6.0,
        '3J_H2b_H3': 3.5,
        '3J_H3_H4': 3.5,
    },
    
    # Phosphate couplings
    'phosphate': {
        '3J_P_H3': 8.0,              # 31P-H3'
        '3J_P_H5a': 5.0,
        '3J_P_H5b': 3.0,
        '4J_P_H4': 2.0,
    },
    
    # Nitrogenous bases
    'adenine': {
        '3J_H2_H8': 0.0,             # No coupling (separated)
    },
    
    'guanine': {
        '1J_N1H': 90.0,              # If detectable
    },
    
    'cytosine': {
        '3J_H5_H6': 7.5,             # Vinyl coupling
    },
    
    'thymine': {
        '4J_CH3_H6': 1.2,            # Allylic coupling
    },
    
    'uracil': {
        '3J_H5_H6': 8.0,
    },
    
    # Glycosidic coupling H1'-base
    'glycosidic': {
        # Depends on chi angle
        '4J_H1_H8_syn': 1.0,         # Syn conformation
        '4J_H1_H8_anti': 0.3,        # Anti conformation
    },
}

# =============================================================================
# EXTENDED GLYCOSIDIC COUPLINGS (4J and 5J long-range through glycosidic bond)
# =============================================================================
# References:
# - Tvaroska I., Carbohydr. Res. 206, 55 (1990)
# - Bock K., Pedersen C., Adv. Carbohydr. Chem. Biochem. 41, 27 (1983)
# =============================================================================

_GLYCOSIDIC_LONG_RANGE_TABLE: Dict[str, Dict[str, float]] = {
    # 4J couplings through glycosidic bond (H-C-O-C-H path)
    '1_2_linkage': {
        '4J_H1_H2_alpha': 0.8,       # α-(1→2) linkage
        '4J_H1_H2_beta': 0.5,        # β-(1→2) linkage
        '5J_H1_H3': 0.3,             # 5-bond coupling (rare)
    },
    '1_3_linkage': {
        '4J_H1_H3_alpha': 0.7,
        '4J_H1_H3_beta': 0.4,
        '5J_H1_H2': 0.2,
        '5J_H1_H4': 0.2,
    },
    '1_4_linkage': {
        '4J_H1_H4_alpha': 0.9,       # Common in maltose, cellobiose
        '4J_H1_H4_beta': 0.6,
        '5J_H1_H3': 0.3,
        '5J_H1_H5': 0.2,
    },
    '1_6_linkage': {
        '4J_H1_H6a_alpha': 0.5,
        '4J_H1_H6b_alpha': 0.5,
        '4J_H1_H6a_beta': 0.3,
        '4J_H1_H6b_beta': 0.3,
        '5J_H1_H5': 0.2,
    },
}

# =============================================================================
# NUCLEOSIDE CHI ANGLE (glycosidic torsion base-sugar)
# =============================================================================
# References:
# - Davies D.B., Prog. NMR Spectrosc. 12, 135 (1978)
# - Wijmenga S.S., van Buuren B.N.M., Prog. NMR Spectrosc. 32, 287 (1998)
# =============================================================================

_CHI_ANGLE_COUPLINGS: Dict[str, Dict[str, Any]] = {
    # Purines (adenine, guanine) - chi = O4'-C1'-N9-C4
    'purine': {
        'chi_anti': (-180, -90),      # Anti conformation (most common)
        'chi_syn': (0, 90),           # Syn conformation (rare for purines)
        # 4J couplings H1'-H8 (purine H8)
        '4J_H1_H8': {
            'anti_high': 0.5,         # χ ~ -150° (typical anti)
            'anti_low': 0.2,          # χ ~ -90° (low anti)
            'syn': 1.5,               # χ ~ +60° (syn, rare)
        },
        # Cross-ring couplings
        '5J_H2_H8': 0.3,              # Through imidazole ring
    },
    
    # Pyrimidines (cytosine, thymine, uracil) - chi = O4'-C1'-N1-C2
    'pyrimidine': {
        'chi_anti': (-180, -90),
        'chi_syn': (0, 90),
        # 4J couplings H1'-H6 (pyrimidine H6)
        '4J_H1_H6': {
            'anti': 0.3,
            'syn': 1.2,
        },
        # H5-H6 coupling within base (vinyl-like)
        '3J_H5_H6': 7.5,              # Standard vinyl coupling in base
    },
    
    # Modified bases
    '8-oxo-guanine': {
        '4J_H1_H8': 0.0,              # No H8 in 8-oxo-G
        'notes': 'H8 replaced by O8'
    },
    
    '5-methyl-cytosine': {
        '4J_H1_H6': {
            'anti': 0.3,
            'syn': 1.0,
        },
        '4J_H6_CH3': 1.2,             # Allylic coupling to methyl
    },
}


def get_glycosidic_long_range_coupling(
    linkage_type: str,
    coupling_key: str,
    is_alpha: bool = True
) -> Optional[float]:
    """
    Get long-range 4J/5J coupling through glycosidic bond.
    
    Args:
        linkage_type: '1_2', '1_3', '1_4', or '1_6'
        coupling_key: Specific coupling (e.g., 'H1_H4')
        is_alpha: True for α-linkage, False for β
    
    Returns:
        Coupling in Hz or None
    """
    table_key = f"{linkage_type}_linkage"
    if table_key not in _GLYCOSIDIC_LONG_RANGE_TABLE:
        return None
    
    table = _GLYCOSIDIC_LONG_RANGE_TABLE[table_key]
    
    # Build full coupling key
    suffix = '_alpha' if is_alpha else '_beta'
    full_key = f"4J_{coupling_key}{suffix}"
    
    if full_key in table:
        return table[full_key]
    
    # Try without alpha/beta
    simple_key = f"4J_{coupling_key}"
    if simple_key in table:
        return table[simple_key]
    
    # Try 5J
    five_j_key = f"5J_{coupling_key}"
    if five_j_key in table:
        return table[five_j_key]
    
    return None


def get_nucleoside_chi_coupling(
    base_type: str,
    chi_angle_deg: float = None,
    conformation: str = 'anti'
) -> Dict[str, float]:
    """
    Get couplings in nucleoside based on chi angle (glycosidic torsion).
    
    Args:
        base_type: 'purine' or 'pyrimidine' (or specific base name)
        chi_angle_deg: Chi angle in degrees (optional)
        conformation: 'anti' or 'syn' (used if chi_angle not given)
    
    Returns:
        Dictionary of couplings
    """
    # Determine base class
    if base_type.lower() in ['adenine', 'guanine', 'purine', 'a', 'g']:
        params = _CHI_ANGLE_COUPLINGS.get('purine', {})
        coupling_key = '4J_H1_H8'
    elif base_type.lower() in ['cytosine', 'thymine', 'uracil', 'pyrimidine', 'c', 't', 'u']:
        params = _CHI_ANGLE_COUPLINGS.get('pyrimidine', {})
        coupling_key = '4J_H1_H6'
    else:
        return {}
    
    result = {}
    
    # Determine conformation from chi angle if provided
    if chi_angle_deg is not None:
        chi_ranges = params.get('chi_anti', (-180, -90))
        if chi_ranges[0] <= chi_angle_deg <= chi_ranges[1]:
            conformation = 'anti'
        else:
            conformation = 'syn'
    
    # Get coupling value
    coupling_table = params.get(coupling_key, {})
    if isinstance(coupling_table, dict):
        if conformation == 'anti':
            result[coupling_key] = coupling_table.get('anti', coupling_table.get('anti_high', 0.3))
        else:
            result[coupling_key] = coupling_table.get('syn', 1.0)
    
    # Add base internal couplings
    if '3J_H5_H6' in params:
        result['3J_H5_H6'] = params['3J_H5_H6']
    
    return result


# Pseudorotation equations
def _calculate_ribose_coupling(P: float, coupling_type: str) -> float:
    """
    Calculate ribose coupling according to pseudorotation phase P.
    
    P = 0 deg -> C3'-endo (N-type, RNA)
    P = 180 deg -> C2'-endo (S-type, DNA)
    
    Equation: J = A + B*cos(P + phase) + C*cos(2*(P + phase))
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
    """Calculate couplings in nucleosides/nucleotides."""
    if sugar_type in _NUCLEOTIDE_COUPLING_TABLE:
        table = _NUCLEOTIDE_COUPLING_TABLE[sugar_type]
        
        # If pseudorotation phase is known, calculate precisely
        if pseudorotation_phase is not None and sugar_type == 'ribose':
            # Determine which coupling
            return _calculate_ribose_coupling(pseudorotation_phase, '3J_H1_H2')
        
        # Otherwise use average value
        return table.get('3J_H1_H2_C2endo', 5.0)
    
    return None


# =============================================================================
# SECTION 12: AUTOMATIC DETECTION AND DISPATCH
# =============================================================================

class MoleculeClassifier:
    """Classify molecules to choose the right coupling table."""
    
    @staticmethod
    def classify(mol: Chem.Mol) -> List[str]:
        """
        Return the structure types present in the molecule.
        
        Returns:
            List of types: ['alkyne', 'aldehyde', 'peptide', ...]
        """
        types = []
        
        # Traverse atoms
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            
            # Phosphorus
            if symbol == 'P':
                types.append('phosphorus')
            
            # Boron
            elif symbol == 'B':
                types.append('borane')
            
            # Silicon
            elif symbol == 'Si':
                types.append('silane')
            
            # Metals
            elif symbol in ['Fe', 'Ru', 'Os', 'Co', 'Rh', 'Ir', 'Ni', 'Pd', 'Pt',
                           'Cu', 'Ag', 'Au', 'Zn', 'Ti', 'Zr', 'Hf', 'Mo', 'W']:
                types.append('organometallic')
        
        # Structural patterns via SMARTS
        smarts_patterns = {
            'alkyne': '[C]#[C]',
            'allene': '[C]=[C]=[C]',
            'aldehyde': '[CH]=O',
            'amide': '[NX3][CX3](=[OX1])',
            'ester': '[OX2][CX3](=[OX1])',
            'macrolactone': '[OR1](=[O])',  # Approximate
        }
        
        for name, smarts in smarts_patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                types.append(name)
        
        # Ring detection
        ring_info = mol.GetRingInfo()
        ring_sizes = [len(r) for r in ring_info.AtomRings()]
        
        if any(s > 8 for s in ring_sizes):
            types.append('macrocycle')
        
        if 'amide' in types and len(types) > 3:
            # Probably a peptide
            types.append('peptide')
        
        # Steroid: 4 fused rings (6-6-6-5 or similar) with stereocenters
        if len(ring_sizes) >= 4:
            fused_count = sum(1 for s in ring_sizes if s in [5, 6])
            if fused_count >= 4:
                types.append('steroid')
        
        # ✅ NEW: Détection améliorée des stéroïdes par stéréocentres
        # Un stéroïde typique a 4+ cycles fusionnés ET 5+ stéréocentres
        if len(ring_sizes) >= 3:
            n_stereocenters = sum(1 for a in mol.GetAtoms() 
                                 if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED)
            if n_stereocenters >= 4 and len(ring_sizes) >= 3:
                if 'steroid' not in types:
                    types.append('steroid')
        
        # Terpene: isoprene patterns OR gem-dimethyl on cyclohexane
        isoprene = Chem.MolFromSmarts('[CH3]C(=C)C')
        if isoprene and mol.HasSubstructMatch(isoprene):
            types.append('terpene')
        
        # ✅ NEW: Détection terpènes par gem-diméthyl (ex: pinène)
        gem_dimethyl = Chem.MolFromSmarts('[CH3]C([CH3])([#6])[#6]')  # C(CH3)2 lié à 2 carbones
        if gem_dimethyl and mol.HasSubstructMatch(gem_dimethyl):
            if 'terpene' not in types:
                types.append('terpene')
        
        # ✅ NEW: Détection terpènes par cyclohexène substitué
        cyclohexene_methyl = Chem.MolFromSmarts('[CH3]C1=CCCCC1')  # Limonène-like
        if cyclohexene_methyl and mol.HasSubstructMatch(cyclohexene_methyl):
            if 'terpene' not in types:
                types.append('terpene')
        
        # ✅ NEW: Monoterpènes type menthane (cyclohexane + isopropyle + méthyle)
        # Menthol pattern: cyclohexane avec CH(CH3)2 et CH3
        # Note: [!a] excludes aromatic carbons to avoid false positives (e.g., cumene)
        isopropyl_on_ring = Chem.MolFromSmarts('[CH;R;!a]([CH3])([CH3])')  # isopropyle lié au cycle non-aromatique
        if isopropyl_on_ring and mol.HasSubstructMatch(isopropyl_on_ring):
            # Vérifier qu'il y a aussi un méthyle sur le cycle
            methyl_on_ring = Chem.MolFromSmarts('[CH3][C;R1;!a]')  # CH3 lié à un C du cycle non-aromatique
            if methyl_on_ring and mol.HasSubstructMatch(methyl_on_ring):
                if 'terpene' not in types:
                    types.append('terpene')

        # ✅ NEW: Détection simple - isopropyle attaché à un cyclohexane
        # Note: [!a] excludes aromatic carbons to avoid false positives (e.g., cumene)
        isopropyl_cyclohexane = Chem.MolFromSmarts('CC(C)[CH;R1;!a]')  # isopropyle sur cycle non-aromatique
        if isopropyl_cyclohexane and mol.HasSubstructMatch(isopropyl_cyclohexane):
            if 'terpene' not in types:
                types.append('terpene')
        
        # ✅ NEW: Détection bicyclo[3.1.1] (pinane) ou [2.2.1] (norbornane) pour terpènes
        for ring in ring_info.AtomRings():
            ring_size = len(ring)
            if ring_size in (4, 5) and len(ring_sizes) >= 2:
                # Petit cycle fusionné → possible système terpénique bicyclique
                if any(s == 6 for s in ring_sizes):
                    if 'terpene' not in types:
                        types.append('terpene')
                    break
        
        # Nucleoside
        nucleobase = Chem.MolFromSmarts('[nR1]1[cR1][nR1][cR1][cR1]1')  # Purine-like
        if nucleobase and mol.HasSubstructMatch(nucleobase):
            types.append('nucleoside')
        
        # ✅ NEW: Pyrimidine base pour nucléosides
        pyrimidine = Chem.MolFromSmarts('[nR1]1[cR1][cR1][nR1][cR1][cR1]1')
        if pyrimidine and mol.HasSubstructMatch(pyrimidine):
            if 'nucleoside' not in types:
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
    Main dispatch function to specialized tables.
    
    Args:
        mol: RDKit molecule
        h1_idx, h2_idx: Proton indices
        distance: Topological distance (number of bonds)
        bond_path: Atom path between h1 and h2
    
    Returns:
        J in Hz or None if not applicable
    """
    # ✅ NEW: Input validation to prevent out-of-bounds errors
    if mol is None or not isinstance(mol, Chem.Mol):
        return None
    
    num_atoms = mol.GetNumAtoms()
    if num_atoms == 0:
        return None
    
    # Check that h1_idx and h2_idx are within bounds
    if not (0 <= h1_idx < num_atoms) or not (0 <= h2_idx < num_atoms):
        return None  # Return None instead of raising to be robust

    # =========================================================================
    # ✅ CRITICAL FIX: Filter out long-range couplings (⁴J and beyond)
    # =========================================================================
    # Only ²J (geminal) and ³J (vicinal) have significant magnitude (2-15 Hz)
    # ⁴J (W-coupling) and ⁵J are typically < 2 Hz and context-specific
    # Returning None here prevents tables from assigning 11-12 Hz to ⁴J/⁵J
    # =========================================================================
    if distance > 3:
        return None  # Long-range couplings handled elsewhere (typically ~0-2 Hz)

    # Classify molecule (with cache for optimization)
    if not hasattr(mol, '_coupling_types_cache'):
        mol._coupling_types_cache = MoleculeClassifier.classify(mol)
    mol_types = mol._coupling_types_cache
    
    # Try each specialized table
    
    # 1. Alkynes
    if 'alkyne' in mol_types and bond_path:
        J = get_alkyne_coupling(mol, h1_idx, h2_idx, bond_path)
        if J is not None:
            return J
    
    # 2. Aldehydes (with high resolution support)
    if 'aldehyde' in mol_types:
        J = get_aldehyde_coupling(mol, h1_idx, h2_idx, distance, high_resolution=True)
        if J is not None:
            return J
    
    # 3. Peptides/Amides
    if 'peptide' in mol_types or 'amide' in mol_types:
        J = get_peptide_coupling(mol, h1_idx, h2_idx)
        if J is not None:
            return J
    
    # 4-6. Heteronuclear elements (optimized)
    hetero_elements = {'phosphorus': 'P', 'borane': 'B', 'silane': 'Si'}
    hetero_functions = {
        'P': get_phosphine_coupling,
        'B': get_borane_coupling, 
        'Si': get_silane_coupling
    }
    
    for mol_type, symbol in hetero_elements.items():
        if mol_type in mol_types:
            # Pre-calculate indices to avoid O(n^2)
            hetero_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == symbol]
            for hetero_idx in hetero_indices:
                J = hetero_functions[symbol](mol, h1_idx, hetero_idx, distance)
                if J is not None:
                    return J
    
    # 7. Steroids (with automatic detection)
    if 'steroid' in mol_types:
        J = get_steroid_coupling(mol, h1_idx, h2_idx)  # Automatic detection
        if J is not None:
            return J
    
    # 8. Terpenes
    if 'terpene' in mol_types:
        J = get_terpene_coupling(mol, h1_idx, h2_idx, 'menthane')
        if J is not None:
            return J
    
    # 9. Macrocycles
    if 'macrocycle' in mol_types:
        # Find the size of the largest ring
        ring_info = mol.GetRingInfo()
        max_ring = max(len(r) for r in ring_info.AtomRings()) if ring_info.AtomRings() else 0
        
        J = get_macrocycle_coupling(mol, h1_idx, h2_idx, max_ring)
        if J is not None:
            return J
    
    # 10. Nucleosides
    if 'nucleoside' in mol_types:
        J = get_nucleotide_coupling(mol, h1_idx, h2_idx, 'ribose')
        if J is not None:
            return J
    
    # No specialized table applicable
    return None


# =============================================================================
# SECTION 13: TESTS
# =============================================================================

def validate_critical_couplings() -> bool:
    """
    Validates critical coupling values against literature references.
    
    Returns True if all validations pass, False otherwise.
    """
    passed = True
    errors = []
    
    # Aldehyde 3J validation (JACS 1963 reference, acetaldehyde ~2.8 Hz)
    aldehyde_3J = _ALDEHYDE_COUPLING_TABLE['3J_CHO_alphaH']['default']
    if not (2.0 <= aldehyde_3J <= 3.5):
        errors.append(f"Aldehyde 3J out of range: {aldehyde_3J} Hz (expected 2.0-3.5)")
        passed = False
    
    # Aldehyde 4J validation
    aldehyde_4J = _ALDEHYDE_COUPLING_TABLE['4J_CHO_betaH']['trans_planar']
    if not (0.8 <= aldehyde_4J <= 1.5):
        errors.append(f"Aldehyde 4J out of range: {aldehyde_4J} Hz (expected 0.8-1.5)")
        passed = False
    
    # Borane 2J validation
    borane_2J = _BORANE_COUPLING_TABLE['2J_BCH']['default']
    if not (2.0 <= borane_2J <= 5.0):
        errors.append(f"Borane 2J out of range: {borane_2J} Hz (expected 2.0-5.0)")
        passed = False
    
    # Silane 2J validation
    silane_2J = _SILANE_COUPLING_TABLE['2J_SiCH']['default']
    if not (3.0 <= silane_2J <= 5.5):
        errors.append(f"Silane 2J out of range: {silane_2J} Hz (expected 3.0-5.5)")
        passed = False
    
    # Alkyne 3J_across_triple validation
    if '3J_across_triple' not in _ALKYNE_COUPLING_TABLE:
        errors.append("Missing 3J_across_triple in alkyne table")
        passed = False
    else:
        alkyne_3J = _ALKYNE_COUPLING_TABLE['3J_across_triple']['default']
        if not (8.0 <= alkyne_3J <= 12.0):
            errors.append(f"Alkyne 3J_across_triple out of range: {alkyne_3J} Hz (expected 8-12)")
            passed = False
    
    # Fluorosilane validation
    if 'fluorosilane' not in _SILANE_COUPLING_TABLE:
        errors.append("Missing fluorosilane section in silane table")
        passed = False
    
    # Polyhedral boranes validation
    if 'polyhedral' not in _BORANE_COUPLING_TABLE:
        errors.append("Missing polyhedral section in borane table")
        passed = False
    
    # Print results
    if passed:
        print("[PASS] All critical coupling validations passed")
    else:
        print("[FAIL] Critical coupling validation errors:")
        for error in errors:
            print(f"  - {error}")
    
    return passed


def _test_coupling_tables():
    """Test coupling tables."""
    print("=" * 60)
    print("COMPLETE COUPLING TABLES TEST")
    print("=" * 60)
    
    # First run critical validations
    print("\n--- Critical Values Validation ---")
    validate_critical_couplings()
    
    test_cases = [
        ("Propyne (alkyne)", "CC#C"),
        ("Allene", "C=C=C"),
        ("Benzaldehyde", "O=Cc1ccccc1"),
        ("Propanal", "CCC=O"),
        ("N-methylacetamide", "CC(=O)NC"),
        ("Alanine", "C[C@H](N)C(=O)O"),
        ("Trimethylphosphine", "CP(C)C"),
        ("Trimethylsilane", "C[SiH](C)C"),
        ("Pinene", "CC1=CCC2CC1C2(C)C"),
        ("12-membered Macrolide", "C1CCCCCCCCCCC(=O)O1"),
        ("Adenosine", "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"),
    ]
    
    for name, smiles in test_cases:
        print(f"\n{name}")
        print("-" * 40)
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print("  ERROR: Invalid SMILES")
            continue
        
        # Classify
        types = MoleculeClassifier.classify(mol)
        print(f"  Detected types: {types}")
        
        # Test dispatch
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        
        # Find 2 protons
        h_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'H']
        
        if len(h_indices) >= 2:
            J = get_specialized_coupling(mol, h_indices[0], h_indices[1], 3)
            print(f"  J test: {J} Hz")
    
    print("\n" + "=" * 60)
    print("TESTS COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    _test_coupling_tables()