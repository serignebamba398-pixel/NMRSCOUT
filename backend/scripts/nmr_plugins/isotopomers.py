# nmr_plugins/isotopomers.py
# Isotopomer coupling plugin (¹³C, ¹⁵N, ²H)

"""Handle couplings involving isotopomers.

Covers ¹³C satellites (0.55% natural abundance per satellite), enriched ¹⁵N
for biomolecular NMR, deuterium isotope effects, ³¹P-¹H in nucleotides and
phospholipids, and ¹⁹F-¹H in fluorinated drugs."""

from typing import Dict, List, Tuple, Optional, Any
import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
except ImportError:
    raise ImportError("RDKit requis pour le plugin isotopomers")


# =============================================================================
# Isotopomer coupling tables
# =============================================================================

# Direct ¹³C-¹H couplings
C13_H_DIRECT_COUPLINGS: Dict[str, float] = {
    '1J_C13_H_sp3': 125.0,      # CH₃, CH₂, CH aliphatiques
    '1J_C13_H_sp2_alkene': 160.0,  # =CH vinyliques
    '1J_C13_H_sp2_arom': 160.0,    # Aromatiques
    '1J_C13_H_sp': 250.0,          # ≡CH alcynes
    '1J_C13_H_aldehyde': 175.0,    # CHO
    '1J_C13_H_formate': 200.0,     # HCOO-
}

# Indirect ¹³C-¹H couplings through bonds
C13_H_INDIRECT_COUPLINGS: Dict[str, float] = {
    '2J_C13_C_H_gem': 5.0,         # ¹³C-C-H (geminal)
    '3J_C13_C_C_H_vicinal': 4.0,   # ¹³C-C-C-H (vicinal)
    '4J_C13_long_range': 1.0,      # Long-range ⁴J
    '2J_C13_O_H': 3.0,             # ¹³C-O-H (alcools)
    '2J_C13_N_H': 2.0,             # ¹³C-N-H (amines)
}

# ¹⁵N-¹H couplings in biomolecular NMR
N15_H_COUPLINGS: Dict[str, float] = {
    '1J_N15_H_amide': 90.0,        # CONH
    '1J_N15_H_amine_primary': 70.0,   # NH₂ aliphatique
    '1J_N15_H_amine_secondary': 75.0, # NHR aliphatique
    '1J_N15_H_imine': 95.0,        # C=NH
    '1J_N15_H_pyrrole': 96.0,      # Pyrrole NH
    '1J_N15_H_imidazole': 92.0,    # Imidazole NH
    '1J_N15_H_indole': 94.0,       # Indole NH
    '2J_N15_C_H': 11.0,            # ¹⁵N-C-H
    '3J_N15_C_C_H': 7.0,           # ¹⁵N-C-C-H
}

# Deuterium isotope effects
DEUTERIUM_EFFECTS: Dict[str, float] = {
    '1J_C_D': 20.0,                # ¹J(C-D) = ¹J(C-H)/6.5
    '2J_H_D_gem': 1.5,             # H-C-D geminal
    '3J_H_D_vicinal': 0.8,         # H-C-C-D vicinal
    'alpha_isotope_shift': 0.15,   # Δδ per alpha D (ppm)
    'beta_isotope_shift': 0.02,    # Δδ per beta D (ppm)
    'gamma_isotope_shift': 0.005,  # Δδ per gamma D (ppm)
}

# ³¹P-¹H couplings in nucleotides and phospholipids
P31_H_COUPLINGS: Dict[str, float] = {
    '2J_P_C_H': 15.0,              # ³¹P-C-H direct
    '3J_P_O_C_H': 8.0,             # ³¹P-O-C-H (phosphoesters)
    '3J_P_C_C_H': 12.0,            # ³¹P-C-C-H
    '4J_P_O_C_C_H': 2.0,           # Long-range ⁴J
    '3J_P_N_C_H': 6.0,             # Phosphoramides
}

# ¹⁹F-¹H couplings in fluorinated drugs
F19_H_COUPLINGS: Dict[str, float] = {
    '2J_F_C_H_gem': 45.0,          # F-C-H geminal
    '3J_F_C_C_H_vicinal': 25.0,    # F-C-C-H vicinal
    '4J_F_aromatic': 8.0,          # ⁴J aromatique F...H
    '5J_F_meta': 3.0,              # ⁵J meta aromatique
    '3J_F_vinyl_cis': 20.0,        # F-C=C-H cis
    '3J_F_vinyl_trans': 35.0,      # F-C=C-H trans
}

# ¹³C satellite patterns for identification
C13_SATELLITE_PATTERNS: Dict[str, Dict[str, Any]] = {
    'methyl_CH3': {
        'intensity_ratio': 0.0055,     # 0.55% on each side
        'separation_hz': 125.0,        # ¹J(¹³C-H) = 125 Hz
        'multiplicity': 'doublet',     # ¹³C coupled to 3H: doublet
    },
    'methylene_CH2': {
        'intensity_ratio': 0.0055,
        'separation_hz': 125.0,
        'multiplicity': 'doublet',     # ¹³C coupled to 2H: doublet
    },
    'methine_CH': {
        'intensity_ratio': 0.0055,
        'separation_hz': 125.0,
        'multiplicity': 'doublet',     # ¹³C coupled to 1H: doublet
    },
    'aromatic_CH': {
        'intensity_ratio': 0.0055,
        'separation_hz': 160.0,        # Aromatique
        'multiplicity': 'doublet',
    },
    'alkene_CH': {
        'intensity_ratio': 0.0055,
        'separation_hz': 160.0,        # Vinylique
        'multiplicity': 'doublet',
    },
}


# =============================================================================
# Detection functions
# =============================================================================

def detect_isotopomer_type(mol: Chem.Mol, atom_idx: int) -> Optional[str]:
    """Detect an enriched isotope at atom_idx in an RDKit molecule. Return the isotope type or None."""
    atom = mol.GetAtomWithIdx(atom_idx)
    
    # Check the explicitly assigned isotope.
    isotope = atom.GetIsotope()
    symbol = atom.GetSymbol()
    
    if isotope != 0:  # Isotope explicite
        if symbol == 'C' and isotope == 13:
            return 'C13'
        elif symbol == 'N' and isotope == 15:
            return 'N15'
        elif symbol == 'H' and isotope == 2:
            return 'D'  # Deuterium
        elif symbol == 'P' and isotope == 31:
            return 'P31'
        elif symbol == 'F' and isotope == 19:
            return 'F19'
    
    # Detect naturally abundant isotopes by element symbol.
    elif symbol == 'P':  # ³¹P naturel
        return 'P31'
    elif symbol == 'F':  # ¹⁹F naturel
        return 'F19'
    
    return None


def detect_c13_satellite_context(mol: Chem.Mol, c_idx: int) -> Optional[str]:
    """Classify carbon for ¹³C satellites. Return 'methyl', 'methylene', 'methine', 'aromatic', 'alkene', or 'quaternary'."""
    carbon = mol.GetAtomWithIdx(c_idx)
    
    if carbon.GetSymbol() != 'C':
        return None
    
    # Count hydrogens.
    h_count = carbon.GetTotalNumHs()
    
    # Determine hybridization.
    if carbon.GetIsAromatic():
        if h_count == 1:
            return 'aromatic_CH'
        elif h_count == 0:
            return 'aromatic_quaternary'
    else:
        # sp² ou sp³ ?
        neighbors = carbon.GetNeighbors()
        double_bonds = sum(1 for n in neighbors 
                          if mol.GetBondBetweenAtoms(c_idx, n.GetIdx()).GetBondTypeAsDouble() == 2)
        
        if double_bonds > 0:  # sp² ou sp
            if h_count == 1:
                return 'alkene_CH'
            elif h_count == 2:
                return 'alkene_CH2'
        else:  # sp³
            if h_count == 3:
                return 'methyl_CH3'
            elif h_count == 2:
                return 'methylene_CH2'
            elif h_count == 1:
                return 'methine_CH'
            elif h_count == 0:
                return 'quaternary'
    
    return None


def get_phosphorus_coupling_context(mol: Chem.Mol, p_idx: int, h_idx: int) -> Optional[str]:
    """Classify the P-H coupling environment, returning a bond-path type such as 'P_C_H', 'P_O_C_H', or 'P_N_C_H'."""
    # Find the P...H path.
    path = Chem.GetShortestPath(mol, p_idx, h_idx)
    if not path or len(path) > 5:  # Max ⁴J
        return None
    
    path_atoms = [mol.GetAtomWithIdx(i) for i in path]
    path_symbols = [a.GetSymbol() for a in path_atoms]
    
    if len(path) == 3:  # ²J: P-X-H
        middle_atom = path_symbols[1]
        if middle_atom == 'C':
            return 'P_C_H'
        elif middle_atom == 'O':
            return 'P_O_H'
        elif middle_atom == 'N':
            return 'P_N_H'
    
    elif len(path) == 4:  # ³J: P-X-Y-H
        if 'O' in path_symbols[1:3]:
            return 'P_O_C_H'  # Phosphoester
        elif 'N' in path_symbols[1:3]:
            return 'P_N_C_H'  # Phosphoramide
        else:
            return 'P_C_C_H'  # Direct
    
    elif len(path) == 5:  # ⁴J: P-X-Y-Z-H
        return 'P_long_range'
    
    return None


# =============================================================================
# Main coupling function
# =============================================================================

def get_coupling(mol: Chem.Mol, h1_idx: int, h2_idx: int) -> Optional[float]:
    """Calculate isotopomer couplings in an RDKit molecule.

    h1_idx and h2_idx identify hydrogen or isotope atoms.
    Return J in Hz, or None if no isotopomer coupling is detected."""
    atom1 = mol.GetAtomWithIdx(h1_idx)
    atom2 = mol.GetAtomWithIdx(h2_idx)
    
    # Require at least one isotope of interest.
    isotope1 = detect_isotopomer_type(mol, h1_idx)
    isotope2 = detect_isotopomer_type(mol, h2_idx)
    
    if not isotope1 and not isotope2:
        return None
    
    # Calculer distance
    path = Chem.GetShortestPath(mol, h1_idx, h2_idx)
    if not path:
        return None
    
    distance = len(path) - 1
    
    # ¹³C-¹H couplings
    if isotope1 == 'C13' and atom2.GetSymbol() == 'H':
        if distance == 1:  # ¹J direct
            context = detect_c13_satellite_context(mol, h1_idx)
            if context:
                if 'aromatic' in context or 'alkene' in context:
                    return C13_H_DIRECT_COUPLINGS['1J_C13_H_sp2_arom']
                else:
                    return C13_H_DIRECT_COUPLINGS['1J_C13_H_sp3']
        elif distance == 2:  # ²J
            return C13_H_INDIRECT_COUPLINGS['2J_C13_C_H_gem']
        elif distance == 3:  # ³J
            return C13_H_INDIRECT_COUPLINGS['3J_C13_C_C_H_vicinal']
    
    # ¹⁵N-¹H couplings
    elif isotope1 == 'N15' and atom2.GetSymbol() == 'H':
        if distance == 1:  # ¹J direct N-H
            # Determine the nitrogen environment.
            nitrogen = mol.GetAtomWithIdx(h1_idx)
            
            # Amide ?
            for neighbor in nitrogen.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    for nn in neighbor.GetNeighbors():
                        if nn.GetSymbol() == 'O':
                            bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nn.GetIdx())
                            if bond and bond.GetBondTypeAsDouble() == 2:
                                return N15_H_COUPLINGS['1J_N15_H_amide']
            
            # Default to amine.
            h_count = nitrogen.GetTotalNumHs()
            if h_count >= 2:
                return N15_H_COUPLINGS['1J_N15_H_amine_primary']
            else:
                return N15_H_COUPLINGS['1J_N15_H_amine_secondary']
        
        elif distance == 2:  # ²J
            return N15_H_COUPLINGS['2J_N15_C_H']
        elif distance == 3:  # ³J
            return N15_H_COUPLINGS['3J_N15_C_C_H']
    
    # ³¹P-¹H couplings
    elif isotope1 == 'P31' and atom2.GetSymbol() == 'H':
        context = get_phosphorus_coupling_context(mol, h1_idx, h2_idx)
        if context and context in P31_H_COUPLINGS:
            return P31_H_COUPLINGS[context]
    
    # ¹⁹F-¹H couplings
    elif isotope1 == 'F19' and atom2.GetSymbol() == 'H':
        if distance == 2:  # ²J geminal
            return F19_H_COUPLINGS['2J_F_C_H_gem']
        elif distance == 3:  # ³J vicinal
            # Check for a vinylic environment.
            path_atoms = [mol.GetAtomWithIdx(i) for i in path]
            has_double_bond = any(
                mol.GetBondBetweenAtoms(path[i], path[i+1]).GetBondTypeAsDouble() == 2
                for i in range(len(path)-1)
            )
            if has_double_bond:
                return F19_H_COUPLINGS['3J_F_vinyl_trans']  # Approximation
            else:
                return F19_H_COUPLINGS['3J_F_C_C_H_vicinal']
        elif distance == 4:  # ⁴J
            return F19_H_COUPLINGS['4J_F_aromatic']
    
    # Deuterium effects
    elif isotope1 == 'D' and atom2.GetSymbol() == 'H':
        if distance == 2:  # Geminal H-C-D
            return DEUTERIUM_EFFECTS['2J_H_D_gem']
        elif distance == 3:  # Vicinal H-C-C-D
            return DEUTERIUM_EFFECTS['3J_H_D_vicinal']
    
    return None


# =============================================================================
# UTILITAIRES
# =============================================================================

def get_c13_satellite_info(separation_hz: float) -> Dict[str, Any]:
    """Identify observed ¹³C satellites from separation_hz. Return a dictionary containing the carbon type and predicted multiplicity."""
    # Determine the type from the separation.
    if 120 <= separation_hz <= 130:
        return {
            'carbon_type': 'sp3',
            'J_value': separation_hz,
            'intensity_percent': 1.1,  # 2 × 0.55%
            'multiplicity': 'doublet'
        }
    elif 155 <= separation_hz <= 165:
        return {
            'carbon_type': 'sp2_aromatic_or_alkene',
            'J_value': separation_hz,
            'intensity_percent': 1.1,
            'multiplicity': 'doublet'
        }
    elif separation_hz >= 240:
        return {
            'carbon_type': 'sp_alkyne',
            'J_value': separation_hz,
            'intensity_percent': 1.1,
            'multiplicity': 'doublet'
        }
    else:
        return {
            'carbon_type': 'unknown',
            'J_value': separation_hz,
            'intensity_percent': 1.1,
            'multiplicity': 'doublet'
        }


def estimate_deuterium_shift_effect(mol: Chem.Mol, d_idx: int, h_idx: int) -> float:
    """Estimate the deuterium isotope shift in ppm; negative values indicate shielding."""
    path = Chem.GetShortestPath(mol, d_idx, h_idx)
    if not path:
        return 0.0
    
    distance = len(path) - 1
    
    if distance == 1:  # α
        return -DEUTERIUM_EFFECTS['alpha_isotope_shift']
    elif distance == 2:  # β
        return -DEUTERIUM_EFFECTS['beta_isotope_shift']
    elif distance == 3:  # γ
        return -DEUTERIUM_EFFECTS['gamma_isotope_shift']
    
    return 0.0


# =============================================================================
# TESTS
# =============================================================================

def _test_isotopomers():
    """Test the isotopomers plugin."""
    print("=" * 60)
    print("TEST DU PLUGIN ISOTOPOMÈRES")
    print("=" * 60)
    
    test_cases = [
        # (name, isotope-labeled SMILES, expected_J, tolerance)
        ("¹³C-Méthanol", "[13CH3]O", 125.0, 5.0),
        ("¹⁵N-Méthylamine", "[15NH2]C", 70.0, 5.0),
        ("Deutéro-éthanol", "[2H]CCO", 1.5, 0.5),
        ("³¹P-Phosphate", "P(=O)(O)(O)CC", 8.0, 2.0),
    ]
    
    passed = 0
    
    for name, smiles, expected_J, tolerance in test_cases:
        print(f"\n{name}: {smiles}")
        print("-" * 40)
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print("  FAIL: Invalid SMILES")
            continue
        
        mol = Chem.AddHs(mol)
        
        # Find atoms of interest.
        isotope_atoms = []
        h_atoms = []
        
        for i, atom in enumerate(mol.GetAtoms()):
            if detect_isotopomer_type(mol, i):
                isotope_atoms.append(i)
            elif atom.GetSymbol() == 'H':
                h_atoms.append(i)
        
        if isotope_atoms and h_atoms:
            J = get_coupling(mol, isotope_atoms[0], h_atoms[0])
            
            if J is not None:
                diff = abs(J - expected_J)
                if diff <= tolerance:
                    print(f"  PASS: J = {J:.1f} Hz (expected: {expected_J} ± {tolerance})")
                    passed += 1
                else:
                    print(f"  FAIL: J = {J:.1f} Hz (difference: {diff:.1f})")
            else:
                print("  WARN: No coupling detected")
        else:
            print("  WARN: No isotope-labeled atoms found")
    
    print(f"\nResult: {passed}/{len(test_cases)} tests passed")
    return passed, len(test_cases)


if __name__ == "__main__":
    _test_isotopomers()