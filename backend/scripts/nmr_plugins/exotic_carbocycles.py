# nmr_plugins/exotic_carbocycles.py
# Unusual carbocycle plugin: annulenes, tropylium, and corannulene

"""Handle nonconventional carbocycles.

Covers annulenes with n ≥ 10; tropylium and cycloheptatrienyl ions;
pentalenes and heptalenes; corannulene and partial fullerenes;
strained cyclobutene/cyclopropene rings; and bridged bicyclic systems."""

from typing import Dict, List, Tuple, Optional, Any, Set
import numpy as np
import math

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors
    from rdkit.Chem.rdmolops import GetFormalCharge
    RDKIT_AVAILABLE = True
except ImportError:
    print("Warning: RDKit non disponible pour le plugin exotic_carbocycles")
    RDKIT_AVAILABLE = False


# =============================================================================
# Unusual carbocycle coupling tables
# =============================================================================

# [n]Annulenes
ANNULENE_COUPLINGS: Dict[int, Dict[str, float]] = {
    10: {  # [10]Annulene
        '3J_vicinal': 11.0,      # Couplage vicinal
        '4J_allylic': 2.5,       # Through the pi system
        '5J_homoallylic': 0.8,   # Long range
        'ring_current_effect': 1.2,  # Ring-current effect
    },
    14: {  # [14]Annulene
        '3J_vicinal': 10.5,
        '4J_allylic': 3.0,       # Greater delocalization
        '5J_homoallylic': 1.2,
        'ring_current_effect': 1.5,
    },
    18: {  # [18]Annulene (aromatic)
        '3J_vicinal': 9.8,       # Slightly reduced
        '4J_allylic': 3.5,       # Strong delocalization
        '5J_homoallylic': 1.5,
        'ring_current_effect': 2.0,  # Strong ring current
    },
    22: {  # [22]Annulene
        '3J_vicinal': 10.2,
        '4J_allylic': 3.2,
        '5J_homoallylic': 1.3,
        'ring_current_effect': 1.7,
    },
}

# Tropylium and cycloheptatriene
TROPYLIUM_COUPLINGS: Dict[str, float] = {
    'tropylium_adjacent': 8.5,      # C⁺-adjacent
    'tropylium_meta': 1.8,          # Meta to the cation
    'tropylium_para': 0.5,          # Para (faible)
    'cycloheptatriene_vicinal': 12.0,  # Neutral cycloheptatriene
    'cycloheptatriene_allylic': 4.0,   # Through the pi system
    'norcaradiene_vicinal': 8.0,    # Norcaradiene form
    'norcaradiene_bridge': 2.5,     # Through the bridge
}

# Pentalenes and heptalenes (antiaromatic)
PENTALENE_COUPLINGS: Dict[str, float] = {
    'pentalene_vicinal': 13.5,      # Increased (antiaromatic)
    'pentalene_across_bridge': 3.5, # Through the central bridge
    'pentalene_long_range': 1.8,    # W couplings
    'heptalene_vicinal': 12.8,      # Heptalene
    'heptalene_across': 2.8,        # Through the system
}

# Corannulene and partial fullerenes
CORANNULENE_COUPLINGS: Dict[str, float] = {
    'rim_vicinal': 9.2,             # Bowl rim
    'rim_meta': 1.5,                # Meta on the rim
    'spoke_hub': 8.8,               # Rayon vers centre (spoke)
    'hub_vicinal': 8.5,             # Centre (hub)
    'across_bowl': 0.8,             # Across the curved structure
    'bowl_inversion': 2.0,          # Bowl inversion
}

# Cycles tendus
STRAINED_CYCLES: Dict[str, float] = {
    'cyclopropene_vicinal': 15.0,   # Highly strained
    'cyclopropene_allylic': 5.5,    # Through the pi system
    'cyclobutene_vicinal': 13.2,    # Tendu
    'cyclobutene_across': 2.8,      # Trans-cyclique
    'bicyclobutane_bridge': 6.0,    # Pont cyclobutane
    'spiropentane_vicinal': 11.5,   # Spiropentane
}

# Unusual bridged bicyclic systems
BRIDGED_BICYCLES: Dict[str, float] = {
    'norbornene_endo_exo': 4.2,     # Endo vs exo
    'norbornene_bridge': 2.8,       # Through the bridge
    'norbornadiene_vicinal': 8.5,   # Norbornadiene
    'quadricyclane_cage': 3.5,      # Cage quadricyclane
    'adamantane_cage': 1.2,         # Adamantane (rigide)
    'cubane_edge': 8.8,             # Cube edge
    'cubane_diagonal': 0.5,         # Diagonale cube
}

# Recognition SMARTS patterns
EXOTIC_PATTERNS: Dict[str, str] = {
    # Approximate annulene patterns
    'cyclooctatetraene': '[CH1]=[CH1][CH1]=[CH1][CH1]=[CH1][CH1]=[CH1]1',
    'cyclodecapentaene': 'C1=CC=CC=CC=CC=C1',  # [10]Annulene
    
    # Tropylium
    'tropylium': '[c+]1[c][c][c][c][c][c]1',
    'cycloheptatriene': 'C1=CC=CC=CC1',
    
    # Cycles tendus
    'cyclopropene': 'C1=CC1',
    'cyclobutene': 'C1=CCC1',
    'bicyclobutane': 'C12CC1C2',
    
    # Bridged systems
    'norbornene': 'C12C=CC(C1)C2',
    'norbornadiene': 'C12C=CC(C1)C=C2',
    'adamantane': 'C12CC3CC(C1)CC(C2)C3',
    'cubane': 'C12C3C4C1C5C4C3C25',  # Simplified
    
    # Approximate corannulene pattern
    'corannulene_fragment': '[c]1[c]2[c][c][c][c]2[c]2[c][c][c][c]2[c]1',
}


# =============================================================================
# Detection functions
# =============================================================================

def detect_annulene(mol: Chem.Mol, atom_idx: int) -> Optional[Dict[str, Any]]:
    """Detect annulenes and their ring sizes. Return a dictionary describing the detected annulene."""
    atom = mol.GetAtomWithIdx(atom_idx)
    
    if not atom.GetIsAromatic() and atom.GetSymbol() != 'C':
        return None
    
    # Find all rings containing this atom.
    ring_info = mol.GetRingInfo()
    atom_rings = [ring for ring in ring_info.AtomRings() if atom_idx in ring]
    
    for ring in atom_rings:
        ring_size = len(ring)
        
        # Check for a potential annulene.
        if ring_size >= 8 and ring_size % 2 == 0:  # Pairs ≥ 8
            
            # Check alternating single and double bonds.
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            all_carbon = all(a.GetSymbol() == 'C' for a in ring_atoms)
            
            if all_carbon:
                # Count double bonds in the ring.
                double_bonds = 0
                for i in range(ring_size):
                    atom1_idx = ring[i]
                    atom2_idx = ring[(i + 1) % ring_size]
                    bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
                    if bond and bond.GetBondTypeAsDouble() == 2.0:
                        double_bonds += 1
                
                # Annulene criterion: approximately 50% double bonds.
                if double_bonds >= ring_size // 3:
                    
                    # Determine aromaticity using the Hückel rule.
                    # 4n+2 pi electrons indicate aromaticity.
                    pi_electrons = double_bonds * 2
                    is_aromatic = (pi_electrons - 2) % 4 == 0
                    
                    return {
                        'type': f'[{ring_size}]annulene',
                        'ring_atoms': ring,
                        'ring_size': ring_size,
                        'double_bonds': double_bonds,
                        'pi_electrons': pi_electrons,
                        'is_aromatic': is_aromatic,
                        'coupling_data': ANNULENE_COUPLINGS.get(ring_size, {}),
                    }
    
    return None


def detect_tropylium_system(mol: Chem.Mol, atom_idx: int) -> Optional[Dict[str, Any]]:
    """Detect tropylium/cycloheptatriene systems and return a dictionary describing the match."""
    atom = mol.GetAtomWithIdx(atom_idx)
    
    # Look for seven-carbon rings.
    ring_info = mol.GetRingInfo()
    atom_rings = [ring for ring in ring_info.AtomRings() if atom_idx in ring]
    
    for ring in atom_rings:
        if len(ring) == 7:  # Seven-membered ring
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            
            # Check whether every atom is carbon.
            if all(a.GetSymbol() == 'C' for a in ring_atoms):
                
                # Chercher charge positive (tropylium)
                total_charge = sum(GetFormalCharge(a) for a in ring_atoms)
                
                if total_charge > 0:
                    return {
                        'type': 'tropylium',
                        'ring_atoms': ring,
                        'formal_charge': total_charge,
                        'aromatic': True,  # Six pi electrons
                        'coupling_type': 'tropylium',
                    }
                
                else:
                    # Count double bonds for neutral cycloheptatriene.
                    double_bonds = 0
                    for i in range(7):
                        atom1_idx = ring[i]
                        atom2_idx = ring[(i + 1) % 7]
                        bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
                        if bond and bond.GetBondTypeAsDouble() == 2.0:
                            double_bonds += 1
                    
                    if double_bonds == 3:  # Cycloheptatriene
                        return {
                            'type': 'cycloheptatriene',
                            'ring_atoms': ring,
                            'double_bonds': double_bonds,
                            'aromatic': False,  # 6 π + 2 non-liants
                            'coupling_type': 'cycloheptatriene',
                        }
    
    return None


def detect_pentalene_system(mol: Chem.Mol, atom_idx: int) -> Optional[Dict[str, Any]]:
    """Detect antiaromatic pentalene/heptalene systems and return a dictionary describing the match."""
    atom = mol.GetAtomWithIdx(atom_idx)
    
    if atom.GetSymbol() != 'C':
        return None
    
    # Look for fused bicyclic systems.
    ring_info = mol.GetRingInfo()
    atom_rings = [ring for ring in ring_info.AtomRings() if atom_idx in ring]
    
    # Pentalene: two fused five-membered rings.
    five_rings = [ring for ring in atom_rings if len(ring) == 5]
    
    if len(five_rings) >= 2:
        # Check for shared ring atoms.
        for i, ring1 in enumerate(five_rings):
            for ring2 in five_rings[i+1:]:
                common_atoms = set(ring1) & set(ring2)
                
                if len(common_atoms) == 2:  # Edge-fused rings
                    # Check for conjugated carbon.
                    all_atoms = set(ring1) | set(ring2)
                    carbons = [mol.GetAtomWithIdx(i) for i in all_atoms 
                              if mol.GetAtomWithIdx(i).GetSymbol() == 'C']
                    
                    if len(carbons) == len(all_atoms):  # All carbon atoms
                        return {
                            'type': 'pentalene',
                            'ring1': ring1,
                            'ring2': ring2,
                            'fusion_atoms': list(common_atoms),
                            'antiaromatic': True,  # Eight pi electrons
                            'coupling_enhancement': 1.15,  # Antiaromatique
                        }
    
    # Heptalene candidate rings: 5+7 or 7+7.
    seven_rings = [ring for ring in atom_rings if len(ring) == 7]
    
    if five_rings and seven_rings:
        for ring5 in five_rings:
            for ring7 in seven_rings:
                common_atoms = set(ring5) & set(ring7)
                if len(common_atoms) == 2:
                    return {
                        'type': 'heptalene_5_7',
                        'ring5': ring5,
                        'ring7': ring7,
                        'fusion_atoms': list(common_atoms),
                        'antiaromatic': True,  # Twelve pi electrons
                        'coupling_enhancement': 1.08,
                    }
    
    return None


def detect_corannulene_fragment(mol: Chem.Mol, atom_idx: int) -> Optional[Dict[str, Any]]:
    """Detect corannulene/fullerene fragments and return a dictionary describing the curved structure."""
    atom = mol.GetAtomWithIdx(atom_idx)
    
    if not atom.GetIsAromatic():
        return None
    
    # Count fused aromatic rings around this atom.
    ring_info = mol.GetRingInfo()
    atom_rings = [ring for ring in ring_info.AtomRings() 
                  if atom_idx in ring and len(ring) <= 6]
    
    # Typical corannulene: five six-membered rings around a central five-membered ring.
    if len(atom_rings) >= 3:  # Multiple fused rings
        
        ring_sizes = [len(ring) for ring in atom_rings]
        five_rings = ring_sizes.count(5)
        six_rings = ring_sizes.count(6)
        
        # Corannulene/fullerene pattern
        if five_rings >= 1 and six_rings >= 2:
            
            # Approximate curvature using ring-fusion density.
            fusion_density = len(atom_rings) / 6.0  # Normalized
            
            if fusion_density > 0.5:  # Highly fused
                return {
                    'type': 'corannulene_fragment',
                    'center_atom': atom_idx,
                    'rings': atom_rings,
                    'five_rings': five_rings,
                    'six_rings': six_rings,
                    'fusion_density': fusion_density,
                    'bowl_shaped': True,
                    'coupling_type': 'curved_aromatic',
                }
    
    return None


def detect_strained_cycle(mol: Chem.Mol, atom_idx: int) -> Optional[Dict[str, Any]]:
    """Detect highly strained unsaturated three- and four-membered rings. Return a dictionary describing the ring strain."""
    atom = mol.GetAtomWithIdx(atom_idx)
    
    ring_info = mol.GetRingInfo()
    atom_rings = [ring for ring in ring_info.AtomRings() if atom_idx in ring]
    
    for ring in atom_rings:
        ring_size = len(ring)
        
        if ring_size in [3, 4]:  # Highly strained rings
            
            # Compter insaturations
            double_bonds = 0
            for i in range(ring_size):
                atom1_idx = ring[i]
                atom2_idx = ring[(i + 1) % ring_size]
                bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
                if bond and bond.GetBondTypeAsDouble() == 2.0:
                    double_bonds += 1
            
            if double_bonds >= 1:  # Unsaturated and strained
                
                if ring_size == 3 and double_bonds == 1:
                    return {
                        'type': 'cyclopropene',
                        'ring_atoms': ring,
                        'strain_energy': 'very_high',
                        'coupling_enhancement': 1.25,
                        'bond_angles': 60,  # Highly distorted
                    }
                
                elif ring_size == 4 and double_bonds == 1:
                    return {
                        'type': 'cyclobutene',
                        'ring_atoms': ring,
                        'strain_energy': 'high',
                        'coupling_enhancement': 1.12,
                        'bond_angles': 90,  # Distorted
                    }
    
    return None


# =============================================================================
# FONCTION PRINCIPALE
# =============================================================================

def get_coupling(mol: Chem.Mol, h1_idx: int, h2_idx: int) -> Optional[float]:
    """Calculate couplings in unusual carbocycles.

    mol is an RDKit molecule; h1_idx and h2_idx are hydrogen indices.
    Return J in Hz for a recognized carbocycle, otherwise None."""
    if not RDKIT_AVAILABLE:
        return None
        
    try:
        atom1 = mol.GetAtomWithIdx(h1_idx)
        atom2 = mol.GetAtomWithIdx(h2_idx)
        
        # Only process hydrogen atoms.
        if atom1.GetSymbol() != 'H' or atom2.GetSymbol() != 'H':
            return None
        
        # Find the atoms bearing the hydrogens.
        parent1_idx = None
        parent2_idx = None
        
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            if idx != h1_idx and idx != h2_idx:
                neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
                if h1_idx in neighbors:
                    parent1_idx = idx
                if h2_idx in neighbors:
                    parent2_idx = idx
        
        if parent1_idx is None or parent2_idx is None:
            return None
        
        # Calculer distance entre porteurs
        try:
            path = Chem.GetShortestPath(mol, parent1_idx, parent2_idx)
            if not path:
                return None
            
            distance = len(path) - 1
        except Exception as e:
            # Handle RDKit invariant errors.
            if "Invariant Violation" in str(e):
                return None  # Ignorer silencieusement
            else:
                raise e
        
        # Detect unusual ring systems.
        
        # Annulenes
        annulene_info = detect_annulene(mol, parent1_idx) or detect_annulene(mol, parent2_idx)
        if annulene_info:
            ring_size = annulene_info['ring_size']
            coupling_data = annulene_info['coupling_data']
            
            if distance == 2:  # Vicinal
                base_J = coupling_data.get('3J_vicinal', 10.0)
                ring_effect = coupling_data.get('ring_current_effect', 1.0)
                return base_J * ring_effect
            elif distance == 3:  # Allylic
                return coupling_data.get('4J_allylic', 2.0)
            elif distance == 4:  # Homoallylic  
                return coupling_data.get('5J_homoallylic', 0.8)
        
        # Strained rings
        strain_info = detect_strained_cycle(mol, parent1_idx) or detect_strained_cycle(mol, parent2_idx)
        if strain_info:
            strain_type = strain_info['type']
            enhancement = strain_info['coupling_enhancement']
            
            if distance == 2:
                if strain_type == 'cyclopropene':
                    return STRAINED_CYCLES['cyclopropene_vicinal'] * enhancement
                elif strain_type == 'cyclobutene':
                    return STRAINED_CYCLES['cyclobutene_vicinal'] * enhancement
        
        # Approximate other unusual rings.
        # Check for an unsaturated ring with at least eight members.
        ring_info = mol.GetRingInfo()
        parent1_rings = [ring for ring in ring_info.AtomRings() if parent1_idx in ring]
        parent2_rings = [ring for ring in ring_info.AtomRings() if parent2_idx in ring]
        
        for ring in parent1_rings:
            if len(ring) >= 8:  # Grand cycle
                if distance == 2:
                    return 10.5  # Modified vicinal coupling
        
        return None
        
    except Exception as e:
        # Handle detection errors.
        print(f"Erreur dans exotic_carbocycles: {e}")
        return None


# =============================================================================
# UTILITAIRES
# =============================================================================

def calculate_ring_strain_effect(ring_size: int, bond_angles: List[float]) -> float:
    """Calculate the ring-strain correction factor for couplings.

    ring_size is the number of ring atoms and bond_angles lists the bond angles.
    Return a multiplicative coupling factor."""
    # Ideal angles
    ideal_angles = {3: 60, 4: 90, 5: 108, 6: 120}
    ideal = ideal_angles.get(ring_size, 120)
    
    # Mean deviation from the ideal angle
    if bond_angles:
        avg_deviation = np.mean([abs(angle - ideal) for angle in bond_angles])
        # Greater distortion produces a larger coupling correction.
        enhancement = 1.0 + (avg_deviation / 100.0)  # Approximation
        return min(enhancement, 2.0)  # Cap the factor at two.
    
    return 1.0


def identify_aromatic_system_type(mol: Chem.Mol, ring_atoms: List[int]) -> str:
    """Classify the system using the Hückel rule. Return 'aromatic', 'antiaromatic', or 'nonaromatic'."""
    # Count pi electrons.
    pi_electrons = 0
    
    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        
        # Atom's pi electrons
        if atom.GetIsAromatic():
            if atom.GetSymbol() == 'C':
                pi_electrons += 1  # One pi electron per aromatic carbon
            elif atom.GetSymbol() in ['N', 'O']:
                # Possible lone pair
                pi_electrons += 2 if atom.GetTotalNumHs() == 0 else 1
    
    # Hückel rule: 4n+2 is aromatic; 4n is antiaromatic.
    if pi_electrons % 4 == 2:
        return 'aromatic'
    elif pi_electrons % 4 == 0 and pi_electrons > 0:
        return 'antiaromatic'
    else:
        return 'nonaromatic'


# =============================================================================
# TESTS
# =============================================================================

def _test_exotic_carbocycles():
    """Test the exotic_carbocycles plugin."""
    print("=" * 60)
    print("TEST DU PLUGIN CARBOCYCLES EXOTIQUES")
    print("=" * 60)
    
    test_cases = [
        # (nom, smiles, J_attendu, tolerance)
        ("Cyclooctatetraène", "C1=CC=CC=CC=C1", 11.0, 2.0),
        ("Tropylium", "[CH+]=1C=CC=CC=C1", 8.5, 1.5),
        ("Cyclopropène", "C1=CC1", 15.0, 3.0),
        ("Cycloheptatriène", "C1=CC=CC=CC1", 12.0, 2.0),
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
        
        # Find hydrogens.
        h_atoms = [i for i, atom in enumerate(mol.GetAtoms()) 
                  if atom.GetSymbol() == 'H']
        
        if len(h_atoms) >= 2:
            J = get_coupling(mol, h_atoms[0], h_atoms[1])
            
            if J is not None:
                diff = abs(J - expected_J)
                if diff <= tolerance:
                    print(f"  PASS: J = {J:.1f} Hz (expected: {expected_J} ± {tolerance})")
                    passed += 1
                else:
                    print(f"  FAIL: J = {J:.1f} Hz (difference: {diff:.1f})")
            else:
                print("  WARN: No unusual carbocycle coupling detected")
        else:
            print("  WARN: Not enough hydrogens")
    
    print(f"\nResult: {passed}/{len(test_cases)} tests passed")
    return passed, len(test_cases)


if __name__ == "__main__":
    _test_exotic_carbocycles()