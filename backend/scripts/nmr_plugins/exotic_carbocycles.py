# nmr_plugins/exotic_carbocycles.py
# Plugin pour carbocycles exotiques (annulènes, tropylium, corannulène)

"""
Plugin pour les carbocycles non-conventionnels.

Couvre:
- Annulènes ([n]annulènes, n ≥ 10)
- Tropylium et ions cycloheptatriényle  
- Pentalènes et heptalènes
- Corannulène et fullerènes partiels
- Cycles tendus (cyclobutène, cyclopropène)
- Systèmes bicycliques pontés
"""

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
# TABLES DE COUPLAGES EXOTIQUES
# =============================================================================

# Annulènes [n]annulènes
ANNULENE_COUPLINGS: Dict[int, Dict[str, float]] = {
    10: {  # [10]annulène
        '3J_vicinal': 11.0,      # Couplage vicinal
        '4J_allylic': 2.5,       # À travers π
        '5J_homoallylic': 0.8,   # Long range
        'ring_current_effect': 1.2,  # Effet de courant de cycle
    },
    14: {  # [14]annulène
        '3J_vicinal': 10.5,
        '4J_allylic': 3.0,       # Plus de délocalisation
        '5J_homoallylic': 1.2,
        'ring_current_effect': 1.5,
    },
    18: {  # [18]annulène (aromatique)
        '3J_vicinal': 9.8,       # Légèrement réduit
        '4J_allylic': 3.5,       # Forte délocalisation
        '5J_homoallylic': 1.5,
        'ring_current_effect': 2.0,  # Fort courant de cycle
    },
    22: {  # [22]annulène
        '3J_vicinal': 10.2,
        '4J_allylic': 3.2,
        '5J_homoallylic': 1.3,
        'ring_current_effect': 1.7,
    },
}

# Tropylium et cycloheptatriène
TROPYLIUM_COUPLINGS: Dict[str, float] = {
    'tropylium_adjacent': 8.5,      # C⁺-adjacent
    'tropylium_meta': 1.8,          # Méta au cation
    'tropylium_para': 0.5,          # Para (faible)
    'cycloheptatriene_vicinal': 12.0,  # Cycloheptatriène neutre
    'cycloheptatriene_allylic': 4.0,   # À travers π
    'norcaradiene_vicinal': 8.0,    # Forme norcaradiène
    'norcaradiene_bridge': 2.5,     # À travers pont
}

# Pentalènes et heptalènes (antiaromatiques)
PENTALENE_COUPLINGS: Dict[str, float] = {
    'pentalene_vicinal': 13.5,      # Augmenté (antiaromatique)
    'pentalene_across_bridge': 3.5, # À travers pont central
    'pentalene_long_range': 1.8,    # Couplages W
    'heptalene_vicinal': 12.8,      # Heptalène
    'heptalene_across': 2.8,        # À travers système
}

# Corannulène et fullerènes partiels
CORANNULENE_COUPLINGS: Dict[str, float] = {
    'rim_vicinal': 9.2,             # Bord de bol (rim)
    'rim_meta': 1.5,                # Méta sur bord
    'spoke_hub': 8.8,               # Rayon vers centre (spoke)
    'hub_vicinal': 8.5,             # Centre (hub)
    'across_bowl': 0.8,             # À travers courbure
    'bowl_inversion': 2.0,          # Inversion de bol
}

# Cycles tendus
STRAINED_CYCLES: Dict[str, float] = {
    'cyclopropene_vicinal': 15.0,   # Très tendu
    'cyclopropene_allylic': 5.5,    # À travers π
    'cyclobutene_vicinal': 13.2,    # Tendu
    'cyclobutene_across': 2.8,      # Trans-cyclique
    'bicyclobutane_bridge': 6.0,    # Pont cyclobutane
    'spiropentane_vicinal': 11.5,   # Spiropentane
}

# Bicycliques pontés exotiques
BRIDGED_BICYCLES: Dict[str, float] = {
    'norbornene_endo_exo': 4.2,     # Endo vs exo
    'norbornene_bridge': 2.8,       # À travers pont
    'norbornadiene_vicinal': 8.5,   # Norbornadiène
    'quadricyclane_cage': 3.5,      # Cage quadricyclane
    'adamantane_cage': 1.2,         # Adamantane (rigide)
    'cubane_edge': 8.8,             # Arête de cube
    'cubane_diagonal': 0.5,         # Diagonale cube
}

# Patterns SMARTS de reconnaissance
EXOTIC_PATTERNS: Dict[str, str] = {
    # Annulènes (approximatifs)
    'cyclooctatetraene': '[CH1]=[CH1][CH1]=[CH1][CH1]=[CH1][CH1]=[CH1]1',
    'cyclodecapentaene': 'C1=CC=CC=CC=CC=C1',  # [10]annulène
    
    # Tropylium
    'tropylium': '[c+]1[c][c][c][c][c][c]1',
    'cycloheptatriene': 'C1=CC=CC=CC1',
    
    # Cycles tendus
    'cyclopropene': 'C1=CC1',
    'cyclobutene': 'C1=CCC1',
    'bicyclobutane': 'C12CC1C2',
    
    # Systèmes pontés
    'norbornene': 'C12C=CC(C1)C2',
    'norbornadiene': 'C12C=CC(C1)C=C2',
    'adamantane': 'C12CC3CC(C1)CC(C2)C3',
    'cubane': 'C12C3C4C1C5C4C3C25',  # Simplifié
    
    # Corannulène (approximatif)
    'corannulene_fragment': '[c]1[c]2[c][c][c][c]2[c]2[c][c][c][c]2[c]1',
}


# =============================================================================
# FONCTIONS DE DÉTECTION
# =============================================================================

def detect_annulene(mol: Chem.Mol, atom_idx: int) -> Optional[Dict[str, Any]]:
    """
    Détecte les annulènes et détermine leur taille.
    
    Returns:
        Dict avec infos sur l'annulène détecté
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    
    if not atom.GetIsAromatic() and atom.GetSymbol() != 'C':
        return None
    
    # Trouver tous les cycles contenant cet atome
    ring_info = mol.GetRingInfo()
    atom_rings = [ring for ring in ring_info.AtomRings() if atom_idx in ring]
    
    for ring in atom_rings:
        ring_size = len(ring)
        
        # Vérifier si c'est potentiellement un annulène
        if ring_size >= 8 and ring_size % 2 == 0:  # Pairs ≥ 8
            
            # Vérifier alternance simple/double liaisons
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            all_carbon = all(a.GetSymbol() == 'C' for a in ring_atoms)
            
            if all_carbon:
                # Compter liaisons doubles dans le cycle
                double_bonds = 0
                for i in range(ring_size):
                    atom1_idx = ring[i]
                    atom2_idx = ring[(i + 1) % ring_size]
                    bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
                    if bond and bond.GetBondTypeAsDouble() == 2.0:
                        double_bonds += 1
                
                # Critère d'annulène: ~50% de liaisons doubles
                if double_bonds >= ring_size // 3:
                    
                    # Déterminer aromaticité selon règle de Hückel
                    # 4n+2 π-électrons = aromatique
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
    """
    Détecte les systèmes tropylium/cycloheptatriène.
    
    Returns:
        Dict avec infos sur le système détecté
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    
    # Chercher cycles à 7 carbones
    ring_info = mol.GetRingInfo()
    atom_rings = [ring for ring in ring_info.AtomRings() if atom_idx in ring]
    
    for ring in atom_rings:
        if len(ring) == 7:  # Cycle à 7
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            
            # Tous carbones ?
            if all(a.GetSymbol() == 'C' for a in ring_atoms):
                
                # Chercher charge positive (tropylium)
                total_charge = sum(GetFormalCharge(a) for a in ring_atoms)
                
                if total_charge > 0:
                    return {
                        'type': 'tropylium',
                        'ring_atoms': ring,
                        'formal_charge': total_charge,
                        'aromatic': True,  # 6 π-électrons
                        'coupling_type': 'tropylium',
                    }
                
                else:
                    # Cycloheptatriène neutre - compter doubles liaisons
                    double_bonds = 0
                    for i in range(7):
                        atom1_idx = ring[i]
                        atom2_idx = ring[(i + 1) % 7]
                        bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
                        if bond and bond.GetBondTypeAsDouble() == 2.0:
                            double_bonds += 1
                    
                    if double_bonds == 3:  # Cycloheptatriène
                        return {
                            'type': 'cycloheptatriene',
                            'ring_atoms': ring,
                            'double_bonds': double_bonds,
                            'aromatic': False,  # 6 π + 2 non-liants
                            'coupling_type': 'cycloheptatriene',
                        }
    
    return None


def detect_pentalene_system(mol: Chem.Mol, atom_idx: int) -> Optional[Dict[str, Any]]:
    """
    Détecte les systèmes pentalène/heptalène (antiaromatiques).
    
    Returns:
        Dict avec infos sur le système antiaromatique
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    
    if atom.GetSymbol() != 'C':
        return None
    
    # Chercher systèmes bicycliques fusionnés
    ring_info = mol.GetRingInfo()
    atom_rings = [ring for ring in ring_info.AtomRings() if atom_idx in ring]
    
    # Pentalène: deux cycles à 5 fusionnés
    five_rings = [ring for ring in atom_rings if len(ring) == 5]
    
    if len(five_rings) >= 2:
        # Vérifier fusion (atomes communs)
        for i, ring1 in enumerate(five_rings):
            for ring2 in five_rings[i+1:]:
                common_atoms = set(ring1) & set(ring2)
                
                if len(common_atoms) == 2:  # Fusion par arête
                    # Vérifier si c'est du carbone conjugué
                    all_atoms = set(ring1) | set(ring2)
                    carbons = [mol.GetAtomWithIdx(i) for i in all_atoms 
                              if mol.GetAtomWithIdx(i).GetSymbol() == 'C']
                    
                    if len(carbons) == len(all_atoms):  # Tous carbones
                        return {
                            'type': 'pentalene',
                            'ring1': ring1,
                            'ring2': ring2,
                            'fusion_atoms': list(common_atoms),
                            'antiaromatic': True,  # 8 π-électrons
                            'coupling_enhancement': 1.15,  # Antiaromatique
                        }
    
    # Heptalène: cycles 5+7 ou 7+7
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
                        'antiaromatic': True,  # 12 π-électrons
                        'coupling_enhancement': 1.08,
                    }
    
    return None


def detect_corannulene_fragment(mol: Chem.Mol, atom_idx: int) -> Optional[Dict[str, Any]]:
    """
    Détecte les fragments corannulène/fullerène.
    
    Returns:
        Dict avec infos sur la structure courbée
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    
    if not atom.GetIsAromatic():
        return None
    
    # Compter cycles aromatiques fusionnés autour de cet atome
    ring_info = mol.GetRingInfo()
    atom_rings = [ring for ring in ring_info.AtomRings() 
                  if atom_idx in ring and len(ring) <= 6]
    
    # Corannulène typique: 5 cycles à 6 + 1 cycle à 5 central
    if len(atom_rings) >= 3:  # Multi-fusionné
        
        ring_sizes = [len(ring) for ring in atom_rings]
        five_rings = ring_sizes.count(5)
        six_rings = ring_sizes.count(6)
        
        # Pattern corannulène/fullerène
        if five_rings >= 1 and six_rings >= 2:
            
            # Vérifier courbure (approximation par densité de fusion)
            fusion_density = len(atom_rings) / 6.0  # Normalisé
            
            if fusion_density > 0.5:  # Fortement fusionné
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
    """
    Détecte les cycles très tendus (3, 4 membresinsaturés).
    
    Returns:
        Dict avec infos sur la contrainte cyclique
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    
    ring_info = mol.GetRingInfo()
    atom_rings = [ring for ring in ring_info.AtomRings() if atom_idx in ring]
    
    for ring in atom_rings:
        ring_size = len(ring)
        
        if ring_size in [3, 4]:  # Cycles très tendus
            
            # Compter insaturations
            double_bonds = 0
            for i in range(ring_size):
                atom1_idx = ring[i]
                atom2_idx = ring[(i + 1) % ring_size]
                bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
                if bond and bond.GetBondTypeAsDouble() == 2.0:
                    double_bonds += 1
            
            if double_bonds >= 1:  # Insaturé ET tendu
                
                if ring_size == 3 and double_bonds == 1:
                    return {
                        'type': 'cyclopropene',
                        'ring_atoms': ring,
                        'strain_energy': 'very_high',
                        'coupling_enhancement': 1.25,
                        'bond_angles': 60,  # Très déformés
                    }
                
                elif ring_size == 4 and double_bonds == 1:
                    return {
                        'type': 'cyclobutene',
                        'ring_atoms': ring,
                        'strain_energy': 'high',
                        'coupling_enhancement': 1.12,
                        'bond_angles': 90,  # Déformés
                    }
    
    return None


# =============================================================================
# FONCTION PRINCIPALE
# =============================================================================

def get_coupling(mol: Chem.Mol, h1_idx: int, h2_idx: int) -> Optional[float]:
    """
    Calcule les couplages dans les carbocycles exotiques.
    
    Args:
        mol: Molécule RDKit
        h1_idx, h2_idx: Indices des hydrogènes
        
    Returns:
        J en Hz pour carbocycles exotiques ou None
    """
    if not RDKIT_AVAILABLE:
        return None
        
    try:
        atom1 = mol.GetAtomWithIdx(h1_idx)
        atom2 = mol.GetAtomWithIdx(h2_idx)
        
        # Ne traiter que les hydrogènes
        if atom1.GetSymbol() != 'H' or atom2.GetSymbol() != 'H':
            return None
        
        # Trouver les atomes porteurs
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
            # Protection contre les erreurs RDKit invariant
            if "Invariant Violation" in str(e):
                return None  # Ignorer silencieusement
            else:
                raise e
        
        # === DÉTECTION SYSTÈMES EXOTIQUES ===
        
        # 1. Annulènes
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
        
        # 2. Cycles tendus (plus simple et robuste)
        strain_info = detect_strained_cycle(mol, parent1_idx) or detect_strained_cycle(mol, parent2_idx)
        if strain_info:
            strain_type = strain_info['type']
            enhancement = strain_info['coupling_enhancement']
            
            if distance == 2:
                if strain_type == 'cyclopropene':
                    return STRAINED_CYCLES['cyclopropene_vicinal'] * enhancement
                elif strain_type == 'cyclobutene':
                    return STRAINED_CYCLES['cyclobutene_vicinal'] * enhancement
        
        # 3. Approximation pour autres cycles exotiques
        # Si c'est un grand cycle (≥8) avec insaturation
        ring_info = mol.GetRingInfo()
        parent1_rings = [ring for ring in ring_info.AtomRings() if parent1_idx in ring]
        parent2_rings = [ring for ring in ring_info.AtomRings() if parent2_idx in ring]
        
        for ring in parent1_rings:
            if len(ring) >= 8:  # Grand cycle
                if distance == 2:
                    return 10.5  # Couplage vicinal modifié
        
        return None
        
    except Exception as e:
        # Gestion d'erreur gracieuse
        print(f"Erreur dans exotic_carbocycles: {e}")
        return None


# =============================================================================
# UTILITAIRES
# =============================================================================

def calculate_ring_strain_effect(ring_size: int, bond_angles: List[float]) -> float:
    """
    Calcule l'effet de la contrainte cyclique sur les couplages.
    
    Args:
        ring_size: Taille du cycle
        bond_angles: Liste des angles de liaison
        
    Returns:
        Facteur multiplicatif pour les couplages
    """
    # Angles idéaux
    ideal_angles = {3: 60, 4: 90, 5: 108, 6: 120}
    ideal = ideal_angles.get(ring_size, 120)
    
    # Écart moyen à l'idéal
    if bond_angles:
        avg_deviation = np.mean([abs(angle - ideal) for angle in bond_angles])
        # Plus la déformation est grande, plus les couplages sont affectés
        enhancement = 1.0 + (avg_deviation / 100.0)  # Approximation
        return min(enhancement, 2.0)  # Limiter à x2
    
    return 1.0


def identify_aromatic_system_type(mol: Chem.Mol, ring_atoms: List[int]) -> str:
    """
    Identifie le type de système aromatique selon Hückel.
    
    Returns:
        'aromatic', 'antiaromatic', ou 'nonaromatic'
    """
    # Compter π-électrons
    pi_electrons = 0
    
    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        
        # π-électrons de l'atome
        if atom.GetIsAromatic():
            if atom.GetSymbol() == 'C':
                pi_electrons += 1  # Un électron π par carbone aromatique
            elif atom.GetSymbol() in ['N', 'O']:
                # Paire libre éventuelle
                pi_electrons += 2 if atom.GetTotalNumHs() == 0 else 1
    
    # Règle de Hückel: 4n+2 = aromatique, 4n = antiaromatique
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
    """Tests du plugin exotic_carbocycles."""
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
            print("  ❌ SMILES invalide")
            continue
        
        mol = Chem.AddHs(mol)
        
        # Trouver hydrogènes
        h_atoms = [i for i, atom in enumerate(mol.GetAtoms()) 
                  if atom.GetSymbol() == 'H']
        
        if len(h_atoms) >= 2:
            J = get_coupling(mol, h_atoms[0], h_atoms[1])
            
            if J is not None:
                diff = abs(J - expected_J)
                if diff <= tolerance:
                    print(f"  ✅ J = {J:.1f} Hz (attendu: {expected_J} ± {tolerance})")
                    passed += 1
                else:
                    print(f"  ❌ J = {J:.1f} Hz (écart: {diff:.1f})")
            else:
                print("  ⚠️ Pas de couplage exotique détecté")
        else:
            print("  ⚠️ Pas assez d'hydrogènes")
    
    print(f"\n📊 Résultat: {passed}/{len(test_cases)} tests réussis")
    return passed, len(test_cases)


if __name__ == "__main__":
    _test_exotic_carbocycles()