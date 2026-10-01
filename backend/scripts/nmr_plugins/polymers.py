# nmr_plugins/polymers.py
# Polymer and monomer plugin

"""Handle polymers and monomers.

Covers vinyl monomers, polymer repeating units, copolymers,
functionalized polymers, and tacticity effects."""

from typing import Dict, List, Tuple, Optional, Any
import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    raise ImportError("RDKit requis pour le plugin polymers")


# Polymer coupling tables
POLYMER_COUPLINGS: Dict[str, float] = {
    # Vinyl monomers
    'vinyl_trans': 17.0,        # Trans vinylique
    'vinyl_cis': 11.0,          # Cis vinylique
    'vinyl_geminal': 1.5,       # Geminal vinylique
    
    # Polymer chains
    'polymer_backbone': 7.0,     # Polymer backbone
    'side_chain': 6.5,          # Side chain
    
    # Tacticity
    'isotactic_CH2': 13.5,      # Isotactique
    'syndiotactic_CH2': 12.8,   # Syndiotactique  
    'atactic_CH2': 13.0,        # Atactique (moyenne)
    
    # Functionalized polymers
    'ester_side_chain': 4.0,    # Ester chain
    'amide_side_chain': 6.0,    # Amide chain
}

# Patterns SMARTS
POLYMER_PATTERNS: Dict[str, str] = {
    'vinyl_monomer': 'C=C',
    'styrene': 'C=Cc1ccccc1',
    'acrylate': 'C=CC(=O)O',
    'methacrylate': 'C=C(C)C(=O)O',
    'vinyl_alcohol': 'C=CO',
}


def detect_vinyl_monomer(mol: Chem.Mol, atom_idx: int) -> Optional[Dict[str, Any]]:
    """Detect vinyl monomers."""
    atom = mol.GetAtomWithIdx(atom_idx)
    
    if atom.GetSymbol() != 'C':
        return None
    
    # Chercher double liaison C=C
    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() == 'C':
            bond = mol.GetBondBetweenAtoms(atom_idx, neighbor.GetIdx())
            if bond and bond.GetBondTypeAsDouble() == 2.0:
                # C=C double bond
                
                # Determine the monomer type.
                monomer_type = 'simple_alkene'
                
                # Check for styrene.
                for nn in atom.GetNeighbors():
                    if nn.GetIsAromatic():
                        monomer_type = 'styrene'
                        break
                
                # Acrylate ?
                for nn in neighbor.GetNeighbors():
                    if nn.GetSymbol() == 'C':
                        for nnn in nn.GetNeighbors():
                            if nnn.GetSymbol() == 'O':
                                bond_co = mol.GetBondBetweenAtoms(nn.GetIdx(), nnn.GetIdx())
                                if bond_co and bond_co.GetBondTypeAsDouble() == 2.0:
                                    monomer_type = 'acrylate'
                
                return {
                    'type': monomer_type,
                    'vinyl_carbons': [atom_idx, neighbor.GetIdx()],
                    'coupling_type': 'vinyl'
                }
    
    return None


def get_coupling(mol: Chem.Mol, h1_idx: int, h2_idx: int) -> Optional[float]:
    """Calculate polymer couplings."""
    atom1 = mol.GetAtomWithIdx(h1_idx)
    atom2 = mol.GetAtomWithIdx(h2_idx)
    
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
    
    # Calculer distance
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
            if hasattr(e, '__str__'):
                print(f"Erreur RDKit dans polymers: {e}")
            return None
    
    # Detect a vinyl monomer.
    vinyl1 = detect_vinyl_monomer(mol, parent1_idx)
    vinyl2 = detect_vinyl_monomer(mol, parent2_idx)
    
    if vinyl1 or vinyl2:
        if distance == 1:  # On the same double bond
            # Approximate cis/trans geometry.
            return POLYMER_COUPLINGS['vinyl_trans']  # Default to trans.
        elif distance == 0:  # Geminal
            return POLYMER_COUPLINGS['vinyl_geminal']
    
    # Generic polymer coupling
    if distance == 2:  # Vicinal
        return POLYMER_COUPLINGS['polymer_backbone']
    
    return None


def _test_polymers():
    """Test the polymers plugin."""
    print("=" * 60)
    print("TEST DU PLUGIN POLYMÈRES")
    print("=" * 60)
    
    test_cases = [
        ("Styrène", "C=Cc1ccccc1", 17.0),
        ("Acrylate", "C=CC(=O)O", 17.0),
        ("Éthylène", "C=C", 17.0),
    ]
    
    passed = 0
    
    for name, smiles, expected in test_cases:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol = Chem.AddHs(mol)
            h_atoms = [i for i, atom in enumerate(mol.GetAtoms()) if atom.GetSymbol() == 'H']
            if len(h_atoms) >= 2:
                J = get_coupling(mol, h_atoms[0], h_atoms[1])
                if J and abs(J - expected) <= 3.0:
                    print(f"PASS: {name}: {J:.1f} Hz")
                    passed += 1
                else:
                    print(f"FAIL: {name}: {J} Hz (expected: {expected})")
            else:
                print(f"WARN: {name}: Not enough H atoms")
        else:
            print(f"FAIL: {name}: Invalid SMILES")
    
    return passed, len(test_cases)


if __name__ == "__main__":
    _test_polymers()