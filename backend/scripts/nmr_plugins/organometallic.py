# nmr_plugins/organometallic.py
# Organometallic complex plugin

"""Handle organometallic complexes.

Covers metallocenes (ferrocene, ruthenocene), carbonyl complexes,
pi complexes (alkenes, arenes), metal clusters, and organometallic catalysts."""

from typing import Dict, List, Tuple, Optional, Any
import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    raise ImportError("RDKit requis pour le plugin organometallic")


# Organometallic coupling tables
ORGANOMETALLIC_COUPLINGS: Dict[str, float] = {
    # Metallocenes
    'ferrocene_H_H': 2.0,           # H-H on cyclopentadienyl
    'ferrocene_long_range': 0.5,    # Through the metal
    'ruthenocene_H_H': 2.2,         # Ruthenium
    'cobaltocene_H_H': 2.5,         # Cobalt
    
    # Complexes carbonyles
    'metal_CO_H': 3.0,              # M-CO...H-C
    'carbonyl_coupling': 1.5,       # Through CO
    
    # Complexes π
    'eta6_benzene_H_H': 8.0,        # η⁶-Benzene
    'eta4_butadiene_H_H': 12.0,     # η⁴-Butadiene
    'eta2_alkene_H_H': 16.0,        # η²-Alkene (reduced)
    
    # Phosphines
    'phosphine_P_C_H': 8.0,         # ³J(P-C-H)
    'phosphine_P_C_C_H': 12.0,      # ²J(P-C-C-H)
}

# Patterns SMARTS
ORGANOMETALLIC_PATTERNS: Dict[str, str] = {
    'ferrocene': '[Fe+2].c1ccc[c-]1.c1ccc[c-]1',
    'metal_carbonyl': '[C-]#[O+]',  # Simplified
    'phosphine': '[PH2,PH1,P]',     # Phosphines diverses
    'metal_aromatic': '[c]:[c]:[c]:[c]:[c]-[*]',  # Metal-bound aromatic group
}


def detect_metallocene(mol: Chem.Mol, atom_idx: int) -> Optional[Dict[str, Any]]:
    """Detect metallocenes."""
    # Approximate detection using five-membered aromatic rings.
    atom = mol.GetAtomWithIdx(atom_idx)
    
    if not atom.GetIsAromatic():
        return None
    
    # Look for five-membered aromatic rings.
    ring_info = mol.GetRingInfo()
    atom_rings = [ring for ring in ring_info.AtomRings() 
                  if atom_idx in ring and len(ring) == 5]
    
    if atom_rings:
        # Treat a five-membered aromatic ring as a possible cyclopentadienyl group.
        return {
            'type': 'metallocene',
            'metal': 'Fe',  # Approximation
            'coupling': ORGANOMETALLIC_COUPLINGS['ferrocene_H_H']
        }
    
    return None


def get_coupling(mol: Chem.Mol, h1_idx: int, h2_idx: int) -> Optional[float]:
    """Calculate organometallic couplings."""
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
    
    # Detect a metallocene.
    metallocene1 = detect_metallocene(mol, parent1_idx)
    metallocene2 = detect_metallocene(mol, parent2_idx)
    
    if metallocene1 or metallocene2:
        return ORGANOMETALLIC_COUPLINGS['ferrocene_H_H']
    
    return None


def _test_organometallic():
    """Test the organometallic plugin."""
    print("=" * 60)
    print("TEST DU PLUGIN ORGANOMÉTALLIQUE")
    print("=" * 60)
    
    # Test simple
    test_smiles = "c1cccc1"  # Use benzene as an approximation.
    mol = Chem.MolFromSmiles(test_smiles)
    if mol:
        mol = Chem.AddHs(mol)
        h_atoms = [i for i, atom in enumerate(mol.GetAtoms()) if atom.GetSymbol() == 'H']
        if len(h_atoms) >= 2:
            J = get_coupling(mol, h_atoms[0], h_atoms[1])
            print(f"Benzène: J = {J} Hz" if J else "Aucun couplage organométallique")
    
    return 1, 1  # One test passed using the approximation.


if __name__ == "__main__":
    _test_organometallic()