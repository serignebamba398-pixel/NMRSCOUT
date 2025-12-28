#!/usr/bin/env python
"""Analyser le mapping des positions SMARTS et les valeurs J"""
from rdkit import Chem
from rdkit.Chem import rdmolops
import numpy as np
from heterocycles_complete import identify_heterocycle, _FUSED_BICYCLIC_MISSING
from heterocycles_complete import estimate_heterocycle_couplings_complete

# Quinoxaline
smiles = 'c1ccc2nccnc2c1'
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)

print('=== Quinoxaline - Test estimate_heterocycle_couplings_complete ===')

# Construire heavy_to_H
heavy_to_H = {}
for i, atom in enumerate(mol.GetAtoms()):
    if atom.GetSymbol() == 'C' and atom.GetIsAromatic():
        h_neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == 'H']
        if h_neighbors:
            heavy_to_H[i] = h_neighbors

print(f'heavy_to_H: {heavy_to_H}')
print()

# Construire D (matrice de distances)
D = rdmolops.GetDistanceMatrix(mol)

# Tester estimate_heterocycle_couplings_complete sur chaque C aromatique
for heavy_idx, h_list in heavy_to_H.items():
    result = estimate_heterocycle_couplings_complete(
        mol=mol,
        heavy_idx=heavy_idx,
        h_list=h_list,
        heavy_to_H=heavy_to_H,
        D=D,
        min_j_threshold=0.3
    )
    if result:
        j_vals = result.get('J_values', [])
        mult = result.get('multiplicity', '?')
        source = result.get('source', '?')
        partners = result.get('coupling_partners', [])
        print(f'C{heavy_idx}: mult={mult}, J={j_vals}, partners={partners}, source={source}')
