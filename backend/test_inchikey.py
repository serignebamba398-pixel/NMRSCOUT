#!/usr/bin/env python3
"""Test InChIKey comparison for L-aspartic acid"""
from rdkit import Chem
from rdkit.Chem import inchi as rdkit_inchi
import pickle
import os

# 1. InChIKey du SMILES utilisateur
smiles = 'C([C@@H](C(=O)O)N)C(=O)O'
mol_user = Chem.MolFromSmiles(smiles)
ik_user = rdkit_inchi.MolToInchiKey(mol_user)
print(f'SMILES utilisateur: {smiles}')
print(f'InChIKey utilisateur: {ik_user}')

# 2. Charger le SDF BMRB et trouver bmse000031
sdf_path = os.path.join(os.path.dirname(__file__), 'ml_dependencies/data/sdf_bmrb_rdkit_v2/bmrb_all_molecules.sdf')
suppl = Chem.SDMolSupplier(sdf_path)
for mol in suppl:
    if mol is None:
        continue
    bmrb_id = mol.GetProp('BMRB_ID') if mol.HasProp('BMRB_ID') else ''
    if bmrb_id == 'bmse000031':
        name = mol.GetProp('Name') if mol.HasProp('Name') else '?'
        mol_no_h = Chem.RemoveHs(mol)
        ik_bmrb = rdkit_inchi.MolToInchiKey(mol_no_h)
        print(f'\nBMRB bmse000031: {name}')
        print(f'InChIKey BMRB:     {ik_bmrb}')
        print(f'Match exact: {ik_bmrb == ik_user}')

        # Comparer les composants
        user_parts = ik_user.split('-')
        bmrb_parts = ik_bmrb.split('-')
        print(f'\nComparaison InChIKey:')
        print(f'  Connectivity (1ere partie): {user_parts[0]}')
        print(f'                              {bmrb_parts[0]}')
        print(f'                              -> {"MATCH" if user_parts[0]==bmrb_parts[0] else "DIFF"}')
        print(f'  Stereo (2e partie): {user_parts[1]} vs {bmrb_parts[1]} -> {"MATCH" if user_parts[1]==bmrb_parts[1] else "DIFF"}')
        break

# 3. Verifier le pickle
pkl_path = os.path.join(os.path.dirname(__file__), 'ml_dependencies/data/fingerprint_database.pkl')
if os.path.exists(pkl_path):
    with open(pkl_path, 'rb') as f:
        data = pickle.load(f)
    exact_db = data.get('exact_db', {})
    print(f'\n--- Fingerprint Database ---')
    print(f'Molecules dans exact_db: {len(exact_db)}')

    if ik_user in exact_db:
        print(f'InChIKey utilisateur TROUVE dans exact_db!')
        print(f'  Solvants: {list(exact_db[ik_user].keys())}')
    else:
        print(f'InChIKey utilisateur PAS dans exact_db')
        # Chercher des InChIKeys similaires (meme connectivity)
        user_conn = ik_user.split('-')[0]
        similar = [ik for ik in exact_db.keys() if ik.split('-')[0] == user_conn]
        if similar:
            print(f'  InChIKeys avec meme connectivity:')
            for ik in similar[:5]:
                print(f'    {ik}')
