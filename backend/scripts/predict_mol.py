#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rdkit import RDLogger, Chem
RDLogger.DisableLog('rdApp.*')
from fingerprint_db_predictor import FingerprintPredictor
from beyond_features_stereo import ensure_gasteiger_charges

smiles = sys.argv[1] if len(sys.argv) > 1 else ""
solvent = sys.argv[2] if len(sys.argv) > 2 else "CDCl3"

mol = Chem.MolFromSmiles(smiles)
if mol is None:
    print("SMILES invalide"); sys.exit(1)

mol_h = Chem.AddHs(mol)
ensure_gasteiger_charges(mol_h)
predictor = FingerprintPredictor()

results_1h, results_13c = {}, {}
for atom_idx in range(mol_h.GetNumAtoms()):
    atom = mol_h.GetAtomWithIdx(atom_idx)
    symbol = atom.GetSymbol()
    if symbol == 'C': nucleus = '13C'
    elif symbol == 'H': nucleus = '1H'
    else: continue
    result = predictor.predict_atom(mol_h, atom_idx, solvent, nucleus)
    if nucleus == '1H': results_1h[atom_idx] = result
    else: results_13c[atom_idx] = result

print(f"\n{'='*65}")
print(f"Prédiction FP — {solvent}")
print(f"{'='*65}")

print(f"\n--- 1H ({sum(1 for r in results_1h.values() if r)}/{len(results_1h)} FP) ---")
print(f"{'Idx':>3} | {'Shift':>7} | {'±Std':>5} | {'n':>4} | R  | Source")
print("-"*50)
for idx, res in sorted(results_1h.items()):
    if res:
        s,std,n,r = res
        src = 'EXACT' if n==99 else ('COMP/FP' if r>=4 else 'FP')
        print(f"{idx:3d} | {s:7.3f} | {std:5.3f} | {n:4d} | R{r} | {src}")
    else:
        print(f"{idx:3d} |     N/A |   --- |  --- | -- | ML")

print(f"\n--- 13C ({sum(1 for r in results_13c.values() if r)}/{len(results_13c)} FP) ---")
print(f"{'Idx':>3} | {'Shift':>7} | {'±Std':>5} | {'n':>4} | R  | Source")
print("-"*50)
for idx, res in sorted(results_13c.items()):
    if res:
        s,std,n,r = res
        src = 'EXACT' if n==99 else ('COMP/FP' if r>=4 else 'FP')
        print(f"{idx:3d} | {s:7.3f} | {std:5.3f} | {n:4d} | R{r} | {src}")
    else:
        print(f"{idx:3d} |     N/A |   --- |  --- | -- | ML")

fp_1h = sum(1 for r in results_1h.values() if r)
ml_1h = sum(1 for r in results_1h.values() if not r)
fp_13c = sum(1 for r in results_13c.values() if r)
ml_13c = sum(1 for r in results_13c.values() if not r)
print(f"\n1H:  {fp_1h}/{fp_1h+ml_1h} FP ({fp_1h/(fp_1h+ml_1h)*100:.1f}%), {ml_1h} ML")
print(f"13C: {fp_13c}/{fp_13c+ml_13c} FP ({fp_13c/(fp_13c+ml_13c)*100:.1f}%), {ml_13c} ML")
print(f"TOTAL: {fp_1h+fp_13c}/{fp_1h+ml_1h+fp_13c+ml_13c} FP ({(fp_1h+fp_13c)/(fp_1h+ml_1h+fp_13c+ml_13c)*100:.1f}%)")
