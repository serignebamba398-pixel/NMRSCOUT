#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tracer exactement comment predict_atom() trouve le shift de 1.241 ppm
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from rdkit import RDLogger, Chem
RDLogger.DisableLog('rdApp.*')

from fingerprint_db_predictor import FingerprintPredictor
from beyond_features_stereo import ensure_gasteiger_charges
import numpy as np

if __name__ == '__main__':
    smiles = "C[C@@H]1C(=O)N/C(=C\\C2=C(NC3=CC=CC=C32)C(C)(C)C=C)/C(=O)N1"
    solvent = "CDCl3"

    print("=" * 80)
    print("TRAÇAGE PREDICTION CH3")
    print("=" * 80)
    print(f"\nMolécule: {smiles}")
    print()

    mol = Chem.MolFromSmiles(smiles)
    mol_h = Chem.AddHs(mol)
    ensure_gasteiger_charges(mol_h)

    predictor = FingerprintPredictor()

    # CH3 H index
    ch3_h_idx = 24

    print(f"Prédiction pour H{ch3_h_idx} (CH3)...")
    print()

    # Activer le logging pour voir ce qui se passe
    import logging
    logging.basicConfig(level=logging.DEBUG)

    result = predictor.predict_atom(mol_h, ch3_h_idx, solvent, '1H')

    if result:
        shift, std, count, radius = result
        print()
        print("=" * 80)
        print("RÉSULTAT")
        print("=" * 80)
        print(f"  Shift: {shift:.3f} ppm")
        print(f"  Std:   {std:.3f} ppm")
        print(f"  Count: {count}")
        print(f"  Radius: R{radius}")
        print()

        # Maintenant, récupérer manuellement les shifts pour analyser la distribution
        print("=" * 80)
        print("RÉCUPÉRATION MANUELLE DES SHIFTS")
        print("=" * 80)
        print()

        # Essayer de récupérer via _predict_with_composite_fp_internal
        # On doit accéder à la méthode interne

        from fingerprint_db_predictor import get_atom_fingerprint
        from beyond_features_stereo import compute_beyond_with_stereo_safe

        beyond = compute_beyond_with_stereo_safe(mol_h, ch3_h_idx, include_chiral=False)

        print(f"Features beyond: {beyond}")
        print()

        # Essayer chaque radius
        for r in [6, 5, 4, 3, 2]:
            fp = get_atom_fingerprint(mol_h, ch3_h_idx, r)
            composite_key = (fp, beyond)

            db_composite = predictor._db_composite

            if composite_key in db_composite:
                entry = db_composite[composite_key]

                # Essayer le solvant exact
                if solvent in entry:
                    if '1H' in entry[solvent]:
                        shifts = entry[solvent]['1H']
                        print(f"✓ Trouvé avec R{r} (solvant exact {solvent})")
                        print(f"  Nombre de shifts: {len(shifts)}")
                        print(f"  Moyenne: {np.mean(shifts):.3f} ppm")
                        print(f"  Std: {np.std(shifts):.3f} ppm")
                        print()

                        # Distribution
                        shifts_array = np.array(shifts)
                        print("  Distribution:")
                        bins = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 10.0]
                        hist, _ = np.histogram(shifts_array, bins=bins)

                        for i in range(len(bins)-1):
                            cnt = hist[i]
                            pct = (cnt / len(shifts)) * 100
                            bar = "█" * int(pct / 2)
                            print(f"    {bins[i]:.1f}-{bins[i+1]:.1f} ppm: {cnt:3d} ({pct:5.1f}%) {bar}")

                        print()

                        # Top 20 et bottom 20
                        sorted_shifts = sorted(shifts)
                        print(f"  20 plus petits shifts:")
                        for i in range(min(20, len(sorted_shifts))):
                            print(f"    [{i:2d}] {sorted_shifts[i]:.3f} ppm")

                        print()
                        print(f"  20 plus grands shifts:")
                        start = max(0, len(sorted_shifts) - 20)
                        for i in range(start, len(sorted_shifts)):
                            print(f"    [{i:2d}] {sorted_shifts[i]:.3f} ppm")

                        print()
                        print(f"  Valeur expérimentale: 1.530 ppm")
                        print(f"  Erreur: {abs(np.mean(shifts) - 1.530):.3f} ppm")

                        # Combien sont proches?
                        close = np.sum((shifts_array >= 1.4) & (shifts_array <= 1.6))
                        print(f"  Shifts dans [1.4-1.6]: {close}/{len(shifts)} ({100*close/len(shifts):.1f}%)")

                        break

                # Essayer fallback multi-solvants
                all_shifts = []
                for sol_key in entry:
                    if '1H' in entry[sol_key]:
                        all_shifts.extend(entry[sol_key]['1H'])

                if all_shifts and solvent not in entry:
                    print(f"✓ Trouvé avec R{r} (fallback multi-solvants)")
                    print(f"  Nombre total de shifts: {len(all_shifts)}")
                    print(f"  Moyenne: {np.mean(all_shifts):.3f} ppm")
                    break

    else:
        print("❌ Pas de résultat")
