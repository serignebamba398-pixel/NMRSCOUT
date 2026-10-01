#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Rebuild the fingerprint database to version 13.0 (14 features). Supports Windows multiprocessing."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def main():
    """Main function protected for Windows multiprocessing"""
    from fingerprint_db_predictor import FingerprintPredictor
    import time

    print("=" * 80)
    print("REBUILD FINGERPRINT DATABASE - Version 13.0")
    print("=" * 80)
    print()
    print("Nouvelles features (11 → 14):")
    print("  [8]  relative_config (Étape 3) - Fusion énantiomères")
    print("  [10] nH (Étape 1) - Nombre d'hydrogènes")
    print("  [11] hybridization_code (Étape 1) - SP/SP2/SP3")
    print("  [12] dist_heteroatom (Étape 4) - Distance BFS → O/N/S/X")
    print("  [13] dist_aromatic (Étape 4) - Distance BFS → aromatique")
    print()
    print("Estimated duration: 30-60 minutes (parallel processing)")
    print()

    # Check for an existing database.
    db_path = os.path.join(os.path.dirname(__file__), '..', 'ml_dependencies', 'data', 'fingerprint_database.pkl')
    if os.path.exists(db_path):
        import time
        size_mb = os.path.getsize(db_path) / (1024 * 1024)
        mtime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(db_path)))
        print(f"   Base existante détectée:")
        print(f"   Taille: {size_mb:.1f} MB")
        print(f"   Date: {mtime}")
        print(f"   → Sera remplacée par version 13.0")
        print()

    print("Lancement du rebuild...")
    print()

    start_time = time.time()

    try:
        # Constructing the predictor triggers an automatic rebuild
        # when the version differs (11.0/12.0 to 13.0).
        predictor = FingerprintPredictor()

        elapsed = time.time() - start_time
        elapsed_min = elapsed / 60

        print()
        print("=" * 80)
        print("REBUILD TERMINÉ AVEC SUCCÈS!")
        print("=" * 80)
        print()
        print(f"  Durée: {elapsed_min:.1f} minutes ({elapsed:.0f} secondes)")
        print()

        # Check the rebuilt database.
        if os.path.exists(db_path):
            size_mb = os.path.getsize(db_path) / (1024 * 1024)
            print(f" Nouvelle base créée:")
            print(f"   Version: 13.0")
            print(f"   Taille: {size_mb:.1f} MB")
            print(f"   Features: 14 (au lieu de 11)")
            print()

        # Statistiques
        if predictor.is_available():
            print(" Statistiques:")
            print(f"   Base de données chargée: OUI")
            print(f"   Prêt pour les prédictions: OUI")
            print()

        print("✓ Le système est maintenant prêt avec les 14 features!")
        print()
        print("Prochaines étapes:")
        print("  1. Lancer les tests de validation")
        print("     → python test_features_only.py")
        print("  2. Tester sur molécules réelles")
        print("     → python test_famciclovir.py")
        print("  3. Démarrer le serveur et tester via l'API")
        print("     → cd backend && npm start")
        print()

    except KeyboardInterrupt:
        print()
        print("Rebuild interrompu par l'utilisateur (Ctrl+C)")
        print("La base est incomplète et doit être reconstruite.")
        sys.exit(1)
    except Exception as e:
        print()
        print("=" * 80)
        print("ERREUR LORS DU REBUILD")
        print("=" * 80)
        print()
        print(f"Erreur: {e}")
        print()
        print("Debugging:")
        print("  1. Vérifier que les fichiers source existent:")
        print("     - backend/ml_dependencies/data/nmrshiftdb2withsignals.sd")
        print("     - backend/ml_dependencies/data/sdf_bmrb_rdkit_v2/")
        print("  2. Vérifier l'espace disque disponible (>500 MB)")
        print("  3. Vérifier les logs ci-dessus pour plus de détails")
        print()
        sys.exit(1)

if __name__ == '__main__':
    # Protection Windows multiprocessing
    import multiprocessing
    multiprocessing.freeze_support()
    main()
