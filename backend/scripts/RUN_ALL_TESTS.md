# Tests de validation - Étapes 1-6

Ce document liste tous les tests disponibles pour valider les améliorations du système de prédiction NMR.

## Structure du système (14 features)

```python
beyond_tuple = (
    carbonyl_type,          # 0  - Type de carbonyle (cétone/ester/amide/acide)
    hetero_ring,            # 1  - Hétérocycle (taille + type)
    charge_bin,             # 2  - Charge formelle binned
    donor_acceptor,         # 3  - Donneur/accepteur H
    alpha_hetero,           # 4  - Hétéroatome en alpha (O/N/S/X)
    aromatic_type,          # 5  - Type d'aromatique (benzène/pyridine/indole...)
    allylic_pos,            # 6  - Position allylique/prényle
    vinyl_conj,             # 7  - Vinyl conjugué
    relative_config_code,   # 8  - ✨ ÉTAPE 3: Config relative (fusion énantiomères)
    ez_code,                # 9  - Stéréochimie E/Z
    nH,                     # 10 - ✨ ÉTAPE 1: Nombre d'H attachés (0-3)
    hybridization_code,     # 11 - ✨ ÉTAPE 1: Hybridation (sp/sp2/sp3)
    dist_heteroatom,        # 12 - ✨ ÉTAPE 4: Distance au plus proche O/N/S/halogène
    dist_aromatic           # 13 - ✨ ÉTAPE 4: Distance au plus proche aromatique
)
```

## Tests disponibles

### 1. test_exact_match.py
**Objectif**: Tester la fusion des énantiomères et la séparation E/Z

**Ce qui est testé**:
- Morgan FP seul fusionne TOUS les stéréoisomères
- Beyond features avec `include_chiral=False` fusionne énantiomères mais sépare E/Z
- Composite FP (Morgan + Beyond) donne 2 groupes (chimiquement correct)

**Molécules de test**:
- Neoechinulin A avec 3 stéréoisomères:
  - `C[C@H]1...` (SDF @ \\)
  - `C[C@@H]1...` (User @@ \\)  ← Énantiomère du premier
  - `C[C@H]1.../` (User @ /)    ← Diastéréomère E/Z

**Résultat attendu**:
```
Morgan fingerprints uniques: 1/3 (tous identiques)
Beyond features uniques: 2/3 (@ et @@ fusionnés, E/Z séparés)
Composite fingerprints uniques: 2/3 (chimiquement correct!)
```

**Comment lancer**:
```bash
# Via launcher Python (recommandé)
cd backend
python launcher.py test_exact_match

# Ou directement (si RDKit disponible)
cd backend/scripts
python test_exact_match.py
```

---

### 2. test_step4_distances.py
**Objectif**: Tester les features de distance BFS (étape 4)

**Ce qui est testé**:
- `dist_heteroatom`: Distance au plus proche O/N/S/halogène
- `dist_aromatic`: Distance au plus proche atome aromatique

**Molécules de test**:
1. **Acétate d'éthyle** `CC(=O)OCC`:
   - CH₃ ester: d_hetero=2 (via C=O)
   - CH₃ alkyle: d_hetero=2 (via O)
   - CH₂: d_hetero=1 (O direct)

2. **Toluène** `Cc1ccccc1`:
   - CH₃ benzylique: d_hetero=8, d_arom=1

3. **Propane** `CCC`:
   - Tous: d_hetero=8, d_arom=8 (aucun hétéro ni aromatique)

4. **Éthanol** `CCO`:
   - CH₂: d_hetero=1 (O direct)
   - CH₃: d_hetero=2 (via C-O)

**Résultat attendu**:
Les distances doivent correspondre aux structures chimiques.

**Comment lancer**:
```bash
cd backend/scripts
python test_step4_distances.py
```

---

### 3. test_step5_bimodal.py
**Objectif**: Tester la détection de distributions bimodales (étape 5)

**Ce qui est testé**:
- **Cas 1**: Distribution unimodale normale → False (pas de faux positif)
- **Cas 2**: Bimodal évident (CH₃ ester 2.0 ppm + CH₃ alkyle 0.9 ppm) → True
- **Cas 3**: Bimodal subtil (13C: 128 ppm + 120 ppm) → True
- **Cas 4**: Distribution uniforme (borderline) → Possiblement True
- **Cas 5**: Avec outliers isolés → False (robuste)
- **Cas 6**: Trimodal → True (généralisation)
- **Cas 7**: n < 10 valeurs → False (pas assez de données)

**Algorithmes utilisés**:
1. Coefficient de bimodalité de Sarle: `b = (skew² + 1) / kurtosis`
2. Gap analysis: chercher gap > 2×std avec ≥3 valeurs de chaque côté

**Résultat attendu**:
Tous les cas doivent retourner les valeurs attendues.

**Comment lancer**:
```bash
cd backend/scripts
python test_step5_bimodal.py
```

---

### 4. test_step6_validation.py ⭐ **BATTERIE COMPLÈTE**
**Objectif**: Validation complète des étapes 1-6

**7 tests inclus**:

#### Test 1: Énantiomères fusionnés
- Molécules: Neoechinulin A @ et @@
- Attendu: Même composite fingerprint
- Validation: `include_chiral=False` fonctionne

#### Test 2: Diastéréomères E/Z séparés
- Molécules: Neoechinulin A (@ \\ vs @ /)
- Attendu: Fingerprints différents
- Validation: `ez_code` capture la stéréochimie E/Z

#### Test 3: CH₃ différents séparés
- Molécule: Acétate d'éthyle `CC(=O)OCC`
- Attendu: CH₃ ester ≠ CH₃ alkyle
- Validation: `nH`, `hybridization`, `dist_heteroatom` discriminent

#### Test 4: CH₃ équivalents fusionnés
- Molécule: DMSO `CS(C)=O`
- Attendu: Les deux CH₃ ont le même fingerprint
- Validation: Symétrie chimique respectée

#### Test 5: CH₂ benzylique vs alkyle
- Molécules: Toluène vs Propane
- Attendu: CH₃ benzylique ≠ CH₂ alkyle
- Validation: `dist_aromatic` fonctionne (1 vs 8)

#### Test 6: Détection bimodale
- Distributions: Unimodale vs bimodale
- Attendu: False pour unimodale, True pour bimodale
- Validation: Pas de faux positifs

#### Test 7: Config relative 2+ centres
- Molécules: (R,R) vs (R,S)
- Attendu: `relative_config_code` = 0 avec `include_chiral=False`
- Validation: Diastéréomères séparés par Morgan + autres features

**Résultat attendu**:
```
✓ Tests réussis: 7/7

✓ ✓ ✓ TOUS LES TESTS SONT RÉUSSIS!

Le système est prêt pour la production:
✓ Étape 1: nH + hybridization
✓ Étape 2: Cache version 13.0
✓ Étape 3: Config relative
✓ Étape 4: dist_hetero + dist_arom
✓ Étape 5: Détection bimodale
✓ Étape 6: TOUS LES TESTS VALIDÉS
```

**Comment lancer**:
```bash
cd backend/scripts
python test_step6_validation.py
```

---

### 5. test_famciclovir.py
**Objectif**: Test sur molécule réelle (Famciclovir)

**Ce qui est testé**:
- Prédictions 1H et 13C avec le système composite
- Vérification que toutes les prédictions retournent des valeurs

**Molécule**: Famciclovir `CC(=O)OCC(CCn1cnc2c(N)nc(nc12)N)COC(=O)C`

**Comment lancer**:
```bash
cd backend/scripts
python test_famciclovir.py
```

---

### 6. check_solvents.py
**Objectif**: Vérifier la normalisation des solvants dans la DB

**Ce qui est testé**:
- Liste de tous les solvants uniques
- Vérification de la normalisation (CDCl3, DMSO-d6, etc.)
- Détection de variantes non-normalisées

**Comment lancer**:
```bash
cd backend/scripts
python check_solvents.py
```

---

### 7. check_inchikey.py
**Objectif**: Vérifier les InChIKey des énantiomères

**Ce qui est testé**:
- InChIKey de deux énantiomères
- Comparaison des couches (connectivité vs stéréochimie)
- Explication du système exact-match

**Comment lancer**:
```bash
cd backend/scripts
python check_inchikey.py
```

---

## Lancer tous les tests

### Option 1: Lancer chaque test individuellement

```bash
cd backend/scripts

# Tests principaux (ordre recommandé)
python test_exact_match.py        # Test 1: Énantiomères/E-Z
python test_step4_distances.py    # Test 2: Distances BFS
python test_step5_bimodal.py      # Test 3: Bimodalité
python test_step6_validation.py   # Test 4: Batterie complète ⭐

# Tests complémentaires
python test_famciclovir.py        # Molécule réelle
python check_solvents.py          # Vérif solvants
python check_inchikey.py          # Vérif InChIKey
```

### Option 2: Script de lancement automatique

```bash
# TODO: Créer run_all_tests.sh ou run_all_tests.bat
```

---

## Interprétation des résultats

### ✓ Tous les tests passent
→ Le système est prêt pour rebuild de la base et production

### ✗ Certains tests échouent

**Test 1 échoue (énantiomères)**:
- Vérifier `encode_relative_chirality()` dans `beyond_features_stereo.py`
- Vérifier que `include_chiral=False` est bien utilisé

**Test 2 échoue (E/Z)**:
- Vérifier `get_stereo_environment()` dans `beyond_features_stereo.py`
- Vérifier que les double liaisons sont détectées

**Test 3 échoue (CH₃ différents)**:
- Vérifier que `nH`, `hybridization_code`, `dist_heteroatom` sont bien calculés
- Vérifier le tuple complet (14 features)

**Test 4 échoue (CH₃ équivalents)**:
- Problème possible avec les fingerprints Morgan (sensibilité à l'ordre des atomes)
- Vérifier la canonicalisation RDKit

**Test 5 échoue (benzylique vs alkyle)**:
- Vérifier `_distance_to_nearest_aromatic()` dans `beyond_features_stereo.py`
- Vérifier BFS sur cycles aromatiques

**Test 6 échoue (bimodalité)**:
- Ajuster les seuils dans `_is_bimodal()` (ligne ~2346)
- Vérifier que scipy est installé

**Test 7 échoue (config relative)**:
- Vérifier `encode_relative_chirality()` avec 2+ centres chiraux
- Vérifier la détection des relations like/unlike

---

## Prochaines étapes après validation

1. **Rebuild de la base de données**:
   ```bash
   cd backend/scripts
   rm ../ml_dependencies/data/fingerprint_database.pkl
   python fingerprint_db_predictor.py  # Force rebuild
   ```

2. **Tests sur molécules réelles**:
   - Famciclovir
   - Neoechinulin A
   - Caféine
   - Acétate d'éthyle

3. **Comparaison avant/après**:
   - Nombre de prédictions fingerprint vs ML
   - Précision (MAE) sur dataset de validation
   - Détection de fusions erronées

4. **Mise en production**:
   - Intégration dans `integrated_nmr_prediction.py`
   - Tests end-to-end via API
   - Documentation utilisateur

---

## Debugging

### Activer les logs détaillés

```python
import logging
logger = logging.getLogger('fingerprint_db_predictor')
logger.setLevel(logging.DEBUG)
```

### Vérifier le cache

```python
from fingerprint_db_predictor import FingerprintPredictor

predictor = FingerprintPredictor()
print(f"Cache loaded: {predictor.is_available()}")
print(f"DB version: {FingerprintPredictor._db_mtime}")
```

### Inspecter un fingerprint

```python
from rdkit import Chem
from beyond_features_stereo import compute_beyond_with_stereo_safe, beyond_stereo_to_string
from fingerprint_db_predictor import get_atom_fingerprint, RADIUS_ULTRA_EXACT

mol = Chem.MolFromSmiles("CC(=O)OCC")
mol_h = Chem.AddHs(mol)

for atom_idx in range(mol_h.GetNumAtoms()):
    beyond = compute_beyond_with_stereo_safe(mol_h, atom_idx, include_chiral=False)
    fp = get_atom_fingerprint(mol_h, atom_idx, RADIUS_ULTRA_EXACT)

    print(f"Atom {atom_idx}: {mol_h.GetAtomWithIdx(atom_idx).GetSymbol()}")
    print(f"  FP: {fp}")
    print(f"  Beyond: {beyond_stereo_to_string(beyond)}")
    print()
```

---

## Contact et support

Pour toute question sur les tests ou l'implémentation, consulter:
- `RESUME_IMPLEMENTATION.md` - Résumé des étapes 1-6
- `backend/scripts/beyond_features_stereo.py` - Implémentation des features
- `backend/scripts/fingerprint_db_predictor.py` - Système de prédiction
