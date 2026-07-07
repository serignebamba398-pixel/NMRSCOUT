# Instructions pour tracer le spectre de Flavoglaucin

## 🎯 Objectif

Tracer le spectre RMN ¹H de Flavoglaucin avec simulation Hamiltonienne exacte et fragmentation automatique pour les systèmes > 10 spins.

## 📋 Prérequis

1. **Environnement conda** `rdkit_env_clean` avec:
   - Python 3.10+
   - RDKit
   - NumPy
   - SciPy
   - Matplotlib

2. **Fichier JSON de multiplicité**:
   - Flavoglaucin_multiplicity.json (déjà présent)

## 🚀 Méthode rapide (recommandée)

### Windows

Double-cliquez sur: **`trace_spectre.bat`**

Le script va automatiquement:
1. ✓ Activer l'environnement rdkit_env_clean
2. ✓ Charger le fichier JSON de Flavoglaucin
3. ✓ Générer le spectre avec fragmentation Hamiltonienne
4. ✓ Tracer 3 vues (complet + 2 zooms)
5. ✓ Sauvegarder en PNG (flavoglaucin_spectrum.png)
6. ✓ Afficher le spectre interactif

### Linux/Mac

```bash
conda activate rdkit_env_clean
python trace_flavoglaucin.py
```

## 📊 Ce que vous allez voir

### Console (logs)

```
================================================================================
SPECTRE FLAVOGLAUCIN - AVEC FRAGMENTATION HAMILTONIENNE
================================================================================

Chargement: C:\Users\...\Flavauglocin_multiplicity.json

Molecule: Flavauglocin
SMILES: CCCCCCCC1=C(C=C(C(=C1C=O)O)CC=C(C)C)O
Frequence: 700 MHz
Solvant: CDCl3
Nombre de protons: 28
Nombre de groupes: 15

Systemes second ordre detectes:
  - Groupe 1: m @ 1.29 ppm (AA'BB')
  - Groupe 2: m @ 1.29 ppm (AA'BB')
  - Groupe 3: m @ 1.29 ppm (AA'BB')
  - Groupe 4: m @ 1.39 ppm (AA'BB')
  - Groupe 5: m @ 1.58 ppm (AA'BB')
  - Groupe 6: m @ 2.88 ppm (AA'BB')

================================================================================
GENERATION DU SPECTRE
================================================================================

Appel de la simulation Hamiltonienne...
(La fragmentation sera appliquee automatiquement si necessaire)

[INFO] Using Version2 exact Hamiltonian for 5 spins
[INFO] Using Version2 exact Hamiltonian for 2 spins
...
[INFO] Système de 12 spins > 10, tentative de fragmentation  ← FRAGMENTATION!
[INFO] Fragmentation: 12 spins, seuil J/Δν = 0.05
[INFO] → 3 fragments détectés
[INFO]   Fragment 1: 5 spins [22, 23, 24, 25, 26]
[INFO]   Fragment 2: 4 spins [27, 28, 29, 30]
[INFO]   Fragment 3: 3 spins [31, 32, 33]
[INFO]   Fragment 1: 48 transitions
[INFO]   Fragment 2: 24 transitions
[INFO]   Fragment 3: 12 transitions
[INFO] Total: 84 transitions combinées

Spectre genere avec succes!
  - Nombre de points: 16384
  - Plage ppm: 0.50 - 12.00
  - Intensite max: 1.000

================================================================================
AFFICHAGE DU SPECTRE
================================================================================

Spectre sauvegarde: flavoglaucin_spectrum.png

Affichage du spectre...
```

### Figure matplotlib (3 panneaux)

**Panneau 1 (haut)**: Spectre complet 0.5-12.0 ppm
- Vue d'ensemble de tous les signaux
- Titre: "Flavauglocin - 700 MHz - CDCl3 - SPECTRE COMPLET"

**Panneau 2 (milieu)**: Zoom 0.5-3.5 ppm
- Région aliphatique
- CH₃ triplet @ 0.88 ppm
- Chaîne alkyle @ 1.29-1.58 ppm (avec fragmentation!)
- CH₂ benzylique @ 2.88 ppm

**Panneau 3 (bas)**: Zoom 5.0-12.0 ppm
- Région aromatique
- =CH vinylique @ 5.28 ppm
- Ar-H @ 6.89 ppm
- CHO aldéhyde @ 10.25 ppm
- OH phénolique @ 11.92 ppm

## 🔍 Points clés à observer

### 1. Fragmentation en action

La **chaîne alkyle C₈** contient potentiellement 16 protons couplés:
- **Sans fragmentation**: Impossible (2¹⁶ = 65536×65536 matrice)
- **Avec fragmentation**: 3 fragments simulés séparément puis combinés

Cherchez dans les logs:
```
[INFO] Système de 12 spins > 10, tentative de fragmentation
[INFO] → 3 fragments détectés
```

### 2. Systèmes AA'BB' simulés exactement

Plusieurs systèmes AA'BB' détectés dans la chaîne:
- Simulation quantique exacte (pas d'approximation)
- Patterns de roof effect caractéristiques
- Espacements J précis

### 3. Résolution à 700 MHz

Haute résolution permettant de:
- Séparer les CH₂ proches (1.29, 1.39, 1.58 ppm)
- Distinguer les couplages vicinaux (³J) et géminaux (²J)
- Voir les multiplets complexes

## 📁 Fichiers générés

Après exécution:

1. **flavoglaucin_spectrum.png** (300 dpi)
   - Image haute qualité du spectre
   - 3 panneaux (complet + 2 zooms)
   - Prêt pour publication/rapport

2. **Logs dans la console**
   - Détails de la fragmentation
   - Nombre de transitions par fragment
   - Statistiques du spectre

## 🛠️ Dépannage

### Erreur: "No module named 'numpy'"

```bash
# Vérifier l'environnement actif
conda info --envs

# Activer le bon environnement
conda activate rdkit_env_clean

# Réessayer
python trace_flavoglaucin.py
```

### Erreur: "Fichier non trouvé"

Vérifier le chemin dans `trace_flavoglaucin.py` ligne 23:
```python
json_path = r"C:\Users\...\Flavauglocin_multiplicity.json"
```

Modifier si nécessaire pour pointer vers votre fichier.

### Fragmentation désactivée

Si vous voyez:
```
[WARN] Fragmentation échouée, fallback approximation premier ordre
```

**Raison possible**: networkx non installé

**Solution**:
```bash
conda activate rdkit_env_clean
pip install networkx
python trace_flavoglaucin.py
```

### Spectre ne s'affiche pas

**Causes possibles**:
1. Backend matplotlib mal configuré
2. Problème de permissions d'affichage

**Solution**:
- Le spectre est sauvegardé dans `flavoglaucin_spectrum.png`
- Ouvrir ce fichier manuellement si l'affichage échoue

## 📊 Analyse du spectre

### Intégration attendue (28H total)

| Région ppm | Signal | Intégration | Attribution |
|------------|--------|-------------|-------------|
| 0.88 | t | 3H | CH₃ terminal |
| 1.29-1.58 | m | 12H | Chaîne (CH₂)₆ |
| 1.70-1.76 | d | 6H | 2×CH₃ allyliques |
| 2.88 | m | 2H | CH₂ benzylique |
| 3.29 | d | 2H | CH₂-C= |
| 5.28 | ddt | 1H | =CH vinylique |
| 6.89 | s | 1H | Ar-H |
| 10.25 | s | 1H | CHO |
| 4.31, 11.92 | s (broad) | 2H | OH échangeables |

### Multiplicités caractéristiques

- **t (triplet)** @ 0.88: CH₃ couplé avec CH₂ adjacent (³J ≈ 7 Hz)
- **m (multiplet)** @ 1.29-1.58: Systèmes AA'BB' de la chaîne
- **d (doublet)** @ 1.70, 1.76: CH₃ allyliques (⁴J ≈ 1.5 Hz)
- **ddt** @ 5.28: Proton vinylique avec 3 couplages différents
- **s (singulet)** @ 6.89, 10.25: Protons isolés sans couplage vicinal

## 📚 Références

- **Fragmentation**: Castillo et al., J. Magn. Reson. 209, 123-130 (2011)
- **Documentation**: backend/scripts/FRAGMENTATION_README.md
- **Code source**: backend/scripts/sim_mult_rulebased.py

## ✅ Checklist avant de commencer

- [ ] Environnement conda créé et activé
- [ ] Fichier Flavoglaucin_multiplicity.json présent
- [ ] Packages installés (numpy, scipy, matplotlib, rdkit)
- [ ] Optionnel: networkx installé pour fragmentation

## 🎯 Résultat attendu

Vous devriez voir:
1. ✓ Console avec logs détaillés de la simulation
2. ✓ Détection et fragmentation automatique des grands systèmes
3. ✓ Fenêtre matplotlib avec 3 panneaux de spectre
4. ✓ Fichier PNG sauvegardé (300 dpi, haute qualité)
5. ✓ Spectre réaliste avec patterns de second ordre corrects

**Bon traçage ! 🎉**
