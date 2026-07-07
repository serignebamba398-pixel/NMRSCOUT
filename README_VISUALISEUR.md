# Visualiseur Interactif de Spectres RMN

## Vue d'ensemble

Visualiseur interactif pour spectres RMN ¹H avec:
- **Simulation Hamiltonienne exacte** pour les systèmes fortement couplés
- **Fragmentation automatique** pour les grands systèmes (> 10 spins)
- **Interface interactive** avec zoom et détection automatique de patterns
- **Analyse de multiplicité** en temps réel

## Installation

### Prérequis

1. **Conda** (Anaconda ou Miniconda)
2. **Environnement Python** avec les packages suivants:
   - Python 3.10+
   - RDKit
   - NumPy
   - SciPy
   - Matplotlib

### Créer l'environnement

```bash
conda create -n rdkit_env_clean python=3.10 rdkit numpy scipy matplotlib -c conda-forge
conda activate rdkit_env_clean
```

### Installation optionnelle de networkx (pour fragmentation)

```bash
pip install networkx
```

## Utilisation

### Méthode 1: Script batch (Windows)

Double-cliquez sur `run_interactive_viewer.bat`

Le script:
1. Active automatiquement l'environnement `rdkit_env_clean`
2. Lance le visualiseur interactif
3. Désactive l'environnement à la fin

### Méthode 2: Ligne de commande

```bash
conda activate rdkit_env_clean
python interactive_spectrum_viewer.py
```

## Fichiers de test

Le visualiseur cherche automatiquement un fichier JSON de multiplicité dans les emplacements suivants:

1. `C:\Users\...\Flavauglocin_multiplicity.json` (chemin personnalisé)
2. `test_output_molecule/Flavauglocin_multiplicity.json`
3. `test_output_steroid/test_steroid_multiplicity.json`
4. `test_output_valine_final/valine_multiplicity.json`

## Contrôles

### Boutons de zoom prédéfini (à gauche)
- **Complet**: Vue complète du spectre
- **CH3 0.88**: Zoom sur le triplet méthyle (0.8-0.96 ppm)
- **Chaîne 1.3**: Zoom sur la chaîne alkyle (1.2-1.6 ppm)
- **CH2 2.88**: Zoom sur le CH2 benzylique (2.75-3.0 ppm)
- **ArH 6.89**: Zoom sur le proton aromatique (6.75-7.0 ppm)

### Souris
- **Clic gauche + glisser**: Sélectionner une région à zoomer
- **Clic droit**: Revenir au spectre complet
- **Molette haut**: Zoom avant (centré sur curseur)
- **Molette bas**: Zoom arrière

### Clavier
- **Flèche gauche**: Déplacer vers la gauche
- **Flèche droite**: Déplacer vers la droite
- **Flèche haut**: Zoom avant
- **Flèche bas**: Zoom arrière

## Informations affichées

### Panneau d'information (coin supérieur gauche)

Pour chaque région zoomée, le visualiseur affiche:

1. **Plage ppm**: Limites de la région visible
2. **Nombre de pics détectés**
3. **Liste des pics** (δ en ppm)
4. **Espacements J** (en Hz) entre pics successifs
5. **Type de pattern**: Singulet, Doublet, Triplet, Quartet, ou Multiplet

### Annotations sur le spectre

- **Points rouges**: Pics détectés automatiquement
- **Étiquettes jaunes**: Valeurs δ précises pour chaque pic
- **Flèches rouges**: Pointent vers les pics annotés

## Exemple: Flavoglaucin

Le fichier de test `Flavauglocin_multiplicity.json` contient:

**Molécule**: Flavoglaucin (produit naturel fongique)
**Formule**: C₂₁H₂₈O₃
**SMILES**: `CCCCCCCC1=C(C=C(C(=C1C=O)O)CC=C(C)C)O`

### Structure du spectre:

| Région ppm | Groupe | Multiplicité | Protons | Notes |
|------------|--------|--------------|---------|-------|
| 0.88 | CH₃ | t | 3H | Triplet terminal |
| 1.29-1.58 | (CH₂)ₙ | m | 10-12H | Chaîne alkyle, **système second ordre AA'BB'** |
| 1.70-1.76 | CH₃ | d | 6H | Groupes méthyle allyliques |
| 2.88 | CH₂ | m | 2H | Système AA'BB' |
| 3.29 | CH₂ | d | 2H | Doublet benzylique |
| 4.31 | OH | s | 1H | Proton échangeable |
| 5.28 | =CH | ddt | 1H | Proton vinylique |
| 6.89 | Ar-H | s | 1H | Proton aromatique |
| 10.25 | CHO | s | 1H | Aldéhyde |
| 11.92 | OH | s | 1H | Phénol (échangeable) |

### Systèmes second ordre détectés:

La molécule contient plusieurs systèmes **AA'BB'** dans la chaîne alkyle:
- Détection automatique par l'algorithme
- Simulation exacte par diagonalisation Hamiltonienne
- **Fragmentation automatique** si > 10 spins

## Fragmentation Hamiltonienne

### Principe

Pour les systèmes > 10 spins, la fragmentation découpe automatiquement le système en fragments faiblement couplés.

**Critère de fragmentation**: |J/Δν| < 0.05

### Exemple avec Flavoglaucin

La chaîne alkyle CCCCCCCC contient ~16 protons:

**Sans fragmentation**: Impossible à simuler (2¹⁶ = 65536×65536 matrice)

**Avec fragmentation**:
1. Fragment 1: CH₃ terminal + CH₂ adjacent (5 spins)
2. Fragment 2: CH₂ centraux (6-8 spins)
3. Fragment 3: CH₂ + CH₂ benzylique (4 spins)

**Résultat**: Simulation exacte de chaque fragment puis combinaison

### Logs de fragmentation

Si la fragmentation est activée, vous verrez dans la console:

```
[INFO] Système de 14 spins > 10, tentative de fragmentation
[INFO] Fragmentation: 14 spins, seuil J/Δν = 0.05
[INFO] → 3 fragments détectés
[INFO]   Fragment 1: 5 spins [0, 1, 2, 3, 4]
[INFO]   Fragment 2: 6 spins [5, 6, 7, 8, 9, 10]
[INFO]   Fragment 3: 3 spins [11, 12, 13]
[INFO]   Fragment 1: 48 transitions
[INFO]   Fragment 2: 96 transitions
[INFO]   Fragment 3: 12 transitions
[INFO] Total: 156 transitions combinées
```

## Analyse de patterns

### Détection automatique

Le visualiseur détecte automatiquement:
- **Singulets** (s): 1 pic
- **Doublets** (d): 2 pics équidistants
- **Triplets** (t): 3 pics (1:2:1)
- **Quartets** (q): 4 pics (1:3:3:1)
- **Multiplets complexes** (m): > 4 pics

### Mesure des couplages J

Pour les multiplets simples (2-10 pics), le visualiseur calcule automatiquement:
- Espacement entre pics successifs (en Hz)
- Permet d'identifier les constantes de couplage ³J et ²J

## Troubleshooting

### Erreur: "No module named 'numpy'"

**Solution**: Activez l'environnement conda:
```bash
conda activate rdkit_env_clean
```

### Erreur: "Fichier non trouvé"

**Solution**: Vérifiez que le fichier JSON de multiplicité existe:
```bash
ls test_output_molecule/Flavauglocin_multiplicity.json
```

### Fragmentation désactivée

Si networkx n'est pas installé, la fragmentation sera désactivée automatiquement.

**Solution**: Installer networkx:
```bash
pip install networkx
```

### Spectre vide ou erreur de simulation

Vérifiez les logs dans la console pour identifier:
- Erreurs de parsing du JSON
- Problèmes de couplages
- Échec de la simulation Hamiltonienne

## Fichiers associés

- `interactive_spectrum_viewer.py` - Code principal du visualiseur
- `backend/scripts/sim_mult_rulebased.py` - Moteur de simulation avec fragmentation
- `backend/scripts/simulate_spectrum_api.py` - API de génération de spectres
- `backend/scripts/FRAGMENTATION_README.md` - Documentation technique de la fragmentation
- `backend/scripts/test_fragmentation_example.py` - Démonstration conceptuelle

## Références

1. **Fragmentation**: Castillo, A. M., Patiny, L., & Wist, J. (2011). *J. Magn. Reson.* **209**, 123-130
2. **Simulation quantique**: Pople, J. A., Schneider, W. G., & Bernstein, H. J. (1959). *High-resolution NMR*

## Contact

Pour toute question ou bug:
- Vérifier les logs dans la console
- Consulter `backend/scripts/FRAGMENTATION_README.md`
- Examiner le code source dans `interactive_spectrum_viewer.py`
