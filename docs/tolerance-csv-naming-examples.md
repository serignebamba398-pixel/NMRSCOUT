# Exemples de Nommage des Fichiers CSV d'Optimisation de Tolérance

## 📋 Vue d'ensemble

Chaque test d'optimisation de tolérance génère maintenant **son propre fichier CSV** avec la valeur de tolérance dans le nom.

## 📁 Exemples concrets

### Optimisation 1H avec 10 tolérances (0.05 → 0.50, step 0.05)

```
outputs/67f0c414-4601-41b7-b513-537d6719673a/optimization/
├── annotation_1H_tol_0_05.csv    ← Test avec tol=0.05 ppm
├── annotation_1H_tol_0_10.csv    ← Test avec tol=0.10 ppm
├── annotation_1H_tol_0_15.csv    ← Test avec tol=0.15 ppm
├── annotation_1H_tol_0_20.csv    ← Test avec tol=0.20 ppm
├── annotation_1H_tol_0_25.csv
├── annotation_1H_tol_0_30.csv
├── annotation_1H_tol_0_35.csv
├── annotation_1H_tol_0_40.csv
├── annotation_1H_tol_0_45.csv
└── annotation_1H_tol_0_50.csv    ← Test avec tol=0.50 ppm
```

**Total : 10 fichiers CSV** (un par tolérance testée)

### Optimisation 13C avec 10 tolérances (0.5 → 5.0, step 0.5)

```
outputs/.../optimization/
├── annotation_13C_tol_0_5.csv    ← Test avec tol=0.5 ppm
├── annotation_13C_tol_1_0.csv    ← Test avec tol=1.0 ppm
├── annotation_13C_tol_1_5.csv
├── annotation_13C_tol_2_0.csv
├── annotation_13C_tol_2_5.csv
├── annotation_13C_tol_3_0.csv
├── annotation_13C_tol_3_5.csv
├── annotation_13C_tol_4_0.csv
├── annotation_13C_tol_4_5.csv
└── annotation_13C_tol_5_0.csv    ← Test avec tol=5.0 ppm
```

### Optimisation COSY avec 25 combinaisons (F1: 0.05→0.5, F2: 0.05→0.5)

```
outputs/.../optimization/
├── annotation_COSY_tol_f1_0_05_f2_0_05.csv   ← F1=0.05, F2=0.05
├── annotation_COSY_tol_f1_0_05_f2_0_1.csv    ← F1=0.05, F2=0.10
├── annotation_COSY_tol_f1_0_05_f2_0_15.csv
├── annotation_COSY_tol_f1_0_05_f2_0_2.csv
├── annotation_COSY_tol_f1_0_05_f2_0_25.csv
├── annotation_COSY_tol_f1_0_1_f2_0_05.csv    ← F1=0.10, F2=0.05
├── annotation_COSY_tol_f1_0_1_f2_0_1.csv     ← F1=0.10, F2=0.10
...
└── annotation_COSY_tol_f1_0_25_f2_0_25.csv   ← F1=0.25, F2=0.25
```

**Total : 25 fichiers CSV** (5×5 combinaisons)

### Optimisation HSQC avec 25 combinaisons (F1: 0.5→5.0, F2: 0.05→0.5)

```
outputs/.../optimization/
├── annotation_HSQC_tol_f1_0_5_f2_0_05.csv    ← F1=0.5, F2=0.05
├── annotation_HSQC_tol_f1_0_5_f2_0_1.csv     ← F1=0.5, F2=0.10
├── annotation_HSQC_tol_f1_1_0_f2_0_05.csv    ← F1=1.0, F2=0.05
├── annotation_HSQC_tol_f1_1_0_f2_0_1.csv
...
└── annotation_HSQC_tol_f1_5_0_f2_0_5.csv     ← F1=5.0, F2=0.5
```

## 🔍 Comment retrouver un fichier CSV spécifique ?

### Méthode 1 : Par nom de fichier

Si tu cherches les résultats pour **1H avec tolérance 0.15 ppm** :
```bash
outputs/{taskId}/optimization/annotation_1H_tol_0_15.csv
```

### Méthode 2 : Par pattern (glob)

Tous les tests 1H :
```bash
outputs/{taskId}/optimization/annotation_1H_tol_*.csv
```

Tous les tests COSY avec F1=0.05 :
```bash
outputs/{taskId}/optimization/annotation_COSY_tol_f1_0_05_*.csv
```

### Méthode 3 : Via le JSON résumé

Le fichier `tolerance_results.json` contient le mapping entre tolérance et résultats :

```json
{
  "toleranceParams": {"tol_correspondence": 0.15},
  "totalSignalsMatched": 45,
  "avgCoverage": 78.2
}
```

Tu peux alors ouvrir le CSV correspondant : `annotation_1H_tol_0_15.csv`

## 📊 Contenu des fichiers CSV

Chaque fichier CSV contient **l'annotation complète** avec cette tolérance :

```csv
Molecule,SMILES,InChIKey,Coverage_RANSY,RANSY_signals,Matched_signals,...
mol001,CCO,LFQSCWFLJHTTHZ-UHFFFAOYSA-N,85.5,"1.2|2.4|3.6","1.23|2.45|3.67",...
mol002,CC(C)O,IPFXBQRQXJSAJR-UHFFFAOYSA-N,72.3,"0.9|1.5","0.92|1.54",...
...
```

## 🆚 Comparaison avant/après

### ❌ Avant (écrasement)

```
optimization/
└── annotation_1H_latest.csv  ← Écrasé à chaque test, on perd les précédents
```

**Problème :** Impossible de comparer les résultats de plusieurs tolérances

### ✅ Après (conservation)

```
optimization/
├── annotation_1H_tol_0_05.csv  ← Tous les tests sont conservés
├── annotation_1H_tol_0_10.csv
├── annotation_1H_tol_0_15.csv
└── ...
```

**Avantage :** Tous les CSV sont disponibles pour analyse/comparaison

## 💡 Cas d'usage

### Analyser un CSV spécifique

```bash
# Ouvrir le CSV avec Excel/LibreOffice
open outputs/{taskId}/optimization/annotation_1H_tol_0_15.csv

# Compter les lignes (molécules)
wc -l outputs/{taskId}/optimization/annotation_1H_tol_0_15.csv

# Extraire les molécules avec coverage > 80%
awk -F',' '$4 > 80' outputs/{taskId}/optimization/annotation_1H_tol_0_15.csv
```

### Comparer deux tolérances

```bash
# Comparer les tailles de fichiers
ls -lh outputs/{taskId}/optimization/annotation_1H_tol_0_05.csv
ls -lh outputs/{taskId}/optimization/annotation_1H_tol_0_15.csv

# Compter les molécules avec coverage > 80% pour chaque tolérance
for f in annotation_1H_tol_*.csv; do
  echo "$f: $(awk -F',' '$4 > 80' $f | wc -l) molécules >80%"
done
```

### Trouver la meilleure tolérance

1. Ouvrir `tolerance_results.json`
2. Chercher le test avec `totalSignalsMatched` maximum
3. Noter la valeur de `tol_correspondence`
4. Ouvrir le CSV correspondant pour analyse détaillée

## 🗑️ Nettoyage

Pour supprimer tous les CSV d'optimisation 1H :
```bash
rm outputs/{taskId}/optimization/annotation_1H_tol_*.csv
```

Pour supprimer tous les CSV d'optimisation :
```bash
rm outputs/{taskId}/optimization/annotation_*_tol_*.csv
```

## 📏 Espace disque

**Estimation :** Un CSV d'annotation typique = 5-50 MB

- 10 tests 1H → 50-500 MB
- 25 tests COSY → 125 MB - 1.25 GB
- Tous les tests (1H + 13C + COSY + HSQC) → ~2-5 GB

**Recommandation :** Nettoyer périodiquement les anciens CSV si l'espace disque est limité.
