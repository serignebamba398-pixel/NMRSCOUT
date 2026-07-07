# Stockage des Résultats d'Optimisation de Tolérance

## Vue d'ensemble

Les résultats d'optimisation de tolérance sont maintenant **persistés sur disque** dans le dossier `outputs/{taskId}/tolerance_results.json`.

Auparavant, les résultats étaient uniquement stockés en mémoire dans le state React et étaient perdus lors d'un refresh de la page.

## Structure de fichier

### Emplacement
```
outputs/
  └── {taskId}/
      ├── tolerance_results.json  ← Résumé JSON de TOUS les tests
      ├── annotations/
      │   ├── annotation_1H_latest.csv      ← Annotations normales (pas optimization)
      │   └── history/                      ← Historique des annotations
      ├── optimization/
      │   ├── annotation_1H_tol_0_05.csv    ← NOUVEAU : Un CSV par tolérance testée
      │   ├── annotation_1H_tol_0_10.csv
      │   ├── annotation_1H_tol_0_15.csv
      │   ├── annotation_13C_tol_0_5.csv
      │   ├── annotation_13C_tol_1_0.csv
      │   ├── annotation_COSY_tol_f1_0_05_f2_0_1.csv   ← 2D : F1 et F2
      │   ├── annotation_HSQC_tol_f1_0_5_f2_0_05.csv
      │   └── temp_1H_2026-02-11T08-04-07/  ← Dossier temporaire
      │       ├── params.json
      │       ├── matrix.csv
      │       └── checkpoint.rds
      └── ... (autres fichiers)
```

## 📝 Schéma de nommage des fichiers CSV

### Méthodes 1D (1H, 13C)

Format : `annotation_{type}_tol_{tolerance}.csv`

**Exemples :**
- `annotation_1H_tol_0_05.csv` → 1H avec tolérance 0.05 ppm
- `annotation_1H_tol_0_15.csv` → 1H avec tolérance 0.15 ppm
- `annotation_13C_tol_0_5.csv` → 13C avec tolérance 0.5 ppm
- `annotation_13C_tol_2_5.csv` → 13C avec tolérance 2.5 ppm

**Note :** Le point décimal est remplacé par underscore (`_`) dans le nom de fichier.

### Méthodes 2D (COSY, HSQC)

Format : `annotation_{type}_tol_f1_{tol_f1}_f2_{tol_f2}.csv`

**Exemples :**
- `annotation_COSY_tol_f1_0_05_f2_0_1.csv` → COSY avec F1=0.05, F2=0.10 ppm
- `annotation_HSQC_tol_f1_0_5_f2_0_05.csv` → HSQC avec F1=0.5, F2=0.05 ppm

### Annotations normales (non-optimization)

Format : `annotation_{type}_latest.csv`

**Exemples :**
- `annotation_1H_latest.csv` → Dernière annotation 1H normale
- `annotation_13C_latest.csv` → Dernière annotation 13C normale

**Comportement :** Les annotations normales sont archivées dans `history/` avant d'être écrasées.

### Format JSON

```json
[
  {
    "id": "1H_1739270647000",
    "type": "1H",
    "timestamp": "2026-02-11T08:04:07.000Z",
    "label": "1H - 2/11/2026, 8:04:07 AM",
    "results": [
      {
        "toleranceParams": {
          "tol_correspondence": 0.05,
          "tol_correspondance": 0.05,
          "tol_match": 0.05
        },
        "coverageDistribution": [10, 20, 30, 25, 15],
        "totalSignalsMatched": 42,
        "totalExperimentalSignals": 50,
        "signalsMatchedPct": 84.0,
        "totalMolecules": 125,
        "avgCoverage": 75.5,
        "matchedSignals": [1.23, 2.45, 3.67, 4.89]
      }
    ]
  },
  {
    "id": "13C_1739280000000",
    "type": "13C",
    "timestamp": "2026-02-11T14:53:20.000Z",
    "label": "13C - 2/11/2026, 2:53:20 PM",
    "results": [ /* ... */ ]
  }
]
```

## API Backend

### Routes ajoutées

#### `POST /api/tasks/:taskId/tolerance-results`
Sauvegarde les résultats d'optimisation de tolérance.

**Request Body:**
```json
{
  "results": [ /* array of result sets */ ]
}
```

**Response:**
```json
{
  "message": "Tolerance results saved successfully",
  "taskId": "67f0c414-4601-41b7-b513-537d6719673a",
  "count": 2,
  "path": "C:\\...\\outputs\\67f0c414...\\tolerance_results.json"
}
```

#### `GET /api/tasks/:taskId/tolerance-results`
Récupère les résultats d'optimisation sauvegardés.

**Response:**
```json
{
  "taskId": "67f0c414-4601-41b7-b513-537d6719673a",
  "results": [ /* array of result sets */ ],
  "count": 2
}
```

Si aucun fichier n'existe :
```json
{
  "taskId": "67f0c414-4601-41b7-b513-537d6719673a",
  "results": [],
  "message": "No saved tolerance results found"
}
```

#### `DELETE /api/tasks/:taskId/tolerance-results`
Supprime les résultats sauvegardés.

**Response:**
```json
{
  "message": "Tolerance results deleted successfully",
  "taskId": "67f0c414-4601-41b7-b513-537d6719673a"
}
```

## API Frontend

### Nouvelles fonctions dans `src/api.js`

```javascript
// Sauvegarder les résultats
await saveToleranceResults(taskId, results);

// Charger les résultats
const { results, count } = await loadToleranceResults(taskId);

// Supprimer les résultats
await deleteToleranceResults(taskId);
```

## Comportement

### Sauvegarde automatique
- **Quand** : Après chaque optimisation de tolérance terminée avec succès
- **Où** : Dans `useToleranceOptimization.js`, ligne ~495
- **Données** : Tous les résultats d'optimisation du task en cours (historique complet)

```javascript
const updatedResults = [...annotationState.toleranceResults, newResult];
await saveToleranceResults(annotationState.selectedTask, updatedResults);
```

### Chargement automatique
- **Quand** : Lors de la sélection d'un task dans l'interface d'annotation
- **Où** : Hook `useToleranceOptimization`, `useEffect` sur `selectedTask`
- **Effet** : Les résultats sauvegardés sont automatiquement chargés dans `annotationState.toleranceResults`

```javascript
useEffect(() => {
  if (annotationState.selectedTask) {
    loadSavedToleranceResults(annotationState.selectedTask);
  }
}, [annotationState.selectedTask]);
```

### Logs
- **Sauvegarde** : `[OK] Tolerance results saved to outputs/{taskId}/tolerance_results.json`
- **Chargement** : `[OK] Loaded {count} saved tolerance optimization result(s) from task {taskId}`
- **Erreur** : `[WARN] Failed to save/load tolerance results: {error}`

## Cas d'usage

### Scénario 1 : Optimisation simple
1. Utilisateur lance une optimisation 1H avec 10 valeurs de tolérance
2. Les résultats sont ajoutés à `toleranceResults[]` en mémoire
3. Après la fin, sauvegarde automatique dans `tolerance_results.json`
4. Utilisateur ferme l'application
5. Au prochain démarrage, sélection du même task → résultats rechargés automatiquement

### Scénario 2 : Optimisations multiples
1. Optimisation 1H → sauvegardée
2. Optimisation 13C → sauvegardée (écrase avec les 2 résultats)
3. Optimisation COSY → sauvegardée (écrase avec les 3 résultats)
4. Le fichier contient maintenant 3 sets de résultats

### Scénario 3 : Changement de task
1. Task A sélectionné → chargement des résultats de A
2. Task B sélectionné → chargement des résultats de B (remplace en mémoire)
3. Chaque task a son propre fichier `tolerance_results.json`

## Fichiers modifiés

### Backend
- `backend/routes/tasks.js` : Ajout des 3 routes POST/GET/DELETE pour `/api/tasks/:taskId/tolerance-results`
- `backend/services/annotationService.js` :
  - Ajout de la fonction `getToleranceStringForFilename()` (ligne ~21)
  - Modification du nom de fichier CSV pour inclure la tolérance (ligne ~720)
  - Les fichiers CSV d'optimisation sont maintenant nommés avec la tolérance
  - Les annotations normales continuent d'utiliser `_latest.csv`

### Frontend
- `frontend/src/api.js` : Ajout des fonctions `saveToleranceResults()`, `loadToleranceResults()`, `deleteToleranceResults()`
- `frontend/src/hooks/useToleranceOptimization.js` :
  - Import des nouvelles fonctions API
  - Ajout de `loadSavedToleranceResults()` callback
  - Sauvegarde automatique après optimisation (ligne ~495)
  - Chargement automatique via `useEffect` sur `selectedTask`
  - Export de `loadSavedToleranceResults` dans le return

## Fonctionnement technique

### Génération du nom de fichier CSV

La fonction `getToleranceStringForFilename()` dans `annotationService.js` :

```javascript
function getToleranceStringForFilename(annotationType, params) {
  if (annotationType === '1H' || annotationType === '13C') {
    const tol = params.tol_correspondence || params.tol_correspondance || 0;
    return `tol_${tol}`.replace(/\./g, '_');  // Remplace . par _
  } else if (annotationType === 'COSY' || annotationType === 'HSQC') {
    const tolF1 = params.tol_correspondence_f1 || params.tol_correspondance_F1 || 0;
    const tolF2 = params.tol_correspondence_f2 || params.tol_correspondance_F2 || 0;
    return `tol_f1_${tolF1}_f2_${tolF2}`.replace(/\./g, '_');
  }
  return 'tol_unknown';
}
```

**Logique :**
1. Extrait la tolérance des paramètres (`tol_correspondence` ou variantes)
2. Construit une chaîne avec la valeur de tolérance
3. Remplace les points décimaux par underscores (pour compatibilité filesystem)
4. Retourne la chaîne à insérer dans le nom de fichier

### Distinction optimization vs annotation normale

Dans `annotationService.js` (ligne ~720) :

```javascript
let latestFilename;
if (isOptimization) {
  // Optimisation : nom avec tolérance
  const tolStr = getToleranceStringForFilename(annotationType, optimizedParams);
  latestFilename = `annotation_${annotationType}_${tolStr}.csv`;
} else {
  // Annotation normale : nom standard
  latestFilename = `annotation_${annotationType}_latest.csv`;
}
```

**Comportement :**
- `isOptimization = true` → Fichier sauvegardé dans `outputs/{taskId}/optimization/annotation_1H_tol_0_05.csv`
- `isOptimization = false` → Fichier sauvegardé dans `outputs/{taskId}/annotations/annotation_1H_latest.csv`

## Test manuel

### Test de sauvegarde
```bash
curl -X POST "http://localhost:5000/api/tasks/67f0c414-4601-41b7-b513-537d6719673a/tolerance-results" \
  -H "Content-Type: application/json" \
  -d '{"results":[{"id":"test","type":"1H","results":[]}]}'
```

### Test de chargement
```bash
curl -X GET "http://localhost:5000/api/tasks/67f0c414-4601-41b7-b513-537d6719673a/tolerance-results"
```

### Test de suppression
```bash
curl -X DELETE "http://localhost:5000/api/tasks/67f0c414-4601-41b7-b513-537d6719673a/tolerance-results"
```

## Notes importantes

1. **Backward compatible** : Si le fichier n'existe pas, l'API retourne `results: []` sans erreur
2. **Écrasement** : Chaque sauvegarde écrase le fichier précédent (pas d'append)
3. **Pas de versioning** : Un seul fichier par task, contenant tous les résultats
4. **Format JSON** : Fichier formaté avec 2 espaces pour lisibilité
5. **Sécurité** : Utilisation de `safeResolve()` pour éviter les attaques path traversal

## Améliorations futures possibles

- [ ] Ajouter un bouton "Export" pour télécharger les résultats en CSV
- [ ] Ajouter un bouton "Clear all results" pour vider l'historique
- [ ] Ajouter une notification toast lors de la sauvegarde/chargement
- [ ] Versioning des résultats (fichiers datés)
- [ ] Compression des fichiers JSON pour économiser l'espace disque
- [ ] Sauvegarde incrémentale (append au lieu d'écraser)
