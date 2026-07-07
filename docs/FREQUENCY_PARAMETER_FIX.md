# Fix: Frequency Parameter Always Defaults to 700 MHz

## Problem

L'utilisateur sélectionne 400 MHz dans l'interface, mais le JSON généré montre toujours `"freq_MHz": 700.0`.

**Symptôme** : Le fichier `Flavoglaucin_multiplicity.json` contient :
```json
{
  "freq_MHz": 700.0,
  "solvent": "CDCl3",
  ...
}
```
Même quand l'utilisateur a sélectionné 400 MHz dans l'interface.

## Root Cause

Les paramètres expérimentaux (`freq_mhz`, `solvent`, `temperature_k`, `ph`) n'étaient **PAS définis** dans `DEFAULT_PREDICTION_OPTIONS`.

**Conséquence** :
- Au premier chargement, `options.freq_mhz` est `undefined`
- Le fallback `|| 700` dans le code frontend s'active systématiquement
- La valeur 400 MHz sélectionnée par l'utilisateur n'est jamais persistée correctement

## Solution

### Fichier modifié : `frontend/src/constants/defaults.js`

**Avant** :
```javascript
export const DEFAULT_PREDICTION_OPTIONS = {
  nb_cores: 2,
  extractors: { ... },
  cosy: { min_bonds: 2, max_bonds: 3 },
  idColumn: "id",
  smilesColumn: "smiles",
  // ❌ Pas de freq_mhz, solvent, temperature_k, ph
};
```

**Après** :
```javascript
export const DEFAULT_PREDICTION_OPTIONS = {
  nb_cores: 2,
  extractors: { ... },
  cosy: { min_bonds: 2, max_bonds: 3 },
  idColumn: "id",
  smilesColumn: "smiles",
  // ✅ Experimental parameters for NMR prediction
  freq_mhz: 700,
  solvent: "CDCl3",
  temperature_k: 298,
  ph: 7.4,
};
```

## Comportement attendu après le fix

### Scénario 1 : Premier chargement
- **Avant** : `freq_mhz` = `undefined` → fallback à 700 MHz
- **Après** : `freq_mhz` = `700` (valeur par défaut explicite)

### Scénario 2 : Utilisateur change à 400 MHz
- **Avant** : Valeur changée en local, mais pas persistée → revient à 700 au reload
- **Après** : Valeur persistée via `updateOption('freq_mhz', 400)` → reste 400

### Scénario 3 : Prédiction lancée avec 400 MHz
- **Avant** : Backend reçoit 700 (fallback)
- **Après** : Backend reçoit 400 (valeur réelle)

## Vérification

Pour vérifier que le fix fonctionne :

1. Ouvrir l'interface de prédiction
2. Changer la fréquence à **400 MHz**
3. Lancer une prédiction
4. Ouvrir le fichier `*_multiplicity.json` généré
5. Vérifier que `"freq_MHz": 400.0` ✓

## Impact sur les seuils J

Avec ce fix, les seuils de fréquence seront correctement appliqués :

| Paramètre | @ 400 MHz | @ 700 MHz |
|-----------|-----------|-----------|
| min_resolvable | 1.00 Hz | 0.80 Hz |
| 4J_visible | **1.00 Hz** | **0.50 Hz** |
| 5J_visible | 2.00 Hz | 0.80 Hz |

**Exemple concret** :
- Un couplage ⁴J de 0.6 Hz sera :
  - **Filtré** @ 400 MHz (sous le seuil 1.0 Hz)
  - **Visible** @ 700 MHz (au-dessus du seuil 0.5 Hz)

## Files Modified

1. `frontend/src/constants/defaults.js` - Ajout des paramètres expérimentaux par défaut

## Related Documentation

- [SOLVENT_FILTERING_FIX.md](./SOLVENT_FILTERING_FIX.md) - Fix pour le filtrage des protons échangeables selon le solvant
- [METHODOLOGY_NMR_MULTIPLICITY.md](./METHODOLOGY_NMR_MULTIPLICITY.md) - Méthodologie complète des seuils de fréquence

## Test Procedure

```bash
# 1. Rebuild frontend with new defaults
cd frontend
npm run build

# 2. Test avec molecule
# Interface: Sélectionner 400 MHz
# Vérifier JSON: "freq_MHz": 400.0

# 3. Comparer 400 vs 700 MHz
python backend/scripts/test_frequency_comparison_simple.py
```

## Notes

- Les valeurs par défaut (700 MHz, CDCl3, 298K, pH 7.4) représentent des conditions standard de laboratoire
- L'utilisateur peut toujours changer ces valeurs via l'interface
- Les changements sont maintenant correctement persistés et transmis au backend
