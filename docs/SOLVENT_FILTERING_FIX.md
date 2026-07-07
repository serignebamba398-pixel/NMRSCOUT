# Fix: Solvent-Based Filtering of Exchangeable Protons

## Problem

Les protons échangeables (OH, NH, SH) apparaissaient toujours dans les spectres prédits, même en solvants protiques (D2O, CD3OD) où ils devraient être invisibles à cause de l'échange rapide avec le solvant.

**Exemple**: En D2O, les protons OH et NH2 d'une molécule comme la glycine (NH2-CH2-COOH) devraient être invisibles, mais ils apparaissaient quand même dans le SDF et le JSON générés.

## Root Causes

1. **`predictionService.js` ne passait pas `--sim_enable_dynamics`**
   - Le flag pour activer le post-processeur dynamique n'était jamais passé au script Python
   - Résultat: Le module `dynamics.py` n'était jamais appelé

2. **`_build_x_h_connectivity_from_mol_with_H()` incluait tous les H**
   - La fonction construisait la connectivité H sans filtrer selon le solvant
   - Tous les protons (CH, OH, NH, SH) étaient systématiquement inclus

3. **Pas de filtrage basé sur le solvant**
   - Aucune logique pour distinguer solvants protiques (D2O, CD3OD) vs aprotiques (CDCl3)
   - Les règles physico-chimiques (échange rapide en D2O) n'étaient pas appliquées

## Solution Implemented

### 1. Backend Service (predictionService.js)

**Fichier**: `backend/services/predictionService.js`

**Changement**: Ajout du flag `--sim_enable_dynamics` pour les solvants protiques

```javascript
// Enable dynamics processor for protic solvents (D2O, CD3OD, CD3COCD3)
// This handles rapid proton exchange (OH, NH, SH become invisible)
const proticSolvents = ['D2O', 'd2o', 'CD3OD', 'cd3od', 'MeOD', 'meod', 'CD3COCD3', 'cd3cocd3'];
if (proticSolvents.includes(solvent)) {
  args.push('--sim_enable_dynamics');
}
```

**Impact**: Le module `dynamics.py` est maintenant activé en solvants protiques.

### 2. Python Integration Script (integrated_nmr_prediction.py)

**Fichier**: `backend/scripts/integrated_nmr_prediction.py`

#### A. Modification de `_build_x_h_connectivity_from_mol_with_H()`

**Avant**:
```python
def _build_x_h_connectivity_from_mol_with_H(molH: Chem.Mol) -> Dict[int, List[int]]:
    # Incluait TOUS les H sans distinction
    for a in molH.GetAtoms():
        if a.GetSymbol() != "H":
            continue
        # ...inclure tous les H
```

**Après**:
```python
def _build_x_h_connectivity_from_mol_with_H(molH: Chem.Mol, solvent: str = "CDCl3") -> Dict[int, List[int]]:
    # Solvants protiques: protons échangeables (OH, NH, SH) ne sont pas visibles
    protic_solvents = {'D2O', 'd2o', 'CD3OD', 'cd3od', 'MeOD', 'meod', 'CD3COCD3', 'cd3cocd3'}
    exclude_exchangeable = solvent in protic_solvents

    for a in molH.GetAtoms():
        if a.GetSymbol() != "H":
            continue
        # ...
        heavy_symbol = heavy_atom.GetSymbol()

        # Filtrer les protons échangeables en solvants protiques
        if exclude_exchangeable and heavy_symbol in ('O', 'N', 'S'):
            # OH, NH, SH ne sont pas visibles en D2O/CD3OD (échange rapide)
            continue

        # ...inclure seulement les H visibles
```

#### B. Ajout du logging

```python
# Log proton filtering based on solvent
total_h = sum(1 for a in molH_for_sim.GetAtoms() if a.GetSymbol() == 'H')
included_h = sum(len(h_list) for h_list in x_h_conn_sim.values())
if total_h > included_h:
    send_log(f"Solvent {solvent_for_sim}: Excluded {total_h - included_h} exchangeable protons (OH/NH/SH)")
```

#### C. Appel avec le solvant

```python
solvent_for_sim = sp.get("solvent", "CDCl3")
x_h_conn_sim = _build_x_h_connectivity_from_mol_with_H(molH_for_sim, solvent=solvent_for_sim)
```

## Behavior After Fix

### CDCl3 (Aprotic Solvent)
- **Tous les protons sont visibles**: CH, OH, NH, SH
- **Exemple (éthanol CH3-CH2-OH)**: 6 protons visibles (3 CH3 + 2 CH2 + 1 OH)
- **enable_dynamics**: `False` (pas nécessaire)

### D2O (Protic Solvent)
- **Protons CH visibles**, protons OH/NH/SH **invisibles** (échange rapide)
- **Exemple (éthanol CH3-CH2-OH)**: 5 protons visibles (3 CH3 + 2 CH2), OH exclu
- **enable_dynamics**: `True` (module dynamics.py activé)

### CD3OD (Protic Solvent)
- **Même comportement que D2O**
- **Exemple (glycine NH2-CH2-COOH)**: 2 protons visibles (2 CH2), NH2 et COOH exclus
- **enable_dynamics**: `True`

## Testing

Un script de test a été créé: `backend/scripts/test_solvent_filtering.py`

**Tests effectués**:
1. Éthanol en CDCl3: 6 H visibles ✓
2. Éthanol en D2O: 5 H visibles (OH exclu) ✓
3. Glycine en CDCl3: 5 H visibles ✓
4. Glycine en D2O: 2 H visibles (NH2 et COOH exclus) ✓

## Files Modified

1. `backend/services/predictionService.js` (lignes 512-522)
2. `backend/scripts/integrated_nmr_prediction.py`:
   - Fonction `_build_x_h_connectivity_from_mol_with_H()` (lignes 1052-1091)
   - Appel avec solvant (lignes 1513-1522)

## Related Files

- `backend/scripts/nmr_plugins/dynamics.py`: Post-processeur pour effets dynamiques (pH, échange)
- `backend/scripts/sim_mult_rulebased.py`: Simulateur de multiplicités (utilise les H filtrés)
- `backend/scripts/fingerprint_db_predictor.py`: Prédicteur basé sur fingerprints (reçoit le solvant)

## Impact on User Experience

**Avant**:
- Spectre en D2O montrait des pics OH/NH qui n'existent pas expérimentalement
- Confusion pour l'utilisateur: "Pourquoi je vois un OH en D2O?"

**Après**:
- Spectre en D2O ne montre QUE les protons non-échangeables (CH)
- Comportement conforme à la réalité expérimentale
- Message de log informatif: "Solvent D2O: Excluded 3 exchangeable protons (OH/NH/SH)"

## Physical Chemistry Background

### Échange rapide en solvants protiques

En D2O ou CD3OD, les protons labiles (OH, NH, SH) échangent rapidement avec le solvant deutéré:

```
R-OH + D2O ⇌ R-OD + HOD (très rapide, k >> 1000 s⁻¹)
R-NH2 + D2O ⇌ R-ND2 + HOD (rapide, k > 100 s⁻¹)
```

**Conséquence RMN**:
- L'échange est plus rapide que la différence de fréquence Δν entre les sites
- Le signal est moyenné et disparaît (ou apparaît comme un pic large du solvant résiduel)
- Seuls les protons C-H restent visibles

### Échange lent en solvants aprotiques

En CDCl3, pas d'échange rapide:
- Les protons OH, NH, SH donnent des signaux nets avec couplages J visibles
- Par exemple: OH d'un alcool en CDCl3 = triplet si voisin d'un CH2

## References

- Claridge, T. D. W. *High-Resolution NMR Techniques in Organic Chemistry*, 3rd ed., Elsevier, 2016.
- Section on "Exchangeable Protons in Protic Solvents"
