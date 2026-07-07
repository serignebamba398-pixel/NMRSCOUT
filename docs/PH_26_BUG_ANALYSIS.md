# Bug Analysis: pH = 26 Error

## Symptôme

```
[Flavoglaucin] ⚠️ Simulation multiplicites echouee: pH 26 hors plage valide [1.0, 14.0]
```

## Investigation

### Test Reproduction

**Test exécuté** :
```bash
cd backend/scripts
python integrated_nmr_prediction.py \
  --smiles "CCCCCCCC1=C(O)C=C(CC=C(C)C)C(O)=C1C=O" \
  --mol_id "Flavoglaucin_test" \
  --model_1h "../models/1H_sparsified_proposed_proposed_set2set_1.pt" \
  --model_13c "../models/13C_sparsified_proposed_proposed_set2set_1.pt" \
  --simulate_multiplicity \
  --simulator_path "sim_mult_rulebased.py" \
  --sim_solvent "CDCl3" \
  --sim_freq_mhz 700 \
  --sim_temperature_k 298 \
  --sim_ph 7.4 \
  --sim_enable_dynamics
```

**Résultat** : Bug reproduit avec la version actuelle de `backend/scripts/`!

### Logs de Debug

```json
"[DEBUG] args.temperature_k = 298.15"
"[DEBUG] args.sim_temperature_k = 298.0"
"[DEBUG] args.ph = 7.4"
"[DEBUG] args.sim_ph = 7.4"
"[DEBUG] sim_params = {'solvent': 'CDCl3', 'freq_mhz': 700.0, 'temperature_k': 298.15, 'use_3d': False, 'debug': False, 'enable_dynamics': True, 'ph': 7.4}"
```

✅ **Les paramètres initiaux sont corrects (pH = 7.4)**

Mais plus tard:

```
"[Flavoglaucin_test] ⚠️ Simulation multiplicites echouee: pH 26 hors plage valide [1.0, 14.0]"
```

❌ **Le pH devient 26 lors de l'appel à `DynamicsPostProcessor`**

## Analyse du Traceback

```
File "sim_mult_rulebased.py", line 14505, in simulate_multiplicity_per_proton
  J_effective = apply_dynamics_to_coupling(
File "nmr_plugins/dynamics.py", line 1316, in apply_dynamics_to_coupling
  processor = DynamicsPostProcessor(ph, solvent, temperature_k)
File "nmr_plugins/dynamics.py", line 225, in __init__
  self._validate_parameters(ph, solvent, temperature_k)
File "nmr_plugins/dynamics.py", line 237, in _validate_parameters
  raise ValueError(f"pH {ph} hors plage valide [1.0, 14.0]")
ValueError: pH 26 hors plage valide [1.0, 14.0]
```

### Appel à `apply_dynamics_to_coupling()` (ligne 14505)

```python
J_effective = apply_dynamics_to_coupling(
    mol=molH,
    h1_idx=h_idx,
    h2_idx=partner_h_idx,
    J_static=J_static,
    ph=ph,  # ← La variable pH locale
    solvent=solv.get("name", solvent),
    temperature_k=temperature_k,
    use_cache=True
)
```

✅ **L'appel utilise des arguments nommés, donc pas d'erreur d'ordre**

### Signature de `apply_dynamics_to_coupling()` (ligne 1261)

```python
def apply_dynamics_to_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    J_static: float,
    ph: float = 7.4,
    solvent: str = "D2O",
    temperature_k: float = 298.15,
    use_cache: bool = True
) -> Optional[float]:
```

✅ **Signature correcte**

### Création de `DynamicsPostProcessor` (ligne 1316)

```python
processor = DynamicsPostProcessor(ph, solvent, temperature_k)
```

✅ **Ordre correct : (ph, solvent, temperature_k)**

### Signature de `DynamicsPostProcessor.__init__` (ligne 208)

```python
def __init__(
    self,
    ph: float = 7.4,
    solvent: str = "D2O",
    temperature_k: float = 298.15,
    validate_params: bool = True
):
```

✅ **Signature correcte**

## Théorie: 298 - 272 = 26

La valeur 26 est suspecte:
- **298 K** (température passée)
- **298 - 272 = 26** (conversion Kelvin→Celsius incorrecte)

La conversion correcte est: **298 K - 273.15 = 24.85°C**

Mais 272 n'apparaît nulle part dans le code source actuel.

## Hypothèses

### Hypothèse 1: Code compilé obsolète

Le bug se produit avec la version compilée dans `dist/win-unpacked/` qui contient peut-être un ancien code bugué.

**Solution**: Rebuild l'application avec le code source actuel.

### Hypothèse 2: Variable `ph` écrasée

Quelque part dans `simulate_multiplicity_per_proton()`, la variable locale `ph` (paramètre de la fonction, valeur 7.4) est écrasée avec une valeur de 26.

**Recherche effectuée**:
- Aucune assignation `ph = ...` trouvée dans la fonction
- Les seuls usages sont `ph=ph` (passage de paramètre)

### Hypothèse 3: Ordre des arguments inversé quelque part

Si quelque part un appel passe les arguments dans le mauvais ordre:
```python
# INCORRECT (hypothétique)
processor = DynamicsPostProcessor(temperature_k, solvent, ph)
# → DynamicsPostProcessor(298.15, "CDCl3", 7.4)
# → mais interprété comme (ph=298.15, solvent="CDCl3", temperature_k=7.4)
```

Mais la recherche n'a trouvé aucun appel positionnel.

### Hypothèse 4: Cache corrompu

Le cache dans `apply_dynamics_to_coupling._cache` pourrait contenir une valeur incorrecte d'un appel précédent.

**Cache key** (ligne 1302):
```python
cache_key = f"{mol_smiles}_{h1_idx}_{h2_idx}_{J_static:.3f}_{ph:.2f}_{solvent}_{temperature_k:.2f}"
```

Si un appel précédent a créé une entrée avec `ph=26`, elle serait réutilisée.

## Recommandations

### Action immédiate

1. **Rebuild l'application** pour utiliser le code source actuel:
   ```bash
   npm run build
   ```

2. **Vider le cache** au démarrage de chaque prédiction:
   ```python
   # Dans apply_dynamics_to_coupling()
   if not hasattr(apply_dynamics_to_coupling, '_cache'):
       apply_dynamics_to_coupling._cache = {}
   # OU
   apply_dynamics_to_coupling._cache.clear()  # Vider au lieu de réutiliser
   ```

### Debug supplémentaire

Ajouter des logs dans `apply_dynamics_to_coupling()` avant la ligne 1316:

```python
# Ligne 1315
print(f"[DEBUG] apply_dynamics_to_coupling called with ph={ph}, solvent={solvent}, temperature_k={temperature_k}")

# Ligne 1316
processor = DynamicsPostProcessor(ph, solvent, temperature_k)
```

Cela permettra de voir la valeur exacte de `ph` juste avant l'erreur.

### Fix définitif

Si le cache est la cause, désactiver le cache temporairement:

```python
def apply_dynamics_to_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    J_static: float,
    ph: float = 7.4,
    solvent: str = "D2O",
    temperature_k: float = 298.15,
    use_cache: bool = False  # ← Désactiver par défaut
) -> Optional[float]:
    # ... sans cache
    processor = DynamicsPostProcessor(ph, solvent, temperature_k)
    return processor.apply_dynamics(mol, h1_idx, h2_idx, J_static)
```

## Conclusion

Le bug est reproductible mais la source exacte n'est pas identifiée dans le code source actuel. Les hypothèses les plus probables sont:

1. **Code compilé obsolète** dans `dist/` (solution: rebuild)
2. **Cache corrompu** (solution: désactiver cache ou le vider au démarrage)

**Prochaine étape**: Rebuild l'application et retester.
