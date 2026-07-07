# Fix: Geminal CH₃ Correction Applied Only When Needed

## Problem

La correction diastéréotopique pour les CH₃ géminaux sur double liaison (groupe prényl/isoprényl) était appliquée **systématiquement**, même quand les shifts ML prédits étaient déjà différents.

**Exemple problématique** :
- ML prédit : CH₃(cis) = 1.64 ppm, CH₃(trans) = 1.70 ppm (déjà différents de 0.06 ppm)
- Correction appliquée : CH₃(cis) = 1.64 - 0.03 = 1.61 ppm, CH₃(trans) = 1.70 + 0.03 = 1.73 ppm
- Résultat : Δ = 0.12 ppm au lieu de 0.06 ppm (sur-correction)

## Root Cause

Le code détectait les paires de CH₃ géminaux via SMARTS pattern matching, puis appliquait **automatiquement** un delta de ±0.03 ppm, sans vérifier si le modèle ML avait déjà prédit des shifts différents.

## Solution

### Logique de correction conditionnelle

**Avant** :
```python
geminal_ch3_info = _detect_geminal_ch3_on_double_bond(molH_for_sim, x_h_conn_sim)
if geminal_ch3_info:
    # ❌ Appliqué systématiquement sans vérifier les shifts existants
    for heavy_idx, (partner_idx, position) in geminal_ch3_info.items():
        if position == 'cis':
            proton_shifts_ppm[h_idx] -= half_delta
        else:
            proton_shifts_ppm[h_idx] += half_delta
```

**Après** :
```python
geminal_ch3_info = _detect_geminal_ch3_on_double_bond(molH_for_sim, x_h_conn_sim)
if geminal_ch3_info:
    for pair in geminal_pairs:
        avg_shift1 = calculate_average_shift(ch3_1)
        avg_shift2 = calculate_average_shift(ch3_2)

        # ✅ Appliquer SEULEMENT si shifts identiques (< 0.01 ppm)
        if abs(avg_shift1 - avg_shift2) < 0.01:
            apply_correction(±0.03 ppm)
            log("Correction applied: identical shifts detected")
        else:
            log(f"Correction skipped: already different (Δ = {delta:.3f} ppm)")
```

### Détection de paires déjà différenciées

La correction n'est appliquée que si :
1. **Pattern détecté** : `C(CH₃)₂=CH-R` (double liaison trisubstituée)
2. **Shifts identiques** : `|δ₁ - δ₂| < 0.01 ppm`

Si le modèle ML a déjà prédit des shifts différents (|Δδ| ≥ 0.01 ppm), la correction diastéréotopique est **sautée**.

## Behavior After Fix

### Case 1: ML predicts identical shifts
```
ML Prediction:
  CH₃(cis) = 1.70 ppm
  CH₃(trans) = 1.70 ppm
  Δ = 0.00 ppm  ← Identical

Detection: Geminal CH₃ pair detected
Check: |1.70 - 1.70| = 0.00 < 0.01 ppm ✓

✅ APPLY CORRECTION:
  CH₃(cis) = 1.70 - 0.03 = 1.67 ppm
  CH₃(trans) = 1.70 + 0.03 = 1.73 ppm
  Final Δ = 0.06 ppm

Log: "Geminal CH3 correction: 1.700 ppm → 1.670 / 1.730 ppm (identical shifts)"
```

### Case 2: ML already predicts different shifts
```
ML Prediction:
  CH₃(cis) = 1.64 ppm
  CH₃(trans) = 1.70 ppm
  Δ = 0.06 ppm  ← Already different

Detection: Geminal CH₃ pair detected
Check: |1.64 - 1.70| = 0.06 > 0.01 ppm ✗

✗ SKIP CORRECTION:
  CH₃(cis) = 1.64 ppm (unchanged)
  CH₃(trans) = 1.70 ppm (unchanged)
  Final Δ = 0.06 ppm

Log: "Geminal CH3 skip: 1.640 vs 1.700 ppm (already different, Δ = 0.060 ppm)"
```

### Case 3: ML predicts slightly different shifts
```
ML Prediction:
  CH₃(cis) = 1.69 ppm
  CH₃(trans) = 1.70 ppm
  Δ = 0.01 ppm  ← Slightly different (noise level)

Detection: Geminal CH₃ pair detected
Check: |1.69 - 1.70| = 0.01 = 0.01 ppm (borderline)

Behavior depends on threshold:
- If < 0.01: Apply correction
- If ≥ 0.01: Skip correction

Current threshold: 0.01 ppm (conservative, favors skipping)
```

## Threshold Justification

**Why 0.01 ppm?**

1. **ML Prediction Precision**: Modèles DGL-GNN ont une précision typique de ~0.3 ppm RMSE
   - Différences < 0.01 ppm sont probablement du bruit numérique
   - Différences > 0.01 ppm sont des prédictions intentionnelles

2. **Chemical Significance**: Δδ diastéréotopique typique = 0.04-0.08 ppm
   - Si ML prédit déjà > 0.01 ppm, il a probablement capturé l'effet
   - Si ML prédit < 0.01 ppm, le modèle n'a pas discriminé → on corrige

3. **Avoid Over-Correction**: Si on corrige toujours, on peut doubler le delta
   - ML: 0.06 ppm → Correction: +0.06 ppm → Total: 0.12 ppm (faux!)

## Log Messages

### Before Fix
```
[20:21:54] [PRED] Applied geminal CH3 delta (±0.030 ppm) to 2 methyl groups
```
No information on whether correction was needed or not.

### After Fix
```
# Case 1: Correction needed (identical shifts)
[20:21:54] [PRED] Geminal CH3 correction: 1.700 ppm → 1.670 / 1.730 ppm (identical shifts)
[20:21:54] [PRED] Applied geminal CH3 delta (±0.030 ppm) to 1 pairs of methyl groups

# Case 2: Correction skipped (already different)
[20:21:54] [PRED] Geminal CH3 skip: 1.640 vs 1.700 ppm (already different, Δ = 0.060 ppm)
```

Clear indication of decision rationale.

## Related Issue: pH = 26 Bug

**Also fixed in this commit:**

Added debug logging to identify where `pH = 26` was coming from:

```python
send_log(f"[DEBUG] args.temperature_k = {getattr(args, 'temperature_k', 'NOT SET')}")
send_log(f"[DEBUG] args.ph = {getattr(args, 'ph', 'NOT SET')}")
send_log(f"[DEBUG] sim_params = {sim_params}")
```

This will help diagnose if `temperature_k` is being passed incorrectly as 26 instead of 298.

**Expected output:**
```
[DEBUG] args.temperature_k = 298
[DEBUG] args.ph = 7.4
[DEBUG] sim_params = {'temperature_k': 298, 'ph': 7.4, ...}
```

**If buggy:**
```
[DEBUG] args.temperature_k = 26  ← BUG! Should be 298
[DEBUG] args.ph = 26  ← BUG! Temperature confused with pH
```

## Testing

### Test Case 1: Identical Shifts (Correction Should Apply)
```bash
# Molecule with geminal CH₃ on C=C
smiles = "CC(C)=CCCC"  # Prenyl-like structure

# Expected ML prediction (if model doesn't discriminate)
CH₃(1) = 1.70 ppm
CH₃(2) = 1.70 ppm

# Expected after correction
CH₃(cis) = 1.67 ppm
CH₃(trans) = 1.73 ppm

# Log should show:
"Geminal CH3 correction: 1.700 ppm → 1.670 / 1.730 ppm (identical shifts)"
```

### Test Case 2: Different Shifts (Correction Should Skip)
```bash
# Molecule: Flavoglaucin (CCCCCCCC1=C(O)C=C(CC=C(C)C)C(O)=C1C=O)

# Expected ML prediction (if model already discriminates)
CH₃(cis) = 1.64 ppm
CH₃(trans) = 1.70 ppm

# Expected: NO correction applied
CH₃(cis) = 1.64 ppm (unchanged)
CH₃(trans) = 1.70 ppm (unchanged)

# Log should show:
"Geminal CH3 skip: 1.640 vs 1.700 ppm (already different, Δ = 0.060 ppm)"
```

## Files Modified

1. `backend/scripts/integrated_nmr_prediction.py` (lines ~1526-1570)
   - Added check for identical shifts before applying correction
   - Added pair deduplication to avoid double-processing
   - Added detailed logging for correction decisions
   - Added debug logging for pH/temperature parameters

## Future Improvements

1. **Adaptive Threshold**: Could use model uncertainty estimates
   - If prediction uncertainty > 0.05 ppm → apply correction
   - If prediction uncertainty < 0.05 ppm → trust ML prediction

2. **ML Training Enhancement**: Include diastereotopic CH₃ pairs in training data
   - Current models may not learn this subtle effect
   - Augment training set with terpenoids, prenylated compounds

3. **Stereochemistry-Aware Correction**: Use 3D conformer to determine cis/trans
   - Current method assigns based on index order (arbitrary)
   - Better: Use spatial orientation relative to substituent
