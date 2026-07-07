# Fix: CH₂ Diastereotopic Correction Applied Only When Needed

## Problem

La correction diastéréotopique pour les CH₂ (méthylènes) était appliquée **systématiquement** dès qu'un pattern diastéréotopique était détecté, même quand les shifts ML prédits étaient déjà différents.

**Exemple problématique** :
- ML prédit : CH₂(Ha) = 2.64 ppm, CH₂(Hb) = 2.70 ppm (déjà différents de 0.06 ppm)
- Correction appliquée : Ha = 2.64 - 0.05 = 2.59 ppm, Hb = 2.70 + 0.05 = 2.75 ppm
- Résultat : Δ = 0.16 ppm au lieu de 0.06 ppm (sur-correction)

## Root Cause

Le code détectait les CH₂ diastéréotopiques via plusieurs méthodes :
1. Centre chiral adjacent
2. Recherche BFS de centres chiraux
3. Prochiralité (benzylique asymétrique, alcène asymétrique, cycles)
4. Analyse de prochiralité individuelle

Puis appliquait **automatiquement** un delta (±0.03-0.15 ppm selon le contexte), sans vérifier si le modèle ML avait déjà prédit des shifts différents.

## Solution

### Logique de correction conditionnelle

**Avant** :
```python
if DIASTEREOTOPIC_TABLE_AVAILABLE:
    delta_ppm, delta_reason = get_diastereotopic_delta(
        mol=molH, h1_idx=h_a, h2_idx=h_b, heavy_idx=heavy_idx,
        spectrometer_mhz=spectrometer_mhz, solvent=solvent
    )

# ❌ Appliqué systématiquement sans vérifier les shifts existants
if delta_ppm > 0:
    diastereotopic_shift_deltas[h_a] = -delta_ppm / 2.0
    diastereotopic_shift_deltas[h_b] = +delta_ppm / 2.0
```

**Après** :
```python
if DIASTEREOTOPIC_TABLE_AVAILABLE:
    delta_ppm, delta_reason = get_diastereotopic_delta(
        mol=molH, h1_idx=h_a, h2_idx=h_b, heavy_idx=heavy_idx,
        spectrometer_mhz=spectrometer_mhz, solvent=solvent
    )

# ✅ Appliquer SEULEMENT si shifts identiques (< 0.01 ppm)
if delta_ppm > 0 and proton_shifts_ppm:
    shift_a = proton_shifts_ppm.get(h_a, 0.0)
    shift_b = proton_shifts_ppm.get(h_b, 0.0)

    if abs(shift_a - shift_b) < 0.01:
        diastereotopic_shift_deltas[h_a] = -delta_ppm / 2.0
        diastereotopic_shift_deltas[h_b] = +delta_ppm / 2.0
        print(f"Applied Δδ={delta_ppm:.3f} ppm (shifts identical)")
    else:
        print(f"Skipped Δδ correction: ML already different (Δ = {abs(shift_a-shift_b):.3f} ppm)")
```

### Détection de paires déjà différenciées

La correction n'est appliquée que si :
1. **Pattern détecté** : CH₂ diastéréotopique (chiral adjacent, prochiral, benzylique asymétrique, etc.)
2. **Shifts identiques** : `|δ₁ - δ₂| < 0.01 ppm`

Si le modèle ML a déjà prédit des shifts différents (|Δδ| ≥ 0.01 ppm), la correction diastéréotopique est **sautée**.

## Behavior After Fix

### Case 1: ML predicts identical shifts for diastereotopic CH₂
```
Context: Benzyl CH₂ adjacent to chiral center
ML Prediction:
  CH₂(Ha) = 2.65 ppm
  CH₂(Hb) = 2.65 ppm
  Δ = 0.00 ppm  ← Identical

Detection: Diastereotopic CH₂ detected (adjacent_chiral)
Check: |2.65 - 2.65| = 0.00 < 0.01 ppm ✓

✅ APPLY CORRECTION (Δδ = 0.08 ppm for adjacent chiral):
  CH₂(Ha) = 2.65 - 0.04 = 2.61 ppm
  CH₂(Hb) = 2.65 + 0.04 = 2.69 ppm
  Final Δ = 0.08 ppm

Log: "[CH2_diast] Applied Δδ=0.080 ppm to H12/H13 (reason: adjacent_chiral, shifts identical: 2.650 ppm)"
```

### Case 2: ML already predicts different shifts
```
Context: Benzyl CH₂ adjacent to aromatic substituent
ML Prediction:
  CH₂(Ha) = 2.64 ppm
  CH₂(Hb) = 2.70 ppm
  Δ = 0.06 ppm  ← Already different

Detection: Diastereotopic CH₂ detected (benzylic_asymmetric_aromatic)
Check: |2.64 - 2.70| = 0.06 > 0.01 ppm ✗

✗ SKIP CORRECTION:
  CH₂(Ha) = 2.64 ppm (unchanged)
  CH₂(Hb) = 2.70 ppm (unchanged)
  Final Δ = 0.06 ppm

Log: "[CH2_diast] Skipped Δδ correction for H12/H13 (reason: benzylic_asymmetric_aromatic): ML already predicts different shifts (2.640 vs 2.700 ppm, Δ=0.060 ppm)"
```

### Case 3: Vinyl CH₂ (terminal =CH₂)
```
Context: Terminal vinyl group (R-CH=CH₂)
ML Prediction:
  CH₂(Hcis) = 5.10 ppm
  CH₂(Htrans) = 5.10 ppm
  Δ = 0.00 ppm  ← Identical

Detection: Terminal vinyl CH₂ detected
Check: |5.10 - 5.10| = 0.00 < 0.01 ppm ✓

✅ APPLY CORRECTION (Δδ = 0.20 ppm for vinyl):
  CH₂(Hcis) = 5.10 - 0.10 = 5.00 ppm
  CH₂(Htrans) = 5.10 + 0.10 = 5.20 ppm
  Final Δ = 0.20 ppm

Log: "[CH2_vinyl] Applied Δδ=0.200 ppm to H8/H9 (shifts identical: 5.100 ppm)"
```

### Case 4: Vinyl CH₂ with ML already different
```
Context: Terminal vinyl group (R-CH=CH₂)
ML Prediction:
  CH₂(Hcis) = 4.95 ppm
  CH₂(Htrans) = 5.18 ppm
  Δ = 0.23 ppm  ← Already different (ML captured E/Z effect)

Detection: Terminal vinyl CH₂ detected
Check: |4.95 - 5.18| = 0.23 > 0.01 ppm ✗

✗ SKIP CORRECTION:
  CH₂(Hcis) = 4.95 ppm (unchanged)
  CH₂(Htrans) = 5.18 ppm (unchanged)
  Final Δ = 0.23 ppm

Log: "[CH2_vinyl] Skipped Δδ correction for H8/H9: ML already predicts different shifts (4.950 vs 5.180 ppm, Δ=0.230 ppm)"
```

## CH₂ Diastereotopic Contexts

The `get_diastereotopic_delta()` function from `diastereotopic_delta_table.py` analyzes various structural contexts:

### 1. **Adjacent Chiral Center** (Δδ = 0.04-0.15 ppm)
- Distance 1 bond: 0.08-0.15 ppm
- Distance 2 bonds: 0.04-0.08 ppm
- Distance 3+ bonds: 0.03-0.05 ppm

### 2. **Benzylic CH₂** (Δδ = 0.05-0.12 ppm)
- Asymmetric aromatic substituents
- Chiral aromatic ring
- Solvent effects (ASIS in C₆D₆)

### 3. **Vinyl CH₂** (Δδ = 0.15-0.25 ppm)
- Terminal =CH₂: 0.20 ppm
- Substituted vinyl: varies by substituent

### 4. **Cyclic CH₂** (Δδ = 0.05-0.20 ppm)
- Cyclohexane axial/equatorial: 0.15 ppm
- Cyclopentane envelope: 0.10 ppm
- Fused rings: up to 0.20 ppm

### 5. **Amino Acid β-CH₂** (Δδ = 0.05-0.30 ppm)
- Aspartate, glutamate: 0.05-0.10 ppm
- Leucine, isoleucine: 0.08-0.15 ppm

### 6. **Anisotropic Systems** (Δδ = 0.10-0.50 ppm)
- Carbonyl shielding cone
- Aromatic ring current effects
- Solvent anisotropy (C₆D₆)

## Threshold Justification

**Why 0.01 ppm?**

1. **ML Prediction Precision**: DGL-GNN models have ~0.3 ppm RMSE
   - Differences < 0.01 ppm are likely numerical noise
   - Differences > 0.01 ppm are intentional predictions

2. **Chemical Significance**: Typical Δδ diastereotopic = 0.04-0.15 ppm
   - If ML predicts > 0.01 ppm, it may have captured the effect
   - If ML predicts < 0.01 ppm, the model didn't discriminate → apply correction

3. **Avoid Over-Correction**:
   - ML: 0.06 ppm → Correction: +0.08 ppm → Total: 0.14 ppm (excessive!)
   - Better: ML: 0.06 ppm → No correction → Total: 0.06 ppm (realistic)

## Log Messages

### Before Fix
```
[CH2] Separated diastereotopic protons H12/H13 (adjacent_chiral)
[CH2] Applied Δδ = 0.080 ppm
```
No information on whether correction was needed or not.

### After Fix
```
# Case 1: Correction needed (identical shifts)
[CH2_diast] Applied Δδ=0.080 ppm to H12/H13 (reason: adjacent_chiral, shifts identical: 2.650 ppm)

# Case 2: Correction skipped (already different)
[CH2_diast] Skipped Δδ correction for H12/H13 (reason: benzylic_asymmetric_aromatic): ML already predicts different shifts (2.640 vs 2.700 ppm, Δ=0.060 ppm)

# Case 3: Vinyl correction applied
[CH2_vinyl] Applied Δδ=0.200 ppm to H8/H9 (shifts identical: 5.100 ppm)

# Case 4: Vinyl correction skipped
[CH2_vinyl] Skipped Δδ correction for H8/H9: ML already predicts different shifts (4.950 vs 5.180 ppm, Δ=0.230 ppm)
```

Clear indication of decision rationale.

## Code Locations

### File: `backend/scripts/sim_mult_rulebased.py`

#### 1. Vinyl CH₂ Correction (lines ~14068-14081)
```python
# ✅ Appliquer le delta SEULEMENT si les shifts ML sont identiques
if delta_ppm > 0 and proton_shifts_ppm:
    shift_a = proton_shifts_ppm.get(h_a, 0.0)
    shift_b = proton_shifts_ppm.get(h_b, 0.0)

    if abs(shift_a - shift_b) < 0.01:
        diastereotopic_shift_deltas[h_a] = -delta_ppm / 2.0
        diastereotopic_shift_deltas[h_b] = +delta_ppm / 2.0
        print(f"[CH2_vinyl] Applied Δδ={delta_ppm:.3f} ppm (identical)")
    else:
        print(f"[CH2_vinyl] Skipped (Δ={abs(shift_a-shift_b):.3f} ppm)")
```

#### 2. General Diastereotopic CH₂ Correction (lines ~14199-14214)
```python
# ✅ Appliquer le delta SEULEMENT si les shifts ML sont identiques
if delta_ppm > 0 and proton_shifts_ppm:
    shift_a = proton_shifts_ppm.get(h_a, 0.0)
    shift_b = proton_shifts_ppm.get(h_b, 0.0)

    if abs(shift_a - shift_b) < 0.01:
        diastereotopic_shift_deltas[h_a] = -delta_ppm / 2.0
        diastereotopic_shift_deltas[h_b] = +delta_ppm / 2.0
        print(f"[CH2_diast] Applied Δδ={delta_ppm:.3f} ppm (reason: {diastereotopic_reason})")
    else:
        print(f"[CH2_diast] Skipped (reason: {diastereotopic_reason}, Δ={abs(shift_a-shift_b):.3f} ppm)")
```

## Testing

### Test Case 1: Adjacent Chiral Center (Identical Shifts)
```python
# Molecule: (2S)-2-methylbutanoic acid
# CH₂ adjacent to chiral center

smiles = "CC[C@@H](C)C(O)=O"

# Expected ML prediction (if model doesn't discriminate)
CH₂(Ha) = 1.45 ppm
CH₂(Hb) = 1.45 ppm

# Expected after correction (adjacent_chiral, distance 1)
CH₂(Ha) = 1.45 - 0.04 = 1.41 ppm (pro-R)
CH₂(Hb) = 1.45 + 0.04 = 1.49 ppm (pro-S)

# Log should show:
"[CH2_diast] Applied Δδ=0.080 ppm to H5/H6 (reason: adjacent_chiral, shifts identical: 1.450 ppm)"
```

### Test Case 2: Benzylic CH₂ (Different Shifts)
```python
# Molecule: Flavoglaucin (CCCCCCCC1=C(O)C=C(CC=C(C)C)C(O)=C1C=O)
# Benzyl CH₂ adjacent to asymmetric aromatic

# Expected ML prediction (if model already discriminates)
CH₂(Ha) = 3.28 ppm
CH₂(Hb) = 3.34 ppm

# Expected: NO correction applied
CH₂(Ha) = 3.28 ppm (unchanged)
CH₂(Hb) = 3.34 ppm (unchanged)

# Log should show:
"[CH2_diast] Skipped Δδ correction for H12/H13 (reason: benzylic_asymmetric_aromatic): ML already predicts different shifts (3.280 vs 3.340 ppm, Δ=0.060 ppm)"
```

### Test Case 3: Terminal Vinyl CH₂ (Identical Shifts)
```python
# Molecule: 1-hexene
smiles = "CCCCC=C"

# Expected ML prediction (if model doesn't discriminate cis/trans)
CH₂(Hcis) = 5.05 ppm
CH₂(Htrans) = 5.05 ppm

# Expected after correction
CH₂(Hcis) = 5.05 - 0.10 = 4.95 ppm
CH₂(Htrans) = 5.05 + 0.10 = 5.15 ppm

# Log should show:
"[CH2_vinyl] Applied Δδ=0.200 ppm to H11/H12 (shifts identical: 5.050 ppm)"
```

## Related Documentation

- [GEMINAL_CH3_CORRECTION_FIX.md](./GEMINAL_CH3_CORRECTION_FIX.md) - Same fix for geminal CH₃ on double bonds
- [METHODOLOGY_NMR_MULTIPLICITY.md](./METHODOLOGY_NMR_MULTIPLICITY.md) - Complete NMR multiplicity methodology
- `backend/scripts/diastereotopic_delta_table.py` - Table of diastereotopic deltas

## Future Improvements

1. **Adaptive Threshold**: Use model uncertainty estimates
   - If prediction uncertainty > 0.05 ppm → apply correction
   - If prediction uncertainty < 0.05 ppm → trust ML prediction

2. **ML Training Enhancement**: Include diastereotopic CH₂ pairs in training data
   - Current models may not learn this subtle effect
   - Augment training set with labeled diastereotopic examples

3. **3D Conformer Analysis**: Use spatial orientation to determine pro-R/pro-S
   - Current method assigns based on index order (arbitrary)
   - Better: Use Cahn-Ingold-Prelog rules for priority assignment

4. **Frequency-Dependent Thresholds**: Adjust Δδ based on spectrometer frequency
   - 400 MHz: minimum resolvable Δδ ~0.04 ppm
   - 900 MHz: minimum resolvable Δδ ~0.02 ppm
