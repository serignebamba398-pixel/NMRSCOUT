# Spectrum Modal: Removal of Editable Experimental Parameters

## Change Summary

Les paramètres expérimentaux (Frequency, Solvent, Temperature, pH) ont été **retirés des contrôles éditables** du modal de détail du spectre. Ces paramètres sont maintenant **définis uniquement lors de la prédiction** et affichés en lecture seule dans le modal.

## Rationale

### Avant
- L'utilisateur pouvait changer la fréquence, le solvant, la température et le pH **dans le modal**
- Confusion : "Est-ce que je modifie la prédiction ou juste la visualisation ?"
- Risque d'incohérence entre les paramètres de prédiction et ceux du modal
- Les données de multiplicité JSON étaient générées avec certains paramètres, puis re-simulées avec d'autres

### Après
- Les paramètres sont **définis une seule fois** lors de la prédiction (interface de prédiction)
- Le modal **lit et affiche** ces paramètres depuis le JSON de multiplicité
- Cohérence garantie : le spectre visualisé correspond exactement à la prédiction
- Interface simplifiée : moins de contrôles, moins de confusion

## Changes Made

### File: `frontend/src/components/viewer/SpectrumDetailModal.js`

#### 1. Removed UI Controls

**Removed** (lines ~908-1001):
```javascript
// ❌ REMOVED: Frequency selector
<div style={MODAL_STYLES.controlGroup}>
  <span style={MODAL_STYLES.label}>Frequency:</span>
  <select value={freqMhz} onChange={(e) => setFreqMhz(...)}>
    <option value={400}>400 MHz</option>
    ...
  </select>
</div>

// ❌ REMOVED: Solvent selector
<div style={MODAL_STYLES.controlGroup}>
  <span style={MODAL_STYLES.label}>Solvent:</span>
  <select value={solvent} onChange={...}>
    <option value="CDCl3">CDCl₃</option>
    ...
  </select>
</div>

// ❌ REMOVED: Temperature input
<div style={MODAL_STYLES.controlGroup}>
  <span style={MODAL_STYLES.label}>Temp:</span>
  <input type="number" value={temperature} onChange={...} />
</div>

// ❌ REMOVED: pH input
<div style={MODAL_STYLES.controlGroup}>
  <span style={MODAL_STYLES.label}>pH:</span>
  <input type="number" value={ph} onChange={...} />
</div>
```

#### 2. Added Parameter Initialization from JSON

**Added** (after line ~271):
```javascript
// ✅ Initialize parameters from multiplicity JSON
useEffect(() => {
  const multData = multiplicityJson || fetchedMultiplicityData;
  if (multData && !userChangedParams) {
    // Load experimental parameters from multiplicity JSON
    if (multData.freq_MHz) setFreqMhzInternal(multData.freq_MHz);
    if (multData.solvent) setSolventInternal(multData.solvent);
    if (multData.temperature_K) setTemperatureInternal(multData.temperature_K);
    console.log('[SpectrumDetailModal] Loaded parameters from multiplicity JSON:', {
      freq: multData.freq_MHz,
      solvent: multData.solvent,
      temp: multData.temperature_K
    });
  }
}, [multiplicityJson, fetchedMultiplicityData, userChangedParams]);
```

#### 3. Added Read-Only Display

**Added** (after Line Width control):
```javascript
// ✅ Display experimental conditions (read-only)
<div style={{ display: 'flex', gap: 12, fontSize: 11, color: '#6b7280',
              alignItems: 'center', marginLeft: 12, paddingLeft: 12,
              borderLeft: '1px solid #e5e7eb' }}>
  <span><strong>Conditions:</strong> {freqMhz} MHz, {solvent}, {temperature}K</span>
</div>
```

## User Experience Changes

### Before
<img src="before-spectrum-modal.png" alt="Before: Many editable parameters" />

**Controls visible:**
- Frequency (dropdown): 400, 500, 600, 700, 800, 900 MHz
- Line Shape: Lorentzian, Gaussian, Voigt
- Line Width: slider + input
- Solvent (dropdown): CDCl₃, D₂O, DMSO-d₆, etc.
- Temp (input): 200-400 K
- pH (input): 1-14
- Refresh button
- Vertical Scale

**Problem:** Too many controls, unclear which affect prediction vs visualization

### After
<img src="after-spectrum-modal.png" alt="After: Simplified controls" />

**Controls visible:**
- Line Shape: Lorentzian, Gaussian, Voigt *(visualization only)*
- Line Width: slider + input *(visualization only)*
- **Conditions:** 700 MHz, CDCl3, 298K *(read-only, from prediction)*
- Refresh button
- Vertical Scale

**Benefit:** Clear separation between prediction parameters (fixed) and visualization parameters (adjustable)

## Data Flow

```
┌─────────────────────────────────────────────────────────────┐
│ Prediction Interface (TabbedPredictionInterface.jsx)       │
│                                                              │
│ User selects:                                                │
│  • Frequency: 400 MHz                                        │
│  • Solvent: D2O                                              │
│  • Temperature: 298 K                                        │
│  • pH: 7.4                                                   │
│                                                              │
│ Launches prediction ──────────────────────────────────────┐ │
└─────────────────────────────────────────────────────────────┘ │
                                                                │
                                                                ▼
┌─────────────────────────────────────────────────────────────────┐
│ Backend (integrated_nmr_prediction.py)                          │
│                                                                  │
│ Generates multiplicity JSON:                                    │
│ {                                                                │
│   "freq_MHz": 400.0,           ◄─── From user input            │
│   "solvent": "D2O",            ◄─── From user input            │
│   "temperature_K": 298,        ◄─── From user input            │
│   "groups": [...],                                              │
│   ...                                                            │
│ }                                                                │
│                                                                  │
│ Saves: Flavoglaucin_multiplicity.json                           │
└─────────────────────────────────────────────────────────────────┘
                                                                │
                                                                ▼
┌─────────────────────────────────────────────────────────────────┐
│ Spectrum Modal (SpectrumDetailModal.js)                         │
│                                                                  │
│ Loads multiplicity JSON                                         │
│ Reads parameters:                                               │
│  • freqMhz = json.freq_MHz        (400 MHz) ◄─── Read-only     │
│  • solvent = json.solvent         (D2O)     ◄─── Read-only     │
│  • temperature = json.temperature_K (298K)  ◄─── Read-only     │
│                                                                  │
│ User can only adjust:                                           │
│  • lineShape (Lorentzian/Gaussian/Voigt)    ◄─── Editable      │
│  • lineWidthHz (0.5-10 Hz)                  ◄─── Editable      │
│                                                                  │
│ Display: "Conditions: 400 MHz, D2O, 298K"  ◄─── Info display   │
└─────────────────────────────────────────────────────────────────┘
```

## Benefits

### 1. **Single Source of Truth**
- Les paramètres expérimentaux sont définis **une seule fois** lors de la prédiction
- Le JSON de multiplicité contient les paramètres exacts utilisés
- Aucune possibilité d'incohérence entre prédiction et visualisation

### 2. **Simplified UI**
- Moins de contrôles dans le modal (4 contrôles retirés)
- Interface plus claire et moins intimidante
- Focus sur les paramètres de visualisation (Line Shape, Line Width)

### 3. **User Clarity**
- **Prédiction** → Définir les conditions expérimentales (Frequency, Solvent, Temp, pH)
- **Visualisation** → Ajuster l'apparence du spectre (Line Shape, Line Width)
- Séparation claire des responsabilités

### 4. **Reproducibility**
- Le fichier `*_multiplicity.json` contient tous les paramètres utilisés
- Un autre utilisateur peut ouvrir le JSON et voir exactement les conditions de prédiction
- Traçabilité complète

## Migration Notes

### For Users
- **Before:** "I changed to 400 MHz in the modal, but the spectrum looks wrong!"
- **After:** "The spectrum shows 700 MHz because that's what I selected during prediction. To get 400 MHz, I need to re-run the prediction."

### For Developers
If you need to allow re-simulation with different parameters:
1. Add a "Re-predict" button in the modal
2. This button should navigate back to the prediction interface
3. Pre-populate the form with current molecule + new parameters
4. Launch a new prediction task

**Do NOT** re-enable inline editing of experimental parameters in the modal.

## Testing

### Test Case 1: Default Parameters
1. Launch prediction with default settings (700 MHz, CDCl3, 298K, pH 7.4)
2. Open spectrum modal
3. **Expected:** Display shows "Conditions: 700 MHz, CDCl3, 298K"

### Test Case 2: Custom Parameters
1. Launch prediction with 400 MHz, D2O, 298K, pH 7.4
2. Open spectrum modal
3. **Expected:** Display shows "Conditions: 400 MHz, D2O, 298K"
4. **Expected:** OH/NH protons are NOT visible (D2O filtering applied)

### Test Case 3: Parameter Persistence
1. Open spectrum modal (shows 700 MHz)
2. Close modal
3. Re-open modal
4. **Expected:** Still shows 700 MHz (parameters persist from JSON)

## Related Documentation

- [FREQUENCY_PARAMETER_FIX.md](./FREQUENCY_PARAMETER_FIX.md) - Fix for frequency parameter defaulting to 700 MHz
- [SOLVENT_FILTERING_FIX.md](./SOLVENT_FILTERING_FIX.md) - Fix for solvent-based proton filtering

## Files Modified

1. `frontend/src/components/viewer/SpectrumDetailModal.js`
   - Removed frequency, solvent, temperature, pH controls (UI)
   - Added useEffect to initialize parameters from JSON
   - Added read-only display of experimental conditions
