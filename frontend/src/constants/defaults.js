/**
 * Valeurs par défaut et constantes globales
 */

/**
 * Couleurs par défaut pour chaque méthode NMR
 */
export const DEFAULT_COLORS = {
  H1: { name: "Turquoise", hex: "#00B8A9", rgb: [0, 0.72, 0.66] },
  C13: { name: "Violet Profond", hex: "#3A0CA3", rgb: [0.23, 0.05, 0.64] },
  COSY: { name: "Orange Vif", hex: "#FFB000", rgb: [1, 0.69, 0] },
  HSQC: { name: "Rose Magenta", hex: "#F72585", rgb: [0.97, 0.15, 0.52] },
};

/**
 * Wrapper pour sessionStorage avec gestion d'erreurs
 */
export const sessionStorage = {
  get(key, fallback = null) {
    try {
      const raw = window.sessionStorage.getItem(key);
      return raw ? JSON.parse(raw) : fallback;
    } catch {
      return fallback;
    }
  },

  set(key, value) {
    try {
      if (value === null || value === undefined) {
        window.sessionStorage.removeItem(key);
      } else {
        window.sessionStorage.setItem(key, JSON.stringify(value));
      }
    } catch (e) {
      console.warn(`sessionStorage.set failed for key "${key}":`, e);
    }
  },

  remove(key) {
    try {
      window.sessionStorage.removeItem(key);
    } catch (e) {
      console.warn(`sessionStorage.remove failed for key "${key}":`, e);
    }
  },

  clear() {
    try {
      window.sessionStorage.clear();
    } catch (e) {
      console.warn('sessionStorage.clear failed:', e);
    }
  },
};

/**
 * Configuration des paramètres de prédiction par défaut
 */
export const DEFAULT_PREDICTION_OPTIONS = {
  nb_cores: 2,
  extractors: {
    h1: { keep_exchangeable: false },
    c13: { carbon_mode: "nonquaternary" },
    hsqc: { fold: true, fold_low: 4.8618, fold_high: 69.8395 },
  },
  cosy: { min_bonds: 2, max_bonds: 3 },
  idColumn: "id",
  smilesColumn: "smiles",
  // Experimental parameters for NMR prediction
  freq_mhz: 700,
  solvent: "CDCl3",
  temperature_k: 298,
  ph: 7.4,
};

/**
 * Plages de tolérance par défaut pour l'optimisation
 */
export const DEFAULT_TOLERANCE_RANGES = {
  "1H": {
    tol_correspondence: { min: 0.05, max: 0.5, step: 0.05 },
  },
  "13C": {
    tol_correspondence: { min: 0.5, max: 5, step: 0.5 },
  },
  COSY: {
    tol_correspondence_f1: { min: 0.05, max: 0.5, step: 0.1 },
    tol_correspondence_f2: { min: 0.05, max: 0.5, step: 0.1 },
  },
  HSQC: {
    tol_correspondence_f1: { min: 0.5, max: 5, step: 1 },
    tol_correspondence_f2: { min: 0.05, max: 0.5, step: 0.1 },
  },
};

/**
 * Options de background par défaut (fenêtres PPM)
 */
export const DEFAULT_BG_OPTIONS = {
  "1H": { ppm_min: -1.0, ppm_max: 12.0 },
  "13C": { ppm_min: 0.0, ppm_max: 220.0 },
  COSY: { f2_ppm_min: 0, f2_ppm_max: 10, f1_ppm_min: 0, f1_ppm_max: 10 },
  HSQC: { f2_ppm_min: 0, f2_ppm_max: 10, f1_ppm_min: 0, f1_ppm_max: 220 },
};