// src/constants/colors.js
// Palette de couleurs et constantes visuelles pour l'application NMR

export const COLOR_OPTIONS = [
  // === BLEUS ===
  { name: "Bleu Acier", hex: "#31688e", rgb: [0.19, 0.41, 0.56] },
  { name: "Bleu Ardoise", hex: "#355C7D", rgb: [0.21, 0.36, 0.49] },
  { name: "Bleu Azur", hex: "#1B98E0", rgb: [0.11, 0.6, 0.88] },
  { name: "Bleu Ciel", hex: "#00B4D8", rgb: [0, 0.71, 0.85] },
  { name: "Bleu Cobalt", hex: "#247BA0", rgb: [0.14, 0.48, 0.63] },
  { name: "Bleu Glacé", hex: "#CAF0F8", rgb: [0.79, 0.94, 0.97] },
  { name: "Bleu Marine", hex: "#006994", rgb: [0, 0.41, 0.58] },
  { name: "Bleu Minuit", hex: "#13293D", rgb: [0.07, 0.16, 0.24] },
  { name: "Bleu Nuit", hex: "#0d0887", rgb: [0.05, 0.03, 0.53] },
  { name: "Bleu Océan", hex: "#0088CC", rgb: [0, 0.53, 0.8] },
  { name: "Bleu Pâle", hex: "#90E0EF", rgb: [0.56, 0.88, 0.94] },
  { name: "Bleu Pur", hex: "#0000FF", rgb: [0, 0, 1] },
  { name: "Bleu Saphir", hex: "#0077B6", rgb: [0, 0.47, 0.71] },
  { name: "Bleu Violet", hex: "#6666FF", rgb: [0.4, 0.4, 1] },

  // === CYANS / TURQUOISES ===
  { name: "Cyan Brillant", hex: "#00FFFF", rgb: [0, 1, 1] },
  { name: "Cyan Clair", hex: "#0099CC", rgb: [0, 0.6, 0.8] },
  { name: "Cyan Électrique", hex: "#00F5FF", rgb: [0, 0.96, 1] },
  { name: "Turquoise", hex: "#00B8A9", rgb: [0, 0.72, 0.66] },

  // === GRIS / NOIR ===
  { name: "Gris", hex: "#808080", rgb: [0.5, 0.5, 0.5] },
  { name: "Gris Clair", hex: "#C0C0C0", rgb: [0.75, 0.75, 0.75] },
  { name: "Gris Foncé", hex: "#404040", rgb: [0.25, 0.25, 0.25] },
  { name: "Noir", hex: "#000000", rgb: [0, 0, 0] },

  // === JAUNES / OR ===
  { name: "Citron Vert", hex: "#fde724", rgb: [0.99, 0.91, 0.14] },
  { name: "Jaune Doré", hex: "#FFCC00", rgb: [1, 0.8, 0] },
  { name: "Jaune Fluo", hex: "#FFBE0B", rgb: [1, 0.75, 0.04] },
  { name: "Jaune Pastel", hex: "#FFDD44", rgb: [1, 0.87, 0.27] },
  { name: "Jaune Pur", hex: "#FFFF00", rgb: [1, 1, 0] },
  { name: "Or", hex: "#FDC830", rgb: [0.99, 0.78, 0.19] },

  // === MARRONS / BRUNS / BEIGE ===
  { name: "Abricot", hex: "#F4A259", rgb: [0.96, 0.64, 0.35] },
  { name: "Beige", hex: "#DEB887", rgb: [0.87, 0.72, 0.53] },
  { name: "Brun Terre", hex: "#996600", rgb: [0.6, 0.4, 0] },
  { name: "Chocolat", hex: "#D2691E", rgb: [0.82, 0.41, 0.12] },
  { name: "Marron", hex: "#8B4513", rgb: [0.55, 0.27, 0.07] },
  { name: "Saumon", hex: "#F4A460", rgb: [0.96, 0.64, 0.38] },
  { name: "Terracotta", hex: "#E07A5F", rgb: [0.88, 0.48, 0.37] },

  // === ORANGES ===
  { name: "Ambre", hex: "#F7931E", rgb: [0.97, 0.58, 0.12] },
  { name: "Orange Feu", hex: "#FF6B35", rgb: [1, 0.42, 0.21] },
  { name: "Orange Foncé", hex: "#CC6600", rgb: [0.8, 0.4, 0] },
  { name: "Orange Pêche", hex: "#f89540", rgb: [0.97, 0.58, 0.25] },
  { name: "Orange Vif", hex: "#FFB000", rgb: [1, 0.69, 0] },
  { name: "Pêche Clair", hex: "#FFD580", rgb: [1, 0.84, 0.5] },

  // === ROSES / MAGENTA ===
  { name: "Magenta Pourpre", hex: "#CC0099", rgb: [0.8, 0, 0.6] },
  { name: "Magenta Pur", hex: "#FF00FF", rgb: [1, 0, 1] },
  { name: "Mauve", hex: "#C06C84", rgb: [0.75, 0.42, 0.52] },
  { name: "Rose Bonbon", hex: "#FF6B9D", rgb: [1, 0.42, 0.62] },
  { name: "Rose Bubble", hex: "#FF69B4", rgb: [1, 0.41, 0.71] },
  { name: "Rose Corail", hex: "#cc4778", rgb: [0.8, 0.28, 0.47] },
  { name: "Rose Fluo", hex: "#FF006E", rgb: [1, 0, 0.43] },
  { name: "Rose Foncé", hex: "#CC3366", rgb: [0.8, 0.2, 0.4] },
  { name: "Rose Magenta", hex: "#F72585", rgb: [0.97, 0.15, 0.52] },
  { name: "Rose Néon", hex: "#FF10F0", rgb: [1, 0.06, 0.94] },
  { name: "Rose Poudré", hex: "#F67280", rgb: [0.96, 0.45, 0.5] },
  { name: "Rose Tendre", hex: "#FF9ECD", rgb: [1, 0.62, 0.8] },

  // === ROUGES ===
  { name: "Bordeaux", hex: "#8E3B46", rgb: [0.56, 0.23, 0.27] },
  { name: "Rouge Brique", hex: "#C84B31", rgb: [0.78, 0.29, 0.19] },
  { name: "Rouge Corail", hex: "#FF6666", rgb: [1, 0.4, 0.4] },
  { name: "Rouge Pur", hex: "#FF0000", rgb: [1, 0, 0] },
  { name: "Rouge Vermillon", hex: "#E63946", rgb: [0.9, 0.22, 0.27] },

  // === VERTS ===
  { name: "Menthe", hex: "#88D8C0", rgb: [0.53, 0.85, 0.75] },
  { name: "Vert Foncé", hex: "#2D5016", rgb: [0.18, 0.31, 0.09] },
  { name: "Vert Forêt", hex: "#4A7C59", rgb: [0.29, 0.49, 0.35] },
  { name: "Vert Lime", hex: "#C5D86D", rgb: [0.77, 0.85, 0.43] },
  { name: "Vert Mousse", hex: "#87A96B", rgb: [0.53, 0.66, 0.42] },
  { name: "Vert Néon", hex: "#39FF14", rgb: [0.22, 1, 0.08] },
  { name: "Vert Pomme", hex: "#35b779", rgb: [0.21, 0.72, 0.47] },
  { name: "Vert Pur", hex: "#00FF00", rgb: [0, 1, 0] },
  { name: "Vert Turquoise", hex: "#00CC99", rgb: [0, 0.8, 0.6] },

  // === VIOLETS / MAUVES ===
  { name: "Indigo", hex: "#440154", rgb: [0.27, 0, 0.33] },
  { name: "Lavande", hex: "#9370DB", rgb: [0.58, 0.44, 0.86] },
  { name: "Lilas", hex: "#A890D3", rgb: [0.66, 0.56, 0.83] },
  { name: "Violet Clair", hex: "#9933CC", rgb: [0.6, 0.2, 0.8] },
  { name: "Violet Électrique", hex: "#7e03a8", rgb: [0.49, 0.01, 0.66] },
  { name: "Violet Fluo", hex: "#8338EC", rgb: [0.51, 0.22, 0.93] },
  { name: "Violet Profond", hex: "#3A0CA3", rgb: [0.23, 0.05, 0.64] },
];

// Couleurs par défaut par méthode (couleurs distinctes et bien visibles)
export const DEFAULT_COLORS = {
  H1: COLOR_OPTIONS.find(c => c.name === "Bleu Océan"),
  C13: COLOR_OPTIONS.find(c => c.name === "Violet Profond"),
  COSY: COLOR_OPTIONS.find(c => c.name === "Orange Vif"),
  HSQC: COLOR_OPTIONS.find(c => c.name === "Rose Magenta"),
};

// Labels affichés pour les méthodes (H1, C13, ...)
export const METHOD_LABELS_SHORT = {
  H1: "¹H",
  C13: "¹³C",
  COSY: "COSY",
  HSQC: "HSQC",
};

// Labels affichés pour les types de spectres (1H, 13C, ...)
export const TYPE_LABELS = {
  "1H": "¹H",
  "13C": "¹³C",
  COSY: "COSY",
  HSQC: "HSQC",
};

// Constantes de style
export const COLOR_THEORY = "#2563eb";
export const COLOR_EXP = "#ef4444";
export const BG = "#ffffff";
export const AXIS = "#6b7280";

// Dimensions
export const PLOT_HEIGHT = 420;
export const MOL_HEIGHT = 300;
export const SIDEBAR_W = 200;
