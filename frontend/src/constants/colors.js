// src/constants/colors.js
// NMR display colors and layout constants.

export const COLOR_OPTIONS = [
  // === BLUES ===
  { name: "Steel Blue", hex: "#31688e", rgb: [0.19, 0.41, 0.56] },
  { name: "Slate Blue", hex: "#355C7D", rgb: [0.21, 0.36, 0.49] },
  { name: "Azure Blue", hex: "#1B98E0", rgb: [0.11, 0.6, 0.88] },
  { name: "Sky Blue", hex: "#00B4D8", rgb: [0, 0.71, 0.85] },
  { name: "Cobalt Blue", hex: "#247BA0", rgb: [0.14, 0.48, 0.63] },
  { name: "Ice Blue", hex: "#CAF0F8", rgb: [0.79, 0.94, 0.97] },
  { name: "Navy Blue", hex: "#006994", rgb: [0, 0.41, 0.58] },
  { name: "Midnight Blue", hex: "#13293D", rgb: [0.07, 0.16, 0.24] },
  { name: "Night Blue", hex: "#0d0887", rgb: [0.05, 0.03, 0.53] },
  { name: "Ocean Blue", hex: "#0088CC", rgb: [0, 0.53, 0.8] },
  { name: "Pale Blue", hex: "#90E0EF", rgb: [0.56, 0.88, 0.94] },
  { name: "Pure Blue", hex: "#0000FF", rgb: [0, 0, 1] },
  { name: "Sapphire Blue", hex: "#0077B6", rgb: [0, 0.47, 0.71] },
  { name: "Blue Violet", hex: "#6666FF", rgb: [0.4, 0.4, 1] },

  // === CYANS / TURQUOISES ===
  { name: "Bright Cyan", hex: "#00FFFF", rgb: [0, 1, 1] },
  { name: "Light Cyan", hex: "#0099CC", rgb: [0, 0.6, 0.8] },
  { name: "Electric Cyan", hex: "#00F5FF", rgb: [0, 0.96, 1] },
  { name: "Turquoise", hex: "#00B8A9", rgb: [0, 0.72, 0.66] },
  { name: "Emerald Green", hex: "#00A878", rgb: [0, 0.66, 0.47] },

  // === GRAYS / BLACK ===
  { name: "Gray", hex: "#808080", rgb: [0.5, 0.5, 0.5] },
  { name: "Light Gray", hex: "#C0C0C0", rgb: [0.75, 0.75, 0.75] },
  { name: "Dark Gray", hex: "#404040", rgb: [0.25, 0.25, 0.25] },
  { name: "Black", hex: "#000000", rgb: [0, 0, 0] },

  // === YELLOWS / GOLD ===
  { name: "Lime", hex: "#fde724", rgb: [0.99, 0.91, 0.14] },
  { name: "Golden Yellow", hex: "#FFCC00", rgb: [1, 0.8, 0] },
  { name: "Neon Yellow", hex: "#FFBE0B", rgb: [1, 0.75, 0.04] },
  { name: "Pastel Yellow", hex: "#FFDD44", rgb: [1, 0.87, 0.27] },
  { name: "Pure Yellow", hex: "#FFFF00", rgb: [1, 1, 0] },
  { name: "Gold", hex: "#FDC830", rgb: [0.99, 0.78, 0.19] },

  // === BROWNS / BEIGE ===
  { name: "Apricot", hex: "#F4A259", rgb: [0.96, 0.64, 0.35] },
  { name: "Beige", hex: "#DEB887", rgb: [0.87, 0.72, 0.53] },
  { name: "Earth Brown", hex: "#996600", rgb: [0.6, 0.4, 0] },
  { name: "Chocolate", hex: "#D2691E", rgb: [0.82, 0.41, 0.12] },
  { name: "Brown", hex: "#8B4513", rgb: [0.55, 0.27, 0.07] },
  { name: "Salmon", hex: "#F4A460", rgb: [0.96, 0.64, 0.38] },
  { name: "Terracotta", hex: "#E07A5F", rgb: [0.88, 0.48, 0.37] },

  // === ORANGES ===
  { name: "Amber", hex: "#F7931E", rgb: [0.97, 0.58, 0.12] },
  { name: "Flame Orange", hex: "#FF6B35", rgb: [1, 0.42, 0.21] },
  { name: "Dark Orange", hex: "#CC6600", rgb: [0.8, 0.4, 0] },
  { name: "Peach Orange", hex: "#f89540", rgb: [0.97, 0.58, 0.25] },
  { name: "Bright Orange", hex: "#FFB000", rgb: [1, 0.69, 0] },
  { name: "Light Peach", hex: "#FFD580", rgb: [1, 0.84, 0.5] },

  // === ROSES / MAGENTA ===
  { name: "Purple Magenta", hex: "#CC0099", rgb: [0.8, 0, 0.6] },
  { name: "Pure Magenta", hex: "#FF00FF", rgb: [1, 0, 1] },
  { name: "Mauve", hex: "#C06C84", rgb: [0.75, 0.42, 0.52] },
  { name: "Bubblegum Pink", hex: "#FF6B9D", rgb: [1, 0.42, 0.62] },
  { name: "Bubble Pink", hex: "#FF69B4", rgb: [1, 0.41, 0.71] },
  { name: "Coral Pink", hex: "#cc4778", rgb: [0.8, 0.28, 0.47] },
  { name: "Neon Pink", hex: "#FF006E", rgb: [1, 0, 0.43] },
  { name: "Dark Pink", hex: "#CC3366", rgb: [0.8, 0.2, 0.4] },
  { name: "Magenta Pink", hex: "#F72585", rgb: [0.97, 0.15, 0.52] },
  { name: "Neon Rose", hex: "#FF10F0", rgb: [1, 0.06, 0.94] },
  { name: "Powder Pink", hex: "#F67280", rgb: [0.96, 0.45, 0.5] },
  { name: "Soft Pink", hex: "#FF9ECD", rgb: [1, 0.62, 0.8] },

  // === REDS ===
  { name: "Burgundy", hex: "#8E3B46", rgb: [0.56, 0.23, 0.27] },
  { name: "Brick Red", hex: "#C84B31", rgb: [0.78, 0.29, 0.19] },
  { name: "Coral Red", hex: "#FF6666", rgb: [1, 0.4, 0.4] },
  { name: "Pure Red", hex: "#FF0000", rgb: [1, 0, 0] },
  { name: "Vermilion Red", hex: "#E63946", rgb: [0.9, 0.22, 0.27] },

  // === GREENS ===
  { name: "Mint", hex: "#88D8C0", rgb: [0.53, 0.85, 0.75] },
  { name: "Dark Green", hex: "#2D5016", rgb: [0.18, 0.31, 0.09] },
  { name: "Forest Green", hex: "#4A7C59", rgb: [0.29, 0.49, 0.35] },
  { name: "Lime Green", hex: "#C5D86D", rgb: [0.77, 0.85, 0.43] },
  { name: "Moss Green", hex: "#87A96B", rgb: [0.53, 0.66, 0.42] },
  { name: "Neon Green", hex: "#39FF14", rgb: [0.22, 1, 0.08] },
  { name: "Apple Green", hex: "#35b779", rgb: [0.21, 0.72, 0.47] },
  { name: "Pure Green", hex: "#00FF00", rgb: [0, 1, 0] },
  { name: "Turquoise Green", hex: "#00CC99", rgb: [0, 0.8, 0.6] },

  // === VIOLETS / MAUVES ===
  { name: "Indigo", hex: "#440154", rgb: [0.27, 0, 0.33] },
  { name: "Lavender", hex: "#9370DB", rgb: [0.58, 0.44, 0.86] },
  { name: "Lilac", hex: "#A890D3", rgb: [0.66, 0.56, 0.83] },
  { name: "Light Violet", hex: "#9933CC", rgb: [0.6, 0.2, 0.8] },
  { name: "Electric Violet", hex: "#7e03a8", rgb: [0.49, 0.01, 0.66] },
  { name: "Neon Violet", hex: "#8338EC", rgb: [0.51, 0.22, 0.93] },
  { name: "Deep Violet", hex: "#3A0CA3", rgb: [0.23, 0.05, 0.64] },
];

export const normalizeColorLabel = (color) => {
  if (!color) return color;
  const paletteColor = COLOR_OPTIONS.find(
    (option) => option.hex.toLowerCase() === color.hex?.toLowerCase()
  );
  if (paletteColor) return { ...color, name: paletteColor.name };

  const legacyNames = { "Vert Émeraude": "Emerald Green" };
  return { ...color, name: legacyNames[color.name] || color.name };
};

// Default colors for each method.
export const DEFAULT_COLORS = {
  H1: COLOR_OPTIONS.find(c => c.name === "Ocean Blue"),
  C13: COLOR_OPTIONS.find(c => c.name === "Deep Violet"),
  COSY: COLOR_OPTIONS.find(c => c.name === "Bright Orange"),
  TOCSY: COLOR_OPTIONS.find(c => c.name === "Apple Green"),
  HSQC: COLOR_OPTIONS.find(c => c.name === "Magenta Pink"),
};

// Display labels for method identifiers (H1, C13, ...).
export const METHOD_LABELS_SHORT = {
  H1: "¹H",
  C13: "¹³C",
  COSY: "COSY",
  TOCSY: "TOCSY",
  HSQC: "HSQC",
};

// Display labels for spectrum types (1H, 13C, ...).
export const TYPE_LABELS = {
  "1H": "¹H",
  "13C": "¹³C",
  COSY: "COSY",
  TOCSY: "TOCSY",
  HSQC: "HSQC",
};

// Style constants
export const COLOR_THEORY = "#2563eb";
export const COLOR_EXP = "#ef4444";
export const BG = "#ffffff";
export const AXIS = "#6b7280";

// Dimensions
export const PLOT_HEIGHT = 420;
export const MOL_HEIGHT = 300;
export const SIDEBAR_W = 200;
