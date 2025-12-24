// src/utils/colorUtils.js
// Utilitaires de conversion et manipulation des couleurs

/**
 * Clamp une valeur entre min et max
 */
export const clamp = (v, min = 0, max = 255) =>
  Math.max(min, Math.min(max, v | 0));

/**
 * Convertit une couleur hexadécimale en RGB (0-255)
 */
export const hexToRgb255 = (hex) => {
  const h = (hex || "").replace("#", "").trim();
  if (h.length === 3) {
    const r = parseInt(h[0] + h[0], 16);
    const g = parseInt(h[1] + h[1], 16);
    const b = parseInt(h[2] + h[2], 16);
    return { r, g, b };
  }
  if (h.length === 6) {
    const r = parseInt(h.slice(0, 2), 16);
    const g = parseInt(h.slice(2, 4), 16);
    const b = parseInt(h.slice(4, 6), 16);
    return { r, g, b };
  }
  return { r: 89, g: 233, b: 69 };
};

/**
 * Convertit RGB (0-255) en hexadécimal
 */
export const rgb255ToHex = ({ r, g, b }) => {
  const to2 = (n) => clamp(n).toString(16).padStart(2, "0");
  return `#${to2(r)}${to2(g)}${to2(b)}`.toUpperCase();
};

/**
 * Convertit RGB (0-255) en RGB float (0-1)
 */
export const rgb255ToFloat = ({ r, g, b }) => [r / 255, g / 255, b / 255];

/**
 * Convertit HSL en Hex
 */
export const hslToHex = (h, s, l) => {
  l /= 100;
  const a = (s * Math.min(l, 1 - l)) / 100;
  const f = (n) => {
    const k = (n + h / 30) % 12;
    const color = l - a * Math.max(Math.min(k - 3, 9 - k, 1), -1);
    return Math.round(255 * color)
      .toString(16)
      .padStart(2, "0");
  };
  return `#${f(0)}${f(8)}${f(4)}`;
};

/**
 * Convertit Hex en HSL
 */
export const hexToHsl = (hex) => {
  const rgb = hexToRgb255(hex);
  const r = rgb.r / 255;
  const g = rgb.g / 255;
  const b = rgb.b / 255;

  const max = Math.max(r, g, b);
  const min = Math.min(r, g, b);
  let h,
    s,
    l = (max + min) / 2;

  if (max === min) {
    h = s = 0;
  } else {
    const d = max - min;
    s = l > 0.5 ? d / (2 - max - min) : d / (max + min);

    switch (max) {
      case r:
        h = ((g - b) / d + (g < b ? 6 : 0)) / 6;
        break;
      case g:
        h = ((b - r) / d + 2) / 6;
        break;
      case b:
        h = ((r - g) / d + 4) / 6;
        break;
      default:
        h = 0;
    }
  }

  return {
    h: Math.round(h * 360),
    s: Math.round(s * 100),
    l: Math.round(l * 100),
  };
};
