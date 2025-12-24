// src/Viewer.js — VERSION CORRIGÉE avec timeouts, gestion des couleurs globale, remount forcé, reset 2D fiable
// + Bouton “Afficher toutes les molécules” avec modale paginée triée par moyenne RANSY
// Vue combinée NMR + Molécules 2D/3D + sidebar
// — Zoom rect., relevel 2D au scroll, persistance d'état
// — Sélecteur de couleurs "style Word" avec Hex/RGB + aperçu
// — Couleurs par méthode pour atomes matchés (H1, C13, COSY, HSQC)
// — CORRECTION : timeouts de sécurité pour éviter les blocages UI
// — CORRECTION : reset 1D/2D fiabilisé via nonce + recompute backend + snapshots par type & source

import React, {
  useEffect,
  useMemo,
  useState,
  useCallback,
  useRef,
} from "react";

import {
  make1DMapper,
  interactive1H,
  interactive13C,
  interactiveCOSY,
  interactiveHSQC,
  relevelByWheelCancelable,
  zoom1DFromPixels,
  API_BASE_URL,
} from "./api";

// Activer le debug en tapant dans la console : window.DEBUG_NMR = true
if (typeof window !== 'undefined' && !window.hasOwnProperty('DEBUG_NMR')) {
  window.DEBUG_NMR = false; // false par défaut en production
}

// Helper pour logs conditionnels
const debugLog = (...args) => {
  if (window.DEBUG_NMR) {
    console.log('[NMR Debug]', ...args);
  }
};

const debugWarn = (...args) => {
  if (window.DEBUG_NMR) {
    console.warn('[NMR Debug]', ...args);
  }
};

const debugError = (...args) => {
  if (window.DEBUG_NMR) {
    console.error('[NMR Debug]', ...args);
  }
};


/* ===========================
   Helper fetch avec retry
   =========================== */

async function fetchWithRetry(url, maxRetries = 3, delayMs = 500) {
  let lastError = null;

  for (let attempt = 1; attempt <= maxRetries; attempt++) {
    try {
      const r = await fetch(url);
      if (!r.ok) throw new Error(`HTTP ${r.status}`);
      return r;
    } catch (e) {
      lastError = e;
      console.warn(
        `Tentative ${attempt}/${maxRetries} échouée:`,
        e?.message || e
      );
      if (attempt < maxRetries) {
        const d = delayMs * attempt;
        await new Promise(res => setTimeout(res, d));
      }
    }
  }
  throw lastError;
}

// Construit une URL d'image "safe". N'ajoute ?cb=... que si l'URL n'est pas signée.
// ⚠️ Ajout : on prend aussi en compte meta.nonce pour casser le cache quand le backend recalcule.
function buildImgSrc(url, meta) {
  try {
    // En Electron, utiliser API_BASE_URL comme base au lieu de window.location.href (file://)
    const baseUrl = typeof API_BASE_URL === 'string' && API_BASE_URL.startsWith('http')
      ? API_BASE_URL
      : window.location.href;
    const u = new URL(url, baseUrl);

    // Tag prioritaire pour le cache-busting
    const tag =
      meta?.etag ?? meta?.version ?? meta?.last_modified ?? meta?.nonce;
    if (!tag) return url;

    const qs = u.search || "";
    const looksSigned =
      /(token|signature|sig|expires|X-Amz-Algorithm|X-Amz-Credential|X-Amz-Signature)/i.test(qs);

    if (looksSigned) {
      // ⚠️ Ne pas modifier la query (sinon signature invalide)
      // → on force le rechargement via le fragment (#cb=...)
      const currentHash = (u.hash || "").replace(/^#/, "");
      const parts = currentHash ? currentHash.split("&").filter(Boolean) : [];
      // supprime cb existant s'il y en a un
      const others = parts.filter(p => !/^cb=/.test(p));
      others.push("cb=" + encodeURIComponent(tag));
      u.hash = "#" + others.join("&");
      return u.toString();
    }

    // URL non signée → on peut utiliser la query
    const hasQuery = u.search && u.search.length > 1;
    const params = new URLSearchParams(u.search);
    params.set("cb", String(tag));
    u.search = (hasQuery ? "?" : "?") + params.toString();
    return u.toString();
  } catch {
    return url;
  }
}

/* ===========================
   Palette de couleurs
   =========================== */

const COLOR_OPTIONS = [
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
const DEFAULT_COLORS = {
  H1:   COLOR_OPTIONS.find(c => c.name === "Bleu Océan"),      // Bleu
  C13:  COLOR_OPTIONS.find(c => c.name === "Violet Profond"),  // Violet
  COSY: COLOR_OPTIONS.find(c => c.name === "Orange Vif"),      // Orange
  HSQC: COLOR_OPTIONS.find(c => c.name === "Rose Magenta"),    // Rose
};

// Labels affichés pour les méthodes (H1, C13, ...)
const METHOD_LABELS_SHORT = {
  H1: "¹H",
  C13: "¹³C",
  COSY: "COSY",
  HSQC: "HSQC",
};

// Labels affichés pour les types de spectres (1H, 13C, ...)
const TYPE_LABELS = {
  "1H": "¹H",
  "13C": "¹³C",
  COSY: "COSY",
  HSQC: "HSQC",
};

/* ===========================
   Styling basique
   =========================== */

const COLOR_THEORY = "#2563eb";
const COLOR_EXP    = "#ef4444";
const BG           = "#ffffff";
const AXIS         = "#6b7280";

const PLOT_HEIGHT  = 420;
const MOL_HEIGHT   = 300;
const SIDEBAR_W    = 200;

/* ===========================
   Petites pastilles métriques
   =========================== */

const Pill = ({ label, value }) => (
  <div
    style={{
      background: "#f2f7ff",
      border: "1px solid #dbeafe",
      borderRadius: 999,
      padding: "2px 6px",
      fontSize: 10,
      display: "inline-flex",
      gap: 3,
      alignItems: "center",
      marginRight: 4,
      marginBottom: 3,
    }}
  >
    <span style={{ color: "#2563eb", fontWeight: 600 }}>
      {label}:
    </span>
    <span style={{ color: "#0f172a", fontWeight: 700 }}>
      {String(value)}
    </span>
  </div>
);

/* ===========================
   Helpers couleur (Hex/RGB)
   =========================== */

const clamp = (v, min = 0, max = 255) =>
  Math.max(min, Math.min(max, v | 0));

const hexToRgb255 = (hex) => {
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

const rgb255ToHex = ({ r, g, b }) => {
  const to2 = (n) => clamp(n).toString(16).padStart(2, "0");
  return `#${to2(r)}${to2(g)}${to2(b)}`.toUpperCase();
};

const rgb255ToFloat = ({ r, g, b }) => [r / 255, g / 255, b / 255];

/* ===========================
   Sélecteur de couleurs "Word"
   =========================== */

/* ===========================
   Sélecteur de couleurs "Word" - VERSION CROSS-BROWSER
   =========================== */

function ColorPicker({
  isOpen,
  onClose,
  customColors,
  onColorChange,
}) {
  const [draft, setDraft] = useState(customColors);
  const [active, setActive] = useState("H1");
  const [showColorPalette, setShowColorPalette] = useState(false);
  const [tempColor, setTempColor] = useState({ h: 180, s: 70, l: 50 });

  useEffect(() => { setDraft(customColors); }, [customColors]);
  if (!isOpen) return null;

  const applyColor = (method, hexOrRgb) => {
    let next;
    if (typeof hexOrRgb === "string") {
      const rgb255 = hexToRgb255(hexOrRgb);
      next = {
        name: "Custom",
        hex: rgb255ToHex(rgb255),
        rgb: rgb255ToFloat(rgb255),
      };
    } else {
      const hex = rgb255ToHex(hexOrRgb);
      next = {
        name: "Custom",
        hex,
        rgb: rgb255ToFloat(hexOrRgb),
      };
    }
    setDraft((d) => ({ ...d, [method]: next }));
  };

  const applyAll = () => {
    ["H1", "C13", "COSY", "HSQC"].forEach((m) => {
      onColorChange(m, draft[m]);
    });
    onClose();
  };

  // ✅ Conversion HSL vers Hex
  const hslToHex = (h, s, l) => {
    l /= 100;
    const a = s * Math.min(l, 1 - l) / 100;
    const f = n => {
      const k = (n + h / 30) % 12;
      const color = l - a * Math.max(Math.min(k - 3, 9 - k, 1), -1);
      return Math.round(255 * color).toString(16).padStart(2, '0');
    };
    return `#${f(0)}${f(8)}${f(4)}`;
  };

  // ✅ Conversion Hex vers HSL
  const hexToHsl = (hex) => {
    const rgb = hexToRgb255(hex);
    const r = rgb.r / 255;
    const g = rgb.g / 255;
    const b = rgb.b / 255;
    
    const max = Math.max(r, g, b);
    const min = Math.min(r, g, b);
    let h, s, l = (max + min) / 2;

    if (max === min) {
      h = s = 0;
    } else {
      const d = max - min;
      s = l > 0.5 ? d / (2 - max - min) : d / (max + min);
      
      switch (max) {
        case r: h = ((g - b) / d + (g < b ? 6 : 0)) / 6; break;
        case g: h = ((b - r) / d + 2) / 6; break;
        case b: h = ((r - g) / d + 4) / 6; break;
        default: h = 0;
      }
    }

    return {
      h: Math.round(h * 360),
      s: Math.round(s * 100),
      l: Math.round(l * 100)
    };
  };

  const MethodRow = ({ method, title }) => {
    const current = customColors[method];
    const d = draft[method];
    const newRGB = hexToRgb255(d.hex);

    return (
      <div
        style={{
          border: "1px solid #e5e7eb",
          borderRadius: 10,
          padding: 14,
          display: "grid",
          gap: 10,
        }}
      >
        <div
          style={{
            display: "flex",
            alignItems: "center",
            justifyContent: "space-between",
            gap: 10,
          }}
        >
          <div style={{ display: "flex", gap: 8, alignItems: "center" }}>
            <div
              style={{
                width: 24,
                height: 24,
                borderRadius: 4,
                background: d.hex,
                border: "2px solid #e5e7eb",
              }}
            />
            <div>
              <div style={{ fontWeight: 700, fontSize: 14 }}>
                {title}
              </div>
              <div style={{ fontSize: 11, color: "#64748b" }}>
                {d.name}
              </div>
            </div>
          </div>

          <select
            value={d.name}
            onChange={(e) => {
              const chosen =
                COLOR_OPTIONS.find(
                  c => c.name === e.target.value
                ) ||
                { name: "Custom", hex: d.hex, rgb: d.rgb };
              setDraft((prev) => ({
                ...prev,
                [method]:
                  chosen.name === "Custom"
                    ? { ...prev[method], name: "Custom" }
                    : chosen,
              }));
            }}
            style={{
              padding: "6px 10px",
              fontSize: 12,
              borderRadius: 6,
              border: "1px solid #cbd5e1",
              background: "#fff",
              cursor: "pointer",
            }}
            title="Palette rapide"
          >
            <option value={d.name} disabled>
              — Palette rapide —
            </option>
            {COLOR_OPTIONS.map((c) => (
              <option key={c.name} value={c.name}>
                {c.name}
              </option>
            ))}
          </select>
        </div>

        <div
          style={{
            display: "grid",
            gridTemplateColumns: "220px 1fr",
            gap: 14,
            alignItems: "start",
          }}
        >
          {/* Aperçu des couleurs */}
          <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 10 }}>
            <div>
              <div
                style={{
                  fontSize: 11,
                  marginBottom: 6,
                  color: "#475569",
                }}
              >
                New
              </div>
              <div
                style={{
                  width: "100%",
                  height: 56,
                  borderRadius: 8,
                  border: "1px solid #cbd5e1",
                  background: d.hex,
                }}
              />
            </div>

            <div>
              <div
                style={{
                  fontSize: 11,
                  marginBottom: 6,
                  color: "#475569",
                }}
              >
                Current
              </div>
              <div
                style={{
                  width: "100%",
                  height: 56,
                  borderRadius: 8,
                  border: "1px solid #cbd5e1",
                  background: current.hex,
                }}
              />
            </div>

            {/* ✅ Bouton pour ouvrir la palette custom */}
            <div style={{ gridColumn: "1 / span 2", marginTop: 8, position: "relative" }}>
              <button
                onClick={() => {
                  const hsl = hexToHsl(d.hex);
                  setTempColor(hsl);
                  setShowColorPalette(!showColorPalette);
                }}
                style={{
                  width: "100%",
                  height: 36,
                  border: "1px solid #cbd5e1",
                  borderRadius: 6,
                  cursor: "pointer",
                  background: `linear-gradient(to right, 
                    ${d.hex} 0%, ${d.hex} 70%, 
                    #fff 70%, #fff 100%)`,
                  display: "flex",
                  alignItems: "center",
                  justifyContent: "flex-end",
                  padding: "0 8px",
                  fontSize: 11,
                  fontWeight: 600,
                  color: "#475569",
                }}
              >
                Choisir ▼
              </button>

              {/* ✅ Palette custom (s'affiche au clic) */}
              {showColorPalette && (
                <div
                  style={{
                    position: "absolute",
                    top: "calc(100% + 4px)",
                    left: 0,
                    zIndex: 10000,
                    background: "#fff",
                    border: "1px solid #cbd5e1",
                    borderRadius: 8,
                    padding: 12,
                    boxShadow: "0 10px 40px rgba(0,0,0,0.15)",
                    width: "100%",
                  }}
                  onClick={(e) => e.stopPropagation()}
                >
                  {/* Gradient 2D (Saturation x Luminosité) */}
                  <div
                    style={{
                      width: "100%",
                      height: 150,
                      borderRadius: 6,
                      background: `linear-gradient(to bottom, 
                        hsl(${tempColor.h}, 100%, 50%) 0%,
                        hsl(${tempColor.h}, 100%, 0%) 100%)`,
                      position: "relative",
                      cursor: "crosshair",
                      marginBottom: 10,
                    }}
                    onClick={(e) => {
                      const rect = e.currentTarget.getBoundingClientRect();
                      const x = e.clientX - rect.left;
                      const y = e.clientY - rect.top;
                      const s = Math.round((x / rect.width) * 100);
                      const l = Math.round(100 - (y / rect.height) * 100);
                      setTempColor({ ...tempColor, s, l });
                      const hex = hslToHex(tempColor.h, s, l);
                      applyColor(method, hex);
                    }}
                  >
                    <div
                      style={{
                        position: "absolute",
                        top: 0,
                        left: 0,
                        right: 0,
                        bottom: 0,
                        background: "linear-gradient(to right, #fff 0%, transparent 100%)",
                        borderRadius: 6,
                      }}
                    />
                  </div>

                  {/* Slider Teinte (Hue) */}
                  <div style={{ marginBottom: 8 }}>
                    <div style={{ fontSize: 10, marginBottom: 4, color: "#475569" }}>
                      Teinte: {tempColor.h}°
                    </div>
                    <input
                      type="range"
                      min="0"
                      max="360"
                      value={tempColor.h}
                      onChange={(e) => {
                        const h = parseInt(e.target.value);
                        setTempColor({ ...tempColor, h });
                        const hex = hslToHex(h, tempColor.s, tempColor.l);
                        applyColor(method, hex);
                      }}
                      style={{
                        width: "100%",
                        height: 20,
                        borderRadius: 10,
                        background: "linear-gradient(to right, #ff0000 0%, #ffff00 17%, #00ff00 33%, #00ffff 50%, #0000ff 67%, #ff00ff 83%, #ff0000 100%)",
                        outline: "none",
                        WebkitAppearance: "none",
                        appearance: "none",
                      }}
                    />
                  </div>

                  <button
                    onClick={() => setShowColorPalette(false)}
                    style={{
                      width: "100%",
                      padding: "6px",
                      borderRadius: 6,
                      border: "1px solid #cbd5e1",
                      background: "#f8fafc",
                      cursor: "pointer",
                      fontSize: 11,
                      fontWeight: 600,
                    }}
                  >
                    Fermer
                  </button>
                </div>
              )}
            </div>
          </div>

          {/* Champs Hex et RGB */}
          <div
            style={{
              display: "grid",
              gridTemplateColumns: "120px 1fr",
              rowGap: 8,
              columnGap: 10,
            }}
          >
            <label
              style={{
                fontSize: 12,
                color: "#334155",
                alignSelf: "center",
              }}
            >
              Hex :
            </label>
            <div style={{ display: "flex", gap: 8 }}>
              <input
                value={d.hex.toUpperCase()}
                onChange={(e) => {
                  const raw = e.target.value.replace("#", "");
                  const ok = /^[0-9a-fA-F]{0,6}$/.test(raw);
                  if (!ok) return;
                  const v = "#" + raw;
                  applyColor(method, v);
                  if (raw.length === 6) {
                    const hsl = hexToHsl(v);
                    setTempColor(hsl);
                  }
                }}
                maxLength={7}
                style={{
                  width: 120,
                  padding: "6px 8px",
                  fontSize: 12,
                  borderRadius: 6,
                  border: "1px solid #cbd5e1",
                }}
              />
              <button
                onClick={() => onColorChange(method, draft[method])}
                style={{
                  padding: "6px 12px",
                  borderRadius: 6,
                  border: "none",
                  background: "#2563eb",
                  color: "#fff",
                  fontSize: 12,
                  fontWeight: 600,
                  cursor: "pointer",
                }}
              >
                Appliquer
              </button>
            </div>

            <label
              style={{
                fontSize: 12,
                color: "#334155",
                alignSelf: "center",
              }}
            >
              RGB :
            </label>
            <div style={{ display: "flex", gap: 8 }}>
              <input
                type="number"
                min={0}
                max={255}
                value={newRGB.r}
                onChange={(e) => {
                  const next = {
                    r: clamp(+e.target.value),
                    g: newRGB.g,
                    b: newRGB.b,
                  };
                  applyColor(method, next);
                  const hex = rgb255ToHex(next);
                  const hsl = hexToHsl(hex);
                  setTempColor(hsl);
                }}
                style={{
                  width: 70,
                  padding: "6px 8px",
                  fontSize: 12,
                  borderRadius: 6,
                  border: "1px solid #cbd5e1",
                }}
              />
              <input
                type="number"
                min={0}
                max={255}
                value={newRGB.g}
                onChange={(e) => {
                  const next = {
                    r: newRGB.r,
                    g: clamp(+e.target.value),
                    b: newRGB.b,
                  };
                  applyColor(method, next);
                  const hex = rgb255ToHex(next);
                  const hsl = hexToHsl(hex);
                  setTempColor(hsl);
                }}
                style={{
                  width: 70,
                  padding: "6px 8px",
                  fontSize: 12,
                  borderRadius: 6,
                  border: "1px solid #cbd5e1",
                }}
              />
              <input
                type="number"
                min={0}
                max={255}
                value={newRGB.b}
                onChange={(e) => {
                  const next = {
                    r: newRGB.r,
                    g: newRGB.g,
                    b: clamp(+e.target.value),
                  };
                  applyColor(method, next);
                  const hex = rgb255ToHex(next);
                  const hsl = hexToHsl(hex);
                  setTempColor(hsl);
                }}
                style={{
                  width: 70,
                  padding: "6px 8px",
                  fontSize: 12,
                  borderRadius: 6,
                  border: "1px solid #cbd5e1",
                }}
              />
            </div>
          </div>
        </div>

        <div style={{ display: "flex", gap: 8, flexWrap: "wrap" }}>
          <button
            onClick={() => onColorChange(method, draft[method])}
            style={{
              padding: "6px 12px",
              borderRadius: 6,
              border: "none",
              background: "#2563eb",
              color: "#fff",
              fontSize: 12,
              fontWeight: 700,
              cursor: "pointer",
            }}
          >
            Appliquer {method}
          </button>
          <button
            onClick={() => onColorChange(method, DEFAULT_COLORS[method])}
            style={{
              padding: "6px 12px",
              borderRadius: 6,
              border: "1px solid #cbd5e1",
              background: "#fff",
              fontSize: 12,
              cursor: "pointer",
            }}
          >
            Réinitialiser {method}
          </button>
        </div>
      </div>
    );
  };

  return (
    <div
      style={{
        position: "fixed",
        inset: 0,
        background: "rgba(0,0,0,0.5)",
        display: "flex",
        alignItems: "center",
        justifyContent: "center",
        zIndex: 9999,
      }}
      onClick={(e) => {
        if (e.target === e.currentTarget) {
          setShowColorPalette(false);
          onClose();
        }
      }}
    >
      <div
        style={{
          background: "#fff",
          borderRadius: 12,
          padding: 18,
          width: "min(880px, 94vw)",
          maxHeight: "86vh",
          overflow: "auto",
          boxShadow: "0 20px 60px rgba(0,0,0,.3)",
        }}
        onClick={(e) => e.stopPropagation()}
      >
        <div
          style={{
            display: "flex",
            justifyContent: "space-between",
            alignItems: "center",
          }}
        >
          <div style={{ fontSize: 18, fontWeight: 700 }}>
            Personnaliser les couleurs (style Word)
          </div>
          <button
            onClick={() => {
              setShowColorPalette(false);
              onClose();
            }}
            style={{
              background: "none",
              border: "none",
              fontSize: 24,
              cursor: "pointer",
              color: "#64748b",
            }}
          >
            ×
          </button>
        </div>

        <div style={{ display: "flex", gap: 8, margin: "12px 0" }}>
          {["H1", "C13", "COSY", "HSQC"].map((m) => (
            <button
              key={m}
              onClick={() => {
                setShowColorPalette(false);
                setActive(m);
              }}
              style={{
                padding: "6px 10px",
                borderRadius: 999,
                border: "1px solid #cbd5e1",
                background: active === m ? "#2563eb" : "#fff",
                color: active === m ? "#fff" : "#0f172a",
                fontSize: 12,
                cursor: "pointer",
              }}
            >
              {METHOD_LABELS_SHORT[m] || m}
            </button>
          ))}
        </div>

        {active === "H1"   && <MethodRow method="H1"   title="¹H NMR" />}
        {active === "C13"  && <MethodRow method="C13"  title="¹³C NMR" />}
        {active === "COSY" && <MethodRow method="COSY" title="COSY" />}
        {active === "HSQC" && <MethodRow method="HSQC" title="HSQC" />}

        <div
          style={{ display: "flex", justifyContent: "flex-end",
                   gap: 8, marginTop: 14 }}
        >
          <button
            onClick={() => {
              ["H1", "C13", "COSY", "HSQC"].forEach((m) =>
                onColorChange(m, DEFAULT_COLORS[m])
              );
            }}
            style={{
              padding: "8px 14px",
              borderRadius: 6,
              border: "1px solid #cbd5e1",
              background: "#fff",
              fontSize: 13,
              cursor: "pointer",
            }}
          >
            Réinitialiser tout
          </button>
          <button
            onClick={() => {
              setShowColorPalette(false);
              applyAll();
            }}
            style={{
              padding: "8px 14px",
              borderRadius: 6,
              border: "none",
              background: "#2563eb",
              color: "#fff",
              fontWeight: 700,
              fontSize: 13,
              cursor: "pointer",
            }}
          >
            Terminé
          </button>
        </div>
      </div>
    </div>
  );
}
/* ===========================
   RDKit (2D) — chargeur CDN
   =========================== */

let rdkitPromise = null;

function loadRDKitFromCDN() {
  if (typeof window === "undefined") {
    return Promise.reject(new Error("window indisponible"));
  }
  if (window.RDKit) return Promise.resolve(window.RDKit);
  if (rdkitPromise) return rdkitPromise;

  rdkitPromise = new Promise((resolve, reject) => {
    const sel = 'script[data-rdkit="RDKit_minimal"]';
    const ex = document.querySelector(sel);

    const onReady = () => {
      if (typeof window.initRDKitModule !== "function") {
        reject(new Error("initRDKitModule introuvable"));
        return;
      }
      window
        .initRDKitModule()
        .then((RDKit) => {
          window.RDKit = RDKit;
          resolve(RDKit);
        })
        .catch(reject);
    };

    if (ex) {
      ex.addEventListener("load", onReady, { once: true });
      ex.addEventListener(
        "error",
        () => reject(new Error("Échec RDKit.js")),
        { once: true }
      );
      return;
    }

    const s = document.createElement("script");
    s.async = true;
    s.defer = true;
    s.dataset.rdkit = "RDKit_minimal";
    s.src = "https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js";
    s.onload = onReady;
    s.onerror = () => reject(new Error("Échec RDKit.js"));
    document.head.appendChild(s);
  });

  return rdkitPromise;
}

/* ===========================
   3Dmol (3D) — chargeur CDN
   =========================== */

function isProblematicBrowser() {
  const ua = navigator.userAgent.toLowerCase();
  const isFirefox = ua.includes("firefox");
  const isSafari =
    /^((?!chrome|android).)*safari/i.test(navigator.userAgent);
  return isFirefox || isSafari;
}

// ✅ NOUVELLE FONCTION : Tester WebGL avant d'essayer 3Dmol
function isWebGLAvailable() {
  try {
    const canvas = document.createElement('canvas');
    const gl = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');
    if (!gl) return false;
    
    // Test basique de rendu
    gl.clearColor(0, 0, 0, 1);
    gl.clear(gl.COLOR_BUFFER_BIT);
    
    return true;
  } catch (e) {
    console.error("Test WebGL échoué:", e);
    return false;
  }
}

function disableOffscreenCanvasCompletely() {
  try {
    if ("OffscreenCanvas" in window) {
      window._originalOffscreenCanvas = window.OffscreenCanvas;
      delete window.OffscreenCanvas;
      Object.defineProperty(window, "OffscreenCanvas", {
        value: undefined,
        writable: false,
        configurable: false,
      });
    }
    if (HTMLCanvasElement.prototype.transferControlToOffscreen) {
      HTMLCanvasElement.prototype._originalTransferControlToOffscreen =
        HTMLCanvasElement.prototype.transferControlToOffscreen;
      HTMLCanvasElement.prototype.transferControlToOffscreen =
        function () { return this; };
    }
    if (isProblematicBrowser() && window.createImageBitmap) {
      window._originalCreateImageBitmap = window.createImageBitmap;
      window.createImageBitmap = function () {
        throw new Error("createImageBitmap désactivé");
      };
    }
    window.OffscreenCanvas = function () {
      throw new Error("OffscreenCanvas désactivé");
    };
  } catch {}
}


function patch3DmolForOffscreenCanvas() {
  if (typeof window.$3Dmol === "undefined") {
    setTimeout(patch3DmolForOffscreenCanvas, 100);
    return;
  }
  try {
    if (window.$3Dmol) {
      window.$3Dmol.global_enable_offscreen = false;
      if (window.$3Dmol.GLViewer) {
        window.$3Dmol.GLViewer.prototype._enable_offscreen = false;
      }
      const orig = window.$3Dmol.createViewer;
      if (orig) {
        window.$3Dmol.createViewer = function (el, cfg = {}) {
          const safe = {
            ...cfg,
            disableOffscreenCanvas: true,
            antialias: isProblematicBrowser()
              ? false
              : cfg.antialias !== false,
            preserveDrawingBuffer: true,
          };
          return orig.call(this, el, safe);
        };
      }
    }
  } catch {}
}

let mol3dPromise = null;

function load3DmolFromCDN() {
  if (typeof window === "undefined") {
    return Promise.reject(new Error("window indisponible"));
  }
  if (window.$3Dmol) return Promise.resolve(window.$3Dmol);
  if (mol3dPromise) return mol3dPromise;

  disableOffscreenCanvasCompletely();

  const urls = [
    "https://3dmol.org/build/3Dmol-min.js",
    "https://cdn.jsdelivr.net/npm/3dmol/build/3Dmol-min.js",
  ];

  mol3dPromise = new Promise(async (resolve, reject) => {
    // ✅ AJOUT : Vérifier si WebGL fonctionne AVANT de charger 3Dmol
    if (!isWebGLAvailable()) {
      reject(new Error("WebGL non disponible sur ce navigateur. Veuillez activer WebGL dans about:config (webgl.disabled = false) ou mettre à jour vos pilotes graphiques."));
      return;
    }

    for (const u of urls) {
      try {
        await new Promise((res, rej) => {
          const ex = document.querySelector(`script[src="${u}"]`);
          if (ex) {
            ex.addEventListener("load", res, { once: true });
            ex.addEventListener("error", rej, { once: true });
            return;
          }
          const s = document.createElement("script");
          s.async = true;
          s.defer = true;
          s.src = u;
          s.onload = res;
          s.onerror = rej;
          document.head.appendChild(s);
        });
        if (window.$3Dmol) {
          patch3DmolForOffscreenCanvas();
          resolve(window.$3Dmol);
          return;
        }
      } catch {}
    }
    reject(new Error("Impossible de charger 3Dmol"));
  });

  return mol3dPromise;
}

/* ===========================
   CSV parsing + petits helpers
   =========================== */

function parseCSVSmart(text) {
  if (!text || typeof text !== "string") {
    return { headers: [], data: [] };
  }
  const first = text.split(/\r?\n/).slice(0, 5).join("\n");
  const cand = [";", ",", "\t"];
  let best = ";";
  let bestCount = -1;

  for (const d of cand) {
    const counts = first
      .split(/\r?\n/)
      .map(l => (l.match(new RegExp(`\\${d}`, "g")) || []).length);
    const cnt = Math.max(...counts, 0);
    if (cnt > bestCount) {
      best = d;
      bestCount = cnt;
    }
  }

  const rows = [];
  let cur = [];
  let val = "";
  let inQuotes = false;

  const pushVal = () => {
    if (val && val.charCodeAt(0) === 0xfeff) val = val.slice(1);
    cur.push(val);
    val = "";
  };
  const pushRow = () => { rows.push(cur); cur = []; };

  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (inQuotes) {
      if (ch === '"') {
        if (text[i + 1] === '"') { val += '"'; i++; }
        else inQuotes = false;
      } else {
        val += ch;
      }
      continue;
    }
    if (ch === '"') inQuotes = true;
    else if (ch === best) pushVal();
    else if (ch === "\n") { pushVal(); pushRow(); }
    else if (ch !== "\r") val += ch;
  }
  if (val.length > 0 || cur.length > 0) {
    pushVal();
    pushRow();
  }
  if (rows.length === 0) return { headers: [], data: [] };

  const headers = rows[0].map(h => (h || "").trim());
  const data = rows
    .slice(1)
    .filter(r => r.some(x => (x || "").trim() !== ""))
    .map(r => {
      const o = {};
      headers.forEach((h, i) => { o[h] = (r[i] ?? "").trim(); });
      return o;
    });

  return { headers, data };
}


/* ===========================
   Graph helpers (bonds)
   =========================== */

function findBondsInPath(molblock, atomPairs) {
  if (!molblock || !atomPairs.length) return [];
  const lines = molblock.split(/\r?\n/);
  const counts = lines[3] || "";
  const nAtoms = parseInt(counts.slice(0, 3).trim(), 10);
  const nBonds = parseInt(counts.slice(3, 6).trim(), 10);

  const adjacency = new Map();
  for (let i = 0; i < nAtoms; i++) adjacency.set(i, []);

  for (let i = 0; i < nBonds; i++) {
    const l = lines[4 + nAtoms + i] || "";
    const a1 = parseInt(l.slice(0, 3).trim(), 10) - 1;
    const a2 = parseInt(l.slice(3, 6).trim(), 10) - 1;
    if (a1 >= 0 && a2 >= 0 && a1 < nAtoms && a2 < nAtoms) {
      adjacency.get(a1).push({ neighbor: a2, bondIndex: i });
      adjacency.get(a2).push({ neighbor: a1, bondIndex: i });
    }
  }

  const findShortestPath = (start, end) => {
    if (start === end) return [];
    const q = [{ atom: start, path: [] }];
    const seen = new Set([start]);
    while (q.length > 0) {
      const { atom, path } = q.shift();
      for (const { neighbor, bondIndex } of adjacency.get(atom) || []) {
        if (neighbor === end) return [...path, bondIndex];
        if (!seen.has(neighbor) && path.length < 4) {
          seen.add(neighbor);
          q.push({ atom: neighbor, path: [...path, bondIndex] });
        }
      }
    }
    return [];
  };

  const hi = new Set();
  for (const [a, b] of atomPairs) {
    const p = findShortestPath(a, b);
    p.forEach(idx => hi.add(idx));
  }
  return Array.from(hi);
}



/* ===========================
   SDF / MOL helpers
   =========================== */

function extractFirstMolBlockPreserving(text) {
  if (!text) return null;
  const i = text.indexOf("M  END");
  if (i < 0) return null;
  const end = i + "M  END".length;
  let out = text.slice(0, end);
  if (!out.endsWith("\n")) out += "\n";
  return out;
}

function extractFirstSDFRecord(text) {
  if (!text) return null;
  const upto = text.indexOf("$$$$");
  const rec = upto >= 0 ? text.slice(0, upto) : text;
  const mend = rec.indexOf("M  END");
  if (mend >= 0) {
    const s = rec.slice(0, mend + "M  END".length);
    return s.endsWith("\n") ? s : s + "\n";
  }
  return rec.endsWith("\n") ? rec : rec + "\n";
}

function extractShiftArrayFromSDFRecord(sdfText) {
  if (!sdfText) return null;
  const first = sdfText.split("$$$$")[0] || sdfText;
  
  // ✅ CORRECTION : chercher <(1)> avec le format exact du backend: >  <(1)>  (1)
  // Le format est: >  <(1)>  (1)\n22.067 54.492 0.000 ...
  const re = />\s*<\(1\)>\s*(?:\(\d+\))?\s*\r?\n([\s\S]*?)(?:\r?\n\r?\n|$$$$|$)/;
  const m = first.match(re);
  if (!m) {
    // Fallback: essayer aussi <1> sans parenthèses
    const re2 = />\s*<1>\s*(?:\(\d+\))?\s*\r?\n([\s\S]*?)(?:\r?\n\r?\n|$$$$|$)/;
    const m2 = first.match(re2);
    if (!m2) return null;
    return parseShiftBlock(m2[1]);
  }

  return parseShiftBlock(m[1]);
}

function parseShiftBlock(block) {
  if (!block) return null;
  const cleaned = block
    .split(/\r?\n/)
    .map(s => s.trim())
    .filter(Boolean)
    .join(" ");

  if (!cleaned) return null;

  // ✅ Les valeurs sont déjà séparées par des espaces
  return cleaned
    .split(/\s+/)
    .map(x => parseFloat(x.replace(",", ".")))
    .map(v => (Number.isFinite(v) ? v : null));
}

function parseAtomSymbolsFromMolblock(molblock) {
  if (!molblock) return [];
  const lines = molblock.split(/\r?\n/);
  const counts = lines[3] || "";
  const nAtoms = parseInt(counts.slice(0, 3).trim(), 10);
  const atoms = [];
  for (let i = 0; i < nAtoms; i++) {
    const l = lines[4 + i] || "";
    const sym = l.slice(31, 34).trim();
    atoms.push(sym || "?");
  }
  return atoms;
}

async function tryFetchText(url) {
  try {
    const r = await fetch(url, { method: "GET" });
    if (!r.ok) return null;
    const txt = await r.text();
    if (/<!doctype|<html/i.test(txt)) return null;
    return txt;
  } catch {
    return null;
  }
}

function dirnameFromUrl(u) {
  try {
    // En Electron, utiliser API_BASE_URL comme base
    const baseUrl = typeof API_BASE_URL === 'string' && API_BASE_URL.startsWith('http')
      ? API_BASE_URL
      : window.location.href;
    const url = new URL(u, baseUrl);
    const p = url.pathname;
    return `${url.origin}${p.slice(0, p.lastIndexOf("/") + 1)}`;
  } catch {
    return null;
  }
}

async function resolveMolForSource(csvUrl, sourceName, taskId = null) {
  if (!csvUrl || !sourceName) {
    debugLog('[resolveMolForSource] Empty csvUrl or sourceName:', { csvUrl, sourceName });
    return null;
  }
  const baseName = String(sourceName)
    .trim()
    .replace(/\.(sdf|mol)$/i, "");

  debugLog('[resolveMolForSource] Searching for:', { csvUrl, sourceName, taskId, baseName });

  if (taskId) {
    try {
      // En Electron, utiliser API_BASE_URL pour construire l'URL de base
      const baseUrl = typeof API_BASE_URL === 'string' && API_BASE_URL.startsWith('http')
        ? API_BASE_URL
        : window.location.href;
      const origin = new URL(csvUrl, baseUrl).origin;
      const base = `${origin}/outputs/${taskId}/`;
      const cand = [
        `${base}${baseName}.sdf`,
        `${base}${baseName}.mol`,
        `${base}${encodeURIComponent(baseName)}.sdf`,
        `${base}${encodeURIComponent(baseName)}.mol`,
      ];
      debugLog('[resolveMolForSource] Trying taskId candidates:', cand);
      for (const u of cand) {
        const txt = await tryFetchText(u);
        if (txt) {
          debugLog('[resolveMolForSource] ✅ Found at:', u);
          return {
            url: u,
            rawText: txt,
            firstRecord: extractFirstSDFRecord(txt),
            molblock: extractFirstMolBlockPreserving(txt) || txt,
          };
        }
      }
      debugLog('[resolveMolForSource] ❌ Not found via taskId');
    } catch (e) {
      debugLog('[resolveMolForSource] Error with taskId:', e?.message);
    }
  }

  const csvDir = dirnameFromUrl(csvUrl);
  if (csvDir) {
    const cand = [
      `${csvDir}${encodeURIComponent(baseName)}.sdf`,
      `${csvDir}${encodeURIComponent(baseName)}.mol`,
      `${csvDir}${baseName}.sdf`,
      `${csvDir}${baseName}.mol`,
      `${csvDir}${baseName}.SDF`,
      `${csvDir}${baseName}.MOL`,
    ];
    for (const u of cand) {
      const txt = await tryFetchText(u);
      if (txt)
        return {
          url: u,
          rawText: txt,
          firstRecord: extractFirstSDFRecord(txt),
          molblock: extractFirstMolBlockPreserving(txt) || txt,
        };
    }
  }

  try {
    // En Electron, utiliser API_BASE_URL pour construire l'URL de base
    const baseUrl = typeof API_BASE_URL === 'string' && API_BASE_URL.startsWith('http')
      ? API_BASE_URL
      : window.location.href;
    const origin = new URL(csvUrl, baseUrl).origin;
    const api = `${origin}/api/prediction-tasks`;
    const r = await fetch(api);
    if (r.ok) {
      const tasks = await r.json();
      for (const t of tasks || []) {
        if (taskId && t.taskId === taskId) continue;
        const base = `${origin}/outputs/${t.taskId}/`;
        const cand = [
          `${base}${baseName}.sdf`,
          `${base}${baseName}.mol`,
          `${base}${encodeURIComponent(baseName)}.sdf`,
          `${base}${encodeURIComponent(baseName)}.mol`,
        ];
        for (const u of cand) {
          const txt = await tryFetchText(u);
          if (txt)
            return {
              url: u,
              rawText: txt,
              firstRecord: extractFirstSDFRecord(txt),
              molblock: extractFirstMolBlockPreserving(txt) || txt,
            };
        }
      }
    }
  } catch {}

  return null;
}

/* ===========================
   Hook de mesure largeur
   =========================== */

function useMeasureWidth(initial = 600) {
  const ref = useRef(null);
  const [width, setWidth] = useState(initial);

  useEffect(() => {
    const el = ref.current;
    if (!el) return;
    const ro = new ResizeObserver((entries) => {
      for (const entry of entries) {
        const w = entry.contentRect?.width || el.clientWidth || initial;
        setWidth(Math.max(100, Math.floor(w)));
      }
    });
    ro.observe(el);
    return () => ro.disconnect();
  }, [initial]);

  return [ref, width];
}

/* ===========================
   Canvas 1D — StickChart
   =========================== */

function StickChart({
  theory = [],
  expPpm = [],
  expInt = [],
  reverse = true,
  height = PLOT_HEIGHT,
  showTheory = true,
  showExp = true,
  backgroundImage = null,
  backgroundMeta = null,
  canInteract = false,
  zoomEnabled = true,
  onZoom,
  onZoomPx,
  onReset,
  busy = false,
}) {
  const [wrapRef, width] = useMeasureWidth(700);
  const canvasRef = useRef(null);
  const [bgLoaded, setBgLoaded] = useState(false);
  const bgImageRef = useRef(null);
  const [drag, setDrag] = useState(null);

  const range = useMemo(() => {
    if (backgroundMeta) {
      if (
        backgroundMeta.ppm_min != null &&
        backgroundMeta.ppm_max != null
      ) {
        const a = Number(backgroundMeta.ppm_min);
        const b = Number(backgroundMeta.ppm_max);
        return { minX: Math.min(a, b), maxX: Math.max(a, b) };
      }
      if (
        backgroundMeta.ppm_left != null &&
        backgroundMeta.ppm_right != null
      ) {
        const a = Number(backgroundMeta.ppm_right);
        const b = Number(backgroundMeta.ppm_left);
        return { minX: Math.min(a, b), maxX: Math.max(a, b) };
      }
    }
    const all = [
      ...(showTheory ? theory : []),
      ...(showExp ? expPpm : []),
    ].filter(Number.isFinite);

    if (all.length === 0) return { minX: 0, maxX: 1 };
    let minX = Math.min(...all);
    let maxX = Math.max(...all);
    if (minX === maxX) { minX -= 0.01; maxX += 0.01; }
    return { minX, maxX };
  }, [theory, expPpm, showTheory, showExp, backgroundMeta]);

  const expSeries = useMemo(() => {
    const norm = (arr) => {
      if (!arr || arr.length === 0) return [];
      const m = Math.max(
        ...arr.map(v => (Number.isFinite(v) ? v : 0)),
        1
      );
      return arr.map(v => (Number.isFinite(v) ? v / m : 0.6));
    };
    const a = norm(expInt);
    return expPpm.map((x, i) => ({ x, a: a[i] ?? 0.6 }));
  }, [expPpm, expInt]);

  useEffect(() => {
    if (!backgroundImage) {
      setBgLoaded(false);
      bgImageRef.current = null;
      return;
    }
    const img = new Image();
    img.crossOrigin = "anonymous";

    const safeSrc = buildImgSrc(backgroundImage, backgroundMeta);

    let retriedWithoutCb = false;
    img.onload = () => { bgImageRef.current = img; setBgLoaded(true); };
    img.onerror = () => {
      if (!retriedWithoutCb && safeSrc !== backgroundImage) {
        retriedWithoutCb = true;
        img.src = backgroundImage;
        return;
      }
      bgImageRef.current = null; setBgLoaded(false);
    };
    img.src = safeSrc;
  }, [backgroundImage, backgroundMeta]);

  useEffect(() => {
    const cvs = canvasRef.current;
    if (!cvs) return;

    const DPR = window.devicePixelRatio || 1;
    const W = width;
    const H = height;
    const PAD = 30;
    const gridW = W - 2 * PAD;
    const gridH = H - 2 * PAD;

    cvs.width = Math.round(W * DPR);
    cvs.height = Math.round(H * DPR);
    const ctx = cvs.getContext("2d");
    if (DPR !== 1) ctx.scale(DPR, DPR);

    ctx.fillStyle = BG;
    ctx.fillRect(0, 0, W, H);

    // ✅ CORRECTION : Vérification image valide avant drawImage
    if (bgLoaded && bgImageRef.current) {
      const imgValid = 
        bgImageRef.current.complete &&
        bgImageRef.current.naturalWidth > 0 &&
        bgImageRef.current.naturalHeight > 0;
      
      if (!imgValid) {
        debugWarn('⚠️ Image 1D pas encore chargée, skip drawImage');
      } else {
        const bb = backgroundMeta?.plot_bbox_px;
        if (bb && [bb.x, bb.y, bb.w, bb.h].every(Number.isFinite)) {
          try {
            ctx.drawImage(
              bgImageRef.current,
              Math.max(0, bb.x),
              Math.max(0, bb.y),
              Math.max(1, bb.w),
              Math.max(1, bb.h),
              PAD,
              PAD,
              gridW,
              gridH
            );
          } catch (e) {
            debugError('❌ Erreur drawImage 1D:', e.message);
          }
        } else {
          try {
            ctx.drawImage(bgImageRef.current, PAD, PAD, gridW, gridH);
          } catch (e) {
            debugError('❌ Erreur drawImage 1D (full):', e.message);
          }
        }
        ctx.globalAlpha = 0.12;
        ctx.fillStyle = "#fff";
        ctx.fillRect(PAD, PAD, gridW, gridH);
        ctx.globalAlpha = 1.0;
      }
    }

    ctx.strokeStyle = "rgba(99,110,131,0.25)";
    ctx.lineWidth = 1;
    const gridY = 5;
    for (let i = 0; i <= gridY; i++) {
      const y = PAD + (i / gridY) * gridH;
      ctx.beginPath();
      ctx.moveTo(PAD, y);
      ctx.lineTo(PAD + gridW, y);
      ctx.stroke();
    }
    ctx.strokeStyle = AXIS;
    ctx.beginPath();
    ctx.moveTo(PAD, PAD + gridH);
    ctx.lineTo(PAD + gridW, PAD + gridH);
    ctx.stroke();

    const denom = range.maxX - range.minX || 1;

    const toX = (() => {
      const canMap =
        bgLoaded &&
        bgImageRef.current &&
        backgroundMeta &&
        backgroundMeta.ppm_left != null &&
        backgroundMeta.ppm_right != null;

      if (canMap) {
        const mapper = make1DMapper?.(backgroundMeta, gridW);
        if (mapper && mapper.axisRect?.w > 0) {
          return (ppmRaw) => {
            const ppm = Number(ppmRaw);
            if (!Number.isFinite(ppm)) return PAD;
            const xpx = mapper(ppm);
            const t =
              (xpx - mapper.axisRect.x) / mapper.axisRect.w;
            return PAD + t * gridW;
          };
        }
      }
      return (ppmRaw) => {
        const ppm = Number(ppmRaw);
        if (!Number.isFinite(ppm)) return PAD;
        const t = (ppm - range.minX) / denom;
        const s = reverse ? 1 - t : t;
        return PAD + s * gridW;
      };
    })();

    if (showTheory && theory?.length) {
      ctx.strokeStyle = COLOR_THEORY;
      ctx.lineWidth = 2.5;
      for (const ppm of theory) {
        const x = Math.round(toX(ppm)) + 0.5;
        ctx.beginPath();
        ctx.moveTo(x, PAD + gridH);
        ctx.lineTo(x, PAD + gridH * 0.95);
        ctx.stroke();
      }
    }

    if (showExp && expSeries?.length) {
      ctx.strokeStyle = COLOR_EXP;
      ctx.lineWidth = 2.5;
      for (const { x: ppm } of expSeries) {
        const xv = Math.round(toX(ppm)) + 0.5;
        ctx.beginPath();
        ctx.moveTo(xv, PAD + gridH);
        ctx.lineTo(xv, PAD + gridH * 0.95);
        ctx.stroke();
      }
    }

    ctx.fillStyle = AXIS;
    ctx.font = "12px ui-sans-serif, system-ui, -apple-system";
    const ticks = 6;
    for (let i = 0; i <= ticks; i++) {
      const t = i / ticks;
      const ppm = range.minX + t * (range.maxX - range.minX);
      const x = PAD + (reverse ? (1 - t) : t) * gridW;
      ctx.fillText(ppm.toFixed(2), x - 14, H - 8);
    }

    if (drag &&
        Number.isFinite(drag.x0) &&
        Number.isFinite(drag.x1)) {
      const x0 = Math.min(drag.x0, drag.x1);
      const x1 = Math.max(drag.x0, drag.x1);
      ctx.fillStyle = "rgba(37,99,235,0.15)";
      ctx.fillRect(x0, PAD, x1 - x0, gridH);
      ctx.strokeStyle = "rgba(37,99,235,0.65)";
      ctx.lineWidth = 2;
      ctx.strokeRect(
        x0 + 0.5, PAD + 0.5, x1 - x0 - 1, gridH - 1
      );
    }

    if (busy) {
      ctx.fillStyle = "rgba(255,255,255,0.85)";
      ctx.fillRect(PAD, PAD, gridW, gridH);
      ctx.fillStyle = "#64748b";
      ctx.font = "13px ui-sans-serif, system-ui, -apple-system";
      ctx.fillText("Mise à jour du fond…", PAD + 10, PAD + 22);
    }
  }, [
    width, height, reverse, theory, expSeries, range,
    showExp, showTheory, bgLoaded, backgroundMeta, drag, busy
  ]);

  const canDragNow = canInteract && zoomEnabled;

  const onMouseDown = useCallback((e) => {
    if (!canDragNow || e.button !== 0) return;
    e.preventDefault();
    const rect = e.currentTarget.getBoundingClientRect();
    const x = e.clientX - rect.left;
    setDrag({ x0: x, x1: x });
  }, [canDragNow]);

  const onMouseMove = useCallback((e) => {
    if (!drag) return;
    e.preventDefault();
    const rect = e.currentTarget.getBoundingClientRect();
    const x = e.clientX - rect.left;
    setDrag(d => ({ ...d, x1: x }));
  }, [drag]);

  const onMouseUp = useCallback((e) => {
    if (!drag) return;
    e.preventDefault();

    const PAD = 30;
    const rect = e.currentTarget.getBoundingClientRect();
    const W = rect.width;
    const H = rect.height;
    const gridW = W - 2 * PAD;
    const gridH = H - 2 * PAD;

    const eps = 0.5;
    const x0 = Math.max(
      PAD + eps,
      Math.min(W - PAD - eps, Math.min(drag.x0, drag.x1))
    );
    const x1 = Math.max(
      PAD + eps,
      Math.min(W - PAD - eps, Math.max(drag.x0, drag.x1))
    );

    setDrag(null);

    if (!canDragNow || x1 - x0 < 4) return;

    if (typeof onZoom === "function") {
      let left, right;
      if (
        Number.isFinite(+backgroundMeta?.ppm_left) &&
        Number.isFinite(+backgroundMeta?.ppm_right)
      ) {
        left = +backgroundMeta.ppm_left;
        right = +backgroundMeta.ppm_right;
      } else if (
        Number.isFinite(+backgroundMeta?.ppm_min) &&
        Number.isFinite(+backgroundMeta?.ppm_max)
      ) {
        const a = +backgroundMeta.ppm_min;
        const b = +backgroundMeta.ppm_max;
        left = Math.max(a, b);
        right = Math.min(a, b);
      }

      if (Number.isFinite(left) && Number.isFinite(right)) {
        const frac = (xc) => Math.max(0, Math.min(1, (xc - PAD) / gridW));
        const ppmAt = (xc) => left - frac(xc) * (left - right);
        const pA = ppmAt(x0);
        const pB = ppmAt(x1);
        onZoom(Math.min(pA, pB), Math.max(pA, pB));
        return;
      }
    }

    if (typeof onZoomPx === "function") {
      const rectPxDOM = {
        x: Math.min(x0, x1) - PAD,
        y: 0,
        w: Math.abs(x1 - x0),
        h: gridH,
      };
      if (rectPxDOM.w >= 1) {
        onZoomPx(rectPxDOM, {
          displayWidth: gridW,
          displayHeight: gridH,
        });
        return;
      }
    }

    if (
      range &&
      Number.isFinite(range.minX) &&
      Number.isFinite(range.maxX) &&
      typeof onZoom === "function"
    ) {
      const t0 = (x0 - PAD) / gridW;
      const t1 = (x1 - PAD) / gridW;
      const width = range.maxX - range.minX;
      const ppmA = reverse
        ? (range.maxX - t0 * width)
        : (range.minX + t0 * width);
      const ppmB = reverse
        ? (range.maxX - t1 * width)
        : (range.minX + t1 * width);
      onZoom(Math.min(ppmA, ppmB), Math.max(ppmA, ppmB));
    }
  }, [drag, canDragNow, backgroundMeta, onZoom, onZoomPx, range, reverse]);

  useEffect(() => {
    const onKey = (e) => { if (e.key === "Escape") setDrag(null); };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, []);

  return (
    <div
      ref={wrapRef}
      style={{ width: "100%", position: "relative", userSelect: "none" }}
    >
      <canvas
        ref={canvasRef}
        style={{
          width: "100%",
          height,
          borderRadius: 8,
          background: BG,
          boxShadow: "0 1px 4px rgba(0,0,0,.08)",
          display: "block",
          cursor: canDragNow
            ? (drag ? "ew-resize" : "crosshair")
            : "default",
        }}
        onMouseDown={onMouseDown}
        onMouseMove={onMouseMove}
        onMouseUp={onMouseUp}
        onMouseLeave={() => drag && setDrag(null)}
      />
      {backgroundImage && !bgLoaded && (
        <div
          style={{
            position: "absolute",
            top: 10,
            right: 10,
            fontSize: 12,
            color: "#666",
            background: "rgba(255,255,255,0.85)",
            padding: 4,
            borderRadius: 4,
          }}
        >
          Chargement fond 1D...
        </div>
      )}
      {canInteract && (
        <div
          style={{
            position: "absolute",
            top: 8,
            left: 8,
            fontSize: 12,
            color: "#475569",
          }}
        >
          Glisser pour zoomer •{" "}
          <button
            onClick={onReset}
            style={{
              marginLeft: 6,
              padding: "2px 6px",
              borderRadius: 6,
              border: "1px solid #cbd5e1",
              background: "#fff",
              cursor: "pointer",
            }}
          >
            Reset Zoom
          </button>
        </div>
      )}
    </div>
  );
}
/* ===========================
   Canvas 2D — Scatter2D
   =========================== */

function Scatter2D({
  theoryPts = [],
  expPts = [],
  expInt = [],
  showTheory = true,
  showExp = true,
  height = PLOT_HEIGHT,
  maxSize = 10,
  backgroundImage = null,
  backgroundMeta = null,
  canInteract = false,
  zoomEnabled = true,
  onZoom2D,
  onWheelDelta,
  onReset,
  busy = false,
}) {
  const [wrapRef, width] = useMeasureWidth(720);
  const canvasRef = useRef(null);
  const [bgLoaded, setBgLoaded] = useState(false);
  const bgImageRef = useRef(null);
  const [drag, setDrag] = useState(null);

  const ranges = useMemo(() => {
    if (
      backgroundMeta?.f1_top !== undefined &&
      backgroundMeta?.f1_bottom !== undefined &&
      backgroundMeta?.f2_left !== undefined &&
      backgroundMeta?.f2_right !== undefined
    ) {
      return {
        f1: [
          Math.min(backgroundMeta.f1_bottom, backgroundMeta.f1_top),
          Math.max(backgroundMeta.f1_bottom, backgroundMeta.f1_top),
        ],
        f2: [
          Math.min(backgroundMeta.f2_left, backgroundMeta.f2_right),
          Math.max(backgroundMeta.f2_left, backgroundMeta.f2_right),
        ],
      };
    }
    if (
      backgroundMeta &&
      backgroundMeta.f1_ppm_min != null &&
      backgroundMeta.f1_ppm_max != null &&
      backgroundMeta.f2_ppm_min != null &&
      backgroundMeta.f2_ppm_max != null
    ) {
      return {
        f1: [
          Math.min(backgroundMeta.f1_ppm_min, backgroundMeta.f1_ppm_max),
          Math.max(backgroundMeta.f1_ppm_min, backgroundMeta.f1_ppm_max),
        ],
        f2: [
          Math.min(backgroundMeta.f2_ppm_min, backgroundMeta.f2_ppm_max),
          Math.max(backgroundMeta.f2_ppm_min, backgroundMeta.f2_ppm_max),
        ],
      };
    }
    const pts = [];
    if (showTheory) pts.push(...theoryPts);
    if (showExp) pts.push(...expPts);
    if (pts.length === 0) return { f1: [0, 1], f2: [0, 1] };
    const f1s = pts.map(p => p[0]).filter(Number.isFinite);
    const f2s = pts.map(p => p[1]).filter(Number.isFinite);
    let minF1 = Math.min(...f1s);
    let maxF1 = Math.max(...f1s);
    let minF2 = Math.min(...f2s);
    let maxF2 = Math.max(...f2s);
    if (minF1 === maxF1) { minF1 -= 0.001; maxF1 += 0.001; }
    if (minF2 === maxF2) { minF2 -= 0.001; maxF2 += 0.001; }
    return { f1: [minF1, maxF1], f2: [minF2, maxF2] };
  }, [theoryPts, expPts, showTheory, showExp, backgroundMeta]);

  const normInt = useMemo(() => {
    if (!showExp || !expInt?.length) return null;
    const m = Math.max(
      ...expInt.map(v => (Number.isFinite(v) ? v : 0)),
      1
    );
    return expInt.map(v => (Number.isFinite(v) ? v / m : 0.4));
  }, [expInt, showExp]);

  useEffect(() => {
    if (!backgroundImage) {
      setBgLoaded(false);
      bgImageRef.current = null;
      return;
    }
    const img = new Image();
    img.crossOrigin = "anonymous";

    const safeSrc = buildImgSrc(backgroundImage, backgroundMeta);

    let retriedWithoutCb = false;
    img.onload = () => { bgImageRef.current = img; setBgLoaded(true); };
    img.onerror = () => {
      if (!retriedWithoutCb && safeSrc !== backgroundImage) {
        retriedWithoutCb = true;
        img.src = backgroundImage;
        return;
      }
      bgImageRef.current = null; setBgLoaded(false);
    };
    img.src = safeSrc;
  }, [backgroundImage, backgroundMeta]);

  useEffect(() => {
    const cvs = canvasRef.current;
    if (!cvs) return;

    const DPR = window.devicePixelRatio || 1;
    const W = width;
    const H = height;
    const PAD = 36;
    const gridW = W - 2 * PAD;
    const gridH = H - 2 * PAD;

    cvs.width = Math.round(W * DPR);
    cvs.height = Math.round(H * DPR);
    const ctx = cvs.getContext("2d");
    if (DPR !== 1) ctx.scale(DPR, DPR);

    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, W, H);

    // ✅ CORRECTION : Vérification image valide avant drawImage
    if (bgLoaded && bgImageRef.current) {
      const imgValid = 
        bgImageRef.current.complete &&
        bgImageRef.current.naturalWidth > 0 &&
        bgImageRef.current.naturalHeight > 0;
      
      if (!imgValid) {
        debugWarn('⚠️ Image 2D pas encore chargée, skip drawImage');
      } else {
        const bb =
          backgroundMeta?.plot_bbox_px ||
          backgroundMeta?.axis_bbox_px ||
          backgroundMeta?.axisRect?.pixel ||
          backgroundMeta?.axisRect?.px ||
          (backgroundMeta?.axisRect &&
           ['x','y','w','h'].every(k => Number.isFinite(backgroundMeta.axisRect[k]))
             ? backgroundMeta.axisRect : null);

        if (bb && [bb.x, bb.y, bb.w, bb.h].every(Number.isFinite)) {
          try {
            ctx.drawImage(
              bgImageRef.current,
              Math.max(0, bb.x),
              Math.max(0, bb.y),
              Math.max(1, bb.w),
              Math.max(1, bb.h),
              PAD, PAD, gridW, gridH
            );
          } catch (e) {
            debugError('❌ Erreur drawImage 2D:', e.message);
          }
        } else {
          try {
            ctx.drawImage(bgImageRef.current, PAD, PAD, gridW, gridH);
          } catch (e) {
            debugError('❌ Erreur drawImage 2D (full):', e.message);
          }
        }
        ctx.globalAlpha = 0.18;
        ctx.fillStyle = "#fff";
        ctx.fillRect(PAD, PAD, gridW, gridH);
        ctx.globalAlpha = 1.0;
      }
    }

    const [minF1, maxF1] = ranges.f1;
    const [minF2, maxF2] = ranges.f2;

    const { toX, toY } = (() => {
      const clamp01 = (t) => Math.max(0, Math.min(1, t));
      const canMap =
        bgLoaded &&
        bgImageRef.current &&
        backgroundMeta &&
        backgroundMeta.f1_top != null &&
        backgroundMeta.f1_bottom != null &&
        backgroundMeta.f2_left != null &&
        backgroundMeta.f2_right != null;

      if (canMap) {
        const f2_left   = Number(backgroundMeta.f2_left);
        const f2_right  = Number(backgroundMeta.f2_right);
        const f1_top    = Number(backgroundMeta.f1_top);
        const f1_bottom = Number(backgroundMeta.f1_bottom);
        const dx = (f2_left - f2_right) || 1e-9;
        const dy = (f1_top  - f1_bottom) || 1e-9;
        return {
          toX: (f2) => {
            const t = (f2_left - Number(f2)) / dx;
            return PAD + clamp01(t) * gridW;
          },
          toY: (f1) => {
            const t = (Number(f1) - f1_bottom) / dy;
            return PAD + clamp01(t) * gridH;
          },
        };
      }

      const dx = (maxF2 - minF2) || 1e-9;
      const dy = (maxF1 - minF1) || 1e-9;
      return {
        toX: (f2) => {
          const t = (f2 - minF2) / dx;
          return PAD + (1 - t) * gridW;
        },
        toY: (f1) => {
          const t = (f1 - minF1) / dy;
          return PAD + t * gridH;
        },
      };
    })();

    // Grille
    ctx.strokeStyle = "rgba(99,110,131,.20)";
    ctx.lineWidth = 1;
    const steps = 10;
    for (let i = 0; i <= steps; i++) {
      const gx = PAD + (i / steps) * gridW;
      const gy = PAD + (i / steps) * gridH;
      ctx.beginPath(); ctx.moveTo(gx, PAD); ctx.lineTo(gx, PAD + gridH);
      ctx.stroke();
      ctx.beginPath(); ctx.moveTo(PAD, gy); ctx.lineTo(PAD + gridW, gy);
      ctx.stroke();
    }

    // Points théorie
    if (showTheory && theoryPts?.length) {
      ctx.fillStyle = "#2563eb";
      ctx.strokeStyle = "rgba(255,255,255,.85)";
      ctx.lineWidth = 1.5;
      for (const [f1, f2] of theoryPts) {
        const x = toX(f2);
        const y = toY(f1);
        ctx.beginPath();
        ctx.arc(x, y, 4.5, 0, Math.PI * 2);
        ctx.fill();
        ctx.stroke();
      }
    }

    // Points expé
    if (showExp && expPts?.length) {
      ctx.strokeStyle = "rgba(255,255,255,.85)";
      ctx.lineWidth = 1.5;
      for (let i = 0; i < expPts.length; i++) {
        const [f1, f2] = expPts[i];
        const a = (normInt ? normInt[i] : 0.8);
        const r = Math.max(3, a * maxSize);
        const x = toX(f2);
        const y = toY(f1);
        ctx.beginPath();
        ctx.fillStyle = "rgba(239, 68, 68, 0.6)";
        ctx.arc(x, y, r, 0, Math.PI * 2);
        ctx.fill();
        ctx.stroke();
      }
    }

    // ✅ Ticks F1 (à gauche) - valeurs uniquement
    ctx.fillStyle = AXIS;
    ctx.font = "12px ui-sans-serif, system-ui, -apple-system";
    const ticks = 6;
    for (let i = 0; i <= ticks; i++) {
      const t = i / ticks;
      const v1 = minF1 + t * (maxF1 - minF1);
      const y1 = toY(v1);
      ctx.fillText(v1.toFixed(2), 6, y1 + 4);
    }

    // ✅ Ticks F2 (en bas) - valeurs uniquement
    for (let i = 0; i <= ticks; i++) {
      const t = i / ticks;
      const v2 = minF2 + t * (maxF2 - minF2);
      const x2 = toX(v2);
      ctx.fillText(v2.toFixed(2), x2 - 14, H - 8);
    }

    // ✅ CORRECTION : Labels repositionnés EN DEHORS de la zone des ticks
    ctx.font = "bold 13px ui-sans-serif, system-ui, -apple-system";

    // Label F2 (en bas, complètement à droite)
    ctx.fillText("F2 (ppm)", PAD + gridW - 55, H - 22);

    // Label F1 (vertical, en haut à gauche)
    ctx.save();
    ctx.translate(8, PAD + 50);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText("F1 (ppm)", 0, 0);
    ctx.restore();

    // Rectangle de zoom
    if (drag) {
      const x = Math.min(drag.x0, drag.x1);
      const y = Math.min(drag.y0, drag.y1);
      const w = Math.abs(drag.x1 - drag.x0);
      const h = Math.abs(drag.y1 - drag.y0);
      ctx.fillStyle = "rgba(34,197,94,0.15)";
      ctx.fillRect(x, y, w, h);
      ctx.strokeStyle = "rgba(34,197,94,0.65)";
      ctx.lineWidth = 2;
      ctx.strokeRect(x + 0.5, y + 0.5, Math.max(0, w - 1), Math.max(0, h - 1));
    }

    // Overlay busy
    if (busy) {
      ctx.fillStyle = "rgba(255,255,255,0.85)";
      ctx.fillRect(PAD, PAD, gridW, gridH);
      ctx.fillStyle = "#64748b";
      ctx.font = "13px ui-sans-serif, system-ui, -apple-system";
      ctx.fillText("Mise à jour du fond…", PAD + 10, PAD + 22);
    }
  }, [
    width, height, ranges, theoryPts, expPts, normInt,
    showTheory, showExp, maxSize, bgLoaded, backgroundMeta, drag, busy
  ]);

  const wheelTimeoutRef = useRef(null);

  useEffect(() => {
    const cvs = canvasRef.current;
    if (!cvs || !canInteract) return;
    
    const handler = (e) => {
      e.preventDefault();
      e.stopPropagation();
      const dy = e.deltaY;
      window.clearTimeout(wheelTimeoutRef.current);
      wheelTimeoutRef.current = window.setTimeout(() => {
        onWheelDelta && onWheelDelta(dy);
      }, 60);
    };
    
    cvs.addEventListener('wheel', handler, { passive: false });
    return () => cvs.removeEventListener('wheel', handler);
  }, [canInteract, onWheelDelta]);

  const onMouseDown = useCallback((e) => {
    if (!canInteract || !zoomEnabled) return;
    const rect = e.currentTarget.getBoundingClientRect();
    setDrag({
      x0: e.clientX - rect.left,
      y0: e.clientY - rect.top,
      x1: e.clientX - rect.left,
      y1: e.clientY - rect.top,
    });
  }, [canInteract, zoomEnabled]);

  const onMouseMove = useCallback((e) => {
    if (!drag) return;
    const rect = e.currentTarget.getBoundingClientRect();
    setDrag((d) => ({
      ...d,
      x1: e.clientX - rect.left,
      y1: e.clientY - rect.top,
    }));
  }, [drag]);

  const onMouseUp = useCallback((e) => {
    if (!drag) return;
    const PAD = 36;
    const rect = e.currentTarget.getBoundingClientRect();
    const W = rect.width;
    const H = rect.height;
    const gridW = W - 2 * PAD;
    const gridH = H - 2 * PAD;

    const x0 = Math.max(PAD, Math.min(W - PAD, Math.min(drag.x0, drag.x1)));
    const x1 = Math.max(PAD, Math.min(W - PAD, Math.max(drag.x0, drag.x1)));
    const y0 = Math.max(PAD, Math.min(H - PAD, Math.min(drag.y0, drag.y1)));
    const y1 = Math.max(PAD, Math.min(H - PAD, Math.max(drag.y0, drag.y1)));

    setDrag(null);

    if (!canInteract || !zoomEnabled || x1 - x0 < 8 || y1 - y0 < 8) {
      return;
    }

    const tX0 = (x0 - PAD) / gridW;
    const tX1 = (x1 - PAD) / gridW;
    const tY0 = (y0 - PAD) / gridH;
    const tY1 = (y1 - PAD) / gridH;

    if (
      backgroundMeta?.f2_left != null &&
      backgroundMeta?.f2_right != null &&
      backgroundMeta?.f1_top != null &&
      backgroundMeta?.f1_bottom != null
    ) {
      const f2_left = Number(backgroundMeta.f2_left);
      const f2_right = Number(backgroundMeta.f2_right);
      const f1_top = Number(backgroundMeta.f1_top);
      const f1_bottom = Number(backgroundMeta.f1_bottom);

      const denomX = (f2_left - f2_right) || 1e-9;
      const denomY = (f1_top - f1_bottom) || 1e-9;

      const f2_A = f2_left - tX0 * denomX;
      const f2_B = f2_left - tX1 * denomX;
      const f1_A = f1_bottom + tY0 * denomY;
      const f1_B = f1_bottom + tY1 * denomY;

      onZoom2D &&
        onZoom2D({
          f2_min: Math.min(f2_A, f2_B),
          f2_max: Math.max(f2_A, f2_B),
          f1_min: Math.min(f1_A, f1_B),
          f1_max: Math.max(f1_A, f1_B),
        });
    }
  }, [drag, canInteract, zoomEnabled, backgroundMeta, onZoom2D]);

  return (
    <div
      ref={wrapRef}
      style={{ width: "100%", position: "relative", userSelect: "none",
               overscrollBehavior: "contain" }}
    >
      <canvas
        ref={canvasRef}
        style={{
          width: "100%",
          height,
          borderRadius: 8,
          background: "#ffffff",
          boxShadow: "0 1px 4px rgba(0,0,0,.08)",
          display: "block",
          cursor: canInteract && zoomEnabled ? "crosshair" : "default",
        }}
        onMouseDown={onMouseDown}
        onMouseMove={onMouseMove}
        onMouseUp={onMouseUp}
        onMouseLeave={() => drag && setDrag(null)}
      />
      {backgroundImage && !bgLoaded && (
        <div
          style={{
            position: "absolute",
            top: 10,
            right: 10,
            fontSize: 12,
            color: "#666",
            background: "rgba(255,255,255,0.85)",
            padding: 4,
            borderRadius: 4,
          }}
        >
          Chargement fond 2D...
        </div>
      )}
      {canInteract && (
        <div
          style={{
            position: "absolute",
            top: 8,
            left: 8,
            fontSize: 12,
            color: "#475569",
          }}
        >
          Molette = niveaux • Glisser = zoom •{" "}
          <button
            onClick={onReset}
            style={{
              marginLeft: 6,
              padding: "2px 6px",
              borderRadius: 6,
              border: "1px solid #cbd5e1",
              background: "#fff",
              cursor: "pointer",
            }}
          >
            Reset Zoom
          </button>
        </div>
      )}
    </div>
  );
}
function MoleculePane({
  csvUrl,
  selectedRow,
  show3D,
  taskId,
  multi = null,
  customColors,
  molZoom = 1,
  molHeight = MOL_HEIGHT,
  flagsOverride = null,            // coches globales (optionnelles)
  showLocalToggles = true,         // masquer/afficher les coches locales
  compact = false,                 // mode miniatures: juste le SVG, sans cadres/titres
}) {
  const [svg, setSvg] = useState(null);
  const [err, setErr] = useState(null);
  const [loading, setLoading] = useState(false);

  // Coches locales (utilisées seulement si flagsOverride n'est pas fourni)
  const [localFlags, setLocalFlags] = useState({
    H1: true, C13: true, COSY: true, HSQC: true,
  });
  const flags = flagsOverride || localFlags;       // source unique de vérité
  const toggle = (k) => setLocalFlags(f => ({ ...f, [k]: !f[k] }));

  const div3DRef = useRef(null);
  const viewerRef = useRef(null);

  // eslint-disable-next-line react-hooks/exhaustive-deps
  const COLORS_RDK = {
    H1:  customColors.H1.rgb,
    C13: customColors.C13.rgb,
    COSY:customColors.COSY.rgb,
    HSQC:customColors.HSQC.rgb,
  };
  // eslint-disable-next-line react-hooks/exhaustive-deps
  const COLORS_3D = {
    H1:  customColors.H1.hex,
    C13: customColors.C13.hex,
    COSY:customColors.COSY.hex,
    HSQC:customColors.HSQC.hex,
  };

  // --- HSQC: colorer la liaison C–H en 2D -----------------------------------
  const computeMultiHighlights = useCallback((resolved, rowsByType) => {
    const atoms  = parseAtomSymbolsFromMolblock(resolved.molblock) || [];
    const shifts = extractShiftArrayFromSDFRecord(resolved.rawText) || [];
    
    // ✅ CORRECTION : Logs conditionnels uniquement en mode debug
    debugLog("🔍 Shifts extraits:", shifts);
    debugLog("🔍 Symboles d'atomes:", atoms);
    debugLog("🔍 Nombre d'atomes:", atoms.length);
    
    // Protection contre les tableaux vides
    if (!atoms.length || !shifts.length) {
      debugWarn("⚠️ Pas d'atomes ou de shifts trouvés, coloration désactivée");
      return { 
        highlightByType: { H1: new Set(), C13: new Set(), COSY: new Set(), HSQC: new Set() }, 
        bondsCOSY: [], 
        bondsHSQC: [] 
      };
    }
    
    const maxIdx = Math.min(shifts.length - 1, atoms.length - 1);

    const hIdx   = atoms.map((s,i)=>s==="H"?i:-1).filter(i=>i>=0);
    const cIdx   = atoms.map((s,i)=>s==="C"?i:-1).filter(i=>i>=0);

    const sets   = { H1:new Set(), C13:new Set(), COSY:new Set(), HSQC:new Set() };
    let bondsCOSY = [];
    let bondsHSQC = [];

    // 1H
    if (flags.H1 && rowsByType["1H"]) {
      const m = (rowsByType["1H"].ppm_molecule_matched||"")
        .split("|").map(v=>parseFloat(String(v).replace(",","."))).filter(Number.isFinite);
      const r = new Set(m.map(v=>Math.round(v*1000)/1000));
      for (const ppm of r) for (const i of hIdx)
        if (i<=maxIdx && (Math.round(shifts[i]*1000)/1000)===ppm) sets.H1.add(i);
    }
    
    // 13C
    if (flags.C13 && rowsByType["13C"]) {
      const m = (rowsByType["13C"].ppm_molecule_matched||"")
        .split("|").map(v=>parseFloat(String(v).replace(",","."))).filter(Number.isFinite);
      const r = new Set(m.map(v=>Math.round(v*1000)/1000));
      for (const ppm of r) for (const i of cIdx)
        if (i<=maxIdx && (Math.round(shifts[i]*1000)/1000)===ppm) sets.C13.add(i);
    }
    
    // COSY (chemins min)
    if (flags.COSY && rowsByType["COSY"]) {
      const pairs = (rowsByType["COSY"].ppm_molecule_matched||"")
        .split("||").map(p=>p.split("|")).filter(a=>a.length===2)
        .map(([f1,f2])=>[+String(f1).replace(",","."),+String(f2).replace(",","." )])
        .filter(([a,b])=>Number.isFinite(a)&&Number.isFinite(b));
      const matched = [];
      for (const [d1,d2] of pairs) {
        const r1=Math.round(d1*1000)/1000, r2=Math.round(d2*1000)/1000;
        const a1=[], a2=[];
        for (const i of hIdx) if (i<=maxIdx&&(Math.round(shifts[i]*1000)/1000)===r1) a1.push(i);
        for (const i of hIdx) if (i<=maxIdx&&(Math.round(shifts[i]*1000)/1000)===r2) a2.push(i);
        for (const x of a1) for (const y of a2) if (x!==y) matched.push([x,y]);
      }
      if (matched.length) {
        bondsCOSY = findBondsInPath(resolved.molblock, matched);
        const lines = resolved.molblock.split(/\r?\n/);
        const nAtoms = parseInt((lines[3]||"").slice(0,3).trim(),10);
        for (const b of bondsCOSY) {
          const l = lines[4+nAtoms+b]||"";
          const a1 = parseInt(l.slice(0,3).trim(),10)-1;
          const a2 = parseInt(l.slice(3,6).trim(),10)-1;
          if (a1>=0) sets.COSY.add(a1);
          if (a2>=0) sets.COSY.add(a2);
        }
      }
    }
    
    // HSQC (liaison directe C–H)
    if (flags.HSQC && rowsByType["HSQC"]) {
      const pairs = (rowsByType["HSQC"].ppm_molecule_matched_unfolded ||
                     rowsByType["HSQC"].ppm_molecule_matched || "")
        .split("||").map(p=>p.split("|")).filter(a=>a.length===2)
        .map(([c,h])=>[+String(c).replace(",","."), +String(h).replace(",","." )])
        .filter(([c,h])=>Number.isFinite(c)&&Number.isFinite(h));

      const lines = resolved.molblock.split(/\r?\n/);
      const counts = lines[3] || "";
      const nAtoms = parseInt(counts.slice(0,3).trim(),10);
      const nBonds = parseInt(counts.slice(3,6).trim(),10);
      const bondMap = new Map();
      for (let bi=0; bi<nBonds; bi++) {
        const l = lines[4+nAtoms+bi] || "";
        const a1 = parseInt(l.slice(0,3).trim(),10)-1;
        const a2 = parseInt(l.slice(3,6).trim(),10)-1;
        if (a1>=0 && a2>=0) {
          const key = `${Math.min(a1,a2)}-${Math.max(a1,a2)}`;
          bondMap.set(key, bi);
        }
      }

      for (const [cppm,hppm] of pairs) {
        const rc = Math.round(cppm*1000)/1000;
        const rh = Math.round(hppm*1000)/1000;

        const cCandidates = [];
        const hCandidates = [];

        for (const ci of cIdx)
          if (ci<=maxIdx && (Math.round(shifts[ci]*1000)/1000)===rc) cCandidates.push(ci);

        for (const hi of hIdx)
          if (hi<=maxIdx && (Math.round(shifts[hi]*1000)/1000)===rh) {
            hCandidates.push(hi);
            sets.HSQC.add(hi);
          }

        for (const ci of cCandidates) for (const hi of hCandidates) {
          const key = `${Math.min(ci,hi)}-${Math.max(ci,hi)}`;
          const b = bondMap.get(key);
          if (b!=null) bondsHSQC.push(b);
        }
      }
    }
    
    return { highlightByType: sets, bondsCOSY, bondsHSQC };
  }, [flags]);
  // ---------------------------------------------------------------------------

  const render2D = useCallback(async (resolved) => {
    setLoading(true);
    try {
      const RDKit = await loadRDKitFromCDN();
      const mol = RDKit.get_mol(resolved.molblock);
      try { mol.generate_aligned_coords(); } catch { try { mol.set_new_coords(); } catch {} }

      let highlightAtoms=[], atomColorsMap={}, highlightBonds=[], bondColorsMap={};

      if (multi?.fromAll && multi?.rowsByType) {
        const { highlightByType, bondsCOSY, bondsHSQC } =
          computeMultiHighlights(resolved, multi.rowsByType);

        // Couleurs atomes (priorité HSQC > COSY > H1 pour les H)
        if (flags.C13)
          for (const i of highlightByType.C13)
            atomColorsMap[String(i)] = COLORS_RDK.C13;

        const h1List = Array.from(highlightByType.H1);
        const cosySet = highlightByType.COSY;
        const hsqcSet = highlightByType.HSQC;

        for (const i of h1List) {
          if (flags.HSQC && hsqcSet.has(i))      atomColorsMap[String(i)] = COLORS_RDK.HSQC;
          else if (flags.COSY && cosySet.has(i)) atomColorsMap[String(i)] = COLORS_RDK.COSY;
          else if (flags.H1)                     atomColorsMap[String(i)] = COLORS_RDK.H1;
        }

        // Atomes à surligner
        const all = new Set([
          ...(flags.H1   ? Array.from(highlightByType.H1)   : []),
          ...(flags.C13  ? Array.from(highlightByType.C13)  : []),
          ...(flags.COSY ? Array.from(highlightByType.COSY) : []),
          ...(flags.HSQC ? Array.from(highlightByType.HSQC) : []),
        ]);
        highlightAtoms = Array.from(all);

        // Liaisons (COSY chemins + HSQC liaison directe C–H)
        if (flags.COSY) {
          highlightBonds.push(...bondsCOSY);
          for (const b of bondsCOSY) bondColorsMap[String(b)] = COLORS_RDK.COSY;
        }
        if (flags.HSQC) {
          highlightBonds.push(...bondsHSQC);
          for (const b of bondsHSQC) bondColorsMap[String(b)] = COLORS_RDK.HSQC;
        }
        // Dédoublonnage
        highlightBonds = Array.from(new Set(highlightBonds));
      }

      const base = show3D ? 350 : 420;
      const sizePx = Math.round(base * Math.max(0.6, molZoom));
      const pad = Math.max(0.04, 0.24 / Math.max(1, molZoom));

      const opts = {
        width: sizePx, height: sizePx, useMolBlockCoords: false,
        padding: pad, bondLineWidth: 1.4, atomLabelFontSize: 10, annotationFontScale: .9,
      };

      let svgStr;
      try {
        svgStr = mol.get_svg_with_highlights(JSON.stringify({
          ...opts,
          atoms: highlightAtoms,
          highlight_atom_colors: atomColorsMap,
          bonds: highlightBonds,
          highlight_bond_colors: bondColorsMap,
          highlightRadius: 0.35,
          explicitMethyl: true,
          addStereoAnnotation: true,
          atomHighlightsAreCircles: true,
          highlightBondWidthMultiplier: 1.5,
        }));
      } catch {
        svgStr = mol.get_svg(JSON.stringify(opts));
      }

      svgStr = svgStr
        .replace(/width="[^"]+"/, 'width="100%"')
        .replace(/height="[^"]+"/, 'height="100%"')
        .replace('<svg ', '<svg preserveAspectRatio="xMidYMid meet" ');

      setSvg(svgStr);
      setErr(null);
      mol.delete?.();
    } catch (e) {
      setErr("Erreur 2D : " + (e?.message || String(e)));
    } finally {
      setLoading(false);
    }
  }, [multi, computeMultiHighlights, COLORS_RDK, molZoom, show3D, flags]);

  const render3D = useCallback(async (resolved) => {
    setLoading(true);
    try {
      // ✅ Vérifier WebGL avant d'essayer de créer le viewer
      if (!isWebGLAvailable()) {
        throw new Error("WebGL n'est pas disponible sur ce navigateur.\n\n" +
          "Solutions :\n" +
          "1. Firefox : Tapez about:config dans la barre d'adresse, puis réglez webgl.disabled sur false\n" +
          "2. Mettez à jour vos pilotes graphiques\n" +
          "3. Essayez Chrome/Edge si le problème persiste");
      }

      const $3Dmol = await load3DmolFromCDN();
      const el = div3DRef.current;
      
      // ✅ VÉRIFICATION : Le conteneur doit exister
      if (!el) {
        throw new Error("Conteneur 3D introuvable");
      }
      
      // ✅ VÉRIFICATION : Le conteneur doit avoir une taille
      const elWidth = el.offsetWidth;
      const elHeight = el.offsetHeight;
      
      if (elWidth === 0 || elHeight === 0) {
        debugWarn(`⚠️ Conteneur 3D pas encore visible (${elWidth}x${elHeight}px), attente...`);
        setTimeout(() => {
          if (div3DRef.current) {
            render3D(resolved);
          }
        }, 100);
        return;
      }
      
      if (!viewerRef.current) {
        viewerRef.current = (typeof $3Dmol.createViewer === "function")
          ? $3Dmol.createViewer(el, { 
              backgroundColor:"white", 
              antialias:false, 
              preserveDrawingBuffer:true, 
              disableOffscreenCanvas:true 
            })
          : new $3Dmol.GLViewer(el, { backgroundColor:"white" });
      }
      
      const viewer = viewerRef.current;
      
      // ✅ VÉRIFICATION : Le viewer doit être valide
      if (!viewer) {
        throw new Error("Échec de création du viewer 3D");
      }
      
      viewer.removeAllModels();
      const sdf = resolved.firstRecord || resolved.rawText || resolved.molblock;
      const model = viewer.addModel(sdf, "sdf");
      
      // ✅ VÉRIFICATION CRITIQUE : Le modèle doit exister
      if (!model) {
        throw new Error("Échec de chargement du modèle 3D");
      }
      
      model.setStyle({}, { 
        stick:{radius:0.15,color:"#808080"}, 
        sphere:{scale:0.2,color:"#808080"} 
      });

      if (multi?.fromAll && multi?.rowsByType) {
        const { highlightByType } = computeMultiHighlights(resolved, multi.rowsByType);
        const CAtoms = flags.C13 ? Array.from(highlightByType.C13) : [];
        const cosySet = flags.COSY ? new Set(highlightByType.COSY) : new Set();
        const hsqcSet = flags.HSQC ? new Set(highlightByType.HSQC) : new Set();
        const h1List  = flags.H1   ? Array.from(highlightByType.H1)  : [];

        // ✅ Colorier les atomes C13
        for (const i of CAtoms) {
          if (typeof i === 'number' && i >= 0) {
            try {
              model.setStyle({ index:i }, { 
                stick:{radius:0.15,color:"#808080"}, 
                sphere:{radius:0.4,color:COLORS_3D.C13} 
              });
            } catch (e) {
              debugWarn(`⚠️ Erreur setStyle C13 pour atome ${i}:`, e.message);
            }
          }
        }

        // ✅ Colorier les atomes COSY
        if (flags.COSY) {
          for (const i of cosySet) {
            if (typeof i === 'number' && i >= 0) {
              try {
                model.setStyle({ index:i }, { 
                  stick:{radius:0.25,color:COLORS_3D.COSY}, 
                  sphere:{radius:0.4,color:CAtoms.includes(i)?COLORS_3D.C13:"#808080"} 
                });
              } catch (e) {
                debugWarn(`⚠️ Erreur setStyle COSY pour atome ${i}:`, e.message);
              }
            }
          }
        }

        // ✅ Colorier les atomes HSQC
        if (flags.HSQC) {
          for (const i of hsqcSet) {
            if (typeof i === 'number' && i >= 0) {
              try {
                model.setStyle({ index:i }, { 
                  stick:{radius:0.2,color:COLORS_3D.HSQC}, 
                  sphere:{radius:0.35,color:COLORS_3D.HSQC} 
                });
              } catch (e) {
                debugWarn(`⚠️ Erreur setStyle HSQC pour atome ${i}:`, e.message);
              }
            }
          }
        }

        // ✅ Colorier les atomes H1
        for (const i of h1List) {
          if (typeof i === 'number' && i >= 0) {
            let stickColor="#808080", stickRadius=0.15;
            if (flags.HSQC && hsqcSet.has(i)) { 
              stickColor=COLORS_3D.HSQC; 
              stickRadius=0.2; 
            } else if (flags.COSY && cosySet.has(i)) { 
              stickColor=COLORS_3D.COSY; 
              stickRadius=0.25; 
            }
            try {
              model.setStyle({ index:i }, { 
                sphere:{radius:0.3,color:COLORS_3D.H1}, 
                stick:{radius:stickRadius,color:stickColor} 
              });
            } catch (e) {
              debugWarn(`⚠️ Erreur setStyle H1 pour atome ${i}:`, e.message);
            }
          }
        }
      }

      // Labels 3D centrés pour tous les hétéroatomes (≠ C et H)
      try {
        if (viewer.removeAllLabels) {
          viewer.removeAllLabels();
        }

        const atoms = (typeof model.selectedAtoms === "function")
          ? model.selectedAtoms({})
          : [];

        for (const a of atoms) {
          if (!a) continue;
          
          const sym = (a?.elem || "").toUpperCase();
          if (!sym || sym === "C" || sym === "H") continue;
          
          if (!Number.isFinite(a.x) || !Number.isFinite(a.y) || !Number.isFinite(a.z)) {
            continue;
          }

          viewer.addLabel(sym, {
            position: { x: a.x, y: a.y, z: a.z },
            alignment: "center",
            showBackground: false,
            fontSize: 16,
            fontFamily: "sans-serif",
            fontColor: "white",
            inFront: true,
          });
        }
      } catch (e) {
        debugWarn("⚠️ Erreur lors de l'ajout des labels 3D:", e?.message || e);
      }

      viewer.zoomTo(); 
      viewer.render();
      setErr(null);
    } catch (e) {
      const msg = e?.message || String(e);
      console.error("❌ Erreur 3D:", msg);
      setErr(msg);
    } finally {
      setLoading(false);
    }
  }, [multi, computeMultiHighlights, COLORS_3D, flags]);



  const load = useCallback(async () => {
    setErr(null);
    if (!(multi?.fromAll ? multi?.pivotSource : selectedRow?.Source) || !csvUrl) {
      setSvg(null); return;
    }
    const srcName = multi?.fromAll ? multi.pivotSource : selectedRow.Source;
    const resolved = await resolveMolForSource(csvUrl, srcName, multi?.taskId || taskId);
    if (!resolved || !resolved.molblock) { setErr("Molécule introuvable pour cette source."); setSvg(null); return; }
    await render2D(resolved);
    if (show3D) await render3D(resolved);
  }, [csvUrl, selectedRow, render2D, render3D, taskId, multi, show3D]);

  useEffect(() => { load(); }, [load]);

  useEffect(() => () => { try { viewerRef.current?.clear?.(); viewerRef.current=null; } catch {} }, []);

  if (compact) {
    return (
      <div style={{ position:"relative", width:"100%", height: molHeight, overflow:"hidden" }}>
        {!err && svg && (
          <div
            style={{ position:"absolute", inset:0, display:"flex", alignItems:"center", justifyContent:"center" }}
            dangerouslySetInnerHTML={{ __html: svg }}
          />
        )}
        {loading && (
          <div style={{ position:"absolute", inset:0, display:"flex", alignItems:"center", justifyContent:"center",
                        color:"#64748b", pointerEvents:"none", fontSize:12 }}>
            Chargement…
          </div>
        )}
        {err && (
          <div style={{ position:"absolute", inset:0, display:"flex", alignItems:"center", justifyContent:"center",
                        color:"#dc2626", background:"#fff", fontSize:12, padding:8, textAlign:"center" }}>
            {err}
          </div>
        )}
      </div>
    );
  }

  return (
    <div style={{ border:"1px solid #e5e7eb", borderRadius:8, background:"#fff", overflow:"hidden", padding:8 }}>
      <div style={{ display:"grid", gridTemplateColumns: show3D?"1fr 1fr":"1fr", gap:12, height: molHeight }}>
        {/* 2D */}
        <div style={{ border:"1px solid #e5e7eb", borderRadius:8, background:"#fff", overflow:"hidden", position:"relative" }}>
          <div style={{ fontWeight:600, padding:8, borderBottom:"1px solid #e5e7eb", fontSize:12 }}>2D</div>
          {!err && svg && (
            <div style={{ position:"absolute", top:40,left:0,right:0,bottom:0, display:"flex",alignItems:"center",justifyContent:"center" }}
                dangerouslySetInnerHTML={{ __html: svg }} />
          )}
          {loading && (
            <div style={{ position:"absolute", top:40,left:0,right:0,bottom:0, display:"flex",alignItems:"center",justifyContent:"center",
                          color:"#64748b", pointerEvents:"none" }}>
              Chargement 2D…
            </div>
          )}
        </div>
        {/* 3D */}
        {show3D && (
          <div style={{ border:"1px solid #e5e7eb", borderRadius:8, background:"#fff", overflow:"hidden", position:"relative" }}>
            <div style={{ fontWeight:600, padding:8, borderBottom:"1px solid #e5e7eb", fontSize:12 }}>3D</div>
            <div ref={div3DRef} style={{ position:"absolute", top:40,left:0,right:0,bottom:0 }} />
            {loading && (
              <div style={{ position:"absolute", top:40,left:0,right:0,bottom:0, display:"flex",alignItems:"center",justifyContent:"center",
                            color:"#64748b", pointerEvents:"none" }}>
                Chargement 3D…
              </div>
            )}
          </div>
        )}
      </div>

      {/* Barre de coches locale */}
      {showLocalToggles && (
        <div style={{ marginTop:8, display:"flex", gap:10, flexWrap:"wrap", alignItems:"center", fontSize:12 }}>
          {[
            ["H1", customColors.H1.hex],
            ["C13", customColors.C13.hex],
            ["COSY", customColors.COSY.hex],
            ["HSQC", customColors.HSQC.hex],
          ].map(([k,hex]) => (
            <label key={k}
              style={{ display:"inline-flex", alignItems:"center", gap:6, padding:"4px 8px",
                      border:"1px solid #e5e7eb", borderRadius:999, background: flags[k]?"#f0f7ff":"#fff", cursor:"pointer" }}>
              <input type="checkbox" checked={!!flags[k]} onChange={()=>toggle(k)} style={{ width:14,height:14 }}/>
              <span style={{ width:12,height:12,borderRadius:3,background:hex,border:"1px solid #cbd5e1" }}/>
              <span style={{ fontWeight:600 }}>
                {METHOD_LABELS_SHORT[k] || k}
              </span>
            </label>
          ))}
        </div>
      )}

      {/* ✅ AFFICHAGE D'ERREUR AMÉLIORÉ */}
      {err && (
        <div style={{ 
          padding:12, 
          marginTop:8,
          color:"#dc2626", 
          background:"#fee2e2", 
          borderRadius:8,
          border:"1px solid #fca5a5",
          fontSize:12 
        }}>
          <div style={{ fontWeight:700, marginBottom:6, fontSize:14, display:"flex", alignItems:"center", gap:6 }}>
            <span style={{ fontSize:20 }}>⚠️</span>
            {err.includes("WebGL") ? "WebGL indisponible" : "Erreur de rendu"}
          </div>
          <div style={{ whiteSpace:"pre-wrap", lineHeight:1.5 }}>{err}</div>
          {err.includes("WebGL") && (
            <div style={{ 
              marginTop:12, 
              padding:10,
              background:"#fff",
              borderRadius:6,
              fontSize:11, 
              color:"#991b1b",
              border:"1px solid #fca5a5"
            }}>
              <div style={{ fontWeight:700, marginBottom:6 }}>💡 Solutions rapides :</div>
              <ol style={{ margin:0, paddingLeft:20, lineHeight:1.8 }}>
                <li>
                  <strong>Firefox :</strong> Tapez <code style={{ background:"#fef2f2", padding:"2px 6px", borderRadius:3 }}>about:config</code> dans la barre d'adresse, 
                  puis réglez <code style={{ background:"#fef2f2", padding:"2px 6px", borderRadius:3 }}>webgl.disabled</code> sur <strong>false</strong>
                </li>
                <li>
                  <strong>Pilotes graphiques :</strong> Mettez à jour les pilotes de votre carte graphique (NVIDIA/AMD/Intel)
                </li>
                <li>
                  <strong>Navigateur alternatif :</strong> Essayez Chrome ou Edge si le problème persiste
                </li>
              </ol>
            </div>
          )}
        </div>
      )}
    </div>
  );
}


/* ===========================
   Métriques d'une ligne CSV
   =========================== */

const extractMetricsFromRow = (row) => {
  if (!row) return {};
  const metrics = {};
  if (row.coverage_ransy !== undefined && row.coverage_ransy !== "") {
    const v = Number(row.coverage_ransy || 0).toFixed(2) + "%";
    metrics.coverage_ransy = v;
  }
  if (row.coverage_sub !== undefined && row.coverage_sub !== "") {
    const v = Number(row.coverage_sub || 0).toFixed(2) + "%";
    metrics.coverage_sub = v;
  }
  if (row.avg_variabilite !== undefined && row.avg_variabilite !== "") {
    metrics.avg_variabilite = row.avg_variabilite;
  }
  return metrics;
};

/* ===========================
   🔎 MODALE — Toutes les molécules (grille visuelle)
   =========================== */

function ShowAllMoleculesModal({
  isOpen,
  onClose,
  sources = [],
  page = 1,
  setPage,
  pageSize = 20,
  setPageSize,
  getRowsByType,
  csvUrl,
  taskId,
  customColors,
}) {
  const stripMethodPrefix = (name) => {
    if (!name) return name;
    // Ex: "1H_MOL", "13C - MOL", "COSY: MOL", "HSQC MOL"
    return String(name).replace(/^(?:1H|13C|COSY|HSQC)[\s_\-:]+/i, "");
  };

  const totalPages = Math.max(1, Math.ceil(sources.length / pageSize));
  const clampedPage = Math.min(Math.max(1, page), totalPages);
  const start = (clampedPage - 1) * pageSize;
  const end = Math.min(sources.length, start + pageSize);
  const pageSources = sources.slice(start, end);

  if (!isOpen) return null;

  return (
    <div
      style={{ position:"fixed", inset:0, zIndex:10000, background:"rgba(0,0,0,0.55)",
               display:"flex", alignItems:"center", justifyContent:"center" }}
      onClick={onClose}
    >
      <div
        style={{ width:"min(1100px, 96vw)", maxHeight:"92vh", background:"#fff", borderRadius:12,
                 boxShadow:"0 24px 80px rgba(0,0,0,.35)", display:"flex", flexDirection:"column" }}
        onClick={(e)=>e.stopPropagation()}
      >
        {/* Header */}
        <div style={{ padding:12, borderBottom:"1px solid #e5e7eb",
                      display:"flex", alignItems:"center", justifyContent:"space-between" }}>
          <div style={{ fontWeight:800, fontSize:14 }}>
            Toutes les molécules (triées par moyenne RANSY)
          </div>
          <div style={{ display:"flex", gap:8, alignItems:"center" }}>
            <label style={{ fontSize:12, color:"#475569" }}>
              Par page:&nbsp;
              <select
                value={pageSize}
                onChange={(e)=>setPageSize(Math.max(10, +e.target.value||20))}
                style={{ padding:"4px 6px", borderRadius:6, border:"1px solid #cbd5e1" }}
              >
                {[20,30,50].map(n=><option key={n} value={n}>{n}</option>)}
              </select>
            </label>
            <button onClick={onClose}
              style={{ padding:"6px 10px", borderRadius:8, border:"1px solid #cbd5e1", background:"#fff", cursor:"pointer", fontWeight:700 }}>
              Fermer
            </button>
          </div>
        </div>

        {/* Pagination top */}
        <div style={{ padding:"8px 12px", display:"flex", gap:8, alignItems:"center", borderBottom:"1px solid #eef2f7" }}>
          <button onClick={()=>setPage(p=>Math.max(1,p-1))} disabled={clampedPage<=1}
                  style={{ padding:"4px 8px", borderRadius:6, border:"1px solid #cbd5e1", background:"#fff", cursor:"pointer" }}>◀</button>
          <div style={{ fontSize:12, color:"#334155" }}>
            Page {clampedPage} / {totalPages} &nbsp;•&nbsp; {sources.length} molécules
          </div>
          <button onClick={()=>setPage(p=>Math.min(totalPages,p+1))} disabled={clampedPage>=totalPages}
                  style={{ padding:"4px 8px", borderRadius:6, border:"1px solid #cbd5e1", background:"#fff", cursor:"pointer" }}>▶</button>
        </div>


        {/* Grille de molécules visuelles */}
        <div style={{ flex:1, overflow:"auto", padding:"12px" }}>
          <div style={{ 
            display:"grid", 
            gridTemplateColumns:"repeat(auto-fill, minmax(180px, 1fr))", 
            gap:12 
          }}>
            {pageSources.map((source, idx) => {
              const globalIdx = start + idx + 1;
              const moleculeName = stripMethodPrefix(source.replace(/_predicted$/i, ""));
              
              return (
                <div key={source} style={{ 
                  border:"1px solid #e5e7eb", 
                  borderRadius:8, 
                  padding:8, 
                  background:"#fff",
                  boxShadow:"0 1px 3px rgba(0,0,0,0.1)",
                  transition:"transform 0.2s ease, box-shadow 0.2s ease",
                  transformOrigin:"center",
                  position:"relative",
                  cursor:"pointer"
                }} 
                onMouseEnter={(e) => {
                  // Important: use currentTarget (the card), not the hovered child
                  e.currentTarget.style.transform = "scale(1.02)";
                  e.currentTarget.style.boxShadow = "0 4px 12px rgba(0,0,0,0.15)";
                  e.currentTarget.style.zIndex = "2";
                }}
                onMouseLeave={(e) => {
                  e.currentTarget.style.transform = "scale(1)";
                  e.currentTarget.style.boxShadow = "0 1px 3px rgba(0,0,0,0.1)";
                  e.currentTarget.style.zIndex = "0";
                }}>
                  
                  {/* Vraie molécule avec MoleculePane */}
                  <div style={{ height:140, marginBottom:8 }}>
                    <MoleculePane
                      csvUrl={csvUrl}
                      selectedRow={{ Source: source }}
                      show3D={false}
                      taskId={taskId}
                      customColors={customColors}
                      molHeight={140}
                      showLocalToggles={false}
                      compact={true}
                      style={{ height:"100%" }}
                      multi={{
                        fromAll: true,
                        pivotSource: source,
                        taskId,
                        rowsByType: getRowsByType(source),
                      }}
                    />
                  </div>
                  
                  {/* Nom de la molécule */}
                  <div style={{ 
                    fontWeight:600, 
                    fontSize:11, 
                    marginBottom:6, 
                    wordBreak:"break-word",
                    textAlign:"center"
                  }}>
                    #{globalIdx} {moleculeName}
                  </div>
                </div>
              );
            })}
          </div>
        </div>

        {/* Pagination bottom */}
        <div style={{ padding:"8px 12px", borderTop:"1px solid #e5e7eb",
                      display:"flex", gap:8, alignItems:"center", justifyContent:"center" }}>
          <div style={{ display:"flex", gap:8, alignItems:"center" }}>
            <button onClick={()=>setPage(p=>Math.max(1,p-1))} disabled={clampedPage<=1}
                    style={{ padding:"6px 10px", borderRadius:8, border:"1px solid #cbd5e1", background:"#fff", cursor:"pointer" }}>◀ Précédent</button>
            <div style={{ fontSize:12, color:"#334155" }}>Page {clampedPage} / {totalPages}</div>
            <button onClick={()=>setPage(p=>Math.min(totalPages,p+1))} disabled={clampedPage>=totalPages}
                    style={{ padding:"6px 10px", borderRadius:8, border:"1px solid #cbd5e1", background:"#fff", cursor:"pointer" }}>Suivant ▶</button>
          </div>
        </div>
      </div>
    </div>
  );
}



/* ===========================
   Viewer principal exporté avec TIMEOUTS DE SÉCURITÉ
   =========================== */

export default function Viewer({
  csvs,
  taskId = null,
  defaultShow3D = true,
  annotationResults = {},
  pdataDirs = {},
  bgOptions = {},
  customColors: propsCustomColors,
  onColorChange,
  // NEW: Props pour contrôle externe de la source (style SIRIUS)
  externalSource = null,
  hideSidebar = false,
}) {
  // Utiliser les couleurs passées en props ou les valeurs par défaut
  const customColors = propsCustomColors || DEFAULT_COLORS;

  const setCustomColor = useCallback((method, color) => {
    if (onColorChange) {
      onColorChange(method, color);
    }
  }, [onColorChange]);

  const [, setRowsByType] = useState({});
  const [mapsByType, setMapsByType] = useState({});
  const [loading, setLoading] = useState(true);
  const [err, setErr] = useState(null);

  // Si externalSource est fourni, l'utiliser ; sinon garder l'état interne
  const [internalSource, setInternalSource] = useState(null);
  const selectedSource = externalSource || internalSource;
  const setSelectedSource = useCallback(
    (next) => {
      if (externalSource) return;
      setInternalSource(next);
    },
    [externalSource]
  );
  
  const [query, setQuery] = useState("");

  const [showTheory, setShowTheory] = useState(true);
  const [showExp, setShowExp] = useState(true);
  const [reverse1D, setReverse1D] = useState(true);
  const [show3D, setShow3D] = useState(defaultShow3D);

  const [zoomRectEnabled, setZoomRectEnabled] = useState(true);
  const [colorPickerOpen, setColorPickerOpen] = useState(false);
  const [relevel2DEnabled, setRelevel2DEnabled] = useState(true); // Option pour désactiver le relevel au scroll

  // NEW: nonce par type pour forcer un remount des canvases
  const [refreshKeyByType, setRefreshKeyByType] = useState({});
  const bumpKey = useCallback((type) => {
    setRefreshKeyByType(k => ({ ...k, [type]: (k[type] ?? 0) + 1 }));
  }, []);

  const wheelRequestRef = useRef({});
  
  const presentTypes = useMemo(
    () => Object.keys(csvs || {}).filter(Boolean),
    [csvs]
  );

  const hasAnyData = useMemo(() => presentTypes.length > 0, [presentTypes]);

  const [bgByType, setBgByType] = useState(
    () => annotationResults || {}
  );
  const [busyByType, setBusyByType] = useState({});

  // État local pour les paramètres de résolution et niveau
  const [resByType, setResByType] = useState({
    "1H":  { width_px: 1200, height_px: 500, dpi: 150 },
    "13C": { width_px: 1200, height_px: 500, dpi: 150 },
    COSY:  { width_px: 1000, height_px: 800, dpi: 150 },
    HSQC:  { width_px: 1000, height_px: 800, dpi: 150 },
  });

  const [levelByType, setLevelByType] = useState({
    COSY: 0.005,
    HSQC: 0.005
  });

  // 🔁 AVANT : const originalBgRef = useRef({});
  // ✅ APRÈS : structure par type ET par source : { [type]: { [source]: snap } }
  const originalBgRef = useRef({});

  useEffect(() => {
    const snap = {};
    for (const k of ["1H", "13C", "COSY", "HSQC"]) {
      if (annotationResults?.[k]) snap[k] = annotationResults[k];
    }
    setBgByType(prev => ({ ...prev, ...snap }));
  }, [annotationResults]);

  // Quand des résultats initiaux arrivent, on snap selon la source courante
  useEffect(() => {
    if (!selectedSource) return;
    for (const k of ["1H", "13C", "COSY", "HSQC"]) {
      if (annotationResults?.[k]) {
        originalBgRef.current[k] = originalBgRef.current[k] || {};
        if (!originalBgRef.current[k][selectedSource]) {
          originalBgRef.current[k][selectedSource] = annotationResults[k];
        }
      }
    }
  }, [annotationResults, selectedSource]);

  // Quand la source change, si on a un fond en mémoire, on le "snap" aussi
  useEffect(() => {
    if (!selectedSource) return;
    for (const k of ["1H", "13C", "COSY", "HSQC"]) {
      const bg = bgByType?.[k];
      if (bg) {
        originalBgRef.current[k] = originalBgRef.current[k] || {};
        if (!originalBgRef.current[k][selectedSource]) {
          originalBgRef.current[k][selectedSource] = bg;
        }
      }
    }
  }, [selectedSource, bgByType]);

  const loadAll = useCallback(async () => {
    if (!csvs) return;
    setLoading(true);
    setErr(null);
    try {
      const entries = await Promise.all(
        presentTypes.map(async (t) => {
          try {
            const r = await fetchWithRetry(csvs[t], 3, 500);
            const text = await r.text();
            const { data } = parseCSVSmart(text);
            const map = new Map();
            for (const row of data || []) {
              if (row?.Source) map.set(row.Source, row);
            }
            return [t, data || [], map, null];
          } catch (error) {
            return [t, [], new Map(), error.message];
          }
        })
      );

      const rb = {};
      const mb = {};
      const errors = {};
      for (const [t, rows, map, error] of entries) {
        rb[t] = rows;
        mb[t] = map;
        if (error) errors[t] = error;
      }

      setRowsByType(rb);
      setMapsByType(mb);

      const failed = Object.keys(errors);
      if (failed.length > 0) {
        setErr(`Spectres non disponibles: ${failed.join(", ")}`);
      }

      const allSets = Object.values(mb).map(m => new Set(m.keys()));
      let common = [];
      if (allSets.length) {
        common = [...allSets[0]].filter(s =>
          allSets.every(S => S.has(s))
        );
        if (common.length === 0) {
          const union = new Set();
          Object.values(mb).forEach(m => {
            for (const s of m.keys()) union.add(s);
          });
          common = [...union];
        }
      }
      common.sort();
      setSelectedSource(common[0] || null);
    } catch (e) {
      setErr(e.message || "Erreur de chargement des CSV");
    } finally {
      setLoading(false);
    }
  }, [csvs, presentTypes, setSelectedSource]);

  useEffect(() => { loadAll(); }, [loadAll]);

  const allSources = useMemo(() => {
    const union = new Set();
    Object.values(mapsByType).forEach(m => {
      for (const s of m.keys()) union.add(s);
    });
    let arr = [...union];

    const parseCoverage = (row) => {
      if (!row) return null;
      // Accept a few possible field names (CSV variations)
      const candidates = [
        row.coverage_ransy,
        row.coverage_ransy_mean,
        row.coverage_ransy_avg,
        row.coverage_ransy_percent,
      ];
      for (const c of candidates) {
        if (c === undefined || c === null || c === "") continue;
        const v = parseFloat(String(c).replace(",", "."));
        if (Number.isFinite(v)) return v;
      }
      return null;
    };

    const withCov = arr.map(source => {
      const coverages = [];
      for (const type of ["1H", "13C", "COSY", "HSQC"]) {
        const row = mapsByType[type]?.get(source);
        const cov = parseCoverage(row);
        if (cov != null) coverages.push(cov);
      }
      const avg = coverages.length > 0
        ? coverages.reduce((s, v) => s + v, 0) / coverages.length
        : -1;
      return { source, avgCoverage: avg };
    });

    withCov.sort((a, b) => (b.avgCoverage - a.avgCoverage) || a.source.localeCompare(b.source));
    arr = withCov.map(x => x.source);

    if (!query) return arr;
    const q = query.toLowerCase();
    return arr.filter(s => s.toLowerCase().includes(q));
  }, [mapsByType, query]);

  // 🔥 Liste globale triée (sans filtre de recherche) pour la modale
  const sortedSourcesAll = useMemo(() => {
    const union = new Set();
    Object.values(mapsByType).forEach(m => {
      for (const s of m.keys()) union.add(s);
    });
    const arr = [...union];

    const parseCoverage = (row) => {
      if (!row) return null;
      const candidates = [
        row.coverage_ransy,
        row.coverage_ransy_mean,
        row.coverage_ransy_avg,
        row.coverage_ransy_percent,
      ];
      for (const c of candidates) {
        if (c === undefined || c === null || c === "") continue;
        const v = parseFloat(String(c).replace(",", "."));
        if (Number.isFinite(v)) return v;
      }
      return null;
    };

    const withCov = arr.map(source => {
      const covs = [];
      for (const type of ["1H", "13C", "COSY", "HSQC"]) {
        const row = mapsByType[type]?.get(source);
        const v = parseCoverage(row);
        if (v != null) covs.push(v);
      }
      const avg = covs.length ? (covs.reduce((s,v)=>s+v,0)/covs.length) : -1;
      return { source, avgCoverage: avg };
    });

    withCov.sort((a,b) => (b.avgCoverage - a.avgCoverage) || a.source.localeCompare(b.source));
    return withCov.map(x => x.source);
  }, [mapsByType]);

  // Helper lignes par méthode pour une source (modale)
  const rowsByTypeForSource = useCallback((source) => ({
    "1H":   mapsByType["1H"]?.get(source),
    "13C":  mapsByType["13C"]?.get(source),
    "COSY": mapsByType["COSY"]?.get(source),
    "HSQC": mapsByType["HSQC"]?.get(source),
  }), [mapsByType]);

  const rowFor = useCallback(
    (t) => mapsByType[t]?.get(selectedSource) || null,
    [mapsByType, selectedSource]
  );

  const oneH = useMemo(() => {
    const r = rowFor("1H");
    return r ? {
      theory: (r.ppm_molecule || "")
        .split("|")
        .map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      expX: (r.RANSY_signals || "")
        .split("|")
        .map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      expI: (r.RANSY_mean_intensities || r.intensity_ratios || "")
        .split(";")
        .map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      matches: (r.ppm_molecule_matched || "")
        .split("|")
        .map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      metrics: extractMetricsFromRow(r),
    } : null;
  }, [rowFor]);

  const c13 = useMemo(() => {
    const r = rowFor("13C");
    return r ? {
      theory: (r.ppm_molecule || "")
        .split("|").map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      expX: (r.RANSY_signals || "")
        .split("|").map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      expI: (r.RANSY_mean_intensities || r.intensity_ratios || "")
        .split(";").map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      matches: (r.ppm_molecule_matched || "")
        .split("|").map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      metrics: extractMetricsFromRow(r),
    } : null;
  }, [rowFor]);

  const cosy = useMemo(() => {
    const r = rowFor("COSY");
    return r ? {
      theory: (r.ppm_molecule || "")
        .split("||")
        .map(p => p.split("|")
        .map(v => parseFloat(String(v).replace(",", "."))))
        .filter(a => a.length === 2 && a.every(Number.isFinite)),
      expPts: (r.RANSY_signals || "")
        .split("||")
        .map(p => p.split("|")
        .map(v => parseFloat(String(v).replace(",", "."))))
        .filter(a => a.length === 2 && a.every(Number.isFinite)),
      expI: (r.RANSY_mean_intensities || r.intensity_ratios || "")
        .split(";")
        .map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      pairs: (r.ppm_molecule_matched || "")
        .split("||")
        .map(p => p.split("|")
        .map(v => parseFloat(String(v).replace(",", "."))))
        .filter(a => a.length === 2 && a.every(Number.isFinite)),
      metrics: extractMetricsFromRow(r),
    } : null;
  }, [rowFor]);

  const hsqc = useMemo(() => {
    const r = rowFor("HSQC");
    return r ? {
      theory: (r.ppm_molecule || "")
        .split("||")
        .map(p => p.split("|")
        .map(v => parseFloat(String(v).replace(",", "."))))
        .filter(a => a.length === 2 && a.every(Number.isFinite)),
      expPts: (r.RANSY_signals || "")
        .split("||")
        .map(p => p.split("|")
        .map(v => parseFloat(String(v).replace(",", "."))))
        .filter(a => a.length === 2 && a.every(Number.isFinite)),
      expI: (r.RANSY_mean_intensities || r.intensity_ratios || "")
        .split(";")
        .map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      pairs: (
        r.ppm_molecule_matched_unfolded ||
        r.ppm_molecule_matched || ""
      )
        .split("||")
        .map(p => p.split("|")
        .map(v => parseFloat(String(v).replace(",", "."))))
        .filter(a => a.length === 2 && a.every(Number.isFinite)),
      metrics: extractMetricsFromRow(r),
    } : null;
  }, [rowFor]);

  const pivotCsvUrl = useMemo(
    () => csvs["1H"] || csvs["13C"] || csvs["COSY"] || csvs["HSQC"] || null,
    [csvs]
  );

  const getMeta = (obj) => obj?.background?.meta || obj?.meta || null;

  // eslint-disable-next-line react-hooks/exhaustive-deps
  const setBusy = useCallback((type, v) => setBusyByType(prev => ({ ...prev, [type]: !!v })));
  
  // ✅ storeResult corrige : ajoute un nonce anti-cache et met à jour le snapshot (type, source)
  // eslint-disable-next-line react-hooks/exhaustive-deps
  const storeResult = useCallback((type, res) => {
    if (!res) return;

    // nonce anti-cache
    const nonce = `${Date.now()}-${Math.random().toString(36).slice(2)}`;
    res.meta = { ...(res.meta || {}), nonce };

    // met à jour l'état courant (affichage)
    setBgByType(prev => ({ ...prev, [type]: res }));

    // si 2D, on synchronise le lvl UI
    if ((type === "COSY" || type === "HSQC") &&
        Number.isFinite(+res?.meta?.level_start)) {
      setLevelByType(p => ({ ...p, [type]: +res.meta.level_start }));
    }

    // 🔒 snapshot initial par (type, source) : ne l'écrit que s'il est absent
    if (selectedSource) {
      originalBgRef.current[type] = originalBgRef.current[type] || {};
      if (!originalBgRef.current[type][selectedSource]) {
        originalBgRef.current[type][selectedSource] = res;
      }
    }
  });

  const extractWindow1D = (meta) => {
    if (!meta) return null;
    if (meta.ppm_min != null && meta.ppm_max != null) {
      return { ppmMin: +meta.ppm_min, ppmMax: +meta.ppm_max };
    }
    if (meta.ppm_left != null && meta.ppm_right != null) {
      const left = +meta.ppm_left;
      const right = +meta.ppm_right;
      return {
        ppmMin: Math.min(left, right),
        ppmMax: Math.max(left, right),
      };
    }
    return null;
  };

  const extractWindow2D = (meta) => {
    if (!meta) return null;
    let f2Min, f2Max, f1Min, f1Max;
    if (meta.f2_ppm_min != null && meta.f2_ppm_max != null) {
      f2Min = +meta.f2_ppm_min; f2Max = +meta.f2_ppm_max;
    } else if (meta.f2_left != null && meta.f2_right != null) {
      const L = +meta.f2_left; const R = +meta.f2_right;
      f2Min = Math.min(L, R); f2Max = Math.max(L, R);
    }
    if (meta.f1_ppm_min != null && meta.f1_ppm_max != null) {
      f1Min = +meta.f1_ppm_min; f1Max = +meta.f1_ppm_max;
    } else if (meta.f1_top != null && meta.f1_bottom != null) {
      const T = +meta.f1_top; const B = +meta.f1_bottom;
      f1Min = Math.min(T, B); f1Max = Math.max(T, B);
    }
    if ([f2Min,f2Max,f1Min,f1Max].every(v => Number.isFinite(v))) {
      return { f2Min, f2Max, f1Min, f1Max };
    }
    return null;
  };

  
  // CORRECTION : Ajout de timeout de sécurité pour applySettings
  const applySettings = useCallback(async (type) => {
    const pd = pdataDirs?.[type];
    if (!pd) return;
    
    try {
      setBusy(type, true);
      
      // Timeout de sécurité
      const timeoutId = setTimeout(() => {
        console.warn(`Timeout pour applySettings ${type}, déblocage forcé`);
        setBusy(type, false);
      }, 30000);
      
      const size = resByType[type] || {};
      const meta = getMeta(bgByType?.[type]);
      let res;
      
      if (type === "1H" || type === "13C") {
        const win = extractWindow1D(meta) || null;
        const line = { width: 1.5 };
        const fn = type === "1H" ? interactive1H : interactive13C;
        if (win) {
          const payload = {
            ppm_min: Math.min(win.ppmMin, win.ppmMax),
            ppm_max: Math.max(win.ppmMin, win.ppmMax),
            ppmMin: Math.min(win.ppmMin, win.ppmMax),
            ppmMax: Math.max(win.ppmMin, win.ppmMax),
            size, line, taskId,
          };
          res = await fn(pd, payload);
        } else {
          res = await fn(pd, { size, line, taskId });
        }
      } else {
        const win = extractWindow2D(meta) || {};
        const contour = {
          level_start: levelByType[type] ?? 0.005,
          n_levels: 20,
          transparent: true,
        };
        const fn = type === "COSY" ? interactiveCOSY : interactiveHSQC;
        res = await fn(pd, { ...win, size, contour, taskId });
      }
      
      clearTimeout(timeoutId);
      storeResult(type, res);
    } catch (e) {
      console.error(`applySettings ${type}:`, e);
    } finally {
      setBusy(type, false);
    }
  }, [pdataDirs, setBusy, resByType, bgByType, storeResult, taskId, levelByType]);

  // CORRECTION : Ajout de timeout de sécurité pour doZoom1DPixels
  const doZoom1DPixels = useCallback(async (type, rectPxDOM, display) => {
    const pd = pdataDirs?.[type];
    const meta = getMeta(bgByType?.[type]);
    if (!pd || !meta || !rectPxDOM || !(rectPxDOM.w > 1)) return;

    try {
      setBusy(type, true);
      
      // Timeout de sécurité
      const timeoutId = setTimeout(() => {
        console.warn(`Timeout pour zoom1DPixels ${type}, déblocage forcé`);
        setBusy(type, false);
      }, 30000);
      
      const size = resByType[type] || {};
      const line = { width: 1.5 };
      const { displayWidth, displayHeight } = display || {};
      const res = await zoom1DFromPixels(
        type, pd, meta, rectPxDOM,
        { taskId, size, line, displayWidth, displayHeight }
      );
      
      clearTimeout(timeoutId);
      storeResult(type, res);
    } catch (e) {
      console.error(`ZoomPixels ${type} error:`, e);
    } finally {
      setBusy(type, false);
    }
  }, [pdataDirs, bgByType, setBusy, resByType, taskId, storeResult]);

  // CORRECTION : Ajout de timeout de sécurité pour doZoom1D
  const doZoom1D = useCallback(async (type, ppmMin, ppmMax) => {
    const pd = pdataDirs?.[type];

    if (!pd) {
      setBgByType(prev => ({
        ...prev,
        [type]: {
          ...prev[type],
          meta: {
            ...prev[type]?.meta,
            ppm_min: Math.min(ppmMin, ppmMax),
            ppm_max: Math.max(ppmMin, ppmMax),
            ppm_left: Math.max(ppmMin, ppmMax),
            ppm_right: Math.min(ppmMin, ppmMax),
            nonce: `${Date.now()}`
          },
        },
      }));
      bumpKey(type);
      return;
    }

    try {
      setBusy(type, true);
      
      // Timeout de sécurité
      const timeoutId = setTimeout(() => {
        console.warn(`Timeout pour zoom1D ${type}, déblocage forcé`);
        setBusy(type, false);
      }, 30000);
      
      const size = resByType[type] || {};
      const line = { width: 1.5 };
      const fn = type === "1H" ? interactive1H : interactive13C;
      const payload = {
        ppm_min: Math.min(ppmMin, ppmMax),
        ppm_max: Math.max(ppmMin, ppmMax),
        ppmMin: Math.min(ppmMin, ppmMax),
        ppmMax: Math.max(ppmMin, ppmMax),
        size, line, taskId,
      };
      const res = await fn(pd, payload);
      
      clearTimeout(timeoutId);
      storeResult(type, res);
    } catch (e) {
      console.error(`Zoom ${type} error:`, e);
    } finally {
      setBusy(type, false);
    }
  }, [pdataDirs, bumpKey, setBusy, resByType, taskId, storeResult]);

  // *** CORRIGÉ *** : doReset1D préfère recalculer au backend si possible, sinon remet le snapshot
  const doReset1D = useCallback(async (type) => {
    const pd = pdataDirs?.[type];
    const snap = originalBgRef.current?.[type]?.[selectedSource];

    // ✅ Récupérer les bornes DÉFINIES PAR L'UTILISATEUR
    const userWindow = bgOptions?.[type] || {};
    const hasUserWindow = 
      (type === "1H" || type === "13C") &&
      userWindow.ppm_min != null && 
      userWindow.ppm_max != null;

    // Si interaction backend dispo → on redemande l'image avec les bornes définies
    if (pd) {
      try {
        setBusy(type, true);
        const timeoutId = setTimeout(() => {
          console.warn(`Timeout pour reset1D ${type}, déblocage forcé`);
          setBusy(type, false);
        }, 30000);

        const size = resByType[type] || {};
        const line = { width: 1.5 };
        const fn = type === "1H" ? interactive1H : interactive13C;
        
        // ✅ UTILISER LES BORNES UTILISATEUR si définies
        const payload = { size, line, taskId };
        
        if (hasUserWindow) {
          // Bornes utilisateur
          payload.ppm_min = Math.min(userWindow.ppm_min, userWindow.ppm_max);
          payload.ppm_max = Math.max(userWindow.ppm_min, userWindow.ppm_max);
          payload.ppmMin = payload.ppm_min;
          payload.ppmMax = payload.ppm_max;
        }
        // Sinon, pas de bornes = plein champ natif du backend

        const res = await fn(pd, payload);

        clearTimeout(timeoutId);
        storeResult(type, res);
        bumpKey(type); // force remount canvas
        return;
      } catch (e) {
        console.error(`Reset ${type} via backend error:`, e);
      } finally {
        setBusy(type, false);
      }
    }

    // Sinon, fallback à l'état snapshot local
    if (snap) {
      const restored = {
        ...snap,
        meta: { ...snap.meta, nonce: `${Date.now()}` },
      };
      setBgByType(prev => ({ ...prev, [type]: restored }));
      originalBgRef.current[type] = originalBgRef.current[type] || {};
      originalBgRef.current[type][selectedSource] = restored;
      bumpKey(type); // force remount
    }
  }, [pdataDirs, selectedSource, bgOptions, setBusy, resByType, taskId, storeResult, bumpKey]);

  // CORRECTION : Ajout de timeout de sécurité pour doZoom2D
  const doZoom2D = useCallback(async (type, { f2_min, f2_max, f1_min, f1_max }) => {
    const pd = pdataDirs?.[type];
    if (!pd) return;
    
    try {
      setBusy(type, true);
      
      // Timeout de sécurité
      const timeoutId = setTimeout(() => {
        console.warn(`Timeout pour zoom2D ${type}, déblocage forcé`);
        setBusy(type, false);
      }, 30000);
      
      const size = resByType[type] || {};
      const contour = {
        level_start: levelByType[type] ?? 0.005,
        n_levels: 20,
        transparent: true,
      };
      const fn = type === "COSY" ? interactiveCOSY : interactiveHSQC;
      const res = await fn(
        pd,
        {
          f2Min: f2_min, f2Max: f2_max,
          f1Min: f1_min, f1Max: f1_max,
          size, contour, taskId,
        }
      );
      
      clearTimeout(timeoutId);
      storeResult(type, res);
    } catch (e) {
      console.error(`Zoom ${type} error:`, e);
    } finally {
      setBusy(type, false);
    }
  }, [pdataDirs, setBusy, resByType, levelByType, taskId, storeResult]);


  const doWheel2D = useCallback(async (type, deltaY) => {
    // Si le relevel 2D est désactivé, ne rien faire
    if (!relevel2DEnabled) return;
    
    const pd = pdataDirs?.[type];
    const meta = getMeta(bgByType?.[type]);
    if (!pd || !meta) return;
    
    if (wheelRequestRef.current[type]) {
      try {
        wheelRequestRef.current[type].cancel();
      } catch {}
    }
    
    try {
      setBusy(type, true);
      
      const timeoutId = setTimeout(() => {
        console.warn(`Timeout wheel2D ${type}, annulation`);
        wheelRequestRef.current[type]?.cancel();
        setBusy(type, false);
      }, 10000);
      
      // ✅ CORRECTION : toujours envoyer size
      const request = relevelByWheelCancelable(
        type, pd, meta, deltaY,
        { 
          nLevels: 20, 
          transparent: true, 
          taskId,
          size: resByType[type] || {} // ← AJOUT
        }
      );
      
      wheelRequestRef.current[type] = request;
      const res = await request;
      
      clearTimeout(timeoutId);
      storeResult(type, res);
      wheelRequestRef.current[type] = null;
    } catch (e) {
      if (e?.name === 'AbortError') {
        console.log(`Requête wheel2D ${type} annulée`);
        return;
      }
      
      console.error(`Wheel relevel ${type}:`, e);
      
      // Fallback local (garder tel quel)
      if (e?.response?.status === 500 || e?.code === "ERR_BAD_RESPONSE") {
        console.warn(`Backend défaillant, ajustement local du niveau pour ${type}`);
        
        const currentLevel = meta?.level_start ?? levelByType[type] ?? 0.005;
        const direction = deltaY > 0 ? 1.2 : 0.8;
        const newLevel = Math.max(1e-6, currentLevel * direction);
        
        setLevelByType(p => ({ ...p, [type]: newLevel }));
        
        if (bgByType[type]?.meta) {
          setBgByType(prev => ({
            ...prev,
            [type]: {
              ...prev[type],
              meta: {
                ...prev[type].meta,
                level_start: newLevel,
                nonce: `${Date.now()}`
              },
            },
          }));
          bumpKey(type);
        }
      }
    } finally {
      setBusy(type, false);
    }
  }, [pdataDirs, bgByType, setBusy, taskId, storeResult, levelByType, bumpKey, resByType, relevel2DEnabled]);
  useEffect(() => {
    return () => {
      Object.values(wheelRequestRef.current).forEach(request => {
        try {
          request?.cancel();
        } catch {}
      });
      wheelRequestRef.current = {};
    };
  }, []);
    
  // Reset 2D fiable : conserve le fond actuel pendant l'attente backend
  // Reset 2D instantané + recalcul conditionnel et léger en arrière-plan
  const doReset2D = useCallback(async (type) => {
    const snap = originalBgRef.current?.[type]?.[selectedSource];
    const pd   = pdataDirs?.[type];

    // ✅ Récupérer les bornes DÉFINIES PAR L'UTILISATEUR
    const userWindow = bgOptions?.[type] || {};
    const hasUserWindow = 
      (type === "COSY" || type === "HSQC") &&
      [userWindow.f1_ppm_min, userWindow.f1_ppm_max,
      userWindow.f2_ppm_min, userWindow.f2_ppm_max].every(v => v != null);

    // Niveau de départ
    const currentLevelStart =
      (snap?.meta?.level_start) ??
      bgByType[type]?.meta?.level_start ??
      levelByType[type] ?? 0.005;

    // Affichage instantané du snapshot
    if (snap) {
      const restored = {
        ...snap,
        meta: {
          ...snap.meta,
          level_start: currentLevelStart,
          nonce: `${Date.now()}`,
        },
      };
      setBgByType(prev => ({ ...prev, [type]: restored }));
      bumpKey(type);
    }

    // Si pas de backend, on s'arrête
    if (!pd) return;

    setBusy(type, true);

    try {
      const timeoutId = setTimeout(() => {
        console.warn(`Timeout pour reset2D ${type}, déblocage forcé`);
        setBusy(type, false);
      }, 30000);

      const lightSize = {
        width_px:  Math.min(900,  (resByType[type]?.width_px  ?? 1000)),
        height_px: Math.min(700,  (resByType[type]?.height_px ?? 800)),
        dpi:       Math.min(120,  (resByType[type]?.dpi       ?? 150)),
      };
      const lightContour = {
        level_start: currentLevelStart,
        n_levels: 12,
        transparent: true,
      };

      const fn = type === "COSY" ? interactiveCOSY : interactiveHSQC;

      // ✅ UTILISER LES BORNES UTILISATEUR au lieu de "full: true"
      const payload = {
        size: lightSize,
        contour: lightContour,
        taskId,
        with_bbox: true,
        return_bbox: true,
      };

      // ✅ Si l'utilisateur a défini des bornes, les utiliser
      if (hasUserWindow) {
        payload.f1Min = userWindow.f1_ppm_min;
        payload.f1Max = userWindow.f1_ppm_max;
        payload.f2Min = userWindow.f2_ppm_min;
        payload.f2Max = userWindow.f2_ppm_max;
      } else {
        // Sinon, plein champ
        payload.reset = true;
        payload.full = true;
      }

      const resLight = await fn(pd, payload);

      clearTimeout(timeoutId);

      // Normaliser bbox
      if (resLight?.meta && !resLight.meta.plot_bbox_px) {
        const bb =
          resLight.meta.axis_bbox_px ||
          resLight.meta.axisRect?.pixel ||
          resLight.meta.axisRect?.px || null;
        if (bb && [bb.x, bb.y, bb.w, bb.h].every(Number.isFinite)) {
          resLight.meta.plot_bbox_px = { x: bb.x, y: bb.y, w: bb.w, h: bb.h };
        }
      }

      storeResult(type, resLight);
      bumpKey(type);

      // Qualité HD en arrière-plan si nécessaire
      const targetSize = resByType[type] || {};
      const wantQuality =
        // eslint-disable-next-line no-mixed-operators
        targetSize.width_px  && targetSize.width_px  > lightSize.width_px ||
        // eslint-disable-next-line no-mixed-operators
        targetSize.height_px && targetSize.height_px > lightSize.height_px ||
        (resByType[type]?.dpi ?? 150) > lightSize.dpi;

      if (wantQuality) {
        const qualityContour = {
          level_start: currentLevelStart,
          n_levels: 20,
          transparent: true,
        };
        
        const payloadHD = {
          size: resByType[type] || {},
          contour: qualityContour,
          taskId,
          with_bbox: true,
          return_bbox: true,
        };

        // ✅ Même logique pour la qualité HD
        if (hasUserWindow) {
          payloadHD.f1Min = userWindow.f1_ppm_min;
          payloadHD.f1Max = userWindow.f1_ppm_max;
          payloadHD.f2Min = userWindow.f2_ppm_min;
          payloadHD.f2Max = userWindow.f2_ppm_max;
        } else {
          payloadHD.reset = true;
          payloadHD.full = true;
        }

        fn(pd, payloadHD)
          .then((resHi) => {
            if (resHi?.meta && !resHi.meta.plot_bbox_px) {
              const bb =
                resHi.meta.axis_bbox_px ||
                resHi.meta.axisRect?.pixel ||
                resHi.meta.axisRect?.px || null;
              if (bb && [bb.x, bb.y, bb.w, bb.h].every(Number.isFinite)) {
                resHi.meta.plot_bbox_px = { x: bb.x, y: bb.y, w: bb.w, h: bb.h };
              }
            }
            storeResult(type, resHi);
            bumpKey(type);
          })
          .catch((e) => {
            console.warn(`Qualité ${type} échouée:`, e?.message || e);
          })
          .finally(() => setBusy(type, false));
      } else {
        setBusy(type, false);
      }
    } catch (e) {
      console.error(`Reset ${type} error:`, e);
      setBusy(type, false);
    }
  }, [
    pdataDirs,
    resByType,
    taskId,
    levelByType,
    bgByType,
    bgOptions,  // ✅ AJOUT dépendance
    bumpKey,
    selectedSource,
    storeResult,
    setBusy,
  ]);

  const canInteract1H   = !!pdataDirs?.["1H"];
  const canInteract13C  = !!pdataDirs?.["13C"];
  const canInteractCOSY = !!pdataDirs?.["COSY"];
  const canInteractHSQC = !!pdataDirs?.["HSQC"];

  // État pour la modale "Afficher toutes les molécules"
  const [showAllOpen, setShowAllOpen] = useState(false);
  const [allPage, setAllPage] = useState(1);
  const [allPageSize, setAllPageSize] = useState(20);

  /* petits panneaux de réglage */
  const renderControls1D = (type) => (
    <div style={{
      display: "flex", gap: 6, alignItems: "center",
      flexWrap: "wrap", marginBottom: 6, fontSize: 11,
    }}>
      <label style={{ fontSize: 10 }}>
        W:
        <input
          type="number" min={300}
          value={resByType[type].width_px}
          onChange={(e) =>
            setResByType(p => ({
              ...p,
              [type]: {
                ...p[type],
                width_px: Math.max(300, +e.target.value || 300),
              },
            }))
          }
          style={{ width: 70, fontSize: 10, padding: "1px 3px" }}
        />
      </label>
      <label style={{ fontSize: 10 }}>
        H:
        <input
          type="number" min={200}
          value={resByType[type].height_px}
          onChange={(e) =>
            setResByType(p => ({
              ...p,
              [type]: {
                ...p[type],
                height_px: Math.max(200, +e.target.value || 200),
              },
            }))
          }
          style={{ width: 70, fontSize: 10, padding: "1px 3px" }}
        />
      </label>
      <label style={{ fontSize: 10 }}>
        DPI:
        <input
          type="number" min={72}
          value={resByType[type].dpi}
          onChange={(e) =>
            setResByType(p => ({
              ...p,
              [type]: {
                ...p[type],
                dpi: Math.max(72, +e.target.value || 72),
              },
            }))
          }
          style={{ width: 60, fontSize: 10, padding: "1px 3px" }}
        />
      </label>
      <button
        onClick={() => applySettings(type)}
        style={{
          padding: "2px 6px", borderRadius: 4,
          border: "1px solid #cbd5e1",
          background: "#fff", cursor: "pointer", fontSize: 10,
        }}
        disabled={!pdataDirs?.[type] || !!busyByType?.[type]}
      >
        Appliquer
      </button>
    </div>
  );

  const renderControls2D = (type) => {
    const currentLevel = bgByType[type]?.meta?.level_start ?? levelByType[type] ?? 0.005;
    
    return (
      <div style={{
        display: "flex", gap: 6, alignItems: "center",
        flexWrap: "wrap", marginBottom: 6, fontSize: 11,
      }}>
        <label style={{ fontSize: 10 }}>
          W:
          <input
            type="number" min={400}
            value={resByType[type].width_px}
            onChange={(e) =>
              setResByType(p => ({
                ...p,
                [type]: {
                  ...p[type],
                  width_px: Math.max(400, +e.target.value || 400),
                },
              }))
            }
            style={{ width: 70, fontSize: 10, padding: "1px 3px" }}
          />
        </label>
        <label style={{ fontSize: 10 }}>
          H:
          <input
            type="number" min={300}
            value={resByType[type].height_px}
            onChange={(e) =>
              setResByType(p => ({
                ...p,
                [type]: {
                  ...p[type],
                  height_px: Math.max(300, +e.target.value || 300),
                },
              }))
            }
            style={{ width: 70, fontSize: 10, padding: "1px 3px" }}
          />
        </label>
        <label style={{ fontSize: 10 }}>
          DPI:
          <input
            type="number" min={72}
            value={resByType[type].dpi}
            onChange={(e) =>
              setResByType(p => ({
                ...p,
                [type]: {
                  ...p[type],
                  dpi: Math.max(72, +e.target.value || 72),
                },
              }))
            }
            style={{ width: 60, fontSize: 10, padding: "1px 3px" }}
          />
        </label>
        <label style={{ fontSize: 10 }}>
          Lvl:
          <input
            type="number" step="0.001" min="0.000001"
            value={currentLevel}
            onChange={(e) => {
              const newVal = Math.max(1e-6, +e.target.value || 0.005);
              setLevelByType(p => ({ ...p, [type]: newVal }));
              if (bgByType[type]?.meta) {
                setBgByType(prev => ({
                  ...prev,
                  [type]: {
                    ...prev[type],
                    meta: {
                      ...prev[type].meta,
                      level_start: newVal,
                      nonce: `${Date.now()}`
                    },
                  },
                }));
                bumpKey(type);
              }
            }}
            style={{ width: 70, fontSize: 10, padding: "1px 3px" }}
          />
        </label>
        <button
          onClick={() => applySettings(type)}
          style={{
            padding: "2px 6px", borderRadius: 4,
            border: "1px solid #cbd5e1",
            background: "#fff", cursor: "pointer", fontSize: 10,
          }}
          disabled={!pdataDirs?.[type] || !!busyByType?.[type]}
        >
          Appliquer
        </button>
      </div>
    );
  };

  /* =======================
     Rendu principal
     ======================= */
  return (
    <div className="viewer" style={{ marginTop: 8 }}>
      <h5 style={{ marginBottom: 8 }}>
        Vue combinée : ¹H · ¹³C · COSY · HSQC
      </h5>

      {loading && (
        <div className="status-item">Chargement…</div>
      )}
      {err && (
        <div className="status-item error">Erreur: {err}</div>
      )}

      {!loading && !err && (
        <div
          style={{
            display: "grid",
            gridTemplateColumns: hideSidebar ? "1fr" : `${SIDEBAR_W}px 1fr`,
            gap: 12, width: "100%", maxWidth: "none",
            alignItems: "start" ,
          }}
        >
          {/* ==== Sidebar ==== */}
          {!hideSidebar && (
          <aside
            style={{
              background: "#f7f7f9", borderRadius: 8, padding: 8,
               display: "flex", flexDirection: "column",
               position: "sticky", top: 8,
               maxHeight: "calc(100vh - 16px)",
               overflow: "auto",
                alignSelf: "start" 

            }}
          >
            <input
              placeholder="Source…"
              value={query}
              onChange={(e) => setQuery(e.target.value)}
              className="parameter-field"
              style={{ marginBottom: 6, fontSize: 10,
                       padding: "4px 6px", flexShrink: 0 }}
            />
            <div
              style={{
                overflow: "auto", flex: 1,
                border: "1px solid #e5e7eb",
                borderRadius: 6, background: "#fff",
              }}
            >
              {allSources.length === 0 && (
                <div style={{ padding: 6, color: "#666", fontSize: 10 }}>
                  Aucune correspondance
                </div>
              )}
              {allSources.map((s) => {
                const stripMethodPrefix = (name) => String(name).replace(/^(?:1H|13C|COSY|HSQC)[\s_\-:]+/i, "");
                const label = stripMethodPrefix(s.replace(/_predicted$/i, ""));
                const active = selectedSource === s;
                
                return (
                  <button
                    key={s}
                    onClick={() => setSelectedSource(s)}
                    className={`type-tab ${active ? "active" : ""}`}
                    style={{
                      width: "100%", textAlign: "left",
                      padding: "6px 8px",
                      borderBottom: "1px solid #f0f0f0",
                      background: active ? "#e6f4ff" : "transparent",
                      cursor: "pointer", fontSize: 10,
                    }}
                    title={label}
                  >
                    <div
                      style={{
                        fontWeight: 600,
                        whiteSpace: "normal",
                        wordBreak: "break-word",
                        lineHeight: 1.3,
                      }}
                    >
                      {label}
                    </div>
                  </button>
                );
              })}
            </div>

            <div style={{ marginTop: 8, display: "grid", gap: 3, fontSize: 9 }}>
              <label style={{ display: "inline-flex", alignItems: "center", gap: 4 }}>
                <input
                  type="checkbox" checked={showTheory}
                  onChange={(e) => setShowTheory(e.target.checked)}
                  style={{ width: 12, height: 12 }}
                /> Théorie
              </label>
              <label style={{ display: "inline-flex", alignItems: "center", gap: 4 }}>
                <input
                  type="checkbox" checked={showExp}
                  onChange={(e) => setShowExp(e.target.checked)}
                  style={{ width: 12, height: 12 }}
                /> Exp.
              </label>
              <label style={{ display: "inline-flex", alignItems: "center", gap: 4 }}>
                <input
                  type="checkbox" checked={reverse1D}
                  onChange={(e) => setReverse1D(e.target.checked)}
                  style={{ width: 12, height: 12 }}
                /> 1D↓
              </label>
              <label style={{ display: "inline-flex", alignItems: "center", gap: 4 }}>
                <input
                  type="checkbox" checked={zoomRectEnabled}
                  onChange={(e) => setZoomRectEnabled(e.target.checked)}
                  style={{ width: 12, height: 12 }}
                /> Zoom
              </label>
              <label style={{ display: "inline-flex", alignItems: "center", gap: 4 }}>
                <input
                  type="checkbox" checked={show3D}
                  onChange={(e) => setShow3D(e.target.checked)}
                  style={{ width: 12, height: 12 }}
                /> 3D
              </label>
              <label style={{ display: "inline-flex", alignItems: "center", gap: 4 }}>
                <input
                  type="checkbox" checked={relevel2DEnabled}
                  onChange={(e) => setRelevel2DEnabled(e.target.checked)}
                  style={{ width: 12, height: 12 }}
                /> Relevel 2D
              </label>
              <div style={{ fontSize: 8, color: "#475569", marginTop: 2, lineHeight: 1.2 }}>
                {(canInteract1H || canInteract13C ||
                  canInteractCOSY || canInteractHSQC)
                  ? "✓ Interactions" : "✗ Pas d'interactions"}
              </div>
            </div>
          </aside>
          )}

          {/* ==== Contenu ==== */}
          <section style={{ minWidth: 0 }}>
            {/* Molécule + boutons */}
            <div style={{ marginBottom: 12 }}>
              <div
                style={{
                  display: "flex", justifyContent: "space-between",
                  alignItems: "center", marginBottom: 6,
                }}
              >
                <div style={{ fontWeight: 700, fontSize: 12 }}>
                  Molécule — colorations combinées
                </div>
                <div style={{ display: "flex", gap: 8 }}>
                  <button
                    onClick={() => setShowAllOpen(true)}
                    style={{
                      padding: "4px 12px", borderRadius: 6,
                      border: "1px solid #0ea5e9",
                      background: "#0ea5e9", color: "#fff",
                      cursor: hasAnyData ? "pointer" : "not-allowed",
                      opacity: hasAnyData ? 1 : 0.5,
                      fontSize: 11, fontWeight: 700,
                    }}
                    title="Afficher toutes les molécules (2D)"
                    disabled={!hasAnyData}
                  >
                    Afficher toutes les molécules
                  </button>
                  <button
                    onClick={() => setColorPickerOpen(true)}
                    style={{
                      padding: "4px 12px", borderRadius: 6,
                      border: "1px solid #2563eb",
                      background: "#2563eb", color: "#fff",
                      cursor: hasAnyData ? "pointer" : "not-allowed",
                      opacity: hasAnyData ? 1 : 0.5,
                      fontSize: 11, fontWeight: 600,
                    }}
                    disabled={!hasAnyData}
                  >
                    Changer les couleurs
                  </button>
                </div>
              </div>

              {selectedSource ? (
                <MoleculePane
                  csvUrl={pivotCsvUrl}
                  selectedRow={{ Source: selectedSource }}
                  show3D={show3D}
                  taskId={taskId}
                  customColors={customColors}
                  multi={{
                    fromAll: true,
                    pivotSource: selectedSource,
                    taskId,
                    rowsByType: {
                      "1H":  rowFor("1H"),
                      "13C": rowFor("13C"),
                      "COSY":rowFor("COSY"),
                      "HSQC":rowFor("HSQC"),
                    },
                  }}
                />
              ) : (
                <div
                  style={{
                    height: MOL_HEIGHT,
                    border: "1px dashed #e5e7eb",
                    borderRadius: 8,
                    display: "flex", alignItems: "center",
                    justifyContent: "center",
                    color: "#64748b", fontSize: 12,
                  }}
                >
                  {hasAnyData ? "Aucune source sélectionnée" : "Aucune molécule — lancez une prédiction puis une annotation"}
                </div>
              )}
            </div>

            {/* Spectres */}
            {(() => {
              const blocks = [];

              if (oneH) {
                const bg1H  = bgByType?.["1H"];
                const img   = bg1H?.pngUrl || null;
                const meta  = getMeta(bg1H);

                // Ajout : ppmPart pour refléter la fenêtre dans la key
                const ppmPart1H = meta
                  ? [meta.ppm_min, meta.ppm_max, meta.ppm_left, meta.ppm_right]
                      .map(v => v ?? "").join(",")
                  : "";

                blocks.push({
                  key: "b-1H",
                  el: (
                    <div
                      style={{
                        border: "1px solid #e5e7eb",
                        borderRadius: 8, padding: 8,
                        overflow: "hidden",
                      }}
                    >
                      <div
                        style={{
                          fontWeight: 700, marginBottom: 4,
                          display: "flex", justifyContent: "space-between",
                          alignItems: "center", fontSize: 12,
                        }}
                      >
                        <span>{TYPE_LABELS["1H"]}</span>
                        {img && (
                          <span
                            style={{
                              fontSize: 9, color: "#16a34a",
                              background: "#dcfce7",
                              padding: "1px 4px", borderRadius: 3,
                            }}
                          >
                            ✓
                          </span>
                        )}
                      </div>

                      {oneH.metrics &&
                        Object.keys(oneH.metrics).length > 0 && (
                        <div style={{ marginBottom: 6, minHeight: 20 }}>
                          {oneH.metrics.coverage_ransy && (
                            <Pill label="RANSY" value={oneH.metrics.coverage_ransy}/>
                          )}
                          {oneH.metrics.coverage_sub && (
                            <Pill label="Sub" value={oneH.metrics.coverage_sub}/>
                          )}
                          {oneH.metrics.avg_variabilite && (
                            <Pill label="Var" value={oneH.metrics.avg_variabilite}/>
                          )}
                        </div>
                      )}

                      {renderControls1D("1H")}

                      <StickChart
                        key={[
                          refreshKeyByType["1H"] ?? 0,
                          img || "",
                          meta?.etag ?? meta?.version ?? meta?.last_modified ?? "",
                          JSON.stringify(meta?.plot_bbox_px ?? null) || "",
                          ppmPart1H,
                          String(meta?.nonce ?? "")
                        ].join("|")}
                        theory={showTheory ? oneH.theory : []}
                        expPpm={showExp ? oneH.expX : []}
                        expInt={showExp ? oneH.expI : []}
                        reverse={reverse1D}
                        showTheory={showTheory}
                        showExp={showExp}
                        height={PLOT_HEIGHT}
                        backgroundImage={img}
                        backgroundMeta={meta}
                        canInteract={canInteract1H}
                        zoomEnabled={zoomRectEnabled}
                        onZoom={(a, b) => doZoom1D("1H", a, b)}
                        onZoomPx={(r, d) => doZoom1DPixels("1H", r, d)}
                        onReset={() => doReset1D("1H")}
                        busy={!!busyByType?.["1H"]}
                      />
                    </div>
                  ),
                });
              }

              if (c13) {
                const bg13  = bgByType?.["13C"];
                const img   = bg13?.pngUrl || null;
                const meta  = getMeta(bg13);

                const ppmPart13 = meta
                  ? [meta.ppm_min, meta.ppm_max, meta.ppm_left, meta.ppm_right]
                      .map(v => v ?? "").join(",")
                  : "";

                blocks.push({
                  key: "b-13C",
                  el: (
                    <div
                      style={{
                        border: "1px solid #e5e7eb",
                        borderRadius: 8, padding: 8,
                        overflow: "hidden"
                      }}
                    >
                      <div
                        style={{
                          fontWeight: 700, marginBottom: 4,
                          display: "flex", justifyContent: "space-between",
                          alignItems: "center", fontSize: 12,
                        }}
                      >
                        <span>{TYPE_LABELS["13C"]}</span>
                        {img && (
                          <span
                            style={{
                              fontSize: 9, color: "#16a34a",
                              background: "#dcfce7",
                              padding: "1px 4px", borderRadius: 3,
                            }}
                          >
                            ✓
                          </span>
                        )}
                      </div>

                      {c13.metrics &&
                        Object.keys(c13.metrics).length > 0 && (
                        <div style={{ marginBottom: 6, minHeight: 20 }}>
                          {c13.metrics.coverage_ransy && (
                            <Pill label="RANSY" value={c13.metrics.coverage_ransy}/>
                          )}
                          {c13.metrics.coverage_sub && (
                            <Pill label="Sub" value={c13.metrics.coverage_sub}/>
                          )}
                          {c13.metrics.avg_variabilite && (
                            <Pill label="Var" value={c13.metrics.avg_variabilite}/>
                          )}
                        </div>
                      )}

                      {renderControls1D("13C")}

                      <StickChart
                        key={[
                          refreshKeyByType["13C"] ?? 0,
                          img || "",
                          meta?.etag ?? meta?.version ?? meta?.last_modified ?? "",
                          JSON.stringify(meta?.plot_bbox_px ?? null) || "",
                          ppmPart13,
                          String(meta?.nonce ?? "")
                        ].join("|")}
                        theory={showTheory ? c13.theory : []}
                        expPpm={showExp ? c13.expX : []}
                        expInt={showExp ? c13.expI : []}
                        reverse={reverse1D}
                        showTheory={showTheory}
                        showExp={showExp}
                        height={PLOT_HEIGHT}
                        backgroundImage={img}
                        backgroundMeta={meta}
                        canInteract={canInteract13C}
                        zoomEnabled={zoomRectEnabled}
                        onZoom={(a, b) => doZoom1D("13C", a, b)}
                        onZoomPx={(r, d) => doZoom1DPixels("13C", r, d)}
                        onReset={() => doReset1D("13C")}
                        busy={!!busyByType?.["13C"]}
                      />
                    </div>
                  ),
                });
              }

              if (cosy) {
                const bgC  = bgByType?.["COSY"];
                const img  = bgC?.pngUrl || null;
                const meta = getMeta(bgC);

                blocks.push({
                  key: "b-COSY",
                  el: (
                    <div
                      style={{
                        border: "1px solid #e5e7eb",
                        borderRadius: 8, padding: 8,
                        overflow: "hidden"
                      }}
                    >
                      <div
                        style={{
                          fontWeight: 700, marginBottom: 4,
                          display: "flex", justifyContent: "space-between",
                          alignItems: "center", fontSize: 12,
                        }}
                      >
                        <span>COSY</span>
                        {img && (
                          <span
                            style={{
                              fontSize: 9, color: "#16a34a",
                              background: "#dcfce7",
                              padding: "1px 4px", borderRadius: 3,
                            }}
                          >
                            ✓
                          </span>
                        )}
                      </div>

                      {cosy.metrics &&
                        Object.keys(cosy.metrics).length > 0 && (
                        <div style={{ marginBottom: 6, minHeight: 20 }}>
                          {cosy.metrics.coverage_ransy && (
                            <Pill label="RANSY" value={cosy.metrics.coverage_ransy}/>
                          )}
                          {cosy.metrics.coverage_sub && (
                            <Pill label="Sub" value={cosy.metrics.coverage_sub}/>
                          )}
                          {cosy.metrics.avg_variabilite && (
                            <Pill label="Var" value={cosy.metrics.avg_variabilite}/>
                          )}
                        </div>
                      )}

                      {renderControls2D("COSY")}

                      <Scatter2D
                        key={[
                          refreshKeyByType["COSY"] ?? 0,
                          img || "",
                          meta?.etag ?? meta?.version ?? meta?.last_modified ?? "",
                          JSON.stringify(meta?.plot_bbox_px ?? null) || "",
                          String(meta?.level_start ?? ""),
                          String(meta?.nonce ?? "")
                        ].join("|")}
                        theoryPts={showTheory ? cosy.theory : []}
                        expPts={showExp ? cosy.expPts : []}
                        expInt={showExp ? cosy.expI : []}
                        showTheory={showTheory}
                        showExp={showExp}
                        height={PLOT_HEIGHT}
                        backgroundImage={img}
                        backgroundMeta={meta}
                        canInteract={canInteractCOSY}
                        zoomEnabled={zoomRectEnabled}
                        onZoom2D={(w) => doZoom2D("COSY", w)}
                        onWheelDelta={(dy) => doWheel2D("COSY", dy)}
                        onReset={() => doReset2D("COSY")}
                        busy={!!busyByType?.["COSY"]}
                      />
                    </div>
                  ),
                });
              }

              if (hsqc) {
                const bgH  = bgByType?.["HSQC"];
                const img  = bgH?.pngUrl || null;
                const meta = getMeta(bgH);

                blocks.push({
                  key: "b-HSQC",
                  el: (
                    <div
                      style={{
                        border: "1px solid #e5e7eb",
                        borderRadius: 8, padding: 8,
                        overflow: "hidden"
                      }}
                    >
                      <div
                        style={{
                          fontWeight: 700, marginBottom: 4,
                          display: "flex", justifyContent: "space-between",
                          alignItems: "center", fontSize: 12,
                        }}
                      >
                        <span>HSQC</span>
                        {img && (
                          <span
                            style={{
                              fontSize: 9, color: "#16a34a",
                              background: "#dcfce7",
                              padding: "1px 4px", borderRadius: 3,
                            }}
                          >
                            ✓
                          </span>
                        )}
                      </div>

                      {hsqc.metrics &&
                        Object.keys(hsqc.metrics).length > 0 && (
                        <div style={{ marginBottom: 6, minHeight: 20 }}>
                          {hsqc.metrics.coverage_ransy && (
                            <Pill label="RANSY" value={hsqc.metrics.coverage_ransy}/>
                          )}
                          {hsqc.metrics.coverage_sub && (
                            <Pill label="Sub" value={hsqc.metrics.coverage_sub}/>
                          )}
                          {hsqc.metrics.avg_variabilite && (
                            <Pill label="Var" value={hsqc.metrics.avg_variabilite}/>
                          )}
                        </div>
                      )}

                      {renderControls2D("HSQC")}

                      <Scatter2D
                        key={[
                          refreshKeyByType["HSQC"] ?? 0,
                          img || "",
                          meta?.etag ?? meta?.version ?? meta?.last_modified ?? "",
                          JSON.stringify(meta?.plot_bbox_px ?? null) || "",
                          String(meta?.level_start ?? ""),
                          String(meta?.nonce ?? "")
                        ].join("|")}
                        theoryPts={showTheory ? hsqc.theory : []}
                        expPts={showExp ? hsqc.expPts : []}
                        expInt={showExp ? hsqc.expI : []}
                        showTheory={showTheory}
                        showExp={showExp}
                        height={PLOT_HEIGHT}
                        backgroundImage={img}
                        backgroundMeta={meta}
                        canInteract={canInteractHSQC}
                        zoomEnabled={zoomRectEnabled}
                        onZoom2D={(w) => doZoom2D("HSQC", w)}
                        onWheelDelta={(dy) => doWheel2D("HSQC", dy)}
                        onReset={() => doReset2D("HSQC")}
                        busy={!!busyByType?.["HSQC"]}
                      />
                    </div>
                  ),
                });
              }

              const n = blocks.length;
              if (n === 0) {
                if (!hasAnyData) {
                  const emptyBox = (title) => (
                    <div
                      style={{
                        border: "1px solid #e5e7eb",
                        borderRadius: 8,
                        padding: 8,
                        overflow: "hidden",
                        background: "#fff",
                      }}
                    >
                      <div
                        style={{
                          fontWeight: 700,
                          marginBottom: 6,
                          display: "flex",
                          justifyContent: "space-between",
                          alignItems: "center",
                          fontSize: 12,
                        }}
                      >
                        <span>{title}</span>
                      </div>
                      <div
                        style={{
                          height: PLOT_HEIGHT,
                          border: "1px dashed #e5e7eb",
                          borderRadius: 8,
                          display: "flex",
                          alignItems: "center",
                          justifyContent: "center",
                          color: "#64748b",
                          fontSize: 12,
                          background: "#fff",
                        }}
                      >
                        Aucun spectre — lancez l'annotation
                      </div>
                    </div>
                  );

                  return (
                    <div
                      style={{
                        display: "grid",
                        gridTemplateColumns: "1fr 1fr",
                        gridTemplateRows: "1fr 1fr",
                        gap: 12,
                        marginBottom: 12,
                        alignItems: "start",
                        gridAutoRows: "min-content",
                      }}
                    >
                      <div>{emptyBox(TYPE_LABELS["1H"])}</div>
                      <div>{emptyBox(TYPE_LABELS["13C"])}</div>
                      <div>{emptyBox("COSY")}</div>
                      <div>{emptyBox("HSQC")}</div>
                    </div>
                  );
                }

                return <div>Aucun spectre disponible</div>;
              }
              if (n === 1) return <div style={{ marginBottom: 12 }}>{blocks[0].el}</div>;
              if (n === 2)
                return (
                  <div
                    style={{
                      display: "grid",
                      gridTemplateColumns: "1fr 1fr",
                      gap: 12, marginBottom: 12,
                      alignItems: "start",
                      gridAutoRows: "min-content"
                    }}
                  >
                    {blocks.map(b => <div key={b.key}>{b.el}</div>)}
                  </div>
                );
              if (n === 3)
                return (
                  <div style={{ marginBottom: 12 }}>
                    <div
                      style={{
                        display: "grid",
                        gridTemplateColumns: "1fr 1fr",
                        gap: 12, marginBottom: 12,
                        alignItems: "start",
                        gridAutoRows: "min-content"
                      }}
                    >
                      {blocks.slice(0, 2)
                        .map(b => <div key={b.key}>{b.el}</div>)}
                    </div>
                    <div style={{ display: "flex", justifyContent: "center" }}>
                      <div style={{ width: "48%" }}>{blocks[2].el}</div>
                    </div>
                  </div>
                );
              return (
                <div
                  style={{
                    display: "grid",
                    gridTemplateColumns: "1fr 1fr",
                    gridTemplateRows: "1fr 1fr",
                    gap: 12, marginBottom: 12,
                    alignItems: "start",
                    gridAutoRows: "min-content"
                  }}
                >
                  {blocks.map(b => <div key={b.key}>{b.el}</div>)}
                </div>
              );
            })()}
          </section>
        </div>
      )}

      {/* Sélecteur de couleurs (ouvert via bouton) */}
      <ColorPicker
        isOpen={colorPickerOpen}
        onClose={() => setColorPickerOpen(false)}
        customColors={customColors}
        onColorChange={setCustomColor}
      />

      {/* Modale — toutes les molécules */}
      <ShowAllMoleculesModal
        isOpen={showAllOpen}
        onClose={() => setShowAllOpen(false)}
        sources={sortedSourcesAll}
        page={allPage}
        setPage={setAllPage}
        pageSize={allPageSize}
        setPageSize={setAllPageSize}
        getRowsByType={rowsByTypeForSource}
        csvUrl={pivotCsvUrl}
        taskId={taskId}
        customColors={customColors}
      />
    </div>
  );
}

if (typeof window !== 'undefined' && !window.__NMR_ERROR_FILTER_INSTALLED__) {
  window.__NMR_ERROR_FILTER_INSTALLED__ = true;
  
  const originalError = console.error;
  console.error = function(...args) {
    const msg = args.join(' ');
    
    // Filtrer les erreurs WebGL connues et non-critiques
    if (msg.includes('GL_INVALID_FRAMEBUFFER_OPERATION') ||
        msg.includes('Framebuffer is incomplete') ||
        msg.includes('Attachment has zero size')) {
      if (window.DEBUG_NMR) {
        originalError.apply(console, ['[WebGL Filtré]', ...args]);
      }
      return;
    }
    
    // Afficher les vraies erreurs
    originalError.apply(console, args);
  };
  
  debugLog('✅ Filtre d\'erreurs WebGL installé');
}

