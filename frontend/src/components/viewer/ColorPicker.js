// src/components/viewer/ColorPicker.js
// Advanced color picker Word-style

import React, { useState, useEffect } from "react";
import {
  COLOR_OPTIONS,
  DEFAULT_COLORS,
  METHOD_LABELS_SHORT,
} from "../../constants/colors";
import {
  hexToRgb255,
  rgb255ToHex,
  rgb255ToFloat,
  hslToHex,
  hexToHsl,
  clamp,
} from "../../utils/colorUtils";

/**
 * Advanced color selection component
 */
function ColorPicker({ isOpen, onClose, customColors, onColorChange }) {
  // Merge with defaults so every method always has hex/rgb fields
  const safeColors = {
    H1:   { ...DEFAULT_COLORS.H1,   ...customColors?.H1 },
    C13:  { ...DEFAULT_COLORS.C13,  ...customColors?.C13 },
    COSY: { ...DEFAULT_COLORS.COSY, ...customColors?.COSY },
    HSQC: { ...DEFAULT_COLORS.HSQC, ...customColors?.HSQC },
  };
  const [draft, setDraft] = useState(safeColors);
  const [active, setActive] = useState("H1");
  const [showColorPalette, setShowColorPalette] = useState(false);
  const [tempColor, setTempColor] = useState({ h: 180, s: 70, l: 50 });

  useEffect(() => {
    setDraft({
      H1:   { ...DEFAULT_COLORS.H1,   ...customColors?.H1 },
      C13:  { ...DEFAULT_COLORS.C13,  ...customColors?.C13 },
      COSY: { ...DEFAULT_COLORS.COSY, ...customColors?.COSY },
      HSQC: { ...DEFAULT_COLORS.HSQC, ...customColors?.HSQC },
    });
  }, [customColors]);

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

  const MethodRow = ({ method, title }) => {
    const current = safeColors[method];
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
              <div style={{ fontWeight: 700, fontSize: 14 }}>{title}</div>
              <div style={{ fontSize: 11, color: "#64748b" }}>{d.name}</div>
            </div>
          </div>

          <select
            value={d.name}
            onChange={(e) => {
              const chosen =
                COLOR_OPTIONS.find((c) => c.name === e.target.value) || {
                  name: "Custom",
                  hex: d.hex,
                  rgb: d.rgb,
                };
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
            title="Quick palette"
          >
            <option value={d.name} disabled>
              — Quick palette —
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
          {/* Color preview */}
          <div
            style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 10 }}
          >
            <div>
              <div style={{ fontSize: 11, marginBottom: 6, color: "#475569" }}>
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
              <div style={{ fontSize: 11, marginBottom: 6, color: "#475569" }}>
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

            {/* Button to open custom palette */}
            <div
              style={{
                gridColumn: "1 / span 2",
                marginTop: 8,
                position: "relative",
              }}
            >
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
                Choose ▼
              </button>

              {/* Custom palette */}
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
                  {/* Gradient 2D */}
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
                        background:
                          "linear-gradient(to right, #fff 0%, transparent 100%)",
                        borderRadius: 6,
                      }}
                    />
                  </div>

                  {/* Slider Teinte */}
                  <div style={{ marginBottom: 8 }}>
                    <div
                      style={{ fontSize: 10, marginBottom: 4, color: "#475569" }}
                    >
                      Hue: {tempColor.h}°
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
                        background:
                          "linear-gradient(to right, #ff0000 0%, #ffff00 17%, #00ff00 33%, #00ffff 50%, #0000ff 67%, #ff00ff 83%, #ff0000 100%)",
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
                    Close
                  </button>
                </div>
              )}
            </div>
          </div>

          {/* Hex and RGB fields */}
          <div
            style={{
              display: "grid",
              gridTemplateColumns: "120px 1fr",
              rowGap: 8,
              columnGap: 10,
            }}
          >
            <label
              style={{ fontSize: 12, color: "#334155", alignSelf: "center" }}
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
                Apply
              </button>
            </div>

            <label
              style={{ fontSize: 12, color: "#334155", alignSelf: "center" }}
            >
              RGB :
            </label>
            <div style={{ display: "flex", gap: 8 }}>
              {["r", "g", "b"].map((channel) => (
                <input
                  key={channel}
                  type="number"
                  min={0}
                  max={255}
                  value={newRGB[channel]}
                  onChange={(e) => {
                    const next = {
                      ...newRGB,
                      [channel]: clamp(+e.target.value),
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
              ))}
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
            Apply {method}
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
            Reset {method}
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
            Customize colors (Word-style)
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

        {active === "H1" && <MethodRow method="H1" title="¹H NMR" />}
        {active === "C13" && <MethodRow method="C13" title="¹³C NMR" />}
        {active === "COSY" && <MethodRow method="COSY" title="COSY" />}
        {active === "HSQC" && <MethodRow method="HSQC" title="HSQC" />}

        <div
          style={{
            display: "flex",
            justifyContent: "flex-end",
            gap: 8,
            marginTop: 14,
          }}
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
            Reset all
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
            Done
          </button>
        </div>
      </div>
    </div>
  );
}

export default ColorPicker;
