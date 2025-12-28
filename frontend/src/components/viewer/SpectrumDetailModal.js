// src/components/viewer/SpectrumDetailModal.js
// Modal pour afficher le spectre 1H simulé avec formes de raie

import React, { useRef, useEffect, useState, useCallback } from "react";
import { simulateSpectrum } from "../../api";

const MODAL_STYLES = {
  overlay: {
    position: "fixed",
    top: 0,
    left: 0,
    right: 0,
    bottom: 0,
    backgroundColor: "rgba(0, 0, 0, 0.7)",
    display: "flex",
    alignItems: "center",
    justifyContent: "center",
    zIndex: 10000,
  },
  modal: {
    backgroundColor: "#fff",
    borderRadius: 12,
    padding: 20,
    width: "90%",
    maxWidth: 1200,
    maxHeight: "90vh",
    overflow: "hidden",
    display: "flex",
    flexDirection: "column",
    boxShadow: "0 20px 60px rgba(0,0,0,0.3)",
  },
  header: {
    display: "flex",
    justifyContent: "space-between",
    alignItems: "center",
    marginBottom: 16,
    paddingBottom: 12,
    borderBottom: "1px solid #e5e7eb",
  },
  title: {
    fontSize: 18,
    fontWeight: 700,
    color: "#1e293b",
  },
  closeBtn: {
    background: "none",
    border: "none",
    fontSize: 24,
    cursor: "pointer",
    color: "#64748b",
    padding: "4px 8px",
    borderRadius: 4,
  },
  controls: {
    display: "flex",
    gap: 16,
    marginBottom: 12,
    flexWrap: "wrap",
    alignItems: "center",
  },
  controlGroup: {
    display: "flex",
    alignItems: "center",
    gap: 6,
  },
  label: {
    fontSize: 12,
    fontWeight: 500,
    color: "#64748b",
  },
  select: {
    padding: "4px 8px",
    borderRadius: 4,
    border: "1px solid #d1d5db",
    fontSize: 12,
  },
  input: {
    padding: "4px 8px",
    borderRadius: 4,
    border: "1px solid #d1d5db",
    fontSize: 12,
    width: 70,
  },
  canvasContainer: {
    flex: 1,
    minHeight: 400,
    position: "relative",
    backgroundColor: "#f8fafc",
    borderRadius: 8,
    overflow: "hidden",
  },
  loading: {
    position: "absolute",
    top: "50%",
    left: "50%",
    transform: "translate(-50%, -50%)",
    display: "flex",
    flexDirection: "column",
    alignItems: "center",
    gap: 12,
  },
  spinner: {
    width: 40,
    height: 40,
    border: "3px solid #e5e7eb",
    borderTop: "3px solid #3b82f6",
    borderRadius: "50%",
    animation: "spin 1s linear infinite",
  },
  error: {
    color: "#dc2626",
    textAlign: "center",
    padding: 20,
  },
  peaksTable: {
    marginTop: 12,
    maxHeight: 150,
    overflow: "auto",
    fontSize: 11,
    border: "1px solid #e5e7eb",
    borderRadius: 6,
  },
  table: {
    width: "100%",
    borderCollapse: "collapse",
  },
  th: {
    padding: "6px 8px",
    textAlign: "left",
    backgroundColor: "#f1f5f9",
    fontWeight: 600,
    position: "sticky",
    top: 0,
  },
  td: {
    padding: "4px 8px",
    borderTop: "1px solid #e5e7eb",
  },
};

// CSS keyframes for spinner
const spinnerKeyframes = `
  @keyframes spin {
    0% { transform: rotate(0deg); }
    100% { transform: rotate(360deg); }
  }
`;

function SpectrumDetailModal({
  isOpen,
  onClose,
  smiles,
  molId = "mol",
  theoryShifts = [],  // Positions théoriques (traits bleus)
  expShifts = [],      // Positions expérimentales (traits rouges)
  multiplicityJson = null,  // JSON de multiplicité pré-calculé (recommandé)
}) {
  const canvasRef = useRef(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [spectrumData, setSpectrumData] = useState(null);
  
  // Options de simulation
  const [freqMhz, setFreqMhz] = useState(700);
  const [lineShape, setLineShape] = useState("lorentzian");
  const [ppmMin, setPpmMin] = useState(-0.5);
  const [ppmMax, setPpmMax] = useState(12);
  const [lineWidthHz, setLineWidthHz] = useState(1.0);
  
  // Zoom state
  const [zoomRange, setZoomRange] = useState(null);

  // Fetch spectrum data
  const fetchSpectrum = useCallback(async () => {
    // On a besoin soit du JSON de multiplicité, soit du SMILES
    if (!multiplicityJson && !smiles) {
      setError("No multiplicity data or SMILES available");
      return;
    }
    
    setLoading(true);
    setError(null);
    
    try {
      const data = await simulateSpectrum({
        // Priorité au JSON de multiplicité (plus précis)
        multiplicity_json: multiplicityJson,
        // Fallback sur SMILES si pas de JSON
        smiles: multiplicityJson ? undefined : smiles,
        mol_id: molId,
        freq_mhz: freqMhz,
        line_shape: lineShape,
        ppm_range: [ppmMin, ppmMax],
        digital_resolution_hz: 0.1,
        line_width_hz: lineWidthHz,
      });
      
      setSpectrumData(data);
      setZoomRange(null);  // Reset zoom on new data
    } catch (err) {
      console.error("[SpectrumDetailModal] Error:", err);
      setError(err.message || "Failed to simulate spectrum");
    } finally {
      setLoading(false);
    }
  }, [smiles, molId, freqMhz, lineShape, ppmMin, ppmMax, multiplicityJson, lineWidthHz]);

  // Fetch on mount or when params change
  useEffect(() => {
    if (isOpen && (smiles || multiplicityJson)) {
      fetchSpectrum();
    }
  }, [isOpen, smiles, multiplicityJson, fetchSpectrum]);

  // Draw spectrum on canvas
  useEffect(() => {
    if (!spectrumData || !canvasRef.current) return;
    
    const canvas = canvasRef.current;
    const ctx = canvas.getContext("2d");
    const DPR = window.devicePixelRatio || 1;
    
    const container = canvas.parentElement;
    const W = container.clientWidth || 800;
    const H = container.clientHeight || 400;
    
    canvas.width = W * DPR;
    canvas.height = H * DPR;
    canvas.style.width = `${W}px`;
    canvas.style.height = `${H}px`;
    ctx.scale(DPR, DPR);
    
    const PAD = { top: 20, right: 30, bottom: 40, left: 50 };
    const plotW = W - PAD.left - PAD.right;
    const plotH = H - PAD.top - PAD.bottom;
    
    // Clear
    ctx.fillStyle = "#fff";
    ctx.fillRect(0, 0, W, H);
    
    // Determine PPM range (use zoom if set)
    const range = zoomRange || { min: ppmMin, max: ppmMax };
    const ppmToX = (ppm) => PAD.left + ((range.max - ppm) / (range.max - range.min)) * plotW;
    
    // Find data indices for current range
    const ppmArr = spectrumData.ppm;
    const intArr = spectrumData.intensity;
    
    // Filter and get max intensity in range
    let maxInt = 0;
    const visibleData = [];
    for (let i = 0; i < ppmArr.length; i++) {
      const ppm = ppmArr[i];
      if (ppm >= range.min && ppm <= range.max) {
        visibleData.push({ ppm, int: intArr[i] });
        if (intArr[i] > maxInt) maxInt = intArr[i];
      }
    }
    if (maxInt === 0) maxInt = 1;
    
    const intToY = (int) => PAD.top + plotH - (int / maxInt) * plotH * 0.9;
    
    // Draw grid
    ctx.strokeStyle = "#e5e7eb";
    ctx.lineWidth = 1;
    
    // Vertical grid lines (ppm)
    const ppmStep = (range.max - range.min) > 6 ? 1 : 0.5;
    for (let ppm = Math.ceil(range.min); ppm <= range.max; ppm += ppmStep) {
      const x = ppmToX(ppm);
      ctx.beginPath();
      ctx.moveTo(x, PAD.top);
      ctx.lineTo(x, PAD.top + plotH);
      ctx.stroke();
    }
    
    // Draw spectrum line
    ctx.strokeStyle = "#1e40af";
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    let started = false;
    for (const pt of visibleData) {
      const x = ppmToX(pt.ppm);
      const y = intToY(pt.int);
      if (!started) {
        ctx.moveTo(x, y);
        started = true;
      } else {
        ctx.lineTo(x, y);
      }
    }
    ctx.stroke();
    
    // Fill under curve (light blue)
    ctx.fillStyle = "rgba(59, 130, 246, 0.1)";
    ctx.beginPath();
    ctx.moveTo(ppmToX(range.max), intToY(0));
    for (const pt of visibleData) {
      ctx.lineTo(ppmToX(pt.ppm), intToY(pt.int));
    }
    ctx.lineTo(ppmToX(range.min), intToY(0));
    ctx.closePath();
    ctx.fill();
    
    // Draw theory markers (blue sticks)
    if (theoryShifts && theoryShifts.length > 0) {
      ctx.strokeStyle = "#3b82f6";
      ctx.lineWidth = 2;
      for (const ppm of theoryShifts) {
        if (ppm >= range.min && ppm <= range.max) {
          const x = ppmToX(ppm);
          ctx.beginPath();
          ctx.moveTo(x, PAD.top + plotH);
          ctx.lineTo(x, PAD.top + plotH - 30);
          ctx.stroke();
        }
      }
    }
    
    // Draw experimental markers (red sticks)
    if (expShifts && expShifts.length > 0) {
      ctx.strokeStyle = "#dc2626";
      ctx.lineWidth = 2;
      for (const ppm of expShifts) {
        if (ppm >= range.min && ppm <= range.max) {
          const x = ppmToX(ppm);
          ctx.beginPath();
          ctx.moveTo(x, PAD.top + plotH);
          ctx.lineTo(x, PAD.top + plotH - 20);
          ctx.stroke();
        }
      }
    }
    
    // Draw axes
    ctx.strokeStyle = "#374151";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(PAD.left, PAD.top);
    ctx.lineTo(PAD.left, PAD.top + plotH);
    ctx.lineTo(PAD.left + plotW, PAD.top + plotH);
    ctx.stroke();
    
    // Draw ppm labels
    ctx.fillStyle = "#374151";
    ctx.font = "11px system-ui, sans-serif";
    ctx.textAlign = "center";
    for (let ppm = Math.ceil(range.min); ppm <= range.max; ppm += ppmStep) {
      const x = ppmToX(ppm);
      ctx.fillText(ppm.toFixed(1), x, PAD.top + plotH + 15);
    }
    
    // X-axis label
    ctx.fillStyle = "#64748b";
    ctx.font = "12px system-ui, sans-serif";
    ctx.fillText("Chemical Shift (ppm)", PAD.left + plotW / 2, H - 8);
    
    // Title
    ctx.fillStyle = "#1e293b";
    ctx.font = "bold 13px system-ui, sans-serif";
    ctx.textAlign = "left";
    ctx.fillText(`¹H NMR Spectrum - ${molId}`, PAD.left, 14);
    
    // Legend
    const legendX = W - 180;
    ctx.font = "11px system-ui, sans-serif";
    ctx.fillStyle = "#3b82f6";
    ctx.fillRect(legendX, 8, 12, 3);
    ctx.fillStyle = "#374151";
    ctx.fillText("Predicted", legendX + 18, 12);
    
    ctx.fillStyle = "#dc2626";
    ctx.fillRect(legendX + 80, 8, 12, 3);
    ctx.fillStyle = "#374151";
    ctx.fillText("Experimental", legendX + 98, 12);
    
  }, [spectrumData, theoryShifts, expShifts, zoomRange, ppmMin, ppmMax, molId]);

  // Handle zoom with mouse wheel
  const handleWheel = useCallback((e) => {
    e.preventDefault();
    if (!spectrumData) return;
    
    const canvas = canvasRef.current;
    const rect = canvas.getBoundingClientRect();
    const mouseX = e.clientX - rect.left;
    
    const PAD = { left: 50, right: 30 };
    const plotW = rect.width - PAD.left - PAD.right;
    
    const range = zoomRange || { min: ppmMin, max: ppmMax };
    const centerPpm = range.max - ((mouseX - PAD.left) / plotW) * (range.max - range.min);
    
    const zoomFactor = e.deltaY < 0 ? 0.8 : 1.25;
    const newWidth = (range.max - range.min) * zoomFactor;
    
    // Limit zoom
    if (newWidth < 0.5 || newWidth > 15) return;
    
    let newMin = centerPpm - newWidth / 2;
    let newMax = centerPpm + newWidth / 2;
    
    // Clamp to valid range
    if (newMin < 0) { newMin = 0; newMax = newWidth; }
    if (newMax > 15) { newMax = 15; newMin = 15 - newWidth; }
    
    setZoomRange({ min: newMin, max: newMax });
  }, [spectrumData, zoomRange, ppmMin, ppmMax]);

  // Reset zoom
  const resetZoom = () => setZoomRange(null);

  if (!isOpen) return null;

  return (
    <>
      <style>{spinnerKeyframes}</style>
      <div style={MODAL_STYLES.overlay} onClick={onClose}>
        <div style={MODAL_STYLES.modal} onClick={(e) => e.stopPropagation()}>
          <div style={MODAL_STYLES.header}>
            <span style={MODAL_STYLES.title}>
              ¹H NMR Spectrum Detail - {molId}
            </span>
            <button style={MODAL_STYLES.closeBtn} onClick={onClose}>×</button>
          </div>
          
          <div style={MODAL_STYLES.controls}>
            <div style={MODAL_STYLES.controlGroup}>
              <span style={MODAL_STYLES.label}>Frequency:</span>
              <select
                style={MODAL_STYLES.select}
                value={freqMhz}
                onChange={(e) => setFreqMhz(Number(e.target.value))}
              >
                <option value={400}>400 MHz</option>
                <option value={500}>500 MHz</option>
                <option value={600}>600 MHz</option>
                <option value={700}>700 MHz</option>
                <option value={800}>800 MHz</option>
                <option value={900}>900 MHz</option>
              </select>
            </div>
            
            <div style={MODAL_STYLES.controlGroup}>
              <span style={MODAL_STYLES.label}>Line Shape:</span>
              <select
                style={MODAL_STYLES.select}
                value={lineShape}
                onChange={(e) => setLineShape(e.target.value)}
              >
                <option value="lorentzian">Lorentzian</option>
                <option value="gaussian">Gaussian</option>
                <option value="voigt">Voigt</option>
              </select>
            </div>
            
            <div style={MODAL_STYLES.controlGroup}>
              <span style={MODAL_STYLES.label}>Line Width:</span>
              <input
                type="number"
                style={{ ...MODAL_STYLES.input, width: 60 }}
                value={lineWidthHz}
                onChange={(e) => setLineWidthHz(Number(e.target.value))}
                min={0.1}
                max={10}
                step={0.1}
              />
              <span style={{ fontSize: 11, color: "#6b7280" }}>Hz</span>
            </div>
            
            <div style={MODAL_STYLES.controlGroup}>
              <span style={MODAL_STYLES.label}>PPM Range:</span>
              <input
                type="number"
                style={MODAL_STYLES.input}
                value={ppmMin}
                onChange={(e) => setPpmMin(Number(e.target.value))}
                step={0.5}
              />
              <span>-</span>
              <input
                type="number"
                style={MODAL_STYLES.input}
                value={ppmMax}
                onChange={(e) => setPpmMax(Number(e.target.value))}
                step={0.5}
              />
            </div>
            
            <button
              onClick={fetchSpectrum}
              disabled={loading}
              style={{
                padding: "6px 12px",
                backgroundColor: "#3b82f6",
                color: "#fff",
                border: "none",
                borderRadius: 4,
                cursor: loading ? "not-allowed" : "pointer",
                fontSize: 12,
                fontWeight: 500,
              }}
            >
              {loading ? "Simulating..." : "Refresh"}
            </button>
            
            <div style={{ display: "flex", gap: 8, alignItems: "center" }}>
              {zoomRange && (
                <button
                  onClick={resetZoom}
                  style={{
                    padding: "6px 12px",
                    backgroundColor: "#f1f5f9",
                    color: "#374151",
                    border: "1px solid #d1d5db",
                    borderRadius: 4,
                    cursor: "pointer",
                    fontSize: 12,
                  }}
                >
                  Reset Zoom
                </button>
              )}
              <button
                onClick={onClose}
                style={{
                  padding: "6px 12px",
                  backgroundColor: "#3b82f6",
                  color: "#fff",
                  border: "none",
                  borderRadius: 4,
                  cursor: "pointer",
                  fontSize: 12,
                  fontWeight: 500,
                }}
              >
                View Detail
              </button>
            </div>
          </div>

          <div style={MODAL_STYLES.canvasContainer}>
            {loading && (
              <div style={MODAL_STYLES.loading}>
                <div style={MODAL_STYLES.spinner}></div>
                <span>Simulating spectrum...</span>
              </div>
            )}
            
            {error && (
              <div style={MODAL_STYLES.error}>
                <p>❌ {error}</p>
                <button
                  onClick={fetchSpectrum}
                  style={{
                    marginTop: 12,
                    padding: "8px 16px",
                    backgroundColor: "#3b82f6",
                    color: "#fff",
                    border: "none",
                    borderRadius: 4,
                    cursor: "pointer",
                  }}
                >
                  Retry
                </button>
              </div>
            )}
            
            {!loading && !error && (
              <canvas
                ref={canvasRef}
                onWheel={handleWheel}
                style={{ display: "block", width: "100%", height: "100%" }}
              />
            )}
          </div>
          
          {spectrumData?.peaks && spectrumData.peaks.length > 0 && (
            <div style={MODAL_STYLES.peaksTable}>
              <table style={MODAL_STYLES.table}>
                <thead>
                  <tr>
                    <th style={MODAL_STYLES.th}>#</th>
                    <th style={MODAL_STYLES.th}>δ (ppm)</th>
                    <th style={MODAL_STYLES.th}>Multiplicity</th>
                    <th style={MODAL_STYLES.th}>J (Hz)</th>
                  </tr>
                </thead>
                <tbody>
                  {spectrumData.peaks.map((peak, i) => (
                    <tr key={i}>
                      <td style={MODAL_STYLES.td}>{i + 1}</td>
                      <td style={MODAL_STYLES.td}>{peak.ppm.toFixed(2)}</td>
                      <td style={MODAL_STYLES.td}>{peak.multiplicity || peak.label}</td>
                      <td style={MODAL_STYLES.td}>
                        {peak.J_hz?.length ? peak.J_hz.map(j => j.toFixed(1)).join(", ") : "-"}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}
        </div>
      </div>
    </>
  );
}

export default SpectrumDetailModal;
