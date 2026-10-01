// src/components/viewer/StickChart.js
// Canvas renderer for 1D spectra.
// Draw raw spectrum data so zooming does not require a backend request.

import React, { useRef, useState, useEffect, useMemo, useCallback } from "react";
import { useMeasureWidth } from "../../hooks/useMeasureWidth";
import {
  COLOR_THEORY,
  COLOR_EXP,
  BG,
  AXIS,
  PLOT_HEIGHT,
} from "../../constants/colors";

/**
 * Render a 1D spectrum and predicted peaks on a canvas.
 * Use raw data to redraw the visible region on zoom.
 * 
 * @param {Object} spectrumData - Raw spectrum data: {ppm: number[], intensity: number[]}.
 * @param {Object} ppmRange - Optional zoom range in ppm: {min, max}.
 */
function StickChart({
  theory = [],
  expPpm = [],
  expInt = [],
  reverse = true,
  height = PLOT_HEIGHT,
  showTheory = true,
  showExp = true,
  spectrumData = null,      // Raw spectrum data: {ppm: [], intensity: []}.
  ppmRange = null,          // Zoom range: {min, max}.
  spectrumColor = "#1e40af",
  canInteract = false,
  zoomEnabled = true,
  onZoom,
  onReset,
  busy = false,
  extraButtons = null,
  intensityScale = 1.0,     // Intensity scale factor.
  onIntensityChange = null, // Intensity adjustment callback: (delta) => void.
}) {
  const [wrapRef, width] = useMeasureWidth(700);
  const canvasRef = useRef(null);
  const [drag, setDrag] = useState(null);

  // Determine the displayed ppm range.
  const range = useMemo(() => {
    // Use the explicit zoom range first.
    if (ppmRange && ppmRange.min != null && ppmRange.max != null) {
      return { minX: Math.min(ppmRange.min, ppmRange.max), maxX: Math.max(ppmRange.min, ppmRange.max) };
    }
    // Otherwise, use the raw spectrum bounds.
    if (spectrumData?.ppm?.length > 0) {
      const minPpm = Math.min(...spectrumData.ppm);
      const maxPpm = Math.max(...spectrumData.ppm);
      return { minX: minPpm, maxX: maxPpm };
    }
    // Fall back to the predicted and experimental peak bounds.
    const all = [
      ...(showTheory ? theory : []),
      ...(showExp ? expPpm : []),
    ].filter(Number.isFinite);

    if (all.length === 0) return { minX: 0, maxX: 14 };
    let minX = Math.min(...all);
    let maxX = Math.max(...all);
    if (minX === maxX) {
      minX -= 0.01;
      maxX += 0.01;
    }
    return { minX, maxX };
  }, [theory, expPpm, showTheory, showExp, spectrumData, ppmRange]);

  const expSeries = useMemo(() => {
    const norm = (arr) => {
      if (!arr || arr.length === 0) return [];
      const m = Math.max(
        ...arr.map((v) => (Number.isFinite(v) ? v : 0)),
        1
      );
      return arr.map((v) => (Number.isFinite(v) ? v / m : 0.6));
    };
    const a = norm(expInt);
    return expPpm.map((x, i) => ({ x, a: a[i] ?? 0.6 }));
  }, [expPpm, expInt]);

  // Render the canvas.
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

    // Map ppm to pixels; NMR axes run from higher ppm on the left to lower ppm on the right.
    const denom = range.maxX - range.minX || 1;
    const ppmToX = (ppm) => {
      const t = (ppm - range.minX) / denom;
      return PAD + (reverse ? (1 - t) : t) * gridW;
    };

    // Draw the background spectrum from raw data.
    if (spectrumData?.ppm?.length > 0 && spectrumData?.intensity?.length > 0) {
      const ppmArr = spectrumData.ppm;
      const intArr = spectrumData.intensity;

      // Keep points inside the visible range.
      const visibleData = [];
      let maxInt = 0;
      for (let i = 0; i < ppmArr.length; i++) {
        const ppm = ppmArr[i];
        if (ppm >= range.minX && ppm <= range.maxX) {
          visibleData.push({ ppm, int: intArr[i] });
          if (Math.abs(intArr[i]) > maxInt) maxInt = Math.abs(intArr[i]);
        }
      }
      if (maxInt === 0) maxInt = 1;

      // Map scaled intensity to canvas y coordinates.
      const intToY = (int) => {
        const scaled = int * intensityScale;
        return PAD + gridH - (scaled / maxInt) * gridH * 0.85;
      };

      // Draw the spectrum curve.
      if (visibleData.length > 1) {
        ctx.strokeStyle = spectrumColor;
        ctx.lineWidth = 1;
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
      }
    }

    // Bottom horizontal axis
    ctx.strokeStyle = AXIS;
    ctx.beginPath();
    ctx.moveTo(PAD, PAD + gridH);
    ctx.lineTo(PAD + gridW, PAD + gridH);
    ctx.stroke();

    // Theoretical points (sticks)
    if (showTheory && theory?.length) {
      ctx.strokeStyle = COLOR_THEORY;
      ctx.lineWidth = 2.5;
      for (const ppm of theory) {
        const x = Math.round(ppmToX(ppm)) + 0.5;
        ctx.beginPath();
        ctx.moveTo(x, PAD + gridH);
        ctx.lineTo(x, PAD + gridH * 0.95);
        ctx.stroke();
      }
    }

    // Experimental points (sticks)
    if (showExp && expSeries?.length) {
      ctx.strokeStyle = COLOR_EXP;
      ctx.lineWidth = 2.5;
      for (const { x: ppm } of expSeries) {
        const xv = Math.round(ppmToX(ppm)) + 0.5;
        ctx.beginPath();
        ctx.moveTo(xv, PAD + gridH);
        ctx.lineTo(xv, PAD + gridH * 0.95);
        ctx.stroke();
      }
    }

    // PPM tick labels.
    ctx.fillStyle = AXIS;
    ctx.font = "12px ui-sans-serif, system-ui, -apple-system";
    const ticks = 6;
    for (let i = 0; i <= ticks; i++) {
      const t = i / ticks;
      const ppm = range.minX + t * (range.maxX - range.minX);
      const x = PAD + (reverse ? 1 - t : t) * gridW;
      ctx.fillText(ppm.toFixed(2), x - 14, H - 8);
    }

    // Zoom selection rectangle
    if (drag && Number.isFinite(drag.x0) && Number.isFinite(drag.x1)) {
      const x0 = Math.min(drag.x0, drag.x1);
      const x1 = Math.max(drag.x0, drag.x1);
      ctx.fillStyle = "rgba(37,99,235,0.15)";
      ctx.fillRect(x0, PAD, x1 - x0, gridH);
      ctx.strokeStyle = "rgba(37,99,235,0.65)";
      ctx.lineWidth = 2;
      ctx.strokeRect(x0 + 0.5, PAD + 0.5, x1 - x0 - 1, gridH - 1);
    }

    // Overlay busy
    if (busy) {
      ctx.fillStyle = "rgba(255,255,255,0.85)";
      ctx.fillRect(PAD, PAD, gridW, gridH);
      ctx.fillStyle = "#64748b";
      ctx.font = "13px ui-sans-serif, system-ui, -apple-system";
      ctx.fillText("Loading...", PAD + 10, PAD + 22);
    }
  }, [
    width,
    height,
    reverse,
    theory,
    expSeries,
    range,
    showExp,
    showTheory,
    spectrumData,
    spectrumColor,
    drag,
    busy,
    intensityScale,
  ]);

  const canDragNow = canInteract && zoomEnabled;

  const onMouseDown = useCallback(
    (e) => {
      if (!canDragNow || e.button !== 0) return;
      e.preventDefault();
      const rect = e.currentTarget.getBoundingClientRect();
      const x = e.clientX - rect.left;
      setDrag({ x0: x, x1: x });
    },
    [canDragNow]
  );

  const onMouseMove = useCallback(
    (e) => {
      if (!drag) return;
      e.preventDefault();
      const rect = e.currentTarget.getBoundingClientRect();
      const x = e.clientX - rect.left;
      setDrag((d) => ({ ...d, x1: x }));
    },
    [drag]
  );

  const onMouseUp = useCallback(
    (e) => {
      if (!drag) return;
      e.preventDefault();

      const PAD = 30;
      const rect = e.currentTarget.getBoundingClientRect();
      const W = rect.width;
      const gridW = W - 2 * PAD;

      const eps = 0.5;
      const x0 = Math.max(PAD + eps, Math.min(W - PAD - eps, Math.min(drag.x0, drag.x1)));
      const x1 = Math.max(PAD + eps, Math.min(W - PAD - eps, Math.max(drag.x0, drag.x1)));

      setDrag(null);

      if (!canDragNow || x1 - x0 < 4) return;

      // Convert the selection to a ppm range.
      if (typeof onZoom === "function" && range) {
        const t0 = (x0 - PAD) / gridW;
        const t1 = (x1 - PAD) / gridW;
        const rangeWidth = range.maxX - range.minX;
        const ppmA = reverse
          ? range.maxX - t0 * rangeWidth
          : range.minX + t0 * rangeWidth;
        const ppmB = reverse
          ? range.maxX - t1 * rangeWidth
          : range.minX + t1 * rangeWidth;
        onZoom(Math.min(ppmA, ppmB), Math.max(ppmA, ppmB));
      }
    },
    [drag, canDragNow, onZoom, range, reverse]
  );

  useEffect(() => {
    const onKey = (e) => {
      if (e.key === "Escape") setDrag(null);
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, []);

  // Register a non-passive wheel listener so intensity adjustments can prevent page scrolling.
  useEffect(() => {
    const cvs = canvasRef.current;
    if (!cvs || typeof onIntensityChange !== "function") return;
    
    const handleWheel = (e) => {
      e.preventDefault();
      onIntensityChange(e.deltaY);
    };
    
    cvs.addEventListener("wheel", handleWheel, { passive: false });
    return () => cvs.removeEventListener("wheel", handleWheel);
  }, [onIntensityChange]);

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
          cursor: canDragNow ? (drag ? "ew-resize" : "crosshair") : "default",
        }}
        onMouseDown={onMouseDown}
        onMouseMove={onMouseMove}
        onMouseUp={onMouseUp}
        onMouseLeave={() => drag && setDrag(null)}
      />
      {canInteract && (
        <div
          style={{
            position: "absolute",
            top: 8,
            left: 8,
            fontSize: 12,
            color: "#475569",
            display: "flex",
            alignItems: "center",
            gap: 8,
          }}
        >
          <span>Drag to zoom •</span>
          <button
            onClick={onReset}
            style={{
              padding: "2px 6px",
              borderRadius: 6,
              border: "1px solid #cbd5e1",
              background: "#fff",
              cursor: "pointer",
            }}
          >
            Reset Zoom
          </button>
          {extraButtons}
        </div>
      )}
    </div>
  );
}

export default StickChart;
