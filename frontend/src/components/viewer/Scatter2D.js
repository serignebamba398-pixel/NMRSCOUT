// src/components/viewer/Scatter2D.js
// Canvas renderer for 2D spectra.
// Render contours directly so zooming does not require a backend request.

import React, { useRef, useState, useEffect, useMemo, useCallback } from "react";
import { useMeasureWidth } from "../../hooks/useMeasureWidth";
import { AXIS, PLOT_HEIGHT } from "../../constants/colors";

const COLOR_THEORY = "#2563eb";
const BG = "#f8fafc";

/**
 * Canvas renderer for 2D spectra (COSY and HSQC).
 * Use contour data to redraw the visible region on zoom.
 * 
 * @param {Object} contourData - Contour data: {contours: [{level, paths}], f1_range, f2_range}.
 * @param {Object} zoomRange - Optional viewport: {f1: {min, max}, f2: {min, max}}.
 */
function Scatter2D({
  theoryPts = [],
  expPts = [],
  expInt = [],
  showTheory = true,
  showExp = true,
  height = PLOT_HEIGHT,
  maxSize = 10,
  contourData = null,
  zoomRange = null,
  contourColor = "#1f2937",
  canInteract = false,
  zoomEnabled = true,
  onZoom2D,
  onReset,
  busy = false,
}) {
  const [wrapRef, width] = useMeasureWidth(720);
  const canvasRef = useRef(null);
  const overlayRef = useRef(null);  // Overlay canvas for zoom rectangle
  const [drag, setDrag] = useState(null);

  const ranges = useMemo(() => {
    // Helper to check if a value is a valid finite number
    const isValid = (v) => Number.isFinite(v);
    
    // Try zoomRange first
    if (zoomRange?.f1 && zoomRange?.f2 &&
        isValid(zoomRange.f1.min) && isValid(zoomRange.f1.max) &&
        isValid(zoomRange.f2.min) && isValid(zoomRange.f2.max)) {
      return {
        f1: [Math.min(zoomRange.f1.min, zoomRange.f1.max), Math.max(zoomRange.f1.min, zoomRange.f1.max)],
        f2: [Math.min(zoomRange.f2.min, zoomRange.f2.max), Math.max(zoomRange.f2.min, zoomRange.f2.max)],
      };
    }
    
    // Try contourData ranges
    if (contourData?.f1_range?.length === 2 && contourData?.f2_range?.length === 2) {
      const f1_0 = contourData.f1_range[0];
      const f1_1 = contourData.f1_range[1];
      const f2_0 = contourData.f2_range[0];
      const f2_1 = contourData.f2_range[1];
      if (isValid(f1_0) && isValid(f1_1) && isValid(f2_0) && isValid(f2_1)) {
        return {
          f1: [Math.min(f1_0, f1_1), Math.max(f1_0, f1_1)],
          f2: [Math.min(f2_0, f2_1), Math.max(f2_0, f2_1)],
        };
      }
    }
    
    // Fallback to theory/exp points
    const pts = [];
    if (showTheory) pts.push(...theoryPts);
    if (showExp) pts.push(...expPts);
    if (pts.length === 0) return { f1: [0, 10], f2: [0, 10] };
    
    const f1s = pts.map((p) => p[0]).filter(Number.isFinite);
    const f2s = pts.map((p) => p[1]).filter(Number.isFinite);
    if (f1s.length === 0 || f2s.length === 0) return { f1: [0, 10], f2: [0, 10] };
    
    let minF1 = f1s.reduce((a, b) => a < b ? a : b, f1s[0]);
    let maxF1 = f1s.reduce((a, b) => a > b ? a : b, f1s[0]);
    let minF2 = f2s.reduce((a, b) => a < b ? a : b, f2s[0]);
    let maxF2 = f2s.reduce((a, b) => a > b ? a : b, f2s[0]);
    if (minF1 === maxF1) { minF1 -= 0.5; maxF1 += 0.5; }
    if (minF2 === maxF2) { minF2 -= 0.5; maxF2 += 0.5; }
    return { f1: [minF1, maxF1], f2: [minF2, maxF2] };
  }, [theoryPts, expPts, showTheory, showExp, contourData, zoomRange]);

  const normInt = useMemo(() => {
    if (!showExp || !expInt?.length) return null;
    const m = Math.max(...expInt.map((v) => (Number.isFinite(v) ? v : 0)), 1);
    return expInt.map((v) => (Number.isFinite(v) ? v / m : 0.4));
  }, [expInt, showExp]);

  useEffect(() => {
    const cvs = canvasRef.current;
    if (!cvs) return;

    const DPR = window.devicePixelRatio || 1;
    const W = width;
    const H = height;
    // Asymmetric padding: more space for spectrum on right and top
    const PAD_LEFT = 35;   // Space for F1 labels
    const PAD_RIGHT = 10;  // Minimal right padding
    const PAD_TOP = 8;     // Minimal top padding
    const PAD_BOTTOM = 28; // Space for F2 labels
    const gridW = W - PAD_LEFT - PAD_RIGHT;
    const gridH = H - PAD_TOP - PAD_BOTTOM;

    cvs.width = Math.round(W * DPR);
    cvs.height = Math.round(H * DPR);
    const ctx = cvs.getContext("2d");
    if (DPR !== 1) ctx.scale(DPR, DPR);

    ctx.fillStyle = BG;
    ctx.fillRect(0, 0, W, H);

    const f2Range = ranges.f2[1] - ranges.f2[0] || 1;
    const f1Range = ranges.f1[1] - ranges.f1[0] || 1;
    
    const f2ToX = (f2) => PAD_LEFT + ((ranges.f2[1] - f2) / f2Range) * gridW;
    const f1ToY = (f1) => PAD_TOP + ((f1 - ranges.f1[0]) / f1Range) * gridH;

    if (contourData?.contours?.length > 0) {
      // Use clipping to handle points outside the visible area
      ctx.save();
      ctx.beginPath();
      ctx.rect(PAD_LEFT, PAD_TOP, gridW, gridH);
      ctx.clip();
      
      ctx.strokeStyle = contourColor;
      ctx.lineWidth = 0.8;
      
      // Show ALL contours - filtering is done by backend via levelStart parameter
      const validContours = contourData.contours.filter(
        c => c.paths?.length > 0
      );
      
      // Draw all contours
      for (const contour of validContours) {
        for (const path of contour.paths || []) {
          if (path.length < 2) continue;
          ctx.beginPath();
          const [firstF2, firstF1] = path[0];
          ctx.moveTo(f2ToX(firstF2), f1ToY(firstF1));
          for (let i = 1; i < path.length; i++) {
            const [f2, f1] = path[i];
            ctx.lineTo(f2ToX(f2), f1ToY(f1));
          }
          ctx.stroke();
        }
      }
      
      ctx.restore();
    }

    // Draw axes without background grid lines.
    const nGridX = 5;
    const nGridY = 5;
    ctx.strokeStyle = AXIS;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(PAD_LEFT, PAD_TOP + gridH);
    ctx.lineTo(PAD_LEFT + gridW, PAD_TOP + gridH);
    ctx.lineTo(PAD_LEFT + gridW, PAD_TOP);
    ctx.stroke();

    // Clip predicted and experimental points to the plot area.
    ctx.save();
    ctx.beginPath();
    ctx.rect(PAD_LEFT, PAD_TOP, gridW, gridH);
    ctx.clip();

    if (showTheory && theoryPts?.length) {
      ctx.fillStyle = COLOR_THEORY;
      for (const [f1, f2] of theoryPts) {
        if (!Number.isFinite(f1) || !Number.isFinite(f2)) continue;
        const x = f2ToX(f2);
        const y = f1ToY(f1);
        ctx.beginPath();
        ctx.arc(x, y, 5, 0, 2 * Math.PI);  // Reduced from 8 to 5
        ctx.fill();
      }
    }

    if (showExp && expPts?.length) {
      for (let i = 0; i < expPts.length; i++) {
        const [f1, f2] = expPts[i];
        if (!Number.isFinite(f1) || !Number.isFinite(f2)) continue;
        const x = f2ToX(f2);
        const y = f1ToY(f1);
        const alpha = Math.max(0.6, normInt?.[i] ?? 0.8);
        const size = 3 + alpha * (maxSize - 3);  // Reduced base size from 5 to 3
        ctx.fillStyle = `rgba(220, 38, 38, ${alpha})`;
        ctx.beginPath();
        ctx.arc(x, y, size, 0, 2 * Math.PI);
        ctx.fill();
      }
    }

    ctx.restore();  // Restore the unclipped drawing context.

    ctx.fillStyle = AXIS;
    ctx.font = "10px ui-sans-serif, system-ui, -apple-system";
    ctx.textAlign = "center";
    for (let i = 0; i <= nGridX; i++) {
      const t = i / nGridX;
      const f2 = ranges.f2[1] - t * f2Range;
      const x = PAD_LEFT + t * gridW;
      ctx.fillText(f2.toFixed(1), x, H - 8);
    }
    ctx.fillText("F2 (ppm)", PAD_LEFT + gridW / 2, H - 1);

    ctx.textAlign = "right";
    for (let i = 0; i <= nGridY; i++) {
      const t = i / nGridY;
      const f1 = ranges.f1[0] + t * f1Range;
      const y = PAD_TOP + t * gridH;
      ctx.fillText(f1.toFixed(1), PAD_LEFT - 3, y + 3);
    }
    ctx.save();
    ctx.translate(8, PAD_TOP + gridH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center";
    ctx.fillText("F1 (ppm)", 0, 0);
    ctx.restore();

    // Don't draw drag rectangle here - it's handled by overlay canvas

    if (busy) {
      ctx.fillStyle = "rgba(255,255,255,0.85)";
      ctx.fillRect(PAD_LEFT, PAD_TOP, gridW, gridH);
      ctx.fillStyle = "#64748b";
      ctx.font = "13px ui-sans-serif, system-ui, -apple-system";
      ctx.textAlign = "left";
      ctx.fillText("Loading...", PAD_LEFT + 10, PAD_TOP + 22);
    }
  }, [width, height, theoryPts, expPts, normInt, ranges, showExp, showTheory, contourData, contourColor, maxSize, busy]);

  // Separate effect for drawing zoom rectangle on overlay (fast, doesn't redraw spectrum)
  useEffect(() => {
    const overlay = overlayRef.current;
    if (!overlay) return;

    const DPR = window.devicePixelRatio || 1;
    const W = width;
    const H = height;

    overlay.width = Math.round(W * DPR);
    overlay.height = Math.round(H * DPR);
    const ctx = overlay.getContext("2d");
    if (DPR !== 1) ctx.scale(DPR, DPR);

    // Clear overlay
    ctx.clearRect(0, 0, W, H);

    // Draw drag rectangle if active
    if (drag && Number.isFinite(drag.x0) && Number.isFinite(drag.x1)) {
      const x0 = Math.min(drag.x0, drag.x1);
      const y0 = Math.min(drag.y0, drag.y1);
      const w = Math.abs(drag.x1 - drag.x0);
      const h = Math.abs(drag.y1 - drag.y0);
      ctx.fillStyle = "rgba(37,99,235,0.15)";
      ctx.fillRect(x0, y0, w, h);
      ctx.strokeStyle = "rgba(37,99,235,0.65)";
      ctx.lineWidth = 2;
      ctx.strokeRect(x0 + 0.5, y0 + 0.5, w - 1, h - 1);
    }
  }, [drag, width, height]);

  const canDragNow = canInteract && zoomEnabled;

  const onMouseDown = useCallback((e) => {
    if (!canDragNow || e.button !== 0) return;
    e.preventDefault();
    const rect = e.currentTarget.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    setDrag({ x0: x, y0: y, x1: x, y1: y });
  }, [canDragNow]);

  const onMouseMove = useCallback((e) => {
    if (!drag) return;
    e.preventDefault();
    const rect = e.currentTarget.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    setDrag((d) => ({ ...d, x1: x, y1: y }));
  }, [drag]);

  const onMouseUp = useCallback((e) => {
    if (!drag) return;
    e.preventDefault();

    // Use same padding values as in useEffect
    const PAD_LEFT = 35;
    const PAD_RIGHT = 10;
    const PAD_TOP = 8;
    const PAD_BOTTOM = 28;
    const rect = e.currentTarget.getBoundingClientRect();
    const W = rect.width;
    const H = rect.height;
    const gridW = W - PAD_LEFT - PAD_RIGHT;
    const gridH = H - PAD_TOP - PAD_BOTTOM;

    const x0 = Math.max(PAD_LEFT, Math.min(drag.x0, drag.x1));
    const x1 = Math.min(W - PAD_RIGHT, Math.max(drag.x0, drag.x1));
    const y0 = Math.max(PAD_TOP, Math.min(drag.y0, drag.y1));
    const y1 = Math.min(H - PAD_BOTTOM, Math.max(drag.y0, drag.y1));

    setDrag(null);

    if (!canDragNow || (x1 - x0 < 4 && y1 - y0 < 4)) return;

    if (typeof onZoom2D === "function" && ranges) {
      const f2Range = ranges.f2[1] - ranges.f2[0];
      const f1Range = ranges.f1[1] - ranges.f1[0];

      const f2Min = ranges.f2[1] - ((x1 - PAD_LEFT) / gridW) * f2Range;
      const f2Max = ranges.f2[1] - ((x0 - PAD_LEFT) / gridW) * f2Range;
      const f1Min = ranges.f1[0] + ((y0 - PAD_TOP) / gridH) * f1Range;
      const f1Max = ranges.f1[0] + ((y1 - PAD_TOP) / gridH) * f1Range;

      onZoom2D({
        f1: { min: f1Min, max: f1Max },
        f2: { min: f2Min, max: f2Max },
      });
    }
  }, [drag, canDragNow, onZoom2D, ranges]);

  useEffect(() => {
    const onKey = (e) => {
      if (e.key === "Escape") setDrag(null);
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, []);

  return (
    <div ref={wrapRef} style={{ width: "100%", position: "relative", userSelect: "none" }}>
      <canvas
        ref={canvasRef}
        style={{
          width: "100%",
          height,
          borderRadius: 8,
          background: BG,
          boxShadow: "0 1px 4px rgba(0,0,0,.08)",
          display: "block",
        }}
      />
      {/* Overlay canvas for zoom rectangle - handles mouse events */}
      <canvas
        ref={overlayRef}
        style={{
          position: "absolute",
          top: 0,
          left: 0,
          width: "100%",
          height,
          borderRadius: 8,
          cursor: canDragNow ? "crosshair" : "default",
          pointerEvents: canDragNow ? "auto" : "none",
        }}
        onMouseDown={onMouseDown}
        onMouseMove={onMouseMove}
        onMouseUp={onMouseUp}
        onMouseLeave={() => drag && setDrag(null)}
      />
      {canInteract && (
        <div style={{ position: "absolute", top: 8, left: 48, fontSize: 12, color: "#475569", display: "flex", alignItems: "center", gap: 8, pointerEvents: "auto" }}>
          <span>Drag to zoom •</span>
          <button onClick={onReset} style={{ padding: "2px 6px", borderRadius: 6, border: "1px solid #cbd5e1", background: "#fff", cursor: "pointer" }}>
            Reset Zoom
          </button>
        </div>
      )}
    </div>
  );
}

export default Scatter2D;
