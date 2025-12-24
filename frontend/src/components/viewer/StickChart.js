// src/components/viewer/StickChart.js
// Composant de visualisation de spectre 1D (Canvas)

import React, { useRef, useState, useEffect, useMemo, useCallback } from "react";
import { useMeasureWidth } from "../../hooks/useMeasureWidth";
import { buildImgSrc } from "../../utils/imageUtils";
import { debugWarn, debugError } from "../../utils/debug";
import { make1DMapper } from "../../api";
import {
  COLOR_THEORY,
  COLOR_EXP,
  BG,
  AXIS,
  PLOT_HEIGHT,
} from "../../constants/colors";

/**
 * Composant Canvas pour afficher un spectre 1D (stick chart)
 */
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
    if (minX === maxX) {
      minX -= 0.01;
      maxX += 0.01;
    }
    return { minX, maxX };
  }, [theory, expPpm, showTheory, showExp, backgroundMeta]);

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

  // Loading background image
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
    img.onload = () => {
      bgImageRef.current = img;
      setBgLoaded(true);
    };
    img.onerror = () => {
      if (!retriedWithoutCb && safeSrc !== backgroundImage) {
        retriedWithoutCb = true;
        img.src = backgroundImage;
        return;
      }
      bgImageRef.current = null;
      setBgLoaded(false);
    };
    img.src = safeSrc;
  }, [backgroundImage, backgroundMeta]);

  // Rendu du canvas
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

    // Dessin de l'image de fond
    if (bgLoaded && bgImageRef.current) {
      const imgValid =
        bgImageRef.current.complete &&
        bgImageRef.current.naturalWidth > 0 &&
        bgImageRef.current.naturalHeight > 0;

      if (!imgValid) {
        debugWarn("[WARN] Image 1D not yet loaded, skip drawImage");
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
            debugError("[ERROR] drawImage 1D:", e.message);
          }
        } else {
          try {
            ctx.drawImage(bgImageRef.current, PAD, PAD, gridW, gridH);
          } catch (e) {
            debugError("[ERROR] drawImage 1D (full):", e.message);
          }
        }
        ctx.globalAlpha = 0.12;
        ctx.fillStyle = "#fff";
        ctx.fillRect(PAD, PAD, gridW, gridH);
        ctx.globalAlpha = 1.0;
      }
    }

    // Grille
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

    // Fonction de mapping ppm -> pixel
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
            const t = (xpx - mapper.axisRect.x) / mapper.axisRect.w;
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

    // Theoretical points
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

    // Experimental points
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

    // Labels des ticks
    ctx.fillStyle = AXIS;
    ctx.font = "12px ui-sans-serif, system-ui, -apple-system";
    const ticks = 6;
    for (let i = 0; i <= ticks; i++) {
      const t = i / ticks;
      const ppm = range.minX + t * (range.maxX - range.minX);
      const x = PAD + (reverse ? 1 - t : t) * gridW;
      ctx.fillText(ppm.toFixed(2), x - 14, H - 8);
    }

    // Rectangle de zoom
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
      ctx.fillText("Updating background…", PAD + 10, PAD + 22);
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
    bgLoaded,
    backgroundMeta,
    drag,
    busy,
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
    [drag, canDragNow, backgroundMeta, onZoom, onZoomPx, range, reverse]
  );

  useEffect(() => {
    const onKey = (e) => {
      if (e.key === "Escape") setDrag(null);
    };
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
          cursor: canDragNow ? (drag ? "ew-resize" : "crosshair") : "default",
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
          Loading 1D background...
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

export default StickChart;
