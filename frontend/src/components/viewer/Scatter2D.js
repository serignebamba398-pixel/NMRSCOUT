// src/components/viewer/Scatter2D.js
// Composant de visualisation de spectre 2D (Canvas)

import React, { useRef, useState, useEffect, useMemo, useCallback } from "react";
import { useMeasureWidth } from "../../hooks/useMeasureWidth";
import { buildImgSrc } from "../../utils/imageUtils";
import { debugWarn, debugError } from "../../utils/debug";
import { AXIS, PLOT_HEIGHT } from "../../constants/colors";

/**
 * Composant Canvas pour afficher un spectre 2D (scatter plot)
 */
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
    const f1s = pts.map((p) => p[0]).filter(Number.isFinite);
    const f2s = pts.map((p) => p[1]).filter(Number.isFinite);
    let minF1 = Math.min(...f1s);
    let maxF1 = Math.max(...f1s);
    let minF2 = Math.min(...f2s);
    let maxF2 = Math.max(...f2s);
    if (minF1 === maxF1) {
      minF1 -= 0.001;
      maxF1 += 0.001;
    }
    if (minF2 === maxF2) {
      minF2 -= 0.001;
      maxF2 += 0.001;
    }
    return { f1: [minF1, maxF1], f2: [minF2, maxF2] };
  }, [theoryPts, expPts, showTheory, showExp, backgroundMeta]);

  const normInt = useMemo(() => {
    if (!showExp || !expInt?.length) return null;
    const m = Math.max(
      ...expInt.map((v) => (Number.isFinite(v) ? v : 0)),
      1
    );
    return expInt.map((v) => (Number.isFinite(v) ? v / m : 0.4));
  }, [expInt, showExp]);

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
    const PAD = 36;
    const gridW = W - 2 * PAD;
    const gridH = H - 2 * PAD;

    cvs.width = Math.round(W * DPR);
    cvs.height = Math.round(H * DPR);
    const ctx = cvs.getContext("2d");
    if (DPR !== 1) ctx.scale(DPR, DPR);

    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, W, H);

    // Dessin de l'image de fond
    if (bgLoaded && bgImageRef.current) {
      const imgValid =
        bgImageRef.current.complete &&
        bgImageRef.current.naturalWidth > 0 &&
        bgImageRef.current.naturalHeight > 0;

      if (!imgValid) {
        debugWarn("[WARN] Image 2D not yet loaded, skip drawImage");
      } else {
        const bb =
          backgroundMeta?.plot_bbox_px ||
          backgroundMeta?.axis_bbox_px ||
          backgroundMeta?.axisRect?.pixel ||
          backgroundMeta?.axisRect?.px ||
          (backgroundMeta?.axisRect &&
          ["x", "y", "w", "h"].every((k) =>
            Number.isFinite(backgroundMeta.axisRect[k])
          )
            ? backgroundMeta.axisRect
            : null);

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
            debugError("[ERROR] drawImage 2D:", e.message);
          }
        } else {
          try {
            ctx.drawImage(bgImageRef.current, PAD, PAD, gridW, gridH);
          } catch (e) {
            debugError("[ERROR] drawImage 2D (full):", e.message);
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
        const f2_left = Number(backgroundMeta.f2_left);
        const f2_right = Number(backgroundMeta.f2_right);
        const f1_top = Number(backgroundMeta.f1_top);
        const f1_bottom = Number(backgroundMeta.f1_bottom);
        const dx = f2_left - f2_right || 1e-9;
        const dy = f1_top - f1_bottom || 1e-9;
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

      const dx = maxF2 - minF2 || 1e-9;
      const dy = maxF1 - minF1 || 1e-9;
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
      ctx.beginPath();
      ctx.moveTo(gx, PAD);
      ctx.lineTo(gx, PAD + gridH);
      ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(PAD, gy);
      ctx.lineTo(PAD + gridW, gy);
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
        const a = normInt ? normInt[i] : 0.8;
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

    // Ticks F1 (à gauche)
    ctx.fillStyle = AXIS;
    ctx.font = "12px ui-sans-serif, system-ui, -apple-system";
    const ticks = 6;
    for (let i = 0; i <= ticks; i++) {
      const t = i / ticks;
      const v1 = minF1 + t * (maxF1 - minF1);
      const y1 = toY(v1);
      ctx.fillText(v1.toFixed(2), 6, y1 + 4);
    }

    // Ticks F2 (en bas)
    for (let i = 0; i <= ticks; i++) {
      const t = i / ticks;
      const v2 = minF2 + t * (maxF2 - minF2);
      const x2 = toX(v2);
      ctx.fillText(v2.toFixed(2), x2 - 14, H - 8);
    }

    // Labels
    ctx.font = "bold 13px ui-sans-serif, system-ui, -apple-system";
    ctx.fillText("F2 (ppm)", PAD + gridW - 55, H - 22);
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
    width,
    height,
    ranges,
    theoryPts,
    expPts,
    normInt,
    showTheory,
    showExp,
    maxSize,
    bgLoaded,
    backgroundMeta,
    drag,
    busy,
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

    cvs.addEventListener("wheel", handler, { passive: false });
    return () => cvs.removeEventListener("wheel", handler);
  }, [canInteract, onWheelDelta]);

  const onMouseDown = useCallback(
    (e) => {
      if (!canInteract || !zoomEnabled) return;
      const rect = e.currentTarget.getBoundingClientRect();
      setDrag({
        x0: e.clientX - rect.left,
        y0: e.clientY - rect.top,
        x1: e.clientX - rect.left,
        y1: e.clientY - rect.top,
      });
    },
    [canInteract, zoomEnabled]
  );

  const onMouseMove = useCallback(
    (e) => {
      if (!drag) return;
      const rect = e.currentTarget.getBoundingClientRect();
      setDrag((d) => ({
        ...d,
        x1: e.clientX - rect.left,
        y1: e.clientY - rect.top,
      }));
    },
    [drag]
  );

  const onMouseUp = useCallback(
    (e) => {
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

        const denomX = f2_left - f2_right || 1e-9;
        const denomY = f1_top - f1_bottom || 1e-9;

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
    },
    [drag, canInteract, zoomEnabled, backgroundMeta, onZoom2D]
  );

  return (
    <div
      ref={wrapRef}
      style={{
        width: "100%",
        position: "relative",
        userSelect: "none",
        overscrollBehavior: "contain",
      }}
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
          Loading 2D background...
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

export default Scatter2D;
