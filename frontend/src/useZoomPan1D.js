import { useRef, useEffect } from "react";

/**
 * Adds zoom (wheel) + pan (drag) on a 1D canvas.
 * @param {HTMLCanvasElement} canvas
 * @param {{min:number,max:number}} view  - current ppm window
 * @param {(v:{ppm_min:number, ppm_max:number, pixelWidth:number})=>void} onViewportChange
 * @param {boolean} reverse - if visual axis decreases left→right (TopSpin-like)
 */
export default function useZoomPan1D(canvas, view, onViewportChange, reverse = true) {
  const drag = useRef(null);

  useEffect(() => {
    if (!canvas) return;
    const el = canvas;

    const clamp = (a, b) => [Math.min(a, b), Math.max(a, b)];
    const emit = (min, max) => {
      const [lo, hi] = clamp(min, max);
      onViewportChange?.({ ppm_min: lo, ppm_max: hi, pixelWidth: el.clientWidth || el.width || 1 });
    };

    const onWheel = (e) => {
      e.preventDefault();
      const rect = el.getBoundingClientRect();
      const mx = (e.clientX - rect.left) / rect.width; // 0..1
      const center = reverse ? 1 - mx : mx;
      const z = e.deltaY < 0 ? 0.85 : 1.15;            // in/out
      const span = (view.max - view.min) * z;
      const c = view.min + center * (view.max - view.min);
      const newMin = c - center * span;
      const newMax = c + (1 - center) * span;
      emit(newMin, newMax);
    };

    const onDown = (e) => { drag.current = { x: e.clientX, min: view.min, max: view.max }; };
    const onUp   = () => { drag.current = null; };
    const onMove = (e) => {
      if (!drag.current) return;
      const dx = (e.clientX - drag.current.x) / (el.clientWidth || el.width || 1);
      const span = drag.current.max - drag.current.min;
      const shift = (reverse ? -dx : dx) * span;
      emit(drag.current.min - shift, drag.current.max - shift);
    };

    el.addEventListener("wheel", onWheel, { passive: false });
    el.addEventListener("mousedown", onDown);
    window.addEventListener("mouseup", onUp);
    window.addEventListener("mousemove", onMove);
    return () => {
      el.removeEventListener("wheel", onWheel);
      el.removeEventListener("mousedown", onDown);
      window.removeEventListener("mouseup", onUp);
      window.removeEventListener("mousemove", onMove);
    };
  }, [canvas, view?.min, view?.max, onViewportChange, reverse]);
}
