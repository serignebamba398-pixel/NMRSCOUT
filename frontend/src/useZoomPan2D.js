import { useRef, useEffect } from "react";

/**
 * Zoom/pan 2D (wheel + drag) on a 2D canvas.
 * @param {HTMLCanvasElement} canvas
 * @param {{f2L:number,f2R:number,f1T:number,f1B:number}} win - current ppm window
 * @param {(w:{f2_ppm_min:number,f2_ppm_max:number,f1_ppm_min:number,f1_ppm_max:number,pixelWidth:number,pixelHeight:number})=>void} onViewportChange
 */
export default function useZoomPan2D(canvas, win, onViewportChange) {
  const drag = useRef(null);

  useEffect(() => {
    if (!canvas) return;
    const el = canvas;

    const clamp = (a, b) => [Math.min(a, b), Math.max(a, b)];
    const emit = (f2L, f2R, f1T, f1B) => {
      const [x0, x1] = clamp(f2L, f2R);
      const [y0, y1] = clamp(f1T, f1B);
      onViewportChange?.({
        f2_ppm_min: x0, f2_ppm_max: x1,
        f1_ppm_min: y0, f1_ppm_max: y1,
        pixelWidth: el.clientWidth || el.width || 1,
        pixelHeight: el.clientHeight || el.height || 1
      });
    };

    const onWheel = (e) => {
      e.preventDefault();
      const rect = el.getBoundingClientRect();
      const mx = (e.clientX - rect.left) / rect.width;   // 0..1
      const my = (e.clientY - rect.top)  / rect.height;  // 0..1 (0 at the top)

      const z = e.deltaY < 0 ? 0.85 : 1.15;
      const spanX = (win.f2R - win.f2L) * z;
      const spanY = (win.f1B - win.f1T) * z;

      const c2 = win.f2L + mx * (win.f2R - win.f2L);
      const c1 = win.f1T + my * (win.f1B - win.f1T);

      emit(c2 - mx * spanX, c2 + (1 - mx) * spanX,
           c1 - my * spanY, c1 + (1 - my) * spanY);
    };

    const onDown = (e) => { drag.current = { x: e.clientX, y: e.clientY, ...win }; };
    const onUp   = () => { drag.current = null; };
    const onMove = (e) => {
      if (!drag.current) return;
      const rect = el.getBoundingClientRect();
      const dx = (e.clientX - drag.current.x) / rect.width;
      const dy = (e.clientY - drag.current.y) / rect.height;

      const spanX = drag.current.f2R - drag.current.f2L;
      const spanY = drag.current.f1B - drag.current.f1T;

      emit(drag.current.f2L - dx * spanX, drag.current.f2R - dx * spanX,
           drag.current.f1T - dy * spanY, drag.current.f1B - dy * spanY);
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
  }, [canvas, win?.f2L, win?.f2R, win?.f1T, win?.f1B, onViewportChange]);
}
