import { useEffect, useRef } from 'react';

const number = value => value == null || String(value).trim() === '' ? NaN : Number(value);
const range = (min, max) => {
  const a = number(min), b = number(max);
  return Number.isFinite(a) && Number.isFinite(b) && a !== b
    ? { min: Math.min(a, b), max: Math.max(a, b) } : null;
};

// Axis settings apply to the canvas even when no Bruker contours are loaded.
export function useSpectrumWindowSettings(options, setZoom) {
  const applied = useRef({});
  useEffect(() => {
    const updates = {};
    for (const type of ['COSY', 'TOCSY', 'HSQC']) {
      const opts = options?.[type] || {};
      const f1 = range(opts.f1_ppm_min, opts.f1_ppm_max);
      const f2 = range(opts.f2_ppm_min, opts.f2_ppm_max);
      const key = JSON.stringify({ f1, f2 });
      if (applied.current[type] === key) continue;
      applied.current[type] = key;
      if (f1 && f2) updates[type] = { f1, f2 };
    }
    if (Object.keys(updates).length) setZoom(prev => ({ ...prev, ...updates }));
  }, [options, setZoom]);
}
