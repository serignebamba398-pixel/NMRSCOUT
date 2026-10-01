// src/hooks/useMeasureWidth.js
// Hook to measure the width of an element

import { useRef, useState, useEffect } from "react";

/**
 * Measure the width of a DOM element.
 * @param {number} initial - Initial width before the element is measured.
 * @returns {[React.RefObject, number]} - Element ref and measured width.
 */
export function useMeasureWidth(initial = 600) {
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

export default useMeasureWidth;
