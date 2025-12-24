// src/hooks/useMeasureWidth.js
// Hook to measure the width of an element

import { useRef, useState, useEffect } from "react";

/**
 * Hook qui mesure la largeur d'un élément DOM
 * @param {number} initial - Largeur initiale par défaut
 * @returns {[React.RefObject, number]} - Ref à attacher et largeur mesurée
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
