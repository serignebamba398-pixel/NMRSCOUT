// src/utils/debug.js
// Conditional debug utilities

// Enable debug by typing in console: window.DEBUG_NMR = true
if (typeof window !== "undefined" && !window.hasOwnProperty("DEBUG_NMR")) {
  window.DEBUG_NMR = false; // false by default in production
}

/**
 * Conditional log (only if DEBUG_NMR is enabled)
 */
export const debugLog = (...args) => {
  if (typeof window !== "undefined" && window.DEBUG_NMR) {
    console.log("[NMR Debug]", ...args);
  }
};

/**
 * Conditional warning
 */
export const debugWarn = (...args) => {
  if (typeof window !== "undefined" && window.DEBUG_NMR) {
    console.warn("[NMR Debug]", ...args);
  }
};

/**
 * Conditional error
 */
export const debugError = (...args) => {
  if (typeof window !== "undefined" && window.DEBUG_NMR) {
    console.error("[NMR Debug]", ...args);
  }
};

/**
 * Installs the WebGL error filter
 */
export function installWebGLErrorFilter() {
  if (
    typeof window !== "undefined" &&
    !window.__NMR_ERROR_FILTER_INSTALLED__
  ) {
    window.__NMR_ERROR_FILTER_INSTALLED__ = true;

    const originalError = console.error;
    console.error = function (...args) {
      const msg = args.join(" ");

      // Filter known non-critical WebGL errors
      if (
        msg.includes("GL_INVALID_FRAMEBUFFER_OPERATION") ||
        msg.includes("Framebuffer is incomplete") ||
        msg.includes("Attachment has zero size")
      ) {
        if (window.DEBUG_NMR) {
          originalError.apply(console, ["[WebGL Filtered]", ...args]);
        }
        return;
      }

      // Display real errors
      originalError.apply(console, args);
    };

    debugLog("[OK] WebGL error filter installed");
  }
}
