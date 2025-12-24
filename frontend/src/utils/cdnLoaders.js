// src/utils/cdnLoaders.js
// Loading CDN libraries (RDKit, 3Dmol)

// eslint-disable-next-line no-unused-vars
import { debugLog, debugWarn, debugError } from "./debug";

let rdkitPromise = null;
let mol3dPromise = null;

/**
 * Loads RDKit from CDN
 */
export function loadRDKitFromCDN() {
  if (typeof window === "undefined") {
    return Promise.reject(new Error("window unavailable"));
  }
  if (window.RDKit) return Promise.resolve(window.RDKit);
  if (rdkitPromise) return rdkitPromise;

  rdkitPromise = new Promise((resolve, reject) => {
    const sel = 'script[data-rdkit="RDKit_minimal"]';
    const ex = document.querySelector(sel);

    const onReady = () => {
      if (typeof window.initRDKitModule !== "function") {
        reject(new Error("initRDKitModule not found"));
        return;
      }
      window
        .initRDKitModule()
        .then((RDKit) => {
          window.RDKit = RDKit;
          resolve(RDKit);
        })
        .catch(reject);
    };

    if (ex) {
      ex.addEventListener("load", onReady, { once: true });
      ex.addEventListener(
        "error",
        () => reject(new Error("RDKit.js failed")),
        { once: true }
      );
      return;
    }

    const s = document.createElement("script");
    s.async = true;
    s.defer = true;
    s.dataset.rdkit = "RDKit_minimal";
    s.src = "https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js";
    s.onload = onReady;
    s.onerror = () => reject(new Error("RDKit.js failed"));
    document.head.appendChild(s);
  });

  return rdkitPromise;
}

/**
 * Detects problematic browsers for WebGL
 */
export function isProblematicBrowser() {
  const ua = navigator.userAgent.toLowerCase();
  const isFirefox = ua.includes("firefox");
  const isSafari = /^((?!chrome|android).)*safari/i.test(navigator.userAgent);
  return isFirefox || isSafari;
}

/**
 * Tests if WebGL is available
 */
export function isWebGLAvailable() {
  try {
    const canvas = document.createElement("canvas");
    const gl =
      canvas.getContext("webgl") || canvas.getContext("experimental-webgl");
    if (!gl) return false;

    // Basic render test
    gl.clearColor(0, 0, 0, 1);
    gl.clear(gl.COLOR_BUFFER_BIT);

    return true;
  } catch (e) {
    console.error("WebGL test failed:", e);
    return false;
  }
}

/**
 * Completely disables OffscreenCanvas
 */
function disableOffscreenCanvasCompletely() {
  try {
    if ("OffscreenCanvas" in window) {
      window._originalOffscreenCanvas = window.OffscreenCanvas;
      delete window.OffscreenCanvas;
      Object.defineProperty(window, "OffscreenCanvas", {
        value: undefined,
        writable: false,
        configurable: false,
      });
    }
    if (HTMLCanvasElement.prototype.transferControlToOffscreen) {
      HTMLCanvasElement.prototype._originalTransferControlToOffscreen =
        HTMLCanvasElement.prototype.transferControlToOffscreen;
      HTMLCanvasElement.prototype.transferControlToOffscreen = function () {
        return this;
      };
    }
    if (isProblematicBrowser() && window.createImageBitmap) {
      window._originalCreateImageBitmap = window.createImageBitmap;
      window.createImageBitmap = function () {
        throw new Error("createImageBitmap disabled");
      };
    }
    window.OffscreenCanvas = function () {
      throw new Error("OffscreenCanvas disabled");
    };
  } catch {
    // Silently ignore
  }
}

/**
 * Patch 3Dmol to disable OffscreenCanvas
 */
function patch3DmolForOffscreenCanvas() {
  if (typeof window.$3Dmol === "undefined") {
    setTimeout(patch3DmolForOffscreenCanvas, 100);
    return;
  }
  try {
    if (window.$3Dmol) {
      window.$3Dmol.global_enable_offscreen = false;
      if (window.$3Dmol.GLViewer) {
        window.$3Dmol.GLViewer.prototype._enable_offscreen = false;
      }
      const orig = window.$3Dmol.createViewer;
      if (orig) {
        window.$3Dmol.createViewer = function (el, cfg = {}) {
          const safe = {
            ...cfg,
            disableOffscreenCanvas: true,
            antialias: isProblematicBrowser() ? false : cfg.antialias !== false,
            preserveDrawingBuffer: true,
          };
          return orig.call(this, el, safe);
        };
      }
    }
  } catch {
    // Silently ignore
  }
}

/**
 * Charge 3Dmol depuis le CDN
 */
export function load3DmolFromCDN() {
  if (typeof window === "undefined") {
    return Promise.reject(new Error("window indisponible"));
  }
  if (window.$3Dmol) return Promise.resolve(window.$3Dmol);
  if (mol3dPromise) return mol3dPromise;

  disableOffscreenCanvasCompletely();

  const urls = [
    "https://3dmol.org/build/3Dmol-min.js",
    "https://cdn.jsdelivr.net/npm/3dmol/build/3Dmol-min.js",
  ];

  mol3dPromise = new Promise(async (resolve, reject) => {
    // Check if WebGL works BEFORE loading 3Dmol
    if (!isWebGLAvailable()) {
      reject(
        new Error(
          "WebGL not available in this browser. Please enable WebGL in about:config (webgl.disabled = false) or update your graphics drivers."
        )
      );
      return;
    }

    for (const u of urls) {
      try {
        await new Promise((res, rej) => {
          const ex = document.querySelector(`script[src="${u}"]`);
          if (ex) {
            ex.addEventListener("load", res, { once: true });
            ex.addEventListener("error", rej, { once: true });
            return;
          }
          const s = document.createElement("script");
          s.async = true;
          s.defer = true;
          s.src = u;
          s.onload = res;
          s.onerror = rej;
          document.head.appendChild(s);
        });
        if (window.$3Dmol) {
          patch3DmolForOffscreenCanvas();
          resolve(window.$3Dmol);
          return;
        }
      } catch {
        // Try next URL
      }
    }
    reject(new Error("Impossible de charger 3Dmol"));
  });

  return mol3dPromise;
}
