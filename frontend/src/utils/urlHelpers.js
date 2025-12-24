/**
 * Helpers to normalize static and server URLs
 */

// Electron mode detection
const IS_ELECTRON = typeof window !== 'undefined' && window.electronAPI?.isElectron === true;

// Server base URL (without /api for static files)
const getServerBase = () => {
  if (IS_ELECTRON) {
    // In Electron, backend is on http://localhost:5000
    return "http://localhost:5000";
  }
  // In web mode, relative URLs work
  return "";
};

/**
 * Returns the API backend base URL
 * In Electron: http://localhost:5000
 * In Docker/Web: uses REACT_APP_API_URL or /api via nginx proxy
 */
export const getApiBaseUrl = () => {
  if (IS_ELECTRON) {
    return "http://localhost:5000";
  }
  // In Docker/Web mode, we can use an env variable or nginx proxy
  // Nginx redirects /api to the backend
  return process.env.REACT_APP_API_URL || "";
};

/**
 * Normalizes a static URL (outputs) for use on client-side
 * In Electron mode, prefixes with backend server URL
 * @param {string} url - Raw server URL
 * @returns {string} Normalized URL
 */
export const normalizeStaticUrl = (url) => {
  if (!url || typeof url !== "string") return url;
  
  // Already a full URL
  if (url.startsWith("http://") || url.startsWith("https://")) {
    return url;
  }

  const serverBase = getServerBase();
  
  // Remove "api/" from start if present (outputs are served directly, not via /api)
  let normalized = url;
  if (url.startsWith("/api/outputs/")) {
    normalized = url.slice(4); // → /outputs/...
  } else if (url.startsWith("api/outputs/")) {
    normalized = "/" + url.slice(4); // → /outputs/...
  } else if (!url.startsWith("/")) {
    normalized = "/" + url;
  }

  // In Electron mode, add server base
  if (serverBase) {
    return serverBase + normalized;
  }

  return normalized;
};

/**
 * Downloads a blob as a file
 * @param {Blob} blob - Data to download
 * @param {string} filename - File name
 */
export const downloadBlob = (blob, filename) => {
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
};