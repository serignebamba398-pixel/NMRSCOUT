// src/api.js
import axios from "axios";

// Electron mode detection via window.electronAPI (exposed by preload.js)
const IS_ELECTRON = window.electronAPI?.isElectron === true;

// Robust normalization of API base with double "/api" correction
// In Electron mode, use http://localhost:5000 (port defined in main.js line 84)
const RAW_BASE = IS_ELECTRON
  ? "http://localhost:5000/api"  // Electron mode - backend on port 5000
  : (process.env.REACT_APP_API_URL ?? "/api").trim();  // Web mode

let API_BASE_URL = (RAW_BASE === "" ? "/api" : RAW_BASE).replace(/\/+$/, "");

// Automatic fix for /api/api duplication
if (API_BASE_URL.includes("/api/api")) {
  API_BASE_URL = API_BASE_URL.replace(/\/api\/api/g, "/api");
  console.warn("[WARN] Double /api/ detected and auto-corrected");
}

// Ensure leading slash for relative URLs (except for complete URLs)
if (!API_BASE_URL) {
  API_BASE_URL = "/api";
} else if (!API_BASE_URL.startsWith("/") && !API_BASE_URL.startsWith("http")) {
  API_BASE_URL = "/" + API_BASE_URL;
}

export { API_BASE_URL };
console.log("[OK] API Base URL:", API_BASE_URL);

// Persistent key
export const DEFAULT_PERSIST_KEY = process.env.REACT_APP_PERSIST_KEY || "frontend-session";

// Axios without timeout (backend SSE streaming)
const api = axios.create({ baseURL: API_BASE_URL, timeout: 0 });


/* ============================================================================
   GENERAL UTILITIES
   ============================================================================ */
export const downloadBlob = (blob, filename) => {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
};

export const formatFileSize = (bytes) => {
  if (!Number.isFinite(bytes) || bytes <= 0) return "0 Bytes";
  const units = ["Bytes", "KB", "MB", "GB", "TB"];
  const i = Math.min(
    units.length - 1,
    Math.floor(Math.log(bytes) / Math.log(1024))
  );
  const val = bytes / Math.pow(1024, i);
  return `${val.toFixed(val >= 10 || i === 0 ? 0 : 2)} ${units[i]}`;
};

export const formatDuration = (ms) => {
  if (!Number.isFinite(ms) || ms < 0) return "0s";
  const s = Math.floor(ms / 1000);
  const m = Math.floor(s / 60);
  const h = Math.floor(m / 60);
  if (h > 0) return `${h}h ${m % 60}m ${s % 60}s`;
  if (m > 0) return `${m}m ${s % 60}s`;
  return `${s}s`;
};

export const outputUrl = (taskId, filename) =>
  `${API_BASE_URL}/outputs/${encodeURIComponent(taskId)}/${encodeURIComponent(
    filename
  )}`;

const absUrl = (u) => u ? (String(u).startsWith("http") ? u : `${API_BASE_URL}${u}`) : u;
const addCacheBust = (u) => (!u ? u : (u.includes("?") ? `${u}&v=${Date.now()}` : `${u}?v=${Date.now()}`));

const coerceTransparent = (v) =>
  typeof v === "boolean" ? (v ? 1 : 0)
  : (Number(v) === 0 || Number(v) === 1) ? Number(v)
  : undefined;

/* ============================================================================
   [NEW] RETRY MECHANISM (to handle 503 errors)
   ============================================================================ */

/**
 * Async wait
 */
const sleep = (ms) => new Promise(resolve => setTimeout(resolve, ms));

/**
 * Automatic retry with exponential backoff on 503 errors
 * @param {Function} fn - Async function to execute
 * @param {Object} options - Retry options
 * @param {number} options.retries - Max number of attempts (default: 3)
 * @param {number} options.base - Base delay in ms (default: 600)
 * @param {number} options.factor - Multiplicative factor (default: 1.8)
 * @returns {Promise} Function result or throw error
 */
export async function withRetry503(fn, { retries = 3, base = 600, factor = 1.8 } = {}) {
  let multiplier = 1;
  
  for (let attempt = 0; attempt <= retries; attempt++) {
    try {
      return await fn();
    } catch (error) {
      const status = error?.response?.status || error?.status;
      
      // If not a 503 OR it's the last attempt -> final failure
      if (status !== 503 || attempt === retries) {
        throw error;
      }
      
      // Calculate delay with jitter (random to avoid collisions)
      const jitter = 0.8 + Math.random() * 0.4; // Between 0.8 and 1.2
      const delay = Math.round(base * multiplier * jitter);
      
      console.warn(
        `[WARN] Error 503 (attempt ${attempt + 1}/${retries + 1}) - ` +
        `retrying in ${delay}ms...`
      );
      
      await sleep(delay);
      multiplier *= factor; // Progressive delay increase
    }
  }
}

// === Path normalization (Variante A) ========================================
export async function normalizeServerPath(rawPath) {
  const r = await fetch(`${API_BASE_URL}/path/normalize`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ rawPath }),
  });
  const data = await r.json().catch(() => ({}));
  if (!r.ok || !data.ok) throw new Error(data.error || `HTTP ${r.status}`);
  return data.path; // ex: "/host/c/Users/.../pdata/1"
}

// === 13C carbon modes =======================================================
// IMPORTANT: backend expects 'all' or 'nonquaternary'
export const CARBON_MODES = {
  NON_QUATERNARY: 'nonquaternary',
  ALL: 'all',
};

// Normalize free value -> 'all' or 'nonquaternary' (default)

// === HSQC fold helpers =======================================================
export const HSQC_FOLD_DEFAULTS = {
  fold: true,           // fold ON by default (backend: without --no-fold)
  fold_low: 4.8618,     // ppm (C) low by default
  fold_high: 69.8395,   // ppm (C) high by default
};

/**
 * Retrieves prediction task details with paths to latest annotation results
 */
export const getPredictionTaskWithAnnotations = async (taskId) => {
  const tasks = await getPredictionTasks();
  const task = tasks.find((t) => t.taskId === taskId);
  if (!task) throw new Error(`Task ${taskId} not found`);

  const validCsv = filterValidCsvFiles(task.csvFiles || {});
  
  // Build paths to "latest" annotation files
  const annotationResults = {};
  for (const type of ["1H", "13C", "COSY", "HSQC"]) {
    const latestPath = `/outputs/${taskId}/annotations/annotation_${type}_latest.csv`;
    annotationResults[type] = latestPath;
  }

  return {
    ...task,
    csvFiles: validCsv,
    csvUrls: Object.fromEntries(
      Object.entries(validCsv).map(([type, filename]) => [
        type,
        `/outputs/${taskId}/${filename}`,
      ])
    ),
    sdfUrls: (task.sdfFiles || []).map(
      (filename) => `/outputs/${taskId}/${filename}`
    ),
    annotationResults,
  };
};

/**
 * Lists annotation archives for a given task and type
 * (not listable via HTTP directory index -> placeholder)
 */
export const getAnnotationHistory = async (taskId, annotationType) => {
  try {
    const historyPath = `/outputs/${taskId}/annotations/history`;
    const response = await fetch(`${API_BASE_URL}${historyPath}`);
    if (!response.ok) return [];
    return [];
  } catch (e) {
    console.error(`Error retrieving history ${taskId}/${annotationType}:`, e);
    return [];
  }
};

// Normalize boolean tolerant: true/false/1/0/"true"/"false"
const coerceBool = (v) => {
  if (typeof v === "boolean") return v;
  if (v === 1 || v === "1" || String(v).toLowerCase() === "true") return true;
  if (v === 0 || v === "0" || String(v).toLowerCase() === "false") return false;
  return undefined;
};

// Quick helpers to compose options
export const withHSQCFold = (options = {}, { fold, fold_low, fold_high } = {}) => ({
  ...options,
  extractors: {
    ...(options.extractors || {}),
    hsqc: {
      ...(options.extractors?.hsqc || {}),
      ...(fold !== undefined ? { fold: !!coerceBool(fold) } : {}),
      ...(Number.isFinite(+fold_low)  ? { fold_low:  +fold_low }  : {}),
      ...(Number.isFinite(+fold_high) ? { fold_high: +fold_high } : {}),
    },
  },
});

export const with13CAllCarbons = (options = {}) => ({
  ...options,
  extractors: {
    ...(options.extractors || {}),
    c13: {
      ...(options.extractors?.c13 || {}),
      carbon_mode: CARBON_MODES.ALL,
    },
  },
});


/* ============================================================================
   PARAMS HELPERS
   ============================================================================ */
const clampInt = (v, min, max) => {
  const n = Number(v);
  if (!Number.isFinite(n)) return null;
  const i = Math.trunc(n);
  if (i < min) return min;
  if (i > max) return max;
  return i;
};

const sanitizeParameters = (params = {}) => {
  const out = {};
  for (const [k, v] of Object.entries(params)) {
    if (v === null || v === undefined) continue;
    if (typeof v === "number" && !Number.isFinite(v)) continue;
    out[k] = v;
  }
  // Avoid nb_cores "null" on server side + bounds [1..3]
  if ("nb_cores" in out) {
    const c = clampInt(out.nb_cores, 1, 3);
    if (c === null) delete out.nb_cores;
    else out.nb_cores = c;
  }
  return out;
};

/* ============================================================================
   [NEW] HELPER: NETWORK ERROR FORMATTING
   ============================================================================ */
const formatNetworkError = (error, context = "") => {
  const status = error?.response?.status;
  const url = error?.config?.url || "";
  
  if (status === 503) {
    return `[WARN] Service temporarily unavailable (503)\n` +
           `The backend cannot execute R/Python scripts.\n` +
           `Check that containers are properly started.\n` +
           `Endpoint: ${url}`;
  }
  
  if (status === 404) {
    if (url.includes("/tasks/")) {
      return `[WARN] Task not found (404)\n` +
             `This task no longer exists or was never created.\n` +
             `URL: ${url}`;
    }
    return `[WARN] Resource not found (404): ${url}`;
  }
  
  if (status === 502 || status === 504) {
    return `[WARN] Backend not responding (${status})\n` +
           `The server may be overloaded or restarting.\n` +
           `Please retry in a few moments.`;
  }
  
  if (error?.code === "ECONNABORTED" || error?.message?.includes("timeout")) {
    return `⏱️ Timeout: Request took too long\n` +
           `Check your connection or retry.`;
  }
  
  const msg = error?.response?.data?.error || 
              error?.response?.data?.message || 
              error?.message || 
              "Unknown network error";
  
  return `[ERROR] Error ${context}: ${msg}`;
};

/* ============================================================================
   FILES (upload / listing / delete)
   ============================================================================ */
export const uploadFile = async (file) => {
  const formData = new FormData();
  formData.append("file", file);
  const { data } = await api.post("/upload", formData, {
    headers: { "Content-Type": "multipart/form-data" },
  });
  return data; // { filename }
};

export const getFiles = async () => {
  const { data } = await api.get("/files");
  return data; // existing backend
};

export const deleteFile = async (filename) => {
  const { data } = await api.delete(
    `/files/${encodeURIComponent(filename)}`
  );
  return data;
};

export const downloadFile = async (taskId, filename) => {
  const { data } = await api.get(
    `/download/${encodeURIComponent(taskId)}/${encodeURIComponent(
      filename
    )}`,
    { responseType: "blob" }
  );
  return data;
};


/* === NEW =====================================================================
   LIST OF ALREADY UPLOADED FILES (to select without re-uploading)
   ========================================================================== */
const _inferKindByName = (name = "") => {
  const n = String(name).toLowerCase();
  if (n.endsWith(".csv")) return "csv";
  if (n.endsWith(".sdf")) return "sdf";
  if (n.endsWith(".json")) return "json";
  return "other";
};
const _filterByKind = (arr, kind) => {
  if (!kind) return arr;
  return arr.filter((f) => (f.kind || _inferKindByName(f.filename)) === kind);
};

const _normalizeUploadsPayload = (raw) => {
  const files = Array.isArray(raw?.files)
    ? raw.files
    : Array.isArray(raw)
    ? raw
    : [];
  return files
    .map((f) => {
      const filename = f.filename || f.name || f.file || "";
      const size = Number.isFinite(f.size)
        ? f.size
        : Number(f.length) || undefined;
      const mtime =
        f.mtime ||
        f.uploadedAt || // (present server-side)
        f.updatedAt ||
        f.modifiedAt ||
        (f.mtimeMs ? new Date(f.mtimeMs) : undefined) ||
        f.date ||
        f.createdAt;
      const kind = f.kind || _inferKindByName(filename);
      return { filename, size, mtime, kind };
    })
    .filter((f) => f.filename);
};

/**
 * Lists files present server-side (uploads).
 * @param {'csv'|'matrix'|'sdf'|'json'|undefined} kind  Optional filter
 */
export async function listUploads(kind) {
  // 1) Dedicated endpoint if available
  try {
    const r = await fetch(
      `${API_BASE_URL}/uploads${
        kind ? `?kind=${encodeURIComponent(kind)}` : ""
      }`
    );
    if (r.ok) {
      const json = await r.json();
      const items = _normalizeUploadsPayload(json);
      return _filterByKind(items, kind).sort(
        (a, b) => new Date(b.mtime || 0) - new Date(a.mtime || 0)
      );
    }
    // [NEW] LOG IF 404 (missing endpoint, normal fallback)
    if (r.status === 404) {
    }
  } catch (e) {
    console.warn("[WARN] /uploads failed, fallback to /files:", e.message);
  }

  // 2) Fallback to /files (already present)
  try {
    const data = await getFiles();
    const items = _normalizeUploadsPayload(data);
    const kindFix = kind === "matrix" ? "csv" : kind;
    return _filterByKind(items, kindFix).sort(
      (a, b) => new Date(b.mtime || 0) - new Date(a.mtime || 0)
    );
  } catch (e) {
    console.error("[ERROR] listUploads error:", e);
    return [];
  }
}

/* ============================================================================
   SSE HELPERS (POST streaming recommended)
   ============================================================================ */
const normalizeProgress = (msg) => {
  let p = Number(msg?.percent);
  if (!Number.isFinite(p)) {
    const cur = Number(msg?.current);
    const tot = Number(msg?.total);
    if (Number.isFinite(cur) && Number.isFinite(tot) && tot > 0) {
      p = (cur / tot) * 100;
    }
  }
  if (Number.isFinite(p)) {
    return { ...msg, percent: Math.max(0, Math.min(100, p)) };
  }
  return msg;
};

// POST (streaming fetch) — cancellable via .cancel() and .cancelTask()
// IMPORTANT: this must NOT be `async`, otherwise the returned Promise would be wrapped
// and we'd lose custom methods/properties (cancel/cancelTask/getTaskId).
const createPostSSEConnection = (
  url,
  body,
  headers = {},
  { onProgress, onLog, onError, onComplete, onTaskId } = {}
) => {
  const controller = new AbortController();
  const { signal } = controller;

  let currentTaskId = null;

  const p = new Promise(async (resolve, reject) => {
    try {
      const response = await fetch(url, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Accept: "text/event-stream",
          ...headers,
        },
        body: JSON.stringify(body),
        signal,
      });

      if (!response.ok) {
        const errorText = await response.text().catch(() => "");
        
        // [NEW] IMPROVED HTTP ERROR HANDLING
        let userMessage = "";
        
        if (response.status === 503) {
          userMessage = "[WARN] Service temporarily unavailable\n" +
                       "The backend cannot start R/Python processes.\n" +
                       "Possible causes:\n" +
                       "• R or Python not installed in the container\n" +
                       "• Insufficient system resources\n" +
                       "• Backend is restarting\n\n" +
                       `Details: ${errorText || response.statusText}`;
        } else if (response.status === 404) {
          userMessage = `[ERROR] Endpoint not found (404)\n` +
                       `The URL ${url} does not exist on the server.\n` +
                       `Check Nginx configuration and backend routes.`;
        } else if (response.status === 502 || response.status === 504) {
          userMessage = `[WARN] Backend not responding (${response.status})\n` +
                       "The server is overloaded, restarted or unreachable.\n" +
                       "Wait a few moments and retry.";
        } else {
          userMessage = `HTTP ${response.status}: ${errorText || response.statusText}`;
        }
        
        onError && onError(userMessage);
        reject(new Error(userMessage));
        return;
      }

      const reader = response.body.getReader();
      const decoder = new TextDecoder();
      let buffer = "";
      let finished = false;

      const pump = async () => {
        const { done, value } = await reader.read();
        if (done) {
          if (!finished) {
            const msg = "[WARN] Connection ended without 'completed' message";
            onError && onError(msg);
            reject(new Error(msg));
          }
          return;
        }

        buffer += decoder.decode(value, { stream: true });
        const lines = buffer.split("\n");
        buffer = lines.pop() || "";

        for (const raw of lines) {
          const line = raw.trim();
          if (!line || line.startsWith(":")) continue;
          if (!line.startsWith("data:")) continue;
          const payload = line.slice(5).trim();

          try {
            const data = JSON.parse(payload);
            switch (data.type) {
              case "taskId":
                currentTaskId = data.taskId;
                onLog && onLog(`[OK] Task started: ${currentTaskId}`);
                onTaskId && onTaskId(currentTaskId);
                // Expose taskId immediately
                p._taskId = currentTaskId;
                break;
              case "log":
                onLog && onLog(data.message);
                break;
              case "progress":
                onProgress && onProgress(normalizeProgress(data));
                break;
              case "error":
                onError && onError(data.message || "Unknown error");
                break;
              case "result":
                onLog && onLog(`[OK] Result received`);
                break;
              case "completed":
                finished = true;
                onComplete && onComplete({ ...data, taskId: currentTaskId });
                if (data.success) resolve({ ...data, taskId: currentTaskId });
                else
                  reject(new Error(data.error || "Annotation/prediction failed"));
                return;
              default:
                onLog && onLog(`[INFO] Message: ${JSON.stringify(data)}`);
            }
          } catch (e) {
            onLog && onLog(`[WARN] Parsing error: ${e.message}`);
          }
        }
        return pump();
      };

      try {
        await pump();
      } finally {
        try {
          reader.releaseLock?.();
        } catch {}
      }
    } catch (err) {
      if (err.name === "AbortError") {
        onLog && onLog("[STOP] Connection cancelled by user");
        reject(new Error("Cancelled by user"));
      } else {
        const formattedError = formatNetworkError(err, "SSE");
        onError && onError(formattedError);
        reject(err);
      }
    }
  });

  // Expose taskId via a property
  p._taskId = null;

  // Client-side cancel method (abort fetch)
  p.cancel = () => {
    try {
      controller.abort();
    } catch {}
  };

  // Server-side cancel method (calls /tasks/:taskId/cancel)
  p.cancelTask = async () => {
    const taskId = p._taskId;
    if (!taskId) {
      console.warn("[WARN] No taskId available for server-side cancellation");
      return;
    }
    try {
      await cancelTask(taskId);
    } catch (e) {
      console.error("[ERROR] Server cancellation error:", e);
    }
  };

  // Method to get taskId
  p.getTaskId = () => p._taskId;

  return p;
};

// Wrap a cancelable SSE promise with a 404 fallback without losing cancel methods.
const _wrapCancelableWith404Fallback = (primary, fallbackFactory, onLog) => {
  let inner = primary;
  const wrapped = new Promise((resolve, reject) => {
    primary
      .then(resolve)
      .catch(async (e) => {
        const msg = String(e?.message || "");
        const is404 = /\b404\b/i.test(msg) || /endpoint\s+not\s+found/i.test(msg);
        if (!is404) {
          reject(e);
          return;
        }

        try {
          onLog && onLog("Route /predict/start unavailable (404) — trying REST fallback...");
          const fallback = await fallbackFactory();
          inner = fallback;
          resolve(await fallback);
        } catch (err) {
          reject(err);
        }
      });
  });

  // Delegate cancel methods to the current inner promise.
  wrapped.cancel = () => {
    try {
      inner?.cancel?.();
    } catch {}
  };
  wrapped.cancelTask = async () => {
    try {
      await inner?.cancelTask?.();
    } catch {}
  };
  wrapped.getTaskId = () => {
    try {
      return inner?.getTaskId?.() || inner?._taskId || null;
    } catch {
      return null;
    }
  };

  return wrapped;
};

/* ============================================================================
   TASK CANCELLATION (NEW)
   ============================================================================ */
/**
 * Cancels a task server-side by killing the R/Python process
 * @param {string} taskId - ID of the task to cancel
 * @returns {Promise<{message: string, taskId: string, canceled: boolean}>}
 */
export const cancelTask = async (taskId) => {
  try {
    const { data } = await api.post(`/tasks/${encodeURIComponent(taskId)}/cancel`);
    return data;
  } catch (e) {
    if (e?.response?.status === 404) {
      return { 
        message: 'Task not found or already completed', 
        taskId,
        canceled: false 
      };
    }
    throw e;
  }
};

/* ============================================================================
   R ANNOTATIONS — MAIN API
   ============================================================================ */
export const runAnnotationPost = (
  tableauFile,
  matrixFile,
  parameters,
  annotationType = "1H",
  onProgress,
  onLog,
  onError,
  onComplete
) => {
  const url = `${API_BASE_URL}/run-annotation`;
  const cleanParams = sanitizeParameters(parameters);
  const body = {
    tableauFile,
    matrixFile,
    parameters: cleanParams,
    annotationType,
  };
  return createPostSSEConnection(url, body, {}, { onProgress, onLog, onError, onComplete });
};

export const runAnnotationValidated = async (
  tableauFile,
  matrixFile,
  parameters = {},
  annotationType = "1H",
  onProgress,
  onLog,
  onError,
  onComplete
) => {
  if (!tableauFile || !matrixFile) {
    throw new Error("Table and matrix files are required");
  }
  if (!["1H", "13C", "COSY", "HSQC"].includes(annotationType)) {
    throw new Error(
      "Unsupported annotation type. Use: 1H, 13C, COSY, HSQC"
    );
  }
  const cleanParams = sanitizeParameters(parameters);
  onLog && onLog(`Starting ${annotationType} annotation...`);
  onLog && onLog(`Table file: ${tableauFile}`);
  onLog && onLog(`Matrix file: ${matrixFile}`);
  onLog && onLog(`Parameters: ${JSON.stringify(cleanParams, null, 2)}`);

  return runAnnotationPost(
    tableauFile,
    matrixFile,
    cleanParams,
    annotationType,
    onProgress,
    onLog,
    onError,
    onComplete
  );
};

export const runAnnotation1H = (
  tableauFile,
  matrixFile,
  parameters = {},
  onProgress,
  onLog,
  onError,
  onComplete
) => {
  const defaults = {
    tol_correspondance: 0.1,
    tol_match: 0.1,
    correlation_threshold: 0.9,
    coverage_threshold_group: 0,
    coverage_threshold_subchain: 0,
    ratio_min: 0.01,
    ransy_cv_threshold: 0.05,
    coverage_threshold_ransy: 0,
    max_missing_ransy: 50,
  };
  return runAnnotationValidated(
    tableauFile,
    matrixFile,
    { ...defaults, ...parameters },
    "1H",
    onProgress,
    onLog,
    onError,
    onComplete
  );
};

export const runAnnotation13C = (
  tableauFile,
  matrixFile,
  parameters = {},
  onProgress,
  onLog,
  onError,
  onComplete
) => {
  const defaults = {
    tol_correspondance: 2,
    tol_match: 2,
    correlation_threshold: 0.9,
    coverage_threshold_group: 0,
    coverage_threshold_subchain: 0,
    ratio_min: 0.1,
    ransy_cv_threshold: 0.1,
    coverage_threshold_ransy: 0,
    max_missing_ransy: 50,
  };
  return runAnnotationValidated(
    tableauFile,
    matrixFile,
    { ...defaults, ...parameters },
    "13C",
    onProgress,
    onLog,
    onError,
    onComplete
  );
};

export const runAnnotationCOSY = (
  tableauFile,
  matrixFile,
  parameters = {},
  onProgress,
  onLog,
  onError,
  onComplete
) => {
  const defaults = {
    tol_correspondance_F1: 0.2,
    tol_correspondance_F2: 0.2,
    tol_match_F1: 0.2,
    tol_match_F2: 0.2,
    correlation_threshold: 0.8,
    coverage_threshold_group: 0,
    coverage_threshold_subchain: 0,
    ratio_min: 0.01,
    ransy_cv_threshold: 0.15,
    coverage_threshold_ransy: 0,
    max_missing_ransy: 50,
  };
  return runAnnotationValidated(
    tableauFile,
    matrixFile,
    { ...defaults, ...parameters },
    "COSY",
    onProgress,
    onLog,
    onError,
    onComplete
  );
};

export const runAnnotationHSQC = (
  tableauFile,
  matrixFile,
  parameters = {},
  onProgress,
  onLog,
  onError,
  onComplete
) => {
  const defaults = {
    tol_correspondance_F1: 4,
    tol_correspondance_F2: 0.2,
    tol_match_F1: 4,
    tol_match_F2: 0.2,
    correlation_threshold: 0.8,
    coverage_threshold_group: 0,
    coverage_threshold_subchain: 0,
    ratio_min: 0.1,
    ransy_cv_threshold: 0.15,
    coverage_threshold_ransy: 0,
    max_missing_ransy: 50,
  };
  return runAnnotationValidated(
    tableauFile,
    matrixFile,
    { ...defaults, ...parameters },
    "HSQC",
    onProgress,
    onLog,
    onError,
    onComplete
  );
};

/* ============================================================================
   TASKS / DIAGNOSTIC / HEALTH
   ============================================================================ */
export const getTasks = async () => {
  try {
    const direct = await api.get("/tasks");
    return direct.data;
  } catch {
    const { data } = await api.get("/diagnostic");
    return data?.tasks || [];
  }
};

export const stopTask = async (taskId) => {
  try {
    const { data } = await api.post(
      `/tasks/${encodeURIComponent(taskId)}/stop`
    );
    return data;
  } catch (e) {
    if (e?.response?.status === 404) {
      return { message: "Stop not supported on this server", id: taskId };
    }
    throw e;
  }
};

export const stopPredictionTask = stopTask;

export const healthCheck = async () => {
  const { data } = await api.get("/health");
  return data;
};

export const getDiagnostic = async () => {
  const { data } = await api.get("/diagnostic");
  return data;
};

export const getMetrics = async (format = "json") => {
  try {
    const { data } = await api.get("/metrics", { params: { format } });
    return data;
  } catch {
    const health = await healthCheck();
    return {
      runningTasks: (await getTasks()).length,
      timestamp: health.timestamp,
      uptime: health.uptime,
      memory: health.memory,
      node: health.node,
      rscript: health.rscript,
      fallback: "health",
    };
  }
};

export const cleanup = async () => {
  const { data } = await api.post("/cleanup");
  return data;
};

/* ============================================================================
   ANNOTATION UTILITIES
   ============================================================================ */
export const validateAnnotationFiles = async (tableauFile, matrixFile) => {
  try {
    const files = await getFiles();
    const names = files.map((f) => f.name || f.filename).filter(Boolean);
    const errors = [];
    if (!names.includes(tableauFile))
      errors.push(`Table file not found: ${tableauFile}`);
    if (!names.includes(matrixFile))
      errors.push(`Matrix file not found: ${matrixFile}`);
    return { valid: errors.length === 0, errors };
  } catch (error) {
    return {
      valid: false,
      errors: [`Error during validation: ${error.message}`],
    };
  }
};

/* === NEW =====================================================================
   Strict detection & short labels for prediction CSV
   ========================================================================== */
const CSV_PATTERNS = {
  "1H": /(^|\/)prediction_1H_.*\.csv$/i,
  "13C": /(^|\/)prediction_13C_.*\.csv$/i,
  COSY: /(^|\/)prediction_COSY_.*\.csv$/i,
  HSQC: /(^|\/)prediction_HSQC_.*\.csv$/i,
};

export const isValidCsvForType = (filename, type) => {
  if (!filename || !type || !CSV_PATTERNS[type]) return false;
  return CSV_PATTERNS[type].test(String(filename));
};

export const shortCsvLabel = (filename) => {
  if (!filename) return "";
  if (CSV_PATTERNS["1H"].test(filename)) return "prediction_1H";
  if (CSV_PATTERNS["13C"].test(filename)) return "prediction_13C";
  if (CSV_PATTERNS.COSY.test(filename)) return "prediction_COSY";
  if (CSV_PATTERNS.HSQC.test(filename)) return "prediction_HSQC";
  return filename; // fallback
};

const filterValidCsvFiles = (csvFiles = {}) => {
  const out = {};
  for (const type of ["1H", "13C", "COSY", "HSQC"]) {
    const f = csvFiles[type];
    if (isValidCsvForType(f, type)) out[type] = f;
  }
  return out;
};

export const getAnnotationTypes = () => [
  {
    key: "1H",
    name: "1H NMR",
    description: "1D proton NMR annotation",
    expectedOutput: "1H_With_Ransy_Julien.csv",
    tolerances: { correspondance: 0.1, match: 0.1 },
    dimensions: 1,
  },
  {
    key: "13C",
    name: "13C NMR",
    description: "1D carbon-13 NMR annotation",
    expectedOutput: "13C_With_Ransy_and_ECFP4.csv",
    tolerances: { correspondance: 2, match: 2 },
    dimensions: 1,
  },
  {
    key: "COSY",
    name: "COSY 2D",
    description: "2D COSY NMR annotation (1H-1H)",
    expectedOutput: "COSY_With_Ransy_and_ECFP4.csv",
    tolerances: { F1: 0.2, F2: 0.2 },
    dimensions: 2,
  },
  {
    key: "HSQC",
    name: "HSQC 2D",
    description: "2D HSQC NMR annotation (13C-1H)",
    expectedOutput: "HSQC_With_Ransy_and_ECFP4.csv",
    tolerances: { F1: 4, F2: 0.2 },
    dimensions: 2,
  },
];

export const getDefaultParameters = (annotationType) =>
  ({
    "1H": {
      tol_correspondance: 0.1,
      tol_match: 0.1,
      correlation_threshold: 0.9,
      ratio_min: 0.01,
      ransy_cv_threshold: 0.05,
      coverage_threshold_group: 0,
      coverage_threshold_subchain: 0,
      coverage_threshold_ransy: 0,
      max_missing_ransy: 50,
    },
    "13C": {
      tol_correspondance: 2,
      tol_match: 2,
      correlation_threshold: 0.9,
      ratio_min: 0.1,
      ransy_cv_threshold: 0.1,
      coverage_threshold_group: 0,
      coverage_threshold_subchain: 0,
      coverage_threshold_ransy: 0,
      max_missing_ransy: 50,
    },
    COSY: {
      tol_correspondance_F1: 0.2,
      tol_correspondance_F2: 0.2,
      tol_match_F1: 0.2,
      tol_match_F2: 0.2,
      correlation_threshold: 0.8,
      ratio_min: 0.01,
      ransy_cv_threshold: 0.15,
      coverage_threshold_group: 0,
      coverage_threshold_subchain: 0,
      coverage_threshold_ransy: 0,
      max_missing_ransy: 50,
    },
    HSQC: {
      tol_correspondance_F1: 4,
      tol_correspondance_F2: 0.2,
      tol_match_F1: 4,
      tol_match_F2: 0.2,
      correlation_threshold: 0.8,
      coverage_threshold_group: 0,
      coverage_threshold_subchain: 0,
      ratio_min: 0.1,
      ransy_cv_threshold: 0.15,
      coverage_threshold_ransy: 0,
      max_missing_ransy: 50,
    },
  }[annotationType] || {});

/* ============================================================================
   PREDICTIONS (SMILES/CSV → SDF + tables 1H/13C/COSY/HSQC)
   ============================================================================ */
const DEFAULT_PREDICT_OPTIONS = { run_extractors: true };

const sanitizePredictOptions = (opts = {}) => {
  // Shallow clone and filter invalid values
  const out = {};
  for (const [k, v] of Object.entries(opts || {})) {
    if (v === null || v === undefined) continue;
    if (typeof v === "number" && !Number.isFinite(v)) continue;
    out[k] = v;
  }

  // Generic bounds/typing
  if ("nb_cores" in out) {
    const c = clampInt(out.nb_cores, 1, 3);
    if (c === null) delete out.nb_cores;
    else out.nb_cores = c;
  }
  if ("min_bonds" in out)
    out.min_bonds = Math.max(1, Math.min(6, Math.trunc(out.min_bonds)));
  if ("max_bonds" in out)
    out.max_bonds = Math.max(1, Math.min(6, Math.trunc(out.max_bonds)));

  // Ensure extractors
  out.extractors = { ...(out.extractors || {}) };
  out.extractors.c13  = { ...(out.extractors.c13  || {}) };
  out.extractors.hsqc = { ...(out.extractors.hsqc || {}) };

  // =========================
  // 13C carbon_mode
  // =========================
  // Known top-level variants
  let rawCarbonMode =
    out.carbon_mode_13c ?? out.carbon_mode ?? out.carbonMode ?? out?.extractors?.c13?.carbon_mode;

  // Nested camelCase variants: options.c13.carbonMode / carbon_mode
  if (rawCarbonMode == null) {
    rawCarbonMode = out?.c13?.carbonMode ?? out?.c13?.carbon_mode ?? out?.c13?.carbonmode;
  }

  if (rawCarbonMode != null) {
    const s = String(rawCarbonMode).toLowerCase().trim();
    const cm = (['all','tout','tous','all_carbons','full'].some(p => s.startsWith(p)))
      ? 'all'
      : 'nonquaternary';
    out.extractors.c13.carbon_mode = cm;
    // Cleanup used aliases
    delete out.carbon_mode_13c;
    delete out.carbon_mode;
    delete out.carbonMode;
    if (out.c13) {
      delete out.c13.carbonMode;
      delete out.c13.carbon_mode;
      delete out.c13.carbonmode;
    }
  }

  // =========================
  // HSQC fold / fold_low / fold_high
  // =========================
  // Top-level variants already sometimes used
  let hsqcFold     = out?.extractors?.hsqc?.fold     ?? out.hsqc_fold     ?? out.hsqcFold;
  let hsqcFoldLow  = out?.extractors?.hsqc?.fold_low ?? out.hsqc_fold_low ?? out.hsqcFoldLow;
  let hsqcFoldHigh = out?.extractors?.hsqc?.fold_high?? out.hsqc_fold_high?? out.hsqcFoldHigh;

  // Nested camelCase variants (the common problematic case)
  if (out.hsqc && (hsqcFold == null || hsqcFoldLow == null || hsqcFoldHigh == null)) {
    if (hsqcFold == null && 'fold' in out.hsqc) {
      hsqcFold = !!out.hsqc.fold;
    }
    if (hsqcFoldLow == null && 'foldLow' in out.hsqc) {
      const v = Number(out.hsqc.foldLow);
      if (Number.isFinite(v)) hsqcFoldLow = v;
    }
    if (hsqcFoldHigh == null && 'foldHigh' in out.hsqc) {
      const v = Number(out.hsqc.foldHigh);
      if (Number.isFinite(v)) hsqcFoldHigh = v;
    }
  }

  const t = (v) => (v === true || v === false) ? v : String(v).toLowerCase() === 'true';

  if (hsqcFold !== undefined) out.extractors.hsqc.fold = t(hsqcFold);
  if (hsqcFoldLow  !== undefined && Number.isFinite(+hsqcFoldLow)) {
    out.extractors.hsqc.fold_low = +hsqcFoldLow;
  }
  if (hsqcFoldHigh !== undefined && Number.isFinite(+hsqcFoldHigh)) {
    out.extractors.hsqc.fold_high = +hsqcFoldHigh;
  }

  // Cleanup aliases
  delete out.hsqc_fold;
  delete out.hsqcFold;
  delete out.hsqc_fold_low;
  delete out.hsqcFoldLow;
  delete out.hsqc_fold_high;
  delete out.hsqcFoldHigh;
  if (out.hsqc) {
    delete out.hsqc.fold;
    delete out.hsqc.foldLow;
    delete out.hsqc.foldHigh;
  }

  return out;
};



export const runPredictionFromSmiles = (
  smiles,
  { molId = "mol", methods = ["1H", "13C", "COSY", "HSQC"], methodOptions = {}, global = {} } = {},
  onProgress,
  onLog,
  onError,
  onComplete
) => {
  const body = {
    mode: "SMILES",
    smiles,
    methods,
    options: sanitizePredictOptions({
      ...DEFAULT_PREDICT_OPTIONS,
      ...global,
      mol_id: molId,
      extractors: {
        ...(methodOptions || {}),
        ...(global?.extractors || {}),
      }
    }),
  };

  const primary = createPostSSEConnection(
    `${API_BASE_URL}/predict/start`,
    body,
    {},
    { onProgress, onLog, onError, onComplete }
  );

  const fallbackFactory = async () => {
    const candidates = ["predict", "prediction/start", "prediction"];
    for (const path of candidates) {
      try {
        const { data } = await api.post(path, body);
        const normalized = {
          ...(data || {}),
          success: !!data?.success,
          taskId: data?.taskId || data?.id || null,
          outputs: data?.outputs || data?.result || {},
        };
        onComplete && onComplete(normalized);

        const p = Promise.resolve(normalized);
        p._taskId = normalized.taskId || null;
        p.cancel = () => {};
        p.cancelTask = async () => {};
        p.getTaskId = () => p._taskId;
        return p;
      } catch (err) {
        const st = err?.response?.status;
        const m = String(err?.message || "");
        if (st === 404 || /404|not\s*found/i.test(m)) {
          onLog && onLog(`Route ${path} unavailable (404) — trying next...`);
          continue;
        }
        throw err;
      }
    }
    throw new Error("No compatible /predict route was found.");
  };

  return _wrapCancelableWith404Fallback(primary, fallbackFactory, onLog);
};

export const runPredictionFromCSV = (
  csvFile,
  { idColumn, smilesColumn, methods = ["1H", "13C", "COSY", "HSQC"], methodOptions = {}, global = {} } = {},
  onProgress,
  onLog,
  onError,
  onComplete
) => {
  const body = {
    mode: "CSV",
    source: "csv",
    csvFile,
    methods,
    options: sanitizePredictOptions({ 
      ...DEFAULT_PREDICT_OPTIONS, 
      ...global,
      extractors: {
        ...(methodOptions || {}),
        ...(global?.extractors || {}),
      }
    }),
  };
  if (idColumn && String(idColumn).trim())
    body.id_column = String(idColumn).trim();
  if (smilesColumn && String(smilesColumn).trim())
    body.smiles_column = String(smilesColumn).trim();

  const primary = createPostSSEConnection(
    `${API_BASE_URL}/predict/start`,
    body,
    {},
    { onProgress, onLog, onError, onComplete }
  );

  const fallbackFactory = async () => {
    const candidates = ["predict", "prediction/start", "prediction"];
    for (const path of candidates) {
      try {
        const { data } = await api.post(path, body);
        const normalized = {
          ...(data || {}),
          success: !!data?.success,
          taskId: data?.taskId || data?.id || null,
          outputs: data?.outputs || data?.result || {},
        };
        onComplete && onComplete(normalized);

        const p = Promise.resolve(normalized);
        p._taskId = normalized.taskId || null;
        p.cancel = () => {};
        p.cancelTask = async () => {};
        p.getTaskId = () => p._taskId;
        return p;
      } catch (err) {
        const st = err?.response?.status;
        const m = String(err?.message || "");
        if (st === 404 || /404|not\s*found/i.test(m)) {
          onLog && onLog(`Route ${path} unavailable (404) — trying next...`);
          continue;
        }
        throw err;
      }
    }
    throw new Error("No compatible /predict route was found.");
  };

  return _wrapCancelableWith404Fallback(primary, fallbackFactory, onLog);
};

/**
 * NOTE: events/result endpoints not provided by current server -> legacy helpers kept
 * but return explicit error if called.
 */
export const connectPredictionEvents = () => {
  throw new Error("connectPredictionEvents not supported: server-side SSE endpoint unavailable");
};

export const getPredictionResult = async () => {
  throw new Error("getPredictionResult not supported: server-side REST endpoint unavailable");
};

/* ============================================================================
   INTEGRATED WORKFLOW: PRED TASK LIST + ANNOTATION FROM PRED CSV
   ============================================================================ */
export const getPredictionTasks = async () => {
  try {
    const res = await fetch(`${API_BASE_URL}/prediction-tasks`);
    if (res.ok) return await res.json();
  } catch {}
  try {
    const { data } = await api.get(`/predict/tasks`);
    return data;
  } catch {
    return [];
  }
};

export const runAnnotationFromPrediction = (
  predictionTaskId,
  matrixFile,
  annotationType,
  parameters,
  onProgress,
  onLog,
  onError,
  onComplete
) => {
  const url = `${API_BASE_URL}/run-annotation-from-prediction`;
  const body = {
    predictionTaskId,
    matrixFile,
    annotationType,
    parameters: sanitizeParameters(parameters || {}),
  };
  return createPostSSEConnection(url, body, {}, { onProgress, onLog, onError, onComplete });
};

export const validatePredictionTask = async (
  taskId,
  requiredTypes = ["1H", "13C", "COSY", "HSQC"]
) => {
  try {
    const tasks = await getPredictionTasks();
    const task = tasks.find((t) => t.taskId === taskId);

    if (!task) {
      return { valid: false, errors: [`Task ${taskId} not found`] };
    }

    const filtered = filterValidCsvFiles(task.csvFiles || {});
    const errors = [];

    requiredTypes.forEach((t) => {
      if (!filtered[t]) {
        errors.push(
          `Missing or invalid CSV for ${t} (must contain "prediction_${t}_")`
        );
      }
    });

    if (!task.canRunAnnotations) {
      errors.push("Task marked as unusable for annotations");
    }

    return {
      valid: errors.length === 0,
      errors,
      task: { ...task, csvFiles: filtered },
      availableTypes: Object.keys(filtered),
    };
  } catch (error) {
    return {
      valid: false,
      errors: [`Error during validation: ${error.message}`],
    };
  }
};

export const getPredictionTaskDetails = async (taskId) => {
  const tasks = await getPredictionTasks();
  const task = tasks.find((t) => t.taskId === taskId);
  if (!task) throw new Error(`Task ${taskId} not found`);

  const validCsv = filterValidCsvFiles(task.csvFiles || {});
  return {
    ...task,
    csvFiles: validCsv,
    csvUrls: Object.fromEntries(
      Object.entries(validCsv).map(([type, filename]) => [
        type,
        `/outputs/${taskId}/${filename}`,
      ])
    ),
    sdfUrls: (task.sdfFiles || []).map(
      (filename) => `/outputs/${taskId}/${filename}`
    ),
  };
};

export const downloadPredictionCsv = async (taskId, type) => {
  const details = await getPredictionTaskDetails(taskId);
  const csvUrl = details.csvUrls?.[type];

  if (!csvUrl) {
    throw new Error(`CSV ${type} not available for task ${taskId}`);
  }

  const response = await fetch(`${API_BASE_URL}${csvUrl}`);
  if (!response.ok) throw new Error(`HTTP ${response.status}`);

  const blob = await response.blob();
  const filename = `prediction_${type}_${taskId}.csv`;

  downloadBlob(blob, filename);
  return { success: true, filename };
};

export const runIntegrated1H = (
  predictionTaskId,
  matrixFile,
  parameters = {},
  onProgress,
  onLog,
  onError,
  onComplete
) => {
  const defaults = {
    tol_correspondance: 0.1,
    tol_match: 0.1,
    correlation_threshold: 0.9,
    coverage_threshold_group: 0,
    coverage_threshold_subchain: 0,
    ratio_min: 0.01,
    ransy_cv_threshold: 0.05,
    coverage_threshold_ransy: 0,
    max_missing_ransy: 50,
  };
  return runAnnotationFromPrediction(
    predictionTaskId,
    matrixFile,
    "1H",
    { ...defaults, ...parameters },
    onProgress,
    onLog,
    onError,
    onComplete
  );
};

export const runIntegrated13C = (
  predictionTaskId,
  matrixFile,
  parameters = {},
  onProgress,
  onLog,
  onError,
  onComplete
) => {
  const defaults = {
    tol_correspondance: 2,
    tol_match: 2,
    correlation_threshold: 0.9,
    coverage_threshold_group: 0,
    coverage_threshold_subchain: 0,
    ratio_min: 0.1,
    ransy_cv_threshold: 0.1,
    coverage_threshold_ransy: 0,
    max_missing_ransy: 50,
  };
  return runAnnotationFromPrediction(
    predictionTaskId,
    matrixFile,
    "13C",
    { ...defaults, ...parameters },
    onProgress,
    onLog,
    onError,
    onComplete
  );
};

export const runIntegratedCOSY = (
  predictionTaskId,
  matrixFile,
  parameters = {},
  onProgress,
  onLog,
  onError,
  onComplete
) => {
  const defaults = {
    tol_correspondance_F1: 0.2,
    tol_correspondance_F2: 0.2,
    tol_match_F1: 0.2,
    tol_match_F2: 0.2,
    correlation_threshold: 0.8,
    coverage_threshold_group: 0,
    coverage_threshold_subchain: 0,
    ratio_min: 0.01,
    ransy_cv_threshold: 0.15,
    coverage_threshold_ransy: 0,
    max_missing_ransy: 50,
  };
  return runAnnotationFromPrediction(
    predictionTaskId,
    matrixFile,
    "COSY",
    { ...defaults, ...parameters },
    onProgress,
    onLog,
    onError,
    onComplete
  );
};

export const runIntegratedHSQC = (
  predictionTaskId,
  matrixFile,
  parameters = {},
  onProgress,
  onLog,
  onError,
  onComplete
) => {
  const defaults = {
    tol_correspondance_F1: 4,
    tol_correspondance_F2: 0.2,
    tol_match_F1: 4,
    tol_match_F2: 0.2,
    correlation_threshold: 0.8,
    coverage_threshold_group: 0,
    coverage_threshold_subchain: 0,
    ratio_min: 0.1,
    ransy_cv_threshold: 0.15,
    coverage_threshold_ransy: 0,
    max_missing_ransy: 50,
  };
  return runAnnotationFromPrediction(
    predictionTaskId,
    matrixFile,
    "HSQC",
    { ...defaults, ...parameters },
    onProgress,
    onLog,
    onError,
    onComplete
  );
};

export const checkPredictionAnnotationCompatibility = (task, annotationType) => {
  if (!task || !annotationType) {
    return { compatible: false, reason: "Missing parameters" };
  }
  const csvFile = task.csvFiles?.[annotationType];
  if (!csvFile) {
    return {
      compatible: false,
      reason: `CSV ${annotationType} not available for task ${task.taskId}`,
    };
  }
  if (!isValidCsvForType(csvFile, annotationType)) {
    return {
      compatible: false,
      reason: `Invalid CSV file format for ${annotationType}`,
    };
  }
  return { compatible: true, csvFile };
};

export const getPredictionTasksSummary = async () => {
  const tasks = await getPredictionTasks();

  const byType = { "1H": 0, "13C": 0, COSY: 0, HSQC: 0 };
  const usableForAnnotations = tasks.filter((t) => {
    const valid = filterValidCsvFiles(t.csvFiles || {});
    Object.keys(valid).forEach((k) => (byType[k] = (byType[k] || 0) + 1));
    return t.canRunAnnotations && Object.keys(valid).length > 0;
  }).length;

  return {
    total: tasks.length,
    usableForAnnotations,
    byType,
    recent: tasks
      .slice()
      .sort((a, b) => new Date(b.createdAt) - new Date(a.createdAt))
      .slice(0, 5)
      .map((t) => {
        const valid = filterValidCsvFiles(t.csvFiles || {});
        return {
          taskId: t.taskId,
          createdAt: t.createdAt,
          availableTypes: Object.keys(valid),
          canRunAnnotations: t.canRunAnnotations,
        };
      }),
  };
};

export const runPredictionFromSmilesWithAnnotations = (
  smiles,
  options = {},
  onProgress,
  onLog,
  onError,
  onComplete
) => {
  // IMPORTANT: pass options in "global" (not "options")
  const enhancedGlobal = { ...options, run_extractors: true };

  return runPredictionFromSmiles(
    smiles,
    { global: enhancedGlobal },
    onProgress,
    (msg) => onLog && onLog(msg),
    (err) => onError && onError(err),
    (result) => {
      if (result?.success && onComplete) {
        const enhancedResult = {
          ...result,
          annotationReady: true,
          availableAnnotationTypes: ["1H", "13C", "COSY", "HSQC"],
        };
        onComplete(enhancedResult);
      } else {
        onComplete && onComplete(result);
      }
    }
  );
};

export const runPredictionFromCSVWithAnnotations = (
  csvFile,
  options = {},
  onProgress,
  onLog,
  onError,
  onComplete
) => {
  // IMPORTANT: pass options in "global" (not "options")
  const enhancedGlobal = { ...options, run_extractors: true };

  return runPredictionFromCSV(
    csvFile,
    { global: enhancedGlobal },
    onProgress,
    (msg) => onLog && onLog(msg),
    (err) => onError && onError(err),
    (result) => {
      if (result?.success && onComplete) {
        const enhancedResult = {
          ...result,
          annotationReady: true,
          availableAnnotationTypes: ["1H", "13C", "COSY", "HSQC"],
        };
        onComplete(enhancedResult);
      } else {
        onComplete && onComplete(result);
      }
    }
  );
};

/* ============================================================================
   >>> NEW: EXPORT BACKGROUNDS (PNG + JSON) FROM BRUKER PDATA
   ============================================================================ */

const withPersist = (body, { persistKey = DEFAULT_PERSIST_KEY, overwrite = true } = {}) => ({
  ...body,
  persistKey: persistKey || DEFAULT_PERSIST_KEY,
  overwrite: !!overwrite,
});

// === Cancellable POST JSON (fetch + AbortController) ============================
export function postJSONCancelable(path, body) {
  const controller = new AbortController();
  const signal = controller.signal;

  const p = (async () => {
    const res = await fetch(`${API_BASE_URL}${path}`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
      signal,
    });
    const data = await res.json().catch(() => ({}));
    if (!res.ok || data?.ok === false) {
      throw new Error(data?.error || `HTTP ${res.status}`);
    }
    return data;
  })();

  // expose cancel()
  p.cancel = () => {
    try { controller.abort(); } catch {}
  };

  return p;
}

const _normalizeBgMeta = (meta = {}, method) => {
  const out = { ...meta, method };
  if (method === "1H" || method === "13C") {
    if (meta.ppm_min != null) out.ppm_min = Number(meta.ppm_min);
    if (meta.ppm_max != null) out.ppm_max = Number(meta.ppm_max);
    if (out.ppm_min != null && out.ppm_max != null) {
      const a = out.ppm_min, b = out.ppm_max;
      out.ppm_left = Math.max(a, b);
      out.ppm_right = Math.min(a, b);
    }
  } else {
    ["f2_ppm_min","f2_ppm_max","f1_ppm_min","f1_ppm_max"].forEach(k=>{
      if (meta[k] != null) out[k] = Number(meta[k]);
    });
    if (out.f2_ppm_min != null && out.f2_ppm_max != null) {
      out.f2_left  = Math.max(out.f2_ppm_min, out.f2_ppm_max);
      out.f2_right = Math.min(out.f2_ppm_min, out.f2_ppm_max);
    }
    if (out.f1_ppm_min != null && out.f1_ppm_max != null) {
      out.f1_top    = Math.max(out.f1_ppm_min, out.f1_ppm_max);
      out.f1_bottom = Math.min(out.f1_ppm_min, out.f1_ppm_max);
    }
  }
  if (!out.fig_px) {
    const w = Number(out.fig_w || out.width_px || out.width || 0);
    const h = Number(out.fig_h || out.height_px || out.height || 0);
    if (w && h) out.fig_px = { w, h };
  }
  if (out.level_start == null) out.level_start = 0.005;
  return out;
};

export const exportSpectralBackgrounds = async (payload, { persistKey = DEFAULT_PERSIST_KEY, overwrite = true } = {}) => {
  const { data } = await api.post("/spectral-bg/export", withPersist(payload, { persistKey, overwrite }));
  const results = { taskId: data.taskId, workDir: data.workDir, outputs: {} };
  const outs = data?.outputs || {};
  for (const [method, obj] of Object.entries(outs)) {
    results.outputs[method] = {
      ...obj,
      pngUrl: addCacheBust(absUrl(obj?.pngUrl)),
      jsonUrl: addCacheBust(absUrl(obj?.jsonUrl)),
      meta: _normalizeBgMeta(obj?.meta || {}, method),
    };
  }
  return results;
};

export const exportBg1H = async (pdataDir1H, opts = {}, extra = {}) =>
  exportSpectralBackgrounds(
    {
      pdata: { "1H": pdataDir1H },
      methods: ["1H"],
      options: { oneH: opts },
    },
    extra
  );

export const exportBg13C = async (pdataDir13C, opts = {}, extra = {}) =>
  exportSpectralBackgrounds(
    {
      pdata: { "13C": pdataDir13C },
      methods: ["13C"],
      options: { c13: opts },
    },
    extra
  );

export const exportBgCOSY = async (pdataDirCOSY, opts = {}, extra = {}) =>
  exportSpectralBackgrounds(
    {
      pdata: { COSY: pdataDirCOSY },
      methods: ["COSY"],
      options: { cosy: { level_start: 0.005, ...opts } },
    },
    extra
  );

export const exportBgHSQC = async (pdataDirHSQC, opts = {}, extra = {}) =>
  exportSpectralBackgrounds(
    {
      pdata: { HSQC: pdataDirHSQC },
      methods: ["HSQC"],
      options: { hsqc: { level_start: 0.005, ...opts } },
    },
    extra
  );

/* ============================================================================
   [NEW] GROUPED EXPORTS WITH AUTOMATIC RETRY
   ============================================================================ */

/**
 * Grouped export of all spectral backgrounds with automatic retry on 503
 * 
 * @param {Object} pdataPaths - Spectrum paths { "1H": "/path/...", "13C": "...", "COSY": "...", "HSQC": "..." }
 * @param {Object} options - Options for each method { oneH: {...}, c13: {...}, cosy: {...}, hsqc: {...} }
 * @param {Object} extra - Additional options { persistKey, overwrite }
 * @returns {Promise<Object>} Results with taskId and outputs
 * 
 * @example
 * const results = await exportAllSpectralBackgroundsWithRetry({
 *   "1H": "/host/.../pdata/1",
 *   "13C": "/host/.../pdata/1",
 *   "COSY": "/host/.../pdata/1",
 *   "HSQC": "/host/.../pdata/1"
 * }, {
 *   oneH: { ppm_min: 0, ppm_max: 10 },
 *   c13: { ppm_min: 0, ppm_max: 200 },
 *   cosy: { level_start: 0.005 },
 *   hsqc: { level_start: 0.005 }
 * });
 */
export async function exportAllSpectralBackgroundsWithRetry(
  pdataPaths = {},
  options = {},
  extra = {}
) {
  // Filter methods that have a valid path
  const methods = Object.keys(pdataPaths).filter(
    method => pdataPaths[method] && typeof pdataPaths[method] === 'string'
  );
  
  if (methods.length === 0) {
    throw new Error('No valid spectrum path provided');
  }
  
  console.log(`[INFO] Grouped export of ${methods.length} spectra: ${methods.join(', ')}`);
  
  // Appel unique avec retry automatique
  return withRetry503(
    () => exportSpectralBackgrounds(
      {
        pdata: pdataPaths,
        methods,
        options: {
          // Options 1H
          ...(pdataPaths["1H"] ? { oneH: options.oneH || {} } : {}),
          
          // Options 13C
          ...(pdataPaths["13C"] ? { c13: options.c13 || {} } : {}),
          
          // Options COSY
          ...(pdataPaths["COSY"] ? { 
            cosy: { 
              level_start: 0.005, 
              ...(options.cosy || {}) 
            } 
          } : {}),
          
          // Options HSQC
          ...(pdataPaths["HSQC"] ? { 
            hsqc: { 
              level_start: 0.005, 
              ...(options.hsqc || {}) 
            } 
          } : {}),
        },
      },
      extra
    ),
    { retries: 3, base: 600 } // 3 attempts, starting at 600ms
  );
}

/**
 * Simplified grouped export (without built-in automatic retry)
 * Use exportAllSpectralBackgroundsWithRetry for more robustness
 */
export async function exportAllSpectralBackgrounds(
  pdataPaths = {},
  options = {},
  extra = {}
) {
  const methods = Object.keys(pdataPaths).filter(
    method => pdataPaths[method] && typeof pdataPaths[method] === 'string'
  );
  
  if (methods.length === 0) {
    throw new Error('No valid spectrum path provided');
  }
  
  return exportSpectralBackgrounds(
    {
      pdata: pdataPaths,
      methods,
      options: {
        ...(pdataPaths["1H"] ? { oneH: options.oneH || {} } : {}),
        ...(pdataPaths["13C"] ? { c13: options.c13 || {} } : {}),
        ...(pdataPaths["COSY"] ? { cosy: { level_start: 0.005, ...(options.cosy || {}) } } : {}),
        ...(pdataPaths["HSQC"] ? { hsqc: { level_start: 0.005, ...(options.hsqc || {}) } } : {}),
      },
    },
    extra
  );
}

export const make1DMapper = (meta, widthPx, heightPx) => {
  const a = Number(meta?.ppm_min), b = Number(meta?.ppm_max);
  const left  = Number.isFinite(Number(meta?.ppm_left))  ? Number(meta.ppm_left)  : Math.max(a, b);
  const right = Number.isFinite(Number(meta?.ppm_right)) ? Number(meta.ppm_right) : Math.min(a, b);
  const figW = Number(meta?.fig_px?.w) || Math.max(1, Number(widthPx || 0));
  const figH = Number(meta?.fig_px?.h) || Math.max(1, Number(heightPx || 0));
  const bbox = meta?.plot_bbox_px || { x: 0, y: 0, w: figW, h: figH };
  const sx = Math.max(1e-6, Number(widthPx  || figW) / figW);
  const sy = Math.max(1e-6, Number(heightPx || figH) / figH);
  const axisRect = { x: bbox.x * sx, y: bbox.y * sy, w: bbox.w * sx, h: bbox.h * sy };
  const clamp01 = (t) => (t < 0 ? 0 : t > 1 ? 1 : t);
  const frac = (ppm) => (left - Number(ppm)) / (left - right);
  const mapX = (ppm) => axisRect.x + clamp01(frac(ppm)) * axisRect.w;
  mapX.axisRect = axisRect;
  return mapX;
};

export const make2DMapper = (meta, widthPx, heightPx, { invertF1 = false } = {}) => {
  const hasF2 =
    Number.isFinite(Number(meta?.f2_ppm_min)) &&
    Number.isFinite(Number(meta?.f2_ppm_max));
  const hasF1 =
    Number.isFinite(Number(meta?.f1_ppm_min)) &&
    Number.isFinite(Number(meta?.f1_ppm_max));

  const f2_left = hasF2
    ? Math.max(Number(meta.f2_ppm_min), Number(meta.f2_ppm_max))
    : Number(meta?.f2_left ?? meta?.f2_ppm_left);
  const f2_right = hasF2
    ? Math.min(Number(meta.f2_ppm_min), Number(meta.f2_ppm_max))
    : Number(meta?.f2_right ?? meta?.f2_ppm_right);

  const f1_top = hasF1
    ? Math.max(Number(meta.f1_ppm_min), Number(meta.f1_ppm_max))
    : Number(meta?.f1_top ?? meta?.f1_ppm_top ?? meta?.f1_up ?? meta?.f1_ppm_up);
  const f1_bottom = hasF1
    ? Math.min(Number(meta.f1_ppm_min), Number(meta.f1_ppm_max))
    : Number(meta?.f1_bottom ?? meta?.f1_ppm_bottom ?? meta?.f1_down ?? meta?.f1_ppm_down);

  const figW = Number(meta?.fig_px?.w) || Math.max(1, Number(widthPx || 0));
  const figH = Number(meta?.fig_px?.h) || Math.max(1, Number(heightPx || 0));
  const bbox = meta?.plot_bbox_px || { x: 0, y: 0, w: figW, h: figH };

  const sx = Math.max(1e-6, Number(widthPx || figW) / figW);
  const sy = Math.max(1e-6, Number(heightPx || figH) / figH);

  const axisX0 = bbox.x * sx, axisY0 = bbox.y * sy;
  const axisW  = bbox.w * sx, axisH  = bbox.h * sy;

  const clamp01 = (t) => (t < 0 ? 0 : t > 1 ? 1 : t);

  const tx = (ppm2) => (f2_left - Number(ppm2)) / (f2_left - f2_right);
  const tyNMR  = (ppm1) => (Number(ppm1) - f1_top) / (f1_bottom - f1_top);
  const tyMath = (ppm1) => (Number(ppm1) - f1_bottom) / (f1_top - f1_bottom);

  const ty = invertF1 ? tyMath : tyNMR;

  const mapX = (ppm2) => axisX0 + clamp01(tx(ppm2)) * axisW;
  const mapY = (ppm1) => axisY0 + clamp01(ty(ppm1)) * axisH;

  return { mapX, mapY, axisRect: { x: axisX0, y: axisY0, w: axisW, h: axisH } };
};

/* ============================================================================
   >>> NEW: INTERACTIVE — zoom / relevel / targeted PNG+JSON background (re)generation
   ============================================================================ */

const _clip = (v, a, b) => Math.max(a, Math.min(b, v));

const _wrapInteractiveResult = (data, method) => ({
  ...data,
  pngUrl: addCacheBust(absUrl(data?.pngUrl)),
  jsonUrl: addCacheBust(absUrl(data?.jsonUrl)),
  meta: _normalizeBgMeta(data?.meta || {}, method),
});

export const spectralBgInteractive = async (method, pdataDir, opts = {}) => {
  const base = { method, pdata_dir: pdataDir };
  if (opts.taskId) base.taskId = String(opts.taskId);

  // Image size
  if (opts.size) {
    const { width_px, height_px, dpi } = opts.size;
    base.size = {
      ...(Number.isFinite(+width_px) ? { width_px: +width_px } : {}),
      ...(Number.isFinite(+height_px) ? { height_px: +height_px } : {}),
      ...(Number.isFinite(+dpi) ? { dpi: +dpi } : {}),
    };
  }

  // Windows/options for 1D/2D
  if (method === "1H" || method === "13C") {
    if (opts.window) {
      const { ppm_min, ppm_max } = opts.window;
      const w = {};
      if (Number.isFinite(+ppm_min)) w.ppm_min = +ppm_min;
      if (Number.isFinite(+ppm_max)) w.ppm_max = +ppm_max;
      if (Object.keys(w).length) base.window = w;
    }
    if (opts.line) {
      const { color, width } = opts.line;
      base.line = {
        ...(color ? { color: String(color) } : {}),
        ...(Number.isFinite(+width) ? { width: +width } : {}),
      };
    }
  } else {
    if (opts.window2d) {
      const { f2_ppm_min, f2_ppm_max, f1_ppm_min, f1_ppm_max } = opts.window2d;
      base.window2d = {
        ...(Number.isFinite(+f2_ppm_min) ? { f2_ppm_min: +f2_ppm_min } : {}),
        ...(Number.isFinite(+f2_ppm_max) ? { f2_ppm_max: +f2_ppm_max } : {}),
        ...(Number.isFinite(+f1_ppm_min) ? { f1_ppm_min: +f1_ppm_min } : {}),
        ...(Number.isFinite(+f1_ppm_max) ? { f1_ppm_max: +f1_ppm_max } : {}),
      };
    }
    const { level_start = 0.005, n_levels, mode, cmap, transparent } = opts.contour || {};
    const t = coerceTransparent(transparent);
    base.contour = {
      ...(Number.isFinite(+level_start) ? { level_start: +level_start } : { level_start: 0.005 }),
      ...(Number.isFinite(+n_levels) ? { n_levels: +n_levels } : {}),
      ...(mode ? { mode: String(mode) } : {}),
      ...(cmap ? { cmap: String(cmap) } : {}),
      ...(t !== undefined ? { transparent: t } : {}),
    };
  }

  const res = await fetch(`${API_BASE_URL}/spectral-bg/interactive`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(withPersist(base, {
      persistKey: opts.taskId || DEFAULT_PERSIST_KEY,
      overwrite: true,
    })),
  });

  const data = await res.json().catch(() => ({}));
  if (!res.ok || data?.ok === false) {
    throw new Error(data?.error || `HTTP ${res.status}`);
  }

  return {
    ...data,
    pngUrl: addCacheBust(absUrl(data?.pngUrl)),
    jsonUrl: addCacheBust(absUrl(data?.jsonUrl)),
    meta: _normalizeBgMeta(data?.meta || {}, method),
  };
};

export const spectralBgZoom = async (method, pdataDir, meta, rectPx, opts = {}) => {
  const base = {
    method,
    pdata_dir: pdataDir,
    meta,
    select: { x: +rectPx.x, y: +rectPx.y, w: +rectPx.w, h: +rectPx.h },
  };
  
  if (opts.taskId) base.taskId = String(opts.taskId);

  if (opts.size) base.size = opts.size;

  const ds = opts.displaySize || {};
  const dW = Number.isFinite(+opts.displayWidth)
    ? +opts.displayWidth
    : Number.isFinite(+ds.width)
    ? +ds.width
    : null;
  const dH = Number.isFinite(+opts.displayHeight)
    ? +opts.displayHeight
    : Number.isFinite(+ds.height)
    ? +ds.height
    : null;
  if (dW || dH) {
    base.display_size = {
      ...(dW ? { width_px: dW } : {}),
      ...(dH ? { height_px: dH } : {}),
    };
  }

  if ((method === "1H" || method === "13C") && opts.line) {
    base.line = {
      ...(opts.line.color ? { color: String(opts.line.color) } : {}),
      ...(Number.isFinite(+opts.line.width) ? { width: +opts.line.width } : {}),
    };
  }
  if (method === "COSY" || method === "HSQC") {
    const { level_start = 0.005, transparent, ...rest } = opts.contour || {};
    const t = coerceTransparent(transparent);
    base.contour = {
      level_start,
      ...rest,
      ...(t !== undefined ? { transparent: t } : {}),
    };
  }

  const { data } = await api.post(
    "/spectral-bg/zoom", 
    withPersist(base, { 
      persistKey: opts.taskId || DEFAULT_PERSIST_KEY, 
      overwrite: true 
    })
  );
  
  if (!data?.ok) {
    throw new Error(data?.error || "zoom failed");
  }
  return _wrapInteractiveResult(data, method);
};

export function compute1DWindowFromRect(meta, rectPx, displayWidth, displayHeight) {
  if (!meta?.fig_px?.w || !meta?.fig_px?.h || !meta?.plot_bbox_px) {
    throw new Error("meta.fig_px/plot_bbox_px missing for 1D zoom");
  }
  const figW = +meta.fig_px.w, figH = +meta.fig_px.h;
  const bbox = meta.plot_bbox_px;
  const dispW = Number(displayWidth)  || figW;
  const dispH = Number(displayHeight) || figH;
  const sx = figW / dispW;
  const sy = figH / dispH;
  const sel = {
    x: (+rectPx.x) * sx,
    y: (+rectPx.y) * sy,
    w: (+rectPx.w) * sx,
    h: (+rectPx.h) * sy,
  };
  const x0 = _clip(sel.x, bbox.x, bbox.x + bbox.w);
  const x1 = _clip(sel.x + sel.w, bbox.x, bbox.x + bbox.w);
  const wEff = Math.max(1, Math.abs(x1 - x0));
  const t0 = (x0 - bbox.x) / bbox.w;
  const t1 = (x1 - bbox.x) / bbox.w;
  let left  = Number(meta.ppm_left);
  let right = Number(meta.ppm_right);
  if (!Number.isFinite(left) || !Number.isFinite(right) || left <= right) {
    const a = Number(meta.ppm_min);
    const b = Number(meta.ppm_max);
    if (Number.isFinite(a) && Number.isFinite(b) && a !== b) {
      left  = Math.max(a, b);
      right = Math.min(a, b);
    } else {
      throw new Error("meta.ppm_left/right invalid and no usable ppm_min/max");
    }
  }
  const pA = left + t0 * (right - left);
  const pB = left + t1 * (right - left);
  let ppm_min = Math.min(pA, pB);
  let ppm_max = Math.max(pA, pB);
  if (wEff < 2 && Number.isFinite(meta.plot_bbox_px?.w)) {
    const mid = (ppm_min + ppm_max) / 2;
    const pad = Math.max(0.001, (left - right) * 0.002);
    ppm_min = mid - pad;
    ppm_max = mid + pad;
  }
  return { ppm_min, ppm_max, width_px_native: wEff };
}

const cleanMetaForBackend = (meta) => {
  if (!meta || typeof meta !== 'object') return {};
  
  const clean = {};
  
  for (const [key, value] of Object.entries(meta)) {
    // Ignore null, undefined, NaN
    if (value === null || value === undefined) continue;
    if (typeof value === 'number' && !Number.isFinite(value)) continue;
    
    // Clean nested objects (e.g., fig_px, plot_bbox_px)
    if (typeof value === 'object' && !Array.isArray(value)) {
      const nested = {};
      
      for (const [k, v] of Object.entries(value)) {
        // [NOTE] STRICT FILTER: ignore null, undefined, NaN, empty string
        if (v === null || v === undefined || v === '') continue;
        if (typeof v === 'number' && !Number.isFinite(v)) continue;
        
        nested[k] = v;
      }
      
      // Only keep object if it contains at least 1 valid value
      if (Object.keys(nested).length > 0) {
        clean[key] = nested;
      }
    } 
    // Keep valid primitive values
    else {
      clean[key] = value;
    }
  }
  
  return clean;
};
export const spectralBgRelevel = async (method, pdataDir, meta, opts = {}) => {
  if (!["COSY", "HSQC"].includes(method)) {
    throw new Error("relevel only supports COSY/HSQC");
  }
  const base = {
    method,
    pdata_dir: pdataDir,
    meta: cleanMetaForBackend(meta),  
  };
  if (opts.taskId) base.taskId = String(opts.taskId);
  if (opts.size) base.size = opts.size;
  base.level_start = Number.isFinite(+opts.levelStart) ? +opts.levelStart : 0.005;
  if (opts.nLevels != null) base.n_levels = +opts.nLevels;
  if (opts.mode) base.mode = String(opts.mode);
  if (opts.cmap) base.cmap = String(opts.cmap);
  if (opts.transparent != null) base.transparent = +!!opts.transparent;

  const { data } = await api.post(
    "/spectral-bg/relevel", 
    withPersist(base, { 
      persistKey: opts.taskId || DEFAULT_PERSIST_KEY, 
      overwrite: true 
    })
  );
  
  if (!data?.ok) throw new Error(data?.error || "relevel failed");
  return _wrapInteractiveResult(data, method);
};

// Wrap to normalize like existing helpers
const _wrapCancelable = (raw, method) => ({
  ...raw,
  pngUrl: addCacheBust(absUrl(raw?.pngUrl)),
  jsonUrl: addCacheBust(absUrl(raw?.jsonUrl)),
  meta: _normalizeBgMeta(raw?.meta || {}, method),
});

/**
 * Cancelable relevel (COSY/HSQC) — identical to spectralBgRelevel but cancelable
 */
export function spectralBgRelevelCancelable(method, pdataDir, meta, opts = {}) {
  if (!["COSY", "HSQC"].includes(method)) {
    throw new Error("cancelable relevel only supports COSY/HSQC");
  }
  
  // [NOTE] CLEAN meta
  const cleanedMeta = cleanMetaForBackend(meta);
  
  // [NOTE] CLEAN opts.size (remove null/undefined)
  let cleanedSize = null;
  if (opts.size && typeof opts.size === 'object') {
    const temp = {};
    
    if (Number.isFinite(+opts.size.width_px) && +opts.size.width_px > 0) {
      temp.width_px = +opts.size.width_px;
    }
    if (Number.isFinite(+opts.size.height_px) && +opts.size.height_px > 0) {
      temp.height_px = +opts.size.height_px;
    }
    if (Number.isFinite(+opts.size.dpi) && +opts.size.dpi > 0) {
      temp.dpi = +opts.size.dpi;
    }
    
    // Only keep size if it contains at least 1 valid value
    if (Object.keys(temp).length > 0) {
      cleanedSize = temp;
    }
  }
  
  const base = {
    method,
    pdata_dir: pdataDir,
    meta: cleanedMeta,  // [NOTE] Cleaned meta
    ...(cleanedSize ? { size: cleanedSize } : {}),  // [NOTE] Cleaned size or absent
    level_start: Number.isFinite(+opts.levelStart) ? +opts.levelStart : 0.005,
    ...(opts.nLevels != null ? { n_levels: +opts.nLevels } : {}),
    ...(opts.mode ? { mode: String(opts.mode) } : {}),
    ...(opts.cmap ? { cmap: String(opts.cmap) } : {}),
    ...(opts.transparent != null ? { transparent: +!!opts.transparent } : {}),
  };
  
  if (opts.taskId) base.taskId = String(opts.taskId);
  
  // [DEBUG] TEMPORARY (remove after test)
  console.log('[DEBUG] [RELEVEL] Payload sent:', JSON.stringify(base, null, 2));
  
  const req = postJSONCancelable(
    "/spectral-bg/relevel",
    withPersist(base, {
      persistKey: opts.taskId || DEFAULT_PERSIST_KEY,
      overwrite: true,
    })
  );
  
  const wrapped = req.then((data) => _wrapCancelable(data, method));
  wrapped.cancel = req.cancel;
  return wrapped;
}

/**
 * Cancelable interactive (1D/2D) — identical to spectralBgInteractive but cancelable
 */
export function spectralBgInteractiveCancelable(method, pdataDir, opts = {}) {
  const base = { method, pdata_dir: pdataDir };
  if (opts.taskId) base.taskId = String(opts.taskId);

  if (opts.size) {
    const { width_px, height_px, dpi } = opts.size;
    base.size = {
      ...(Number.isFinite(+width_px) ? { width_px: +width_px } : {}),
      ...(Number.isFinite(+height_px) ? { height_px: +height_px } : {}),
      ...(Number.isFinite(+dpi) ? { dpi: +dpi } : {}),
    };
  }

  if (method === "1H" || method === "13C") {
    if (opts.window) {
      const { ppm_min, ppm_max } = opts.window;
      const w = {};
      if (Number.isFinite(+ppm_min)) w.ppm_min = +ppm_min;
      if (Number.isFinite(+ppm_max)) w.ppm_max = +ppm_max;
      if (Object.keys(w).length) base.window = w;
    }
    if (opts.line) {
      const { color, width } = opts.line;
      base.line = {
        ...(color ? { color: String(color) } : {}),
        ...(Number.isFinite(+width) ? { width: +width } : {}),
      };
    }
  } else {
    if (opts.window2d) {
      const { f2_ppm_min, f2_ppm_max, f1_ppm_min, f1_ppm_max } = opts.window2d;
      base.window2d = {
        ...(Number.isFinite(+f2_ppm_min) ? { f2_ppm_min: +f2_ppm_min } : {}),
        ...(Number.isFinite(+f2_ppm_max) ? { f2_ppm_max: +f2_ppm_max } : {}),
        ...(Number.isFinite(+f1_ppm_min) ? { f1_ppm_min: +f1_ppm_min } : {}),
        ...(Number.isFinite(+f1_ppm_max) ? { f1_ppm_max: +f1_ppm_max } : {}),
      };
    }
    const { level_start = 0.005, n_levels, mode, cmap, transparent } = opts.contour || {};
    const t = coerceTransparent(transparent);
    base.contour = {
      ...(Number.isFinite(+level_start) ? { level_start: +level_start } : { level_start: 0.005 }),
      ...(Number.isFinite(+n_levels) ? { n_levels: +n_levels } : {}),
      ...(mode ? { mode: String(mode) } : {}),
      ...(cmap ? { cmap: String(cmap) } : {}),
      ...(t !== undefined ? { transparent: t } : {}),
    };
  }

  const req = postJSONCancelable(
    "/spectral-bg/interactive",
    withPersist(base, {
      persistKey: opts.taskId || DEFAULT_PERSIST_KEY,
      overwrite: true,
    })
  );

  const wrapped = req.then((data) => _wrapCancelable(data, method));
  wrapped.cancel = req.cancel;
  return wrapped;
}

export const interactive1H = (pdataDir, { ppmMin, ppmMax, taskId, size, line } = {}) =>
  spectralBgInteractive("1H", pdataDir, {
    taskId,
    size,
    line,
    window: { 
      ppm_min: Number.isFinite(+ppmMin) ? +ppmMin : undefined, 
      ppm_max: Number.isFinite(+ppmMax) ? +ppmMax : undefined 
    },
  });

export const interactive13C = (pdataDir, { ppmMin, ppmMax, taskId, size, line } = {}) =>
  spectralBgInteractive("13C", pdataDir, {
    taskId,
    size,
    line,
    window: { 
      ppm_min: Number.isFinite(+ppmMin) ? +ppmMin : undefined, 
      ppm_max: Number.isFinite(+ppmMax) ? +ppmMax : undefined 
    },
  });

export const interactiveCOSY = (pdataDir, { f2Min, f2Max, f1Min, f1Max, taskId, size, contour } = {}) =>
  spectralBgInteractive("COSY", pdataDir, {
    taskId,
    size,
    contour: { level_start: 0.005, ...(contour || {}) },
    window2d: { f2_ppm_min: f2Min, f2_ppm_max: f2Max, f1_ppm_min: f1Min, f1_ppm_max: f1Max },
  });

export const interactiveHSQC = (pdataDir, { f2Min, f2Max, f1Min, f1Max, taskId, size, contour } = {}) =>
  spectralBgInteractive("HSQC", pdataDir, {
    taskId,
    size,
    contour: { level_start: 0.005, ...(contour || {}) },
    window2d: { f2_ppm_min: f2Min, f2_ppm_max: f2Max, f1_ppm_min: f1Min, f1_ppm_max: f1Max },
  });

export const zoom1DFromPixels = (
  method,
  pdataDir,
  meta,
  rectPx,
  { taskId, size, line, displayWidth, displayHeight, displaySize } = {}
) =>
  spectralBgZoom(method, pdataDir, meta, rectPx, {
    taskId,
    size,
    line,
    displayWidth,
    displayHeight,
    displaySize,
  });

export const zoom2DFromPixels = (
  method,
  pdataDir,
  meta,
  rectPx,
  { taskId, size, contour, displayWidth, displayHeight, displaySize } = {}
) =>
  spectralBgZoom(method, pdataDir, meta, rectPx, {
    taskId,
    size,
    contour,
    displayWidth,
    displayHeight,
    displaySize,
  });

export async function zoom1DStrict(method, pdataDir, meta, rectPx, {
  taskId,
  line,
  displayWidth,
  displayHeight,
  size,
} = {}) {
  if (!["1H","13C"].includes(method)) {
    throw new Error("zoom1DStrict: method must be '1H' or '13C'");
  }

  const win = compute1DWindowFromRect(meta, rectPx, displayWidth, displayHeight);

  if (win.width_px_native < 2 && Number.isFinite(meta.plot_bbox_px?.w)) {
    const mid = (win.ppm_min + win.ppm_max) / 2;
    const pad = (meta.ppm_left - meta.ppm_right) * 0.002;
    win.ppm_min = mid - pad;
    win.ppm_max = mid + pad;
  }

  return spectralBgInteractive(method, pdataDir, {
    taskId,
    size,
    line,
    window: { ppm_min: win.ppm_min, ppm_max: win.ppm_max },
  });
}

export const computeWheelLevelStart = (
  current = 0.005,
  deltaY = 0,
  step = 1.2,
  min = 1e-6,
  max = 1
) => {
  const cur = Number.isFinite(+current) && +current > 0 ? +current : 0.005;
  const factor = deltaY < 0 ? 1 / step : step;
  const v = cur * factor;
  return Math.max(min, Math.min(max, v));
};

export const relevelByWheel = (
  method,
  pdataDir,
  meta,
  deltaY,
  { taskId, mode, cmap, transparent, nLevels = 20, step = 1.2 } = {}
) => {
  const current = Number.isFinite(+meta?.level_start) ? +meta.level_start : 0.005;
  const levelStart = computeWheelLevelStart(current, deltaY, step);
  return spectralBgRelevel(method, pdataDir, meta, {
    taskId,
    mode,
    cmap,
    transparent,
    nLevels,
    levelStart,
  });
};

/**
 * Relevel by wheel (cancelable) — returns a promise with .cancel()
 */
export function relevelByWheelCancelable(
  method,
  pdataDir,
  meta,
  deltaY,
  { taskId, mode, cmap, transparent, nLevels = 20, step = 1.2 } = {}
) {
  const current = Number.isFinite(+meta?.level_start) && +meta.level_start > 0
    ? +meta.level_start
    : 0.005;
  const levelStart = computeWheelLevelStart(current, deltaY, step);

  const p = spectralBgRelevelCancelable(method, pdataDir, meta, {
    taskId,
    mode,
    cmap,
    transparent,
    nLevels,
    levelStart,
  });
  return p; // p.cancel() available
}

export const startNmrInteract = (
  method,
  pdataDir,
  { mode2d = "contour" } = {},
  onLog,
  onError,
  onComplete,
  onTaskId
) => {
  const url = `${API_BASE_URL}/spectral-bg/interact`;
  const body = withPersist({ method, pdata_dir: pdataDir, ...(method === "COSY" || method === "HSQC" ? { mode2d } : {}) });

  return createPostSSEConnection(url, body, {}, {
    onProgress: null,
    onLog,
    onError,
    onComplete,
    onTaskId,
  });
};

export const startInteract1H = (pdataDir, handlers = {}) =>
  startNmrInteract(
    "1H",
    pdataDir,
    {},
    handlers.onLog,
    handlers.onError,
    handlers.onComplete,
    handlers.onTaskId
  );

export const startInteract13C = (pdataDir, handlers = {}) =>
  startNmrInteract(
    "13C",
    pdataDir,
    {},
    handlers.onLog,
    handlers.onError,
    handlers.onComplete,
    handlers.onTaskId
  );

export const startInteractCOSY = (pdataDir, { mode2d = "contour", ...handlers } = {}) =>
  startNmrInteract(
    "COSY",
    pdataDir,
    { mode2d },
    handlers.onLog,
    handlers.onError,
    handlers.onComplete,
    handlers.onTaskId
  );

export const startInteractHSQC = (pdataDir, { mode2d = "contour", ...handlers } = {}) =>
  startNmrInteract(
    "HSQC",
    pdataDir,
    { mode2d },
    handlers.onLog,
    handlers.onError,
    handlers.onComplete,
    handlers.onTaskId
  );

export const stopNmrInteract = (taskId) => stopTask(taskId);

export const debugSpectralMeta = (meta, method = 'unknown') => {
  console.group(`[DEBUG] Debug metadata ${method}`);
  
  if (!meta) {
    console.warn('[ERROR] No metadata provided');
    console.groupEnd();
    return false;
  }
  
  console.log('[INFO] Complete metadata:', meta);
  console.log('[INFO] Figure dimensions:', meta.fig_px);
  console.log('[INFO] Bbox plot:', meta.plot_bbox_px);
  
  if (method === '1H' || method === '13C') {
    console.log('[INFO] PPM 1D:', {
      ppm_left: meta.ppm_left,
      ppm_right: meta.ppm_right,
      ppm_min: meta.ppm_min,
      ppm_max: meta.ppm_max,
    });
    
    const hasValidBbox = meta.plot_bbox_px && 
      typeof meta.plot_bbox_px === 'object' &&
      [meta.plot_bbox_px.x, meta.plot_bbox_px.y, meta.plot_bbox_px.w, meta.plot_bbox_px.h]
        .every(Number.isFinite);

    const hasValidPpm = (Number.isFinite(+meta.ppm_left) && Number.isFinite(+meta.ppm_right)) ||
                       (Number.isFinite(+meta.ppm_min) && Number.isFinite(+meta.ppm_max));
    
    console.log('[OK] Validations:', { hasValidBbox, hasValidPpm });
    
    if (!hasValidBbox) console.warn('[WARN] plot_bbox_px invalid or missing');
    if (!hasValidPpm) console.warn('[WARN] No valid PPM coordinates');
    
    console.groupEnd();
    return hasValidBbox && hasValidPpm;
  }
  
  console.groupEnd();
  return true;
};

export default api;