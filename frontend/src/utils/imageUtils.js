// src/utils/imageUtils.js
// Utilities for image management and cache-busting

import { API_BASE_URL } from "../api";

/**
 * Builds an image URL with cache-busting
 * Only adds ?cb=... if URL is not signed.
 * Also takes into account meta.nonce to bust cache when backend recalculates.
 */
export function buildImgSrc(url, meta) {
  try {
    // En Electron, utiliser API_BASE_URL comme base au lieu de window.location.href (file://)
    const baseUrl =
      typeof API_BASE_URL === "string" && API_BASE_URL.startsWith("http")
        ? API_BASE_URL
        : window.location.href;
    const u = new URL(url, baseUrl);

    // Tag prioritaire pour le cache-busting
    const tag =
      meta?.etag ?? meta?.version ?? meta?.last_modified ?? meta?.nonce;
    if (!tag) return url;

    const qs = u.search || "";
    const looksSigned =
      /(token|signature|sig|expires|X-Amz-Algorithm|X-Amz-Credential|X-Amz-Signature)/i.test(
        qs
      );

    if (looksSigned) {
      // Do not modify query (otherwise signature invalid)
      // -> force reload via fragment (#cb=...)
      const currentHash = (u.hash || "").replace(/^#/, "");
      const parts = currentHash ? currentHash.split("&").filter(Boolean) : [];
      // remove existing cb if present
      const others = parts.filter((p) => !/^cb=/.test(p));
      others.push("cb=" + encodeURIComponent(tag));
      u.hash = "#" + others.join("&");
      return u.toString();
    }

    // Unsigned URL → we can use the query
    const hasQuery = u.search && u.search.length > 1;
    const params = new URLSearchParams(u.search);
    params.set("cb", String(tag));
    u.search = (hasQuery ? "?" : "?") + params.toString();
    return u.toString();
  } catch {
    return url;
  }
}

/**
 * Fetch with automatic retry
 */
export async function fetchWithRetry(url, maxRetries = 3, delayMs = 500) {
  let lastError = null;

  for (let attempt = 1; attempt <= maxRetries; attempt++) {
    try {
      const r = await fetch(url);
      if (!r.ok) throw new Error(`HTTP ${r.status}`);
      return r;
    } catch (e) {
      lastError = e;
      console.warn(
        `Attempt ${attempt}/${maxRetries} failed:`,
        e?.message || e
      );
      if (attempt < maxRetries) {
        const d = delayMs * attempt;
        await new Promise((res) => setTimeout(res, d));
      }
    }
  }
  throw lastError;
}
