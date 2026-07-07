// src/utils/molParser.js
// Parsing utilities for MOL/SDF files and molecule manipulation

import { API_BASE_URL } from "../api";

/**
 * Extracts the first MOL block from an SDF text
 */
export function extractFirstMolBlockPreserving(text) {
  if (!text) return null;
  const i = text.indexOf("M  END");
  if (i < 0) return null;
  const end = i + "M  END".length;
  let out = text.slice(0, end);
  if (!out.endsWith("\n")) out += "\n";
  return out;
}

/**
 * Extracts the first SDF record
 */
export function extractFirstSDFRecord(text) {
  if (!text) return null;
  const upto = text.indexOf("$$$$");
  const rec = upto >= 0 ? text.slice(0, upto) : text;
  const mend = rec.indexOf("M  END");
  if (mend >= 0) {
    const s = rec.slice(0, mend + "M  END".length);
    return s.endsWith("\n") ? s : s + "\n";
  }
  return rec.endsWith("\n") ? rec : rec + "\n";
}

/**
 * Extracts chemical shifts from an SDF record
 */
export function extractShiftArrayFromSDFRecord(sdfText) {
  if (!sdfText) return null;
  const first = sdfText.split("$$$$")[0] || sdfText;

  // Format: >  <(1)>  (1)\n22.067 54.492 0.000 ...
  const re = />\s*<\(1\)>\s*(?:\(\d+\))?\s*\r?\n([\s\S]*?)(?:\r?\n\r?\n|\$\$\$\$|$)/;
  const m = first.match(re);
  if (!m) {
    // Fallback: also try <1> without parentheses
    const re2 = />\s*<1>\s*(?:\(\d+\))?\s*\r?\n([\s\S]*?)(?:\r?\n\r?\n|$$$$|$)/;
    const m2 = first.match(re2);
    if (!m2) return null;
    return parseShiftBlock(m2[1]);
  }

  return parseShiftBlock(m[1]);
}

/**
 * Parse un bloc de valeurs de shift
 */
function parseShiftBlock(block) {
  if (!block) return null;
  const cleaned = block
    .split(/\r?\n/)
    .map((s) => s.trim())
    .filter(Boolean)
    .join(" ");

  if (!cleaned) return null;

  return cleaned
    .split(/\s+/)
    .map((x) => parseFloat(x.replace(",", ".")))
    .map((v) => (Number.isFinite(v) ? v : null));
}

/**
 * Parse les symboles atomiques depuis un bloc MOL
 */
/**
 * Supprime les H explicites d'un molblock V2000 pour un rendu 2D propre.
 * Retourne aussi le mapping oldIndex → newIndex pour remapper les highlights.
 */
export function stripExplicitHFromMolblock(molblock) {
  if (!molblock) return { molblock, indexMap: {} };
  const lines = molblock.split(/\r?\n/);
  if (lines.length < 4) return { molblock, indexMap: {} };

  const countsLine = lines[3];
  const nAtoms = parseInt(countsLine.slice(0, 3).trim(), 10);
  const nBonds = parseInt(countsLine.slice(3, 6).trim(), 10);
  if (isNaN(nAtoms) || isNaN(nBonds)) return { molblock, indexMap: {} };

  const hIndices = new Set();
  for (let i = 0; i < nAtoms; i++) {
    const sym = (lines[4 + i] || '').slice(31, 34).trim();
    if (sym === 'H') hIndices.add(i);
  }
  if (hIndices.size === 0) {
    const indexMap = {};
    for (let i = 0; i < nAtoms; i++) indexMap[i] = i;
    return { molblock, indexMap };
  }

  const indexMap = {};
  let newIdx = 0;
  for (let i = 0; i < nAtoms; i++) {
    if (!hIndices.has(i)) indexMap[i] = newIdx++;
  }

  const newAtomLines = [];
  for (let i = 0; i < nAtoms; i++) {
    if (!hIndices.has(i)) newAtomLines.push(lines[4 + i]);
  }

  const newBondLines = [];
  const bondIndexMap = {}; // oldBondIdx → newBondIdx
  let newBondIdx = 0;
  for (let i = 0; i < nBonds; i++) {
    const line = lines[4 + nAtoms + i] || '';
    const a1 = parseInt(line.slice(0, 3).trim(), 10) - 1;
    const a2 = parseInt(line.slice(3, 6).trim(), 10) - 1;
    if (hIndices.has(a1) || hIndices.has(a2)) continue;
    bondIndexMap[i] = newBondIdx++;
    const na1 = String(indexMap[a1] + 1).padStart(3);
    const na2 = String(indexMap[a2] + 1).padStart(3);
    newBondLines.push(na1 + na2 + line.slice(6));
  }

  const newCounts = String(newAtomLines.length).padStart(3) +
    String(newBondLines.length).padStart(3) +
    countsLine.slice(6);

  const footer = lines.slice(4 + nAtoms + nBonds);
  const newMolblock = [...lines.slice(0, 3), newCounts, ...newAtomLines, ...newBondLines, ...footer].join('\n');
  return { molblock: newMolblock, indexMap, bondIndexMap };
}

export function parseAtomSymbolsFromMolblock(molblock) {
  if (!molblock) return [];
  const lines = molblock.split(/\r?\n/);
  const counts = lines[3] || "";
  const nAtoms = parseInt(counts.slice(0, 3).trim(), 10);
  const atoms = [];
  for (let i = 0; i < nAtoms; i++) {
    const l = lines[4 + i] || "";
    const sym = l.slice(31, 34).trim();
    atoms.push(sym || "?");
  }
  return atoms;
}

/**
 * Trouve les liaisons dans un chemin entre paires d'atomes
 */
export function findBondsInPath(molblock, atomPairs) {
  if (!molblock || !atomPairs.length) return [];
  const lines = molblock.split(/\r?\n/);
  const counts = lines[3] || "";
  const nAtoms = parseInt(counts.slice(0, 3).trim(), 10);
  const nBonds = parseInt(counts.slice(3, 6).trim(), 10);

  const adjacency = new Map();
  for (let i = 0; i < nAtoms; i++) adjacency.set(i, []);

  for (let i = 0; i < nBonds; i++) {
    const l = lines[4 + nAtoms + i] || "";
    const a1 = parseInt(l.slice(0, 3).trim(), 10) - 1;
    const a2 = parseInt(l.slice(3, 6).trim(), 10) - 1;
    if (a1 >= 0 && a2 >= 0 && a1 < nAtoms && a2 < nAtoms) {
      adjacency.get(a1).push({ neighbor: a2, bondIndex: i });
      adjacency.get(a2).push({ neighbor: a1, bondIndex: i });
    }
  }

  const findShortestPath = (start, end) => {
    // Pour les diagonaux COSY (même proton) : retourner les bonds connectés (liaison C-H)
    if (start === end) {
      const bonds = adjacency.get(start) || [];
      return bonds.map(b => b.bondIndex);
    }
    const q = [{ atom: start, path: [] }];
    const seen = new Set([start]);
    while (q.length > 0) {
      const { atom, path } = q.shift();
      for (const { neighbor, bondIndex } of adjacency.get(atom) || []) {
        if (neighbor === end) return [...path, bondIndex];
        if (!seen.has(neighbor) && path.length < 4) {
          seen.add(neighbor);
          q.push({ atom: neighbor, path: [...path, bondIndex] });
        }
      }
    }
    return [];
  };

  const hi = new Set();
  for (const [a, b] of atomPairs) {
    const p = findShortestPath(a, b);
    p.forEach((idx) => hi.add(idx));
  }
  return Array.from(hi);
}

/**
 * Attempts to fetch a text file from a URL
 */
export async function tryFetchText(url) {
  try {
    const r = await fetch(url, { method: "GET" });
    if (!r.ok) return null;
    const txt = await r.text();
    if (/<!doctype|<html/i.test(txt)) return null;
    return txt;
  } catch {
    return null;
  }
}

/**
 * Extracts the directory from a URL
 */
export function dirnameFromUrl(u) {
  try {
    const baseUrl =
      typeof API_BASE_URL === "string" && API_BASE_URL.startsWith("http")
        ? API_BASE_URL
        : window.location.href;
    const url = new URL(u, baseUrl);
    const p = url.pathname;
    return `${url.origin}${p.slice(0, p.lastIndexOf("/") + 1)}`;
  } catch {
    return null;
  }
}

/**
 * Resolves a MOL file for a given source
 */
export async function resolveMolForSource(csvUrl, sourceName, taskId = null) {
  if (!csvUrl || !sourceName) {
    return null;
  }
  const baseName = String(sourceName)
    .trim()
    .replace(/\.(sdf|mol)$/i, "");

  if (taskId) {
    try {
      const baseUrl =
        typeof API_BASE_URL === "string" && API_BASE_URL.startsWith("http")
          ? API_BASE_URL
          : window.location.href;
      const origin = new URL(csvUrl, baseUrl).origin;
      const base = `${origin}/outputs/${taskId}/`;
      const cand = [
        `${base}${baseName}.sdf`,
        `${base}${baseName}.mol`,
        `${base}${encodeURIComponent(baseName)}.sdf`,
        `${base}${encodeURIComponent(baseName)}.mol`,
      ];
      for (const u of cand) {
        const txt = await tryFetchText(u);
        if (txt) {
          return {
            url: u,
            rawText: txt,
            firstRecord: extractFirstSDFRecord(txt),
            molblock: extractFirstMolBlockPreserving(txt) || txt,
          };
        }
      }
    } catch (e) {
      // Silently ignore
    }
  }

  const csvDir = dirnameFromUrl(csvUrl);
  if (csvDir) {
    const cand = [
      `${csvDir}${encodeURIComponent(baseName)}.sdf`,
      `${csvDir}${encodeURIComponent(baseName)}.mol`,
      `${csvDir}${baseName}.sdf`,
      `${csvDir}${baseName}.mol`,
      `${csvDir}${baseName}.SDF`,
      `${csvDir}${baseName}.MOL`,
    ];
    for (const u of cand) {
      const txt = await tryFetchText(u);
      if (txt)
        return {
          url: u,
          rawText: txt,
          firstRecord: extractFirstSDFRecord(txt),
          molblock: extractFirstMolBlockPreserving(txt) || txt,
        };
    }
  }

  try {
    const baseUrl =
      typeof API_BASE_URL === "string" && API_BASE_URL.startsWith("http")
        ? API_BASE_URL
        : window.location.href;
    const origin = new URL(csvUrl, baseUrl).origin;
    const api = `${origin}/api/prediction-tasks`;
    const r = await fetch(api);
    if (r.ok) {
      const tasks = await r.json();
      for (const t of tasks || []) {
        if (taskId && t.taskId === taskId) continue;
        const base = `${origin}/outputs/${t.taskId}/`;
        const cand = [
          `${base}${baseName}.sdf`,
          `${base}${baseName}.mol`,
          `${base}${encodeURIComponent(baseName)}.sdf`,
          `${base}${encodeURIComponent(baseName)}.mol`,
        ];
        for (const u of cand) {
          const txt = await tryFetchText(u);
          if (txt)
            return {
              url: u,
              rawText: txt,
              firstRecord: extractFirstSDFRecord(txt),
              molblock: extractFirstMolBlockPreserving(txt) || txt,
            };
        }
      }
    }
  } catch {
    // Silently ignore
  }

  return null;
}
