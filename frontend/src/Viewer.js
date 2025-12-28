// src/Viewer.js - CORRECTED VERSION with timeouts, global color management, forced remount, reliable 2D reset
// + "Show all molecules" button with paginated modal sorted by RANSY average
// Combined view NMR + 2D/3D Molecules + sidebar
// - Rect zoom, 2D relevel on scroll, state persistence
// - "Word style" color picker with Hex/RGB + preview
// - Colors by method for matched atoms (H1, C13, COSY, HSQC)
// - FIX: safety timeouts to avoid UI blocking
// - FIX: reliable 1D/2D reset via nonce + backend recompute + snapshots by type and source

import React, {
  useEffect,
  useMemo,
  useState,
  useCallback,
  useRef,
} from "react";

import {
  interactive1H,
  interactive13C,
  interactiveCOSY,
  interactiveHSQC,
  relevelByWheelCancelable,
  zoom1DFromPixels,
  API_BASE_URL,
} from "./api";

// ========================================
// MODULARIZED IMPORTS
// ========================================

// Color constants and styles
import {
  DEFAULT_COLORS,
  METHOD_LABELS_SHORT,
  TYPE_LABELS,
  PLOT_HEIGHT,
  MOL_HEIGHT,
  SIDEBAR_W,
} from "./constants/colors";

// Utilitaires
import { debugLog, debugWarn } from "./utils/debug";
import { fetchWithRetry } from "./utils/imageUtils";
import { 
  loadRDKitFromCDN, 
  load3DmolFromCDN, 
  isWebGLAvailable 
} from "./utils/cdnLoaders";
// colorUtils used only in sub-components

// Composants viewer
import { Pill, ColorPicker, StickChart, Scatter2D, SpectrumDetailModal } from "./components/viewer";

// Hook de mesure
// useMeasureWidth used only in sub-components

// Enable debug by typing in console: window.DEBUG_NMR = true
if (typeof window !== 'undefined' && !window.hasOwnProperty('DEBUG_NMR')) {
  window.DEBUG_NMR = false; // false by default in production
}

/* ===========================
   CDN Loaders (RDKit, 3Dmol) - imported from utils/cdnLoaders.js
   Helpers isProblematicBrowser, isWebGLAvailable also imported
   =========================== */

/* ===========================
   CSV parsing + petits helpers
   =========================== */

function parseCSVSmart(text) {
  if (!text || typeof text !== "string") {
    return { headers: [], data: [] };
  }
  const first = text.split(/\r?\n/).slice(0, 5).join("\n");
  const cand = [";", ",", "\t"];
  let best = ";";
  let bestCount = -1;

  for (const d of cand) {
    const counts = first
      .split(/\r?\n/)
      .map(l => (l.match(new RegExp(`\\${d}`, "g")) || []).length);
    const cnt = Math.max(...counts, 0);
    if (cnt > bestCount) {
      best = d;
      bestCount = cnt;
    }
  }

  const rows = [];
  let cur = [];
  let val = "";
  let inQuotes = false;

  const pushVal = () => {
    if (val && val.charCodeAt(0) === 0xfeff) val = val.slice(1);
    cur.push(val);
    val = "";
  };
  const pushRow = () => { rows.push(cur); cur = []; };

  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (inQuotes) {
      if (ch === '"') {
        if (text[i + 1] === '"') { val += '"'; i++; }
        else inQuotes = false;
      } else {
        val += ch;
      }
      continue;
    }
    if (ch === '"') inQuotes = true;
    else if (ch === best) pushVal();
    else if (ch === "\n") { pushVal(); pushRow(); }
    else if (ch !== "\r") val += ch;
  }
  if (val.length > 0 || cur.length > 0) {
    pushVal();
    pushRow();
  }
  if (rows.length === 0) return { headers: [], data: [] };

  const headers = rows[0].map(h => (h || "").trim());
  const data = rows
    .slice(1)
    .filter(r => r.some(x => (x || "").trim() !== ""))
    .map(r => {
      const o = {};
      headers.forEach((h, i) => { o[h] = (r[i] ?? "").trim(); });
      return o;
    });

  return { headers, data };
}


/* ===========================
   Graph helpers (bonds)
   =========================== */

function findBondsInPath(molblock, atomPairs) {
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
    if (start === end) return [];
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
    p.forEach(idx => hi.add(idx));
  }
  return Array.from(hi);
}



/* ===========================
   SDF / MOL helpers
   =========================== */

function extractFirstMolBlockPreserving(text) {
  if (!text) return null;
  const i = text.indexOf("M  END");
  if (i < 0) return null;
  const end = i + "M  END".length;
  let out = text.slice(0, end);
  if (!out.endsWith("\n")) out += "\n";
  return out;
}

function extractFirstSDFRecord(text) {
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

function extractShiftArrayFromSDFRecord(sdfText) {
  if (!sdfText) return null;
  const first = sdfText.split("$$$$")[0] || sdfText;
  
  // [FIX] CORRECTION: search for <(1)> with backend exact format: >  <(1)>  (1)
  // Format is: >  <(1)>  (1)\n22.067 54.492 0.000 ...
  const re = />\s*<\(1\)>\s*(?:\(\d+\))?\s*\r?\n([\s\S]*?)(?:\r?\n\r?\n|$$$$|$)/;
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

function parseShiftBlock(block) {
  if (!block) return null;
  const cleaned = block
    .split(/\r?\n/)
    .map(s => s.trim())
    .filter(Boolean)
    .join(" ");

  if (!cleaned) return null;

  // [NOTE] Values are already space-separated
  return cleaned
    .split(/\s+/)
    .map(x => parseFloat(x.replace(",", ".")))
    .map(v => (Number.isFinite(v) ? v : null));
}

function parseAtomSymbolsFromMolblock(molblock) {
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

async function tryFetchText(url) {
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

function dirnameFromUrl(u) {
  try {
    // En Electron, utiliser API_BASE_URL comme base
    const baseUrl = typeof API_BASE_URL === 'string' && API_BASE_URL.startsWith('http')
      ? API_BASE_URL
      : window.location.href;
    const url = new URL(u, baseUrl);
    const p = url.pathname;
    return `${url.origin}${p.slice(0, p.lastIndexOf("/") + 1)}`;
  } catch {
    return null;
  }
}

async function resolveMolForSource(csvUrl, sourceName, taskId = null) {
  if (!csvUrl || !sourceName) {
    debugLog('[resolveMolForSource] Empty csvUrl or sourceName:', { csvUrl, sourceName });
    return null;
  }
  const baseName = String(sourceName)
    .trim()
    .replace(/\.(sdf|mol)$/i, "");

  debugLog('[resolveMolForSource] Searching for:', { csvUrl, sourceName, taskId, baseName });

  if (taskId) {
    try {
      // En Electron, utiliser API_BASE_URL pour construire l'URL de base
      const baseUrl = typeof API_BASE_URL === 'string' && API_BASE_URL.startsWith('http')
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
      debugLog('[resolveMolForSource] Trying taskId candidates:', cand);
      for (const u of cand) {
        const txt = await tryFetchText(u);
        if (txt) {
          debugLog('[resolveMolForSource] [OK] Found at:', u);
          return {
            url: u,
            rawText: txt,
            firstRecord: extractFirstSDFRecord(txt),
            molblock: extractFirstMolBlockPreserving(txt) || txt,
          };
        }
      }
      debugLog('[resolveMolForSource] [ERROR] Not found via taskId');
    } catch (e) {
      debugLog('[resolveMolForSource] Error with taskId:', e?.message);
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
    // En Electron, utiliser API_BASE_URL pour construire l'URL de base
    const baseUrl = typeof API_BASE_URL === 'string' && API_BASE_URL.startsWith('http')
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
  } catch {}

  return null;
}

/* useMeasureWidth importe depuis hooks/useMeasureWidth.js */

/* StickChart importe depuis components/viewer/StickChart.jsx */

/* Scatter2D importe depuis components/viewer/Scatter2D.jsx */

function MoleculePane({
  csvUrl,
  selectedRow,
  show3D,
  taskId,
  multi = null,
  customColors,
  molZoom = 1,
  molHeight = MOL_HEIGHT,
  flagsOverride = null,            // coches globales (optionnelles)
  showLocalToggles = true,         // masquer/afficher les coches locales
  compact = false,                 // mode miniatures: juste le SVG, sans cadres/titres
}) {
  const [svg, setSvg] = useState(null);
  const [err, setErr] = useState(null);
  const [loading, setLoading] = useState(false);

  // Local checkboxes (used only if flagsOverride is not provided)
  const [localFlags, setLocalFlags] = useState({
    H1: true, C13: true, COSY: true, HSQC: true,
  });
  const flags = flagsOverride || localFlags;       // single source of truth
  const toggle = (k) => setLocalFlags(f => ({ ...f, [k]: !f[k] }));

  const div3DRef = useRef(null);
  const viewerRef = useRef(null);

  // eslint-disable-next-line react-hooks/exhaustive-deps
  const COLORS_RDK = {
    H1:  customColors.H1.rgb,
    C13: customColors.C13.rgb,
    COSY:customColors.COSY.rgb,
    HSQC:customColors.HSQC.rgb,
  };
  // eslint-disable-next-line react-hooks/exhaustive-deps
  const COLORS_3D = {
    H1:  customColors.H1.hex,
    C13: customColors.C13.hex,
    COSY:customColors.COSY.hex,
    HSQC:customColors.HSQC.hex,
  };

  // --- HSQC: colorer la liaison C–H en 2D -----------------------------------
  const computeMultiHighlights = useCallback((resolved, rowsByType) => {
    const atoms  = parseAtomSymbolsFromMolblock(resolved.molblock) || [];
    const shifts = extractShiftArrayFromSDFRecord(resolved.rawText) || [];
    
    // [FIX] CORRECTION : Logs conditionnels uniquement en mode debug
    debugLog("[DEBUG] Shifts extraits:", shifts);
    debugLog("[DEBUG] Symboles d'atomes:", atoms);
    debugLog("[DEBUG] Nombre d'atomes:", atoms.length);
    
    // Protection contre les tableaux vides
    if (!atoms.length || !shifts.length) {
      debugWarn("[WARN] No atoms or shifts found, coloring disabled");
      return { 
        highlightByType: { H1: new Set(), C13: new Set(), COSY: new Set(), HSQC: new Set() }, 
        bondsCOSY: [], 
        bondsHSQC: [] 
      };
    }
    
    const maxIdx = Math.min(shifts.length - 1, atoms.length - 1);

    const hIdx   = atoms.map((s,i)=>s==="H"?i:-1).filter(i=>i>=0);
    const cIdx   = atoms.map((s,i)=>s==="C"?i:-1).filter(i=>i>=0);

    const sets   = { H1:new Set(), C13:new Set(), COSY:new Set(), HSQC:new Set() };
    let bondsCOSY = [];
    let bondsHSQC = [];

    // 1H
    if (flags.H1 && rowsByType["1H"]) {
      const m = (rowsByType["1H"].ppm_molecule_matched||"")
        .split("|").map(v=>parseFloat(String(v).replace(",","."))).filter(Number.isFinite);
      const r = new Set(m.map(v=>Math.round(v*1000)/1000));
      for (const ppm of r) for (const i of hIdx)
        if (i<=maxIdx && (Math.round(shifts[i]*1000)/1000)===ppm) sets.H1.add(i);
    }
    
    // 13C
    if (flags.C13 && rowsByType["13C"]) {
      const m = (rowsByType["13C"].ppm_molecule_matched||"")
        .split("|").map(v=>parseFloat(String(v).replace(",","."))).filter(Number.isFinite);
      const r = new Set(m.map(v=>Math.round(v*1000)/1000));
      for (const ppm of r) for (const i of cIdx)
        if (i<=maxIdx && (Math.round(shifts[i]*1000)/1000)===ppm) sets.C13.add(i);
    }
    
    // COSY (chemins min)
    if (flags.COSY && rowsByType["COSY"]) {
      const pairs = (rowsByType["COSY"].ppm_molecule_matched||"")
        .split("||").map(p=>p.split("|")).filter(a=>a.length===2)
        .map(([f1,f2])=>[+String(f1).replace(",","."),+String(f2).replace(",","." )])
        .filter(([a,b])=>Number.isFinite(a)&&Number.isFinite(b));
      const matched = [];
      for (const [d1,d2] of pairs) {
        const r1=Math.round(d1*1000)/1000, r2=Math.round(d2*1000)/1000;
        const a1=[], a2=[];
        for (const i of hIdx) if (i<=maxIdx&&(Math.round(shifts[i]*1000)/1000)===r1) a1.push(i);
        for (const i of hIdx) if (i<=maxIdx&&(Math.round(shifts[i]*1000)/1000)===r2) a2.push(i);
        for (const x of a1) for (const y of a2) if (x!==y) matched.push([x,y]);
      }
      if (matched.length) {
        bondsCOSY = findBondsInPath(resolved.molblock, matched);
        const lines = resolved.molblock.split(/\r?\n/);
        const nAtoms = parseInt((lines[3]||"").slice(0,3).trim(),10);
        for (const b of bondsCOSY) {
          const l = lines[4+nAtoms+b]||"";
          const a1 = parseInt(l.slice(0,3).trim(),10)-1;
          const a2 = parseInt(l.slice(3,6).trim(),10)-1;
          if (a1>=0) sets.COSY.add(a1);
          if (a2>=0) sets.COSY.add(a2);
        }
      }
    }
    
    // HSQC (liaison directe C–H)
    if (flags.HSQC && rowsByType["HSQC"]) {
      const pairs = (rowsByType["HSQC"].ppm_molecule_matched_unfolded ||
                     rowsByType["HSQC"].ppm_molecule_matched || "")
        .split("||").map(p=>p.split("|")).filter(a=>a.length===2)
        .map(([c,h])=>[+String(c).replace(",","."), +String(h).replace(",","." )])
        .filter(([c,h])=>Number.isFinite(c)&&Number.isFinite(h));

      const lines = resolved.molblock.split(/\r?\n/);
      const counts = lines[3] || "";
      const nAtoms = parseInt(counts.slice(0,3).trim(),10);
      const nBonds = parseInt(counts.slice(3,6).trim(),10);
      const bondMap = new Map();
      for (let bi=0; bi<nBonds; bi++) {
        const l = lines[4+nAtoms+bi] || "";
        const a1 = parseInt(l.slice(0,3).trim(),10)-1;
        const a2 = parseInt(l.slice(3,6).trim(),10)-1;
        if (a1>=0 && a2>=0) {
          const key = `${Math.min(a1,a2)}-${Math.max(a1,a2)}`;
          bondMap.set(key, bi);
        }
      }

      for (const [cppm,hppm] of pairs) {
        const rc = Math.round(cppm*1000)/1000;
        const rh = Math.round(hppm*1000)/1000;

        const cCandidates = [];
        const hCandidates = [];

        for (const ci of cIdx)
          if (ci<=maxIdx && (Math.round(shifts[ci]*1000)/1000)===rc) cCandidates.push(ci);

        for (const hi of hIdx)
          if (hi<=maxIdx && (Math.round(shifts[hi]*1000)/1000)===rh) {
            hCandidates.push(hi);
            sets.HSQC.add(hi);
          }

        for (const ci of cCandidates) for (const hi of hCandidates) {
          const key = `${Math.min(ci,hi)}-${Math.max(ci,hi)}`;
          const b = bondMap.get(key);
          if (b!=null) bondsHSQC.push(b);
        }
      }
    }
    
    return { highlightByType: sets, bondsCOSY, bondsHSQC };
  }, [flags]);
  // ---------------------------------------------------------------------------

  const render2D = useCallback(async (resolved) => {
    setLoading(true);
    try {
      const RDKit = await loadRDKitFromCDN();
      const mol = RDKit.get_mol(resolved.molblock);
      try { mol.generate_aligned_coords(); } catch { try { mol.set_new_coords(); } catch {} }

      let highlightAtoms=[], atomColorsMap={}, highlightBonds=[], bondColorsMap={};

      if (multi?.fromAll && multi?.rowsByType) {
        const { highlightByType, bondsCOSY, bondsHSQC } =
          computeMultiHighlights(resolved, multi.rowsByType);

        // Atom colors (priority HSQC > COSY > H1 for H atoms)
        if (flags.C13)
          for (const i of highlightByType.C13)
            atomColorsMap[String(i)] = COLORS_RDK.C13;

        const h1List = Array.from(highlightByType.H1);
        const cosySet = highlightByType.COSY;
        const hsqcSet = highlightByType.HSQC;

        for (const i of h1List) {
          if (flags.HSQC && hsqcSet.has(i))      atomColorsMap[String(i)] = COLORS_RDK.HSQC;
          else if (flags.COSY && cosySet.has(i)) atomColorsMap[String(i)] = COLORS_RDK.COSY;
          else if (flags.H1)                     atomColorsMap[String(i)] = COLORS_RDK.H1;
        }

        // Atoms to highlight
        const all = new Set([
          ...(flags.H1   ? Array.from(highlightByType.H1)   : []),
          ...(flags.C13  ? Array.from(highlightByType.C13)  : []),
          ...(flags.COSY ? Array.from(highlightByType.COSY) : []),
          ...(flags.HSQC ? Array.from(highlightByType.HSQC) : []),
        ]);
        highlightAtoms = Array.from(all);

        // Liaisons (COSY chemins + HSQC liaison directe C–H)
        if (flags.COSY) {
          highlightBonds.push(...bondsCOSY);
          for (const b of bondsCOSY) bondColorsMap[String(b)] = COLORS_RDK.COSY;
        }
        if (flags.HSQC) {
          highlightBonds.push(...bondsHSQC);
          for (const b of bondsHSQC) bondColorsMap[String(b)] = COLORS_RDK.HSQC;
        }
        // Deduplication
        highlightBonds = Array.from(new Set(highlightBonds));
      }

      const base = show3D ? 350 : 420;
      const sizePx = Math.round(base * Math.max(0.6, molZoom));
      const pad = Math.max(0.04, 0.24 / Math.max(1, molZoom));

      const opts = {
        width: sizePx, height: sizePx, useMolBlockCoords: false,
        padding: pad, bondLineWidth: 1.4, atomLabelFontSize: 10, annotationFontScale: .9,
      };

      let svgStr;
      try {
        svgStr = mol.get_svg_with_highlights(JSON.stringify({
          ...opts,
          atoms: highlightAtoms,
          highlight_atom_colors: atomColorsMap,
          bonds: highlightBonds,
          highlight_bond_colors: bondColorsMap,
          highlightRadius: 0.35,
          explicitMethyl: true,
          addStereoAnnotation: true,
          atomHighlightsAreCircles: true,
          highlightBondWidthMultiplier: 1.5,
        }));
      } catch {
        svgStr = mol.get_svg(JSON.stringify(opts));
      }

      svgStr = svgStr
        .replace(/width="[^"]+"/, 'width="100%"')
        .replace(/height="[^"]+"/, 'height="100%"')
        .replace('<svg ', '<svg preserveAspectRatio="xMidYMid meet" ');

      setSvg(svgStr);
      setErr(null);
      mol.delete?.();
    } catch (e) {
      setErr("2D Error: " + (e?.message || String(e)));
    } finally {
      setLoading(false);
    }
  }, [multi, computeMultiHighlights, COLORS_RDK, molZoom, show3D, flags]);

  const render3D = useCallback(async (resolved) => {
    setLoading(true);
    try {
      // [FIX] Check WebGL before trying to create the viewer
      if (!isWebGLAvailable()) {
        throw new Error("WebGL is not available on this browser.\n\n" +
          "Solutions:\n" +
          "1. Firefox: Type about:config in the address bar, then set webgl.disabled to false\n" +
          "2. Update your graphics drivers\n" +
          "3. Try Chrome/Edge if the problem persists");
      }

      const $3Dmol = await load3DmolFromCDN();
      const el = div3DRef.current;
      
      // [FIX] CHECK: Container must exist
      if (!el) {
        throw new Error("3D container not found");
      }
      
      // [FIX] CHECK: Container must have size
      const elWidth = el.offsetWidth;
      const elHeight = el.offsetHeight;
      
      if (elWidth === 0 || elHeight === 0) {
        debugWarn(`[WARN] 3D container not yet visible (${elWidth}x${elHeight}px), waiting...`);
        setTimeout(() => {
          if (div3DRef.current) {
            render3D(resolved);
          }
        }, 100);
        return;
      }
      
      if (!viewerRef.current) {
        viewerRef.current = (typeof $3Dmol.createViewer === "function")
          ? $3Dmol.createViewer(el, { 
              backgroundColor:"white", 
              antialias:false, 
              preserveDrawingBuffer:true, 
              disableOffscreenCanvas:true 
            })
          : new $3Dmol.GLViewer(el, { backgroundColor:"white" });
      }
      
      const viewer = viewerRef.current;
      
      // [FIX] CHECK: Viewer must be valid
      if (!viewer) {
        throw new Error("Failed to create 3D viewer");
      }
      
      viewer.removeAllModels();
      const sdf = resolved.firstRecord || resolved.rawText || resolved.molblock;
      const model = viewer.addModel(sdf, "sdf");
      
      // [FIX] CRITICAL CHECK: Model must exist
      if (!model) {
        throw new Error("Failed to load 3D model");
      }
      
      model.setStyle({}, { 
        stick:{radius:0.15,color:"#808080"}, 
        sphere:{scale:0.2,color:"#808080"} 
      });

      if (multi?.fromAll && multi?.rowsByType) {
        const { highlightByType } = computeMultiHighlights(resolved, multi.rowsByType);
        const CAtoms = flags.C13 ? Array.from(highlightByType.C13) : [];
        const cosySet = flags.COSY ? new Set(highlightByType.COSY) : new Set();
        const hsqcSet = flags.HSQC ? new Set(highlightByType.HSQC) : new Set();
        const h1List  = flags.H1   ? Array.from(highlightByType.H1)  : [];

        // [FIX] Colorier les atomes C13
        for (const i of CAtoms) {
          if (typeof i === 'number' && i >= 0) {
            try {
              model.setStyle({ index:i }, { 
                stick:{radius:0.15,color:"#808080"}, 
                sphere:{radius:0.4,color:COLORS_3D.C13} 
              });
            } catch (e) {
              debugWarn(`[WARN] C13 setStyle error for atom ${i}:`, e.message);
            }
          }
        }

        // [FIX] Colorier les atomes COSY
        if (flags.COSY) {
          for (const i of cosySet) {
            if (typeof i === 'number' && i >= 0) {
              try {
                model.setStyle({ index:i }, { 
                  stick:{radius:0.25,color:COLORS_3D.COSY}, 
                  sphere:{radius:0.4,color:CAtoms.includes(i)?COLORS_3D.C13:"#808080"} 
                });
              } catch (e) {
                debugWarn(`[WARN] COSY setStyle error for atom ${i}:`, e.message);
              }
            }
          }
        }

        // [FIX] Colorier les atomes HSQC
        if (flags.HSQC) {
          for (const i of hsqcSet) {
            if (typeof i === 'number' && i >= 0) {
              try {
                model.setStyle({ index:i }, { 
                  stick:{radius:0.2,color:COLORS_3D.HSQC}, 
                  sphere:{radius:0.35,color:COLORS_3D.HSQC} 
                });
              } catch (e) {
                debugWarn(`[WARN] HSQC setStyle error for atom ${i}:`, e.message);
              }
            }
          }
        }

        // [FIX] Colorier les atomes H1
        for (const i of h1List) {
          if (typeof i === 'number' && i >= 0) {
            let stickColor="#808080", stickRadius=0.15;
            if (flags.HSQC && hsqcSet.has(i)) { 
              stickColor=COLORS_3D.HSQC; 
              stickRadius=0.2; 
            } else if (flags.COSY && cosySet.has(i)) { 
              stickColor=COLORS_3D.COSY; 
              stickRadius=0.25; 
            }
            try {
              model.setStyle({ index:i }, { 
                sphere:{radius:0.3,color:COLORS_3D.H1}, 
                stick:{radius:stickRadius,color:stickColor} 
              });
            } catch (e) {
              debugWarn(`[WARN] H1 setStyle error for atom ${i}:`, e.message);
            }
          }
        }
      }

      // Centered 3D labels for all heteroatoms (≠ C and H)
      try {
        if (viewer.removeAllLabels) {
          viewer.removeAllLabels();
        }

        const atoms = (typeof model.selectedAtoms === "function")
          ? model.selectedAtoms({})
          : [];

        for (const a of atoms) {
          if (!a) continue;
          
          const sym = (a?.elem || "").toUpperCase();
          if (!sym || sym === "C" || sym === "H") continue;
          
          if (!Number.isFinite(a.x) || !Number.isFinite(a.y) || !Number.isFinite(a.z)) {
            continue;
          }

          viewer.addLabel(sym, {
            position: { x: a.x, y: a.y, z: a.z },
            alignment: "center",
            showBackground: false,
            fontSize: 16,
            fontFamily: "sans-serif",
            fontColor: "white",
            inFront: true,
          });
        }
      } catch (e) {
        debugWarn("[WARN] Error adding 3D labels:", e?.message || e);
      }

      viewer.zoomTo(); 
      viewer.render();
      setErr(null);
    } catch (e) {
      const msg = e?.message || String(e);
      console.error("[ERROR] 3D error:", msg);
      setErr(msg);
    } finally {
      setLoading(false);
    }
  }, [multi, computeMultiHighlights, COLORS_3D, flags]);



  const load = useCallback(async () => {
    setErr(null);
    if (!(multi?.fromAll ? multi?.pivotSource : selectedRow?.Source) || !csvUrl) {
      setSvg(null); return;
    }
    const srcName = multi?.fromAll ? multi.pivotSource : selectedRow.Source;
    const resolved = await resolveMolForSource(csvUrl, srcName, multi?.taskId || taskId);
    if (!resolved || !resolved.molblock) { setErr("Molecule not found for this source."); setSvg(null); return; }
    await render2D(resolved);
    if (show3D) await render3D(resolved);
  }, [csvUrl, selectedRow, render2D, render3D, taskId, multi, show3D]);

  useEffect(() => { load(); }, [load]);

  useEffect(() => () => { try { viewerRef.current?.clear?.(); viewerRef.current=null; } catch {} }, []);

  if (compact) {
    return (
      <div style={{ position:"relative", width:"100%", height: molHeight, overflow:"hidden" }}>
        {!err && svg && (
          <div
            style={{ position:"absolute", inset:0, display:"flex", alignItems:"center", justifyContent:"center" }}
            dangerouslySetInnerHTML={{ __html: svg }}
          />
        )}
        {loading && (
          <div style={{ position:"absolute", inset:0, display:"flex", alignItems:"center", justifyContent:"center",
                        color:"#64748b", pointerEvents:"none", fontSize:12 }}>
            Loading…
          </div>
        )}
        {err && (
          <div style={{ position:"absolute", inset:0, display:"flex", alignItems:"center", justifyContent:"center",
                        color:"#dc2626", background:"#fff", fontSize:12, padding:8, textAlign:"center" }}>
            {err}
          </div>
        )}
      </div>
    );
  }

  return (
    <div style={{ border:"1px solid #e5e7eb", borderRadius:8, background:"#fff", overflow:"hidden", padding:8 }}>
      <div style={{ display:"grid", gridTemplateColumns: show3D?"1fr 1fr":"1fr", gap:12, height: molHeight }}>
        {/* 2D */}
        <div style={{ border:"1px solid #e5e7eb", borderRadius:8, background:"#fff", overflow:"hidden", position:"relative" }}>
          <div style={{ fontWeight:600, padding:8, borderBottom:"1px solid #e5e7eb", fontSize:12 }}>2D</div>
          {!err && svg && (
            <div style={{ position:"absolute", top:40,left:0,right:0,bottom:0, display:"flex",alignItems:"center",justifyContent:"center" }}
                dangerouslySetInnerHTML={{ __html: svg }} />
          )}
          {loading && (
            <div style={{ position:"absolute", top:40,left:0,right:0,bottom:0, display:"flex",alignItems:"center",justifyContent:"center",
                          color:"#64748b", pointerEvents:"none" }}>
              Loading 2D…
            </div>
          )}
        </div>
        {/* 3D */}
        {show3D && (
          <div style={{ border:"1px solid #e5e7eb", borderRadius:8, background:"#fff", overflow:"hidden", position:"relative" }}>
            <div style={{ fontWeight:600, padding:8, borderBottom:"1px solid #e5e7eb", fontSize:12 }}>3D</div>
            <div ref={div3DRef} style={{ position:"absolute", top:40,left:0,right:0,bottom:0 }} />
            {loading && (
              <div style={{ position:"absolute", top:40,left:0,right:0,bottom:0, display:"flex",alignItems:"center",justifyContent:"center",
                            color:"#64748b", pointerEvents:"none" }}>
                Loading 3D…
              </div>
            )}
          </div>
        )}
      </div>

      {/* Local toggle bar */}
      {showLocalToggles && (
        <div style={{ marginTop:8, display:"flex", gap:10, flexWrap:"wrap", alignItems:"center", fontSize:12 }}>
          {[
            ["H1", customColors.H1.hex],
            ["C13", customColors.C13.hex],
            ["COSY", customColors.COSY.hex],
            ["HSQC", customColors.HSQC.hex],
          ].map(([k,hex]) => (
            <label key={k}
              style={{ display:"inline-flex", alignItems:"center", gap:6, padding:"4px 8px",
                      border:"1px solid #e5e7eb", borderRadius:999, background: flags[k]?"#f0f7ff":"#fff", cursor:"pointer" }}>
              <input type="checkbox" checked={!!flags[k]} onChange={()=>toggle(k)} style={{ width:14,height:14 }}/>
              <span style={{ width:12,height:12,borderRadius:3,background:hex,border:"1px solid #cbd5e1" }}/>
              <span style={{ fontWeight:600 }}>
                {METHOD_LABELS_SHORT[k] || k}
              </span>
            </label>
          ))}
        </div>
      )}

      {/* [FIX] IMPROVED ERROR DISPLAY */}
      {err && (
        <div style={{ 
          padding:12, 
          marginTop:8,
          color:"#dc2626", 
          background:"#fee2e2", 
          borderRadius:8,
          border:"1px solid #fca5a5",
          fontSize:12 
        }}>
          <div style={{ fontWeight:700, marginBottom:6, fontSize:14, display:"flex", alignItems:"center", gap:6 }}>
            <span style={{ fontSize:20 }}>[!]</span>
            {err.includes("WebGL") ? "WebGL unavailable" : "Rendering error"}
          </div>
          <div style={{ whiteSpace:"pre-wrap", lineHeight:1.5 }}>{err}</div>
          {err.includes("WebGL") && (
            <div style={{ 
              marginTop:12, 
              padding:10,
              background:"#fff",
              borderRadius:6,
              fontSize:11, 
              color:"#991b1b",
              border:"1px solid #fca5a5"
            }}>
              <div style={{ fontWeight:700, marginBottom:6 }}>💡 Quick solutions:</div>
              <ol style={{ margin:0, paddingLeft:20, lineHeight:1.8 }}>
                <li>
                  <strong>Firefox:</strong> Type <code style={{ background:"#fef2f2", padding:"2px 6px", borderRadius:3 }}>about:config</code> in the address bar, 
                  then set <code style={{ background:"#fef2f2", padding:"2px 6px", borderRadius:3 }}>webgl.disabled</code> to <strong>false</strong>
                </li>
                <li>
                  <strong>Graphics drivers:</strong> Update your graphics card drivers (NVIDIA/AMD/Intel)
                </li>
                <li>
                  <strong>Alternative browser:</strong> Try Chrome or Edge if the problem persists
                </li>
              </ol>
            </div>
          )}
        </div>
      )}
    </div>
  );
}


/* ===========================
   Metrics from a CSV row
   =========================== */

const extractMetricsFromRow = (row) => {
  if (!row) return {};
  const metrics = {};
  if (row.coverage_ransy !== undefined && row.coverage_ransy !== "") {
    const v = Number(row.coverage_ransy || 0).toFixed(2) + "%";
    metrics.coverage_ransy = v;
  }
  if (row.coverage_sub !== undefined && row.coverage_sub !== "") {
    const v = Number(row.coverage_sub || 0).toFixed(2) + "%";
    metrics.coverage_sub = v;
  }
  if (row.avg_variabilite !== undefined && row.avg_variabilite !== "") {
    metrics.avg_variabilite = row.avg_variabilite;
  }
  return metrics;
};

/* ===========================
    MODAL — All molecules (visual grid)
   =========================== */

function ShowAllMoleculesModal({
  isOpen,
  onClose,
  sources = [],
  page = 1,
  setPage,
  pageSize = 20,
  setPageSize,
  getRowsByType,
  csvUrl,
  taskId,
  customColors,
}) {
  const stripMethodPrefix = (name) => {
    if (!name) return name;
    // Ex: "1H_MOL", "13C - MOL", "COSY: MOL", "HSQC MOL"
    return String(name).replace(/^(?:1H|13C|COSY|HSQC)[\s_\-:]+/i, "");
  };

  const totalPages = Math.max(1, Math.ceil(sources.length / pageSize));
  const clampedPage = Math.min(Math.max(1, page), totalPages);
  const start = (clampedPage - 1) * pageSize;
  const end = Math.min(sources.length, start + pageSize);
  const pageSources = sources.slice(start, end);

  if (!isOpen) return null;

  return (
    <div
      style={{ position:"fixed", inset:0, zIndex:10000, background:"rgba(0,0,0,0.55)",
               display:"flex", alignItems:"center", justifyContent:"center" }}
      onClick={onClose}
    >
      <div
        style={{ width:"min(1100px, 96vw)", maxHeight:"92vh", background:"#fff", borderRadius:12,
                 boxShadow:"0 24px 80px rgba(0,0,0,.35)", display:"flex", flexDirection:"column" }}
        onClick={(e)=>e.stopPropagation()}
      >
        {/* Header */}
        <div style={{ padding:12, borderBottom:"1px solid #e5e7eb",
                      display:"flex", alignItems:"center", justifyContent:"space-between" }}>
          <div style={{ fontWeight:800, fontSize:14 }}>
            All molecules (sorted by RANSY average)
          </div>
          <div style={{ display:"flex", gap:8, alignItems:"center" }}>
            <label style={{ fontSize:12, color:"#475569" }}>
              Per page:&nbsp;
              <select
                value={pageSize}
                onChange={(e)=>setPageSize(Math.max(10, +e.target.value||20))}
                style={{ padding:"4px 6px", borderRadius:6, border:"1px solid #cbd5e1" }}
              >
                {[20,30,50].map(n=><option key={n} value={n}>{n}</option>)}
              </select>
            </label>
            <button onClick={onClose}
              style={{ padding:"6px 10px", borderRadius:8, border:"1px solid #cbd5e1", background:"#fff", cursor:"pointer", fontWeight:700 }}>
              Close
            </button>
          </div>
        </div>

        {/* Pagination top */}
        <div style={{ padding:"8px 12px", display:"flex", gap:8, alignItems:"center", borderBottom:"1px solid #eef2f7" }}>
          <button onClick={()=>setPage(p=>Math.max(1,p-1))} disabled={clampedPage<=1}
                  style={{ padding:"4px 8px", borderRadius:6, border:"1px solid #cbd5e1", background:"#fff", cursor:"pointer" }}>&lt;</button>
          <div style={{ fontSize:12, color:"#334155" }}>
            Page {clampedPage} / {totalPages} &nbsp;•&nbsp; {sources.length} molecules
          </div>
          <button onClick={()=>setPage(p=>Math.min(totalPages,p+1))} disabled={clampedPage>=totalPages}
                  style={{ padding:"4px 8px", borderRadius:6, border:"1px solid #cbd5e1", background:"#fff", cursor:"pointer" }}></button>
        </div>


        {/* Visual molecules grid */}
        <div style={{ flex:1, overflow:"auto", padding:"12px" }}>
          <div style={{ 
            display:"grid", 
            gridTemplateColumns:"repeat(auto-fill, minmax(180px, 1fr))", 
            gap:12 
          }}>
            {pageSources.map((source, idx) => {
              const globalIdx = start + idx + 1;
              const moleculeName = stripMethodPrefix(source.replace(/_predicted$/i, ""));
              
              return (
                <div key={source} style={{ 
                  border:"1px solid #e5e7eb", 
                  borderRadius:8, 
                  padding:8, 
                  background:"#fff",
                  boxShadow:"0 1px 3px rgba(0,0,0,0.1)",
                  transition:"transform 0.2s ease, box-shadow 0.2s ease",
                  transformOrigin:"center",
                  position:"relative",
                  cursor:"pointer"
                }} 
                onMouseEnter={(e) => {
                  // Important: use currentTarget (the card), not the hovered child
                  e.currentTarget.style.transform = "scale(1.02)";
                  e.currentTarget.style.boxShadow = "0 4px 12px rgba(0,0,0,0.15)";
                  e.currentTarget.style.zIndex = "2";
                }}
                onMouseLeave={(e) => {
                  e.currentTarget.style.transform = "scale(1)";
                  e.currentTarget.style.boxShadow = "0 1px 3px rgba(0,0,0,0.1)";
                  e.currentTarget.style.zIndex = "0";
                }}>
                  
                  {/* Real molecule with MoleculePane */}
                  <div style={{ height:140, marginBottom:8 }}>
                    <MoleculePane
                      csvUrl={csvUrl}
                      selectedRow={{ Source: source }}
                      show3D={false}
                      taskId={taskId}
                      customColors={customColors}
                      molHeight={140}
                      showLocalToggles={false}
                      compact={true}
                      style={{ height:"100%" }}
                      multi={{
                        fromAll: true,
                        pivotSource: source,
                        taskId,
                        rowsByType: getRowsByType(source),
                      }}
                    />
                  </div>
                  
                  {/* Molecule name */}
                  <div style={{ 
                    fontWeight:600, 
                    fontSize:11, 
                    marginBottom:6, 
                    wordBreak:"break-word",
                    textAlign:"center"
                  }}>
                    #{globalIdx} {moleculeName}
                  </div>
                </div>
              );
            })}
          </div>
        </div>

        {/* Pagination bottom */}
        <div style={{ padding:"8px 12px", borderTop:"1px solid #e5e7eb",
                      display:"flex", gap:8, alignItems:"center", justifyContent:"center" }}>
          <div style={{ display:"flex", gap:8, alignItems:"center" }}>
            <button onClick={()=>setPage(p=>Math.max(1,p-1))} disabled={clampedPage<=1}
                    style={{ padding:"6px 10px", borderRadius:8, border:"1px solid #cbd5e1", background:"#fff", cursor:"pointer" }}>&lt; Previous</button>
            <div style={{ fontSize:12, color:"#334155" }}>Page {clampedPage} / {totalPages}</div>
            <button onClick={()=>setPage(p=>Math.min(totalPages,p+1))} disabled={clampedPage>=totalPages}
                    style={{ padding:"6px 10px", borderRadius:8, border:"1px solid #cbd5e1", background:"#fff", cursor:"pointer" }}>Next </button>
          </div>
        </div>
      </div>
    </div>
  );
}



/* ===========================
   Main exported Viewer with SAFETY TIMEOUTS
   =========================== */

export default function Viewer({
  csvs,
  taskId = null,
  defaultShow3D = true,
  annotationResults = {},
  pdataDirs = {},
  bgOptions = {},
  customColors: propsCustomColors,
  onColorChange,
  // NEW: Props for external source control (SIRIUS-style)
  externalSource = null,
  hideSidebar = false,
}) {
  // Use colors passed as props or default values
  const customColors = propsCustomColors || DEFAULT_COLORS;

  const setCustomColor = useCallback((method, color) => {
    if (onColorChange) {
      onColorChange(method, color);
    }
  }, [onColorChange]);

  const [, setRowsByType] = useState({});
  const [mapsByType, setMapsByType] = useState({});
  const [loading, setLoading] = useState(true);
  const [err, setErr] = useState(null);

  // If externalSource is provided, use it; otherwise keep internal state
  const [internalSource, setInternalSource] = useState(null);
  const selectedSource = externalSource || internalSource;
  const setSelectedSource = useCallback(
    (next) => {
      if (externalSource) return;
      setInternalSource(next);
    },
    [externalSource]
  );
  
  const [query, setQuery] = useState("");

  const [showTheory, setShowTheory] = useState(true);
  const [showExp, setShowExp] = useState(true);
  const [reverse1D, setReverse1D] = useState(true);
  const [show3D, setShow3D] = useState(defaultShow3D);

  const [zoomRectEnabled, setZoomRectEnabled] = useState(true);
  const [colorPickerOpen, setColorPickerOpen] = useState(false);
  const [relevel2DEnabled, setRelevel2DEnabled] = useState(false); // Disabled by default - only zoom updates background

  // NEW: nonce par type pour forcer un remount des canvases
  const [refreshKeyByType, setRefreshKeyByType] = useState({});
  const bumpKey = useCallback((type) => {
    setRefreshKeyByType(k => ({ ...k, [type]: (k[type] ?? 0) + 1 }));
  }, []);

  const wheelRequestRef = useRef({});
  
  const presentTypes = useMemo(
    () => Object.keys(csvs || {}).filter(Boolean),
    [csvs]
  );

  const hasAnyData = useMemo(() => presentTypes.length > 0, [presentTypes]);

  const [bgByType, setBgByType] = useState(
    () => annotationResults || {}
  );
  const [busyByType, setBusyByType] = useState({});

  // Local state for resolution and level parameters
  const [resByType, setResByType] = useState({
    "1H":  { width_px: 1200, height_px: 500, dpi: 150 },
    "13C": { width_px: 1200, height_px: 500, dpi: 150 },
    COSY:  { width_px: 1000, height_px: 800, dpi: 150 },
    HSQC:  { width_px: 1000, height_px: 800, dpi: 150 },
  });

  const [levelByType, setLevelByType] = useState({
    COSY: 0.005,
    HSQC: 0.005
  });

  // 🔁 AVANT : const originalBgRef = useRef({});
  // [FIX] AFTER: structure by type AND by source: { [type]: { [source]: snap } }
  const originalBgRef = useRef({});

  useEffect(() => {
    const snap = {};
    for (const k of ["1H", "13C", "COSY", "HSQC"]) {
      if (annotationResults?.[k]) snap[k] = annotationResults[k];
    }
    setBgByType(prev => ({ ...prev, ...snap }));
  }, [annotationResults]);

  // When initial results arrive, snap according to current source
  useEffect(() => {
    if (!selectedSource) return;
    for (const k of ["1H", "13C", "COSY", "HSQC"]) {
      if (annotationResults?.[k]) {
        originalBgRef.current[k] = originalBgRef.current[k] || {};
        if (!originalBgRef.current[k][selectedSource]) {
          originalBgRef.current[k][selectedSource] = annotationResults[k];
        }
      }
    }
  }, [annotationResults, selectedSource]);

  // When source changes, if we have a background in memory, snap it too
  useEffect(() => {
    if (!selectedSource) return;
    for (const k of ["1H", "13C", "COSY", "HSQC"]) {
      const bg = bgByType?.[k];
      if (bg) {
        originalBgRef.current[k] = originalBgRef.current[k] || {};
        if (!originalBgRef.current[k][selectedSource]) {
          originalBgRef.current[k][selectedSource] = bg;
        }
      }
    }
  }, [selectedSource, bgByType]);

  const loadAll = useCallback(async () => {
    if (!csvs) return;
    setLoading(true);
    setErr(null);
    try {
      const entries = await Promise.all(
        presentTypes.map(async (t) => {
          try {
            const r = await fetchWithRetry(csvs[t], 3, 500);
            const text = await r.text();
            const { data } = parseCSVSmart(text);
            const map = new Map();
            for (const row of data || []) {
              if (row?.Source) map.set(row.Source, row);
            }
            return [t, data || [], map, null];
          } catch (error) {
            return [t, [], new Map(), error.message];
          }
        })
      );

      const rb = {};
      const mb = {};
      const errors = {};
      for (const [t, rows, map, error] of entries) {
        rb[t] = rows;
        mb[t] = map;
        if (error) errors[t] = error;
      }

      setRowsByType(rb);
      setMapsByType(mb);

      const failed = Object.keys(errors);
      if (failed.length > 0) {
        setErr(`Spectra unavailable: ${failed.join(", ")}`);
      }

      const allSets = Object.values(mb).map(m => new Set(m.keys()));
      let common = [];
      if (allSets.length) {
        common = [...allSets[0]].filter(s =>
          allSets.every(S => S.has(s))
        );
        if (common.length === 0) {
          const union = new Set();
          Object.values(mb).forEach(m => {
            for (const s of m.keys()) union.add(s);
          });
          common = [...union];
        }
      }
      common.sort();
      setSelectedSource(common[0] || null);
    } catch (e) {
      setErr(e.message || "Error loading CSV files");
    } finally {
      setLoading(false);
    }
  }, [csvs, presentTypes, setSelectedSource]);

  useEffect(() => { loadAll(); }, [loadAll]);

  const allSources = useMemo(() => {
    const union = new Set();
    Object.values(mapsByType).forEach(m => {
      for (const s of m.keys()) union.add(s);
    });
    let arr = [...union];

    const parseCoverage = (row) => {
      if (!row) return null;
      // Accept a few possible field names (CSV variations)
      const candidates = [
        row.coverage_ransy,
        row.coverage_ransy_mean,
        row.coverage_ransy_avg,
        row.coverage_ransy_percent,
      ];
      for (const c of candidates) {
        if (c === undefined || c === null || c === "") continue;
        const v = parseFloat(String(c).replace(",", "."));
        if (Number.isFinite(v)) return v;
      }
      return null;
    };

    const withCov = arr.map(source => {
      const coverages = [];
      for (const type of ["1H", "13C", "COSY", "HSQC"]) {
        const row = mapsByType[type]?.get(source);
        const cov = parseCoverage(row);
        if (cov != null) coverages.push(cov);
      }
      const avg = coverages.length > 0
        ? coverages.reduce((s, v) => s + v, 0) / coverages.length
        : -1;
      return { source, avgCoverage: avg };
    });

    withCov.sort((a, b) => (b.avgCoverage - a.avgCoverage) || a.source.localeCompare(b.source));
    arr = withCov.map(x => x.source);

    if (!query) return arr;
    const q = query.toLowerCase();
    return arr.filter(s => s.toLowerCase().includes(q));
  }, [mapsByType, query]);

  // 🔥 Global sorted list (without search filter) for the modal
  const sortedSourcesAll = useMemo(() => {
    const union = new Set();
    Object.values(mapsByType).forEach(m => {
      for (const s of m.keys()) union.add(s);
    });
    const arr = [...union];

    const parseCoverage = (row) => {
      if (!row) return null;
      const candidates = [
        row.coverage_ransy,
        row.coverage_ransy_mean,
        row.coverage_ransy_avg,
        row.coverage_ransy_percent,
      ];
      for (const c of candidates) {
        if (c === undefined || c === null || c === "") continue;
        const v = parseFloat(String(c).replace(",", "."));
        if (Number.isFinite(v)) return v;
      }
      return null;
    };

    const withCov = arr.map(source => {
      const covs = [];
      for (const type of ["1H", "13C", "COSY", "HSQC"]) {
        const row = mapsByType[type]?.get(source);
        const v = parseCoverage(row);
        if (v != null) covs.push(v);
      }
      const avg = covs.length ? (covs.reduce((s,v)=>s+v,0)/covs.length) : -1;
      return { source, avgCoverage: avg };
    });

    withCov.sort((a,b) => (b.avgCoverage - a.avgCoverage) || a.source.localeCompare(b.source));
    return withCov.map(x => x.source);
  }, [mapsByType]);

  // Helper lines by method for a source (modal)
  const rowsByTypeForSource = useCallback((source) => ({
    "1H":   mapsByType["1H"]?.get(source),
    "13C":  mapsByType["13C"]?.get(source),
    "COSY": mapsByType["COSY"]?.get(source),
    "HSQC": mapsByType["HSQC"]?.get(source),
  }), [mapsByType]);

  const rowFor = useCallback(
    (t) => mapsByType[t]?.get(selectedSource) || null,
    [mapsByType, selectedSource]
  );

  const oneH = useMemo(() => {
    const r = rowFor("1H");
    return r ? {
      theory: (r.ppm_molecule || "")
        .split("|")
        .map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      expX: (r.RANSY_signals || "")
        .split("|")
        .map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      expI: (r.RANSY_mean_intensities || r.intensity_ratios || "")
        .split(";")
        .map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      matches: (r.ppm_molecule_matched || "")
        .split("|")
        .map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      metrics: extractMetricsFromRow(r),
      smiles: r.smiles || r.SMILES || r.Smiles || r.smi || "",  // SMILES for spectrum simulation
      source: r.Source || selectedSource || "",  // Source name
    } : null;
  }, [rowFor, selectedSource]);

  const c13 = useMemo(() => {
    const r = rowFor("13C");
    return r ? {
      theory: (r.ppm_molecule || "")
        .split("|").map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      expX: (r.RANSY_signals || "")
        .split("|").map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      expI: (r.RANSY_mean_intensities || r.intensity_ratios || "")
        .split(";").map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      matches: (r.ppm_molecule_matched || "")
        .split("|").map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      metrics: extractMetricsFromRow(r),
    } : null;
  }, [rowFor]);

  const cosy = useMemo(() => {
    const r = rowFor("COSY");
    return r ? {
      theory: (r.ppm_molecule || "")
        .split("||")
        .map(p => p.split("|")
        .map(v => parseFloat(String(v).replace(",", "."))))
        .filter(a => a.length === 2 && a.every(Number.isFinite)),
      expPts: (r.RANSY_signals || "")
        .split("||")
        .map(p => p.split("|")
        .map(v => parseFloat(String(v).replace(",", "."))))
        .filter(a => a.length === 2 && a.every(Number.isFinite)),
      expI: (r.RANSY_mean_intensities || r.intensity_ratios || "")
        .split(";")
        .map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      pairs: (r.ppm_molecule_matched || "")
        .split("||")
        .map(p => p.split("|")
        .map(v => parseFloat(String(v).replace(",", "."))))
        .filter(a => a.length === 2 && a.every(Number.isFinite)),
      metrics: extractMetricsFromRow(r),
    } : null;
  }, [rowFor]);

  const hsqc = useMemo(() => {
    const r = rowFor("HSQC");
    return r ? {
      theory: (r.ppm_molecule || "")
        .split("||")
        .map(p => p.split("|")
        .map(v => parseFloat(String(v).replace(",", "."))))
        .filter(a => a.length === 2 && a.every(Number.isFinite)),
      expPts: (r.RANSY_signals || "")
        .split("||")
        .map(p => p.split("|")
        .map(v => parseFloat(String(v).replace(",", "."))))
        .filter(a => a.length === 2 && a.every(Number.isFinite)),
      expI: (r.RANSY_mean_intensities || r.intensity_ratios || "")
        .split(";")
        .map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      pairs: (
        r.ppm_molecule_matched_unfolded ||
        r.ppm_molecule_matched || ""
      )
        .split("||")
        .map(p => p.split("|")
        .map(v => parseFloat(String(v).replace(",", "."))))
        .filter(a => a.length === 2 && a.every(Number.isFinite)),
      metrics: extractMetricsFromRow(r),
    } : null;
  }, [rowFor]);

  const pivotCsvUrl = useMemo(
    () => csvs["1H"] || csvs["13C"] || csvs["COSY"] || csvs["HSQC"] || null,
    [csvs]
  );

  const getMeta = (obj) => obj?.background?.meta || obj?.meta || null;

  // eslint-disable-next-line react-hooks/exhaustive-deps
  const setBusy = useCallback((type, v) => setBusyByType(prev => ({ ...prev, [type]: !!v })));
  
  // [FIX] storeResult fix: adds anti-cache nonce and updates snapshot (type, source)
  // eslint-disable-next-line react-hooks/exhaustive-deps
  const storeResult = useCallback((type, res) => {
    if (!res) return;

    // nonce anti-cache
    const nonce = `${Date.now()}-${Math.random().toString(36).slice(2)}`;
    res.meta = { ...(res.meta || {}), nonce };

    // updates current state (display)
    setBgByType(prev => ({ ...prev, [type]: res }));

    // si 2D, on synchronise le lvl UI
    if ((type === "COSY" || type === "HSQC") &&
        Number.isFinite(+res?.meta?.level_start)) {
      setLevelByType(p => ({ ...p, [type]: +res.meta.level_start }));
    }

    // 🔒 initial snapshot by (type, source): only write if absent
    if (selectedSource) {
      originalBgRef.current[type] = originalBgRef.current[type] || {};
      if (!originalBgRef.current[type][selectedSource]) {
        originalBgRef.current[type][selectedSource] = res;
      }
    }
  });

  const extractWindow1D = (meta) => {
    if (!meta) return null;
    if (meta.ppm_min != null && meta.ppm_max != null) {
      return { ppmMin: +meta.ppm_min, ppmMax: +meta.ppm_max };
    }
    if (meta.ppm_left != null && meta.ppm_right != null) {
      const left = +meta.ppm_left;
      const right = +meta.ppm_right;
      return {
        ppmMin: Math.min(left, right),
        ppmMax: Math.max(left, right),
      };
    }
    return null;
  };

  const extractWindow2D = (meta) => {
    if (!meta) return null;
    let f2Min, f2Max, f1Min, f1Max;
    if (meta.f2_ppm_min != null && meta.f2_ppm_max != null) {
      f2Min = +meta.f2_ppm_min; f2Max = +meta.f2_ppm_max;
    } else if (meta.f2_left != null && meta.f2_right != null) {
      const L = +meta.f2_left; const R = +meta.f2_right;
      f2Min = Math.min(L, R); f2Max = Math.max(L, R);
    }
    if (meta.f1_ppm_min != null && meta.f1_ppm_max != null) {
      f1Min = +meta.f1_ppm_min; f1Max = +meta.f1_ppm_max;
    } else if (meta.f1_top != null && meta.f1_bottom != null) {
      const T = +meta.f1_top; const B = +meta.f1_bottom;
      f1Min = Math.min(T, B); f1Max = Math.max(T, B);
    }
    if ([f2Min,f2Max,f1Min,f1Max].every(v => Number.isFinite(v))) {
      return { f2Min, f2Max, f1Min, f1Max };
    }
    return null;
  };

  
  // FIX: Added safety timeout for applySettings
  const applySettings = useCallback(async (type) => {
    const pd = pdataDirs?.[type];
    if (!pd) return;
    
    try {
      setBusy(type, true);
      
      // Safety timeout
      const timeoutId = setTimeout(() => {
        console.warn(`Timeout for applySettings ${type}, forced unlock`);
        setBusy(type, false);
      }, 30000);
      
      const size = resByType[type] || {};
      const meta = getMeta(bgByType?.[type]);
      let res;
      
      if (type === "1H" || type === "13C") {
        const win = extractWindow1D(meta) || null;
        const line = { width: 1.5 };
        const fn = type === "1H" ? interactive1H : interactive13C;
        if (win) {
          const payload = {
            ppm_min: Math.min(win.ppmMin, win.ppmMax),
            ppm_max: Math.max(win.ppmMin, win.ppmMax),
            ppmMin: Math.min(win.ppmMin, win.ppmMax),
            ppmMax: Math.max(win.ppmMin, win.ppmMax),
            size, line, taskId,
          };
          res = await fn(pd, payload);
        } else {
          res = await fn(pd, { size, line, taskId });
        }
      } else {
        const win = extractWindow2D(meta) || {};
        const contour = {
          level_start: levelByType[type] ?? 0.005,
          n_levels: 20,
          transparent: true,
        };
        const fn = type === "COSY" ? interactiveCOSY : interactiveHSQC;
        res = await fn(pd, { ...win, size, contour, taskId });
      }
      
      clearTimeout(timeoutId);
      storeResult(type, res);
    } catch (e) {
      console.error(`applySettings ${type}:`, e);
    } finally {
      setBusy(type, false);
    }
  }, [pdataDirs, setBusy, resByType, bgByType, storeResult, taskId, levelByType]);

  // FIX: Added safety timeout for doZoom1DPixels
  const doZoom1DPixels = useCallback(async (type, rectPxDOM, display) => {
    const pd = pdataDirs?.[type];
    const meta = getMeta(bgByType?.[type]);
    if (!pd || !meta || !rectPxDOM || !(rectPxDOM.w > 1)) return;

    try {
      setBusy(type, true);
      
      // Safety timeout
      const timeoutId = setTimeout(() => {
        console.warn(`Timeout for zoom1DPixels ${type}, forced unlock`);
        setBusy(type, false);
      }, 30000);
      
      const size = resByType[type] || {};
      const line = { width: 1.5 };
      const { displayWidth, displayHeight } = display || {};
      const res = await zoom1DFromPixels(
        type, pd, meta, rectPxDOM,
        { taskId, size, line, displayWidth, displayHeight }
      );
      
      clearTimeout(timeoutId);
      storeResult(type, res);
    } catch (e) {
      console.error(`ZoomPixels ${type} error:`, e);
    } finally {
      setBusy(type, false);
    }
  }, [pdataDirs, bgByType, setBusy, resByType, taskId, storeResult]);

  // FIX: Added safety timeout for doZoom1D
  const doZoom1D = useCallback(async (type, ppmMin, ppmMax) => {
    const pd = pdataDirs?.[type];

    if (!pd) {
      setBgByType(prev => ({
        ...prev,
        [type]: {
          ...prev[type],
          meta: {
            ...prev[type]?.meta,
            ppm_min: Math.min(ppmMin, ppmMax),
            ppm_max: Math.max(ppmMin, ppmMax),
            ppm_left: Math.max(ppmMin, ppmMax),
            ppm_right: Math.min(ppmMin, ppmMax),
            nonce: `${Date.now()}`
          },
        },
      }));
      bumpKey(type);
      return;
    }

    try {
      setBusy(type, true);
      
      // Safety timeout
      const timeoutId = setTimeout(() => {
        console.warn(`Timeout for zoom1D ${type}, forced unlock`);
        setBusy(type, false);
      }, 30000);
      
      const size = resByType[type] || {};
      const line = { width: 1.5 };
      const fn = type === "1H" ? interactive1H : interactive13C;
      const payload = {
        ppm_min: Math.min(ppmMin, ppmMax),
        ppm_max: Math.max(ppmMin, ppmMax),
        ppmMin: Math.min(ppmMin, ppmMax),
        ppmMax: Math.max(ppmMin, ppmMax),
        size, line, taskId,
      };
      const res = await fn(pd, payload);
      
      clearTimeout(timeoutId);
      storeResult(type, res);
    } catch (e) {
      console.error(`Zoom ${type} error:`, e);
    } finally {
      setBusy(type, false);
    }
  }, [pdataDirs, bumpKey, setBusy, resByType, taskId, storeResult]);

  // *** FIXED *** : doReset1D prefers backend recalculation if possible, otherwise restores snapshot
  const doReset1D = useCallback(async (type) => {
    const pd = pdataDirs?.[type];
    const snap = originalBgRef.current?.[type]?.[selectedSource];

    // [FIX] Retrieve USER-DEFINED bounds
    const userWindow = bgOptions?.[type] || {};
    const hasUserWindow = 
      (type === "1H" || type === "13C") &&
      userWindow.ppm_min != null && 
      userWindow.ppm_max != null;

    // If backend interaction available → request image with defined bounds
    if (pd) {
      try {
        setBusy(type, true);
        const timeoutId = setTimeout(() => {
          console.warn(`Timeout for reset1D ${type}, forced unlock`);
          setBusy(type, false);
        }, 30000);

        const size = resByType[type] || {};
        const line = { width: 1.5 };
        const fn = type === "1H" ? interactive1H : interactive13C;
        
        // [FIX] USE USER BOUNDS if defined
        const payload = { size, line, taskId };
        
        if (hasUserWindow) {
          // Bornes utilisateur
          payload.ppm_min = Math.min(userWindow.ppm_min, userWindow.ppm_max);
          payload.ppm_max = Math.max(userWindow.ppm_min, userWindow.ppm_max);
          payload.ppmMin = payload.ppm_min;
          payload.ppmMax = payload.ppm_max;
        }
        // Sinon, pas de bornes = plein champ natif du backend

        const res = await fn(pd, payload);

        clearTimeout(timeoutId);
        storeResult(type, res);
        bumpKey(type); // force remount canvas
        return;
      } catch (e) {
        console.error(`Reset ${type} via backend error:`, e);
      } finally {
        setBusy(type, false);
      }
    }

    // Otherwise, fallback to local snapshot state
    if (snap) {
      const restored = {
        ...snap,
        meta: { ...snap.meta, nonce: `${Date.now()}` },
      };
      setBgByType(prev => ({ ...prev, [type]: restored }));
      originalBgRef.current[type] = originalBgRef.current[type] || {};
      originalBgRef.current[type][selectedSource] = restored;
      bumpKey(type); // force remount
    }
  }, [pdataDirs, selectedSource, bgOptions, setBusy, resByType, taskId, storeResult, bumpKey]);

  // FIX: Added safety timeout for doZoom2D
  const doZoom2D = useCallback(async (type, { f2_min, f2_max, f1_min, f1_max }) => {
    const pd = pdataDirs?.[type];
    if (!pd) return;
    
    try {
      setBusy(type, true);
      
      // Safety timeout
      const timeoutId = setTimeout(() => {
        console.warn(`Timeout for zoom2D ${type}, forced unlock`);
        setBusy(type, false);
      }, 30000);
      
      const size = resByType[type] || {};
      const contour = {
        level_start: levelByType[type] ?? 0.005,
        n_levels: 20,
        transparent: true,
      };
      const fn = type === "COSY" ? interactiveCOSY : interactiveHSQC;
      const res = await fn(
        pd,
        {
          f2Min: f2_min, f2Max: f2_max,
          f1Min: f1_min, f1Max: f1_max,
          size, contour, taskId,
        }
      );
      
      clearTimeout(timeoutId);
      storeResult(type, res);
    } catch (e) {
      console.error(`Zoom ${type} error:`, e);
    } finally {
      setBusy(type, false);
    }
  }, [pdataDirs, setBusy, resByType, levelByType, taskId, storeResult]);


  const doWheel2D = useCallback(async (type, deltaY) => {
    // If 2D relevel is disabled, do nothing
    if (!relevel2DEnabled) return;
    
    const pd = pdataDirs?.[type];
    const meta = getMeta(bgByType?.[type]);
    if (!pd || !meta) return;
    
    if (wheelRequestRef.current[type]) {
      try {
        wheelRequestRef.current[type].cancel();
      } catch {}
    }
    
    try {
      setBusy(type, true);
      
      const timeoutId = setTimeout(() => {
        console.warn(`Timeout wheel2D ${type}, annulation`);
        wheelRequestRef.current[type]?.cancel();
        setBusy(type, false);
      }, 10000);
      
      // [FIX] CORRECTION : toujours envoyer size
      const request = relevelByWheelCancelable(
        type, pd, meta, deltaY,
        { 
          nLevels: 20, 
          transparent: true, 
          taskId,
          size: resByType[type] || {} // ← AJOUT
        }
      );
      
      wheelRequestRef.current[type] = request;
      const res = await request;
      
      clearTimeout(timeoutId);
      storeResult(type, res);
      wheelRequestRef.current[type] = null;
    } catch (e) {
      if (e?.name === 'AbortError') {
        console.log(`Request wheel2D ${type} cancelled`);
        return;
      }
      
      console.error(`Wheel relevel ${type}:`, e);
      
      // Fallback local (garder tel quel)
      if (e?.response?.status === 500 || e?.code === "ERR_BAD_RESPONSE") {
        console.warn(`Backend failed, local level adjustment for ${type}`);
        
        const currentLevel = meta?.level_start ?? levelByType[type] ?? 0.005;
        const direction = deltaY > 0 ? 1.2 : 0.8;
        const newLevel = Math.max(1e-6, currentLevel * direction);
        
        setLevelByType(p => ({ ...p, [type]: newLevel }));
        
        if (bgByType[type]?.meta) {
          setBgByType(prev => ({
            ...prev,
            [type]: {
              ...prev[type],
              meta: {
                ...prev[type].meta,
                level_start: newLevel,
                nonce: `${Date.now()}`
              },
            },
          }));
          bumpKey(type);
        }
      }
    } finally {
      setBusy(type, false);
    }
  }, [pdataDirs, bgByType, setBusy, taskId, storeResult, levelByType, bumpKey, resByType, relevel2DEnabled]);
  useEffect(() => {
    return () => {
      Object.values(wheelRequestRef.current).forEach(request => {
        try {
          request?.cancel();
        } catch {}
      });
      wheelRequestRef.current = {};
    };
  }, []);
    
  // Reliable 2D reset: keeps current background while waiting for backend
  // Instant 2D reset + conditional lightweight recalculation in background
  const doReset2D = useCallback(async (type) => {
    const snap = originalBgRef.current?.[type]?.[selectedSource];
    const pd   = pdataDirs?.[type];

    // [FIX] Retrieve USER-DEFINED bounds
    const userWindow = bgOptions?.[type] || {};
    const hasUserWindow = 
      (type === "COSY" || type === "HSQC") &&
      [userWindow.f1_ppm_min, userWindow.f1_ppm_max,
      userWindow.f2_ppm_min, userWindow.f2_ppm_max].every(v => v != null);

    // Starting level
    const currentLevelStart =
      (snap?.meta?.level_start) ??
      bgByType[type]?.meta?.level_start ??
      levelByType[type] ?? 0.005;

    // Instant snapshot display
    if (snap) {
      const restored = {
        ...snap,
        meta: {
          ...snap.meta,
          level_start: currentLevelStart,
          nonce: `${Date.now()}`,
        },
      };
      setBgByType(prev => ({ ...prev, [type]: restored }));
      bumpKey(type);
    }

    // If no backend, stop here
    if (!pd) return;

    setBusy(type, true);

    try {
      const timeoutId = setTimeout(() => {
        console.warn(`Timeout for reset2D ${type}, forced unlock`);
        setBusy(type, false);
      }, 30000);

      const lightSize = {
        width_px:  Math.min(900,  (resByType[type]?.width_px  ?? 1000)),
        height_px: Math.min(700,  (resByType[type]?.height_px ?? 800)),
        dpi:       Math.min(120,  (resByType[type]?.dpi       ?? 150)),
      };
      const lightContour = {
        level_start: currentLevelStart,
        n_levels: 12,
        transparent: true,
      };

      const fn = type === "COSY" ? interactiveCOSY : interactiveHSQC;

      // [FIX] UTILISER LES BORNES UTILISATEUR au lieu de "full: true"
      const payload = {
        size: lightSize,
        contour: lightContour,
        taskId,
        with_bbox: true,
        return_bbox: true,
      };

      // [FIX] If user has defined bounds, use them
      if (hasUserWindow) {
        payload.f1Min = userWindow.f1_ppm_min;
        payload.f1Max = userWindow.f1_ppm_max;
        payload.f2Min = userWindow.f2_ppm_min;
        payload.f2Max = userWindow.f2_ppm_max;
      } else {
        // Sinon, plein champ
        payload.reset = true;
        payload.full = true;
      }

      const resLight = await fn(pd, payload);

      clearTimeout(timeoutId);

      // Normaliser bbox
      if (resLight?.meta && !resLight.meta.plot_bbox_px) {
        const bb =
          resLight.meta.axis_bbox_px ||
          resLight.meta.axisRect?.pixel ||
          resLight.meta.axisRect?.px || null;
        if (bb && [bb.x, bb.y, bb.w, bb.h].every(Number.isFinite)) {
          resLight.meta.plot_bbox_px = { x: bb.x, y: bb.y, w: bb.w, h: bb.h };
        }
      }

      storeResult(type, resLight);
      bumpKey(type);

      // HD quality in background if needed
      const targetSize = resByType[type] || {};
      const wantQuality =
        // eslint-disable-next-line no-mixed-operators
        targetSize.width_px  && targetSize.width_px  > lightSize.width_px ||
        // eslint-disable-next-line no-mixed-operators
        targetSize.height_px && targetSize.height_px > lightSize.height_px ||
        (resByType[type]?.dpi ?? 150) > lightSize.dpi;

      if (wantQuality) {
        const qualityContour = {
          level_start: currentLevelStart,
          n_levels: 20,
          transparent: true,
        };
        
        const payloadHD = {
          size: resByType[type] || {},
          contour: qualityContour,
          taskId,
          with_bbox: true,
          return_bbox: true,
        };

        // [FIX] Même logique pour la qualité HD
        if (hasUserWindow) {
          payloadHD.f1Min = userWindow.f1_ppm_min;
          payloadHD.f1Max = userWindow.f1_ppm_max;
          payloadHD.f2Min = userWindow.f2_ppm_min;
          payloadHD.f2Max = userWindow.f2_ppm_max;
        } else {
          payloadHD.reset = true;
          payloadHD.full = true;
        }

        fn(pd, payloadHD)
          .then((resHi) => {
            if (resHi?.meta && !resHi.meta.plot_bbox_px) {
              const bb =
                resHi.meta.axis_bbox_px ||
                resHi.meta.axisRect?.pixel ||
                resHi.meta.axisRect?.px || null;
              if (bb && [bb.x, bb.y, bb.w, bb.h].every(Number.isFinite)) {
                resHi.meta.plot_bbox_px = { x: bb.x, y: bb.y, w: bb.w, h: bb.h };
              }
            }
            storeResult(type, resHi);
            bumpKey(type);
          })
          .catch((e) => {
            console.warn(`Qualité ${type} échouée:`, e?.message || e);
          })
          .finally(() => setBusy(type, false));
      } else {
        setBusy(type, false);
      }
    } catch (e) {
      console.error(`Reset ${type} error:`, e);
      setBusy(type, false);
    }
  }, [
    pdataDirs,
    resByType,
    taskId,
    levelByType,
    bgByType,
    bgOptions,  // [FIX] AJOUT dépendance
    bumpKey,
    selectedSource,
    storeResult,
    setBusy,
  ]);

  const canInteract1H   = !!pdataDirs?.["1H"];
  const canInteract13C  = !!pdataDirs?.["13C"];
  const canInteractCOSY = !!pdataDirs?.["COSY"];
  const canInteractHSQC = !!pdataDirs?.["HSQC"];

  // État pour la modale "Afficher toutes les molécules"
  const [showAllOpen, setShowAllOpen] = useState(false);
  const [allPage, setAllPage] = useState(1);
  const [allPageSize, setAllPageSize] = useState(20);

  // État pour le modal de spectre 1H détaillé
  const [spectrumDetailOpen, setSpectrumDetailOpen] = useState(false);

  /* petits panneaux de réglage */
  const renderControls1D = (type) => (
    <div style={{
      display: "flex", gap: 6, alignItems: "center",
      flexWrap: "wrap", marginBottom: 6, fontSize: 11,
    }}>
      <label style={{ fontSize: 10 }}>
        W:
        <input
          type="number" min={300}
          value={resByType[type].width_px}
          onChange={(e) =>
            setResByType(p => ({
              ...p,
              [type]: {
                ...p[type],
                width_px: Math.max(300, +e.target.value || 300),
              },
            }))
          }
          style={{ width: 70, fontSize: 10, padding: "1px 3px" }}
        />
      </label>
      <label style={{ fontSize: 10 }}>
        H:
        <input
          type="number" min={200}
          value={resByType[type].height_px}
          onChange={(e) =>
            setResByType(p => ({
              ...p,
              [type]: {
                ...p[type],
                height_px: Math.max(200, +e.target.value || 200),
              },
            }))
          }
          style={{ width: 70, fontSize: 10, padding: "1px 3px" }}
        />
      </label>
      <label style={{ fontSize: 10 }}>
        DPI:
        <input
          type="number" min={72}
          value={resByType[type].dpi}
          onChange={(e) =>
            setResByType(p => ({
              ...p,
              [type]: {
                ...p[type],
                dpi: Math.max(72, +e.target.value || 72),
              },
            }))
          }
          style={{ width: 60, fontSize: 10, padding: "1px 3px" }}
        />
      </label>
      <button
        onClick={() => applySettings(type)}
        style={{
          padding: "2px 6px", borderRadius: 4,
          border: "1px solid #cbd5e1",
          background: "#fff", cursor: "pointer", fontSize: 10,
        }}
        disabled={!pdataDirs?.[type] || !!busyByType?.[type]}
      >
        Appliquer
      </button>
    </div>
  );

  const renderControls2D = (type) => {
    const currentLevel = bgByType[type]?.meta?.level_start ?? levelByType[type] ?? 0.005;
    
    return (
      <div style={{
        display: "flex", gap: 6, alignItems: "center",
        flexWrap: "wrap", marginBottom: 6, fontSize: 11,
      }}>
        <label style={{ fontSize: 10 }}>
          W:
          <input
            type="number" min={400}
            value={resByType[type].width_px}
            onChange={(e) =>
              setResByType(p => ({
                ...p,
                [type]: {
                  ...p[type],
                  width_px: Math.max(400, +e.target.value || 400),
                },
              }))
            }
            style={{ width: 70, fontSize: 10, padding: "1px 3px" }}
          />
        </label>
        <label style={{ fontSize: 10 }}>
          H:
          <input
            type="number" min={300}
            value={resByType[type].height_px}
            onChange={(e) =>
              setResByType(p => ({
                ...p,
                [type]: {
                  ...p[type],
                  height_px: Math.max(300, +e.target.value || 300),
                },
              }))
            }
            style={{ width: 70, fontSize: 10, padding: "1px 3px" }}
          />
        </label>
        <label style={{ fontSize: 10 }}>
          DPI:
          <input
            type="number" min={72}
            value={resByType[type].dpi}
            onChange={(e) =>
              setResByType(p => ({
                ...p,
                [type]: {
                  ...p[type],
                  dpi: Math.max(72, +e.target.value || 72),
                },
              }))
            }
            style={{ width: 60, fontSize: 10, padding: "1px 3px" }}
          />
        </label>
        <label style={{ fontSize: 10 }}>
          Lvl:
          <input
            type="number" step="0.001" min="0.000001"
            value={currentLevel}
            onChange={(e) => {
              const newVal = Math.max(1e-6, +e.target.value || 0.005);
              setLevelByType(p => ({ ...p, [type]: newVal }));
              if (bgByType[type]?.meta) {
                setBgByType(prev => ({
                  ...prev,
                  [type]: {
                    ...prev[type],
                    meta: {
                      ...prev[type].meta,
                      level_start: newVal,
                      nonce: `${Date.now()}`
                    },
                  },
                }));
                bumpKey(type);
              }
            }}
            style={{ width: 70, fontSize: 10, padding: "1px 3px" }}
          />
        </label>
        <button
          onClick={() => applySettings(type)}
          style={{
            padding: "2px 6px", borderRadius: 4,
            border: "1px solid #cbd5e1",
            background: "#fff", cursor: "pointer", fontSize: 10,
          }}
          disabled={!pdataDirs?.[type] || !!busyByType?.[type]}
        >
          Appliquer
        </button>
      </div>
    );
  };

  /* =======================
     Rendu principal
     ======================= */
  return (
    <div className="viewer" style={{ marginTop: 8 }}>
      <h5 style={{ marginBottom: 8 }}>
        Vue combinée : ¹H · ¹³C · COSY · HSQC
      </h5>

      {loading && (
        <div className="status-item">Loading…</div>
      )}
      {err && (
        <div className="status-item error">Error: {err}</div>
      )}

      {!loading && !err && (
        <div
          style={{
            display: "grid",
            gridTemplateColumns: hideSidebar ? "1fr" : `${SIDEBAR_W}px 1fr`,
            gap: 12, width: "100%", maxWidth: "none",
            alignItems: "start" ,
          }}
        >
          {/* ==== Sidebar ==== */}
          {!hideSidebar && (
          <aside
            style={{
              background: "#f7f7f9", borderRadius: 8, padding: 8,
               display: "flex", flexDirection: "column",
               position: "sticky", top: 8,
               maxHeight: "calc(100vh - 16px)",
               overflow: "auto",
                alignSelf: "start" 

            }}
          >
            <input
              placeholder="Source…"
              value={query}
              onChange={(e) => setQuery(e.target.value)}
              className="parameter-field"
              style={{ marginBottom: 6, fontSize: 10,
                       padding: "4px 6px", flexShrink: 0 }}
            />
            <div
              style={{
                overflow: "auto", flex: 1,
                border: "1px solid #e5e7eb",
                borderRadius: 6, background: "#fff",
              }}
            >
              {allSources.length === 0 && (
                <div style={{ padding: 6, color: "#666", fontSize: 10 }}>
                  No match
                </div>
              )}
              {allSources.map((s) => {
                const stripMethodPrefix = (name) => String(name).replace(/^(?:1H|13C|COSY|HSQC)[\s_\-:]+/i, "");
                const label = stripMethodPrefix(s.replace(/_predicted$/i, ""));
                const active = selectedSource === s;
                
                return (
                  <button
                    key={s}
                    onClick={() => setSelectedSource(s)}
                    className={`type-tab ${active ? "active" : ""}`}
                    style={{
                      width: "100%", textAlign: "left",
                      padding: "6px 8px",
                      borderBottom: "1px solid #f0f0f0",
                      background: active ? "#e6f4ff" : "transparent",
                      cursor: "pointer", fontSize: 10,
                    }}
                    title={label}
                  >
                    <div
                      style={{
                        fontWeight: 600,
                        whiteSpace: "normal",
                        wordBreak: "break-word",
                        lineHeight: 1.3,
                      }}
                    >
                      {label}
                    </div>
                  </button>
                );
              })}
            </div>

            <div style={{ marginTop: 8, display: "grid", gap: 3, fontSize: 9 }}>
              <label style={{ display: "inline-flex", alignItems: "center", gap: 4 }}>
                <input
                  type="checkbox" checked={showTheory}
                  onChange={(e) => setShowTheory(e.target.checked)}
                  style={{ width: 12, height: 12 }}
                /> Théorie
              </label>
              <label style={{ display: "inline-flex", alignItems: "center", gap: 4 }}>
                <input
                  type="checkbox" checked={showExp}
                  onChange={(e) => setShowExp(e.target.checked)}
                  style={{ width: 12, height: 12 }}
                /> Exp.
              </label>
              <label style={{ display: "inline-flex", alignItems: "center", gap: 4 }}>
                <input
                  type="checkbox" checked={reverse1D}
                  onChange={(e) => setReverse1D(e.target.checked)}
                  style={{ width: 12, height: 12 }}
                /> 1D↓
              </label>
              <label style={{ display: "inline-flex", alignItems: "center", gap: 4 }}>
                <input
                  type="checkbox" checked={zoomRectEnabled}
                  onChange={(e) => setZoomRectEnabled(e.target.checked)}
                  style={{ width: 12, height: 12 }}
                /> Zoom
              </label>
              <label style={{ display: "inline-flex", alignItems: "center", gap: 4 }}>
                <input
                  type="checkbox" checked={show3D}
                  onChange={(e) => setShow3D(e.target.checked)}
                  style={{ width: 12, height: 12 }}
                /> 3D
              </label>
              <label style={{ display: "inline-flex", alignItems: "center", gap: 4 }}>
                <input
                  type="checkbox" checked={relevel2DEnabled}
                  onChange={(e) => setRelevel2DEnabled(e.target.checked)}
                  style={{ width: 12, height: 12 }}
                /> Relevel 2D
              </label>
              <div style={{ fontSize: 8, color: "#475569", marginTop: 2, lineHeight: 1.2 }}>
                {(canInteract1H || canInteract13C ||
                  canInteractCOSY || canInteractHSQC)
                  ? "✓ Interactions" : "✗ Pas d'interactions"}
              </div>
            </div>
          </aside>
          )}

          {/* ==== Contenu ==== */}
          <section style={{ minWidth: 0 }}>
            {/* Molécule + boutons */}
            <div style={{ marginBottom: 12 }}>
              <div
                style={{
                  display: "flex", justifyContent: "space-between",
                  alignItems: "center", marginBottom: 6,
                }}
              >
                <div style={{ fontWeight: 700, fontSize: 12 }}>
                  Molécule — colorations combinées
                </div>
                <div style={{ display: "flex", gap: 8 }}>
                  <button
                    onClick={() => setShowAllOpen(true)}
                    style={{
                      padding: "4px 12px", borderRadius: 6,
                      border: "1px solid #0ea5e9",
                      background: "#0ea5e9", color: "#fff",
                      cursor: hasAnyData ? "pointer" : "not-allowed",
                      opacity: hasAnyData ? 1 : 0.5,
                      fontSize: 11, fontWeight: 700,
                    }}
                    title="Show all molecules (2D)"
                    disabled={!hasAnyData}
                  >
                    Show all molecules
                  </button>
                  <button
                    onClick={() => setColorPickerOpen(true)}
                    style={{
                      padding: "4px 12px", borderRadius: 6,
                      border: "1px solid #2563eb",
                      background: "#2563eb", color: "#fff",
                      cursor: hasAnyData ? "pointer" : "not-allowed",
                      opacity: hasAnyData ? 1 : 0.5,
                      fontSize: 11, fontWeight: 600,
                    }}
                    disabled={!hasAnyData}
                  >
                    Change colors
                  </button>
                </div>
              </div>

              {selectedSource ? (
                <MoleculePane
                  csvUrl={pivotCsvUrl}
                  selectedRow={{ Source: selectedSource }}
                  show3D={show3D}
                  taskId={taskId}
                  customColors={customColors}
                  multi={{
                    fromAll: true,
                    pivotSource: selectedSource,
                    taskId,
                    rowsByType: {
                      "1H":  rowFor("1H"),
                      "13C": rowFor("13C"),
                      "COSY":rowFor("COSY"),
                      "HSQC":rowFor("HSQC"),
                    },
                  }}
                />
              ) : (
                <div
                  style={{
                    height: MOL_HEIGHT,
                    border: "1px dashed #e5e7eb",
                    borderRadius: 8,
                    display: "flex", alignItems: "center",
                    justifyContent: "center",
                    color: "#64748b", fontSize: 12,
                  }}
                >
                  {hasAnyData ? "No source selected" : "No molecule — run a prediction then an annotation"}
                </div>
              )}
            </div>

            {/* Spectres */}
            {(() => {
              const blocks = [];

              if (oneH) {
                const bg1H  = bgByType?.["1H"];
                const img   = bg1H?.pngUrl || null;
                const meta  = getMeta(bg1H);

                // Ajout : ppmPart pour refléter la fenêtre dans la key
                const ppmPart1H = meta
                  ? [meta.ppm_min, meta.ppm_max, meta.ppm_left, meta.ppm_right]
                      .map(v => v ?? "").join(",")
                  : "";

                blocks.push({
                  key: "b-1H",
                  el: (
                    <div
                      style={{
                        border: "1px solid #e5e7eb",
                        borderRadius: 8, padding: 8,
                        overflow: "hidden",
                      }}
                    >
                      <div
                        style={{
                          fontWeight: 700, marginBottom: 4,
                          display: "flex", justifyContent: "space-between",
                          alignItems: "center", fontSize: 12,
                        }}
                      >
                        <span>{TYPE_LABELS["1H"]}</span>
                        {img && (
                          <span
                            style={{
                              fontSize: 9, color: "#16a34a",
                              background: "#dcfce7",
                              padding: "1px 4px", borderRadius: 3,
                            }}
                          >
                            ✓
                          </span>
                        )}
                      </div>

                      {oneH.metrics &&
                        Object.keys(oneH.metrics).length > 0 && (
                        <div style={{ marginBottom: 6, minHeight: 20 }}>
                          {oneH.metrics.coverage_ransy && (
                            <Pill label="RANSY" value={oneH.metrics.coverage_ransy}/>
                          )}
                          {oneH.metrics.coverage_sub && (
                            <Pill label="Sub" value={oneH.metrics.coverage_sub}/>
                          )}
                          {oneH.metrics.avg_variabilite && (
                            <Pill label="Var" value={oneH.metrics.avg_variabilite}/>
                          )}
                        </div>
                      )}

                      {renderControls1D("1H")}

                      {/* Bouton View Detail pour le spectre 1H simulé */}
                      <div style={{ marginBottom: 8, display: "flex", gap: 8 }}>
                        <button
                          onClick={() => setSpectrumDetailOpen(true)}
                          disabled={!oneH?.smiles}
                          style={{
                            padding: "4px 10px",
                            fontSize: 11,
                            fontWeight: 500,
                            backgroundColor: oneH?.smiles ? "#3b82f6" : "#e5e7eb",
                            color: oneH?.smiles ? "#fff" : "#9ca3af",
                            border: "none",
                            borderRadius: 4,
                            cursor: oneH?.smiles ? "pointer" : "not-allowed",
                            display: "flex",
                            alignItems: "center",
                            gap: 4,
                          }}
                          title={oneH?.smiles ? "View simulated spectrum with peak shapes" : "No SMILES available"}
                        >
                          View Detail
                        </button>
                      </div>

                      <StickChart
                        key={[
                          refreshKeyByType["1H"] ?? 0,
                          img || "",
                          meta?.etag ?? meta?.version ?? meta?.last_modified ?? "",
                          JSON.stringify(meta?.plot_bbox_px ?? null) || "",
                          ppmPart1H,
                          String(meta?.nonce ?? "")
                        ].join("|")}
                        theory={showTheory ? oneH.theory : []}
                        expPpm={showExp ? oneH.expX : []}
                        expInt={showExp ? oneH.expI : []}
                        reverse={reverse1D}
                        showTheory={showTheory}
                        showExp={showExp}
                        height={PLOT_HEIGHT}
                        backgroundImage={img}
                        backgroundMeta={meta}
                        canInteract={canInteract1H}
                        zoomEnabled={zoomRectEnabled}
                        onZoom={(a, b) => doZoom1D("1H", a, b)}
                        onZoomPx={(r, d) => doZoom1DPixels("1H", r, d)}
                        onReset={() => doReset1D("1H")}
                        busy={!!busyByType?.["1H"]}
                      />
                    </div>
                  ),
                });
              }

              if (c13) {
                const bg13  = bgByType?.["13C"];
                const img   = bg13?.pngUrl || null;
                const meta  = getMeta(bg13);

                const ppmPart13 = meta
                  ? [meta.ppm_min, meta.ppm_max, meta.ppm_left, meta.ppm_right]
                      .map(v => v ?? "").join(",")
                  : "";

                blocks.push({
                  key: "b-13C",
                  el: (
                    <div
                      style={{
                        border: "1px solid #e5e7eb",
                        borderRadius: 8, padding: 8,
                        overflow: "hidden"
                      }}
                    >
                      <div
                        style={{
                          fontWeight: 700, marginBottom: 4,
                          display: "flex", justifyContent: "space-between",
                          alignItems: "center", fontSize: 12,
                        }}
                      >
                        <span>{TYPE_LABELS["13C"]}</span>
                        {img && (
                          <span
                            style={{
                              fontSize: 9, color: "#16a34a",
                              background: "#dcfce7",
                              padding: "1px 4px", borderRadius: 3,
                            }}
                          >
                            ✓
                          </span>
                        )}
                      </div>

                      {c13.metrics &&
                        Object.keys(c13.metrics).length > 0 && (
                        <div style={{ marginBottom: 6, minHeight: 20 }}>
                          {c13.metrics.coverage_ransy && (
                            <Pill label="RANSY" value={c13.metrics.coverage_ransy}/>
                          )}
                          {c13.metrics.coverage_sub && (
                            <Pill label="Sub" value={c13.metrics.coverage_sub}/>
                          )}
                          {c13.metrics.avg_variabilite && (
                            <Pill label="Var" value={c13.metrics.avg_variabilite}/>
                          )}
                        </div>
                      )}

                      {renderControls1D("13C")}

                      <StickChart
                        key={[
                          refreshKeyByType["13C"] ?? 0,
                          img || "",
                          meta?.etag ?? meta?.version ?? meta?.last_modified ?? "",
                          JSON.stringify(meta?.plot_bbox_px ?? null) || "",
                          ppmPart13,
                          String(meta?.nonce ?? "")
                        ].join("|")}
                        theory={showTheory ? c13.theory : []}
                        expPpm={showExp ? c13.expX : []}
                        expInt={showExp ? c13.expI : []}
                        reverse={reverse1D}
                        showTheory={showTheory}
                        showExp={showExp}
                        height={PLOT_HEIGHT}
                        backgroundImage={img}
                        backgroundMeta={meta}
                        canInteract={canInteract13C}
                        zoomEnabled={zoomRectEnabled}
                        onZoom={(a, b) => doZoom1D("13C", a, b)}
                        onZoomPx={(r, d) => doZoom1DPixels("13C", r, d)}
                        onReset={() => doReset1D("13C")}
                        busy={!!busyByType?.["13C"]}
                      />
                    </div>
                  ),
                });
              }

              if (cosy) {
                const bgC  = bgByType?.["COSY"];
                const img  = bgC?.pngUrl || null;
                const meta = getMeta(bgC);

                blocks.push({
                  key: "b-COSY",
                  el: (
                    <div
                      style={{
                        border: "1px solid #e5e7eb",
                        borderRadius: 8, padding: 8,
                        overflow: "hidden"
                      }}
                    >
                      <div
                        style={{
                          fontWeight: 700, marginBottom: 4,
                          display: "flex", justifyContent: "space-between",
                          alignItems: "center", fontSize: 12,
                        }}
                      >
                        <span>COSY</span>
                        {img && (
                          <span
                            style={{
                              fontSize: 9, color: "#16a34a",
                              background: "#dcfce7",
                              padding: "1px 4px", borderRadius: 3,
                            }}
                          >
                            ✓
                          </span>
                        )}
                      </div>

                      {cosy.metrics &&
                        Object.keys(cosy.metrics).length > 0 && (
                        <div style={{ marginBottom: 6, minHeight: 20 }}>
                          {cosy.metrics.coverage_ransy && (
                            <Pill label="RANSY" value={cosy.metrics.coverage_ransy}/>
                          )}
                          {cosy.metrics.coverage_sub && (
                            <Pill label="Sub" value={cosy.metrics.coverage_sub}/>
                          )}
                          {cosy.metrics.avg_variabilite && (
                            <Pill label="Var" value={cosy.metrics.avg_variabilite}/>
                          )}
                        </div>
                      )}

                      {renderControls2D("COSY")}

                      <Scatter2D
                        key={[
                          refreshKeyByType["COSY"] ?? 0,
                          img || "",
                          meta?.etag ?? meta?.version ?? meta?.last_modified ?? "",
                          JSON.stringify(meta?.plot_bbox_px ?? null) || "",
                          String(meta?.level_start ?? ""),
                          String(meta?.nonce ?? "")
                        ].join("|")}
                        theoryPts={showTheory ? cosy.theory : []}
                        expPts={showExp ? cosy.expPts : []}
                        expInt={showExp ? cosy.expI : []}
                        showTheory={showTheory}
                        showExp={showExp}
                        height={PLOT_HEIGHT}
                        backgroundImage={img}
                        backgroundMeta={meta}
                        canInteract={canInteractCOSY}
                        zoomEnabled={zoomRectEnabled}
                        onZoom2D={(w) => doZoom2D("COSY", w)}
                        onWheelDelta={(dy) => doWheel2D("COSY", dy)}
                        onReset={() => doReset2D("COSY")}
                        busy={!!busyByType?.["COSY"]}
                      />
                    </div>
                  ),
                });
              }

              if (hsqc) {
                const bgH  = bgByType?.["HSQC"];
                const img  = bgH?.pngUrl || null;
                const meta = getMeta(bgH);

                blocks.push({
                  key: "b-HSQC",
                  el: (
                    <div
                      style={{
                        border: "1px solid #e5e7eb",
                        borderRadius: 8, padding: 8,
                        overflow: "hidden"
                      }}
                    >
                      <div
                        style={{
                          fontWeight: 700, marginBottom: 4,
                          display: "flex", justifyContent: "space-between",
                          alignItems: "center", fontSize: 12,
                        }}
                      >
                        <span>HSQC</span>
                        {img && (
                          <span
                            style={{
                              fontSize: 9, color: "#16a34a",
                              background: "#dcfce7",
                              padding: "1px 4px", borderRadius: 3,
                            }}
                          >
                            ✓
                          </span>
                        )}
                      </div>

                      {hsqc.metrics &&
                        Object.keys(hsqc.metrics).length > 0 && (
                        <div style={{ marginBottom: 6, minHeight: 20 }}>
                          {hsqc.metrics.coverage_ransy && (
                            <Pill label="RANSY" value={hsqc.metrics.coverage_ransy}/>
                          )}
                          {hsqc.metrics.coverage_sub && (
                            <Pill label="Sub" value={hsqc.metrics.coverage_sub}/>
                          )}
                          {hsqc.metrics.avg_variabilite && (
                            <Pill label="Var" value={hsqc.metrics.avg_variabilite}/>
                          )}
                        </div>
                      )}

                      {renderControls2D("HSQC")}

                      <Scatter2D
                        key={[
                          refreshKeyByType["HSQC"] ?? 0,
                          img || "",
                          meta?.etag ?? meta?.version ?? meta?.last_modified ?? "",
                          JSON.stringify(meta?.plot_bbox_px ?? null) || "",
                          String(meta?.level_start ?? ""),
                          String(meta?.nonce ?? "")
                        ].join("|")}
                        theoryPts={showTheory ? hsqc.theory : []}
                        expPts={showExp ? hsqc.expPts : []}
                        expInt={showExp ? hsqc.expI : []}
                        showTheory={showTheory}
                        showExp={showExp}
                        height={PLOT_HEIGHT}
                        backgroundImage={img}
                        backgroundMeta={meta}
                        canInteract={canInteractHSQC}
                        zoomEnabled={zoomRectEnabled}
                        onZoom2D={(w) => doZoom2D("HSQC", w)}
                        onWheelDelta={(dy) => doWheel2D("HSQC", dy)}
                        onReset={() => doReset2D("HSQC")}
                        busy={!!busyByType?.["HSQC"]}
                      />
                    </div>
                  ),
                });
              }

              const n = blocks.length;
              if (n === 0) {
                if (!hasAnyData) {
                  const emptyBox = (title) => (
                    <div
                      style={{
                        border: "1px solid #e5e7eb",
                        borderRadius: 8,
                        padding: 8,
                        overflow: "hidden",
                        background: "#fff",
                      }}
                    >
                      <div
                        style={{
                          fontWeight: 700,
                          marginBottom: 6,
                          display: "flex",
                          justifyContent: "space-between",
                          alignItems: "center",
                          fontSize: 12,
                        }}
                      >
                        <span>{title}</span>
                      </div>
                      <div
                        style={{
                          height: PLOT_HEIGHT,
                          border: "1px dashed #e5e7eb",
                          borderRadius: 8,
                          display: "flex",
                          alignItems: "center",
                          justifyContent: "center",
                          color: "#64748b",
                          fontSize: 12,
                          background: "#fff",
                        }}
                      >
                        No spectrum — run the annotation
                      </div>
                    </div>
                  );

                  return (
                    <div
                      style={{
                        display: "grid",
                        gridTemplateColumns: "1fr 1fr",
                        gridTemplateRows: "1fr 1fr",
                        gap: 12,
                        marginBottom: 12,
                        alignItems: "start",
                        gridAutoRows: "min-content",
                      }}
                    >
                      <div>{emptyBox(TYPE_LABELS["1H"])}</div>
                      <div>{emptyBox(TYPE_LABELS["13C"])}</div>
                      <div>{emptyBox("COSY")}</div>
                      <div>{emptyBox("HSQC")}</div>
                    </div>
                  );
                }

                return <div>No spectrum available</div>;
              }
              if (n === 1) return <div style={{ marginBottom: 12 }}>{blocks[0].el}</div>;
              if (n === 2)
                return (
                  <div
                    style={{
                      display: "grid",
                      gridTemplateColumns: "1fr 1fr",
                      gap: 12, marginBottom: 12,
                      alignItems: "start",
                      gridAutoRows: "min-content"
                    }}
                  >
                    {blocks.map(b => <div key={b.key}>{b.el}</div>)}
                  </div>
                );
              if (n === 3)
                return (
                  <div style={{ marginBottom: 12 }}>
                    <div
                      style={{
                        display: "grid",
                        gridTemplateColumns: "1fr 1fr",
                        gap: 12, marginBottom: 12,
                        alignItems: "start",
                        gridAutoRows: "min-content"
                      }}
                    >
                      {blocks.slice(0, 2)
                        .map(b => <div key={b.key}>{b.el}</div>)}
                    </div>
                    <div style={{ display: "flex", justifyContent: "center" }}>
                      <div style={{ width: "48%" }}>{blocks[2].el}</div>
                    </div>
                  </div>
                );
              return (
                <div
                  style={{
                    display: "grid",
                    gridTemplateColumns: "1fr 1fr",
                    gridTemplateRows: "1fr 1fr",
                    gap: 12, marginBottom: 12,
                    alignItems: "start",
                    gridAutoRows: "min-content"
                  }}
                >
                  {blocks.map(b => <div key={b.key}>{b.el}</div>)}
                </div>
              );
            })()}
          </section>
        </div>
      )}

      {/* Sélecteur de couleurs (ouvert via bouton) */}
      <ColorPicker
        isOpen={colorPickerOpen}
        onClose={() => setColorPickerOpen(false)}
        customColors={customColors}
        onColorChange={setCustomColor}
      />

      {/* Modale — toutes les molécules */}
      <ShowAllMoleculesModal
        isOpen={showAllOpen}
        onClose={() => setShowAllOpen(false)}
        sources={sortedSourcesAll}
        page={allPage}
        setPage={setAllPage}
        pageSize={allPageSize}
        setPageSize={setAllPageSize}
        getRowsByType={rowsByTypeForSource}
        csvUrl={pivotCsvUrl}
        taskId={taskId}
        customColors={customColors}
      />

      {/* Modale — spectre 1H simulé détaillé */}
      <SpectrumDetailModal
        isOpen={spectrumDetailOpen}
        onClose={() => setSpectrumDetailOpen(false)}
        smiles={oneH?.smiles || ""}
        molId={oneH?.source || selectedSource || "mol"}
        theoryShifts={oneH?.theory || []}
        expShifts={oneH?.expX || []}
      />
    </div>
  );
}

if (typeof window !== 'undefined' && !window.__NMR_ERROR_FILTER_INSTALLED__) {
  window.__NMR_ERROR_FILTER_INSTALLED__ = true;
  
  const originalError = console.error;
  console.error = function(...args) {
    const msg = args.join(' ');
    
    // Filter known non-critical WebGL errors
    if (msg.includes('GL_INVALID_FRAMEBUFFER_OPERATION') ||
        msg.includes('Framebuffer is incomplete') ||
        msg.includes('Attachment has zero size')) {
      if (window.DEBUG_NMR) {
        originalError.apply(console, ['[WebGL Filtré]', ...args]);
      }
      return;
    }
    
    // Display real errors
    originalError.apply(console, args);
  };
  
  debugLog('[OK] WebGL error filter installed');
}


