import { isMethodAvailable, availableMethodEntries } from './constants/releaseFeatures';
// NMR results viewer with molecule rendering and spectrum controls.








import React, {
  useEffect,
  useMemo,
  useState,
  useCallback,
  useRef,
} from "react";
import { jsPDF } from "jspdf";
import MoleculeScoreSummary from "./components/viewer/MoleculeScoreSummary";

import {
  interactive1H,
  interactive13C,
  interactiveCOSY,
  interactiveHSQC,
  API_BASE_URL,
  extractBrukerSpectrum,
  extractBruker2D,
  getMultiplicityData,
  simulateSpectrum,
  simulateSpectrumBatch,
} from "./api";





// Color constants and styles
import {
  DEFAULT_COLORS,
  METHOD_LABELS_SHORT,
  TYPE_LABELS,
  PLOT_HEIGHT,
  MOL_HEIGHT,
  SIDEBAR_W,
} from "./constants/colors";

// Utilities
import { debugLog, debugWarn } from "./utils/debug";
import { fetchWithRetry } from "./utils/imageUtils";
import { parseCSVSmart } from "./utils/csvParser";
import { homonuclearSpectrumRows } from "./utils/homonuclearSpectrumRows";
import { findBondsInPath } from "./utils/molParser";
import { tocsyHighlights } from "./utils/tocsyHighlights";
import { useSpectrumWindowSettings } from "./hooks/useSpectrumWindowSettings";
import { useMeasureWidth } from "./hooks/useMeasureWidth";
import { 
  loadRDKitFromCDN, 
  load3DmolFromCDN, 
  isWebGLAvailable 
} from "./utils/cdnLoaders";


// Viewer components
import { Pill, ColorPicker, StickChart, Scatter2D, SpectrumDetailModal } from "./components/viewer";




// Enable debug by typing in console: window.DEBUG_NMR = true
if (typeof window !== 'undefined' && !window.hasOwnProperty('DEBUG_NMR')) {
  window.DEBUG_NMR = false; // false by default in production
}

/*   CDN Loaders (RDKit, 3Dmol) - imported from utils/cdnLoaders.js
   Helpers isProblematicBrowser, isWebGLAvailable also imported*/

/*   CSV helpers*/



/*   SDF / MOL helpers*/

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
  
  // Match the SDF property header emitted by the backend: >  <(1)>  (1).
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

  // Values are already space-separated
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
    // Use API_BASE_URL as the base URL in Electron.
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
      // Use API_BASE_URL for requests from Electron.
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
    // Use API_BASE_URL for requests from Electron.
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

/*  */

/*  */

/*  */

function MoleculePane({
  csvUrl,
  selectedRow,
  show3D,
  taskId,
  multi = null,
  customColors,
  molZoom = 1,
  molHeight = MOL_HEIGHT,
  flagsOverride = null,            // Optional shared method toggles.
  showLocalToggles = true,         // Show or hide the local method toggles.
  compact = false,                 // Thumbnail mode renders only the SVG.
}) {
  const [svg, setSvg] = useState(null);
  const [drawingRef, drawingWidth] = useMeasureWidth(compact ? 180 : 600);
  const [err, setErr] = useState(null);
  const [loading, setLoading] = useState(false);

  // Local checkboxes (used only if flagsOverride is not provided)
  const [localFlags, setLocalFlags] = useState({
    H1: true, C13: true, COSY: true, TOCSY: isMethodAvailable('TOCSY'), HSQC: true,
  });
  const flags = flagsOverride || localFlags;       // single source of truth
  const toggle = (k) => setLocalFlags(f => ({ ...f, [k]: !f[k] }));

  const div3DRef = useRef(null);
  const viewerRef = useRef(null);

  // Safe accessors — customColors may come from sessionStorage with missing fields
  const defaultRgb = { H1: [0, 0.72, 0.66], C13: [0.23, 0.05, 0.64], COSY: [1, 0.69, 0], TOCSY: [0, 0.66, 0.47], HSQC: [0.97, 0.15, 0.52] };
  const defaultHex = { H1: '#00B8A9', C13: '#3A0CA3', COSY: '#FFB000', TOCSY: '#00A878', HSQC: '#F72585' };
  const safeRgb = (key) => customColors?.[key]?.rgb || defaultRgb[key] || [0, 0.5, 0.5];
  const safeHex = (key) => customColors?.[key]?.hex || defaultHex[key] || '#008080';

  // eslint-disable-next-line react-hooks/exhaustive-deps
  const COLORS_RDK = {
    H1:  safeRgb('H1'),
    C13: safeRgb('C13'),
    COSY:safeRgb('COSY'),
    TOCSY:safeRgb('TOCSY'),
    HSQC:safeRgb('HSQC'),
  };
  // eslint-disable-next-line react-hooks/exhaustive-deps
  const COLORS_3D = {
    H1:  safeHex('H1'),
    C13: safeHex('C13'),
    COSY:safeHex('COSY'),
    TOCSY:safeHex('TOCSY'),
    HSQC:safeHex('HSQC'),
  };

  // HSQC C–H bond highlighting in 2D.
  const computeMultiHighlights = useCallback((resolved, rowsByType) => {
    const atoms  = parseAtomSymbolsFromMolblock(resolved.molblock) || [];
    const shifts = extractShiftArrayFromSDFRecord(resolved.rawText) || [];
    
    // Log diagnostics only in debug mode.
    debugLog("[DEBUG] Shifts extraits:", shifts);
    debugLog("[DEBUG] Symboles d'atomes:", atoms);
    debugLog("[DEBUG] Nombre d'atomes:", atoms.length);
    
    // Handle empty arrays.
    if (!atoms.length || !shifts.length) {
      debugWarn("[WARN] No atoms or shifts found, coloring disabled");
      return { 
        highlightByType: { H1: new Set(), C13: new Set(), COSY: new Set(), TOCSY: new Set(), HSQC: new Set() },
        bondsCOSY: [],
        bondsTOCSY: [],
        tocsyStickAtoms: new Set(),
        tocsyCrossCarbons: new Set(),
        cosyDiagonalProtons: new Set(),
        cosyDiagonalCarbons: new Set(),
        cosyCrosspeakCarbons: new Set(),
        bondsHSQC: [] 
      };
    }
    
    const maxIdx = Math.min(shifts.length - 1, atoms.length - 1);

    const hIdx   = atoms.map((s,i)=>s==="H"?i:-1).filter(i=>i>=0);
    const cIdx   = atoms.map((s,i)=>s==="C"?i:-1).filter(i=>i>=0);

    const sets   = { H1:new Set(), C13:new Set(), COSY:new Set(), TOCSY:new Set(), HSQC:new Set() };
    const cosyDiagonalProtons = new Set();  // Protons matched to COSY diagonal peaks.
    const cosyDiagonalCarbons = new Set();  // Parent carbons of COSY diagonal protons.
    const cosyCrosspeakCarbons = new Set(); // Parent carbons of COSY cross-peak protons.
    let bondsCOSY = [];
    let bondsCOSYDiagonal = [];  // C–H bonds matched to COSY diagonal peaks.
    let bondsHSQC = [];

    // Helper: round to 3 decimals for exact matching
    const r3 = v => Math.round(v * 1000) / 1000;
    
    // Find a proton's parent carbon in the molblock.
    const findCarbonParent = (hIdx) => {
      const lines = resolved.molblock.split(/\r?\n/);
      const counts = lines[3] || "";
      const nAtoms = parseInt(counts.slice(0,3).trim(),10);
      const nBonds = parseInt(counts.slice(3,6).trim(),10);
      
      for (let bi=0; bi<nBonds; bi++) {
        const l = lines[4+nAtoms+bi] || "";
        const a1 = parseInt(l.slice(0,3).trim(),10)-1;
        const a2 = parseInt(l.slice(3,6).trim(),10)-1;
        if (a1 === hIdx || a2 === hIdx) {
          const other = a1 === hIdx ? a2 : a1;
          // Check whether the other atom is carbon.
          const atomLine = lines[4+other] || "";
          const symbol = atomLine.slice(31,34).trim();
          if (symbol === "C") {
            return other;
          }
        }
      }
      return -1;
    };

    // 1H - Try matched_atom_indices first
    if (flags.H1 && rowsByType["1H"]) {
      const atomIndicesStr = rowsByType["1H"].matched_atom_indices || "";
      if (atomIndicesStr) {
        // Format: "idx1|idx2|idx3|..." where indices are separated by pipes
        const parts = atomIndicesStr.split("|").filter(Boolean);
        for (const part of parts) {
          const indices = [parseInt(part, 10)].filter(Number.isFinite);
          for (const idx of indices) if (idx >= 0) sets.H1.add(idx);
        }
      } else {
        // Fallback: match by shift
        const m = (rowsByType["1H"].ppm_molecule_matched||"")
          .split("|").map(v=>parseFloat(String(v).replace(",","."))).filter(Number.isFinite);
        const r = new Set(m.map(r3));
        for (const ppm of r) for (const i of hIdx)
          if (i<=maxIdx && r3(shifts[i])===ppm) sets.H1.add(i);
      }
    }
    
    // 13C - Try matched_atom_indices first
    if (flags.C13 && rowsByType["13C"]) {
      const atomIndicesStr = rowsByType["13C"].matched_atom_indices || "";
      if (atomIndicesStr) {
        // Format: "idx1~idx2|idx3~idx4|..." where ~ separates equivalent carbons, | separates groups
        const parts = atomIndicesStr.split("|").filter(Boolean);
        for (const part of parts) {
          // Parse tilde-separated indices for equivalent carbons
          const indices = part.split("~").map(s => parseInt(s.trim(), 10)).filter(Number.isFinite);
          for (const idx of indices) if (idx >= 0) sets.C13.add(idx);
        }
      } else {
        // Fallback: match by shift
        const m = (rowsByType["13C"].ppm_molecule_matched||"")
          .split("|").map(v=>parseFloat(String(v).replace(",","."))).filter(Number.isFinite);
        const r = new Set(m.map(r3));
        for (const ppm of r) for (const i of cIdx)
          if (i<=maxIdx && r3(shifts[i])===ppm) sets.C13.add(i);
      }
    }
    
    // COSY - Try matched_atom_indices first
    if (flags.COSY && rowsByType["COSY"]) {
      const atomIndicesStr = rowsByType["COSY"].matched_atom_indices || "";
      const lines = resolved.molblock.split(/\r?\n/);
      const nAtoms = parseInt((lines[3]||"").slice(0,3).trim(),10);

      if (atomIndicesStr) {
        // Format: "H1_all:H2_all||..." where H1_all/H2_all can be "34" or "34~35~36" for equivalent protons
        const entries = atomIndicesStr.split("||").filter(Boolean);
        const matchedNormal = [];
        const matchedDiagonal = [];
        for (const entry of entries) {
          const [h1Str, h2Str] = entry.split(":");
          if (!h1Str || !h2Str) continue;
          
          // Parse equivalent proton indices separated by ~.
          const h1List = h1Str.split("~").map(s => parseInt(s.trim(), 10)).filter(Number.isFinite);
          const h2List = h2Str.split("~").map(s => parseInt(s.trim(), 10)).filter(Number.isFinite);
          
          if (h1List.length === 0 || h2List.length === 0) continue;
          
          // Check if diagonal (all H1 same as all H2, or same group)
          const isDiagonal = (h1List.length === h2List.length && 
                            h1List.every((h, i) => h === h2List[i])) ||
                            (h1List.length === 1 && h2List.length === 1 && h1List[0] === h2List[0]);
          
          if (isDiagonal) {
            // Highlight every proton in the matched diagonal group.
            
            for (const h of h1List) {
              if (h >= 0) {
                cosyDiagonalProtons.add(h);
                const carbonIdx = findCarbonParent(h);
                if (carbonIdx >= 0) {
                  cosyDiagonalCarbons.add(carbonIdx);
                }
              }
            }
            // Use first pair for path finding
            matchedDiagonal.push([h1List[0], h2List[0]]);
          } else {
            // Normal correlation (cross-peak): add all protons to COSY set
            // Also add their parent carbons for full stick coloring
            for (const h of h1List) {
              if (h >= 0) {
                sets.COSY.add(h);
                // Add parent carbon for full C-H stick coloring
                const carbonIdx = findCarbonParent(h);
                if (carbonIdx >= 0) {
                  cosyCrosspeakCarbons.add(carbonIdx);
                }
              }
            }
            for (const h of h2List) {
              if (h >= 0) {
                sets.COSY.add(h);
                // Add parent carbon for full C-H stick coloring
                const carbonIdx = findCarbonParent(h);
                if (carbonIdx >= 0) {
                  cosyCrosspeakCarbons.add(carbonIdx);
                }
              }
            }
            // Use the first proton pair to find the bond path.
            matchedNormal.push([h1List[0], h2List[0]]);
          }
        }
        if (matchedNormal.length) {
          bondsCOSY = findBondsInPath(resolved.molblock, matchedNormal);
        }
        if (matchedDiagonal.length) {
          // Resolve C–H bonds separately for diagonal peaks.
          bondsCOSYDiagonal = findBondsInPath(resolved.molblock, matchedDiagonal);
        }
      } else {
        // Fallback: match by shift
        const pairs = (rowsByType["COSY"].ppm_molecule_matched||"")
          .split("||").map(p=>p.split("|")).filter(a=>a.length===2)
          .map(([f1,f2])=>[+String(f1).replace(",","."),+String(f2).replace(",","." )])
          .filter(([a,b])=>Number.isFinite(a)&&Number.isFinite(b));
        const matched = [];
        for (const [d1,d2] of pairs) {
          const r1=r3(d1), r2=r3(d2);
          const a1=[], a2=[];
          for (const i of hIdx) if (i<=maxIdx && r3(shifts[i])===r1) a1.push(i);
          for (const i of hIdx) if (i<=maxIdx && r3(shifts[i])===r2) a2.push(i);
          for (const x of a1) for (const y of a2) {
            matched.push([x,y]);
            // Include protons in the COSY highlight set.
            sets.COSY.add(x);
            sets.COSY.add(y);
          }
        }
        if (matched.length) {
          bondsCOSY = findBondsInPath(resolved.molblock, matched);
        }
      }
    }

    const tocsy = tocsyHighlights(
      resolved.molblock, shifts, flags.TOCSY ? rowsByType.TOCSY : undefined
    );
    sets.TOCSY = tocsy.protons;

    // HSQC direct C–H correlations.
    if (flags.HSQC && rowsByType["HSQC"]) {
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

      // Try matched_atom_indices first (format: "C_idx~C_idx:H_idx1|H_idx2||..." for equivalent atoms)
      const atomIndicesStr = rowsByType["HSQC"].matched_atom_indices || "";
      if (atomIndicesStr) {
        const entries = atomIndicesStr.split("||").filter(Boolean);
        for (const entry of entries) {
          const [cPart, hPart] = entry.split(":");
          // Support multiple equivalent carbons separated by ~ 
          const cIndices = (cPart || "").split("~").map(s => parseInt(s, 10)).filter(Number.isFinite);
          const hIndices = (hPart || "").split("|").map(s => parseInt(s, 10)).filter(Number.isFinite);
          // Add all C-H pairs for equivalent atoms
          for (const ci of cIndices) {
            if (ci < 0) continue;
            // Don't add carbon to HSQC set - only color proton stick, not carbon sphere/stick
            for (const hi of hIndices) {
              if (hi >= 0) {
                sets.HSQC.add(hi);  // Only color the proton
                const key = `${Math.min(ci,hi)}-${Math.max(ci,hi)}`;
                const b = bondMap.get(key);
                if (b != null) bondsHSQC.push(b);
              }
            }
          }
        }
      } else {
        // Fallback: match by shift values
        const pairs = (rowsByType["HSQC"].ppm_molecule_matched_unfolded ||
                       rowsByType["HSQC"].ppm_molecule_matched || "")
          .split("||").map(p=>p.split("|")).filter(a=>a.length===2)
          .map(([c,h])=>[+String(c).replace(",","."), +String(h).replace(",","." )])
          .filter(([c,h])=>Number.isFinite(c)&&Number.isFinite(h));

        for (const [cppm,hppm] of pairs) {
          const rc = r3(cppm);
          const rh = r3(hppm);

          const cCandidates = [];
          const hCandidates = [];

          for (const ci of cIdx)
            if (ci<=maxIdx && r3(shifts[ci])===rc) cCandidates.push(ci);

          for (const hi of hIdx)
            if (hi<=maxIdx && r3(shifts[hi])===rh) {
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
    }
    
    return { highlightByType: sets, bondsCOSY, bondsCOSYDiagonal, bondsHSQC, cosyDiagonalProtons, cosyDiagonalCarbons, cosyCrosspeakCarbons,
      bondsTOCSY: tocsy.bonds, tocsyStickAtoms: tocsy.stickAtoms, tocsyCrossCarbons: tocsy.crossCarbons };
  }, [flags]);
  

  const render2D = useCallback(async (resolved) => {
    setLoading(true);
    try {
      const RDKit = await loadRDKitFromCDN();

      // Keep explicit hydrogens: backend matched_atom_indices refer to the original molblock.
      const mol = RDKit.get_mol(resolved.molblock);

      let highlightAtoms=[], atomColorsMap={}, highlightBonds=[], bondColorsMap={};

      if (multi?.fromAll && multi?.rowsByType) {
        const { highlightByType, bondsCOSY, bondsTOCSY, bondsHSQC } =
          computeMultiHighlights(resolved, multi.rowsByType);

        const cosySet = highlightByType.COSY;
        const tocsySet = highlightByType.TOCSY;
        const hsqcSet = highlightByType.HSQC;

        // Indices already refer to the original molblock.
        if (flags.C13)
          for (const i of highlightByType.C13)
            atomColorsMap[String(i)] = COLORS_RDK.C13;
        if (flags.H1)
          for (const i of highlightByType.H1)
            atomColorsMap[String(i)] = COLORS_RDK.H1;
        if (flags.HSQC)
          for (const i of hsqcSet)
            atomColorsMap[String(i)] = COLORS_RDK.HSQC;
        if (flags.COSY)
          for (const i of cosySet)
            atomColorsMap[String(i)] = COLORS_RDK.COSY;
        if (flags.TOCSY)
          for (const i of tocsySet)
            atomColorsMap[String(i)] = COLORS_RDK.TOCSY;

        highlightAtoms = Array.from(new Set(Object.keys(atomColorsMap).map(Number)));

        // Combine matched atoms using the highlight color priority.
        const atomPriorityColor = {};
        const priority = [COLORS_RDK.COSY, COLORS_RDK.TOCSY, COLORS_RDK.HSQC, COLORS_RDK.H1, COLORS_RDK.C13];
        for (const [idxStr, color] of Object.entries(atomColorsMap)) {
          const idx = Number(idxStr);
          if (!atomPriorityColor[idx] ||
              priority.indexOf(color) < priority.indexOf(atomPriorityColor[idx])) {
            atomPriorityColor[idx] = color;
          }
        }
        const allMatched = new Set(Object.keys(atomPriorityColor).map(Number));

        // Highlight direct bonds between matched atoms in the original molblock.
        const origLines = resolved.molblock.split(/\r?\n/);
        const nAtO = parseInt((origLines[3]||'').slice(0,3).trim(), 10) || 0;
        const nBdO = parseInt((origLines[3]||'').slice(3,6).trim(), 10) || 0;
        for (let bi = 0; bi < nBdO; bi++) {
          const l = origLines[4 + nAtO + bi] || '';
          const a1 = parseInt(l.slice(0,3).trim(), 10) - 1;
          const a2 = parseInt(l.slice(3,6).trim(), 10) - 1;
          if (allMatched.has(a1) && allMatched.has(a2)) {
            highlightBonds.push(bi);
            const c1 = atomPriorityColor[a1];
            const c2 = atomPriorityColor[a2];
            const chosen = priority.find(p => p === c1 || p === c2) || c1 || c2;
            bondColorsMap[String(bi)] = chosen;
          }
        }

        if (flags.HSQC)
          for (const b of bondsHSQC) { highlightBonds.push(b); bondColorsMap[String(b)] = COLORS_RDK.HSQC; }
        if (flags.TOCSY)
          for (const b of bondsTOCSY) { highlightBonds.push(b); bondColorsMap[String(b)] = COLORS_RDK.TOCSY; }
        if (flags.COSY)
          for (const b of bondsCOSY) { highlightBonds.push(b); bondColorsMap[String(b)] = COLORS_RDK.COSY; }

        highlightBonds = Array.from(new Set(highlightBonds));
      }

      const renderScale = 3 * Math.max(0.6, molZoom);
      const drawingHeight = compact ? molHeight : Math.max(1, molHeight - 42);
      const pad = compact ? 0.035 : Math.max(0.02, 0.04 / Math.max(1, molZoom));

      const opts = {
        width: Math.round(drawingWidth * renderScale),
        height: Math.round(drawingHeight * renderScale),
        useMolBlockCoords: false,
        padding: pad, bondLineWidth: 1.5, atomLabelFontSize: 11,
        addStereoAnnotation: false,
        explicitMethyl: false,
        annotationFontScale: 0.9,
        backgroundColour: [1, 1, 1, 1],
      };

      let svgStr;
      try {
        svgStr = mol.get_svg_with_highlights(JSON.stringify({
          ...opts,
          atoms: highlightAtoms,
          highlight_atom_colors: atomColorsMap,
          bonds: highlightBonds,
          highlight_bond_colors: bondColorsMap,
          highlightRadius: 0.2,
          atomHighlightsAreCircles: true,
          highlightBondWidthMultiplier: 8,
        }));
      } catch {
        svgStr = mol.get_svg(JSON.stringify(opts));
      }

      svgStr = svgStr
        .replace(/\bwidth=(["'])[^"']+\1/, 'width="100%"')
        .replace(/\bheight=(["'])[^"']+\1/, 'height="100%"')
        .replace('<svg ', '<svg preserveAspectRatio="xMidYMid meet" ');

      setSvg(svgStr);
      setErr(null);
      mol.delete?.();
    } catch (e) {
      setErr("2D Error: " + (e?.message || String(e)));
    } finally {
      setLoading(false);
    }
  }, [multi, computeMultiHighlights, COLORS_RDK, molZoom, flags, compact, drawingWidth, molHeight]);

  const render3D = useCallback(async (resolved) => {
    setLoading(true);
    try {
      // Check WebGL availability before creating the viewer.
      if (!isWebGLAvailable()) {
        throw new Error("WebGL is not available on this browser.\n\n" +
          "Solutions:\n" +
          "1. Firefox: Type about:config in the address bar, then set webgl.disabled to false\n" +
          "2. Update your graphics drivers\n" +
          "3. Try Chrome/Edge if the problem persists");
      }

      const $3Dmol = await load3DmolFromCDN();
      const el = div3DRef.current;
      
      // Wait for the viewer container.
      if (!el) {
        throw new Error("3D container not found");
      }
      
      // Wait until the container has a nonzero size.
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
      
      // Check that the viewer was created.
      if (!viewer) {
        throw new Error("Failed to create 3D viewer");
      }
      
      viewer.removeAllModels();
      const sdf = resolved.firstRecord || resolved.rawText || resolved.molblock;
      const model = viewer.addModel(sdf, "sdf");
      
      // Stop if the molecule model could not be created.
      if (!model) {
        throw new Error("Failed to load 3D model");
      }
      
      model.setStyle({}, {
        stick:{radius:0.15,color:"#808080"},
        sphere:{scale:0.2,color:"#808080"}
      });

      if (multi?.fromAll && multi?.rowsByType) {
        const { highlightByType, tocsyStickAtoms, tocsyCrossCarbons, cosyDiagonalCarbons: diagCarbons, cosyDiagonalProtons: diagProtons, cosyCrosspeakCarbons: crossCarbons } = computeMultiHighlights(resolved, multi.rowsByType);
        const CAtoms = flags.C13 ? Array.from(highlightByType.C13) : [];
        const cosySet = flags.COSY ? new Set(highlightByType.COSY) : new Set();
        const tocsySet = flags.TOCSY ? new Set(highlightByType.TOCSY) : new Set();
        const hsqcSet = flags.HSQC ? new Set(highlightByType.HSQC) : new Set();
        const h1List  = flags.H1   ? Array.from(highlightByType.H1)  : [];
        
        // Parent carbons of COSY diagonal protons.
        const cosyDiagonalCarbons = flags.COSY ? diagCarbons : new Set();
        // COSY diagonal protons whose bonds use the COSY color.
        const cosyDiagonalProtons = flags.COSY ? diagProtons : new Set();
        // Parent carbons of COSY cross-peak protons. (for C–H bond highlighting)
        const cosyCrosspeakCarbons = flags.COSY ? crossCarbons : new Set();

        // Highlight priority
        // Atom spheres: H1 > HSQC > C13 > COSY > TOCSY.
        // Stick: preserve the experiment-specific color when available.
        

        // Use sets for highlight membership checks.
        const c13Set = new Set(CAtoms);
        const h1Set = new Set(h1List);

        // Collect atoms that need highlighting.
        const allAtoms = new Set([
          ...CAtoms,
          ...cosySet,
          ...tocsySet,
          ...tocsyStickAtoms,
          ...hsqcSet,
          ...h1List,
          ...cosyDiagonalCarbons,
          ...cosyDiagonalProtons,
          ...cosyCrosspeakCarbons  // Parent carbons of cross-peak protons.
        ]);

        // Apply highlight styles in priority order.
        for (const i of allAtoms) {
          if (typeof i !== 'number' || i < 0) continue;
          
          const inCOSY = cosySet.has(i);  // Protons in off-diagonal COSY peaks.
          const inTOCSY = tocsySet.has(i) || tocsyStickAtoms.has(i);
          const inCOSYCrossC = cosyCrosspeakCarbons.has(i);  // Parent carbon of a cross-peak proton.
          const inCOSYDiagC = cosyDiagonalCarbons.has(i);  // Parent carbon of a diagonal proton.
          const inCOSYDiagH = cosyDiagonalProtons.has(i);  // Diagonal proton.
          const inC13 = c13Set.has(i);
          // hsqcSet already excludes atoms when the HSQC toggle is off.
          const inHSQC = hsqcSet.has(i);
          const inH1 = h1Set.has(i);
          const inAnyCOSY = inCOSY || inCOSYCrossC || inCOSYDiagC || inCOSYDiagH;
          
          // Check whether the proton belongs to a COSY cross-peak or diagonal peak.
          const isHydrogenCOSY = inCOSYDiagH || (inCOSY && !inC13 && !inCOSYDiagC);
          // Draw spheres for carbons linked to COSY cross-peaks.
          // Cross-peak highlighting takes priority when a carbon also belongs to a diagonal match.
          const isCarbonCOSYCross = inCOSYCrossC;  
          // Diagonal-only carbons have no sphere.
          const isCarbonCOSYDiagOnly = inCOSYDiagC && !inCOSYCrossC;

          // Choose the sphere color: H1 > HSQC > C13 > COSY cross-peak carbon.
          // For diagonal-only COSY matches, highlight the C–H bond without a carbon sphere.
          let sphereColor = null, sphereRadius = 0;
          if (inH1) {
            sphereColor = COLORS_3D.H1;
            sphereRadius = 0.35;
          } else if (inHSQC) {
            // HSQC takes priority over COSY for atom spheres.
            sphereColor = COLORS_3D.HSQC;
            sphereRadius = 0.35;
          } else if (inC13) {
            sphereColor = COLORS_3D.C13;
            sphereRadius = 0.4;
          } else if (isCarbonCOSYCross) {
            // COSY cross-peaks take priority over diagonal matches for atom spheres.
            sphereColor = COLORS_3D.COSY;
            sphereRadius = 0.4;
          } else if (isCarbonCOSYDiagOnly) {
            // Diagonal-only COSY carbons use a C–H bond highlight without a sphere.
            sphereColor = null;
            sphereRadius = 0;
          } else if (tocsyCrossCarbons.has(i)) {
            sphereColor = COLORS_3D.TOCSY;
            sphereRadius = 0.4;
          } else if (inTOCSY) {
            // Like COSY: protons and diagonal carbons use sticks only.
            sphereColor = null;
            sphereRadius = 0;
          } else if (inAnyCOSY && isHydrogenCOSY) {
            // COSY protons use bond highlights without spheres; continue to apply the bond style.
            sphereColor = null;
            sphereRadius = 0;
          } else {
            continue; // Skip atoms outside the highlight sets.
          }

          // Choose the bond color: HSQC, then COSY, then the default gray.
          // Apply the COSY bond color to both protons and their parent carbons.
          // Styling both endpoints colors the full C–H bond.
          let stickColor, stickRadius;
          if (inHSQC) {
            // HSQC also takes priority over COSY for bond colors.
            stickColor = COLORS_3D.HSQC;
            stickRadius = 0.2;
          } else if (inCOSY || inCOSYCrossC) {
            // For COSY cross-peaks, color the bonds of both protons and their parent carbons.
            stickColor = COLORS_3D.COSY;
            stickRadius = 0.25;
          } else if (inCOSYDiagH || inCOSYDiagC) {
            // For COSY diagonal matches, color the bonds of both the proton and its parent carbon.
            stickColor = COLORS_3D.COSY;
            stickRadius = 0.2;
          } else if (inTOCSY) {
            stickColor = COLORS_3D.TOCSY;
            stickRadius = tocsyCrossCarbons.has(i) || tocsySet.has(i) ? 0.25 : 0.2;
          } else {
            stickColor = "#808080";
            stickRadius = 0.15;
          }

          try {
            // Apply only the bond style when no atom sphere is needed.
            if (sphereColor === null) {
              model.setStyle({ index:i }, { 
                stick:{radius:stickRadius,color:stickColor} 
              });
            } else {
              model.setStyle({ index:i }, { 
                sphere:{radius:sphereRadius,color:sphereColor}, 
                stick:{radius:stickRadius,color:stickColor} 
              });
            }
          } catch (e) {
            debugWarn(`[WARN] setStyle error for atom ${i}:`, e.message);
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
            fontSize: 22,
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
      <div ref={drawingRef} style={{ position:"relative", width:"100%", height: molHeight, overflow:"hidden" }}>
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
    <div style={{ border:"1px solid #e5e7eb", borderRadius:8, background:"#fff", overflow:"visible", padding:8 }}>
      <div style={{ display:"grid", gridTemplateColumns: show3D?"1fr 1fr":"1fr", gap:12, height: molHeight }}>
        {/* 2D */}
        <div ref={drawingRef} style={{ border:"1px solid #e5e7eb", borderRadius:8, background:"#fff", overflow:"hidden", position:"relative", minWidth:0 }}>
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
            ["H1", safeHex('H1')],
            ["C13", safeHex('C13')],
            ["COSY", safeHex('COSY')],
            ["TOCSY", safeHex('TOCSY')],
            ["HSQC", safeHex('HSQC')],
          ].filter(([k]) => isMethodAvailable(k) && multi?.rowsByType?.[k === "H1" ? "1H" : k === "C13" ? "13C" : k]).map(([k,hex]) => (
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

      {/* Molecule rendering errors. */}
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
              <div style={{ fontWeight:700, marginBottom:6 }}>Quick solutions:</div>
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


/*   Metrics from a CSV row*/

const extractMetricsFromRow = (row, specType = null) => {
  if (!row) return {};
  const metrics = {};
  // For 1H and 13C, display ransy_score.
  // For COSY and HSQC, display coverage_ransy.
  const use1D = specType === "1H" || specType === "13C";
  if (use1D && row.ransy_score !== undefined && row.ransy_score !== "") {
    const v = Number(row.ransy_score || 0).toFixed(1);
    metrics.ransy_score = v;
  } else if (row.coverage_ransy !== undefined && row.coverage_ransy !== "") {
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

/*    Molecule grid modal*/

function stripMethodPrefix(name) {
  if (!name) return name;
  return String(name).replace(/^(?:1H|13C|COSY|TOCSY|HSQC)[\s_\-:]+/i, "");
}

function ShowAllMoleculesModal({
  isOpen,
  onClose,
  sources = [],
  sourcesWithScore = [],
  page = 1,
  setPage,
  pageSize = 20,
  setPageSize,
  getRowsByType,
  csvUrl,
  taskId,
  customColors,
}) {
  const [pdfExporting, setPdfExporting] = useState(false);
  const [pdfProgress, setPdfProgress] = useState(0);
  const [showPdfOptions, setShowPdfOptions] = useState(false);
  const [pdfCount, setPdfCount] = useState(50);
  const [maximized, setMaximized] = useState(false);
  const gridRef = useRef(null);

  // Safe accessors — customColors may come from sessionStorage with missing fields
  const defaultRgb = { H1: [0, 0.72, 0.66], C13: [0.23, 0.05, 0.64], COSY: [1, 0.69, 0], HSQC: [0.97, 0.15, 0.52] };
  const safeRgb = (key) => customColors?.[key]?.rgb || defaultRgb[key] || [0, 0.5, 0.5];


  const totalPages = Math.max(1, Math.ceil(sources.length / pageSize));
  const clampedPage = Math.min(Math.max(1, page), totalPages);
  const start = (clampedPage - 1) * pageSize;
  const end = Math.min(sources.length, start + pageSize);
  const pageSources = sources.slice(start, end);

  const scoreMap = useMemo(() => {
    const m = new Map();
    for (const { source, avgCoverage } of sourcesWithScore) {
      m.set(source, avgCoverage);
    }
    return m;
  }, [sourcesWithScore]);

  const exportPDF = useCallback(async () => {
    setPdfExporting(true);
    setPdfProgress(0);
    setShowPdfOptions(false);

    const savedPage = page; // remember current page to restore later

    try {
      const toExport = sourcesWithScore.slice(0, pdfCount);
      const COLS = 5;
      const ROWS_PER_PAGE = 5;
      const PER_PAGE = COLS * ROWS_PER_PAGE;
      const PAGE_W = 210; // A4 mm
      const PAGE_H = 297;
      const MARGIN = 6;
      const GAP = 2; // mm gap between cells
      const HEADER_H = 12; // mm for page title
      const NAME_H = 5;  // mm for molecule name above structure
      const SCORE_H = 8; // mm for score lines below structure
      const CELL_W = (PAGE_W - MARGIN * 2 - GAP * (COLS - 1)) / COLS;
      const CELL_H = (PAGE_H - MARGIN * 2 - HEADER_H - GAP * (ROWS_PER_PAGE - 1)) / ROWS_PER_PAGE;
      const MOL_H = CELL_H - NAME_H - SCORE_H;

      const doc = new jsPDF({ unit: "mm", format: "a4" });

      // Title page header
      doc.setFontSize(13);
      doc.setFont("helvetica", "bold");
      doc.text("NMRScout - Molecule Report", PAGE_W / 2, MARGIN + 4, { align: "center" });
      doc.setFontSize(8);
      doc.setFont("helvetica", "normal");
      doc.setTextColor(80);
      doc.text(`${toExport.length} molecules • sorted by RANSY average score`, PAGE_W / 2, MARGIN + 8, { align: "center" });
      doc.setTextColor(0);

      // Generate SVGs from RDKit for PDF export.
      const RDKit = await loadRDKitFromCDN();
      const SVG_PX = 400; // RDKit SVG dimensions.

      let posInPdfPage = 0;

      for (let globalIdx = 0; globalIdx < toExport.length; globalIdx++) {
        const { source, avgCoverage } = toExport[globalIdx];
        const name = stripMethodPrefix(source.replace(/_predicted$/i, ""));

        // Start a PDF page when the grid is full.
        if (posInPdfPage === PER_PAGE) {
          doc.addPage();
          posInPdfPage = 0;
          doc.setFontSize(7);
          doc.setFont("helvetica", "normal");
          doc.setTextColor(150);
          doc.text(`NMRScout - Page ${doc.getNumberOfPages()}`, PAGE_W - MARGIN, MARGIN - 1, { align: "right" });
          doc.setTextColor(0);
        }

        const col = posInPdfPage % COLS;
        const row = Math.floor(posInPdfPage / COLS);
        const x = MARGIN + col * (CELL_W + GAP);
        const cellY = MARGIN + HEADER_H + row * (CELL_H + GAP);

        const molY = cellY;

        // Render the molecule and highlights with RDKit.
        try {
          const resolved = await resolveMolForSource(csvUrl, source, taskId);
          const molBlock = resolved?.molblock;
          if (molBlock) {
            // Preserve explicit hydrogens to keep matched_atom_indices aligned with the molblock.
            const rdkMol = RDKit.get_mol(molBlock);
            if (rdkMol) {
              const PDF_COLORS = {
                H1:   safeRgb('H1'),
                C13:  safeRgb('C13'),
                COSY: safeRgb('COSY'),
                HSQC: safeRgb('HSQC'),
              };
              let svgStr;
              try {
                const rowsByType = getRowsByType(source);

                // rowsByType keys: "1H", "13C", "COSY", "HSQC"
                // matched_atom_indices are H-inclusive (backend uses AddHs)
                const getIndices = (key) => {
                  const str = rowsByType[key]?.matched_atom_indices || '';
                  if (!str) return [];
                  return str.split(/\|\||\|/).filter(Boolean)
                    .flatMap(p => p.split(':')[0].split('~').map(s => parseInt(s.trim(),10)).filter(Number.isFinite));
                };

                const atomColors = {};
                const order = [['13C', PDF_COLORS.C13], ['1H', PDF_COLORS.H1], ['HSQC', PDF_COLORS.HSQC], ['COSY', PDF_COLORS.COSY]];
                for (const [key, color] of order)
                  for (const i of getIndices(key)) if (!atomColors[String(i)]) atomColors[String(i)] = color;

                const bondColors = {};
                const allMatched = new Set(Object.keys(atomColors).map(Number));
                const origLines = molBlock.split(/\r?\n/);
                const nAtO = parseInt((origLines[3]||'').slice(0,3).trim(),10)||0;
                const nBdO = parseInt((origLines[3]||'').slice(3,6).trim(),10)||0;
                for (let bi=0; bi<nBdO; bi++) {
                  const l = origLines[4+nAtO+bi]||'';
                  const a1 = parseInt(l.slice(0,3).trim(),10)-1;
                  const a2 = parseInt(l.slice(3,6).trim(),10)-1;
                  if (allMatched.has(a1) && allMatched.has(a2))
                    bondColors[String(bi)] = atomColors[String(a1)] || atomColors[String(a2)];
                }

                const highlightAtoms = Object.keys(atomColors).map(Number);
                const highlightBonds = Object.keys(bondColors).map(Number);
                svgStr = rdkMol.get_svg_with_highlights(JSON.stringify({
                  width: SVG_PX, height: SVG_PX,
                  atoms: highlightAtoms,
                  highlight_atom_colors: atomColors,
                  bonds: highlightBonds,
                  highlight_bond_colors: bondColors,
                  highlightRadius: 0.2,
                  atomHighlightsAreCircles: true,
                  highlightBondWidthMultiplier: 8,
                }));
              } catch {
                svgStr = rdkMol.get_svg(SVG_PX, SVG_PX);
              }
              rdkMol.delete();

              // Rasterize the SVG on a canvas before adding it to the PDF.
              await new Promise((resolve) => {
                const img = new Image();
                const blob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
                const url = URL.createObjectURL(blob);
                img.onload = () => {
                  const canvas = document.createElement('canvas');
                  canvas.width = SVG_PX;
                  canvas.height = SVG_PX;
                  const ctx = canvas.getContext('2d');
                  ctx.fillStyle = '#ffffff';
                  ctx.fillRect(0, 0, SVG_PX, SVG_PX);
                  ctx.drawImage(img, 0, 0, SVG_PX, SVG_PX);
                  URL.revokeObjectURL(url);
                  const png = canvas.toDataURL('image/png');
                  doc.addImage(png, 'PNG', x, molY, CELL_W, MOL_H);
                  resolve();
                };
                img.onerror = () => { URL.revokeObjectURL(url); resolve(); };
                img.src = url;
              });
            }
          }
        } catch (pdfErr) { console.error('[PDF cell error]', source, pdfErr); }

        // Molecule name below the structure.
        const label = name.length > 20 ? name.slice(0, 19) + "…" : name;
        const textY = molY + MOL_H;
        doc.setFontSize(6);
        doc.setFont("helvetica", "bolditalic");
        doc.setTextColor(20);
        doc.text(label, x + CELL_W / 2, textY + 3.5, { align: "center" });

        // NMR score.
        doc.setFontSize(5.5);
        doc.setFont("helvetica", "normal");
        doc.setTextColor(60);
        doc.text(`NMR score: ${avgCoverage >= 0 ? avgCoverage.toFixed(0) : "N/A"}`, x + CELL_W / 2, textY + 7, { align: "center" });
        doc.setTextColor(0);

        posInPdfPage++;
        setPdfProgress(Math.round(((globalIdx + 1) / toExport.length) * 100));
      }

      // Page numbers footer
      const totalPdfPages = doc.getNumberOfPages();
      for (let p = 1; p <= totalPdfPages; p++) {
        doc.setPage(p);
        doc.setFontSize(7);
        doc.setTextColor(150);
        doc.text(`${p} / ${totalPdfPages}`, PAGE_W / 2, PAGE_H - 4, { align: "center" });
      }

      doc.save(`NMRScout_molecules_${toExport.length}.pdf`);
    } catch (e) {
      console.error("[PDF export error]", e);
      alert("PDF export failed: " + (e?.message || String(e)));
    } finally {
      // Restore original modal page
      setPage(savedPage);
      setPdfExporting(false);
      setPdfProgress(0);
    }
  }, [sourcesWithScore, pdfCount, page, setPage, csvUrl, taskId, getRowsByType, customColors]);

  if (!isOpen) return null;

  return (
    <div
      style={{ position:"fixed", inset:0, zIndex:10000, background:"rgba(0,0,0,0.55)",
               display:"flex", alignItems:"center", justifyContent:"center" }}
      onClick={onClose}
    >
      <div
        style={{
          width: maximized ? "100vw" : "min(1100px, 96vw)",
          maxHeight: maximized ? "100vh" : "92vh",
          height: maximized ? "100vh" : undefined,
          background:"#fff",
          borderRadius: maximized ? 0 : 12,
          boxShadow: maximized ? "none" : "0 24px 80px rgba(0,0,0,.35)",
          display:"flex", flexDirection:"column",
          transition: "all 0.2s ease",
        }}
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

            {/* PDF Export */}
            <div style={{ position:"relative" }}>
              <button
                onClick={() => setShowPdfOptions(v => !v)}
                disabled={pdfExporting}
                style={{ padding:"6px 12px", borderRadius:8, border:"1px solid #3b82f6",
                         background: pdfExporting ? "#e2e8f0" : "#3b82f6", color:"#fff",
                         cursor: pdfExporting ? "not-allowed" : "pointer", fontWeight:600, fontSize:12 }}
              >
                {pdfExporting ? `Exporting… ${pdfProgress}%` : "⬇ PDF"}
              </button>
              {showPdfOptions && !pdfExporting && (
                <div style={{ position:"absolute", right:0, top:"110%", background:"#fff",
                              border:"1px solid #cbd5e1", borderRadius:10, padding:14,
                              boxShadow:"0 8px 24px rgba(0,0,0,0.15)", zIndex:100, minWidth:220 }}>
                  <div style={{ fontSize:12, fontWeight:600, marginBottom:8 }}>Export PDF</div>
                  <label style={{ fontSize:12, color:"#475569", display:"flex", flexDirection:"column", gap:4 }}>
                    Number of molecules to export:
                    <select
                      value={pdfCount}
                      onChange={e => setPdfCount(+e.target.value)}
                      style={{ padding:"4px 8px", borderRadius:6, border:"1px solid #cbd5e1", marginTop:2 }}
                    >
                      {[25, 50, 100, 200, sources.length].filter((v,i,a) => a.indexOf(v)===i && v > 0).map(n => (
                        <option key={n} value={n}>{n === sources.length ? `All (${n})` : n}</option>
                      ))}
                    </select>
                  </label>
                  <div style={{ fontSize:11, color:"#94a3b8", margin:"6px 0" }}>
                    {Math.ceil(pdfCount / 25)} PDF page(s) • 5×5 grid
                  </div>
                  <button
                    onClick={exportPDF}
                    style={{ width:"100%", padding:"7px 0", borderRadius:8, border:"none",
                             background:"#3b82f6", color:"#fff", fontWeight:700, fontSize:13, cursor:"pointer" }}
                  >
                    Generate PDF
                  </button>
                </div>
              )}
            </div>

            <button onClick={() => setMaximized(m => !m)}
              title={maximized ? "Restore" : "Fullscreen"}
              style={{ padding:"6px 10px", borderRadius:8, border:"1px solid #cbd5e1", background:"#fff", cursor:"pointer", fontWeight:700, fontSize:14 }}>
              {maximized ? "⧉" : "⛶"}
            </button>
            <button onClick={onClose}
              style={{ padding:"6px 10px", borderRadius:8, border:"1px solid #cbd5e1", background:"#fff", cursor:"pointer", fontWeight:700 }}>
              Close
            </button>
          </div>
        </div>

        {/* Pagination top */}
        <div style={{ padding:"8px 12px", display:"flex", gap:8, alignItems:"center", borderBottom:"1px solid #eef2f7" }}>
          <button aria-label="Previous page" title="Previous page" onClick={()=>setPage(p=>Math.max(1,p-1))} disabled={clampedPage<=1}
                  style={{ padding:"4px 8px", borderRadius:6, border:"1px solid #cbd5e1", background:"#fff", cursor:"pointer" }}>&lt;</button>
          <div style={{ fontSize:12, color:"#334155" }}>
            Page {clampedPage} / {totalPages} &nbsp;•&nbsp; {sources.length} molecules
          </div>
          <button aria-label="Next page" title="Next page" onClick={()=>setPage(p=>Math.min(totalPages,p+1))} disabled={clampedPage>=totalPages}
                  style={{ padding:"4px 8px", borderRadius:6, border:"1px solid #cbd5e1", background:"#fff", cursor:"pointer" }}>&gt;</button>
        </div>


        {/* Visual molecules grid */}
        <div style={{ flex:1, overflow:"auto", padding:"12px" }}>
          <div ref={gridRef} style={{
            display:"grid",
            gridTemplateColumns:"repeat(auto-fill, minmax(180px, 1fr))",
            gap:12
          }}>
            {pageSources.map((source, idx) => {
              return (
                <div key={source} data-mol-source={source} className="nmr-molecule-card" style={{
                  borderRadius:8,
                  padding:8,
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
                  <div className="nmr-molecule-card__structure" style={{ height:140, marginBottom:8 }}>
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
                  
                  <MoleculeScoreSummary
                    source={source}
                    rowsByType={getRowsByType(source)}
                    averageScore={scoreMap.get(source)}
                  />
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



/*   Results viewer*/

export default function Viewer({
  csvs: allCsvs,
  taskId = null,
  defaultShow3D = true,
  annotationResults = {},
  pdataDirs: allPdataDirs = {},
  bgOptions = {},
  customColors: propsCustomColors,
  onColorChange,
  // Optional external control of the selected source.
  externalSource = null,
  hideSidebar = false,
  coverageMethods: allCoverageMethods = ["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable),
}) {
  const csvs = useMemo(() => availableMethodEntries(allCsvs), [allCsvs]);
  const pdataDirs = useMemo(() => availableMethodEntries(allPdataDirs), [allPdataDirs]);
  const coverageMethods = useMemo(() => allCoverageMethods.filter(isMethodAvailable), [allCoverageMethods]);

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
  const [activeSpectrum, setActiveSpectrum] = useState("all");
  const [showExp, setShowExp] = useState(true);
  const [reverse1D, setReverse1D] = useState(true);
  const [show3D, setShow3D] = useState(defaultShow3D);

  const [zoomRectEnabled, setZoomRectEnabled] = useState(true);
  const [colorPickerOpen, setColorPickerOpen] = useState(false);
  const [relevel2DEnabled, setRelevel2DEnabled] = useState(false); // Disabled by default - only zoom updates background

  // A nonce per spectrum type forces canvas remounts.
  const [refreshKeyByType, setRefreshKeyByType] = useState({});
  const bumpKey = useCallback((type) => {
    setRefreshKeyByType(k => ({ ...k, [type]: (k[type] ?? 0) + 1 }));
  }, []);

  const presentTypes = useMemo(
    () => Object.keys(csvs || {}).filter(Boolean),
    [csvs]
  );

  const hasAnyData = useMemo(() => presentTypes.length > 0, [presentTypes]);

  const [bgByType, setBgByType] = useState(
    () => annotationResults || {}
  );
  const [busyByType, setBusyByType] = useState({});

  // State for raw spectrum data (for direct canvas rendering without PNG)
  const [rawDataByType, setRawDataByType] = useState({
    "1H": null,
    "13C": null,
    COSY: null,
    TOCSY: null,
    HSQC: null,
  });
  
  // State for current zoom range (ppm) per type
  const [zoomRangeByType, setZoomRangeByType] = useState({
    "1H": null,   // {min, max} for 1D
    "13C": null,  // {min, max} for 1D
    COSY: null,   // {f1: {min, max}, f2: {min, max}} for 2D
    HSQC: null,   // {f1: {min, max}, f2: {min, max}} for 2D
  });

  // State for intensity scale (controlled by mouse wheel)
  useSpectrumWindowSettings(bgOptions, setZoomRangeByType);

  // TopSpin-style: scale controls how many contour levels are visible
  // Uses logarithmic mapping: log10(scale) determines threshold
  // Full range from 0.0001 to 10000:
  // - Scale = 0.0001 → 0% of contours (nothing visible)
  // - Scale = 0.01 → ~25% of contours
  // - Scale = 1.0 → ~50% of contours (center)
  // - Scale = 100 → ~75% of contours  
  // - Scale = 10000 → 100% of contours (all noise visible)
  const [intensityScaleByType, setIntensityScaleByType] = useState({
    "1H": 1.0,
    "13C": 1.0,
    COSY: 1.0,   // Start at center - scroll up for more, down for less
    TOCSY: 1.0,
    HSQC: 1.0,   // Start at center - scroll up for more, down for less
  });

  // Local state for resolution and level parameters
  const [resByType, setResByType] = useState({
    "1H":  { width_px: 1200, height_px: 500, dpi: 150 },
    "13C": { width_px: 1200, height_px: 500, dpi: 150 },
    COSY:  { width_px: 1000, height_px: 800, dpi: 150 },
    TOCSY: { width_px: 1000, height_px: 800, dpi: 150 },
    HSQC:  { width_px: 1000, height_px: 800, dpi: 150 },
  });

  const [levelByType, setLevelByType] = useState({
    COSY: 0.005,
    TOCSY: 0.005,
    HSQC: 0.005
  });

  
  // Snapshots are keyed by spectrum type and source: { [type]: { [source]: snap } }.
  const originalBgRef = useRef({});

  useEffect(() => {
    const snap = {};
    for (const k of ["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable)) {
      if (annotationResults?.[k]) snap[k] = annotationResults[k];
    }
    setBgByType(prev => ({ ...prev, ...snap }));
  }, [annotationResults]);

  // When initial results arrive, snap according to current source
  useEffect(() => {
    if (!selectedSource) return;
    for (const k of ["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable)) {
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
    for (const k of ["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable)) {
      const bg = bgByType?.[k];
      if (bg) {
        originalBgRef.current[k] = originalBgRef.current[k] || {};
        if (!originalBgRef.current[k][selectedSource]) {
          originalBgRef.current[k][selectedSource] = bg;
        }
      }
    }
  }, [selectedSource, bgByType]);

  // Ref to track which pdata paths have been loaded (prevents infinite loop)
  // Stores the path that was loaded, not just a boolean
  const loadedRawDataRef = useRef({});

  // Effect to load raw spectrum data when pdataDirs change
  useEffect(() => {
    const loadRawData = async () => {
      const types1D = ["1H", "13C"];
      const types2D = ["COSY", "TOCSY", "HSQC"].filter(isMethodAvailable);

      // Load 1D spectra
      for (const type of types1D) {
        const pd = pdataDirs?.[type];
        // Reload if path changed or if never loaded
        const pathChanged = pd && loadedRawDataRef.current[type] !== pd;
        if (pd && pathChanged) {
          loadedRawDataRef.current[type] = pd; // Store the loaded path
          try {
            setBusyByType(prev => ({ ...prev, [type]: true }));
            const data = await extractBrukerSpectrum(pd);
            if (data?.status === "ok") {
              setRawDataByType(prev => ({
                ...prev,
                [type]: {
                  ppm: data.ppm,
                  intensity: data.intensity,
                  sf: data.sf
                }
              }));
              // Initialize zoom: use bgOptions bounds if available, otherwise full data range
              const opts = bgOptions?.[type] || {};
              const userPpmMin = Number(opts.ppm_min);
              const userPpmMax = Number(opts.ppm_max);
              
              if (Number.isFinite(userPpmMin) && Number.isFinite(userPpmMax)) {
                // Use user-defined bounds from annotation parameters
                setZoomRangeByType(prev => ({
                  ...prev,
                  [type]: {
                    min: Math.min(userPpmMin, userPpmMax),
                    max: Math.max(userPpmMin, userPpmMax)
                  }
                }));
              } else if (data.ppm?.length) {
                // Fall back to full data range
                const ppmMin = data.ppm.reduce((a, b) => a < b ? a : b, data.ppm[0]);
                const ppmMax = data.ppm.reduce((a, b) => a > b ? a : b, data.ppm[0]);
                setZoomRangeByType(prev => ({
                  ...prev,
                  [type]: {
                    min: ppmMin,
                    max: ppmMax
                  }
                }));
              }
            }
          } catch (err) {
            console.error(`Error loading raw ${type} data:`, err);
          } finally {
            setBusyByType(prev => ({ ...prev, [type]: false }));
          }
        }
      }
      
      // Load 2D spectra
      for (const type of types2D) {
        const pd = pdataDirs?.[type];
        const opts = bgOptions?.[type] || {};

        // Create cache key including both path AND bgOptions
        // This ensures reload when axis windows change, not just when path changes
        const cacheKey = JSON.stringify({
          path: pd,
          f1_min: opts.f1_ppm_min,
          f1_max: opts.f1_ppm_max,
          f2_min: opts.f2_ppm_min,
          f2_max: opts.f2_ppm_max
        });

        const needsReload = pd && loadedRawDataRef.current[type] !== cacheKey;
        if (pd && needsReload) {
          loadedRawDataRef.current[type] = cacheKey; // Store cache key with path + options
          try {
            setBusyByType(prev => ({ ...prev, [type]: true }));
            // Use default lvl of 0.005 for initial load
            const lvl = 0.005;
            const data = await extractBruker2D(pd, {
              levelStart: lvl,   // Use lvl as levelStart (0.005 = 0.5% of max)
              nLevels: 60,
              decimate: 2
            });
            if (data?.status === "ok") {
              setRawDataByType(prev => ({
                ...prev,
                [type]: {
                  contours: data.contours,
                  f1_range: data.f1_range,
                  f2_range: data.f2_range
                }
              }));
              // Initialize zoom: use bgOptions bounds if available, otherwise full data range
              // opts already defined above in the loop scope
              const userF1Min = Number(opts.f1_ppm_min);
              const userF1Max = Number(opts.f1_ppm_max);
              const userF2Min = Number(opts.f2_ppm_min);
              const userF2Max = Number(opts.f2_ppm_max);
              
              const hasUserF1 = Number.isFinite(userF1Min) && Number.isFinite(userF1Max);
              const hasUserF2 = Number.isFinite(userF2Min) && Number.isFinite(userF2Max);
              
              if (hasUserF1 || hasUserF2) {
                // Use user-defined bounds from annotation parameters
                setZoomRangeByType(prev => ({
                  ...prev,
                  [type]: {
                    f1: hasUserF1 
                      ? { min: Math.min(userF1Min, userF1Max), max: Math.max(userF1Min, userF1Max) }
                      : { min: data.f1_range[0], max: data.f1_range[1] },
                    f2: hasUserF2
                      ? { min: Math.min(userF2Min, userF2Max), max: Math.max(userF2Min, userF2Max) }
                      : { min: data.f2_range[0], max: data.f2_range[1] }
                  }
                }));
              } else if (data.f1_range && data.f2_range) {
                // Fall back to full data range
                setZoomRangeByType(prev => ({
                  ...prev,
                  [type]: {
                    f1: { min: data.f1_range[0], max: data.f1_range[1] },
                    f2: { min: data.f2_range[0], max: data.f2_range[1] }
                  }
                }));
              }
            }
          } catch (err) {
            console.error(`Error loading raw ${type} data:`, err);
          } finally {
            setBusyByType(prev => ({ ...prev, [type]: false }));
          }
        }
      }
    };
    
    loadRawData();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [pdataDirs, bgOptions]);

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
            const parsed = parseCSVSmart(text);
            const data = homonuclearSpectrumRows(parsed.data || [], t);
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

    const parseCoverage = (row, type) => {
      if (!row) return null;
      // Rank 1H and 13C by intersection_coverage, the STOCSY group intersection score.
      const use1D = type === "1H" || type === "13C";
      if (use1D) {
        if (row.intersection_coverage !== undefined && row.intersection_coverage !== null && row.intersection_coverage !== "") {
          const v = parseFloat(String(row.intersection_coverage).replace(",", "."));
          if (Number.isFinite(v)) return v;
        }
        return 0; // Missing intersection coverage contributes a zero score.
      }
      // Rank COSY and HSQC by coverage_ransy.
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
      let weightedSum = 0;
      let totalWeight = 0;

      for (const type of ["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable)) {
        const row = mapsByType[type]?.get(source);
        const cov = parseCoverage(row, type);
        if (cov != null) {
          weightedSum += cov;
          totalWeight += 1;
        }
      }
      const avg = totalWeight > 0 ? weightedSum / totalWeight : -1;
      return { source, avgCoverage: avg };
    });

    withCov.sort((a, b) => (b.avgCoverage - a.avgCoverage) || a.source.localeCompare(b.source));
    arr = withCov.map(x => x.source);

    if (!query) return arr;
    const q = query.toLowerCase();
    return arr.filter(s => s.toLowerCase().includes(q));
  }, [mapsByType, query]);

  // Select the highest-ranked molecule after loading.
  const hasAutoSelectedRef = useRef(false);
  useEffect(() => {
    // Select the first molecule once after loading.
    if (!loading && allSources.length > 0 && !selectedSource && !hasAutoSelectedRef.current) {
      debugLog('[Viewer] Auto-selecting first molecule:', allSources[0]);
      setSelectedSource(allSources[0]);
      hasAutoSelectedRef.current = true;
    }
    // Reset selection tracking when a new source set is loaded.
    if (allSources.length === 0) {
      hasAutoSelectedRef.current = false;
    }
  }, [loading, allSources, selectedSource, setSelectedSource]);

  //  Global sorted list (without search filter) for the modal
  const sortedSourcesWithScore = useMemo(() => {
    const union = new Set();
    Object.values(mapsByType).forEach(m => {
      for (const s of m.keys()) union.add(s);
    });
    const arr = [...union];

    const parseCoverage = (row, type) => {
      if (!row) return null;
      // Rank 1H and 13C by intersection_coverage, the STOCSY group intersection score.
      const use1D = type === "1H" || type === "13C";
      if (use1D) {
        if (row.intersection_coverage !== undefined && row.intersection_coverage !== null && row.intersection_coverage !== "") {
          const v = parseFloat(String(row.intersection_coverage).replace(",", "."));
          if (Number.isFinite(v)) return v;
        }
        return 0; // Missing intersection coverage contributes a zero score.
      }
      // Rank COSY and HSQC by coverage_ransy.
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
      let weightedSum = 0;
      let totalWeight = 0;

      const methodsToUse = coverageMethods.length > 0 ? coverageMethods : ["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable);
      for (const type of methodsToUse) {
        const row = mapsByType[type]?.get(source);
        const v = parseCoverage(row, type);
        if (v != null) {
          weightedSum += v;
          totalWeight += 1;
        }
      }
      const avg = totalWeight > 0 ? weightedSum / totalWeight : -1;
      return { source, avgCoverage: avg };
    });

    withCov.sort((a,b) => (b.avgCoverage - a.avgCoverage) || a.source.localeCompare(b.source));
    return withCov;
  }, [mapsByType, coverageMethods]);

  const sortedSourcesAll = useMemo(
    () => sortedSourcesWithScore.map(x => x.source),
    [sortedSourcesWithScore]
  );

  // Rows grouped by spectrum type for the molecule grid.
  const rowsByTypeForSource = useCallback((source) => ({
    "1H":   mapsByType["1H"]?.get(source),
    "13C":  mapsByType["13C"]?.get(source),
    "COSY": mapsByType["COSY"]?.get(source),
    "TOCSY": mapsByType["TOCSY"]?.get(source),
    "HSQC": mapsByType["HSQC"]?.get(source),
  }), [mapsByType]);

  const rowFor = useCallback(
    (t) => mapsByType[t]?.get(selectedSource) || null,
    [mapsByType, selectedSource]
  );

  const availableSpectrumTypes = useMemo(
    () => ["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable).filter(type => rowFor(type) != null),
    [rowFor]
  );
  const displayedSpectrum = availableSpectrumTypes.includes(activeSpectrum) ? activeSpectrum : "all";

  useEffect(() => {
    if (activeSpectrum !== "all" && !availableSpectrumTypes.includes(activeSpectrum)) {
      setActiveSpectrum("all");
    }
  }, [activeSpectrum, availableSpectrumTypes]);

  const oneH = useMemo(() => {
    const r = rowFor("1H");
    return r ? {
      theory: (r.ppm_molecule || "")
        .split("|")
        .map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      expX: (r.intersection_ppms || r.RANSY_signals || "")
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
      metrics: extractMetricsFromRow(r, "1H"),
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
      expX: (r.intersection_ppms || r.RANSY_signals || "")
        .split("|").map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      expI: (r.RANSY_mean_intensities || r.intensity_ratios || "")
        .split(";").map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      matches: (r.ppm_molecule_matched || "")
        .split("|").map(v => parseFloat(String(v).replace(",", ".")))
        .filter(Number.isFinite),
      metrics: extractMetricsFromRow(r, "13C"),
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
      metrics: extractMetricsFromRow(r, "COSY"),
    } : null;
  }, [rowFor]);

  const tocsy = useMemo(() => {
    const r = rowFor("TOCSY");
    return r ? {
      theory: (r.ppm_molecule || "")
        .split("||").map(p => p.split("|").map(v => parseFloat(String(v).replace(",", "."))))
        .filter(a => a.length === 2 && a.every(Number.isFinite)),
      expPts: (r.RANSY_signals || "")
        .split("||").map(p => p.split("|").map(v => parseFloat(String(v).replace(",", "."))))
        .filter(a => a.length === 2 && a.every(Number.isFinite)),
      expI: (r.RANSY_mean_intensities || r.intensity_ratios || "")
        .split(";").map(v => parseFloat(String(v).replace(",", "."))).filter(Number.isFinite),
      pairs: (r.ppm_molecule_matched || "")
        .split("||").map(p => p.split("|").map(v => parseFloat(String(v).replace(",", "."))))
        .filter(a => a.length === 2 && a.every(Number.isFinite)),
      metrics: extractMetricsFromRow(r, "TOCSY"),
      maxBonds: r.MaxBonds,
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
      metrics: extractMetricsFromRow(r, "HSQC"),
    } : null;
  }, [rowFor]);

  const pivotCsvUrl = useMemo(
    () => csvs["1H"] || csvs["13C"] || csvs["COSY"] || csvs["TOCSY"] || csvs["HSQC"] || null,
    [csvs]
  );

  const getMeta = (obj) => obj?.background?.meta || obj?.meta || null;

  // eslint-disable-next-line react-hooks/exhaustive-deps
  const setBusy = useCallback((type, v) => setBusyByType(prev => ({ ...prev, [type]: !!v })));
  
  // Invalidate cached images and update the snapshot for this type and source.
  // eslint-disable-next-line react-hooks/exhaustive-deps
  const storeResult = useCallback((type, res) => {
    if (!res) return;

    // Cache invalidation token.
    const nonce = `${Date.now()}-${Math.random().toString(36).slice(2)}`;
    res.meta = { ...(res.meta || {}), nonce };

    // updates current state (display)
    setBgByType(prev => ({ ...prev, [type]: res }));

    // Synchronize the level control for 2D spectra.
    if ((type === "COSY" || type === "TOCSY" || type === "HSQC") &&
        Number.isFinite(+res?.meta?.level_start)) {
      setLevelByType(p => ({ ...p, [type]: +res.meta.level_start }));
    }

    //  initial snapshot by (type, source): only write if absent
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

  
  // Bound the time spent applying settings.
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
        const fn = type === "COSY" || type === "TOCSY" ? interactiveCOSY : interactiveHSQC;
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

  const canInteract1H   = !!pdataDirs?.["1H"];
  const canInteract13C  = !!pdataDirs?.["13C"];
  const canInteractCOSY = !!pdataDirs?.["COSY"];
  const canInteractTOCSY = !!pdataDirs?.["TOCSY"];
  const canInteractHSQC = !!pdataDirs?.["HSQC"];

  // Molecule grid modal state.
  const [showAllOpen, setShowAllOpen] = useState(false);
  const [allPage, setAllPage] = useState(1);
  const [allPageSize, setAllPageSize] = useState(20);

  // Detailed 1H spectrum modal state.
  const [spectrumDetailOpen, setSpectrumDetailOpen] = useState(false);

  // ====
  // Batch loading of simulated spectra
  // ====
  // Per-molecule cache populated by batch simulation.
  const simulatedSpectrumCache = useRef({});
  // Simulation data for the selected molecule.
  const [preloadedMultiplicityData, setPreloadedMultiplicityData] = useState(null);
  const [preloadedSpectrumData, setPreloadedSpectrumData] = useState(null);
  const [preloadingSpectrum, setPreloadingSpectrum] = useState(false);
  const [preloadingMultiplicity, setPreloadingMultiplicity] = useState(false);
  // Track whether batch simulation has started for this task.
  const batchLoadedForTask = useRef(null);

  // Load simulated spectra in a batch when taskId changes.
  useEffect(() => {
    if (!taskId || batchLoadedForTask.current === taskId) {
      return;
    }

    debugLog('[Viewer] Starting batch spectrum simulation for task:', taskId);
    batchLoadedForTask.current = taskId;
    setPreloadingSpectrum(true);

    // Clear spectra cached for the previous task.
    simulatedSpectrumCache.current = {};

    simulateSpectrumBatch({
      task_id: taskId,
      freq_mhz: 700,
      line_shape: 'lorentzian',
      line_width_hz: 1.5,
      enable_dynamics: true,
      solvent: 'CDCl3',
      temperature_k: 298,
      ph: 7.4,
    })
      .then(batchResult => {
        debugLog('[Viewer] Batch simulation complete:', batchResult.total_molecules, 'ok:', batchResult.successful);

        // Cache the simulated spectra.
        if (batchResult.spectra) {
          for (const spec of batchResult.spectra) {
            if (spec.mol_id && spec.status === 'ok') {
              simulatedSpectrumCache.current[spec.mol_id] = spec;
            }
          }
        }

        // Load the spectrum for the current selection, if any.
        if (selectedSource && simulatedSpectrumCache.current[selectedSource]) {
          setPreloadedSpectrumData(simulatedSpectrumCache.current[selectedSource]);
        }
      })
      .catch(err => {
        console.warn('[Viewer] Batch simulation failed:', err);
        // Fall back to individual requests if batch loading fails.
        batchLoadedForTask.current = null;
      })
      .finally(() => {
        setPreloadingSpectrum(false);
      });
  }, [taskId, selectedSource]);

  // Use cached data for the selected molecule, or request it separately.
  useEffect(() => {
    if (!taskId || !selectedSource) {
      setPreloadedMultiplicityData(null);
      setPreloadedSpectrumData(null);
      return;
    }

    const molId = selectedSource;
    
    // Source name variants used for cache lookup.
    const molIdVariants = [
      molId,
      molId.replace(/_predicted$/i, ''),
      molId.replace(/_annotated$/i, ''),
    ];

    // Try each source name variant in the cache.
    let cachedSpectrum = null;
    for (const variant of molIdVariants) {
      if (simulatedSpectrumCache.current[variant]) {
        cachedSpectrum = simulatedSpectrumCache.current[variant];
        debugLog('[Viewer] Using cached spectrum for:', variant);
        break;
      }
    }

    if (cachedSpectrum) {
      setPreloadedSpectrumData(cachedSpectrum);
      setPreloadedMultiplicityData(null); // A cached spectrum does not require a multiplicity request.
      return;
    }

    // Fetch this molecule separately on a cache miss.
    debugLog('[Viewer] Spectrum not in cache, loading individually:', molId);
    setPreloadingMultiplicity(true);
    setPreloadedSpectrumData(null);
    setPreloadedMultiplicityData(null);

    getMultiplicityData(taskId, molId)
      .then(data => {
        const multData = data.multiplicity_data;
        // Ignore responses for molecules that are no longer selected.
        if (selectedSource === molId) {
          setPreloadedMultiplicityData(multData);
        }

        return simulateSpectrum({
          multiplicity_json: multData,
          mol_id: molId,
          freq_mhz: 700,
          line_shape: 'lorentzian',
          line_width_hz: 1.5,
          enable_dynamics: true,
          solvent: 'CDCl3',
          temperature_k: 298,
          ph: 7.4,
        });
      })
      .then(specData => {
        debugLog('[Viewer] Loaded spectrum for:', molId);
        simulatedSpectrumCache.current[molId] = specData;
        // Ignore responses for molecules that are no longer selected.
        if (selectedSource === molId) {
          setPreloadedSpectrumData(specData);
        }
      })
      .catch(err => {
        console.warn('[Viewer] Failed to load spectrum:', err);
      })
      .finally(() => {
        setPreloadingMultiplicity(false);
      });
  }, [taskId, selectedSource]);

  // Adjust spectrum intensity with the mouse wheel.
  // A factor of 1.1 gives gradual scaling: about 100 wheel steps from 1 to 10,000.
  const handleIntensityWheel = useCallback((type, delta) => {
    setIntensityScaleByType(prev => {
      const current = prev[type] || 1.0;
      // delta > 0 = scroll down = decrease (fewer signals), delta < 0 = scroll up = increase (more signals)
      // Factor 1.1 for smoother control (~100 scrolls for full range)
      const factor = delta > 0 ? (1 / 1.1) : 1.1;  // ~0.909 or 1.1
      // Clamp the scale to 0.0001–10000.
      const newScale = Math.max(0.0001, Math.min(10000, current * factor));
      return { ...prev, [type]: newScale };
    });
  }, []);

  /* Spectrum scale and reset controls. */
  const renderControls1D = (type) => {
    const scale = intensityScaleByType[type] || 1.0;
    return (
      <div style={{
        display: "flex", gap: 8, alignItems: "center",
        marginBottom: 6, fontSize: 10, color: "#64748b",
      }}>
        <span>Scale: {scale.toFixed(2)}×</span>
        <button
          onClick={() => setIntensityScaleByType(prev => ({ ...prev, [type]: 1.0 }))}
          style={{
            padding: "2px 6px", borderRadius: 4,
            border: "1px solid #cbd5e1",
            background: "#fff", cursor: "pointer", fontSize: 9,
          }}
        >
          Reset Scale
        </button>
      </div>
    );
  };

  // Force reload of 2D data with specific level
  const forceReload2D = useCallback(async (type, lvlOverride = null) => {
    const pd = pdataDirs?.[type];
    if (!pd) return;
    
    // Use provided level or current levelByType
    const lvl = lvlOverride ?? levelByType[type] ?? 0.005;
    
    // Clear cache for this type
    loadedRawDataRef.current[type] = false;
    
    // Trigger reload with lvl as levelStart
    setBusyByType(prev => ({ ...prev, [type]: true }));
    try {
      const data = await extractBruker2D(pd, {
        levelStart: lvl,  // Use lvl directly as levelStart
        nLevels: 60,
        decimate: 2
      });
      if (data?.status === "ok") {
        setRawDataByType(prev => ({
          ...prev,
          [type]: {
            contours: data.contours,
            f1_range: data.f1_range,
            f2_range: data.f2_range
          }
        }));
        loadedRawDataRef.current[type] = true;
      }
    } catch (err) {
      console.error(`[forceReload2D] Error reloading ${type}:`, err);
    } finally {
      setBusyByType(prev => ({ ...prev, [type]: false }));
    }
  }, [pdataDirs, levelByType]);

  const renderControls2D = (type) => {
    const lvl = levelByType[type] ?? 0.005;
    return (
      <div style={{
        display: "flex", gap: 8, alignItems: "center",
        marginBottom: 6, fontSize: 10, color: "#64748b",
      }}>
        <label style={{ display: "flex", alignItems: "center", gap: 4 }}>
          <span>Lvl:</span>
          <input
            type="number"
            value={lvl}
            onChange={(e) => {
              const val = parseFloat(e.target.value);
              if (!isNaN(val) && val >= 0 && val <= 1) {
                setLevelByType(prev => ({ ...prev, [type]: val }));
              }
            }}
            step="0.001"
            min="0"
            max="1"
            style={{
              width: 60, padding: "2px 4px", fontSize: 10,
              border: "1px solid #cbd5e1", borderRadius: 4,
            }}
            title="Seuil de contour: 0.001 = bruit visible, 0.1 = signaux forts seulement"
          />
        </label>
        <button
          onClick={() => forceReload2D(type, lvl)}
          style={{
            padding: "2px 6px", borderRadius: 4,
            border: "1px solid #2563eb",
            background: "#3b82f6", color: "#fff",
            cursor: "pointer", fontSize: 9, fontWeight: 600,
          }}
          title="Appliquer le niveau et recharger les contours"
        >
          Appliquer
        </button>
        <button
          onClick={() => {
            setLevelByType(prev => ({ ...prev, [type]: 0.005 }));
            forceReload2D(type, 0.005);
          }}
          style={{
            padding: "2px 6px", borderRadius: 4,
            border: "1px solid #cbd5e1",
            background: "#fff", cursor: "pointer", fontSize: 9,
          }}
        >
          Reset
        </button>
      </div>
    );
  };

  /*     Main layout*/
  return (
    <div className="viewer" style={{ marginTop: 0 }}>

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
          {/* Sidebar */}
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
                const stripMethodPrefix = (name) => String(name).replace(/^(?:1H|13C|COSY|TOCSY|HSQC)[\s_\-:]+/i, "");
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
                  canInteractCOSY || canInteractTOCSY || canInteractHSQC)
                  ? "✓ Interactions" : "✗ Pas d'interactions"}
              </div>
            </div>
          </aside>
          )}

          {/* Results */}
          <section style={{ minWidth: 0 }}>
            {/* Molecule and display controls */}
            <div style={{ marginBottom: 12 }}>
              <div
                style={{
                  display: "flex", justifyContent: "flex-end",
                  alignItems: "center", marginBottom: 4,
                }}
              >
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
                      "TOCSY":rowFor("TOCSY"),
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
                  {hasAnyData ? "No source selected" : "No molecule. Run a prediction then an annotation"}
                </div>
              )}
            </div>

            {/* Spectra */}
            <div role="group" aria-label="Displayed spectrum" style={{ display: "flex", flexWrap: "wrap", gap: 8, marginBottom: 12 }}>
              {(availableSpectrumTypes.length ? ["all", ...availableSpectrumTypes] : []).map(type => (
                <button
                  key={type}
                  type="button"
                  aria-pressed={displayedSpectrum === type}
                  onClick={() => setActiveSpectrum(type)}
                  style={{ padding: "6px 12px", borderRadius: 6, cursor: "pointer", border: "1px solid #cbd5e1", background: displayedSpectrum === type ? "#2563eb" : "#fff", color: displayedSpectrum === type ? "#fff" : "#0f172a" }}
                >
                  {type === "all" ? "All spectra" : TYPE_LABELS[type] || type}
                </button>
              ))}
            </div>
            {(() => {
              let blocks = [];

              if (oneH) {
                const raw1H = rawDataByType?.["1H"];
                const zoom1H = zoomRangeByType?.["1H"];
                const hasData1H = raw1H && raw1H.ppm?.length > 0;

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
                        {hasData1H && (
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
                          {oneH.metrics.ransy_score && (
                            <Pill label="Score" value={oneH.metrics.ransy_score}/>
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

                      <StickChart
                        key={[
                          refreshKeyByType["1H"] ?? 0,
                          hasData1H ? "loaded" : "",
                          JSON.stringify(zoom1H || null),
                        ].join("|")}
                        theory={showTheory ? oneH.theory : []}
                        expPpm={showExp ? oneH.expX : []}
                        expInt={showExp ? oneH.expI : []}
                        reverse={reverse1D}
                        showTheory={showTheory}
                        showExp={showExp}
                        height={PLOT_HEIGHT}
                        spectrumData={raw1H}
                        ppmRange={zoom1H}
                        spectrumColor="#000000"
                        canInteract={canInteract1H}
                        zoomEnabled={zoomRectEnabled}
                        intensityScale={intensityScaleByType["1H"] || 1.0}
                        onIntensityChange={(delta) => handleIntensityWheel("1H", delta)}
                        onZoom={(ppmMin, ppmMax) => {
                          setZoomRangeByType(prev => ({
                            ...prev,
                            "1H": { min: Math.min(ppmMin, ppmMax), max: Math.max(ppmMin, ppmMax) }
                          }));
                        }}
                        onReset={() => {
                          if (raw1H?.ppm?.length) {
                            // Use bgOptions bounds if defined, otherwise full range
                            const opts = bgOptions?.["1H"] || {};
                            const fullMin = raw1H.ppm.reduce((a, b) => a < b ? a : b, raw1H.ppm[0]);
                            const fullMax = raw1H.ppm.reduce((a, b) => a > b ? a : b, raw1H.ppm[0]);
                            const ppmMin = Number.isFinite(opts.ppm_min) ? opts.ppm_min : fullMin;
                            const ppmMax = Number.isFinite(opts.ppm_max) ? opts.ppm_max : fullMax;
                            setZoomRangeByType(prev => ({
                              ...prev,
                              "1H": { min: ppmMin, max: ppmMax }
                            }));
                          }
                        }}
                        busy={!!busyByType?.["1H"]}
                        extraButtons={
                          <button
                            onClick={() => setSpectrumDetailOpen(true)}
                            disabled={!taskId}
                            style={{
                              padding: "2px 6px",
                              fontSize: 11,
                              fontWeight: 500,
                              backgroundColor: taskId ? "#3b82f6" : "#e5e7eb",
                              color: taskId ? "#fff" : "#9ca3af",
                              border: "none",
                              borderRadius: 6,
                              cursor: taskId ? "pointer" : "not-allowed",
                            }}
                            title={taskId ? "View simulated spectrum with multiplicity data" : "No prediction task loaded"}
                          >
                            View Detail
                          </button>
                        }
                      />
                    </div>
                  ),
                });
              }

              if (c13) {
                const raw13C = rawDataByType?.["13C"];
                const zoom13C = zoomRangeByType?.["13C"];
                const hasData13C = raw13C && raw13C.ppm?.length > 0;

                blocks.push({
                  key: "b-13C",
                  el: (
                    <div
                      style={{
                        border: "1px solid #e5e7eb",
                        borderRadius: 6, padding: 4,
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
                        {hasData13C && (
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
                          {c13.metrics.ransy_score && (
                            <Pill label="Score" value={c13.metrics.ransy_score}/>
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
                          hasData13C ? "loaded" : "",
                          JSON.stringify(zoom13C || null),
                        ].join("|")}
                        theory={showTheory ? c13.theory : []}
                        expPpm={showExp ? c13.expX : []}
                        expInt={showExp ? c13.expI : []}
                        reverse={reverse1D}
                        showTheory={showTheory}
                        showExp={showExp}
                        height={PLOT_HEIGHT}
                        spectrumData={raw13C}
                        ppmRange={zoom13C}
                        spectrumColor="#000000"
                        canInteract={canInteract13C}
                        zoomEnabled={zoomRectEnabled}
                        intensityScale={intensityScaleByType["13C"] || 1.0}
                        onIntensityChange={(delta) => handleIntensityWheel("13C", delta)}
                        onZoom={(ppmMin, ppmMax) => {
                          setZoomRangeByType(prev => ({
                            ...prev,
                            "13C": { min: Math.min(ppmMin, ppmMax), max: Math.max(ppmMin, ppmMax) }
                          }));
                        }}
                        onReset={() => {
                          if (raw13C?.ppm?.length) {
                            // Use bgOptions bounds if defined, otherwise full range
                            const opts = bgOptions?.["13C"] || {};
                            const fullMin = raw13C.ppm.reduce((a, b) => a < b ? a : b, raw13C.ppm[0]);
                            const fullMax = raw13C.ppm.reduce((a, b) => a > b ? a : b, raw13C.ppm[0]);
                            const ppmMin = Number.isFinite(opts.ppm_min) ? opts.ppm_min : fullMin;
                            const ppmMax = Number.isFinite(opts.ppm_max) ? opts.ppm_max : fullMax;
                            setZoomRangeByType(prev => ({
                              ...prev,
                              "13C": { min: ppmMin, max: ppmMax }
                            }));
                          }
                        }}
                        busy={!!busyByType?.["13C"]}
                      />
                    </div>
                  ),
                });
              }

              if (cosy) {
                const rawCOSY = rawDataByType?.["COSY"];
                const zoomCOSY = zoomRangeByType?.["COSY"];
                const hasDataCOSY = rawCOSY && rawCOSY.contours?.length > 0;

                blocks.push({
                  key: "b-COSY",
                  el: (
                    <div
                      style={{
                        border: "1px solid #e5e7eb",
                        borderRadius: 6, padding: 4,
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
                        {hasDataCOSY && (
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
                          hasDataCOSY ? "loaded" : "",
                          JSON.stringify(zoomCOSY || null),
                        ].join("|")}
                        theoryPts={showTheory ? cosy.theory : []}
                        expPts={showExp ? cosy.expPts : []}
                        expInt={showExp ? cosy.expI : []}
                        showTheory={showTheory}
                        showExp={showExp}
                        height={PLOT_HEIGHT}
                        contourData={rawCOSY}
                        zoomRange={zoomCOSY}
                        contourColor="#000000"
                        canInteract={canInteractCOSY}
                        zoomEnabled={zoomRectEnabled}
                        onZoom2D={(zoomWindow) => {
                          setZoomRangeByType(prev => ({
                            ...prev,
                            COSY: {
                              f1: { min: zoomWindow.f1.min, max: zoomWindow.f1.max },
                              f2: { min: zoomWindow.f2.min, max: zoomWindow.f2.max }
                            }
                          }));
                        }}
                        onReset={() => {
                          if (rawCOSY) {
                            // Use bgOptions bounds if defined, otherwise full range
                            const opts = bgOptions?.["COSY"] || {};
                            const f1Min = Number.isFinite(opts.f1_ppm_min) ? opts.f1_ppm_min : rawCOSY.f1_range?.[0];
                            const f1Max = Number.isFinite(opts.f1_ppm_max) ? opts.f1_ppm_max : rawCOSY.f1_range?.[1];
                            const f2Min = Number.isFinite(opts.f2_ppm_min) ? opts.f2_ppm_min : rawCOSY.f2_range?.[0];
                            const f2Max = Number.isFinite(opts.f2_ppm_max) ? opts.f2_ppm_max : rawCOSY.f2_range?.[1];
                            setZoomRangeByType(prev => ({
                              ...prev,
                              COSY: {
                                f1: { min: f1Min, max: f1Max },
                                f2: { min: f2Min, max: f2Max }
                              }
                            }));
                          }
                        }}
                        busy={!!busyByType?.["COSY"]}
                      />
                    </div>
                  ),
                });
              }

              if (isMethodAvailable('TOCSY') && tocsy) {
                const rawTOCSY = rawDataByType?.["TOCSY"];
                const zoomTOCSY = zoomRangeByType?.["TOCSY"];
                const hasDataTOCSY = rawTOCSY && rawTOCSY.contours?.length > 0;

                blocks.push({
                  key: "b-TOCSY",
                  el: (
                    <div style={{ border: "1px solid #e5e7eb", borderRadius: 6, padding: 4, overflow: "hidden" }}>
                      <div style={{ fontWeight: 700, marginBottom: 4, display: "flex", justifyContent: "space-between", alignItems: "center", fontSize: 12 }}>
                        <span>TOCSY</span>
                        {tocsy.maxBonds && <span>Max bonds: {tocsy.maxBonds}</span>}
                        {hasDataTOCSY && <span style={{ fontSize: 9, color: "#16a34a", background: "#dcfce7", padding: "1px 4px", borderRadius: 3 }}>OK</span>}
                      </div>
                      {tocsy.metrics && Object.keys(tocsy.metrics).length > 0 && (
                        <div style={{ marginBottom: 6, minHeight: 20 }}>
                          {tocsy.metrics.coverage_ransy && <Pill label="RANSY" value={tocsy.metrics.coverage_ransy} />}
                          {tocsy.metrics.coverage_sub && <Pill label="Sub" value={tocsy.metrics.coverage_sub} />}
                        </div>
                      )}
                      {renderControls2D("TOCSY")}
                      <Scatter2D
                        key={[refreshKeyByType["TOCSY"] ?? 0, hasDataTOCSY ? "loaded" : "", JSON.stringify(zoomTOCSY || null)].join("|")}
                        theoryPts={showTheory ? tocsy.theory : []}
                        expPts={showExp ? tocsy.expPts : []}
                        expInt={showExp ? tocsy.expI : []}
                        showTheory={showTheory}
                        showExp={showExp}
                        height={PLOT_HEIGHT}
                        contourData={rawTOCSY}
                        zoomRange={zoomTOCSY}
                        contourColor={customColors?.TOCSY?.hex || "#00A878"}
                        canInteract={canInteractTOCSY}
                        zoomEnabled={zoomRectEnabled}
                        onZoom2D={(zoomWindow) => setZoomRangeByType(prev => ({
                          ...prev,
                          TOCSY: {
                            f1: { min: zoomWindow.f1.min, max: zoomWindow.f1.max },
                            f2: { min: zoomWindow.f2.min, max: zoomWindow.f2.max },
                          },
                        }))}
                        onReset={() => {
                          if (rawTOCSY) {
                            const options = bgOptions?.["TOCSY"] || {};
                            const f1Min = Number.isFinite(options.f1_ppm_min) ? options.f1_ppm_min : rawTOCSY.f1_range?.[0];
                            const f1Max = Number.isFinite(options.f1_ppm_max) ? options.f1_ppm_max : rawTOCSY.f1_range?.[1];
                            const f2Min = Number.isFinite(options.f2_ppm_min) ? options.f2_ppm_min : rawTOCSY.f2_range?.[0];
                            const f2Max = Number.isFinite(options.f2_ppm_max) ? options.f2_ppm_max : rawTOCSY.f2_range?.[1];
                            setZoomRangeByType(prev => ({ ...prev, TOCSY: { f1: { min: f1Min, max: f1Max }, f2: { min: f2Min, max: f2Max } } }));
                          }
                        }}
                        busy={!!busyByType?.["TOCSY"]}
                      />
                    </div>
                  ),
                });
              }

              if (hsqc) {
                const rawHSQC = rawDataByType?.["HSQC"];
                const zoomHSQC = zoomRangeByType?.["HSQC"];
                const hasDataHSQC = rawHSQC && rawHSQC.contours?.length > 0;

                blocks.push({
                  key: "b-HSQC",
                  el: (
                    <div
                      style={{
                        border: "1px solid #e5e7eb",
                        borderRadius: 6, padding: 4,
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
                        {hasDataHSQC && (
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
                          hasDataHSQC ? "loaded" : "",
                          JSON.stringify(zoomHSQC || null),
                        ].join("|")}
                        theoryPts={showTheory ? hsqc.theory : []}
                        expPts={showExp ? hsqc.expPts : []}
                        expInt={showExp ? hsqc.expI : []}
                        showTheory={showTheory}
                        showExp={showExp}
                        height={PLOT_HEIGHT}
                        contourData={rawHSQC}
                        zoomRange={zoomHSQC}
                        contourColor="#000000"
                        canInteract={canInteractHSQC}
                        zoomEnabled={zoomRectEnabled}
                        onZoom2D={(zoomWindow) => {
                          setZoomRangeByType(prev => ({
                            ...prev,
                            HSQC: {
                              f1: { min: zoomWindow.f1.min, max: zoomWindow.f1.max },
                              f2: { min: zoomWindow.f2.min, max: zoomWindow.f2.max }
                            }
                          }));
                        }}
                        onReset={() => {
                          if (rawHSQC) {
                            // Use bgOptions bounds if defined, otherwise full range
                            const opts = bgOptions?.["HSQC"] || {};
                            const f1Min = Number.isFinite(opts.f1_ppm_min) ? opts.f1_ppm_min : rawHSQC.f1_range?.[0];
                            const f1Max = Number.isFinite(opts.f1_ppm_max) ? opts.f1_ppm_max : rawHSQC.f1_range?.[1];
                            const f2Min = Number.isFinite(opts.f2_ppm_min) ? opts.f2_ppm_min : rawHSQC.f2_range?.[0];
                            const f2Max = Number.isFinite(opts.f2_ppm_max) ? opts.f2_ppm_max : rawHSQC.f2_range?.[1];
                            setZoomRangeByType(prev => ({
                              ...prev,
                              HSQC: {
                                f1: { min: f1Min, max: f1Max },
                                f2: { min: f2Min, max: f2Max }
                              }
                            }));
                          }
                        }}
                        busy={!!busyByType?.["HSQC"]}
                      />
                    </div>
                  ),
                });
              }

              if (displayedSpectrum !== "all") {
                blocks = blocks.filter(block => block.key === `b-${displayedSpectrum}`);
              }
              const n = blocks.length;
              if (n === 0) return <div>No spectrum available</div>;
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

      {/* Color picker */}
      <ColorPicker
        isOpen={colorPickerOpen}
        onClose={() => setColorPickerOpen(false)}
        customColors={customColors}
        onColorChange={setCustomColor}
      />

      {/* Molecule grid modal. */}
      <ShowAllMoleculesModal
        isOpen={showAllOpen}
        onClose={() => setShowAllOpen(false)}
        sources={sortedSourcesAll}
        sourcesWithScore={sortedSourcesWithScore}
        page={allPage}
        setPage={setAllPage}
        pageSize={allPageSize}
        setPageSize={setAllPageSize}
        getRowsByType={rowsByTypeForSource}
        csvUrl={pivotCsvUrl}
        taskId={taskId}
        customColors={customColors}
      />

      {/* Detailed simulated 1H spectrum modal */}
      <SpectrumDetailModal
        isOpen={spectrumDetailOpen}
        onClose={() => setSpectrumDetailOpen(false)}
        molId={selectedSource || oneH?.source || "mol"}
        taskId={taskId}
        theoryShifts={oneH?.theory || []}
        expShifts={oneH?.expX || []}
        pdataDir1H={pdataDirs?.["1H"] || null}
        matchTolerancePpm={0.05}
        // Reuse loaded data when opening the modal.
        multiplicityJson={preloadedMultiplicityData}
        preloadedSpectrumData={preloadedSpectrumData}
        isPreloading={preloadingMultiplicity || preloadingSpectrum}
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
