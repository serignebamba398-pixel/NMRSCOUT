const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const { BG_SCRIPT } = require('../config/paths');
const { ensureDir } = require('../utils/filesystem');
const { getPythonExe, buildPythonEnv } = require('../utils/spawn');
const {
  DEFAULT_BG_OPTIONS_1H,
  DEFAULT_BG_OPTIONS_13C,
  DEFAULT_BG_OPTIONS_COSY,
  DEFAULT_BG_OPTIONS_HSQC,
  BG_OUTPUT_NAMES,
} = require('../constants/defaults');

async function runOneBgExport({ method, pdataDir, workDir, options = {} }) {
  if (!['1H', '13C', 'COSY', 'HSQC'].includes(method)) throw new Error('invalid method');
  if (!pdataDir) throw new Error('pdata_dir is required');
  
  const PY = getPythonExe();
  const env = buildPythonEnv();
  
  // Detect if we're using bundled executable
  const isBundledExe = PY.toLowerCase().endsWith('nmr_engine.exe') || 
                       PY.toLowerCase().endsWith('nmr_engine') ||
                       PY.toLowerCase().endsWith('nmr_engine_production.exe') ||
                       PY.toLowerCase().endsWith('nmr_engine_production');
  
  // Only check for script file if NOT using bundled exe
  if (!isBundledExe && !fs.existsSync(BG_SCRIPT)) {
    throw new Error('nmr_export_bg.py not found in /scripts');
  }

  const out_png = path.join(workDir, BG_OUTPUT_NAMES[method].png);
  const out_json = path.join(workDir, BG_OUTPUT_NAMES[method].json);

  // Build args differently for bundled exe vs script
  let args;
  if (isBundledExe) {
    // Bundled: nmr_engine_production.exe nmr_export_bg --method 1H ...
    args = [
      'nmr_export_bg',
      '--method',
      method,
      '--pdata_dir',
      pdataDir,
      '--out_png',
      out_png,
      '--out_json',
      out_json,
    ];
  } else {
    // Script mode: python nmr_export_bg.py --method 1H ...
    args = [
      BG_SCRIPT,
      '--method',
      method,
      '--pdata_dir',
      pdataDir,
      '--out_png',
      out_png,
      '--out_json',
      out_json,
    ];
  }

  const defaults = {
    oneH: DEFAULT_BG_OPTIONS_1H,
    c13: DEFAULT_BG_OPTIONS_13C,
    cosy: DEFAULT_BG_OPTIONS_COSY,
    hsqc: DEFAULT_BG_OPTIONS_HSQC,
  };

  const bool01 = v =>
    typeof v === 'boolean' ? (v ? 1 : 0) : Number(v) === 0 || Number(v) === 1 ? Number(v) : null;
  const asNum = v => (Number.isFinite(+v) ? +v : null);
  const asStr = v => (typeof v === 'string' && v.trim() ? v.trim() : null);

  if (method === '1H') {
    const o = { ...defaults.oneH, ...(options || {}) };
    args.push(
      '--ppm_min',
      String(o.ppm_min),
      '--ppm_max',
      String(o.ppm_max),
      '--width_px',
      String(o.width_px),
      '--height_px',
      String(o.height_px),
      '--dpi',
      String(o.dpi)
    );
    const t = bool01(o.transparent);
    if (t !== null) args.push('--transparent', String(t));
    const lc = asStr(o.line_color);
    if (lc) args.push('--line_color', lc);
    const lw = asNum(o.line_width);
    if (lw) args.push('--line_width', String(lw));
  } else if (method === '13C') {
    const o = { ...defaults.c13, ...(options || {}) };
    args.push(
      '--ppm_min',
      String(o.ppm_min),
      '--ppm_max',
      String(o.ppm_max),
      '--width_px',
      String(o.width_px),
      '--height_px',
      String(o.height_px),
      '--dpi',
      String(o.dpi)
    );
    const t = bool01(o.transparent);
    if (t !== null) args.push('--transparent', String(t));
    const lc = asStr(o.line_color);
    if (lc) args.push('--line_color', lc);
    const lw = asNum(o.line_width);
    if (lw) args.push('--line_width', String(lw));
  } else if (method === 'COSY') {
    const o = { ...defaults.cosy, ...(options || {}) };
    args.push(
      '--f2_ppm_min',
      String(o.f2_ppm_min),
      '--f2_ppm_max',
      String(o.f2_ppm_max),
      '--f1_ppm_min',
      String(o.f1_ppm_min),
      '--f1_ppm_max',
      String(o.f1_ppm_max),
      '--width_px',
      String(o.width_px),
      '--height_px',
      String(o.height_px),
      '--dpi',
      String(o.dpi),
      '--mode',
      asStr(o.mode) || 'contour'
    );
    const t = bool01(o.transparent);
    if (t !== null) args.push('--transparent', String(t));
    args.push('--n_levels', String(o.n_levels), '--level_start', String(o.level_start));
    const cmap = asStr(o.cmap);
    if (cmap) args.push('--cmap', cmap);
  } else if (method === 'HSQC') {
    const o = { ...defaults.hsqc, ...(options || {}) };
    args.push(
      '--f2_ppm_min',
      String(o.f2_ppm_min),
      '--f2_ppm_max',
      String(o.f2_ppm_max),
      '--f1_ppm_min',
      String(o.f1_ppm_min),
      '--f1_ppm_max',
      String(o.f1_ppm_max),
      '--width_px',
      String(o.width_px),
      '--height_px',
      String(o.height_px),
      '--dpi',
      String(o.dpi),
      '--mode',
      asStr(o.mode) || 'contour'
    );
    const t = bool01(o.transparent);
    if (t !== null) args.push('--transparent', String(t));
    args.push('--n_levels', String(o.n_levels), '--level_start', String(o.level_start));
    const cmap = asStr(o.cmap);
    if (cmap) args.push('--cmap', cmap);
  }

  const r = await new Promise(resolve => {
    const child = spawn(PY, args, {
      cwd: workDir,
      stdio: ['ignore', 'pipe', 'pipe'],
      shell: false,
      env,
    });
    let err = '';
    child.stderr.on('data', d => (err += d.toString()));
    child.on('close', code => resolve({ code, err }));
    child.on('error', e => resolve({ code: -1, err: e.message }));
  });

  if (r.code !== 0) throw new Error(r.err || `python exit ${r.code}`);

  const jsonOk = fs.existsSync(out_json) && fs.statSync(out_json).size > 0;
  const pngOk = fs.existsSync(out_png) && fs.statSync(out_png).size > 0;
  if (!jsonOk || !pngOk) throw new Error('PNG/JSON not produced');

  const meta = JSON.parse(fs.readFileSync(out_json, 'utf8'));
  return {
    ok: true,
    pngUrl: `/outputs/${path.basename(workDir)}/${path.basename(out_png)}`,
    jsonUrl: `/outputs/${path.basename(workDir)}/${path.basename(out_json)}`,
    meta,
  };
}

function computeWindow1DFromRect(meta, rectPx) {
  const bbox = meta?.plot_bbox_px;
  if (
    !bbox ||
    !Number.isFinite(bbox.x) ||
    !Number.isFinite(bbox.y) ||
    !Number.isFinite(bbox.w) ||
    !Number.isFinite(bbox.h)
  ) {
    throw new Error('meta.plot_bbox_px missing or invalid for 1D zoom');
  }

  let left, right;

  if (Number.isFinite(Number(meta.ppm_left)) && Number.isFinite(Number(meta.ppm_right))) {
    left = Number(meta.ppm_left);
    right = Number(meta.ppm_right);
  } else if (Number.isFinite(Number(meta.ppm_min)) && Number.isFinite(Number(meta.ppm_max))) {
    const a = Number(meta.ppm_min);
    const b = Number(meta.ppm_max);
    left = Math.max(a, b);
    right = Math.min(a, b);
  } else if (Number.isFinite(Number(meta.ppm_max)) && Number.isFinite(Number(meta.ppm_min))) {
    left = Number(meta.ppm_max);
    right = Number(meta.ppm_min);
  } else {
    throw new Error('No valid ppm information in meta for 1D zoom');
  }

  if (!Number.isFinite(left) || !Number.isFinite(right) || left <= right) {
    throw new Error(`Invalid ppm values: left=${left}, right=${right}`);
  }

  const x0 = Math.max(bbox.x, Math.min(bbox.x + bbox.w, rectPx.x));
  const x1 = Math.max(bbox.x, Math.min(bbox.x + bbox.w, rectPx.x + rectPx.w));
  const t0 = (x0 - bbox.x) / bbox.w;
  const t1 = (x1 - bbox.x) / bbox.w;
  const p0 = left + t0 * (right - left);
  const p1 = left + t1 * (right - left);

  const ppm_min = Math.min(p0, p1);
  const ppm_max = Math.max(p0, p1);

  if (!Number.isFinite(ppm_min) || !Number.isFinite(ppm_max) || Math.abs(ppm_max - ppm_min) < 1e-9) {
    const mid = (ppm_min + ppm_max) / 2;
    const pad = Math.max(0.001, (left - right) * 0.001);
    return {
      ppm_min: mid - pad,
      ppm_max: mid + pad,
      width_px_native: Math.abs(x1 - x0),
    };
  }

  return {
    ppm_min,
    ppm_max,
    width_px_native: Math.abs(x1 - x0),
  };
}

function computeWindow2DFromRect(meta, rectPx) {
  const bbox = meta?.plot_bbox_px;
  if (
    !bbox ||
    !Number.isFinite(bbox.x) ||
    !Number.isFinite(bbox.y) ||
    !Number.isFinite(bbox.w) ||
    !Number.isFinite(bbox.h)
  ) {
    throw new Error('meta.plot_bbox_px missing or invalid for 2D zoom');
  }

  const f2L = Number(meta.f2_left);
  const f2R = Number(meta.f2_right);
  const f1T = Number(meta.f1_top);
  const f1B = Number(meta.f1_bottom);

  if (![f2L, f2R, f1T, f1B].every(Number.isFinite)) {
    throw new Error('meta f1/f2 incomplete for 2D zoom');
  }

  const x0 = Math.max(bbox.x, Math.min(bbox.x + bbox.w, rectPx.x));
  const x1 = Math.max(bbox.x, Math.min(bbox.x + bbox.w, rectPx.x + rectPx.w));
  const y0 = Math.max(bbox.y, Math.min(bbox.y + bbox.h, rectPx.y));
  const y1 = Math.max(bbox.y, Math.min(bbox.y + bbox.h, rectPx.y + rectPx.h));

  const tx0 = (x0 - bbox.x) / bbox.w;
  const tx1 = (x1 - bbox.x) / bbox.w;
  const ty0 = (y0 - bbox.y) / bbox.h;
  const ty1 = (y1 - bbox.y) / bbox.h;

  const f2a = f2L + tx0 * (f2R - f2L);
  const f2b = f2L + tx1 * (f2R - f2L);
  const f1a = f1T + ty0 * (f1B - f1T);
  const f1b = f1T + ty1 * (f1B - f1T);

  const f2_ppm_min = Math.min(f2a, f2b);
  const f2_ppm_max = Math.max(f2a, f2b);
  const f1_ppm_min = Math.min(f1a, f1b);
  const f1_ppm_max = Math.max(f1a, f1b);

  if ([f2_ppm_min, f2_ppm_max, f1_ppm_min, f1_ppm_max].some(v => !Number.isFinite(v))) {
    throw new Error('2D zoom: non-finite ppm values');
  }
  if (f2_ppm_max - f2_ppm_min <= 0 || f1_ppm_max - f1_ppm_min <= 0) {
    throw new Error('2D zoom: empty window');
  }

  return { f2_ppm_min, f2_ppm_max, f1_ppm_min, f1_ppm_max };
}

module.exports = {
  runOneBgExport,
  computeWindow1DFromRect,
  computeWindow2DFromRect,
};