const express = require('express');
const router = express.Router();
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const { v4: uuidv4 } = require('uuid');
const { outputsDir, getActiveOutputsDir } = require('../config/paths');
const { ensureDir } = require('../utils/filesystem');
const { normalizeHostPathToContainer } = require('../utils/pathNormalizer');
const { logWithMemory } = require('../utils/logger');
const { runOneBgExport, computeWindow1DFromRect, computeWindow2DFromRect } = require('../services/backgroundService');
const {
  DEFAULT_BG_OPTIONS_1H,
  DEFAULT_BG_OPTIONS_13C,
  DEFAULT_BG_OPTIONS_COSY,
  DEFAULT_BG_OPTIONS_HSQC,
  BG_OUTPUT_NAMES,
} = require('../constants/defaults');
const { getPythonExe, buildPythonEnv } = require('../utils/spawn');

// Export batch backgrounds
router.post('/api/spectral-bg/export', async (req, res) => {
  try {
    const { BG_SCRIPT } = require('../config/paths');
    const PY = getPythonExe();
    
    // Detect if we're using bundled executable
    const isBundledExe = PY.toLowerCase().endsWith('nmr_engine.exe') || 
                         PY.toLowerCase().endsWith('nmr_engine') ||
                         PY.toLowerCase().endsWith('nmr_engine_production.exe') ||
                         PY.toLowerCase().endsWith('nmr_engine_production');
    
    // Only check for script file if NOT using bundled exe
    if (!isBundledExe && !fs.existsSync(BG_SCRIPT)) {
      return res.status(400).json({ error: 'nmr_export_bg.py not found in /scripts' });
    }

    const pdataIn = req.body?.pdata || {};
    if (!pdataIn || typeof pdataIn !== 'object' || Object.keys(pdataIn).length === 0) {
      return res
        .status(400)
        .json({ error: 'Body.pdata is required (paths to pdata folders)' });
    }
    
    // In Electron bundled mode, don't normalize paths (they're already native Windows paths)
    // Only normalize for Docker mode - use isBundledExe as it's more reliable
    const pdata = Object.fromEntries(
      Object.entries(pdataIn).map(([k, v]) => [k, isBundledExe ? String(v || '') : normalizeHostPathToContainer(String(v || ''))])
    );
    if (!isBundledExe) {
      for (const [k, raw] of Object.entries(pdataIn)) {
        const norm = normalizeHostPathToContainer(String(raw || ''));
        if (raw && norm !== raw) logWithMemory(`[bg/export] ${k}: ${raw} -> ${norm}`);
      }
    } else {
      logWithMemory('[bg/export] Bundled mode, paths not modified');
    }

    const methods =
      Array.isArray(req.body?.methods) && req.body.methods.length
        ? req.body.methods
        : Object.keys(pdata).filter(k => ['1H', '13C', 'COSY', 'TOCSY', 'HSQC'].includes(k));

    if (!methods.length) {
      return res.status(400).json({ error: 'No methods to process' });
    }

    const taskId = req.body?.persistKey || req.body?.taskId || uuidv4();
    const workDir = path.join(getActiveOutputsDir(), taskId);
    ensureDir(workDir);

    // PY already declared at top of function
    const env = buildPythonEnv();

    const opts = req.body?.options || {};
    const defaults = {
      oneH: DEFAULT_BG_OPTIONS_1H,
      c13: DEFAULT_BG_OPTIONS_13C,
      cosy: DEFAULT_BG_OPTIONS_COSY,
      hsqc: DEFAULT_BG_OPTIONS_HSQC,
    };

    const validateTransparent = value => {
      if (typeof value === 'boolean') {
        return value ? 1 : 0;
      }
      if (typeof value === 'number' && (value === 0 || value === 1)) {
        return value;
      }
      return null;
    };

    const runOne = method =>
      new Promise(resolve => {
        const pdataDir = pdata[method];
        if (!pdataDir) return resolve({ method, ok: false, error: 'Missing pdata' });
        if (!BG_OUTPUT_NAMES[method]) return resolve({ method, ok: false, error: 'Unsupported spectrum method' });
        const requiredFiles = method === '1H' || method === '13C'
          ? ['1r', 'procs'] : ['2rr', 'procs', 'proc2s'];
        const missing = requiredFiles.filter(name => !fs.existsSync(path.join(pdataDir, name)));
        if (missing.length) {
          return resolve({ method, ok: false, error:
            `Incomplete Bruker pdata folder: missing ${missing.join(', ')} in ${pdataDir}. ` +
            'A processed Bruker spectrum is required. NMRPipe .ft2 files are not supported by this reader.' });
        }

        const out_png = path.join(workDir, BG_OUTPUT_NAMES[method].png);
        const out_json = path.join(workDir, BG_OUTPUT_NAMES[method].json);

        // Build args differently for bundled exe vs script
        let baseArgs;
        if (isBundledExe) {
          // Bundled: nmr_engine_production.exe nmr_export_bg --method 1H ...
          baseArgs = [
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
          baseArgs = [
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

        if (method === '1H') {
          const o = { ...defaults.oneH, ...(opts.oneH || {}) };
          baseArgs.push(
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

          const transparentValue = validateTransparent(o.transparent);
          if (transparentValue !== null) {
            baseArgs.push('--transparent', String(transparentValue));
          }

          baseArgs.push('--line_color', String(o.line_color), '--line_width', String(o.line_width));
        } else if (method === '13C') {
          const o = { ...defaults.c13, ...(opts.c13 || {}) };
          baseArgs.push(
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

          const transparentValue = validateTransparent(o.transparent);
          if (transparentValue !== null) {
            baseArgs.push('--transparent', String(transparentValue));
          }

          baseArgs.push('--line_color', String(o.line_color), '--line_width', String(o.line_width));
        } else if (method === 'COSY' || method === 'TOCSY') {
          const o = { ...defaults.cosy, ...(opts[method === 'TOCSY' ? 'tocsy' : 'cosy'] || {}) };
          baseArgs.push(
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
            o.mode
          );

          const transparentValue = validateTransparent(o.transparent);
          if (transparentValue !== null) {
            baseArgs.push('--transparent', String(transparentValue));
          }

          baseArgs.push(
            '--n_levels',
            String(o.n_levels),
            '--level_start',
            String(o.level_start),
            '--cmap',
            String(o.cmap)
          );
        } else if (method === 'HSQC') {
          const o = { ...defaults.hsqc, ...(opts.hsqc || {}) };
          baseArgs.push(
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
            o.mode
          );

          const transparentValue = validateTransparent(o.transparent);
          if (transparentValue !== null) {
            baseArgs.push('--transparent', String(transparentValue));
          }

          baseArgs.push(
            '--n_levels',
            String(o.n_levels),
            '--level_start',
            String(o.level_start),
            '--cmap',
            String(o.cmap)
          );
        }

        const child = spawn(PY, baseArgs, {
          cwd: workDir,
          stdio: ['ignore', 'pipe', 'pipe'],
          shell: false,
          env,
        });
        let out = '',
          err = '';
        child.stdout.on('data', d => (out += d.toString()));
        child.stderr.on('data', d => (err += d.toString()));
        child.on('close', code => {
          if (code !== 0) {
            return resolve({ method, ok: false, error: err || `python exit ${code}` });
          }
          try {
            const metaPath = out_json;
            const pngOk = fs.existsSync(out_png) && fs.statSync(out_png).size > 0;
            const jsonOk = fs.existsSync(metaPath) && fs.statSync(metaPath).size > 0;
            let meta = null;
            if (jsonOk) meta = JSON.parse(fs.readFileSync(metaPath, 'utf8'));
            resolve({
              method,
              ok: pngOk && jsonOk,
              pngUrl: `/outputs/${path.basename(workDir)}/${BG_OUTPUT_NAMES[method].png}`,
              jsonUrl: `/outputs/${path.basename(workDir)}/${BG_OUTPUT_NAMES[method].json}`,
              meta,
              stderr: err || undefined,
            });
          } catch (e) {
            resolve({
              method,
              ok: false,
              error: e.message || 'parse meta failed',
              stderr: err || undefined,
            });
          }
        });
        child.on('error', e => resolve({ method, ok: false, error: e.message }));
      });

    const results = {};
    for (const m of methods) {
      results[m] = await runOne(m);
    }

    res.json({ taskId, workDir, outputs: results });
  } catch (e) {
    res.status(500).json({ error: e.message });
  }
});

// Zoom background
router.post('/api/spectral-bg/zoom', async (req, res) => {
  try {
    const { BG_SCRIPT } = require('../config/paths');
    const PY = getPythonExe();
    const isBundledExe = PY.toLowerCase().endsWith('nmr_engine.exe') || 
                         PY.toLowerCase().endsWith('nmr_engine') ||
                         PY.toLowerCase().endsWith('nmr_engine_production.exe') ||
                         PY.toLowerCase().endsWith('nmr_engine_production');
    
    const method = String(req.body?.method || '');
    const rawPdataDir = String(req.body?.pdata_dir || '');
    const pdataDir = isBundledExe ? rawPdataDir : normalizeHostPathToContainer(rawPdataDir);
    const meta = req.body?.meta;
    const select = req.body?.select;
    const display = req.body?.display_size || {};

    console.log('[ZOOM] Request received:', { method, rawPdataDir, pdataDir, isBundledExe });
    console.log('[ZOOM] Meta:', JSON.stringify(meta?.fig_px), JSON.stringify(meta?.plot_bbox_px));
    console.log('[ZOOM] Select:', JSON.stringify(select));

    if (!['1H', '13C', 'COSY', 'HSQC'].includes(method)) {
      console.log('[ZOOM] Error: invalid method');
      return res.status(400).json({ error: 'Invalid method' });
    }
    if (!pdataDir || !meta || !select) {
      console.log('[ZOOM] Error: missing params', { pdataDir: !!pdataDir, meta: !!meta, select: !!select });
      return res.status(400).json({ error: 'pdata_dir, meta and select are required' });
    }

    const taskId = req.body?.persistKey || req.body?.taskId || uuidv4();
    const workDir = path.join(getActiveOutputsDir(), taskId);
    ensureDir(workDir);

    const figW = Number(meta?.fig_px?.w) || 0;
    const figH = Number(meta?.fig_px?.h) || 0;
    const dispW = Number(display?.width_px) || figW;
    const dispH = Number(display?.height_px) || figH;
    const sx = figW && dispW ? figW / dispW : 1;
    const sy = figH && dispH ? figH / dispH : 1;
    const nativeRect = {
      x: +select.x * sx,
      y: +select.y * sy,
      w: +select.w * sx,
      h: +select.h * sy,
    };

    const size = req.body?.size || {};
    const line = req.body?.line || {};
    const contour = req.body?.contour || {};

    let options = {};
    if (method === '1H' || method === '13C') {
      try {
        const win1d = computeWindow1DFromRect(meta, nativeRect);

        options = {
          ppm_min: win1d.ppm_min,
          ppm_max: win1d.ppm_max,
          width_px: size.width_px,
          height_px: size.height_px,
          dpi: size.dpi,
          line_color: line.color,
          line_width: line.width,
        };
      } catch (error) {
        return res.status(500).json({ ok: false, error: error.message });
      }
    } else {
      const win2d = computeWindow2DFromRect(meta, nativeRect);
      options = {
        ...win2d,
        width_px: size.width_px,
        height_px: size.height_px,
        dpi: size.dpi,
        mode: contour.mode || 'contour',
        transparent: contour.transparent,
        n_levels: contour.n_levels,
        level_start: contour.level_start || 0.005,
        cmap: contour.cmap,
      };
    }

    const result = await runOneBgExport({ method, pdataDir, workDir, options });
    return res.json({ ok: true, taskId, ...result });
  } catch (e) {
    return res.status(500).json({ ok: false, error: e.message });
  }
});

// Relevel background (2D)
router.post('/api/spectral-bg/relevel', async (req, res) => {
  try {
    const { BG_SCRIPT } = require('../config/paths');
    const PY = getPythonExe();
    const isBundledExe = PY.toLowerCase().endsWith('nmr_engine.exe') || 
                         PY.toLowerCase().endsWith('nmr_engine') ||
                         PY.toLowerCase().endsWith('nmr_engine_production.exe') ||
                         PY.toLowerCase().endsWith('nmr_engine_production');
    
    const method = String(req.body?.method || '');
    const rawPdataDir = String(req.body?.pdata_dir || '');
    const pdataDir = isBundledExe ? rawPdataDir : normalizeHostPathToContainer(rawPdataDir);
    const meta = req.body?.meta;
    if (!['COSY', 'HSQC'].includes(method)) {
      return res.status(400).json({ error: 'Method must be COSY or HSQC' });
    }
    if (!pdataDir || !meta) return res.status(400).json({ error: 'pdata_dir and meta are required' });

    const taskId = req.body?.persistKey || req.body?.taskId || uuidv4();
    const workDir = path.join(getActiveOutputsDir(), taskId);
    ensureDir(workDir);

    const _num = v => {
      const n = Number(v);
      return Number.isFinite(n) ? n : null;
    };

    const cleanOptions = {};

    const f2_min = _num(meta.f2_left);
    const f2_max = _num(meta.f2_right);
    const f1_min = _num(meta.f1_top);
    const f1_max = _num(meta.f1_bottom);

    if (f2_min !== null) cleanOptions.f2_ppm_min = f2_min;
    if (f2_max !== null) cleanOptions.f2_ppm_max = f2_max;
    if (f1_min !== null) cleanOptions.f1_ppm_min = f1_min;
    if (f1_max !== null) cleanOptions.f1_ppm_max = f1_max;

    const w = _num(req.body?.size?.width_px);
    const h = _num(req.body?.size?.height_px);
    const d = _num(req.body?.size?.dpi);

    if (w !== null && w > 0) cleanOptions.width_px = w;
    if (h !== null && h > 0) cleanOptions.height_px = h;
    if (d !== null && d > 0) cleanOptions.dpi = d;

    cleanOptions.mode = req.body?.mode || req.body?.contour?.mode || 'contour';

    const trans = _num(req.body?.transparent) ?? _num(req.body?.contour?.transparent);
    if (trans !== null) cleanOptions.transparent = trans;

    const nlev = _num(req.body?.n_levels) ?? _num(req.body?.contour?.n_levels);
    if (nlev !== null && nlev > 0) cleanOptions.n_levels = nlev;

    const lvlStart = _num(req.body?.level_start) ?? _num(req.body?.contour?.level_start);
    if (lvlStart !== null && lvlStart > 0) cleanOptions.level_start = lvlStart;

    const cmap = req.body?.cmap || req.body?.contour?.cmap;
    if (cmap) cleanOptions.cmap = String(cmap);

    const result = await runOneBgExport({ method, pdataDir, workDir, options: cleanOptions });
    res.json({ ok: true, taskId, ...result });
  } catch (e) {
    res.status(500).json({ ok: false, error: e.message });
  }
});

// Interactive background
async function handleBgInteractive(req, res, method) {
  try {
    const { BG_SCRIPT } = require('../config/paths');
    const PY = getPythonExe();
    const isBundledExe = PY.toLowerCase().endsWith('nmr_engine.exe') || 
                         PY.toLowerCase().endsWith('nmr_engine') ||
                         PY.toLowerCase().endsWith('nmr_engine_production.exe') ||
                         PY.toLowerCase().endsWith('nmr_engine_production');
    
    const isGet = method === 'GET';
    const raw = isGet ? req.query : req.body;

    const methodStr = String(raw?.method || '').toUpperCase();
    const rawPdataDir = String(raw?.pdata_dir || '');
    const pdataDir = isBundledExe ? rawPdataDir : normalizeHostPathToContainer(rawPdataDir);

    if (!['1H', '13C', 'COSY', 'HSQC'].includes(methodStr)) {
      return res.status(400).json({ ok: false, error: 'Invalid method' });
    }
    if (!pdataDir) {
      return res.status(400).json({ ok: false, error: 'pdata_dir is required' });
    }

    const taskId = raw?.persistKey || raw?.taskId || uuidv4();
    const workDir = path.join(getActiveOutputsDir(), taskId);
    ensureDir(workDir);

    const options = {};
    const size = raw?.size || {};
    if (Number.isFinite(+size.width_px)) options.width_px = +size.width_px;
    if (Number.isFinite(+size.height_px)) options.height_px = +size.height_px;
    if (Number.isFinite(+size.dpi)) options.dpi = +size.dpi;

    if (methodStr === '1H' || methodStr === '13C') {
      const window = raw?.window || {};
      if (Number.isFinite(+window.ppm_min)) options.ppm_min = +window.ppm_min;
      if (Number.isFinite(+window.ppm_max)) options.ppm_max = +window.ppm_max;

      const line = raw?.line || {};
      if (line.color) options.line_color = String(line.color);
      if (Number.isFinite(+line.width)) options.line_width = +line.width;
    } else {
      const w2 = raw?.window2d || {};
      if (Number.isFinite(+w2.f2_ppm_min)) options.f2_ppm_min = +w2.f2_ppm_min;
      if (Number.isFinite(+w2.f2_ppm_max)) options.f2_ppm_max = +w2.f2_ppm_max;
      if (Number.isFinite(+w2.f1_ppm_min)) options.f1_ppm_min = +w2.f1_ppm_min;
      if (Number.isFinite(+w2.f1_ppm_max)) options.f1_ppm_max = +w2.f1_ppm_max;

      const c = raw?.contour || {};
      if (c.mode) options.mode = String(c.mode);
      if (c.cmap) options.cmap = String(c.cmap);
      if (Number.isFinite(+c.n_levels)) options.n_levels = +c.n_levels;
      if (Number.isFinite(+c.level_start)) options.level_start = +c.level_start;
      if (c.transparent !== undefined) {
        const t =
          typeof c.transparent === 'boolean'
            ? c.transparent
              ? 1
              : 0
            : Number(c.transparent) === 0 || Number(c.transparent) === 1
            ? Number(c.transparent)
            : undefined;
        if (t !== undefined) options.transparent = t;
      }
    }

    const result = await runOneBgExport({
      method: methodStr,
      pdataDir,
      workDir,
      options,
    });

    return res.json({
      ok: true,
      taskId,
      pngUrl: `/outputs/${path.basename(workDir)}/${path.basename(result.pngUrl)}`,
      jsonUrl: `/outputs/${path.basename(workDir)}/${path.basename(result.jsonUrl)}`,
      meta: result.meta,
    });
  } catch (e) {
    return res.status(500).json({ ok: false, error: e.message });
  }
}

router.get('/api/spectral-bg/interactive', (req, res) => handleBgInteractive(req, res, 'GET'));
router.post('/api/spectral-bg/interactive', (req, res) => handleBgInteractive(req, res, 'POST'));

router.get('/api/spectral-bg/interact', (req, res) => handleBgInteractive(req, res, 'GET'));
router.post('/api/spectral-bg/interact', (req, res) => handleBgInteractive(req, res, 'POST'));

module.exports = router;
