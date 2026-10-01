// backend/routes/spectrum.js
// Detailed 1H spectrum simulation routes.

const express = require('express');
const router = express.Router();
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const { getPythonExe } = require('../utils/spawn');
const { normalizeHostPathToContainer } = require('../utils/pathNormalizer');

/**
 * POST /api/simulate-spectrum
 * 
 * Data sources:
 * 
 * 1. JSON: use precomputed multiplicity data.
 *   - multiplicity_json: object (complete multiplicity data)
 *   - or multiplicity_json_path: string (path to the JSON file)
 * 
 * 2. SMILES fallback: calculate multiplicities on demand.
 *   - smiles: string (molecule SMILES)
 * 
 * Shared options:
 *   - mol_id: string (optional identifier)
 *   - freq_mhz: number (frequency in MHz, default 700)
 *   - solvent: string (solvent, default "CDCl3")
 *   - line_shape: string ("lorentzian" | "gaussian" | "voigt", default "lorentzian")
 *   - digital_resolution_hz: number (resolution in Hz, default 0.1)
 *   - ppm_range: [number, number] (ppm range, default [-0.5, 12])
 * 
 * Response:
 *   - ppm: number[] (ppm axis)
 *   - intensity: number[] (normalized intensity)
 *   - peaks: { ppm, intensity, label, multiplicity, J_hz, n_protons }[]
 *   - mol_id: string
 *   - freq_mhz: number
 *   - status: "ok" | "error"
 */
router.post('/api/simulate-spectrum', async (req, res) => {
  const {
    smiles,
    multiplicity_json,
    multiplicity_json_path,
    mol_id = 'mol',
    freq_mhz = 700,
    solvent = 'CDCl3',
    line_shape = 'lorentzian',
    digital_resolution_hz = 0.1,
    ppm_range = [-0.5, 12],
    line_width_hz, // Optional; the Python backend calculates this when omitted.
    // Dynamics postprocessing settings.
    enable_dynamics = false,
    temperature_k = 298,
    ph: ph_raw = 7.4,
  } = req.body;

  // Validate pH against the supported chemical range.
  let ph = Number(ph_raw) || 7.4;
  if (isNaN(ph) || ph < 0.0 || ph > 14.0) {
    console.warn(`[simulate-spectrum] pH invalide: ${ph_raw}, fallback vers 7.4`);
    ph = 7.4;
  }

  console.log('[simulate-spectrum] Request received:', {
    hasSmiles: !!smiles,
    hasMultiplicityJson: !!multiplicity_json,
    hasMultiplicityJsonPath: !!multiplicity_json_path,
    mol_id,
    freq_mhz,
    multiplicity_json_size: multiplicity_json ? JSON.stringify(multiplicity_json).length : 0
  });

  // Require at least one simulation data source.
  if (!smiles && !multiplicity_json && !multiplicity_json_path) {
    return res.status(400).json({ 
      error: 'Either smiles, multiplicity_json, or multiplicity_json_path is required' 
    });
  }

  // Resolve the Python script path.
  const scriptsDir = path.resolve(__dirname, '../scripts');
  const simulatorScript = path.join(scriptsDir, 'simulate_spectrum_api.py');
  
  // Check that the script exists.
  if (!fs.existsSync(simulatorScript)) {
    return res.status(500).json({ error: 'Simulator script not found' });
  }

  // Use the bundled Python runtime so the installed app does not require system Python.
  console.log('[simulate-spectrum] Environment:', {
    IS_ELECTRON: process.env.IS_ELECTRON,
    IS_ELECTRON_BUNDLED: process.env.IS_ELECTRON_BUNDLED,
    BUNDLED_PYTHON_PATH: process.env.BUNDLED_PYTHON_PATH,
    RESOURCES_PATH: process.env.RESOURCES_PATH,
  });

  let pythonPath;
  try {
    pythonPath = getPythonExe();
    console.log('[simulate-spectrum] Python path resolved to:', pythonPath);
    console.log('[simulate-spectrum] Python exists?', fs.existsSync(pythonPath));
  } catch (err) {
    console.error('[simulate-spectrum] Bundled Python not found:', err.message);
    return res.status(500).json({
      status: 'error',
      error: 'Bundled Python not found. The application requires the bundled Python to run.',
      details: err.message,
      env: {
        BUNDLED_PYTHON_PATH: process.env.BUNDLED_PYTHON_PATH,
        RESOURCES_PATH: process.env.RESOURCES_PATH,
      }
    });
  }

  // The bundled Python executable dispatches commands internally.
  const isBundledEnv = process.env.IS_ELECTRON_BUNDLED === 'true';
  const isActuallyBundled = isBundledEnv && pythonPath.endsWith('.exe') && fs.existsSync(pythonPath);
  
  let args;
  if (isActuallyBundled) {
    // Mode bundled: exe simulate_spectrum [args]
    args = ['simulate_spectrum'];
  } else {
    // Mode dev: python script.py [args] 
    args = [simulatorScript];
  }
  
  // Inline multiplicity JSON.
  if (multiplicity_json) {
    args.push('--json_data', JSON.stringify(multiplicity_json));
  }
  // Mode 2: Chemin vers fichier JSON
  else if (multiplicity_json_path && fs.existsSync(multiplicity_json_path)) {
    args.push('--json_path', multiplicity_json_path);
  }
  // Mode 3: SMILES (fallback)
  else if (smiles) {
    args.push('--smiles', smiles);
  }
  
  // Options communes
  args.push(
    '--mol_id', String(mol_id),
    '--freq_mhz', String(freq_mhz),
    '--solvent', String(solvent),
    '--line_shape', String(line_shape),
    '--digital_resolution_hz', String(digital_resolution_hz),
    '--ppm_min', String(ppm_range[0] ?? -0.5),
    '--ppm_max', String(ppm_range[1] ?? 12),
    // Dynamics postprocessing settings.
    '--temperature_k', String(temperature_k),
    '--ph', String(ph),
  );

  // Use the supplied line_width_hz, or let Python calculate it when omitted.
  if (line_width_hz !== null && line_width_hz !== undefined) {
    args.push('--line_width_hz', String(line_width_hz));
  }

  // Flag optionnel enable_dynamics
  if (enable_dynamics) {
    args.push('--enable_dynamics');
  }

  console.log('[simulate-spectrum] Executing Python command:', pythonPath, args.slice(0, 5).join(' '), '...');

  let stdout = '';
  let stderr = '';
  let responseSent = false;

  const sendResponse = (statusCode, data) => {
    if (responseSent) return;
    responseSent = true;
    res.status(statusCode).json(data);
  };

  // For bundled Python exe, use its directory as cwd so it can find _internal/
  // For dev Python, use scriptsDir
  const pythonDir = path.dirname(pythonPath);
  const workingDir = isActuallyBundled ? pythonDir : scriptsDir;

  console.log('[simulate-spectrum] Working directory:', workingDir);

  const child = spawn(pythonPath, args, {
    cwd: workingDir,
    env: { ...process.env, PYTHONIOENCODING: 'utf-8' },
  });

  child.stdout.on('data', (data) => {
    stdout += data.toString();
  });

  child.stderr.on('data', (data) => {
    stderr += data.toString();
  });

  child.on('close', (code) => {
    console.log('[simulate-spectrum] Python process exited with code:', code);
    console.log('[simulate-spectrum] stdout length:', stdout.length, 'stderr length:', stderr.length);

    if (code !== 0) {
      console.error('[simulate-spectrum] Python error:', stderr);
      console.error('[simulate-spectrum] Full stdout:', stdout);
      // Return detailed error to frontend for debugging
      sendResponse(500, {
        status: 'error',
        error: 'Simulation failed',
        details: stderr.substring(0, 2000),
        stdout: stdout.substring(0, 1000),
        exitCode: code,
        pythonPath,
        args: args.join(' ')
      });
      return;
    }

    try {
      // The script writes JSON to stdout.
      const result = JSON.parse(stdout);
      console.log('[simulate-spectrum] Simulation successful, result keys:', Object.keys(result));
      sendResponse(200, result);
    } catch (e) {
      console.error('[simulate-spectrum] JSON parse error:', e.message);
      console.error('[simulate-spectrum] Raw stdout (first 500 chars):', stdout.substring(0, 500));
      sendResponse(500, { 
        status: 'error',
        error: 'Failed to parse simulation result',
        details: stdout.substring(0, 500)
      });
    }
  });

  child.on('error', (err) => {
    console.error('[simulate-spectrum] Spawn error:', err);
    sendResponse(500, { error: 'Failed to start Python process', details: err.message });
  });
});

/**
 * GET /api/extract-bruker-spectrum
 *
 * Extract raw Bruker spectrum data (ppm, intensity arrays) from pdata folder.
 *
 * Query params:
 *   - pdata_dir: string (path to pdata/1 folder)
 *
 * Response:
 *   - ppm: number[] (ppm scale)
 *   - intensity: number[] (normalized intensities)
 *   - n_points: number
 *   - sf: number (spectrometer frequency MHz)
 *   - status: "ok" | "error"
 */
router.get('/api/extract-bruker-spectrum', async (req, res) => {
  const { pdata_dir } = req.query;

  if (!pdata_dir) {
    return res.status(400).json({
      status: 'error',
      error: 'pdata_dir query parameter is required'
    });
  }

  console.log('[extract-bruker] Request for pdata_dir:', pdata_dir);

  const scriptsDir = path.resolve(__dirname, '../scripts');
  const extractScript = path.join(scriptsDir, 'extract_bruker_spectrum.py');

  let pythonPath;
  try {
    pythonPath = getPythonExe();
    console.log('[extract-bruker] Python path resolved to:', pythonPath);
  } catch (err) {
    console.error('[extract-bruker] Python not found:', err);
    return res.status(500).json({
      status: 'error',
      error: 'Python not found',
      details: err.message
    });
  }

  // Detect bundled vs dev mode (same logic as simulate-spectrum)
  const isBundledEnv = process.env.IS_ELECTRON_BUNDLED === 'true';
  const isActuallyBundled = isBundledEnv && pythonPath.endsWith('.exe') && fs.existsSync(pythonPath);

  // Normalize host path to container path in Docker mode
  const resolvedPdataDir = isActuallyBundled ? pdata_dir : normalizeHostPathToContainer(pdata_dir);
  console.log('[extract-bruker] Resolved pdata_dir:', resolvedPdataDir);

  let args;
  if (isActuallyBundled) {
    // Mode bundled: exe extract_bruker_spectrum [args]
    console.log('[extract-bruker] Running in BUNDLED mode');
    args = ['extract_bruker_spectrum'];
  } else {
    // Mode dev: python script.py [args]
    console.log('[extract-bruker] Running in DEV mode');
    if (!fs.existsSync(extractScript)) {
      return res.status(500).json({
        status: 'error',
        error: 'Extract script not found in dev mode',
        path: extractScript
      });
    }
    args = [extractScript];
  }

  args.push(
    '--pdata_dir', resolvedPdataDir,
    '--output', 'json'
  );

  console.log('[extract-bruker] Executing:', pythonPath, args.join(' '));

  // Use same working directory logic as simulate-spectrum
  const pythonDir = path.dirname(pythonPath);
  const workingDir = isActuallyBundled ? pythonDir : scriptsDir;

  console.log('[extract-bruker] Working directory:', workingDir);

  const child = spawn(pythonPath, args, {
    cwd: workingDir,
    env: { ...process.env, PYTHONIOENCODING: 'utf-8' }
  });

  let stdout = '';
  let stderr = '';

  child.stdout.on('data', (data) => {
    stdout += data.toString();
  });

  child.stderr.on('data', (data) => {
    stderr += data.toString();
    console.error('[extract-bruker] stderr:', data.toString());
  });

  child.on('close', (code) => {
    console.log('[extract-bruker] Process exited with code:', code);

    if (code !== 0) {
      console.error('[extract-bruker] Error output:', stderr);
      return res.status(500).json({
        status: 'error',
        error: 'Extraction failed',
        stderr: stderr,
        exit_code: code
      });
    }

    try {
      const result = JSON.parse(stdout);
      res.json(result);
    } catch (e) {
      console.error('[extract-bruker] JSON parse error:', e.message);
      res.status(500).json({
        status: 'error',
        error: 'Failed to parse result',
        details: stdout.substring(0, 500)
      });
    }
  });

  child.on('error', (err) => {
    console.error('[extract-bruker] Spawn error:', err);
    res.status(500).json({
      status: 'error',
      error: 'Failed to start Python',
      details: err.message
    });
  });
});

/**
 * GET /api/extract-bruker-2d
 *
 * Extract 2D Bruker spectrum contours (COSY, HSQC) from pdata folder.
 * Returns polyline contours for direct canvas rendering.
 *
 * Query params:
 *   - pdata_dir: string (path to pdata/1 folder)
 *   - level_start: number (optional, default 0.08)
 *   - n_levels: number (optional, default 15)
 *   - decimate: number (optional, default 4)
 *
 * Response:
 *   - contours: [{level, paths: [[[x,y],...]]}]
 *   - f2_range: [min, max]
 *   - f1_range: [min, max]
 *   - status: "ok" | "error"
 */
router.get('/api/extract-bruker-2d', async (req, res) => {
  const { pdata_dir, level_start = 0.08, n_levels = 15, decimate = 4 } = req.query;

  if (!pdata_dir) {
    return res.status(400).json({
      status: 'error',
      error: 'pdata_dir query parameter is required'
    });
  }

  console.log('[extract-bruker-2d] Request for pdata_dir:', pdata_dir);

  const scriptsDir = path.resolve(__dirname, '../scripts');
  const extractScript = path.join(scriptsDir, 'extract_bruker_2d.py');

  let pythonPath;
  try {
    pythonPath = getPythonExe();
    console.log('[extract-bruker-2d] Python path resolved to:', pythonPath);
  } catch (err) {
    console.error('[extract-bruker-2d] Python not found:', err);
    return res.status(500).json({
      status: 'error',
      error: 'Python not found',
      details: err.message
    });
  }

  // Detect bundled vs dev mode
  const isBundledEnv = process.env.IS_ELECTRON_BUNDLED === 'true';
  const isActuallyBundled = isBundledEnv && pythonPath.endsWith('.exe') && fs.existsSync(pythonPath);

  // Normalize host path to container path in Docker mode
  const resolvedPdataDir = isActuallyBundled ? pdata_dir : normalizeHostPathToContainer(pdata_dir);
  console.log('[extract-bruker-2d] Resolved pdata_dir:', resolvedPdataDir);

  let args;
  if (isActuallyBundled) {
    console.log('[extract-bruker-2d] Running in BUNDLED mode');
    args = ['extract_bruker_2d'];
  } else {
    console.log('[extract-bruker-2d] Running in DEV mode');
    if (!fs.existsSync(extractScript)) {
      return res.status(500).json({
        status: 'error',
        error: 'Extract 2D script not found in dev mode',
        path: extractScript
      });
    }
    args = [extractScript];
  }

  args.push(
    '--pdata_dir', resolvedPdataDir,
    '--level_start', String(level_start),
    '--n_levels', String(n_levels),
    '--decimate', String(decimate),
    '--output', 'json'
  );

  console.log('[extract-bruker-2d] Executing:', pythonPath, args.join(' '));

  const pythonDir = path.dirname(pythonPath);
  const workingDir = isActuallyBundled ? pythonDir : scriptsDir;

  const child = spawn(pythonPath, args, {
    cwd: workingDir,
    env: { ...process.env, PYTHONIOENCODING: 'utf-8' }
  });

  let stdout = '';
  let stderr = '';

  child.stdout.on('data', (data) => {
    stdout += data.toString();
  });

  child.stderr.on('data', (data) => {
    stderr += data.toString();
    console.error('[extract-bruker-2d] stderr:', data.toString());
  });

  child.on('close', (code) => {
    console.log('[extract-bruker-2d] Process exited with code:', code);

    if (code !== 0) {
      console.error('[extract-bruker-2d] Error output:', stderr);
      return res.status(500).json({
        status: 'error',
        error: 'Extraction failed',
        stderr: stderr,
        exit_code: code
      });
    }

    try {
      const result = JSON.parse(stdout);
      res.json(result);
    } catch (e) {
      console.error('[extract-bruker-2d] JSON parse error:', e.message);
      res.status(500).json({
        status: 'error',
        error: 'Failed to parse result',
        details: stdout.substring(0, 500)
      });
    }
  });

  child.on('error', (err) => {
    console.error('[extract-bruker-2d] Spawn error:', err);
    res.status(500).json({
      status: 'error',
      error: 'Failed to start Python',
      details: err.message
    });
  });
});

/**
 * POST /api/simulate-spectrum-batch
 * 
 * Simulate spectra for all molecules in all_multiplicities.json.
 * 
 * Body:
 *   - task_id: string (task UUID used to locate the JSON file)
 *   - freq_mhz: number (frequency in MHz, default 700)
 *   - solvent: string (solvent, default "CDCl3")
 *   - line_shape: string ("lorentzian" | "gaussian" | "voigt")
 *   - line_width_hz: number (line width in Hz)
 *   - enable_dynamics: boolean
 *   - temperature_k: number
 *   - ph: number
 * 
 * Response:
 *   - status: "ok" | "error"
 *   - mode: "batch"
 *   - total_molecules: number
 *   - successful: number
 *   - failed: number
 *   - spectra: array of spectrum objects (one per molecule)
 */
router.post('/api/simulate-spectrum-batch', async (req, res) => {
  const {
    task_id,
    freq_mhz = 700,
    solvent = 'CDCl3',
    line_shape = 'lorentzian',
    digital_resolution_hz = 0.1,
    ppm_range = [-0.5, 12],
    line_width_hz, // Optional; the Python backend calculates this when omitted.
    enable_dynamics = false,
    temperature_k = 298,
    ph: ph_raw = 7.4,
  } = req.body;

  // Validate pH before passing it to the simulator. 
  let ph = ph_raw;
  console.log(`[simulate-spectrum-batch] pH reçu: ${ph_raw} (type: ${typeof ph_raw})`);
  
  if (ph_raw === 26) {
    console.warn(`[simulate-spectrum-batch] DANGER: pH=26 détecté - probablement confusion avec année 2026 ou température! Fallback vers 7.4`);
    ph = 7.4;
  } else if (ph_raw === 2026 || ph_raw === 2025 || ph_raw === 2024 || ph_raw === 2023 || ph_raw === 2022 || ph_raw === 2021 || ph_raw === 2020) {
    console.warn(`[simulate-spectrum-batch] DANGER: pH=${ph_raw} détecté - probablement confusion avec année! Fallback vers 7.4`);
    ph = 7.4;
  } else if (ph_raw === 21 || ph_raw === 20 || ph_raw === 25 || ph_raw === 15 || ph_raw === 16 || ph_raw === 18 || ph_raw === 19 || ph_raw === 22 || ph_raw === 23 || ph_raw === 24) {
    console.warn(`[simulate-spectrum-batch] DANGER: pH=${ph_raw} détecté - probablement confusion avec température! Fallback vers 7.4`);
    ph = 7.4;
  } else if (ph_raw >= 15 && ph_raw <= 30) {
    console.warn(`[simulate-spectrum-batch] DANGER: pH=${ph_raw} dans plage suspecte [15-30] - probablement confusion avec température/année! Fallback vers 7.4`);
    ph = 7.4;
  } else if (typeof ph_raw !== 'number' || ph_raw < 1.0 || ph_raw > 14.0) {
    console.warn(`[simulate-spectrum-batch] pH invalide reçu: ${ph_raw} (type: ${typeof ph_raw}), fallback vers 7.4`);
    ph = 7.4;  // Use the default when validation fails.
  }

  console.log('[simulate-spectrum-batch] Request received:', { task_id, freq_mhz });

  if (!task_id) {
    return res.status(400).json({ 
      status: 'error',
      error: 'task_id is required' 
    });
  }

  // Locate the task directory.
  const projectService = require('../services/projectService');
  const currentProject = projectService.getCurrent();
  const outputsDir = path.resolve(__dirname, '../../outputs');
  const baseDir = currentProject?.workDir
    ? currentProject.workDir
    : outputsDir;

  const taskDir = path.join(baseDir, task_id);
  const allMultiplicitiesPath = path.join(taskDir, 'all_multiplicities.json');

  console.log('[simulate-spectrum-batch] Looking for:', allMultiplicitiesPath);

  // Locate the multiplicity data file.
  let multiplicityFilePath = null;
  let isSingleMode = false;

  if (fs.existsSync(allMultiplicitiesPath)) {
    // Mode batch: all_multiplicities.json existe
    multiplicityFilePath = allMultiplicitiesPath;
    console.log('[simulate-spectrum-batch] Found all_multiplicities.json');
  } else {
    // In single mode, find a *_multiplicity.json file.
    console.log('[simulate-spectrum-batch] all_multiplicities.json not found, looking for single molecule file...');
    
    if (fs.existsSync(taskDir)) {
      const files = fs.readdirSync(taskDir);
      const singleMultFile = files.find(f => f.endsWith('_multiplicity.json'));
      
      if (singleMultFile) {
        multiplicityFilePath = path.join(taskDir, singleMultFile);
        isSingleMode = true;
        console.log('[simulate-spectrum-batch] Found single molecule file:', singleMultFile);
      }
    }
  }

  if (!multiplicityFilePath) {
    return res.status(404).json({
      status: 'error',
      error: 'No multiplicity JSON found (neither all_multiplicities.json nor *_multiplicity.json)',
      path: taskDir
    });
  }

  // Resolve the Python executable.
  const scriptsDir = path.resolve(__dirname, '../scripts');
  const simulatorScript = path.join(scriptsDir, 'simulate_spectrum_api.py');

  let pythonPath;
  try {
    pythonPath = getPythonExe();
  } catch (err) {
    return res.status(500).json({
      status: 'error',
      error: 'Bundled Python not found',
      details: err.message
    });
  }

  const isBundledEnv = process.env.IS_ELECTRON_BUNDLED === 'true';
  const isActuallyBundled = isBundledEnv && pythonPath.endsWith('.exe') && fs.existsSync(pythonPath);

  let args;
  if (isActuallyBundled) {
    args = ['simulate_spectrum'];
  } else {
    args = [simulatorScript];
  }

  // Build arguments for the selected execution mode.
  args.push(
    '--json_path', multiplicityFilePath,
    '--freq_mhz', String(freq_mhz),
    '--solvent', String(solvent),
    '--line_shape', String(line_shape),
    '--digital_resolution_hz', String(digital_resolution_hz),
    '--ppm_min', String(ppm_range[0] ?? -0.5),
    '--ppm_max', String(ppm_range[1] ?? 12),
    '--temperature_k', String(temperature_k),
    '--ph', String(ph),
  );

  // Use the supplied line_width_hz, or let Python calculate it when omitted.
  if (line_width_hz !== null && line_width_hz !== undefined) {
    args.push('--line_width_hz', String(line_width_hz));
  }

  // Mode batch seulement si c'est all_multiplicities.json
  if (!isSingleMode) {
    args.push('--all');
  }

  if (enable_dynamics) {
    args.push('--enable_dynamics');
  }

  console.log('[simulate-spectrum-batch] Executing simulation...', { isSingleMode });

  const pythonDir = path.dirname(pythonPath);
  const workingDir = isActuallyBundled ? pythonDir : scriptsDir;

  let stdout = '';
  let stderr = '';

  const child = spawn(pythonPath, args, {
    cwd: workingDir,
    env: { ...process.env, PYTHONIOENCODING: 'utf-8' },
  });

  child.stdout.on('data', (data) => {
    stdout += data.toString();
  });

  child.stderr.on('data', (data) => {
    stderr += data.toString();
    // Log stderr without failing the request; it can contain warnings.
    console.log('[simulate-spectrum-batch] stderr:', data.toString().substring(0, 200));
  });

  child.on('close', (code) => {
    console.log('[simulate-spectrum-batch] Process exited with code:', code);

    if (code !== 0) {
      return res.status(500).json({
        status: 'error',
        error: 'Simulation failed',
        stderr: stderr.substring(0, 1000),
        exit_code: code
      });
    }

    try {
      const result = JSON.parse(stdout);
      
      // Wrap a single-molecule result in the batch response format.
      if (isSingleMode) {
        console.log('[simulate-spectrum-batch] Single mode result for:', result.mol_id);
        const batchResult = {
          status: 'ok',
          mode: 'single',
          total_molecules: 1,
          successful: result.status === 'ok' ? 1 : 0,
          failed: result.status === 'ok' ? 0 : 1,
          spectra: [result]
        };
        res.json(batchResult);
      } else {
        // Mode batch normal
        console.log('[simulate-spectrum-batch] Batch mode result:', {
          total: result.total_molecules,
          successful: result.successful,
          failed: result.failed
        });
        res.json(result);
      }
    } catch (e) {
      console.error('[simulate-spectrum-batch] JSON parse error:', e.message);
      res.status(500).json({
        status: 'error',
        error: 'Failed to parse result',
        details: stdout.substring(0, 500)
      });
    }
  });

  child.on('error', (err) => {
    console.error('[simulate-spectrum-batch] Spawn error:', err);
    res.status(500).json({
      status: 'error',
      error: 'Failed to start Python',
      details: err.message
    });
  });
});

module.exports = router;
