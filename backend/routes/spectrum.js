// backend/routes/spectrum.js
// Route pour la simulation de spectre 1H détaillé

const express = require('express');
const router = express.Router();
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const { getPythonExe } = require('../utils/spawn');

/**
 * POST /api/simulate-spectrum
 * 
 * DEUX MODES:
 * 
 * 1. Mode JSON (RECOMMANDÉ - utilise les données de multiplicité pré-calculées):
 *   - multiplicity_json: object (données de multiplicité complètes)
 *   - OU multiplicity_json_path: string (chemin vers le fichier JSON)
 * 
 * 2. Mode SMILES (fallback - calcule les multiplicités à la volée):
 *   - smiles: string (SMILES de la molécule)
 * 
 * Options communes:
 *   - mol_id: string (optionnel, identifiant)
 *   - freq_mhz: number (fréquence en MHz, défaut 700)
 *   - solvent: string (solvant, défaut "CDCl3")
 *   - line_shape: string ("lorentzian" | "gaussian" | "voigt", défaut "lorentzian")
 *   - digital_resolution_hz: number (résolution en Hz, défaut 0.1)
 *   - ppm_range: [number, number] (intervalle ppm, défaut [-0.5, 12])
 * 
 * Response:
 *   - ppm: number[] (axe des ppm)
 *   - intensity: number[] (intensité normalisée)
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
    line_width_hz, // Optionnel : si non fourni, calculé automatiquement par le backend Python
    // Paramètres dynamics post-processeur
    enable_dynamics = false,
    temperature_k = 298,
    ph: ph_raw = 7.4,
  } = req.body;

  // Validation basique du pH (plage chimique valide)
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

  // Vérifier qu'on a au moins une source de données
  if (!smiles && !multiplicity_json && !multiplicity_json_path) {
    return res.status(400).json({ 
      error: 'Either smiles, multiplicity_json, or multiplicity_json_path is required' 
    });
  }

  // Détermine le chemin vers le script Python
  const scriptsDir = path.resolve(__dirname, '../scripts');
  const simulatorScript = path.join(scriptsDir, 'simulate_spectrum_api.py');
  
  // Vérifier que le script existe
  if (!fs.existsSync(simulatorScript)) {
    return res.status(500).json({ error: 'Simulator script not found' });
  }

  // IMPORTANT: L'app doit être 100% autonome - toujours utiliser le Python bundled
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

  // Le Python bundled utilise un système de commandes intégré
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
  
  // Mode 1: JSON de multiplicité inline
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
    // Paramètres dynamics post-processeur
    '--temperature_k', String(temperature_k),
    '--ph', String(ph),
  );

  // line_width_hz optionnel : si fourni par le frontend, l'utiliser ; sinon auto-calculé par Python
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
      // Le script retourne du JSON sur stdout
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
    '--pdata_dir', pdata_dir,
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
    '--pdata_dir', pdata_dir,
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
 * Simule les spectres de TOUTES les molécules d'un fichier all_multiplicities.json
 * 
 * Body:
 *   - task_id: string (UUID de la tâche, pour trouver le fichier JSON)
 *   - freq_mhz: number (fréquence en MHz, défaut 700)
 *   - solvent: string (solvant, défaut "CDCl3")
 *   - line_shape: string ("lorentzian" | "gaussian" | "voigt")
 *   - line_width_hz: number (largeur de ligne en Hz)
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
    line_width_hz, // Optionnel : si non fourni, calculé automatiquement par le backend Python
    enable_dynamics = false,
    temperature_k = 298,
    ph: ph_raw = 7.4,
  } = req.body;

  // ✅ VALIDATION pH ULTRA-ROBUSTE - Prévention bug pH erronés 
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
    ph = 7.4;  // Fallback sécurisé
  }

  console.log('[simulate-spectrum-batch] Request received:', { task_id, freq_mhz });

  if (!task_id) {
    return res.status(400).json({ 
      status: 'error',
      error: 'task_id is required' 
    });
  }

  // Trouve le dossier de la tâche
  const projectService = require('../services/projectService');
  const currentProject = projectService.getCurrent();
  const outputsDir = path.resolve(__dirname, '../../outputs');
  const baseDir = currentProject?.workDir
    ? currentProject.workDir
    : outputsDir;

  const taskDir = path.join(baseDir, task_id);
  const allMultiplicitiesPath = path.join(taskDir, 'all_multiplicities.json');

  console.log('[simulate-spectrum-batch] Looking for:', allMultiplicitiesPath);

  // Chercher le fichier de multiplicité
  let multiplicityFilePath = null;
  let isSingleMode = false;

  if (fs.existsSync(allMultiplicitiesPath)) {
    // Mode batch: all_multiplicities.json existe
    multiplicityFilePath = allMultiplicitiesPath;
    console.log('[simulate-spectrum-batch] Found all_multiplicities.json');
  } else {
    // Mode single: chercher un fichier *_multiplicity.json
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

  // Détermine le chemin Python
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

  // Construire les arguments selon le mode
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

  // line_width_hz optionnel : si fourni par le frontend, l'utiliser ; sinon auto-calculé par Python
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
    // Log mais ne pas échouer (peut contenir des warnings)
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
      
      // Si mode single, wrapper le résultat dans un format batch-like
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
