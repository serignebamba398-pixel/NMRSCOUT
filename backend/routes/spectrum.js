// backend/routes/spectrum.js
// Route pour la simulation de spectre 1H détaillé

const express = require('express');
const router = express.Router();
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');

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
    line_width_hz = 1.0,
  } = req.body;

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

  // Trouver l'interpréteur Python
  const pythonPath = process.env.PYTHON_PATH || 'python';

  // Préparer les arguments selon le mode
  const args = [simulatorScript];
  
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
    '--line_width_hz', String(line_width_hz),
  );

  let stdout = '';
  let stderr = '';

  const child = spawn(pythonPath, args, {
    cwd: scriptsDir,
    env: { ...process.env, PYTHONIOENCODING: 'utf-8' },
  });

  child.stdout.on('data', (data) => {
    stdout += data.toString();
  });

  child.stderr.on('data', (data) => {
    stderr += data.toString();
  });

  child.on('close', (code) => {
    if (code !== 0) {
      console.error('[simulate-spectrum] Python error:', stderr);
      return res.status(500).json({ 
        status: 'error',
        error: 'Simulation failed', 
        details: stderr.substring(0, 1000) 
      });
    }

    try {
      // Le script retourne du JSON sur stdout
      const result = JSON.parse(stdout);
      res.json(result);
    } catch (e) {
      console.error('[simulate-spectrum] JSON parse error:', e.message);
      res.status(500).json({ 
        status: 'error',
        error: 'Failed to parse simulation result',
        details: stdout.substring(0, 500)
      });
    }
  });

  child.on('error', (err) => {
    console.error('[simulate-spectrum] Spawn error:', err);
    res.status(500).json({ error: 'Failed to start Python process', details: err.message });
  });
});

module.exports = router;
