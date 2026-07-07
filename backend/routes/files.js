const express = require('express');
const router = express.Router();
const fs = require('fs');
const path = require('path');
const { uploadsDir, outputsDir } = require('../config/paths');
const { safeResolve } = require('../utils/filesystem');
const { logWithMemory } = require('../utils/logger');
const upload = require('../middleware/upload');
const projectService = require('../services/projectService');

router.get('/api/files', (req, res) => {
  try {
    const files = fs
      .readdirSync(uploadsDir)
      .filter(fn => fn.endsWith('.csv') || fn.endsWith('.sdf'))
      .map(fn => {
        const p = path.join(uploadsDir, fn);
        const st = fs.statSync(p);
        return {
          name: fn,
          size: st.size,
          sizeMB: Math.round((st.size / 1024 / 1024) * 100) / 100,
          uploadedAt: st.mtime,
          originalName: fn.split('_').slice(1).join('_'),
        };
      })
      .sort((a, b) => b.uploadedAt - a.uploadedAt);
    res.json(files);
  } catch (e) {
    res.status(500).json({ error: e.message });
  }
});

router.post('/api/upload', upload.single('file'), (req, res) => {
  if (!req.file) return res.status(400).json({ error: 'No file provided' });
  logWithMemory(
    `Upload: ${req.file.filename} (${Math.round((req.file.size / 1024 / 1024) * 100) / 100}MB)`
  );
  res.json({
    message: 'OK',
    filename: req.file.filename,
    originalName: req.file.originalname,
    size: req.file.size,
    sizeMB: Math.round((req.file.size / 1024 / 1024) * 100) / 100,
  });
});

router.delete('/api/files/:filename', (req, res) => {
  try {
    const p = safeResolve(uploadsDir, req.params.filename);
    if (!fs.existsSync(p)) return res.status(404).json({ error: 'File not found' });
    const stats = fs.statSync(p);
    fs.unlinkSync(p);
    logWithMemory(
      `File deleted: ${req.params.filename} (${Math.round((stats.size / 1024 / 1024) * 100) / 100}MB)`
    );
    res.json({ message: 'Deleted' });
  } catch (e) {
    res.status(400).json({ error: e.message });
  }
});

/**
 * GET /api/multiplicity/:taskId/:molId
 * Retrieve multiplicity data for a specific molecule from multiplicity JSON files
 * Supports both formats:
 * - all_multiplicities.json (CSV batch predictions)
 * - {mol_id}_multiplicity.json (single SMILES predictions)
 *
 * @param taskId - UUID of the prediction task
 * @param molId - Molecule ID to find in the JSON
 * @returns JSON object with multiplicity data for the molecule
 */
router.get('/api/multiplicity/:taskId/:molId', (req, res) => {
  console.log('========================================');
  console.log('[MULTIPLICITY] NEW VERSION v2.0 - ENDPOINT HIT!');
  console.log('[MULTIPLICITY] Params:', req.params);
  console.log('========================================');

  try {
    const { taskId, molId } = req.params;

    // Get the current project's outputs directory
    const currentProject = projectService.getCurrent();
    console.log('[MULTIPLICITY] Current project:', currentProject?.name || 'none');
    const baseDir = currentProject?.workDir
      ? path.join(currentProject.workDir, 'outputs')
      : outputsDir;

    console.log('[MULTIPLICITY] baseDir:', baseDir);

    const taskDir = path.join(baseDir, taskId);
    console.log('[MULTIPLICITY] taskDir:', taskDir);
    console.log('[MULTIPLICITY] taskDir exists?', fs.existsSync(taskDir));

    if (!fs.existsSync(taskDir)) {
      console.log('[MULTIPLICITY] ERROR: Task directory not found!');
      return res.status(404).json({
        error: 'Task directory not found',
        path: taskDir,
        baseDir,
        currentProject: currentProject?.name || 'none'
      });
    }

    // Try 1: Look for specific molecule file {molId}_multiplicity.json
    // Try multiple name variations (with/without _predicted suffix)
    const molIdVariants = [
      molId,
      molId.replace(/_predicted$/i, ''),  // Remove _predicted suffix if present
      molId.replace(/_annotated$/i, ''),  // Remove _annotated suffix if present
    ];

    console.log('[MULTIPLICITY] Trying variants:', molIdVariants);

    for (const variant of molIdVariants) {
      const specificPath = path.join(taskDir, `${variant}_multiplicity.json`);
      console.log('[MULTIPLICITY] Trying file:', specificPath, '- exists?', fs.existsSync(specificPath));
      if (fs.existsSync(specificPath)) {
        console.log('[MULTIPLICITY] SUCCESS! Found file:', specificPath);
        const jsonContent = fs.readFileSync(specificPath, 'utf-8');
        const moleculeData = JSON.parse(jsonContent);

        return res.json({
          status: 'ok',
          mol_id: moleculeData.mol_id || molId,
          smiles: moleculeData.smiles || '',
          multiplicity_data: moleculeData
        });
      }
    }

    // Try 2: Look for all_multiplicities.json (batch predictions)
    const allPath = path.join(taskDir, 'all_multiplicities.json');
    if (fs.existsSync(allPath)) {
      const jsonContent = fs.readFileSync(allPath, 'utf-8');
      const allMultiplicities = JSON.parse(jsonContent);

      // Find the molecule by ID (also try without _predicted suffix)
      const molIdClean = molId.replace(/_predicted$/i, '');
      const molecule = allMultiplicities.molecules?.find(m => 
        m.mol_id === molId || m.mol_id === molIdClean
      );

      if (!molecule) {
        return res.status(404).json({
          error: `Molecule ${molId} not found in all_multiplicities.json`,
          available: allMultiplicities.molecules?.map(m => m.mol_id) || []
        });
      }

      return res.json({
        status: 'ok',
        mol_id: molecule.mol_id,
        smiles: molecule.smiles,
        multiplicity_data: molecule
      });
    }

    // Try 3: Search for any *_multiplicity.json or *multiplicities.json file
    const files = fs.readdirSync(taskDir);
    const multFiles = files.filter(f =>
      f.endsWith('_multiplicity.json') || f.endsWith('multiplicities.json')
    );

    if (multFiles.length > 0) {
      // Try each multiplicity file
      for (const file of multFiles) {
        const filePath = path.join(taskDir, file);
        const jsonContent = fs.readFileSync(filePath, 'utf-8');
        const data = JSON.parse(jsonContent);

        // Check if this is a single molecule file
        if (data.mol_id === molId) {
          return res.json({
            status: 'ok',
            mol_id: data.mol_id,
            smiles: data.smiles || '',
            multiplicity_data: data
          });
        }

        // Check if this is a batch file with molecules array
        if (data.molecules && Array.isArray(data.molecules)) {
          const molecule = data.molecules.find(m => m.mol_id === molId);
          if (molecule) {
            return res.json({
              status: 'ok',
              mol_id: molecule.mol_id,
              smiles: molecule.smiles,
              multiplicity_data: molecule
            });
          }
        }
      }

      return res.status(404).json({
        error: `Molecule ${molId} not found in any multiplicity files`,
        files_checked: multFiles
      });
    }

    // No multiplicity files found
    return res.status(404).json({
      error: 'No multiplicity files found in task directory',
      taskDir,
      files: files
    });

  } catch (e) {
    console.error('[multiplicity] Error:', e);
    res.status(500).json({ error: e.message });
  }
});

module.exports = router;