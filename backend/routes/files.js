const express = require('express');
const router = express.Router();
const fs = require('fs');
const path = require('path');
const { uploadsDir, outputsDir } = require('../config/paths');
const { safeResolve } = require('../utils/filesystem');
const { logWithMemory } = require('../utils/logger');
const upload = require('../middleware/upload');

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
 * Retrieve multiplicity data for a specific molecule from all_multiplicities.json
 *
 * @param taskId - UUID of the prediction task
 * @param molId - Molecule ID to find in the JSON
 * @returns JSON object with multiplicity data for the molecule
 */
router.get('/api/multiplicity/:taskId/:molId', (req, res) => {
  try {
    const { taskId, molId } = req.params;

    // Construct path to all_multiplicities.json
    const jsonPath = path.join(outputsDir, taskId, 'all_multiplicities.json');

    if (!fs.existsSync(jsonPath)) {
      return res.status(404).json({
        error: 'Multiplicity file not found',
        path: jsonPath
      });
    }

    // Read and parse the JSON file
    const jsonContent = fs.readFileSync(jsonPath, 'utf-8');
    const allMultiplicities = JSON.parse(jsonContent);

    // Find the molecule by ID
    const molecule = allMultiplicities.molecules?.find(m => m.mol_id === molId);

    if (!molecule) {
      return res.status(404).json({
        error: `Molecule ${molId} not found in multiplicity data`,
        available: allMultiplicities.molecules?.map(m => m.mol_id) || []
      });
    }

    // Return the multiplicity data for this molecule
    res.json({
      status: 'ok',
      mol_id: molecule.mol_id,
      smiles: molecule.smiles,
      multiplicity_data: molecule.multiplicity_data
    });

  } catch (e) {
    console.error('[multiplicity] Error:', e);
    res.status(500).json({ error: e.message });
  }
});

module.exports = router;