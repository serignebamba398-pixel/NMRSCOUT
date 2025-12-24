const express = require('express');
const router = express.Router();
const { normalizeHostPathToContainer } = require('../utils/pathNormalizer');
const { getPythonExe } = require('../utils/spawn');

const normalizePathHandler = (req, res) => {
  try {
    const rawPath = String(req.body?.rawPath || '');
    if (!rawPath) return res.status(400).json({ ok: false, error: 'rawPath is required' });
    
    // In bundled mode (Electron), don't transform paths - keep native Windows paths
    const PY = getPythonExe();
    const isBundledExe = PY.toLowerCase().endsWith('nmr_engine.exe') || 
                         PY.toLowerCase().endsWith('nmr_engine') ||
                         PY.toLowerCase().endsWith('nmr_engine_production.exe') ||
                         PY.toLowerCase().endsWith('nmr_engine_production');
    
    const out = isBundledExe ? rawPath : normalizeHostPathToContainer(rawPath);
    return res.json({ ok: true, path: out });
  } catch (e) {
    return res.status(500).json({ ok: false, error: String(e?.message || e) });
  }
};

router.post('/api/path/normalize', express.json(), normalizePathHandler);
router.post('/path/normalize', express.json(), normalizePathHandler);

module.exports = router;