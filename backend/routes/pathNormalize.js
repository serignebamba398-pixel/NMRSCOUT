const express = require('express');
const router = express.Router();
const { normalizeHostPathToContainer } = require('../utils/pathNormalizer');

const normalizePathHandler = (req, res) => {
  try {
    const rawPath = String(req.body?.rawPath || '');
    if (!rawPath) return res.status(400).json({ ok: false, error: 'rawPath is required' });
    
    // In Electron mode, keep native paths (Windows) to avoid transforming user paths.
    // In non-Electron (e.g., Docker) mode, normalize host paths to container paths.
    const out = process.env.IS_ELECTRON === 'true' ? rawPath : normalizeHostPathToContainer(rawPath);
    return res.json({ ok: true, path: out });
  } catch (e) {
    return res.status(500).json({ ok: false, error: String(e?.message || e) });
  }
};

router.post('/api/path/normalize', express.json(), normalizePathHandler);
router.post('/path/normalize', express.json(), normalizePathHandler);

module.exports = router;