const express = require('express');
const router = express.Router();
const fs = require('fs');
const path = require('path');
const { uploadsDir } = require('../config/paths');
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

module.exports = router;