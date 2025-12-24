const express = require('express');
const router = express.Router();
const fs = require('fs');
const path = require('path');
const { uploadsDir, outputsDir } = require('../config/paths');
const { TASK_TTL_HOURS } = require('../config');
const { getDirSize } = require('../utils/filesystem');
const { logWithMemory } = require('../utils/logger');

router.post('/api/cleanup', (req, res) => {
  try {
    const maxAge = TASK_TTL_HOURS * 3600 * 1000;
    const now = Date.now();
    let cleanedFiles = 0,
      cleanedDirs = 0,
      freedMB = 0;

    fs.readdirSync(uploadsDir).forEach(filename => {
      const p = path.join(uploadsDir, filename);
      const st = fs.statSync(p);
      if (st.isFile() && now - st.mtime.getTime() > maxAge) {
        freedMB += Math.round((st.size / 1024 / 1024) * 100) / 100;
        fs.unlinkSync(p);
        cleanedFiles++;
      }
    });

    fs.readdirSync(outputsDir).forEach(dirname => {
      const dirpath = path.join(outputsDir, dirname);
      const st = fs.statSync(dirpath);
      if (st.isDirectory() && now - st.mtime.getTime() > maxAge) {
        const sz = getDirSize(dirpath);
        freedMB += sz;
        fs.rmSync(dirpath, { recursive: true, force: true });
        cleanedDirs++;
      }
    });

    logWithMemory(`Cleanup: ${cleanedFiles} files + ${cleanedDirs} folders (${freedMB}MB)`);
    res.json({
      message: 'Cleanup completed',
      cleanedFiles,
      cleanedDirs,
      freedMB: Math.round(freedMB * 100) / 100,
    });
  } catch (e) {
    res.status(500).json({ error: e.message });
  }
});

module.exports = router;