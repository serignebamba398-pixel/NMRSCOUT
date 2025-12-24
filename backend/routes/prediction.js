const express = require('express');
const router = express.Router();
const fs = require('fs');
const path = require('path');
const { outputsDir } = require('../config/paths');
const { handlePrediction } = require('../services/predictionService');

router.post('/api/predict/start', handlePrediction);

router.get('/api/prediction-tasks', (req, res) => {
  try {
    const tasks = [];
    const outputDirs = fs
      .readdirSync(outputsDir)
      .filter(dir => {
        const dirPath = path.join(outputsDir, dir);
        return fs.statSync(dirPath).isDirectory();
      });

    const findCsv = (files, token) => {
      const rx = new RegExp(`^prediction_${token}_.+\\.csv$`, 'i');
      const exact = files.find(f => rx.test(f));
      if (exact) return exact;
      const tok = String(token).toLowerCase();
      return (
        files.find(f => f.toLowerCase().endsWith('.csv') && f.toLowerCase().includes(`prediction_${tok}_`)) ||
        null
      );
    };

    for (const taskId of outputDirs) {
      const taskDir = path.join(outputsDir, taskId);
      const files = fs.readdirSync(taskDir);

      const csvFiles = {
        '1H': findCsv(files, '1H'),
        '13C': findCsv(files, '13C'),
        COSY: findCsv(files, 'COSY'),
        HSQC: findCsv(files, 'HSQC'),
      };

      const sdfFiles = files.filter(f => f.toLowerCase().endsWith('.sdf'));
      const stats = fs.statSync(taskDir);

      tasks.push({
        taskId,
        createdAt: stats.mtime,
        csvFiles,
        sdfFiles,
        hasSpectralMatrix: false,
        canRunAnnotations: Object.values(csvFiles).some(Boolean),
      });
    }

    tasks.sort((a, b) => new Date(b.createdAt) - new Date(a.createdAt));
    res.json(tasks);
  } catch (e) {
    res.status(500).json({ error: e.message });
  }
});

module.exports = router;