const express = require('express');
const router = express.Router();
const fs = require('fs');
const path = require('path');
const { getSystemMemoryUsage } = require('../utils/logger');
const { getRscript, getPythonExe } = require('../utils/spawn');
const { getTasksCount } = require('../services/processManager');
const {
  PROCESS_TIMEOUT_HOURS,
  MEMORY_LIMIT_MB,
  ENABLE_CHECKPOINTS,
  KILL_ON_CLIENT_CLOSE,
} = require('../config');
const { scriptsDir, uploadsDir, outputsDir, scriptMap, PRED_SCRIPT } = require('../config/paths');

// Check if Python is available
function checkPythonAvailable() {
  const pythonPath = getPythonExe();
  const isBundled = process.env.IS_ELECTRON_BUNDLED === 'true';
  
  // If bundled, check if exe exists
  if (isBundled) {
    return fs.existsSync(pythonPath);
  }
  
  // If system Python, check if script exists
  if (!fs.existsSync(PRED_SCRIPT)) {
    return false;
  }
  
  // Could also try to spawn python --version but that's slow
  return true;
}

router.get('/health', (req, res) => {
  const memory = getSystemMemoryUsage();
  let pythonPath = null;
  let pythonAvailable = false;
  let pythonError = null;
  try {
    pythonPath = getPythonExe();
    pythonAvailable = checkPythonAvailable();
  } catch (e) {
    pythonPath = null;
    pythonAvailable = false;
    pythonError = String(e?.message || e);
  }
  
  let rPath = null;
  let rAvailable = false;
  try {
    rPath = getRscript();
    rAvailable = fs.existsSync(rPath);
  } catch (e) {
    // R-Portable not found - this is an error, app must be autonomous
    rPath = null;
    rAvailable = false;
  }
  
  res.json({
    status: 'OK',
    timestamp: new Date().toISOString(),
    uptime: process.uptime(),
    memory,
    capabilities: {
      predictions: pythonAvailable,
      annotations: rAvailable,
    },
    python: {
      path: pythonPath,
      available: pythonAvailable,
      bundled: process.env.IS_ELECTRON_BUNDLED === 'true',
      scriptExists: fs.existsSync(PRED_SCRIPT),
      error: pythonError,
    },
    r: {
      path: rPath,
      available: rAvailable,
      bundled: true, // App only uses bundled R-Portable
    },
    node: process.version,
    platform: process.platform,
    isElectron: process.env.IS_ELECTRON === 'true',
    runningTasks: getTasksCount(),
    config: {
      processTimeoutHours: PROCESS_TIMEOUT_HOURS,
      memoryLimitMB: MEMORY_LIMIT_MB,
      enableCheckpoints: ENABLE_CHECKPOINTS,
      killOnClientClose: KILL_ON_CLIENT_CLOSE,
    },
  });
});

router.get('/api/diagnostic', (req, res) => {
  const rStatus = {};
  for (const [k, v] of Object.entries(scriptMap)) {
    rStatus[k] = fs.existsSync(path.join(scriptsDir, v));
  }

  const { getAllTasks } = require('../services/processManager');
  const taskDetails = getAllTasks().map(([id, t]) => ({
    id,
    type: t.annotationType || t.mode || 'unknown',
    start: t.startTime,
    duration: Date.now() - t.startTime.getTime(),
    workDir: t.workDir,
    pid: t.process ? t.process.pid : null,
    killed: t.process ? t.process.killed : null,
  }));

  res.json({
    status: 'OK',
    timestamp: new Date().toISOString(),
    system: {
      platform: process.platform,
      node: process.version,
      rscript: getRscript(),
      uploadsDir,
      outputsDir,
      scriptsDir,
      memory: getSystemMemoryUsage(),
      uptime: process.uptime(),
    },
    scriptsR: rStatus,
    tasks: taskDetails,
  });
});

module.exports = router;