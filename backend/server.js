require('dotenv').config();
'use strict';

const fs = require('fs');
const path = require('path');

const { logWithMemory, logEvent, redactEnv } = require('./utils/logger');

function getStartupLogPath() {
  const dir = process.env.LOG_DIR;
  if (dir) {
    try { fs.mkdirSync(dir, { recursive: true }); } catch {}
    return path.join(dir, 'backend-startup.log.txt');
  }
  return path.join(process.cwd(), 'backend-startup.log.txt');
}

const debugLogPath = getStartupLogPath();
function debugLog(msg) {
  try {
    fs.appendFileSync(debugLogPath, `[${new Date().toISOString()}] ${msg}\n`);
  } catch (e) {}
}

debugLog('Backend starting...');
debugLog(`__dirname: ${__dirname}`);
debugLog(`process.cwd(): ${process.cwd()}`);
debugLog(`ENV (redacted): ${JSON.stringify(redactEnv(process.env), null, 2)}`);

// Global error handlers to catch startup crashes
process.on('uncaughtException', (err) => {
  console.error('[Backend Critical] Uncaught Exception:', err);
  debugLog(`[CRITICAL] Uncaught Exception: ${err.message}\n${err.stack}`);
  logEvent('ERROR', 'Uncaught exception', { error: err?.message, stack: err?.stack });
  process.exit(1);
});

process.on('unhandledRejection', (reason, promise) => {
  console.error('[Backend Critical] Unhandled Rejection at:', promise, 'reason:', reason);
  debugLog(`[CRITICAL] Unhandled Rejection: ${reason}`);
  logEvent('ERROR', 'Unhandled promise rejection', { reason: String(reason) });
  process.exit(1);
});

console.log('[Backend] Starting up...');

let express;
try {
  debugLog('Attempting to require express...');
  express = require('express');
  debugLog('Express loaded successfully.');
} catch (e) {
  console.error('[Backend Critical] Failed to load express:', e);
  debugLog(`[CRITICAL] Failed to load express: ${e.message}\n${e.stack}`);
  debugLog(`require.resolve paths: ${JSON.stringify(require.resolve.paths('express'))}`);
  process.exit(1);
}

const { PORT, BODY_LIMIT_MB, AUTO_CLEANUP_INTERVAL_HOURS, TASK_TTL_HOURS } = require('./config');
const { corsMiddleware } = require('./config/cors');
const helmetConfig = require('./config/helmet');
const { outputsDir, uploadsDir } = require('./config/paths');
const { notFoundHandler, errorHandler } = require('./middleware/errorHandler');
const { getSystemMemoryUsage } = require('./utils/logger');
const { getTasksCount } = require('./services/processManager');

// Import routes
const healthRoutes = require('./routes/health');
const filesRoutes = require('./routes/files');
const predictionRoutes = require('./routes/prediction');
const annotationRoutes = require('./routes/annotation');
const tasksRoutes = require('./routes/tasks');
const spectralBgRoutes = require('./routes/spectralBg');
const cleanupRoutes = require('./routes/cleanup');
const pathNormalizeRoutes = require('./routes/pathNormalize');
const projectRoutes = require('./routes/projects');
const spectrumRoutes = require('./routes/spectrum');  // Simulation de spectre 1H

const shouldServeFrontend = process.env.SERVE_FRONTEND === 'true';
const frontendBuildDir =
  process.env.FRONTEND_BUILD_DIR || path.join(__dirname, '..', 'frontend', 'build');

// Init app
const app = express();
process.setMaxListeners(0);

// Security & CORS
app.use(helmetConfig);
app.use(corsMiddleware);
app.options(/.*/, corsMiddleware);

// Body parsers
app.use(express.json({ limit: `${BODY_LIMIT_MB}mb` }));
app.use(express.urlencoded({ extended: true, limit: `${BODY_LIMIT_MB}mb`, parameterLimit: 200000 }));

// Dynamic static files - serve from current project's folder or default
const projectService = require('./services/projectService');

// Helper to get current outputs directory (project-specific or default)
// Note: outputsDir is already set correctly via OUTPUTS_DIR environment variable
// when Electron starts the backend with a project
function getCurrentOutputsDir() {
  // First priority: environment variable set by Electron
  if (process.env.OUTPUTS_DIR && fs.existsSync(process.env.OUTPUTS_DIR)) {
    return process.env.OUTPUTS_DIR;
  }
  
  // Second priority: projectService current project (for API-based project switching)
  const current = projectService.getCurrent();
  if (current && current.workDir) {
    const projectOutputs = path.join(current.workDir, 'outputs');
    if (fs.existsSync(projectOutputs)) {
      return projectOutputs;
    }
  }
  
  // Fallback to default outputsDir
  return outputsDir;
}

// Dynamic outputs route - serves files from current project
app.use('/api/outputs', (req, res, next) => {
  const dir = getCurrentOutputsDir();
  express.static(dir, { maxAge: '1h', etag: true, lastModified: true })(req, res, next);
});
app.use('/outputs', (req, res, next) => {
  const dir = getCurrentOutputsDir();
  express.static(dir, { maxAge: '1h', etag: true, lastModified: true })(req, res, next);
});

// Static uploads (these can stay static as they're app-level)
app.use('/api/uploads', express.static(uploadsDir, { maxAge: '1h', etag: true, lastModified: true }));
app.use('/uploads', express.static(uploadsDir, { maxAge: '1h', etag: true, lastModified: true }));

// Routes
app.use(healthRoutes);
app.use(filesRoutes);
app.use(predictionRoutes);
app.use(annotationRoutes);
app.use(tasksRoutes);
app.use(spectralBgRoutes);
app.use(cleanupRoutes);
app.use(pathNormalizeRoutes);
app.use(projectRoutes);
app.use(spectrumRoutes);  // Simulation de spectre 1H

// In packaged/Electron mode, serve the React build from the backend so the renderer can load via http://
if (shouldServeFrontend && fs.existsSync(frontendBuildDir)) {
  app.use(express.static(frontendBuildDir, { etag: true, lastModified: true }));
  app.get(/^\/(?!api\/).*/, (req, res) => {
    res.sendFile(path.join(frontendBuildDir, 'index.html'));
  });
}

// Error handlers
app.use(notFoundHandler);
app.use(errorHandler);

// Start server
const server = app.listen(PORT, () => {
  logWithMemory(`Server started on port ${PORT}`);
  logWithMemory(
    `Endpoints: /api/run-annotation, /api/run-annotation-from-prediction, /api/predict/start, /api/prediction-tasks, /api/tasks/:taskId/cancel, /api/diagnostic, /health, /api/files, /api/tasks, /api/cleanup, /api/download/:taskId/:filename, /api/inspect/:taskId, /api/spectral-bg/export, /api/spectral-bg/zoom, /api/spectral-bg/relevel, /api/spectral-bg/interact, /api/path/normalize (+ alias /path/normalize)`
  );
  
  // ✅ PERSISTENCE: Initialiser le système de persistance des tâches au démarrage
  try {
    const { initializeOnStartup } = require('./services/taskPersistence');
    initializeOnStartup();
  } catch (e) {
    logWithMemory(`Task persistence init error: ${e.message}`, 'WARN');
  }
  
  // Notify parent process (Electron) that backend is ready
  if (process.send) {
    process.send({ type: 'ready', port: PORT });
    logWithMemory('IPC "ready" message sent to parent process');
  }
});

if (typeof server.setTimeout === 'function') {
  server.setTimeout(0);
  logWithMemory('server.setTimeout(0)');
}
if ('headersTimeout' in server) {
  server.headersTimeout = 0;
  logWithMemory('server.headersTimeout=0');
}
if ('keepAliveTimeout' in server) {
  server.keepAliveTimeout = 0;
  logWithMemory('server.keepAliveTimeout=0');
}
if ('requestTimeout' in server) {
  server.requestTimeout = 0;
  logWithMemory('server.requestTimeout=0');
}

// Process monitoring
const { ENABLE_PROCESS_MONITORING, PROCESS_TIMEOUT_HOURS } = require('./config');
const { getAllTasks, killProcessGracefully, deleteTask } = require('./services/processManager');

if (ENABLE_PROCESS_MONITORING) {
  setInterval(() => {
    for (const [taskId, task] of getAllTasks()) {
      if (task.process && task.process.pid) {
        try {
          process.kill(task.process.pid, 0);
          const duration = Date.now() - task.startTime.getTime();
          const durationHours = duration / 3_600_000;
          if (durationHours > PROCESS_TIMEOUT_HOURS) {
            logWithMemory(`Security timeout reached ${taskId}`, 'WARN');
            killProcessGracefully(task.process, taskId);
            deleteTask(taskId);
          }
          const mins10 = Math.floor(duration / 600000);
          if ((task.lastLogMinutes || 0) < mins10) {
            logWithMemory(`Task ${taskId} active for ${(duration / 60000).toFixed(1)} min`, 'INFO');
            task.lastLogMinutes = mins10;
          }
        } catch {
          logWithMemory(`Process ${taskId} disparu, cleanup`, 'WARN');
          deleteTask(taskId);
        }
      }
    }
  }, 30000);
}

// Auto cleanup
setInterval(() => {
  try {
    const req = { body: {}, method: 'POST', url: '/api/cleanup' };
    const res = {
      json: d => logWithMemory(`Auto-cleanup: ${JSON.stringify(d)}`),
      status: () => ({ json: d => logWithMemory(`Auto-cleanup error: ${JSON.stringify(d)}`) }),
    };
    app._router.handle(req, res, () => {});
  } catch (e) {
    logWithMemory(`Auto-cleanup error: ${e.message}`, 'WARN');
  }
}, AUTO_CLEANUP_INTERVAL_HOURS * 3600 * 1000);

// Status monitoring
setInterval(() => {
  const memory = getSystemMemoryUsage();
  const uptime = process.uptime();
  const tasksCount = getTasksCount();
  if (tasksCount > 0 || memory.heapUsed > 1000) {
    logWithMemory(
      `Status: ${tasksCount} active tasks, uptime ${Math.floor(uptime / 3600)}h${Math.floor((uptime % 3600) / 60)}m`
    );
  }
}, 5 * 60 * 1000);

module.exports = app;