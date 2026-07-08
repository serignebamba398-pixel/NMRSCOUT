const express = require('express');
const router = express.Router();
const fs = require('fs');
const path = require('path');
const { outputsDir, getActiveOutputsDir } = require('../config/paths');
const { safeResolve } = require('../utils/filesystem');
const { logWithMemory } = require('../utils/logger');
const { getAllTasks, getTask, deleteTask, killProcessGracefully } = require('../services/processManager');
const { sseSend } = require('../middleware/sse');
const { isWin } = require('../config');
const { getResumableTasks, loadCheckpoint, canResumeTask, getResumeData, removeFromPendingTasks, loadPendingTasks } = require('../services/taskPersistence');

router.get('/api/tasks', (req, res) => {
  const tasks = getAllTasks().map(([id, t]) => ({
    id,
    type: t.annotationType || t.mode || 'unknown',
    start: t.startTime,
    duration: Date.now() - t.startTime.getTime(),
    durationMin: Math.round((Date.now() - t.startTime.getTime()) / 60000 * 10) / 10,
    workDir: t.workDir,
    pid: t.process ? t.process.pid : null,
    parameters: t.parameters,
  }));
  res.json(tasks);
});

// GET /api/tasks/:taskId - minimal status (useful for frontend diagnostics)
router.get('/api/tasks/:taskId', (req, res) => {
  const taskId = req.params.taskId;
  const task = getTask(taskId);

  if (!task) {
    return res.status(404).json({ taskId, status: 'not_found' });
  }

  const pid = task.process ? task.process.pid : null;
  const status = pid ? 'running' : 'unknown';

  return res.json({
    taskId,
    status,
    pid,
    type: task.annotationType || task.mode || 'unknown',
    start: task.startTime,
    workDir: task.workDir,
  });
});

router.post('/api/tasks/:taskId/cancel', async (req, res) => {
  const taskId = req.params.taskId;
  const task = getTask(taskId);

  if (!task) {
    return res.status(404).json({
      error: `Task ${taskId} not found or already finished`,
      taskId,
      canceled: false,
    });
  }

  try {
    logWithMemory(`[${taskId}] Cancellation requested by client`, 'INFO');
    task.canceled = true;

    if (task.process && task.process.pid) {
      await new Promise(resolve => {
        const pid = task.process.pid;
        try {
          process.kill(pid, 0);
        } catch {
          logWithMemory(`[${taskId}] Process ${pid} already exited`, 'INFO');
          resolve();
          return;
        }

        logWithMemory(`[${taskId}] Sending SIGTERM to PID ${pid}`, 'INFO');

        const onExit = () => {
          logWithMemory(`[${taskId}] Process ${pid} exited`, 'INFO');
          resolve();
        };

        task.process.once('exit', onExit);
        task.process.once('close', onExit);

        try {
          if (isWin) {
            const { killProcessTree } = require('../services/processManager');
            killProcessTree(pid, 'SIGTERM', taskId);
          } else {
            task.process.kill('SIGTERM');
          }
        } catch (e) {
          logWithMemory(`[${taskId}] SIGTERM error: ${e.message}`, 'WARN');
          resolve();
          return;
        }

        setTimeout(() => {
          try {
            process.kill(pid, 0);
            logWithMemory(`[${taskId}] Process still running, escalating`, 'WARN');
            if (isWin) {
              const { killProcessTree } = require('../services/processManager');
              killProcessTree(pid, 'SIGKILL', taskId);
            } else {
              task.process.kill('SIGKILL');
            }
            setTimeout(resolve, 1000);
          } catch {
            resolve();
          }
        }, 3000);
      });
    }

    if (task.safetyTimeout) {
      clearTimeout(task.safetyTimeout);
      task.safetyTimeout = null;
    }

    if (task.heartbeat) {
      try {
        clearInterval(task.heartbeat);
        task.heartbeat = null;
      } catch (e) {
        logWithMemory(`[${taskId}] Heartbeat cleanup error: ${e.message}`, 'WARN');
      }
    }

    if (task.sseResponse && !task.sseResponse.writableEnded && !task.sseResponse.destroyed) {
      try {
        sseSend(task.sseResponse, {
          type: 'completed',
          success: false,
          error: 'Canceled by user',
          canceled: true,
          taskId,
        });
        task.sseResponse.end();
      } catch (e) {
        logWithMemory(`[${taskId}] SSE end error: ${e.message}`, 'WARN');
      }
    }

    deleteTask(taskId);
    logWithMemory(`[${taskId}] Task canceled successfully`, 'INFO');

    res.json({ message: 'Task canceled successfully', taskId, canceled: true });
  } catch (e) {
    logWithMemory(`[${taskId}] Cancellation error: ${e.message}`, 'ERROR');
    res.status(500).json({ error: e.message, taskId, canceled: false });
  }
});

router.post('/api/tasks/:taskId/stop', (req, res) => {
  const taskId = req.params.taskId;
  const task = getTask(taskId);
  if (!task) return res.status(404).json({ error: 'Task not found' });
  try {
    if (task.safetyTimeout) clearTimeout(task.safetyTimeout);
    if (task.process) {
      killProcessGracefully(task.process, taskId);
    }
    deleteTask(taskId);
    res.json({ message: 'Stop requested', taskId });
  } catch (e) {
    res.status(500).json({ error: e.message });
  }
});

router.get('/api/download/:taskId/:filename', (req, res) => {
  try {
    const dir = safeResolve(getActiveOutputsDir(), req.params.taskId);
    const filepath = safeResolve(dir, req.params.filename);
    if (!fs.existsSync(filepath)) return res.status(404).json({ error: 'File not found' });
    res.download(filepath, req.params.filename);
  } catch (e) {
    res.status(400).json({ error: e.message });
  }
});

router.get('/api/inspect/:taskId', (req, res) => {
  try {
    const workDir = safeResolve(getActiveOutputsDir(), req.params.taskId);
    if (!fs.existsSync(workDir)) return res.status(404).json({ error: 'Task not found' });
    const items = fs
      .readdirSync(workDir)
      .map(name => {
        const p = path.join(workDir, name);
        const st = fs.statSync(p);
        return {
          name,
          size: st.isFile() ? st.size : null,
          sizeMB: st.isFile() ? Math.round((st.size / 1024 / 1024) * 100) / 100 : null,
          mtime: st.mtime,
          isDir: st.isDirectory(),
          isCSV: name.endsWith('.csv'),
        };
      })
      .sort((a, b) => b.mtime - a.mtime);
    res.json({ taskId: req.params.taskId, workDir, items });
  } catch (e) {
    res.status(400).json({ error: e.message });
  }
});

// NEW: Get annotation parameters for a task
router.get('/api/tasks/:taskId/params', (req, res) => {
  try {
    const taskId = req.params.taskId;
    const workDir = safeResolve(getActiveOutputsDir(), taskId);
    const paramsPath = path.join(workDir, 'params.json');
    
    if (!fs.existsSync(paramsPath)) {
      return res.status(404).json({ error: 'Parameters file not found', taskId });
    }
    
    const params = JSON.parse(fs.readFileSync(paramsPath, 'utf8'));
    res.json({ taskId, parameters: params });
  } catch (e) {
    res.status(400).json({ error: e.message });
  }
});

// ============================================================================
// TASK PERSISTENCE ROUTES
// ============================================================================

/**
 * GET /api/tasks/pending
 * Liste toutes les tâches qui peuvent être reprises après une déconnexion
 */
router.get('/api/tasks/pending', (req, res) => {
  try {
    const pending = getResumableTasks();
    res.json({
      count: pending.length,
      tasks: pending,
      message: pending.length > 0 
        ? `${pending.length} tâche(s) en attente peuvent être reprises` 
        : 'Aucune tâche en attente'
    });
  } catch (e) {
    logWithMemory(`Error listing pending tasks: ${e.message}`, 'ERROR');
    res.status(500).json({ error: e.message });
  }
});

/**
 * GET /api/tasks/pending/count
 * Retourne simplement le nombre de tâches en attente (léger pour polling)
 */
router.get('/api/tasks/pending/count', (req, res) => {
  try {
    const pending = loadPendingTasks();
    res.json({ count: Object.keys(pending).length });
  } catch (e) {
    res.json({ count: 0 });
  }
});

/**
 * GET /api/tasks/:taskId/checkpoint
 * Récupère le checkpoint d'une tâche spécifique
 */
router.get('/api/tasks/:taskId/checkpoint', (req, res) => {
  try {
    const taskId = req.params.taskId;
    const checkpoint = loadCheckpoint(taskId);
    
    if (!checkpoint) {
      return res.status(404).json({ 
        error: 'No checkpoint found for this task',
        taskId 
      });
    }
    
    res.json({
      taskId,
      checkpoint,
      canResume: canResumeTask(taskId),
    });
  } catch (e) {
    res.status(400).json({ error: e.message });
  }
});

/**
 * GET /api/tasks/:taskId/status
 * Vérifie le statut d'une tâche (en cours, terminée, interrompue)
 */
router.get('/api/tasks/:taskId/status', (req, res) => {
  try {
    const taskId = req.params.taskId;
    
    // Vérifier si la tâche est en cours d'exécution
    const runningTask = getTask(taskId);
    if (runningTask) {
      return res.json({
        taskId,
        status: 'running',
        progress: runningTask.progress || null,
        startTime: runningTask.startTime,
      });
    }
    
    // Vérifier le checkpoint
    const checkpoint = loadCheckpoint(taskId);
    if (checkpoint) {
      return res.json({
        taskId,
        status: checkpoint.status,
        progress: checkpoint.progress,
        completedAt: checkpoint.completedAt,
        outputs: checkpoint.outputs,
      });
    }
    
    // Vérifier si le dossier existe (tâche terminée sans checkpoint)
    const workDir = path.join(getActiveOutputsDir(), taskId);
    if (fs.existsSync(workDir)) {
      return res.json({
        taskId,
        status: 'completed_legacy',
        message: 'Task completed (no checkpoint available)',
      });
    }
    
    return res.status(404).json({ 
      error: 'Task not found',
      taskId 
    });
  } catch (e) {
    res.status(400).json({ error: e.message });
  }
});

/**
 * DELETE /api/tasks/:taskId/pending
 * Supprime une tâche de la liste des tâches en attente (sans la reprendre)
 */
router.delete('/api/tasks/:taskId/pending', (req, res) => {
  try {
    const taskId = req.params.taskId;
    removeFromPendingTasks(taskId);
    res.json({ 
      message: 'Task removed from pending list',
      taskId 
    });
  } catch (e) {
    res.status(400).json({ error: e.message });
  }
});

/**
 * POST /api/tasks/:taskId/dismiss
 * Marque une tâche pendante comme "vue" (l'utilisateur ne veut pas la reprendre)
 */
router.post('/api/tasks/:taskId/dismiss', (req, res) => {
  try {
    const taskId = req.params.taskId;
    removeFromPendingTasks(taskId);
    logWithMemory(`[${taskId}] Task dismissed by user`, 'INFO');
    res.json({
      message: 'Task dismissed',
      taskId
    });
  } catch (e) {
    res.status(400).json({ error: e.message });
  }
});

// ============================================================================
// TOLERANCE OPTIMIZATION RESULTS PERSISTENCE
// ============================================================================

/**
 * POST /api/tasks/:taskId/tolerance-results
 * Sauvegarde les résultats d'optimisation de tolérance
 */
router.post('/api/tasks/:taskId/tolerance-results', (req, res) => {
  try {
    const taskId = req.params.taskId;
    const results = req.body.results;

    if (!results || !Array.isArray(results)) {
      return res.status(400).json({ error: 'Invalid results format (expected array)' });
    }

    const taskDir = safeResolve(getActiveOutputsDir(), taskId);
    if (!fs.existsSync(taskDir)) {
      fs.mkdirSync(taskDir, { recursive: true });
    }

    const resultsPath = path.join(taskDir, 'tolerance_results.json');
    fs.writeFileSync(resultsPath, JSON.stringify(results, null, 2), 'utf8');

    logWithMemory(`[${taskId}] Saved ${results.length} tolerance optimization result(s)`, 'INFO');

    res.json({
      message: 'Tolerance results saved successfully',
      taskId,
      count: results.length,
      path: resultsPath
    });
  } catch (e) {
    logWithMemory(`Error saving tolerance results: ${e.message}`, 'ERROR');
    res.status(500).json({ error: e.message });
  }
});

/**
 * GET /api/tasks/:taskId/tolerance-results
 * Récupère les résultats d'optimisation de tolérance sauvegardés
 */
router.get('/api/tasks/:taskId/tolerance-results', (req, res) => {
  try {
    const taskId = req.params.taskId;
    const taskDir = safeResolve(getActiveOutputsDir(), taskId);
    const resultsPath = path.join(taskDir, 'tolerance_results.json');

    if (!fs.existsSync(resultsPath)) {
      return res.json({
        taskId,
        results: [],
        message: 'No saved tolerance results found'
      });
    }

    const results = JSON.parse(fs.readFileSync(resultsPath, 'utf8'));

    res.json({
      taskId,
      results,
      count: results.length
    });
  } catch (e) {
    logWithMemory(`Error loading tolerance results: ${e.message}`, 'ERROR');
    res.status(500).json({ error: e.message });
  }
});

/**
 * DELETE /api/tasks/:taskId/tolerance-results
 * Supprime les résultats d'optimisation de tolérance sauvegardés
 */
router.delete('/api/tasks/:taskId/tolerance-results', (req, res) => {
  try {
    const taskId = req.params.taskId;
    const taskDir = safeResolve(getActiveOutputsDir(), taskId);
    const resultsPath = path.join(taskDir, 'tolerance_results.json');

    if (!fs.existsSync(resultsPath)) {
      return res.json({
        message: 'No tolerance results to delete',
        taskId
      });
    }

    fs.unlinkSync(resultsPath);
    logWithMemory(`[${taskId}] Deleted tolerance results`, 'INFO');

    res.json({
      message: 'Tolerance results deleted successfully',
      taskId
    });
  } catch (e) {
    logWithMemory(`Error deleting tolerance results: ${e.message}`, 'ERROR');
    res.status(500).json({ error: e.message });
  }
});

module.exports = router;