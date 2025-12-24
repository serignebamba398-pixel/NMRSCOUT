const express = require('express');
const router = express.Router();
const fs = require('fs');
const path = require('path');
const { outputsDir } = require('../config/paths');
const { safeResolve } = require('../utils/filesystem');
const { logWithMemory } = require('../utils/logger');
const { getAllTasks, getTask, deleteTask, killProcessGracefully } = require('../services/processManager');
const { sseSend } = require('../middleware/sse');
const { isWin } = require('../config');

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
    const dir = safeResolve(outputsDir, req.params.taskId);
    const filepath = safeResolve(dir, req.params.filename);
    if (!fs.existsSync(filepath)) return res.status(404).json({ error: 'File not found' });
    res.download(filepath, req.params.filename);
  } catch (e) {
    res.status(400).json({ error: e.message });
  }
});

router.get('/api/inspect/:taskId', (req, res) => {
  try {
    const workDir = safeResolve(outputsDir, req.params.taskId);
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

module.exports = router;