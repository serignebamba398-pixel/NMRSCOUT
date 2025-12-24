const { exec } = require('child_process');
const { isWin } = require('../config');
const { logWithMemory } = require('../utils/logger');

const runningTasks = new Map();

function killProcessTree(pid, signal = 'SIGTERM', taskId = 'unknown') {
  if (!pid) return;
  try {
    if (isWin) {
      const force = signal === 'SIGKILL' ? '/F' : '';
      const cmd = `taskkill /PID ${pid} /T ${force}`.trim();
      exec(cmd, err => {
        if (err) logWithMemory(`[${taskId}] taskkill error: ${err.message}`, 'WARN');
      });
    } else {
      try {
        process.kill(pid, 0);
      } catch {
        return;
      }
      try {
        process.kill(-pid, signal);
      } catch {
        try {
          process.kill(pid, signal);
        } catch {}
      }
    }
  } catch (e) {
    logWithMemory(`[${taskId}] kill tree error: ${e.message}`, 'WARN');
  }
}

function killProcessGracefully(proc, taskId) {
  if (!proc || !proc.pid) {
    logWithMemory(`[${taskId}] No PID to kill`, 'WARN');
    return;
  }
  const pid = proc.pid;

  try {
    process.kill(pid, 0);
  } catch {
    logWithMemory(`[${taskId}] Process ${pid} already exited`, 'INFO');
    return;
  }

  logWithMemory(`[${taskId}] Sending SIGTERM to PID ${pid}`, 'INFO');
  try {
    if (!isWin) {
      try {
        process.kill(-pid, 'SIGTERM');
      } catch {
        proc.kill('SIGTERM');
      }
    } else {
      setTimeout(() => killProcessTree(pid, 'SIGTERM', taskId), 200);
    }
  } catch (e) {
    logWithMemory(`[${taskId}] SIGTERM error: ${e.message}`, 'WARN');
  }

  setTimeout(() => {
    try {
      process.kill(pid, 0);
      logWithMemory(`[${taskId}] Process still running, escalating to KILL`, 'WARN');
      if (!isWin) {
        try {
          process.kill(-pid, 'SIGKILL');
        } catch {
          proc.kill('SIGKILL');
        }
      } else {
        killProcessTree(pid, 'SIGKILL', taskId);
      }
    } catch {}
  }, 5000);
}

function addTask(taskId, taskInfo) {
  runningTasks.set(taskId, taskInfo);
}

function getTask(taskId) {
  return runningTasks.get(taskId);
}

function deleteTask(taskId) {
  const task = runningTasks.get(taskId);
  runningTasks.delete(taskId);
  return task; // Return the task so callers can check its state (canceled, etc.)
}

function getAllTasks() {
  return Array.from(runningTasks.entries());
}

function getTasksCount() {
  return runningTasks.size;
}

module.exports = {
  runningTasks,
  killProcessTree,
  killProcessGracefully,
  addTask,
  getTask,
  deleteTask,
  getAllTasks,
  getTasksCount,
};