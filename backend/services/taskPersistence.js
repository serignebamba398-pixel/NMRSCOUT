/**
 * Task Persistence Service
 * 
 * Permet de sauvegarder l'état des tâches sur disque pour
 * qu'elles puissent être reprises après un redémarrage du serveur
 * ou une fermeture de l'ordinateur.
 */

const fs = require('fs');
const path = require('path');
const { outputsDir, getActiveOutputsDir } = require('../config/paths');
const { logWithMemory } = require('../utils/logger');
const { ENABLE_CHECKPOINTS } = require('../config');

const CHECKPOINT_FILE = 'checkpoint.json';
const PENDING_TASKS_FILE = path.join(outputsDir, 'pending_tasks.json');

/**
 * États possibles d'une tâche persistée
 */
const TaskStatus = {
  PENDING: 'pending',       // En attente de reprise
  RUNNING: 'running',       // En cours d'exécution
  COMPLETED: 'completed',   // Terminée avec succès
  FAILED: 'failed',         // Échouée
  CANCELLED: 'cancelled',   // Annulée par l'utilisateur
  INTERRUPTED: 'interrupted' // Interrompue (fermeture ordi, etc.)
};

/**
 * Sauvegarde un checkpoint pour une tâche
 */
function saveCheckpoint(taskId, checkpointData) {
  if (!ENABLE_CHECKPOINTS) return false;
  
  try {
    const workDir = path.join(getActiveOutputsDir(), taskId);
    if (!fs.existsSync(workDir)) {
      fs.mkdirSync(workDir, { recursive: true });
    }
    
    const checkpointPath = path.join(workDir, CHECKPOINT_FILE);
    const checkpoint = {
      taskId,
      timestamp: new Date().toISOString(),
      version: 1,
      ...checkpointData
    };
    
    fs.writeFileSync(checkpointPath, JSON.stringify(checkpoint, null, 2), 'utf8');
    logWithMemory(`[${taskId}] Checkpoint saved: ${checkpointData.progress?.current || 0}/${checkpointData.progress?.total || '?'}`, 'DEBUG');
    return true;
  } catch (e) {
    logWithMemory(`[${taskId}] Failed to save checkpoint: ${e.message}`, 'WARN');
    return false;
  }
}

/**
 * Charge un checkpoint pour une tâche
 */
function loadCheckpoint(taskId) {
  try {
    const workDir = path.join(getActiveOutputsDir(), taskId);
    const checkpointPath = path.join(workDir, CHECKPOINT_FILE);
    
    if (!fs.existsSync(checkpointPath)) {
      return null;
    }
    
    const data = fs.readFileSync(checkpointPath, 'utf8');
    return JSON.parse(data);
  } catch (e) {
    logWithMemory(`[${taskId}] Failed to load checkpoint: ${e.message}`, 'WARN');
    return null;
  }
}

/**
 * Marque une tâche comme interrompue (pour reprise ultérieure)
 */
function markTaskInterrupted(taskId, reason = 'connection_lost') {
  try {
    const checkpoint = loadCheckpoint(taskId);
    if (!checkpoint) return false;
    
    checkpoint.status = TaskStatus.INTERRUPTED;
    checkpoint.interruptedAt = new Date().toISOString();
    checkpoint.interruptReason = reason;
    
    saveCheckpoint(taskId, checkpoint);
    addToPendingTasks(taskId, checkpoint);
    return true;
  } catch (e) {
    logWithMemory(`[${taskId}] Failed to mark as interrupted: ${e.message}`, 'WARN');
    return false;
  }
}

/**
 * Marque une tâche comme complétée
 */
function markTaskCompleted(taskId, outputs = {}) {
  try {
    const workDir = path.join(getActiveOutputsDir(), taskId);
    const checkpointPath = path.join(workDir, CHECKPOINT_FILE);
    
    if (fs.existsSync(checkpointPath)) {
      const checkpoint = loadCheckpoint(taskId);
      if (checkpoint) {
        checkpoint.status = TaskStatus.COMPLETED;
        checkpoint.completedAt = new Date().toISOString();
        checkpoint.outputs = outputs;
        fs.writeFileSync(checkpointPath, JSON.stringify(checkpoint, null, 2), 'utf8');
      }
    }
    
    removeFromPendingTasks(taskId);
    return true;
  } catch (e) {
    logWithMemory(`[${taskId}] Failed to mark as completed: ${e.message}`, 'WARN');
    return false;
  }
}

/**
 * Marque une tâche comme échouée
 */
function markTaskFailed(taskId, error) {
  try {
    const checkpoint = loadCheckpoint(taskId);
    if (checkpoint) {
      checkpoint.status = TaskStatus.FAILED;
      checkpoint.failedAt = new Date().toISOString();
      checkpoint.error = error;
      saveCheckpoint(taskId, checkpoint);
    }
    
    removeFromPendingTasks(taskId);
    return true;
  } catch (e) {
    logWithMemory(`[${taskId}] Failed to mark as failed: ${e.message}`, 'WARN');
    return false;
  }
}

/**
 * Ajoute une tâche à la liste des tâches en attente
 */
function addToPendingTasks(taskId, checkpointData) {
  try {
    const pending = loadPendingTasks();
    pending[taskId] = {
      taskId,
      type: checkpointData.type || 'unknown',
      addedAt: new Date().toISOString(),
      progress: checkpointData.progress,
      options: checkpointData.options,
      workDir: checkpointData.workDir,
    };
    
    savePendingTasks(pending);
    logWithMemory(`[${taskId}] Added to pending tasks`, 'INFO');
    return true;
  } catch (e) {
    logWithMemory(`Failed to add pending task: ${e.message}`, 'WARN');
    return false;
  }
}

/**
 * Retire une tâche de la liste des tâches en attente
 */
function removeFromPendingTasks(taskId) {
  try {
    const pending = loadPendingTasks();
    if (pending[taskId]) {
      delete pending[taskId];
      savePendingTasks(pending);
      logWithMemory(`[${taskId}] Removed from pending tasks`, 'DEBUG');
    }
    return true;
  } catch (e) {
    logWithMemory(`Failed to remove pending task: ${e.message}`, 'WARN');
    return false;
  }
}

/**
 * Charge la liste des tâches en attente
 */
function loadPendingTasks() {
  try {
    if (!fs.existsSync(PENDING_TASKS_FILE)) {
      return {};
    }
    const data = fs.readFileSync(PENDING_TASKS_FILE, 'utf8');
    return JSON.parse(data);
  } catch (e) {
    logWithMemory(`Failed to load pending tasks: ${e.message}`, 'WARN');
    return {};
  }
}

/**
 * Sauvegarde la liste des tâches en attente
 */
function savePendingTasks(pending) {
  try {
    // Ensure outputsDir exists
    if (!fs.existsSync(outputsDir)) {
      fs.mkdirSync(outputsDir, { recursive: true });
    }
    fs.writeFileSync(PENDING_TASKS_FILE, JSON.stringify(pending, null, 2), 'utf8');
  } catch (e) {
    logWithMemory(`Failed to save pending tasks: ${e.message}`, 'WARN');
  }
}

/**
 * Obtient toutes les tâches qui peuvent être reprises
 */
function getResumableTasks() {
  const pending = loadPendingTasks();
  const resumable = [];
  
  for (const [taskId, taskInfo] of Object.entries(pending)) {
    const checkpoint = loadCheckpoint(taskId);
    if (checkpoint && checkpoint.status !== TaskStatus.COMPLETED && checkpoint.status !== TaskStatus.CANCELLED) {
      resumable.push({
        taskId,
        type: taskInfo.type,
        progress: checkpoint.progress,
        interruptedAt: checkpoint.interruptedAt || taskInfo.addedAt,
        canResume: Boolean(checkpoint.resumeData),
        options: checkpoint.options,
      });
    }
  }
  
  return resumable;
}

/**
 * Vérifie si une tâche peut être reprise
 */
function canResumeTask(taskId) {
  const checkpoint = loadCheckpoint(taskId);
  if (!checkpoint) return false;
  
  return checkpoint.status === TaskStatus.INTERRUPTED && 
         checkpoint.resumeData && 
         Object.keys(checkpoint.resumeData).length > 0;
}

/**
 * Obtient les données de reprise pour une tâche
 */
function getResumeData(taskId) {
  const checkpoint = loadCheckpoint(taskId);
  if (!checkpoint) return null;
  
  return {
    taskId,
    type: checkpoint.type,
    progress: checkpoint.progress,
    resumeData: checkpoint.resumeData,
    options: checkpoint.options,
    workDir: checkpoint.workDir,
  };
}

/**
 * Nettoie les vieilles tâches pendantes (plus de X heures)
 */
function cleanupOldPendingTasks(maxAgeHours = 48) {
  const pending = loadPendingTasks();
  const now = Date.now();
  const maxAge = maxAgeHours * 60 * 60 * 1000;
  let cleaned = 0;
  
  for (const [taskId, taskInfo] of Object.entries(pending)) {
    const addedAt = new Date(taskInfo.addedAt).getTime();
    if (now - addedAt > maxAge) {
      delete pending[taskId];
      cleaned++;
      logWithMemory(`[${taskId}] Removed old pending task (${maxAgeHours}h+)`, 'INFO');
    }
  }
  
  if (cleaned > 0) {
    savePendingTasks(pending);
    logWithMemory(`Cleaned ${cleaned} old pending tasks`, 'INFO');
  }
  
  return cleaned;
}

/**
 * Initialisation au démarrage du serveur
 * Marque toutes les tâches "running" comme "interrupted"
 */
function initializeOnStartup() {
  logWithMemory('Task persistence: checking for interrupted tasks...', 'INFO');
  
  try {
    // Check outputs directory for checkpoint files
    if (!fs.existsSync(outputsDir)) {
      return;
    }
    
    const dirs = fs.readdirSync(outputsDir).filter(d => {
      const p = path.join(outputsDir, d);
      return fs.statSync(p).isDirectory();
    });
    
    let interrupted = 0;
    for (const dir of dirs) {
      const checkpointPath = path.join(outputsDir, dir, CHECKPOINT_FILE);
      if (fs.existsSync(checkpointPath)) {
        try {
          const checkpoint = JSON.parse(fs.readFileSync(checkpointPath, 'utf8'));
          if (checkpoint.status === TaskStatus.RUNNING) {
            // Task was running when server stopped - mark as interrupted
            checkpoint.status = TaskStatus.INTERRUPTED;
            checkpoint.interruptedAt = new Date().toISOString();
            checkpoint.interruptReason = 'server_restart';
            fs.writeFileSync(checkpointPath, JSON.stringify(checkpoint, null, 2), 'utf8');
            addToPendingTasks(dir, checkpoint);
            interrupted++;
            logWithMemory(`[${dir}] Marked as interrupted (was running)`, 'INFO');
          }
        } catch (e) {
          logWithMemory(`[${dir}] Failed to process checkpoint: ${e.message}`, 'WARN');
        }
      }
    }
    
    if (interrupted > 0) {
      logWithMemory(`Found ${interrupted} interrupted tasks that can be resumed`, 'INFO');
    }
    
    // Cleanup old tasks
    cleanupOldPendingTasks();
    
  } catch (e) {
    logWithMemory(`Task persistence init error: ${e.message}`, 'ERROR');
  }
}

module.exports = {
  TaskStatus,
  saveCheckpoint,
  loadCheckpoint,
  markTaskInterrupted,
  markTaskCompleted,
  markTaskFailed,
  addToPendingTasks,
  removeFromPendingTasks,
  loadPendingTasks,
  getResumableTasks,
  canResumeTask,
  getResumeData,
  cleanupOldPendingTasks,
  initializeOnStartup,
};
