require('dotenv').config();

const PORT = Number(process.env.PORT || 4000);
const MEMORY_LIMIT_MB = Number(process.env.MEMORY_LIMIT_MB || 6144);
const BODY_LIMIT_MB = Number(process.env.BODY_LIMIT_MB || 500);
const MAX_UPLOAD_MB = Number(process.env.MAX_UPLOAD_MB || 100);

// Default to false for task persistence - tasks continue running in background even if browser closes.
// Cancel button uses the /api/tasks/:id/cancel endpoint to explicitly stop tasks.
// Set KILL_ON_CLIENT_CLOSE=true if you want the old behavior (SSE close kills task).
const KILL_ON_CLIENT_CLOSE = String(process.env.KILL_ON_CLIENT_CLOSE || 'false').toLowerCase() === 'true';

// Enable task persistence - saves checkpoints to disk for resume after restart
const ENABLE_TASK_PERSISTENCE = String(process.env.ENABLE_TASK_PERSISTENCE || 'true').toLowerCase() === 'true';
const ENABLE_PROCESS_MONITORING = String(process.env.ENABLE_PROCESS_MONITORING || 'true').toLowerCase() === 'true';
const PROCESS_TIMEOUT_HOURS = Number(process.env.PROCESS_TIMEOUT_HOURS || 168); // 1 week default
const ENABLE_CHECKPOINTS = String(process.env.ENABLE_CHECKPOINTS || 'true').toLowerCase() === 'true';
const TASK_TTL_HOURS = Number(process.env.TASK_TTL_HOURS || 48);
const AUTO_CLEANUP_INTERVAL_HOURS = Number(process.env.AUTO_CLEANUP_INTERVAL_HOURS || TASK_TTL_HOURS);

const DEFAULT_PERSIST_KEY = process.env.REACT_APP_PERSIST_KEY || 'frontend-session';

const isWin = process.platform === 'win32';

const DEFAULT_LIGHT_PARAMS = {
  nb_cores: 2,
  batch_size: 100,
  memory_cleanup_interval: 25,
  enable_checkpoints: true,
  max_memory_mb: MEMORY_LIMIT_MB,
};

module.exports = {
  PORT,
  MEMORY_LIMIT_MB,
  BODY_LIMIT_MB,
  MAX_UPLOAD_MB,
  KILL_ON_CLIENT_CLOSE,
  ENABLE_PROCESS_MONITORING,
  PROCESS_TIMEOUT_HOURS,
  ENABLE_CHECKPOINTS,
  ENABLE_TASK_PERSISTENCE,
  TASK_TTL_HOURS,
  AUTO_CLEANUP_INTERVAL_HOURS,
  DEFAULT_PERSIST_KEY,
  DEFAULT_LIGHT_PARAMS,
  isWin,
};