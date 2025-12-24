const fs = require('fs');
const path = require('path');

function getSystemMemoryUsage() {
  const used = process.memoryUsage();
  return {
    heapUsed: Math.round(used.heapUsed / 1024 / 1024),
    heapTotal: Math.round(used.heapTotal / 1024 / 1024),
    external: Math.round(used.external / 1024 / 1024),
    rss: Math.round(used.rss / 1024 / 1024),
  };
}

function ensureDir(dirPath) {
  if (!dirPath) return;
  if (!fs.existsSync(dirPath)) {
    fs.mkdirSync(dirPath, { recursive: true });
  }
}

function redactEnv(envObj) {
  const out = {};
  const sensitive = /pass|pwd|secret|token|key/i;
  for (const [k, v] of Object.entries(envObj || {})) {
    if (sensitive.test(k)) {
      out[k] = '[REDACTED]';
    } else {
      out[k] = v;
    }
  }
  return out;
}

function getLogFilePath() {
  const dir = process.env.LOG_DIR;
  if (!dir) return null;
  ensureDir(dir);
  return path.join(dir, 'backend.log.jsonl');
}

function writeJsonl(entry) {
  const filePath = getLogFilePath();
  if (!filePath) return;
  try {
    fs.appendFileSync(filePath, `${JSON.stringify(entry)}\n`);
  } catch (e) {
    // best effort
  }
}

function logEvent(level, message, meta = {}) {
  const entry = {
    ts: new Date().toISOString(),
    level,
    source: 'backend',
    message,
    ...meta,
  };
  // Human readable console line
  console.log(`[${entry.ts}] [${level}] ${message}`);
  // Structured line for support
  writeJsonl(entry);
}

function logWithMemory(message, level = 'INFO', meta = {}) {
  const memory = getSystemMemoryUsage();
  logEvent(level, message, {
    memory,
    ...meta,
  });
}

module.exports = {
  getSystemMemoryUsage,
  logWithMemory,
  logEvent,
  redactEnv,
};