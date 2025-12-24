const fs = require('fs');
const path = require('path');
const { logWithMemory } = require('./logger');

function ensureDir(d) {
  if (!fs.existsSync(d)) fs.mkdirSync(d, { recursive: true });
}

function safeResolve(baseDir, p) {
  const base = fs.realpathSync(baseDir);
  const target = path.resolve(base, p);

  // Symbolic link verification
  try {
    const lst = fs.lstatSync(target);
    if (lst.isSymbolicLink()) {
      logWithMemory(`Symbolic link access refused: ${target}`, 'WARN');
      throw new Error('Symbolic link refused');
    }
  } catch (e) {
    if (e.code !== 'ENOENT') {
      logWithMemory(`File verification error: ${e.message}`, 'WARN');
    }
  }

  const real = fs.existsSync(target) ? fs.realpathSync(target) : target;

  // Path traversal protection
  const normalizedBase = path.normalize(base + path.sep);
  const normalizedReal = path.normalize(real + path.sep);

  if (real !== base && !normalizedReal.startsWith(normalizedBase)) {
    logWithMemory(`Path traversal detected: ${real} outside ${base}`, 'ERROR');
    throw new Error('Path traversal detected');
  }

  return real;
}

function getDirSize(dirPath) {
  let size = 0;
  try {
    const files = fs.readdirSync(dirPath);
    files.forEach(file => {
      const p = path.join(dirPath, file);
      const st = fs.statSync(p);
      if (st.isFile()) size += st.size;
      else if (st.isDirectory()) size += getDirSize(p);
    });
  } catch {}
  return Math.round((size / 1024 / 1024) * 100) / 100;
}

module.exports = { ensureDir, safeResolve, getDirSize };