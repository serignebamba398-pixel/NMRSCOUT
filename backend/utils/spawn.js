const { spawn, exec } = require('child_process');
const path = require('path');
const fs = require('fs');
const { MEMORY_LIMIT_MB } = require('../config');
const { logWithMemory } = require('./logger');

function preferWindowsX64Rscript(rscriptPath) {
  if (!rscriptPath || process.platform !== 'win32') return rscriptPath;
  if (process.arch !== 'x64') return rscriptPath;

  const normalized = path.normalize(rscriptPath);
  const lower = normalized.toLowerCase();
  if (!lower.endsWith(`${path.sep}rscript.exe`)) return rscriptPath;

  // If caller points at ...\bin\Rscript.exe and ...\bin\x64\Rscript.exe exists,
  // prefer x64 build to avoid 32/64-bit runtime mismatches.
  const parts = lower.split(path.sep);
  const binIdx = parts.lastIndexOf('bin');
  if (binIdx === -1) return rscriptPath;

  const afterBin = parts[binIdx + 1];
  if (afterBin === 'x64') return rscriptPath;

  const candidate = path.join(path.dirname(normalized), 'x64', 'Rscript.exe');
  if (fs.existsSync(candidate)) return candidate;
  return rscriptPath;
}

function computeRPortableRootFromRscript(rscriptPath) {
  if (!rscriptPath) return null;
  try {
    const p = path.normalize(rscriptPath);
    // ...\R-Portable\bin\Rscript.exe  -> root = ...\R-Portable
    // ...\R-Portable\bin\x64\Rscript.exe -> root = ...\R-Portable
    const dir = path.dirname(p);
    const base = path.basename(dir).toLowerCase();
    if (base === 'x64') return path.resolve(dir, '..', '..');
    if (base === 'bin') return path.resolve(dir, '..');
    // Best effort: walk up a little if needed
    const up1 = path.resolve(dir, '..');
    const up2 = path.resolve(dir, '..', '..');
    if (fs.existsSync(path.join(up1, 'library'))) return up1;
    if (fs.existsSync(path.join(up2, 'library'))) return up2;
  } catch {}
  return null;
}

function getRscript() {
  // PRIORITY 1: BUNDLED_R_PATH explicitly set (Electron production)
  if (process.env.BUNDLED_R_PATH && fs.existsSync(process.env.BUNDLED_R_PATH)) {
    const chosen = preferWindowsX64Rscript(process.env.BUNDLED_R_PATH);
    console.log('[getRscript] Using BUNDLED_R_PATH:', chosen);
    return chosen;
  }
  
  // PRIORITY 2: Check RESOURCES_PATH for R-Portable (Electron packaged)
  if (process.env.RESOURCES_PATH) {
    if (process.platform === 'win32') {
      const rPathX64 = path.join(process.env.RESOURCES_PATH, 'R-Portable', 'bin', 'x64', 'Rscript.exe');
      const rPath = path.join(process.env.RESOURCES_PATH, 'R-Portable', 'bin', 'Rscript.exe');
      if (fs.existsSync(rPathX64)) {
        console.log('[getRscript] Found R-Portable (x64) at:', rPathX64);
        return rPathX64;
      }
      if (fs.existsSync(rPath)) {
        console.log('[getRscript] Found R-Portable at:', rPath);
        return rPath;
      }
    } else {
      const rPath = path.join(process.env.RESOURCES_PATH, 'R-Portable', 'bin', 'Rscript');
      if (fs.existsSync(rPath)) {
        console.log('[getRscript] Found R-Portable at:', rPath);
        return rPath;
      }
    }
  }

  // PRIORITY 3: Dev mode - look for R-Portable in project root
  const devRPortable = path.join(__dirname, '..', '..', 'R-Portable');
  if (fs.existsSync(devRPortable)) {
    if (process.platform === 'win32') {
      const rPathX64 = path.join(devRPortable, 'bin', 'x64', 'Rscript.exe');
      const rPath = path.join(devRPortable, 'bin', 'Rscript.exe');
      if (fs.existsSync(rPathX64)) {
        console.log('[getRscript] Found dev R-Portable (x64) at:', rPathX64);
        return rPathX64;
      }
      if (fs.existsSync(rPath)) {
        console.log('[getRscript] Found dev R-Portable at:', rPath);
        return rPath;
      }
    } else {
      const rPath = path.join(devRPortable, 'bin', 'Rscript');
      if (fs.existsSync(rPath)) {
        console.log('[getRscript] Found dev R-Portable at:', rPath);
        return rPath;
      }
    }
  }

  // PRIORITY 4: Docker/server mode - use R_PATH env var (system R)
  if (process.env.R_PATH && fs.existsSync(process.env.R_PATH)) {
    console.log('[getRscript] Using R_PATH (Docker/server mode):', process.env.R_PATH);
    return process.env.R_PATH;
  }

  // NO FALLBACK - App must have R available
  console.error('[getRscript] ERROR: Rscript not found!');
  console.error('[getRscript] Checked BUNDLED_R_PATH:', process.env.BUNDLED_R_PATH);
  console.error('[getRscript] Checked RESOURCES_PATH:', process.env.RESOURCES_PATH);
  console.error('[getRscript] Checked R_PATH:', process.env.R_PATH);
  console.error('[getRscript] Checked dev path:', devRPortable);
  throw new Error('Rscript not found. Set R_PATH or provide bundled R-Portable.');
}

function buildREnv(memoryLimitMB, rscriptPath = null) {
  const env = { ...process.env };
  const lim = Number.isFinite(+memoryLimitMB) ? +memoryLimitMB : MEMORY_LIMIT_MB;
  const rMaxMB = Math.max(256, Math.floor(lim * 0.8));
  env.R_MAX_VSIZE = `${rMaxMB}M`;

  // On Windows, R can crash (0xC0000005) if it picks up incompatible DLLs
  // from a different toolchain earlier in PATH. When using R-Portable,
  // force a clean-ish environment.
  if (process.platform === 'win32' && rscriptPath) {
    const root = computeRPortableRootFromRscript(rscriptPath);
    if (root && fs.existsSync(root)) {
      try {
        env.R_HOME = root;

        // Prevent mixing user/system libraries with bundled portable R.
        // Mismatched binary packages (e.g. data.table/jsonlite) can crash R
        // with 0xC0000005 before our scripts can catch an error.
        const portableLib = path.join(root, 'library');
        if (fs.existsSync(portableLib)) {
          env.R_LIBS_USER = portableLib;
          env.R_LIBS_SITE = portableLib;
          env.R_LIBS = portableLib;
        }

        // Disable user-level profile/env which can alter library paths.
        env.R_PROFILE_USER = 'NUL';
        env.R_ENVIRON_USER = 'NUL';

        // CRITICAL: Remove ANY other R installation from PATH to prevent DLL conflicts.
        // This is the main cause of 0xC0000005 crashes when user has R system installed.
        const currentPath = env.PATH || process.env.PATH || '';
        const cleanedPathParts = currentPath.split(path.delimiter).filter((p) => {
          const lower = p.toLowerCase();
          // Remove paths containing other R installations
          if (lower.includes('\\r\\') || lower.includes('/r/')) return false;
          if (lower.includes('\\r-') || lower.includes('/r-')) return false;
          if (lower.match(/[\\\/]r[\\\/]r-\d/i)) return false;
          if (lower.includes('rtools')) return false;
          // Keep our own R-Portable path (will be added at front)
          if (lower.includes(root.toLowerCase())) return false;
          return true;
        });

        const bin = path.join(root, 'bin');
        const binX64 = path.join(root, 'bin', 'x64');
        const rPaths = [];
        if (process.arch === 'x64' && fs.existsSync(binX64)) rPaths.push(binX64);
        if (fs.existsSync(bin)) rPaths.push(bin);
        
        // Build clean PATH: R-Portable first, then cleaned system PATH
        env.PATH = [...rPaths, ...cleanedPathParts].join(path.delimiter);
        
        console.log('[buildREnv] R_HOME:', env.R_HOME);
        console.log('[buildREnv] R_LIBS_USER:', env.R_LIBS_USER);
        console.log('[buildREnv] PATH starts with:', env.PATH.substring(0, 200));
      } catch (e) {
        console.error('[buildREnv] Error setting up R env:', e);
      }
    }
  }
  return env;
}

function getPythonExe(forML = false) {
  const isBundled = process.env.IS_ELECTRON_BUNDLED === 'true';
  const exeNames =
    process.platform === 'win32'
      ? ['nmr_engine_production.exe', 'nmr_engine.exe']
      : ['nmr_engine_production', 'nmr_engine'];
  
  // [DEBUG] Log Python resolution
  console.log('[getPythonExe] IS_ELECTRON_BUNDLED:', isBundled, 'forML:', forML);
  console.log('[getPythonExe] BUNDLED_PYTHON_PATH:', process.env.BUNDLED_PYTHON_PATH);
  
  // PRIORITY 1: BUNDLED_PYTHON_PATH explicitly set (Electron production)
  if (process.env.BUNDLED_PYTHON_PATH && fs.existsSync(process.env.BUNDLED_PYTHON_PATH)) {
    console.log('[getPythonExe] Using BUNDLED_PYTHON_PATH:', process.env.BUNDLED_PYTHON_PATH);
    return process.env.BUNDLED_PYTHON_PATH;
  }
  
  // PRIORITY 2: Check RESOURCES_PATH for bundled Python (Electron packaged)
  if (process.env.RESOURCES_PATH) {
    for (const exeName of exeNames) {
      // Try python/nmr_engine_production/nmr_engine_production.exe (new structure)
      const candidate1 = path.join(process.env.RESOURCES_PATH, 'python', 'nmr_engine_production', exeName);
      if (fs.existsSync(candidate1)) {
        console.log('[getPythonExe] Found bundled Python at:', candidate1);
        return candidate1;
      }
      // Fallback to python/nmr_engine_production.exe (old structure)
      const candidate2 = path.join(process.env.RESOURCES_PATH, 'python', exeName);
      if (fs.existsSync(candidate2)) {
        console.log('[getPythonExe] Found bundled Python at:', candidate2);
        return candidate2;
      }
    }
  }
  
  // PRIORITY 3: Dev mode - look for bundled Python in backend/dist
  for (const exeName of exeNames) {
    const devBundledPython = path.join(__dirname, '..', 'dist', 'nmr_engine_production', exeName);
    if (fs.existsSync(devBundledPython)) {
      console.log('[getPythonExe] Found dev bundled Python at:', devBundledPython);
      return devBundledPython;
    }
  }

  // PRIORITY 4: Dev fallback - prebuilt binary under repo resources/python
  for (const exeName of exeNames) {
    const devPrebuiltPython = path.join(__dirname, '..', '..', 'resources', 'python', exeName);
    if (fs.existsSync(devPrebuiltPython)) {
      console.log('[getPythonExe] Found dev prebuilt Python at:', devPrebuiltPython);
      return devPrebuiltPython;
    }
  }
  
  // PRIORITY 5: Docker/server mode - use PRED_PYTHON env var (conda python)
  if (process.env.PRED_PYTHON && fs.existsSync(process.env.PRED_PYTHON)) {
    console.log('[getPythonExe] Using PRED_PYTHON (Docker/server mode):', process.env.PRED_PYTHON);
    return process.env.PRED_PYTHON;
  }

  // NO FALLBACK TO SYSTEM PYTHON - App must be autonomous
  console.error('[getPythonExe] ERROR: Python not found!');
  console.error('[getPythonExe] Checked BUNDLED_PYTHON_PATH:', process.env.BUNDLED_PYTHON_PATH);
  console.error('[getPythonExe] Checked RESOURCES_PATH:', process.env.RESOURCES_PATH);
  console.error('[getPythonExe] Checked PRED_PYTHON:', process.env.PRED_PYTHON);
  console.error('[getPythonExe] Checked dev dist:', path.join(__dirname, '..', 'dist', 'nmr_engine_production'));
  throw new Error('Python not found. Set PRED_PYTHON or provide bundled nmr_engine.');
}

function detectDatasetDir() {
  const candidates = [
    process.env.NMR_DATA_DIR,
    path.join(__dirname, '..', 'ml_dependencies', 'data'),
    path.join(__dirname, '..', '..', 'ml_dependencies', 'data'),
  ].filter(Boolean);

  for (const d of candidates) {
    try {
      if (fs.existsSync(d)) return d;
    } catch {}
  }
  return null;
}

function buildPythonEnv() {
  const env = { ...process.env };

  const envRoot = process.env.PRED_ENV_ROOT;
  if (envRoot && fs.existsSync(envRoot)) {
    const scripts = path.join(envRoot, 'Scripts');
    const bin = path.join(envRoot, 'bin');
    const parts = [];

    if (process.platform === 'win32') {
      if (fs.existsSync(scripts)) parts.push(scripts);
      parts.push(envRoot);
    } else {
      if (fs.existsSync(bin)) parts.push(bin);
      parts.push(envRoot);
    }

    env.PATH = `${parts.join(path.delimiter)}${path.delimiter}${process.env.PATH || ''}`;
    env.CONDA_PREFIX = envRoot;
  }

  env.PYTHONIOENCODING = 'utf-8';
  env.PYTHONUTF8 = '1';
  env.DGLBACKEND = env.DGLBACKEND || 'pytorch';

  const datasetDir = detectDatasetDir();
  if (datasetDir) env.NMR_DATA_DIR = datasetDir;

  try {
    const scriptsDir = path.join(__dirname, '..', 'scripts');
    const scriptsAbs = path.resolve(scriptsDir);
    env.PYTHONPATH = env.PYTHONPATH
      ? `${scriptsAbs}${path.delimiter}${env.PYTHONPATH}`
      : scriptsAbs;
    
    // Allow hot-patching bundled scripts without recompilation
    // NMR_SCRIPTS_OVERRIDE tells the bundled exe to load scripts from external path
    if (fs.existsSync(scriptsAbs)) {
      env.NMR_SCRIPTS_OVERRIDE = scriptsAbs;
      console.log('[buildPythonEnv] NMR_SCRIPTS_OVERRIDE:', scriptsAbs);
    }
  } catch {}

  return env;
}

module.exports = {
  getRscript,
  buildREnv,
  getPythonExe,
  detectDatasetDir,
  buildPythonEnv,
};
