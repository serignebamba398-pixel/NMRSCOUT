const { app, BrowserWindow, ipcMain, dialog, Menu } = require('electron');
const path = require('path');
const { fork } = require('child_process');
const fs = require('fs');
const archiver = require('archiver');

let mainWindow;
let backendProcess;
let currentProjectPath = null;

// Determine if we are in development or production
const isDev = !app.isPackaged;

// Fichier de configuration pour stocker le dernier projet ouvert
const getConfigPath = () => path.join(app.getPath('userData'), 'config.json');

const loadConfig = () => {
  try {
    const configPath = getConfigPath();
    if (fs.existsSync(configPath)) {
      return JSON.parse(fs.readFileSync(configPath, 'utf8'));
    }
  } catch (e) {
    console.error('Erreur chargement config:', e);
  }
  return {};
};

const ensureDir = (dirPath) => {
  if (!fs.existsSync(dirPath)) {
    fs.mkdirSync(dirPath, { recursive: true });
  }
};

const tryEnforceLogRetention = (logsDir, maxBytes) => {
  try {
    const entries = fs
      .readdirSync(logsDir)
      .map((name) => {
        const fullPath = path.join(logsDir, name);
        let stat;
        try {
          stat = fs.statSync(fullPath);
        } catch {
          return null;
        }
        if (!stat.isFile()) return null;
        return { name, fullPath, size: stat.size, mtimeMs: stat.mtimeMs };
      })
      .filter(Boolean);

    let totalBytes = 0;
    for (const f of entries) totalBytes += f.size;
    if (totalBytes <= maxBytes) return;

    entries.sort((a, b) => a.mtimeMs - b.mtimeMs);
    for (const f of entries) {
      try {
        fs.unlinkSync(f.fullPath);
        totalBytes -= f.size;
      } catch {
        // ignore
      }
      if (totalBytes <= maxBytes) break;
    }
  } catch {
    // ignore
  }
};

// Logging helpers (stored under userData/logs)
const getLogsDir = () => {
  const dir = path.join(app.getPath('userData'), 'logs');
  ensureDir(dir);

  const maxMbRaw = process.env.LOG_RETENTION_MAX_MB;
  const parsed = Number(maxMbRaw);
  const maxMb = Number.isFinite(parsed) ? parsed : 25;
  tryEnforceLogRetention(dir, Math.max(1, maxMb) * 1024 * 1024);
  return dir;
};

const writeJsonl = (filePath, entry) => {
  try {
    ensureDir(path.dirname(filePath));
    fs.appendFileSync(filePath, `${JSON.stringify(entry)}\n`);
  } catch (e) {
    console.error('Log write error:', e);
  }
};

const mainLogPath = () => path.join(getLogsDir(), 'main.log.jsonl');
const backendLogPath = () => path.join(getLogsDir(), 'backend.log.txt');

const logMain = (level, message, meta = {}) => {
  writeJsonl(mainLogPath(), {
    ts: new Date().toISOString(),
    level,
    source: 'electron-main',
    message,
    ...meta,
  });
};

const appendTextLog = (filePath, prefix, text) => {
  try {
    ensureDir(path.dirname(filePath));
    fs.appendFileSync(filePath, `[${new Date().toISOString()}] ${prefix} ${text}\n`);
  } catch (e) {
    console.error('Log write error:', e);
  }
};

const saveConfig = (config) => {
  try {
    fs.writeFileSync(getConfigPath(), JSON.stringify(config, null, 2));
  } catch (e) {
    console.error('Erreur sauvegarde config:', e);
  }
};

function createWindow() {
  // Enlever la barre de menu (File, Edit, View, Window, Help)
  Menu.setApplicationMenu(null);
  
  mainWindow = new BrowserWindow({
    width: 1280,
    height: 800,
    webPreferences: {
      nodeIntegration: false,
      contextIsolation: true,
      preload: path.join(__dirname, 'preload.js'),
      webSecurity: false, // ⚠️ Required for SSE streaming from file:// to http://localhost
    },
  });

  if (isDev) {
    // In dev, load the React dev server
    mainWindow.loadURL('http://localhost:3000');
    mainWindow.webContents.openDevTools();
  } else {
    // In prod, load the built React app
    mainWindow.loadFile(path.join(__dirname, 'frontend', 'build', 'index.html'));
    // Ne pas ouvrir DevTools en production
    // mainWindow.webContents.openDevTools();
  }

  mainWindow.on('closed', () => {
    mainWindow = null;
  });
}

function startBackend(projectPath = null) {
  return new Promise((resolve) => {
    ensureDir(getLogsDir());
    logMain('INFO', 'Starting backend process', { projectPath, isDev });

    // Arrêter le backend existant si présent
    if (backendProcess) {
      backendProcess.kill();
      backendProcess = null;
    }

    const backendPath = path.join(__dirname, 'backend', 'server.js');
    
    // Set environment variables for the backend
    const env = { ...process.env };
    env.PORT = 5000; // Force port
    env.IS_ELECTRON = 'true';
    env.LOG_DIR = getLogsDir();
    env.LOG_LEVEL = env.LOG_LEVEL || 'INFO';
    
    if (!isDev) {
      env.RESOURCES_PATH = process.resourcesPath;
      
      // Check if bundled Python executable exists
      const exeName = process.platform === 'win32' ? 'nmr_engine_production.exe' : 'nmr_engine_production';
      // In package.json, we map backend/dist/nmr_engine_production to "python" folder in resources
      const bundledPythonPath = path.join(process.resourcesPath, 'python', exeName);
      
      if (fs.existsSync(bundledPythonPath)) {
        env.IS_ELECTRON_BUNDLED = 'true';
        env.BUNDLED_PYTHON_PATH = bundledPythonPath;
        console.log('[Main] Using bundled Python:', bundledPythonPath);
      } else {
        // No bundled Python - will use system Python
        env.IS_ELECTRON_BUNDLED = 'false';
        console.log('[Main] No bundled Python found, using system Python');
        console.log('[Main] Looked for:', bundledPythonPath);
      }
      
      // Path to bundled R
      const rExeName = process.platform === 'win32' ? 'bin/Rscript.exe' : 'bin/Rscript';
      const bundledRPath = path.join(process.resourcesPath, 'R-Portable', rExeName);
      if (fs.existsSync(bundledRPath)) {
        env.BUNDLED_R_PATH = bundledRPath;
        console.log('[Main] Using bundled R:', bundledRPath);
      }

      // Set Model Paths
      const modelsDir = path.join(process.resourcesPath, 'models');
      const model1H = path.join(modelsDir, '1H_sparsified_proposed_proposed_set2set_1.pt');
      const model13C = path.join(modelsDir, '13C_sparsified_proposed_proposed_set2set_1.pt');
      
      if (fs.existsSync(model1H)) {
        env.MODEL_1H_PATH = model1H;
        console.log('[Main] Found 1H Model:', model1H);
      }
      if (fs.existsSync(model13C)) {
        env.MODEL_13C_PATH = model13C;
        console.log('[Main] Found 13C Model:', model13C);
      }
    } else {
      // In Dev, try to find models in backend/models
      const devModelsDir = path.join(__dirname, 'backend', 'models');
      const model1H = path.join(devModelsDir, '1H_sparsified_proposed_proposed_set2set_1.pt');
      const model13C = path.join(devModelsDir, '13C_sparsified_proposed_proposed_set2set_1.pt');
      
      if (fs.existsSync(model1H)) env.MODEL_1H_PATH = model1H;
      if (fs.existsSync(model13C)) env.MODEL_13C_PATH = model13C;
    }

    // Define persistent paths for cache and projects
    const userDataPath = app.getPath('userData');
    env.CACHE_DIR = path.join(userDataPath, 'cache');
    env.PROJECTS_FILE = path.join(userDataPath, 'projects.json');
    
    // Ensure cache dir exists
    if (!fs.existsSync(env.CACHE_DIR)) {
      fs.mkdirSync(env.CACHE_DIR, { recursive: true });
    }

    // Si un projet est spécifié, utiliser son dossier pour uploads/outputs
    if (projectPath) {
      currentProjectPath = projectPath;
      env.PROJECT_PATH = projectPath;
      env.UPLOADS_DIR = path.join(projectPath, 'uploads');
      env.OUTPUTS_DIR = path.join(projectPath, 'outputs');
      env.NMR_DATA_DIR = path.join(projectPath, 'data');
      
      // Créer les dossiers du projet
      fs.mkdirSync(env.UPLOADS_DIR, { recursive: true });
      fs.mkdirSync(env.OUTPUTS_DIR, { recursive: true });
      fs.mkdirSync(env.NMR_DATA_DIR, { recursive: true });
    } else {
      // Mode sans projet - dossier temporaire
      const tempDir = path.join(app.getPath('temp'), 'nmr-annotation-temp');
      env.UPLOADS_DIR = path.join(tempDir, 'uploads');
      env.OUTPUTS_DIR = path.join(tempDir, 'outputs');
      fs.mkdirSync(env.UPLOADS_DIR, { recursive: true });
      fs.mkdirSync(env.OUTPUTS_DIR, { recursive: true });
    }

    backendProcess = fork(backendPath, [], {
      env,
      stdio: ['pipe', 'pipe', 'pipe', 'ipc'],
    });

    logMain('INFO', 'Backend started', { pid: backendProcess.pid, backendPath });
    appendTextLog(backendLogPath(), '[INFO]', `Backend started. PID: ${backendProcess.pid}`);
    appendTextLog(backendLogPath(), '[INFO]', `Backend path: ${backendPath}`);

    backendProcess.stdout.on('data', (data) => {
      const msg = data.toString();
      console.log(`[Backend]: ${msg}`);
      appendTextLog(backendLogPath(), '[STDOUT]', msg.trim());
    });

    backendProcess.stderr.on('data', (data) => {
      const msg = data.toString();
      console.error(`[Backend Error]: ${msg}`);
      appendTextLog(backendLogPath(), '[STDERR]', msg.trim());
    });

    // Timeout: si pas de "ready", on considère que le backend n'a pas démarré.
    const timeout = setTimeout(() => {
      console.error('[Main] Backend ready timeout');
      try {
        if (backendProcess) backendProcess.kill();
      } catch {}
      if (mainWindow) {
        mainWindow.webContents.send('backend-error', {
          projectPath: currentProjectPath,
          error: 'Backend did not become ready (timeout)',
        });
      }
      resolve({ ready: false, timeout: true, error: 'Backend did not become ready (timeout)' });
    }, 15000);

    // Communiquer avec le backend via IPC
    backendProcess.on('message', (msg) => {
      console.log('[Main] Message from backend:', msg);
      logMain('INFO', 'Message from backend', { msg });
      if (msg.type === 'ready') {
        clearTimeout(timeout);
        if (mainWindow) {
          mainWindow.webContents.send('backend-ready', { projectPath: currentProjectPath });
        }
        resolve({ ready: true, projectPath: currentProjectPath });
      }
    });

    backendProcess.on('error', (err) => {
      console.error('[Main] Backend error:', err);
      logMain('ERROR', 'Backend process error', { error: err?.message, stack: err?.stack });
      clearTimeout(timeout);
      if (mainWindow) {
        mainWindow.webContents.send('backend-error', {
          projectPath: currentProjectPath,
          error: err.message,
        });
      }
      resolve({ ready: false, error: err.message });
    });

    backendProcess.on('exit', (code, signal) => {
      console.error(`[Main] Backend exited with code ${code} and signal ${signal}`);
      logMain('WARN', 'Backend exited', { code, signal });
      clearTimeout(timeout);
      if (code !== 0) {
        if (mainWindow) {
          mainWindow.webContents.send('backend-error', {
            projectPath: currentProjectPath,
            error: `Backend exited with code ${code}`,
          });
        }
        resolve({ ready: false, error: `Backend exited with code ${code}` });
      }
    });
  });
}

// Export logs as a zip file for support
ipcMain.handle('export-logs', async (event, payload = {}) => {
  try {
    ensureDir(getLogsDir());

    const ts = new Date().toISOString().replace(/[:.]/g, '-');
    const defaultPath = path.join(app.getPath('documents'), `nmr-annotation-logs-${ts}.zip`);

    const { canceled, filePath } = await dialog.showSaveDialog(mainWindow, {
      title: 'Exporter les logs',
      defaultPath,
      filters: [{ name: 'Zip', extensions: ['zip'] }],
    });

    if (canceled || !filePath) {
      return { success: false, canceled: true };
    }

    await new Promise((resolve, reject) => {
      const output = fs.createWriteStream(filePath);
      const archive = archiver('zip', { zlib: { level: 9 } });

      output.on('close', resolve);
      output.on('error', reject);
      archive.on('error', reject);

      archive.pipe(output);

      // Include log directory files
      archive.directory(getLogsDir(), 'logs');

      // Include renderer logs (in-memory)
      if (payload && Array.isArray(payload.rendererLogs)) {
        archive.append(JSON.stringify(payload.rendererLogs, null, 2), { name: 'renderer-logs.json' });
      }

      // Include a small manifest
      const manifest = {
        exportedAt: new Date().toISOString(),
        appVersion: app.getVersion(),
        platform: process.platform,
        arch: process.arch,
        isPackaged: app.isPackaged,
        projectPath: currentProjectPath,
      };
      archive.append(JSON.stringify(manifest, null, 2), { name: 'manifest.json' });

      archive.finalize();
    });

    logMain('INFO', 'Logs exported', { filePath });
    return { success: true, filePath };
  } catch (e) {
    logMain('ERROR', 'Failed to export logs', { error: e?.message, stack: e?.stack });
    return { success: false, error: e?.message || String(e) };
  }
});

// IPC Handlers pour la gestion des projets

// Créer un nouveau projet (choisir un dossier vide ou nouveau)
ipcMain.handle('create-project', async () => {
  console.log('[Main] IPC create-project called');
  const result = await dialog.showOpenDialog(mainWindow, {
    properties: ['openDirectory', 'createDirectory'],
    title: 'Créer un nouveau projet - Sélectionner un dossier',
    buttonLabel: 'Créer le projet ici',
  });
  
  if (result.canceled || !result.filePaths[0]) return null;
  
  const projectDir = result.filePaths[0];
  
  // Créer la structure du projet
  const projectFile = path.join(projectDir, 'project.nmrproj');
  let projectData;
  if (fs.existsSync(projectFile)) {
    // Ne jamais écraser un projet existant à la création.
    // L'écrasement doit arriver uniquement via Save / Save As.
    projectData = JSON.parse(fs.readFileSync(projectFile, 'utf8'));
  } else {
    projectData = {
      name: path.basename(projectDir),
      createdAt: new Date().toISOString(),
      updatedAt: new Date().toISOString(),
      version: '1.0.0',
      state: {},
    };
    fs.writeFileSync(projectFile, JSON.stringify(projectData, null, 2));
  }
  fs.mkdirSync(path.join(projectDir, 'uploads'), { recursive: true });
  fs.mkdirSync(path.join(projectDir, 'outputs'), { recursive: true });
  fs.mkdirSync(path.join(projectDir, 'data'), { recursive: true });
  
  // Sauvegarder comme dernier projet
  saveConfig({ lastProject: projectDir });
  
  // Redémarrer le backend avec le nouveau projet et attendre qu'il soit prêt
  const r = await startBackend(projectDir);
  if (!r?.ready) {
    throw new Error(r?.error || 'Backend failed to start for new project');
  }
  
  return { path: projectDir, ...projectData };
});

// Ouvrir un projet existant (dossier contenant project.nmrproj)
ipcMain.handle('open-project', async () => {
  console.log('[Main] IPC open-project called');
  const result = await dialog.showOpenDialog(mainWindow, {
    properties: ['openDirectory'],
    title: 'Ouvrir un projet NMR',
    buttonLabel: 'Ouvrir le projet',
  });
  
  if (result.canceled || !result.filePaths[0]) return null;
  
  const projectDir = result.filePaths[0];
  const projectFile = path.join(projectDir, 'project.nmrproj');
  
  // Vérifier si c'est un projet valide ou créer un nouveau fichier projet
  let projectData;
  if (fs.existsSync(projectFile)) {
    projectData = JSON.parse(fs.readFileSync(projectFile, 'utf8'));
  } else {
    // Créer un fichier projet pour un dossier existant
    projectData = {
      name: path.basename(projectDir),
      createdAt: new Date().toISOString(),
      updatedAt: new Date().toISOString(),
      version: '1.0.0',
      state: {}
    };
    fs.writeFileSync(projectFile, JSON.stringify(projectData, null, 2));
  }
  
  // Créer les sous-dossiers si nécessaire
  fs.mkdirSync(path.join(projectDir, 'uploads'), { recursive: true });
  fs.mkdirSync(path.join(projectDir, 'outputs'), { recursive: true });
  fs.mkdirSync(path.join(projectDir, 'data'), { recursive: true });
  
  // Sauvegarder comme dernier projet
  saveConfig({ lastProject: projectDir });
  
  // Redémarrer le backend avec ce projet et attendre qu'il soit prêt
  const r = await startBackend(projectDir);
  if (!r?.ready) {
    throw new Error(r?.error || 'Backend failed to start for opened project');
  }
  
  return { path: projectDir, ...projectData };
});

// Sauvegarder l'état du projet
ipcMain.handle('save-project-state', async (event, state) => {
  if (!currentProjectPath) return { success: false, error: 'Aucun projet ouvert' };
  
  try {
    const projectFile = path.join(currentProjectPath, 'project.nmrproj');
    let projectData = {};
    
    if (fs.existsSync(projectFile)) {
      projectData = JSON.parse(fs.readFileSync(projectFile, 'utf8'));
    }
    
    projectData.state = state;
    projectData.updatedAt = new Date().toISOString();
    
    fs.writeFileSync(projectFile, JSON.stringify(projectData, null, 2));
    
    return { success: true };
  } catch (e) {
    return { success: false, error: e.message };
  }
});

// Obtenir le projet courant
ipcMain.handle('get-current-project', () => {
  if (!currentProjectPath) return null;
  
  try {
    const projectFile = path.join(currentProjectPath, 'project.nmrproj');
    if (fs.existsSync(projectFile)) {
      const projectData = JSON.parse(fs.readFileSync(projectFile, 'utf8'));
      return { path: currentProjectPath, ...projectData };
    }
  } catch (e) {
    console.error('Erreur lecture projet:', e);
  }
  return null;
});

// Obtenir le dernier projet ouvert
ipcMain.handle('get-last-project', () => {
  const config = loadConfig();
  if (config.lastProject && fs.existsSync(config.lastProject)) {
    return config.lastProject;
  }
  return null;
});

// Fermer le projet courant
ipcMain.handle('close-project', async () => {
  currentProjectPath = null;
  const r = await startBackend(null); // Redémarrer sans projet
  if (!r?.ready) {
    throw new Error(r?.error || 'Backend failed to start after closing project');
  }
  return { success: true };
});

// Sauvegarder sous (Save As) - copier vers un nouveau dossier
ipcMain.handle('save-project-as', async (event, state) => {
  const result = await dialog.showOpenDialog(mainWindow, {
    properties: ['openDirectory', 'createDirectory'],
    title: 'Sauvegarder le projet sous - Sélectionner un dossier',
    buttonLabel: 'Sauvegarder ici',
  });
  
  if (result.canceled || !result.filePaths[0]) return null;
  
  const newProjectDir = result.filePaths[0];
  
  try {
    // Copier les fichiers du projet actuel si existant
    if (currentProjectPath && fs.existsSync(currentProjectPath)) {
      const copyRecursive = (src, dest) => {
        if (!fs.existsSync(src)) return;
        const stat = fs.statSync(src);
        if (stat.isDirectory()) {
          fs.mkdirSync(dest, { recursive: true });
          fs.readdirSync(src).forEach(child => {
            copyRecursive(path.join(src, child), path.join(dest, child));
          });
        } else {
          fs.copyFileSync(src, dest);
        }
      };
      
      // Copier uploads, outputs, data
      ['uploads', 'outputs', 'data'].forEach(dir => {
        const srcDir = path.join(currentProjectPath, dir);
        const destDir = path.join(newProjectDir, dir);
        if (fs.existsSync(srcDir)) {
          copyRecursive(srcDir, destDir);
        } else {
          fs.mkdirSync(destDir, { recursive: true });
        }
      });
    } else {
      // Créer les dossiers vides
      fs.mkdirSync(path.join(newProjectDir, 'uploads'), { recursive: true });
      fs.mkdirSync(path.join(newProjectDir, 'outputs'), { recursive: true });
      fs.mkdirSync(path.join(newProjectDir, 'data'), { recursive: true });
    }
    
    // Créer le fichier projet
    const projectFile = path.join(newProjectDir, 'project.nmrproj');
    const projectData = {
      name: path.basename(newProjectDir),
      createdAt: new Date().toISOString(),
      updatedAt: new Date().toISOString(),
      version: '1.0.0',
      state: state || {}
    };
    
    fs.writeFileSync(projectFile, JSON.stringify(projectData, null, 2));
    
    // Sauvegarder comme dernier projet
    saveConfig({ lastProject: newProjectDir });
    
    // Redémarrer le backend avec le nouveau projet et attendre qu'il soit prêt
    const r = await startBackend(newProjectDir);
    if (!r?.ready) {
      throw new Error(r?.error || 'Backend failed to start after Save As');
    }
    
    return { path: newProjectDir, ...projectData };
  } catch (e) {
    console.error('Erreur Save As:', e);
    return null;
  }
});

ipcMain.handle('select-directory', async () => {
  const result = await dialog.showOpenDialog(mainWindow, {
    properties: ['openDirectory', 'createDirectory'],
    title: 'Sélectionner un dossier',
  });
  return result.canceled ? null : result.filePaths[0];
});

ipcMain.handle('get-user-data-path', () => {
  return app.getPath('userData');
});

ipcMain.handle('get-documents-path', () => {
  return app.getPath('documents');
});

app.on('ready', () => {
  // Charger le dernier projet si disponible
  const config = loadConfig();
  const lastProject = config.lastProject && fs.existsSync(config.lastProject) ? config.lastProject : null;
  
  startBackend(lastProject);
  createWindow();
});

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app.on('activate', () => {
  if (mainWindow === null) {
    createWindow();
  }
});

app.on('before-quit', () => {
  if (backendProcess) {
    backendProcess.kill();
  }
});
