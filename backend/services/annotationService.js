const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const os = require('os');
const { v4: uuidv4 } = require('uuid');
const { Semaphore } = require('async-mutex');

const { DEFAULT_LIGHT_PARAMS, MEMORY_LIMIT_MB, PROCESS_TIMEOUT_HOURS, KILL_ON_CLIENT_CLOSE, ENABLE_TASK_PERSISTENCE } = require('../config');
const { scriptsDir, outputsDir, uploadsDir, scriptMap, legacyOut, getActiveOutputsDir } = require('../config/paths');
const { ensureDir, safeResolve } = require('../utils/filesystem');
const { getRscript, buildREnv } = require('../utils/spawn');
const { logWithMemory } = require('../utils/logger');
const { sseInit, sseSend } = require('../middleware/sse');
const { addTask, deleteTask, killProcessGracefully } = require('./processManager');
const projectService = require('./projectService');
const { saveCheckpoint, markTaskCompleted, markTaskFailed, TaskStatus } = require('./taskPersistence');
const { startPreventSleep, stopPreventSleep } = require('../utils/sleepPrevention');

const ANNOTATION_SEM = new Semaphore(2);

/**
 * Generate tolerance string for filename
 * Examples:
 *   1H/13C: "tol_0.05" or "tol_0.5"
 *   COSY/HSQC: "tol_f1_0.05_f2_0.1"
 */
function getToleranceStringForFilename(annotationType, params) {
  if (annotationType === '1H' || annotationType === '13C') {
    const tol = params.tol_correspondence || params.tol_correspondance || 0;
    return `tol_${tol}`.replace(/\./g, '_');
  } else if (annotationType === 'COSY' || annotationType === 'HSQC') {
    const tolF1 = params.tol_correspondence_f1 || params.tol_correspondance_F1 || 0;
    const tolF2 = params.tol_correspondence_f2 || params.tol_correspondance_F2 || 0;
    return `tol_f1_${tolF1}_f2_${tolF2}`.replace(/\./g, '_');
  }
  return 'tol_unknown';
}

async function handleAnnotation(req, res, method) {
  const [val, release] = await ANNOTATION_SEM.acquire();
  try {
    const startTime = Date.now();
    const isGet = method === 'GET';
    const tableauFile = isGet ? req.query.tableauFile : req.body.tableauFile;
    const matrixFile = isGet ? req.query.matrixFile : req.body.matrixFile;
    const annotationType = (isGet ? req.query.annotationType : req.body.annotationType) || '1H';
    const rawParams = isGet ? req.query.parameters : req.body.parameters;
    const projectId = isGet ? req.query.projectId : req.body.projectId;
    let parameters = {};
    if (rawParams) {
      try {
        parameters = typeof rawParams === 'string' ? JSON.parse(rawParams) : rawParams;
      } catch {}
    }

    if (!tableauFile || !matrixFile) {
      return res.status(400).json({ error: 'tableauFile and matrixFile are required' });
    }

    let tableauPath, matrixPath;
    try {
      tableauPath = safeResolve(uploadsDir, path.basename(tableauFile));
      matrixPath = safeResolve(uploadsDir, path.basename(matrixFile));
    } catch (e) {
      return res.status(400).json({ error: `Invalid file path: ${e.message}` });
    }

    if (!fs.existsSync(tableauPath) || !fs.existsSync(matrixPath)) {
      return res.status(400).json({ error: 'Files not found on server' });
    }

    const taskId = uuidv4();
    const workDir = path.join(getActiveOutputsDir(), taskId);
    ensureDir(workDir);

    const optimizedParams = { ...DEFAULT_LIGHT_PARAMS, ...parameters };

    try {
      const totalSizeMB = (fs.statSync(tableauPath).size + fs.statSync(matrixPath).size) / 1024 / 1024;
      if (totalSizeMB > 1000) {
        optimizedParams.nb_cores = 1;
        optimizedParams.batch_size = 50;
      }
      fs.copyFileSync(tableauPath, path.join(workDir, 'tableau.csv'));
      fs.copyFileSync(matrixPath, path.join(workDir, 'matrix.csv'));
      fs.writeFileSync(path.join(workDir, 'params.json'), JSON.stringify(optimizedParams, null, 2));
    } catch (e) {
      return res.status(400).json({ error: `File preparation error: ${e.message}` });
    }

    const cpuCount = Math.max(1, os.cpus()?.length || 1);
    const parsed = parseInt(optimizedParams.nb_cores, 10);
    const nCores = Number.isFinite(parsed) && parsed >= 1 ? Math.min(parsed, cpuCount) : 1;
    optimizedParams.nb_cores = nCores;
    try {
      fs.writeFileSync(path.join(workDir, 'params.json'), JSON.stringify(optimizedParams, null, 2));
    } catch {}

    const hb = sseInit(res, req.headers.origin);
    sseSend(res, { type: 'taskId', taskId });

    if (projectId) {
      projectService.addTask(projectId, {
        id: taskId,
        type: 'annotation',
        annotationType,
        workDir: path.basename(workDir)
      });
      sseSend(res, { type: 'log', message: `Project: ${projectId}` });
    }

    sseSend(res, { type: 'log', message: `Annotation ${annotationType} started (sem=${val})` });
    sseSend(res, { type: 'log', message: `WorkDir: ${workDir}` });
    sseSend(res, { type: 'log', message: `Workers: ${nCores}, Batch: ${optimizedParams.batch_size}` });
    sseSend(res, { type: 'progress', percent: 1, message: 'Initializing...' });
    sseSend(res, { type: 'progress', percent: 2, message: 'Loading data...' });
    sseSend(res, { type: 'progress', percent: 3, message: 'Preparing R environment...' });

    const rscript = getRscript();
    logWithMemory(`[${taskId}] Rscript selected: ${rscript}`, 'INFO');
    const args = [
      path.join(scriptsDir, scriptMap[annotationType]),
      workDir,
      'tableau.csv',
      'matrix.csv',
      'params.json',
    ];
    const scriptPath = args[0];
    if (!scriptPath || !fs.existsSync(scriptPath)) {
      clearInterval(hb);
      return res.status(400).json({ error: `R script not found for ${annotationType}` });
    }

    const childEnv = buildREnv(MEMORY_LIMIT_MB, rscript);
    logWithMemory(
      `[${taskId}] R env: R_HOME=${childEnv.R_HOME || ''} R_LIBS_USER=${childEnv.R_LIBS_USER || ''} R_LIBS_SITE=${childEnv.R_LIBS_SITE || ''} R_PROFILE_USER=${childEnv.R_PROFILE_USER || ''}`,
      'INFO'
    );
    if (nCores >= 1) {
      childEnv.OMP_NUM_THREADS = String(nCores);
      childEnv.R_NTHREADS = String(nCores);
    }

    const child = spawn(rscript, args, {
      cwd: workDir,
      stdio: ['ignore', 'pipe', 'pipe'],
      shell: false,
      detached: false,
      env: childEnv,
    });

    let outputBuffer = '';
    let lastPercent = 0;
    // No artificial heartbeat - rely on the R script real progress

    if (process.platform === 'win32' && child.pid) {
      const pid = parseInt(child.pid, 10);
      if (Number.isInteger(pid) && pid > 0) {
        try {
          require('child_process').exec(`wmic process where ProcessId=${pid} CALL setpriority "normal"`, () => {});
        } catch {}
      }
    }

    const taskInfo = {
      process: child,
      startTime: new Date(),
      workDir,
      annotationType,
      tableauFile,
      matrixFile,
      parameters: optimizedParams,
      lastLogMinutes: 0,
      safetyTimeout: null,
      canceled: false,
      sseResponse: res,
      heartbeat: hb,
    };
    addTask(taskId, taskInfo);

    // ✅ SLEEP PREVENTION: Empêcher la mise en veille pendant l'annotation
    startPreventSleep(taskId);

    // ✅ PERSISTENCE: Sauvegarder checkpoint initial
    if (ENABLE_TASK_PERSISTENCE) {
      saveCheckpoint(taskId, {
        type: 'annotation',
        status: TaskStatus.RUNNING,
        annotationType,
        progress: { current: 0, total: 100, percent: 0 },
        options: optimizedParams,
        workDir,
        projectId,
        resumeData: {
          tableauFile,
          matrixFile,
          annotationType,
          parameters: optimizedParams,
        },
      });
    }

    child.stdout.on('data', buf => {
      outputBuffer += buf.toString();
      const lines = outputBuffer.split('\n');
      outputBuffer = lines.pop() || '';
      lines.forEach(line => {
        if (!line.trim()) return;
        try {
          const parsed = JSON.parse(line);
          // Update lastPercent when the script emits progress
          if (parsed.type === 'progress' && typeof parsed.percent === 'number') {
            lastPercent = parsed.percent;
            // ✅ PERSISTENCE: Mise à jour du checkpoint à chaque progression
            if (ENABLE_TASK_PERSISTENCE && parsed.percent % 10 === 0) {
              saveCheckpoint(taskId, {
                type: 'annotation',
                status: TaskStatus.RUNNING,
                annotationType,
                progress: { current: parsed.percent, total: 100, percent: parsed.percent },
                options: optimizedParams,
                workDir,
                projectId,
                resumeData: {
                  tableauFile,
                  matrixFile,
                  annotationType,
                  parameters: optimizedParams,
                },
              });
            }
          }
          sseSend(res, parsed);
        } catch {
          sseSend(res, { type: 'log', message: line });
        }
      });
    });

    let stderrBuffer = '';
    let inWarningBlock = false;

    const shouldSilenceRWarningLine = line => {
      const s = (line || '').trim();
      if (!s) return true;
      if (/^warning messages?:/i.test(s)) return true;
      if (/^\d+\s*:/i.test(s) && inWarningBlock) return true;
      if (/was built under r version/i.test(s)) return true;
      if (/\bwarning\b|\bnote\b/i.test(s)) return true;
      return false;
    };

    child.stderr.on('data', buf => {
      stderrBuffer += buf.toString();
      const lines = stderrBuffer.split('\n');
      stderrBuffer = lines.pop() || '';

      for (const rawLine of lines) {
        const line = rawLine.trim();
        if (!line) {
          inWarningBlock = false;
          continue;
        }

        if (/^warning messages?:/i.test(line)) {
          inWarningBlock = true;
        }

        if (inWarningBlock && !/^\d+\s*:/i.test(line) && !/^warning messages?:/i.test(line)) {
          inWarningBlock = false;
        }

        if (shouldSilenceRWarningLine(line)) continue;

        // Anything else on stderr is treated as a real error
        sseSend(res, { type: 'error', message: line });
        logWithMemory(`[${taskId}] STDERR: ${line}`, 'WARN');
      }
    });

    child.on('close', (code, signal) => {
      const task = deleteTask(taskId);
      if (task?.safetyTimeout) clearTimeout(task.safetyTimeout);
      clearInterval(hb);
      
      // ✅ SLEEP PREVENTION: Réactiver la veille
      stopPreventSleep(taskId);
      
      const durationMin = Math.round((Date.now() - startTime) / 60000 * 10) / 10;

      if (task?.canceled) {
        sseSend(res, {
          type: 'completed',
          success: false,
          error: 'Cancelled by user',
          duration: durationMin,
          canceled: true,
        });
        res.end();
        return;
      }

      if (signal) {
        sseSend(res, { type: 'completed', success: false, error: `Signal ${signal}`, duration: durationMin });
        res.end();
        return;
      }
      if (code !== 0) {
        // ✅ PERSISTENCE: Marquer comme échouée
        if (ENABLE_TASK_PERSISTENCE) {
          markTaskFailed(taskId, `R script exited with code ${code}`);
        }
        sseSend(res, {
          type: 'completed',
          success: false,
          error: `R script exited with code ${code}`,
          duration: durationMin,
        });
        res.end();
        return;
      }

      const standardOutput = path.join(workDir, 'annotations.csv');
      const legacyOutput = path.join(workDir, legacyOut[annotationType] || '');
      let outputFile = null;
      if (fs.existsSync(standardOutput)) outputFile = path.basename(standardOutput);
      else if (legacyOut[annotationType] && fs.existsSync(legacyOutput))
        outputFile = path.basename(legacyOutput);

      if (outputFile) {
        const outputPath = path.join(workDir, outputFile);
        const sz = fs.existsSync(outputPath) ? fs.statSync(outputPath).size : 0;
        
        // ✅ PERSISTENCE: Marquer comme complétée
        if (ENABLE_TASK_PERSISTENCE) {
          markTaskCompleted(taskId, { outputFile, downloadUrl: `/outputs/${taskId}/${outputFile}` });
        }
        
        sseSend(res, {
          type: 'completed',
          success: true,
          downloadUrl: `/outputs/${taskId}/${outputFile}`,
          taskId,
          annotationType,
          duration: durationMin,
          outputFile,
          outputSizeMB: Math.round((sz / 1024 / 1024) * 100) / 100,
        });
      } else {
        const csvFiles = fs
          .readdirSync(workDir)
          .filter(f => f.endsWith('.csv') && !['tableau.csv', 'matrix.csv'].includes(f));
        
        // ✅ PERSISTENCE: Marquer comme échouée (pas de fichier output)
        if (ENABLE_TASK_PERSISTENCE) {
          markTaskFailed(taskId, 'No output file found');
        }
        
        sseSend(res, {
          type: 'completed',
          success: false,
          error: 'No output file found',
          diagnostics: { csvFiles, workDir },
          duration: durationMin,
        });
      }
      res.end();
    });

    child.on('error', err => {
      // ✅ PERSISTENCE: Marquer comme échouée
      if (ENABLE_TASK_PERSISTENCE) {
        markTaskFailed(taskId, `R spawn error: ${err.message}`);
      }
      deleteTask(taskId);
      if (taskInfo.safetyTimeout) clearTimeout(taskInfo.safetyTimeout);
      clearInterval(hb);
      sseSend(res, { type: 'completed', success: false, error: `R spawn error: ${err.message}` });
      res.end();
    });

    taskInfo.safetyTimeout = setTimeout(() => {
      if (deleteTask(taskId)) {
        killProcessGracefully(child, taskId);
      }
    }, PROCESS_TIMEOUT_HOURS * 60 * 60 * 1000);

    req.on('close', () => {
      clearInterval(hb);
      if (KILL_ON_CLIENT_CLOSE) {
        killProcessGracefully(child, taskId);
        deleteTask(taskId);
      } else if (ENABLE_TASK_PERSISTENCE) {
        // ✅ PERSISTENCE: La tâche continue en arrière-plan
        logWithMemory(`[${taskId}] Client disconnected - annotation continues in background`, 'INFO');
      }
    });
  } catch (e) {
    logWithMemory(`handleAnnotation error: ${e.message}`, 'ERROR');
    res.status(500).json({ error: e.message });
  } finally {
    release();
  }
}

async function handleAnnotationFromPrediction(req, res, method) {
  const [val, release] = await ANNOTATION_SEM.acquire();
  try {
    const startTime = Date.now();
    const isGet = method === 'GET';
    const predictionTaskId = isGet ? req.query.predictionTaskId : req.body.predictionTaskId;
    const matrixFile = isGet ? req.query.matrixFile : req.body.matrixFile;
    const annotationType = (isGet ? req.query.annotationType : req.body.annotationType) || '1H';
    const isOptimization = isGet ? req.query.isOptimization === 'true' : !!req.body.isOptimization;
    const rawParams = isGet ? req.query.parameters : req.body.parameters;
    let parameters = {};
    if (rawParams) {
      try {
        parameters = typeof rawParams === 'string' ? JSON.parse(rawParams) : rawParams;
      } catch {}
    }

    if (!predictionTaskId || !matrixFile) {
      return res.status(400).json({ error: 'predictionTaskId and matrixFile are required' });
    }

    const predictionDir = path.join(getActiveOutputsDir(), predictionTaskId);
    if (!fs.existsSync(predictionDir)) {
      return res.status(400).json({ error: 'Prediction task not found' });
    }

    const findCsv = (files, token) => {
      const exact = files.find(f => f.startsWith(`prediction_${token}_`) && f.toLowerCase().endsWith('.csv'));
      if (exact) return exact;
      const tok = String(token).toLowerCase();
      return files.find(f => f.toLowerCase().endsWith('.csv') && f.toLowerCase().includes(tok)) || null;
    };

    const files = fs.readdirSync(predictionDir);
    const predictionCsvMap = {
      '1H': findCsv(files, '1H'),
      '13C': findCsv(files, '13C'),
      COSY: findCsv(files, 'COSY'),
      HSQC: findCsv(files, 'HSQC'),
    };

    const predictionCsv = predictionCsvMap[annotationType];
    if (!predictionCsv) {
      return res.status(400).json({
        error: `Prediction CSV ${annotationType} not found in task ${predictionTaskId}`,
      });
    }

    const predictionCsvPath = path.join(predictionDir, predictionCsv);
    const matrixPath = path.join(uploadsDir, path.basename(matrixFile));
    if (!fs.existsSync(matrixPath)) {
      return res.status(400).json({ error: 'Matrix file not found on server' });
    }

    const taskId = uuidv4();
    const workDir = predictionDir;

    // Use separate folder for optimization to avoid overwriting annotation results
    const outputSubfolder = isOptimization ? 'optimization' : 'annotations';
    const annotationsDir = path.join(workDir, outputSubfolder);
    ensureDir(annotationsDir);

    const timestamp = new Date()
      .toISOString()
      .replace(/[:.]/g, '-')
      .slice(0, 19);
    const tempAnnotDir = path.join(annotationsDir, `temp_${annotationType}_${timestamp}`);
    ensureDir(tempAnnotDir);

    const optimizedParams = { ...DEFAULT_LIGHT_PARAMS, ...parameters };

    try {
      fs.copyFileSync(predictionCsvPath, path.join(tempAnnotDir, `prediction_${annotationType}.csv`));
      fs.copyFileSync(matrixPath, path.join(tempAnnotDir, 'matrix.csv'));
      fs.writeFileSync(path.join(tempAnnotDir, 'params.json'), JSON.stringify(optimizedParams, null, 2));
    } catch (e) {
      return res.status(400).json({ error: `File preparation error: ${e.message}` });
    }

    const cpuCount = Math.max(1, os.cpus()?.length || 1);
    const parsed = parseInt(optimizedParams.nb_cores, 10);
    const nCores = Number.isFinite(parsed) && parsed >= 1 ? Math.min(parsed, cpuCount) : 1;
    optimizedParams.nb_cores = nCores;

    const hb = sseInit(res, req.headers.origin);
    sseSend(res, { type: 'taskId', taskId });
    sseSend(res, { type: 'log', message: `Annotation ${annotationType} from prediction ${predictionTaskId}` });
    sseSend(res, { type: 'log', message: `Prediction CSV: ${predictionCsv}` });
    sseSend(res, { type: 'log', message: `Main folder: ${workDir}` });
    sseSend(res, { type: 'log', message: `Temporary folder: ${path.basename(tempAnnotDir)}` });
    sseSend(res, { type: 'progress', percent: 1, message: 'Initializing...' });

    const rscript = getRscript();
    const args = [
      path.join(scriptsDir, scriptMap[annotationType]),
      tempAnnotDir,
      `prediction_${annotationType}.csv`,
      'matrix.csv',
      'params.json',
    ];
    const scriptPath = args[0];
    if (!scriptPath || !fs.existsSync(scriptPath)) {
      clearInterval(hb);
      try {
        fs.rmSync(tempAnnotDir, { recursive: true, force: true });
      } catch {}
      return res.status(400).json({ error: `R script not found for ${annotationType}` });
    }

      const childEnv = buildREnv(MEMORY_LIMIT_MB, rscript);
    if (nCores >= 1) {
      childEnv.OMP_NUM_THREADS = String(nCores);
      childEnv.R_NTHREADS = String(nCores);
    }

    const child = spawn(rscript, args, {
      cwd: tempAnnotDir,
      stdio: ['ignore', 'pipe', 'pipe'],
      shell: false,
      detached: false,
      env: childEnv,
    });

    let outputBuffer = '';
    let lastPercent = 0;
    // No artificial heartbeat - rely on the R script real progress

    const taskInfo = {
      process: child,
      startTime: new Date(),
      workDir: workDir,
      tempDir: tempAnnotDir,
      annotationType,
      predictionTaskId,
      predictionCsv,
      matrixFile,
      parameters: optimizedParams,
      lastLogMinutes: 0,
      safetyTimeout: null,
      canceled: false,
      sseResponse: res,
      heartbeat: hb,
    };
    addTask(taskId, taskInfo);

    // ✅ SLEEP PREVENTION: Empêcher la mise en veille
    startPreventSleep(taskId);

    child.stdout.on('data', buf => {
      outputBuffer += buf.toString();
      const lines = outputBuffer.split('\n');
      outputBuffer = lines.pop() || '';

      lines.forEach(line => {
        if (!line.trim()) return;

        try {
          const parsed = JSON.parse(line);

          if (parsed.type === 'progress') {
            let p = Number(parsed.percent);

            if (p >= 0 && p <= 1) {
              p = p * 100;
            }

            p = Math.max(lastPercent, Math.min(100, Math.floor(p)));

            if (p > lastPercent) {
              lastPercent = p;
              sseSend(res, {
                type: 'progress',
                percent: p,
                message: parsed.message || `${p}%`,
                current: parsed.current,
                total: parsed.total,
              });
            }
          } else {
            sseSend(res, parsed);
          }
        } catch {
          const msg = line.toLowerCase();
          let estimatedPercent = lastPercent;

          if (msg.includes('loading')) {
            estimatedPercent = Math.max(lastPercent, 10);
          } else if (msg.includes('computing') || msg.includes('processing')) {
            estimatedPercent = Math.max(lastPercent, 30);
          } else if (msg.includes('annotation') || msg.includes('matching')) {
            estimatedPercent = Math.max(lastPercent, 50);
          } else if (msg.includes('writing') || msg.includes('saving')) {
            estimatedPercent = Math.max(lastPercent, 80);
          } else if (msg.includes('done') || msg.includes('finished')) {
            estimatedPercent = Math.max(lastPercent, 95);
          }

          if (estimatedPercent > lastPercent) {
            lastPercent = estimatedPercent;
            sseSend(res, {
              type: 'progress',
              percent: estimatedPercent,
              message: line.substring(0, 100),
            });
          }

          sseSend(res, { type: 'log', message: line });
        }
      });
    });

    let stderrBuffer = '';
    let inWarningBlock = false;

    const shouldSilenceRWarningLine = line => {
      const s = (line || '').trim();
      if (!s) return true;
      if (/^warning messages?:/i.test(s)) return true;
      if (/^\d+\s*:/i.test(s) && inWarningBlock) return true;
      if (/was built under r version/i.test(s)) return true;
      if (/\bwarning\b|\bnote\b/i.test(s)) return true;
      return false;
    };

    child.stderr.on('data', buf => {
      stderrBuffer += buf.toString();
      const lines = stderrBuffer.split('\n');
      stderrBuffer = lines.pop() || '';

      for (const rawLine of lines) {
        const line = rawLine.trim();
        if (!line) {
          inWarningBlock = false;
          continue;
        }

        if (/^warning messages?:/i.test(line)) {
          inWarningBlock = true;
        }

        if (inWarningBlock && !/^\d+\s*:/i.test(line) && !/^warning messages?:/i.test(line)) {
          inWarningBlock = false;
        }

        if (shouldSilenceRWarningLine(line)) continue;

        sseSend(res, { type: 'error', message: line });
        logWithMemory(`[${taskId}] STDERR: ${line}`, 'WARN');
      }
    });

    child.on('close', async (code, signal) => {
      const task = deleteTask(taskId);
      if (task?.safetyTimeout) clearTimeout(task.safetyTimeout);
      clearInterval(hb);
      
      // ✅ SLEEP PREVENTION: Réactiver la veille
      stopPreventSleep(taskId);
      
      const durationMin = Math.round((Date.now() - startTime) / 60000 * 10) / 10;

      if (task?.canceled) {
        try {
          fs.rmSync(tempAnnotDir, { recursive: true, force: true });
        } catch {}
        sseSend(res, {
          type: 'completed',
          success: false,
          error: 'Cancelled by user',
          duration: durationMin,
          canceled: true,
        });
        res.end();
        return;
      }

      if (signal) {
        try {
          fs.rmSync(tempAnnotDir, { recursive: true, force: true });
        } catch {}
        sseSend(res, { type: 'completed', success: false, error: `Signal ${signal}`, duration: durationMin });
        res.end();
        return;
      }

      if (code !== 0) {
        try {
          fs.rmSync(tempAnnotDir, { recursive: true, force: true });
        } catch {}
        sseSend(res, {
          type: 'completed',
          success: false,
          error: `R script exited with code ${code}`,
          duration: durationMin,
        });
        res.end();
        return;
      }

      const standardOutput = path.join(tempAnnotDir, 'annotations.csv');
      const legacyOutput = path.join(tempAnnotDir, legacyOut[annotationType] || '');
      let outputFile = null;
      let sourcePath = null;

      if (fs.existsSync(standardOutput)) {
        outputFile = 'annotations.csv';
        sourcePath = standardOutput;
      } else if (legacyOut[annotationType] && fs.existsSync(legacyOutput)) {
        outputFile = legacyOut[annotationType];
        sourcePath = legacyOutput;
      }

      if (outputFile && sourcePath) {
        try {
          // Generate filename with tolerance info for optimization runs
          let latestFilename;
          if (isOptimization) {
            // Extract tolerance from parameters for filename
            const tolStr = getToleranceStringForFilename(annotationType, optimizedParams);
            latestFilename = `annotation_${annotationType}_${tolStr}.csv`;
          } else {
            latestFilename = `annotation_${annotationType}_latest.csv`;
          }
          const latestPath = path.join(annotationsDir, latestFilename);

          // Only create history for normal annotations, not for optimization runs
          if (!isOptimization && fs.existsSync(latestPath)) {
            const historyDir = path.join(annotationsDir, 'history');
            ensureDir(historyDir);

            const stats = fs.statSync(latestPath);
            const archiveTimestamp = stats.mtime
              .toISOString()
              .replace(/[:.]/g, '-')
              .slice(0, 19);
            const archivedFilename = `annotation_${annotationType}_${archiveTimestamp}.csv`;
            const archivedPath = path.join(historyDir, archivedFilename);

            fs.renameSync(latestPath, archivedPath);
            sseSend(res, { type: 'log', message: `Archive created: ${archivedFilename}` });
          }

          fs.copyFileSync(sourcePath, latestPath);

          let fileReady = false;
          let attempts = 0;
          const maxAttempts = 10;

          while (!fileReady && attempts < maxAttempts) {
            try {
              const stats = fs.statSync(latestPath);
              if (stats.size > 0) {
                const fd = fs.openSync(latestPath, 'r');
                const buffer = Buffer.alloc(100);
                fs.readSync(fd, buffer, 0, 100, 0);
                fs.closeSync(fd);
                fileReady = true;
                sseSend(res, { type: 'log', message: `File verified and accessible (${stats.size} bytes)` });
              }
            } catch (e) {
              attempts++;
              if (attempts < maxAttempts) {
                await new Promise(resolve => setTimeout(resolve, 100));
              }
            }
          }

          if (!fileReady) {
            throw new Error('Timeout: file not accessible after write');
          }

          sseSend(res, { type: 'log', message: `Result saved: ${latestFilename}` });

          try {
            fs.rmSync(tempAnnotDir, { recursive: true, force: true });
            sseSend(res, { type: 'log', message: 'Temporary folder cleaned up' });
          } catch (cleanupErr) {
            sseSend(res, { type: 'log', message: `Cleanup warning: ${cleanupErr.message}` });
          }

          const sz = fs.existsSync(latestPath) ? fs.statSync(latestPath).size : 0;
          const fileTimestamp = Date.now();

          sseSend(res, {
            type: 'completed',
            success: true,
            downloadUrl: `/outputs/${predictionTaskId}/${outputSubfolder}/${latestFilename}?v=${fileTimestamp}`,
            taskId,
            annotationType,
            predictionTaskId,
            predictionCsv,
            duration: durationMin,
            outputFile: latestFilename,
            outputSizeMB: Math.round((sz / 1024 / 1024) * 100) / 100,
            workDir: workDir,
            annotationsDir: annotationsDir,
            isOptimization: isOptimization,
            timestamp: fileTimestamp,
          });
        } catch (copyErr) {
          sseSend(res, {
            type: 'completed',
            success: false,
            error: `Result save error: ${copyErr.message}`,
            duration: durationMin,
          });
        }
      } else {
        const csvFiles = fs
          .readdirSync(tempAnnotDir)
          .filter(f => f.endsWith('.csv') && !['matrix.csv', `prediction_${annotationType}.csv`].includes(f));

        try {
          fs.rmSync(tempAnnotDir, { recursive: true, force: true });
        } catch {}

        sseSend(res, {
          type: 'completed',
          success: false,
          error: 'No output file found',
          diagnostics: { csvFiles, tempAnnotDir: path.basename(tempAnnotDir) },
          duration: durationMin,
        });
      }
      res.end();
    });

    child.on('error', err => {
      deleteTask(taskId);
      if (taskInfo.safetyTimeout) clearTimeout(taskInfo.safetyTimeout);
      clearInterval(hb);
      try {
        fs.rmSync(tempAnnotDir, { recursive: true, force: true });
      } catch {}
      sseSend(res, { type: 'completed', success: false, error: `R spawn error: ${err.message}` });
      res.end();
    });

    taskInfo.safetyTimeout = setTimeout(() => {
      if (deleteTask(taskId)) {
        killProcessGracefully(child, taskId);
        try {
          fs.rmSync(tempAnnotDir, { recursive: true, force: true });
        } catch {}
      }
    }, PROCESS_TIMEOUT_HOURS * 60 * 60 * 1000);

    req.on('close', () => {
      clearInterval(hb);
      if (KILL_ON_CLIENT_CLOSE) {
        killProcessGracefully(child, taskId);
        deleteTask(taskId);
      } else if (ENABLE_TASK_PERSISTENCE) {
        // ✅ PERSISTENCE: La tâche continue en arrière-plan
        logWithMemory(`[${taskId}] Client disconnected - annotation continues in background`, 'INFO');
      }
    });
  } catch (e) {
    logWithMemory(`handleAnnotationFromPrediction error: ${e.message}`, 'ERROR');
    res.status(500).json({ error: e.message });
  } finally {
    release();
  }
}

module.exports = {
  handleAnnotation,
  handleAnnotationFromPrediction,
};