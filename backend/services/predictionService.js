const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const { v4: uuidv4 } = require('uuid');

const { PROCESS_TIMEOUT_HOURS, KILL_ON_CLIENT_CLOSE } = require('../config');
const { PRED_SCRIPT, outputsDir, uploadsDir } = require('../config/paths');
const { ensureDir } = require('../utils/filesystem');
const { getPythonExe, buildPythonEnv, detectDatasetDir } = require('../utils/spawn');
const { logWithMemory } = require('../utils/logger');
const { sseInit, sseSend } = require('../middleware/sse');
const { addTask, deleteTask, getTask, killProcessGracefully } = require('./processManager');
const { runExtractorsOnSDF } = require('./extractorService');
const projectService = require('./projectService');

function normalizePredictOptions(raw = {}) {
  const o = { ...(raw || {}) };

  o.extractors = { ...(o.extractors || {}) };
  o.extractors.c13 = { ...(o.extractors.c13 || {}) };
  o.extractors.hsqc = { ...(o.extractors.hsqc || {}) };

  const c13NestedCarbonMode = o?.c13?.carbonMode ?? o?.c13?.carbon_mode ?? o?.c13?.carbonmode;

  if (c13NestedCarbonMode != null && !o.extractors.c13.carbon_mode) {
    o.extractors.c13.carbon_mode = String(c13NestedCarbonMode)
      .toLowerCase()
      .startsWith('all')
      ? 'all'
      : 'nonquaternary';
  }

  const hasHsqcNested =
    o?.hsqc && ('fold' in o.hsqc || 'foldLow' in o.hsqc || 'foldHigh' in o.hsqc);

  if (hasHsqcNested) {
    if ('fold' in o.hsqc && o.extractors.hsqc.fold == null) {
      o.extractors.hsqc.fold = !!o.hsqc.fold;
    }
    if ('foldLow' in o.hsqc && o.extractors.hsqc.fold_low == null) {
      o.extractors.hsqc.fold_low = Number(o.hsqc.foldLow);
    }
    if ('foldHigh' in o.hsqc && o.extractors.hsqc.fold_high == null) {
      o.extractors.hsqc.fold_high = Number(o.hsqc.foldHigh);
    }
  }

  if (o.hsqc_fold !== undefined) {
    o.extractors.hsqc.fold = !!o.hsqc_fold;
  }
  if (o.hsqc_fold_low !== undefined) {
    o.extractors.hsqc.fold_low = Number(o.hsqc_fold_low);
  }
  if (o.hsqc_fold_high !== undefined) {
    o.extractors.hsqc.fold_high = Number(o.hsqc_fold_high);
  }

  ['fold_low', 'fold_high'].forEach(k => {
    if (o.extractors.hsqc[k] != null && !Number.isFinite(+o.extractors.hsqc[k])) {
      delete o.extractors.hsqc[k];
    }
  });

  return o;
}

function extractSmilesFromSdf(sdfPath, pythonExe, env) {
  return new Promise((resolve, reject) => {
    // Use the new script instead of inline code
    // If we are in a bundled environment (Electron), pythonExe will be the bundled exe
    // and we should pass the command name as the first argument.
    // If we are in dev (Node), pythonExe is 'python' and we need to pass the script path.
    
    const isBundledEnv = process.env.IS_ELECTRON_BUNDLED === 'true';
    // Check if we're really using a bundled exe (not just the env var)
    const isActuallyBundled = isBundledEnv && pythonExe.endsWith('.exe') && fs.existsSync(pythonExe);
    
    let args;
    if (isActuallyBundled) {
        // Bundled: exe extract_smiles <sdfPath>
        args = ['extract_smiles', sdfPath];
    } else {
        // Dev: python scripts/extract_smiles.py <sdfPath>
        const scriptPath = path.join(__dirname, '..', 'scripts', 'extract_smiles.py');
        args = [scriptPath, sdfPath];
    }

    const p = spawn(pythonExe, args, { env });
    p.on('error', (err) => {
      reject(new Error(`Failed to spawn Python for SDF conversion: ${err.message}`));
    });
    let out = '',
      err = '';
    p.stdout.on('data', d => (out += d.toString()));
    p.stderr.on('data', d => (err += d.toString()));
    p.on('close', code => {
      if (code !== 0) return reject(new Error(err || `convert SDF->SMILES code ${code}`));
      try {
        const j = JSON.parse(out.trim());
        if (!j.ok) return reject(new Error(j.error || 'conversion failed'));
        resolve(j.smiles);
      } catch (e) {
        reject(new Error('parse conversion output: ' + e.message));
      }
    });
  });
}

function handlePrediction(req, res) {
  const startTime = Date.now();
  try {
    const hb = sseInit(res, req.headers.origin);
    const taskId = uuidv4();
    const workDir = path.join(outputsDir, taskId);
    ensureDir(workDir);

    const SDF_WEIGHT = 95; // 0-95% for SDF, 95-100% for CSV extractors
    const EXTRACTOR_WEIGHT = 5; // 5% distributed across extractors
    let totalMolecules = 1;
    let doneMolecules = 0;
    let lastPercent = 0;

    const bodyMode = ((req.body?.mode || req.body?.source || '') + '').toUpperCase();
    let mode = bodyMode;
    const smiles = req.body?.smiles || '';
    const sdfFile = req.body?.sdfFile ? path.basename(req.body.sdfFile) : null;
    const csvFile = req.body?.csvFile ? path.basename(req.body.csvFile) : null;

    const selectedMethods = req.body?.methods || req.body?.kinds || ['1H', '13C', 'COSY', 'HSQC'];
    const projectId = req.body?.projectId;

    if (!mode) {
      if (csvFile) mode = 'CSV';
      else if (sdfFile) mode = 'SDF';
      else mode = 'SMILES';
    }

    const optionsIn = req.body?.options || {};
    const options = normalizePredictOptions({ ...optionsIn, run_extractors: true });
    const model1H = options.model_1h_path || process.env.MODEL_1H_PATH;
    const model13C = options.model_13c_path || process.env.MODEL_13C_PATH;
    if (!model1H || !model13C) {
      clearInterval(hb);
      sseSend(res, { type: 'error', message: 'MODEL_1H_PATH and MODEL_13C_PATH must be set.' });
      sseSend(res, { type: 'completed', success: false, error: 'Missing model paths' });
      res.end();
      return;
    }
    const propertyKey = options.property_key || process.env.PROPERTY_KEY || '(1)';

    sseSend(res, { type: 'taskId', taskId });
  sseSend(res, { type: 'log', message: `Prediction started: ${mode}` });
    sseSend(res, { type: 'log', message: `WorkDir: ${workDir}` });
  sseSend(res, { type: 'log', message: `Methods: ${selectedMethods.join(', ')}` });

    if (projectId) {
      projectService.addTask(projectId, {
        id: taskId,
        type: 'prediction',
        mode,
        methods: selectedMethods,
        workDir: path.basename(workDir)
      });
      sseSend(res, { type: 'log', message: `Project: ${projectId}` });
    }

    sseSend(res, { type: 'progress', percent: 0, phase: 'prediction', message: 'Initializing...' });
    lastPercent = 0;

    // Determine whether we're using the PyInstaller bundle or the Python script
    const PYTHON = getPythonExe(true); // forML = true
    
    console.log('[Prediction] PYTHON (for ML):', PYTHON);

    const isBundledExe = PYTHON.toLowerCase().endsWith('nmr_engine_production.exe') || 
                         PYTHON.toLowerCase().endsWith('nmr_engine_production') ||
                         PYTHON.toLowerCase().endsWith('nmr_engine.exe') || 
                         PYTHON.toLowerCase().endsWith('nmr_engine');

    if (!isBundledExe && !fs.existsSync(PRED_SCRIPT)) {
      clearInterval(hb);
      sseSend(res, { type: 'error', message: `Python script not found: ${PRED_SCRIPT}` });
      sseSend(res, { type: 'completed', success: false, error: 'integrated_nmr_prediction.py not found' });
      res.end();
      return;
    }

    const childEnv = buildPythonEnv();

    const datasetDir = childEnv.NMR_DATA_DIR || detectDatasetDir();
    if (datasetDir) {
      sseSend(res, { type: 'log', message: `Dataset dir (NMR_DATA_DIR): ${datasetDir}` });
    } else {
      sseSend(res, { type: 'log', message: '[WARN] NMR_DATA_DIR not found - fallback mean/std.' });
    }

    // Construct args based on whether we use bundled exe or python script
    const args = [];
    if (isBundledExe) {
        args.push('integrated_prediction');
    } else {
        args.push(PRED_SCRIPT);
    }

    args.push(
      '--output_dir',
      workDir,
      '--task_id',
      taskId,
      '--model_1h_path',
      model1H,
      '--model_13c_path',
      model13C,
      '--property_key',
      propertyKey,
    );

    const addSmilesArg = smi => args.push('--smiles', smi);

      const useCSV = () => {
      if (!csvFile) {
        clearInterval(hb);
        sseSend(res, { type: 'error', message: 'csvFile is required in CSV mode' });
        sseSend(res, { type: 'completed', success: false, error: 'csvFile required' });
        res.end();
        return false;
      }
      const src = path.join(uploadsDir, csvFile);
      if (!fs.existsSync(src)) {
        clearInterval(hb);
        sseSend(res, { type: 'error', message: 'CSV not found on server' });
        sseSend(res, { type: 'completed', success: false, error: 'CSV not found' });
        res.end();
        return false;
      }      const stripBOM = s => s.replace(/^\uFEFF/, '');
      const removeDiacritics = s =>
        s
          .normalize('NFD')
          .replace(/[\u0300-\u036f]/g, '');
      const normalizeHeaderName = s => {
        s = stripBOM(String(s || ''))
          .replace(/^"(.*)"$/, '$1')
          .trim();
        s = removeDiacritics(s).toLowerCase();
        s = s
          .replace(/[\s\-./\\|;:,]+/g, '_')
          .replace(/[^a-z0-9_]/g, '')
          .replace(/_+/g, '_')
          .replace(/^_+|_+$/g, '');
        return s;
      };
      const splitCSVLine = (line, delim) => {
        const out = [];
        let cur = '',
          inQ = false;
        for (let i = 0; i < line.length; i++) {
          const ch = line[i];
          if (ch === '"') {
            if (inQ && line[i + 1] === '"') {
              cur += '"';
              i++;
            } else inQ = !inQ;
          } else if (ch === delim && !inQ) {
            out.push(cur);
            cur = '';
          } else {
            cur += ch;
          }
        }
        out.push(cur);
        return out;
      };
      const chooseDelimiter = header => {
        const cands = [',', ';', '\t', '|'];
        let best = ',',
          bestScore = -1;
        for (const d of cands) {
          const n = splitCSVLine(header, d).length;
          if (n > bestScore) {
            best = d;
            bestScore = n;
          }
        }
        return best;
      };
      const needsQuote = (cell, delim) => cell.includes(delim) || /["\r\n]/.test(cell);
      const quoteCell = (cell, delim) =>
        needsQuote(cell, delim) ? `"${cell.replace(/"/g, '""')}"` : cell;

      let raw = '';
      try {
        raw = fs.readFileSync(src, { encoding: 'utf8' });
      } catch (e) {
        clearInterval(hb);
        sseSend(res, { type: 'error', message: `Unable to read CSV: ${e.message}` });
        sseSend(res, { type: 'completed', success: false, error: e.message });
        res.end();
        return false;
      }
      raw = stripBOM(raw);
      const nl = raw.indexOf('\n');
      const headerLineRaw = nl >= 0 ? raw.slice(0, nl).replace(/\r$/, '') : raw;
      const rest = nl >= 0 ? raw.slice(nl + 1) : '';

      const dataLines = rest
        .split(/\r?\n/)
        .map(l => l.trim())
        .filter(l => {
          if (!l) return false;
          if (l.startsWith('#')) return false;
          if (l.startsWith('//')) return false;
          return true;
        });

      totalMolecules = Math.max(1, dataLines.length);
  sseSend(res, { type: 'log', message: `CSV detected: ${totalMolecules} molecules` });

      const delim = chooseDelimiter(headerLineRaw);
      let headerCells = splitCSVLine(headerLineRaw, delim).map(c => c.replace(/\r$/, ''));
      const headerPretty = headerCells.map(c =>
        c
          .replace(/^"(.*)"$/, '$1')
          .trim()
      );
      sseSend(res, {
        type: 'log',
        message: `CSV detected — delim='${delim === '\t' ? '\\t' : delim}', columns=[${headerPretty.join(', ')}]`,
      });

      const idWantedRaw = options.id_column
        ? String(options.id_column)
        : req.body?.id_column
        ? String(req.body.id_column)
        : null;
      const smiWantedRaw = options.smiles_column
        ? String(options.smiles_column)
        : req.body?.smiles_column
        ? String(req.body.smiles_column)
        : null;

      const normalizedMap = new Map();
      headerCells.forEach((name, idx) => normalizedMap.set(normalizeHeaderName(name), idx));

      const findByHint = hint => {
        const norm = normalizeHeaderName(hint);
        if (normalizedMap.has(norm)) return normalizedMap.get(norm);
        return null;
      };

      let idIdx = idWantedRaw ? findByHint(idWantedRaw) : null;
      let smiIdx = smiWantedRaw ? findByHint(smiWantedRaw) : null;

      const candidatesId = [
        'id',
        'sample_id',
        'sampleid',
        'record_id',
        'mol_id',
        'molid',
        'name',
        'sample',
        'index',
      ];
      const candidatesSmiles = [
        'smiles',
        'smile',
        'smiles_string',
        'smilestring',
        'smiles_col',
        'smilescode',
      ];
      if (idIdx == null) {
        for (const c of candidatesId) {
          if (normalizedMap.has(c)) {
            idIdx = normalizedMap.get(c);
            break;
          }
        }
      }
      if (smiIdx == null) {
        for (const c of candidatesSmiles) {
          if (normalizedMap.has(c)) {
            smiIdx = normalizedMap.get(c);
            break;
          }
        }
      }

      if (headerCells.length === 2 && (idIdx == null || smiIdx == null)) {
        const a = normalizeHeaderName(headerCells[0]),
          b = normalizeHeaderName(headerCells[1]);
        if (a.includes('smile') || b.includes('smile')) {
          smiIdx = a.includes('smile') ? 0 : 1;
          idIdx = smiIdx === 0 ? 1 : 0;
        }
      }

      if (idIdx == null || smiIdx == null) {
        clearInterval(hb);
        const msg =
          `Unable to find the id/smiles columns.\nHeaders: [${headerPretty.join(', ')}]\n` +
          `Set "id column" and "smiles column" in the UI or rename the columns.`;
        sseSend(res, { type: 'error', message: msg });
        sseSend(res, { type: 'completed', success: false, error: msg });
        res.end();
        return false;
      }

      const usedId = headerPretty[idIdx];
      const usedSmi = headerPretty[smiIdx];
      sseSend(res, { type: 'log', message: `Using columns — id='${usedId}', smiles='${usedSmi}'` });

      const normalizedHeader = headerCells.map((c, i) => {
        if (i === idIdx) return 'id';
        if (i === smiIdx) return 'smiles';
        return c.replace(/\r$/, '');
      });
      const newHeaderLine = normalizedHeader.map(c => quoteCell(c, delim)).join(delim);

      const outPath = path.join(workDir, 'input.norm.csv');
      fs.writeFileSync(outPath, `${newHeaderLine}\n${rest}`, 'utf8');

      args.push('--csv', outPath);
      args.push('--id_column', 'id');
      args.push('--smiles_column', 'smiles');

      return true;
    };

    const prepareInput = async () => {
      if (mode === 'SMILES') {
        if (!smiles || typeof smiles !== 'string') {
          clearInterval(hb);
          sseSend(res, { type: 'error', message: 'smiles is required in SMILES mode' });
          sseSend(res, { type: 'completed', success: false, error: 'smiles required' });
          res.end();
          return false;
        }
        totalMolecules = 1;
        addSmilesArg(smiles);
        return true;
      }
      if (mode === 'SDF') {
        if (!sdfFile) {
          clearInterval(hb);
          sseSend(res, { type: 'error', message: 'sdfFile is required in SDF mode' });
          sseSend(res, { type: 'completed', success: false, error: 'sdfFile required' });
          res.end();
          return false;
        }
        const src = path.join(uploadsDir, sdfFile);
        if (!fs.existsSync(src)) {
          clearInterval(hb);
          sseSend(res, { type: 'error', message: 'SDF not found on server' });
          sseSend(res, { type: 'completed', success: false, error: 'SDF not found' });
          res.end();
          return false;
        }
        const localSdf = path.join(workDir, 'input.sdf');
        fs.copyFileSync(src, localSdf);
        sseSend(res, { type: 'log', message: 'Converting SDF to SMILES...' });
        try {
          const smi = await extractSmilesFromSdf(localSdf, PYTHON, childEnv);
          totalMolecules = 1;
          addSmilesArg(smi);
          return true;
        } catch (e) {
          clearInterval(hb);
          sseSend(res, { type: 'error', message: `SDF to SMILES conversion failed: ${e.message}` });
          sseSend(res, { type: 'completed', success: false, error: e.message });
          res.end();
          return false;
        }
      }
      if (mode === 'CSV') {
        return useCSV();
      }
      clearInterval(hb);
      sseSend(res, { type: 'error', message: `Unsupported mode: ${mode}` });
      sseSend(res, { type: 'completed', success: false, error: `Unsupported mode: ${mode}` });
      res.end();
      return false;
    };

    (async () => {
      sseSend(res, {
        type: 'progress',
        percent: 0,
        phase: 'prediction',
        message: 'Initializing Python...',
      });
      lastPercent = 0;

      const ok = await prepareInput();
      if (!ok) return;

      if (options.device) args.push('--device', String(options.device));
      if (options.batch_size) args.push('--batch_size', String(options.batch_size));
      if (options.mol_id) args.push('--mol_id', String(options.mol_id));

      sseSend(res, { type: 'log', message: `Python: ${PYTHON}` });
      sseSend(res, { type: 'log', message: `Script: ${PRED_SCRIPT}` });
      sseSend(res, {
        type: 'log',
        message: `Args: ${args.filter(a => a !== PRED_SCRIPT).join(' ')}`,
      });
      sseSend(res, { type: 'log', message: `Total molecules expected: ${totalMolecules}` });

      console.log('[Prediction] Spawning:', PYTHON, args.slice(0, 3).join(' '), '...');
      
      let child;
      try {
        child = spawn(PYTHON, args, {
          cwd: workDir,
          stdio: ['ignore', 'pipe', 'pipe'],
          shell: false,
          detached: false,
          env: childEnv,
        });

        child.on('error', (err) => {
          console.error('[Prediction] Spawn error:', err);
          sseSend(res, { type: 'error', message: `Python launch error: ${err.message}` });
          clearInterval(hb);
        });
      } catch (spawnError) {
        console.error('[Prediction] Spawn failed:', spawnError.message);
        sseSend(res, { type: 'error', message: `Python launch failed: ${spawnError.message}` });
        sseSend(res, { type: 'completed', success: false, error: spawnError.message });
        clearInterval(hb);
        return;
      }
      
      // Handle spawn error (e.g., ENOENT if python not found)
      child.on('error', (err) => {
        console.error('[Prediction] Child process error:', err.message);
        sseSend(res, { type: 'error', message: `Process error: ${err.message}` });
        sseSend(res, { type: 'completed', success: false, error: err.message });
        clearInterval(hb);
        deleteTask(taskId);
      });

      let buffer = '';

      const taskInfo = {
        process: child,
        startTime: new Date(),
        workDir,
        mode,
        lastLogMinutes: 0,
        safetyTimeout: null,
        canceled: false,
        sseResponse: res,
        heartbeat: hb,
      };
      addTask(taskId, taskInfo);

      child.stdout.on('data', buf => {
        buffer += buf.toString();
        const lines = buffer.split('\n');
        buffer = lines.pop() || '';

        lines.forEach(line => {
          const s = line.trim();
          if (!s) return;

          let parsed = null;
          try {
            parsed = JSON.parse(s);
          } catch (_) {}

          if (parsed && typeof parsed === 'object' && parsed.type) {
            if (parsed.type === 'progress') {
              // Weight the Python script percentage into the 0-SDF_WEIGHT range (0-80%)
              // SDF phase is 0%->80%, extractors are 80%->100%
              if (typeof parsed.percent === 'number') {
                // Map 0-100 script percent to 0-SDF_WEIGHT
                const weightedPercent = Math.min(SDF_WEIGHT, (parsed.percent / 100) * SDF_WEIGHT);
                lastPercent = Math.max(lastPercent, weightedPercent);
                
                // Forward with weighted percentage
                sseSend(res, {
                  ...parsed,
                  percent: weightedPercent,
                  message: parsed.message || `Progress: ${Math.round(weightedPercent)}%`
                });
              } else {
                sseSend(res, parsed);
              }
              return;
            }
            sseSend(res, parsed);
            return;
          }

          const isMolDone =
            /\[OK\]\s+molecule\s+\d+\/\d+\s+processed\s*:/i.test(s) ||
            /\[OK\]\s+molecule\s+[a-z0-9_-]+\s*:\s*sdf\s+generated/i.test(s) ||
            /sdf\s+(written|saved)\s*:\s*[\w/\\.-]+\.sdf/i.test(s) ||
            /generation\s+sdf\s+(finished|complete|done)/i.test(s) ||
            /successfully\s+(wrote|saved|generated).*\.sdf/i.test(s) ||
            /finished\s+processing\s+molecule/i.test(s) ||
            /saved\s+output.*\.sdf/i.test(s);

          if (isMolDone) {
            console.log(`[BACKEND-TRIGGER] [OK] End message detected: "${s}"`);
            console.log(
              `[BACKEND-BEFORE] doneMolecules=${doneMolecules}, lastPercent=${lastPercent.toFixed(2)}%`
            );

            doneMolecules = Math.min(doneMolecules + 1, totalMolecules);

            const frac = totalMolecules > 0 ? doneMolecules / totalMolecules : 1;
            let percent = Math.floor(frac * SDF_WEIGHT * 100) / 100;

            if (percent <= lastPercent) {
              percent = lastPercent + 0.1;
            }

            percent = Math.min(SDF_WEIGHT, percent);
            lastPercent = percent;

            console.log(
              `[BACKEND-AFTER] doneMolecules=${doneMolecules}/${totalMolecules}, percent=${percent.toFixed(2)}%, frac=${frac.toFixed(4)}`
            );

            sseSend(res, {
              type: 'progress',
              phase: 'prediction',
              current: doneMolecules,
              total: totalMolecules,
              percent: percent,
              message: `Molecule ${doneMolecules}/${totalMolecules} predicted (SDF)`,
            });

            return;
          }

          sseSend(res, { type: 'log', message: s });
        });
      });

      child.stderr.on('data', buf => {
        const msg = buf.toString();
        sseSend(res, { type: 'log', message: msg.replace(/[^\x00-\x7F]/g, '?') });
        if (/error|traceback/i.test(msg)) sseSend(res, { type: 'error', message: msg });
      });

      child.on('close', async (code, signal) => {
        // Retrieve the task BEFORE deleting it to check whether it was canceled
        const task = getTask(taskId);
        deleteTask(taskId);
        if (task?.safetyTimeout) clearTimeout(task.safetyTimeout);
        clearInterval(hb);
        const durationMin = Math.round((Date.now() - startTime) / 60000 * 10) / 10;

        console.log(
          `[BACKEND-CLOSE] code=${code}, signal=${signal}, doneMolecules=${doneMolecules}/${totalMolecules}`
        );

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
          sseSend(res, {
            type: 'completed',
            success: false,
            error: `Python exited with code ${code}`,
            duration: durationMin,
          });
          res.end();
          return;
        }

        const files = fs.readdirSync(workDir);
        const allSdfs = files.filter(
          f => f.toLowerCase().endsWith('.sdf') && f.toLowerCase() !== 'input.sdf'
        );
        const sdfOut = allSdfs[0] || null;
        const csvOut = files.find(f => f.toLowerCase().endsWith('.csv'));
        const jsonOut = files.find(f => f.toLowerCase().endsWith('.json'));
        const payload = { type: 'completed', success: true, taskId, duration: durationMin, outputs: {} };

        if (sdfOut) payload.outputs.sdf = `/outputs/${taskId}/${sdfOut}`;
        if (csvOut) payload.outputs.csv = `/outputs/${taskId}/${csvOut}`;
        if (jsonOut) payload.outputs.json = `/outputs/${taskId}/${jsonOut}`;

        let sdfForExtractors = null;
        if (allSdfs.length > 1) {
          sdfForExtractors = workDir;
        } else if (allSdfs.length === 1) {
          sdfForExtractors = path.join(workDir, allSdfs[0]);
        } else if (mode === 'SDF') {
          const candidate = path.join(workDir, 'input.sdf');
          if (fs.existsSync(candidate)) sdfForExtractors = candidate;
        }

        const wantExtract = options.run_extractors === true;
        if (wantExtract && (sdfForExtractors || mode === 'CSV')) {
          sseSend(res, { type: 'log', message: 'Post-processing: generating CSV 1H/13C/COSY/HSQC...' });

          // Progress for extractor phase (95% to 100%)
          const extractorMethods = selectedMethods.filter(m => ['1H', '13C', 'COSY', 'HSQC'].includes(m));
          const extractorCount = extractorMethods.length || 4;
          const percentPerExtractor = EXTRACTOR_WEIGHT / extractorCount;
          let extractorProgress = SDF_WEIGHT;
          let currentExtractorIndex = 0;

          sseSend(res, {
            type: 'progress',
            phase: 'extractors',
            percent: SDF_WEIGHT,
            message: 'Starting CSV extraction...',
          });

          try {
            const cm =
              options?.extractors?.c13?.carbon_mode ?? options?.carbon_mode_13c ?? options?.carbon_mode;
            if (cm) sseSend(res, { type: 'log', message: `13C carbon_mode=${cm}` });

            const extra = await runExtractorsOnSDF({
              sdfPath: sdfForExtractors || workDir,
              workDir,
              taskId,
              pythonExe: PYTHON,
              env: childEnv,
              sse: msg => {
                if (!msg) return;
                if (msg.type === 'progress') {
                  // Map extractor progress (0-100%) into global progress (80-100%)
                  const extractorLocalPercent = msg.percent || 0;
                  const globalPercent = SDF_WEIGHT + (currentExtractorIndex * percentPerExtractor) + 
                                        (extractorLocalPercent / 100 * percentPerExtractor);
                  sseSend(res, {
                    type: 'progress',
                    phase: 'extractors',
                    percent: Math.min(99, Math.round(globalPercent * 10) / 10),
                    message: msg.message || 'Extraction in progress...',
                  });
                  return;
                }
                if (msg.type === 'log' && msg.message) {
                  // Detect start of a new extractor
                  if (
                    msg.message.toLowerCase().includes('starting') ||
                    msg.message.toLowerCase().includes('started') ||
                    msg.message.toLowerCase().includes('start')
                  ) {
                    currentExtractorIndex++;
                  }
                }
                sseSend(res, msg);
              },
              extractorOptions: options.extractors || options || {},
              selectedMethods: selectedMethods,
            });
            Object.assign(payload.outputs, extra);
          } catch (e) {
            sseSend(res, { type: 'log', message: `[WARN] Post-processing failed: ${e.message}` });
          }
        }

        if (lastPercent < 100) {
          lastPercent = 100;
          sseSend(res, {
            type: 'progress',
            phase: 'extractors',
            percent: 100,
            message: 'Prediction and extractor CSV generation finished',
          });
          console.log('[BACKEND-FINAL] Forced 100%');
        }

        sseSend(res, payload);
        res.end();
      });

      child.on('error', err => {
        deleteTask(taskId);
        if (taskInfo.safetyTimeout) clearTimeout(taskInfo.safetyTimeout);
        clearInterval(hb);
        sseSend(res, { type: 'completed', success: false, error: `Python spawn error: ${err.message}` });
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
        }
      });
    })();
  } catch (e) {
    logWithMemory(`handlePrediction error: ${e.message}`, 'ERROR');
    if (res.headersSent) {
      sseSend(res, { type: 'error', message: `Internal error: ${e.message}` });
      sseSend(res, { type: 'completed', success: false, error: e.message });
      res.end();
    } else {
      res.status(500).json({ error: e.message });
    }
  }
}

module.exports = {
  handlePrediction,
  normalizePredictOptions,
  extractSmilesFromSdf,
};