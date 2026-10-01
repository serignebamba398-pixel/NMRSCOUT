const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const { EXTRACTORS } = require('../config/paths');
const { homonuclearArgs } = require('../utils/homonuclearPrediction');

function runExtractorsOnSDF({
  sdfPath,
  workDir,
  taskId,
  pythonExe,
  env,
  sse,
  extractorOptions = {},
  selectedMethods = ['1H', '13C', 'COSY', 'HSQC'],
}) {
  // Validate before entering the async executor so invalid options cannot hang a job.
  const homonuclear = {};
  for (const experiment of ['COSY', 'TOCSY']) {
    if (selectedMethods.includes(experiment)) {
      homonuclear[experiment] = homonuclearArgs(experiment, {
        min_bonds: extractorOptions.min_bonds,
        max_bonds: extractorOptions.max_bonds,
        include_diagonal: extractorOptions.include_diagonal,
        ...(extractorOptions[experiment.toLowerCase()] || {}),
      });
    }
  }
  return new Promise(async resolve => {
    const outputs = {};

    const findSdfContainer = () => {
      const listSdf = dir => {
        try {
          return fs
            .readdirSync(dir)
            .filter(f => f.toLowerCase().endsWith('.sdf') && f.toLowerCase() !== 'input.sdf');
        } catch {
          return [];
        }
      };
      const subdirs = dir => {
        try {
          return fs
            .readdirSync(dir)
            .map(n => path.join(dir, n))
            .filter(p => fs.existsSync(p) && fs.statSync(p).isDirectory());
        } catch {
          return [];
        }
      };

      if (
        sdfPath &&
        fs.existsSync(sdfPath) &&
        fs.statSync(sdfPath).isFile() &&
        sdfPath.toLowerCase().endsWith('.sdf')
      ) {
        const dir = path.dirname(sdfPath);
        const here = listSdf(dir);
        if (here.length > 1) return { mode: 'folder', folder: dir };
        const sd = subdirs(dir).find(
          d => path.basename(d).toLowerCase().startsWith('sdf') || listSdf(d).length > 0
        );
        if (sd) return { mode: 'folder', folder: sd };
        return { mode: 'single', single: sdfPath };
      }

      if (sdfPath && fs.existsSync(sdfPath) && fs.statSync(sdfPath).isDirectory()) {
        return { mode: 'folder', folder: sdfPath };
      }

      const here = listSdf(workDir);
      if (here.length > 1) return { mode: 'folder', folder: workDir };
      if (here.length === 1) return { mode: 'single', single: path.join(workDir, here[0]) };

      const sd = subdirs(workDir).find(
        d => path.basename(d).toLowerCase().startsWith('sdf') || listSdf(d).length > 0
      );
      if (sd) return { mode: 'folder', folder: sd };

      return { mode: 'folder', folder: workDir };
    };

    const container = findSdfContainer();
    const logContainer =
      container.mode === 'single'
        ? `single file: ${path.basename(container.single)}`
        : `folder: ${path.relative(workDir, container.folder) || '.'}`;
    sse({ type: 'log', message: `Post-processing — SDF detected (${logContainer})` });
    sse({ type: 'log', message: `Methods to run: ${selectedMethods.join(', ')}` });

    try {
      const _hsqcFold =
        extractorOptions?.hsqc?.fold ?? extractorOptions?.hsqc_fold ?? extractorOptions?.fold_hsqc;
      const _hsqcFoldLow = extractorOptions?.hsqc?.fold_low ?? extractorOptions?.hsqc_fold_low;
      const _hsqcFoldHigh = extractorOptions?.hsqc?.fold_high ?? extractorOptions?.hsqc_fold_high;
      if (_hsqcFold !== undefined || _hsqcFoldLow !== undefined || _hsqcFoldHigh !== undefined) {
        sse({
          type: 'log',
          message: `HSQC options → fold=${_hsqcFold ?? '(def:true)'} low=${_hsqcFoldLow ?? '(def:4.8618)'} high=${_hsqcFoldHigh ?? '(def:69.8395)'}`,
        });
      }
    } catch {}

    const ensureOut = (key, absPath, token) => {
      try {
        const tok = String(token).toLowerCase();
        const candidates = [absPath, path.join(workDir, path.basename(absPath))];

        const anyCsv = fs
          .readdirSync(workDir)
          .filter(f => f.toLowerCase().endsWith('.csv') && f.toLowerCase().includes(tok))
          .map(f => path.join(workDir, f));

        for (const p of [...candidates, ...anyCsv]) {
          try {
            if (fs.existsSync(p) && fs.statSync(p).size > 0) {
              outputs[key] = `/outputs/${path.basename(workDir)}/${path.basename(p)}`;
              sse({ type: 'log', message: `[OK] ${path.basename(p)} detected` });
              return true;
            }
          } catch {}
        }
      } catch {}
      sse({ type: 'log', message: `[WARN] CSV ${token} not found` });
      return false;
    };

    const spawnOnce = (title, script, args) =>
      new Promise(res => {
        const isBundledEnv = process.env.IS_ELECTRON_BUNDLED === 'true';
        const isBundledExe = pythonExe.toLowerCase().endsWith('nmr_engine_production.exe') || 
                             pythonExe.toLowerCase().endsWith('nmr_engine_production') ||
                             pythonExe.toLowerCase().endsWith('nmr_engine.exe') || 
                             pythonExe.toLowerCase().endsWith('nmr_engine');
        const useCommandMode = isBundledExe;

        if (!useCommandMode && !fs.existsSync(script)) {
          sse({ type: 'log', message: `[WARN] ${title}: script not found (${path.basename(script)})` });
          return res({ code: -1 });
        }

        sse({ type: 'log', message: `[START] ${title} starting...` });
        sse({ type: 'progress', percent: 0, message: `Extraction ${title}...` });

        let spawnArgs;
        if (useCommandMode) {
            // Map script path to command name
            const scriptName = path.basename(script);
            let commandName;
            if (scriptName === 'predict_1H_real.py') commandName = 'predict_1H';
            else if (scriptName === 'predict_13C_real.py') commandName = 'predict_13C';
            else if (scriptName === 'predict_COSY_real.py') commandName = 'predict_COSY';
            else if (scriptName === 'predict_HSQC_real.py') commandName = 'predict_HSQC';
            else commandName = scriptName.replace('.py', ''); // Fallback
            
            spawnArgs = [commandName, ...args];
        } else {
            spawnArgs = [script, ...args];
        }

        const p = spawn(pythonExe, spawnArgs, {
          cwd: workDir,
          stdio: ['ignore', 'pipe', 'pipe'],
          shell: false,
          env,
        });

        p.on('error', (err) => {
          sse({ type: 'error', message: `Failed to start extractor ${title}: ${err.message}` });
          res({ code: -1, error: err });
        });

        let buf = '';
        let lineCount = 0;
        const estimatedLines = 150;
        let lastReportedPercent = 0;

        p.stdout.on('data', d => {
          buf += d.toString();
          const lines = buf.split('\n');
          buf = lines.pop() || '';

          lines.forEach(l => {
            if (!l || !l.trim()) return;
            lineCount++;

            // Try to parse JSON progress messages from Python scripts
            let jsonProgress = null;
            try {
              if (l.startsWith('{') && l.includes('"type"')) {
                const parsed = JSON.parse(l);
                if (parsed.type === 'progress' && typeof parsed.percent === 'number') {
                  jsonProgress = parsed;
                }
              }
            } catch (e) {
              // Not valid JSON, continue with line-based estimation
            }

            if (jsonProgress) {
              // Use the exact progress from Python script
              const currentPercent = Math.min(95, Math.floor(jsonProgress.percent));
              if (currentPercent > lastReportedPercent) {
                sse({
                  type: 'progress',
                  percent: currentPercent,
                  message: jsonProgress.message || l.substring(0, 100),
                });
                lastReportedPercent = currentPercent;
              }
              // Also log the progress message
              sse({ type: 'log', message: `[Progress ${jsonProgress.percent}%] ${jsonProgress.message || ''}` });
            } else {
              // Fallback to line-count based estimation
              const rawPercent = Math.min(95, (lineCount / estimatedLines) * 100);
              const currentPercent = Math.floor(rawPercent);

              if (currentPercent > lastReportedPercent) {
                sse({
                  type: 'progress',
                  percent: currentPercent,
                  message: l.substring(0, 100),
                });
                lastReportedPercent = currentPercent;
              }

              sse({ type: 'log', message: l });
            }
          });
        });

        p.stderr.on('data', d => {
          const errMsg = d.toString();
          sse({ type: 'log', message: errMsg });
        });

        p.on('close', code => {
          sse({ type: 'progress', percent: 100, message: `${title} completed` });
          sse({ type: 'log', message: `[OK] ${title} exit code: ${code}` });
          res({ code });
        });

        p.on('error', e => {
          sse({ type: 'log', message: `[ERROR] ${title} error: ${e.message}` });
          res({ code: -2 });
        });
      });

    const buildArgVariantsForSingle = single => [
      ['--sdf_path', single],
      ['--sdf', single],
      ['--input_sdf', single],
      ['--input', single],
      ['-i', single],
    ];

    const buildArgVariantsForFolder = folder => [
      ['--sdf_path', folder],
      ['--input_folder', folder],
      ['--input_dir', folder],
      ['--folder', folder],
      ['--dir', folder],
      ['-i', folder],
    ];

    const buildBaseArgs = () => {
      if (container.mode === 'single') return buildArgVariantsForSingle(container.single);
      return buildArgVariantsForFolder(container.folder);
    };

    const extraArgs13C = () => {
      const mode =
        extractorOptions?.carbon_mode_13c ??
        extractorOptions?.c13?.carbon_mode ??
        extractorOptions?.carbon_mode ??
        'nonquaternary';
      if (mode !== 'all' && mode !== 'nonquaternary') return [];
      return ['--carbon_mode', mode];
    };

    const extraArgs1H = () => {
      const ex = extractorOptions || {};
      const h1 = ex['1H'] || ex.h1 || {};
      const keepExchangeable = h1.keep_exchangeable ?? ex.keep_exchangeable ?? false;
      
      const args = [];
      if (keepExchangeable === true || String(keepExchangeable).toLowerCase() === 'true') {
        args.push('--keep_exchangeable');
      }
      return args;
    };

    const extraArgsHSQC = () => {
      const ex = extractorOptions || {};
      const h = ex.hsqc || ex.HSQC || {};
      const fold = h.fold ?? ex.hsqc_fold ?? ex.fold_hsqc ?? ex.hsqcFold;
      const foldLow = h.fold_low ?? ex.hsqc_fold_low ?? ex.hsqcFoldLow;
      const foldHigh = h.fold_high ?? ex.hsqc_fold_high ?? ex.hsqcFoldHigh;

      const args = [];
      if (fold === false || String(fold).toLowerCase() === 'false' || fold === 0) {
        args.push('--no-fold');
      }
      if (foldLow != null && foldLow !== '') args.push('--fold-low', String(foldLow));
      if (foldHigh != null && foldHigh !== '') args.push('--fold-high', String(foldHigh));
      return args;
    };

    const extraArgsCOSY = experiment => homonuclear[experiment];

    const jobs = [];
    const planExtractor = (label, script, outKey, outFilename, extraArgBuilder = () => [], options = {}) => {
      if (!selectedMethods.includes(label)) {
        sse({ type: 'log', message: `${label} not selected, skipping` });
        return;
      }

      const outAbs = path.join(workDir, outFilename);
      const baseArgSets = options.outputCsv
        ? [container.mode === 'single'
          ? ['--sdf_path', container.single]
          : ['--input_folder', container.folder]]
        : buildBaseArgs();

      const makeOutVariants = base => {
        if (options.outputCsv) {
          return [[...base, '--output_csv', outAbs, ...extraArgBuilder()]];
        }
        const canonical = [
          ...base,
          '--output_dir',
          workDir,
          '--task_id',
          taskId,
          ...extraArgBuilder(),
        ];
        const variants = [canonical];

        if (!options.canonicalOnly) {
          variants.push([
            ...base,
            '--input_dir',
            container.folder || workDir,
            '--output_dir',
            workDir,
            '--task_id',
            taskId,
            ...extraArgBuilder(),
          ]);
        }

        return variants;
      };

      jobs.push(async () => {
        for (const base of baseArgSets) {
          const argSets = makeOutVariants(base);
          for (let i = 0; i < argSets.length; i++) {
            const args = argSets[i];
            const title = `${label}${i === 0 ? '' : ` (alt${i})`}`;
            await spawnOnce(title, script, args);
            if (ensureOut(outKey, outAbs, label)) return;
          }
        }
        if (container.mode === 'folder') {
          sse({ type: 'log', message: `[WARN] ${label}: no CSV produced after all variants.` });
        }
      });
    };

    planExtractor('1H', EXTRACTORS['1H'], 'csv_1H', `prediction_1H_${taskId}.csv`, extraArgs1H, { outputCsv: true });
    planExtractor('13C', EXTRACTORS['13C'], 'csv_13C', `prediction_13C_${taskId}.csv`, extraArgs13C);
    planExtractor(
      'COSY',
      EXTRACTORS['COSY'],
      'csv_COSY',
      `prediction_COSY_${taskId}.csv`,
      () => extraArgsCOSY('COSY'),
      { canonicalOnly: true },
    );
    planExtractor(
      'TOCSY',
      EXTRACTORS['COSY'],
      'csv_TOCSY',
      `prediction_TOCSY_${taskId}.csv`,
      () => extraArgsCOSY('TOCSY'),
      { canonicalOnly: true },
    );
    planExtractor('HSQC', EXTRACTORS['HSQC'], 'csv_HSQC', `prediction_HSQC_${taskId}.csv`, extraArgsHSQC);

    const totalJobs = jobs.length;
    let completedJobs = 0;
    
    for (const job of jobs) {
      await job();
      completedJobs++;
      // Emit overall extractor progress
      const jobProgress = Math.round((completedJobs / totalJobs) * 100);
      sse({ 
        type: 'progress', 
        percent: jobProgress, 
        message: `Extractors: ${completedJobs}/${totalJobs} completed` 
      });
    }

    try {
      const listing = fs.readdirSync(workDir).filter(f => f.toLowerCase().endsWith('.csv'));
      sse({ type: 'log', message: `CSV found: ${listing.join(', ') || '(none)'}` });
    } catch {}

    resolve(outputs);
  });
}

module.exports = { runExtractorsOnSDF };
