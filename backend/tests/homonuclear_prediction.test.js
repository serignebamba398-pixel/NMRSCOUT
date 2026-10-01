const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const { EventEmitter } = require('node:events');
const helpers = require('../utils/homonuclearPrediction');

function loadService(name, mocks) {
  const result = { exports: {} };
  vm.runInNewContext(fs.readFileSync(path.join(__dirname, '../services', name), 'utf8'), {
    module: result, process, console,
    require: key => key in mocks ? mocks[key] : key.startsWith('.') ? {} : require(key),
  });
  return result.exports;
}

for (const bundled of [false, true]) {
  for (const folder of [false, true]) {
    test(`1H uses its own CLI (bundled=${bundled}, folder=${folder})`, async () => {
      const workDir = path.resolve('test-output');
      const sdf = path.join(workDir, 'molecule.sdf');
      const script = path.resolve('backend/scripts/predict_1H_real.py');
      const outputCsv = path.join(workDir, 'prediction_1H_test.csv');
      const files = new Set([sdf, script]);
      const calls = [];
      const service = loadService('extractorService.js', {
        '../config/paths': { EXTRACTORS: { '1H': script } },
        '../utils/homonuclearPrediction': helpers,
        fs: {
          existsSync: name => name === workDir || files.has(name),
          statSync: name => ({ isFile: () => files.has(name), isDirectory: () => name === workDir, size: 100 }),
          readdirSync: dir => [...files].filter(name => path.dirname(name) === dir).map(name => path.basename(name)),
        },
        child_process: { spawn(executable, args) {
          calls.push(args);
          const child = new EventEmitter();
          child.stdout = new EventEmitter(); child.stderr = new EventEmitter();
          queueMicrotask(() => {
            // Model the required output argument: incompatible calls cannot produce a CSV.
            const index = args.indexOf('--output_csv');
            if (index >= 0) files.add(args[index + 1]);
            child.emit('close', index >= 0 ? 0 : 2);
          });
          return child;
        } },
      });
      const output = await service.runExtractorsOnSDF({
        sdfPath: folder ? workDir : sdf, workDir, taskId: 'test',
        pythonExe: bundled ? 'nmr_engine_production.exe' : '/opt/conda/envs/nmr/bin/python',
        env: {}, sse() {}, selectedMethods: ['1H'],
        extractorOptions: { h1: { keep_exchangeable: true } },
      });
      assert.equal(calls.length, 1);
      assert.deepEqual(Array.from(calls[0]), [
        bundled ? 'predict_1H' : script,
        folder ? '--input_folder' : '--sdf_path', folder ? workDir : sdf,
        '--output_csv', outputCsv, '--keep_exchangeable',
      ]);
      assert.equal(output.csv_1H, '/outputs/test-output/prediction_1H_test.csv');
    });
  }
}

test('request validation preserves independent TOCSY bounds without mixing time', () => {
  const service = loadService('predictionService.js', {
    '../utils/homonuclearPrediction': helpers,
    '../utils/logger': { logWithMemory() {} },
  });
  const normalized = service.normalizePredictOptions({ extractors: {
    cosy: { min_bonds: 2, max_bonds: 3 },
    tocsy: { min_bonds: 4, max_bonds: 6, mixing_time: '80', include_diagonal: false },
  } });
  assert.deepEqual(normalized.extractors.tocsy, { min_bonds: 4, max_bonds: 6, include_diagonal: false });
  assert.equal(normalized.extractors.cosy.max_bonds, 3);
  for (const config of [{ max_bonds: 0 }, { max_bonds: 2.5 }, { min_bonds: 7, max_bonds: 6 },
    { include_diagonal: 'no' }]) {
    let status, body;
    service.handlePrediction({ body: { options: { extractors: { tocsy: config } } } }, {
      status(value) { status = value; return this; }, json(value) { body = value; },
    });
    assert.equal(status, 400);
    assert.match(body.error, /TOCSY:/);
  }
  assert.equal(helpers.normalizeHomonuclearOptions('TOCSY', { mixing_time: '' }).mixing_time, undefined);
});

for (const pythonExe of ['/opt/conda/envs/nmr/bin/python', 'C:\\app\\nmr_engine_production.exe']) {
  test(`extractor forwards COSY and TOCSY independently: ${pythonExe}`, async () => {
    const workDir = path.resolve('test-output');
    const sdf = path.join(workDir, 'molecule.sdf');
    const script = path.resolve('backend/scripts/predict_COSY_real.py');
    const files = new Set([sdf, script]);
    const calls = [];
    const service = loadService('extractorService.js', {
      '../config/paths': { EXTRACTORS: { COSY: script } },
      '../utils/homonuclearPrediction': helpers,
      fs: {
        existsSync: name => name === workDir || files.has(name),
        statSync: name => ({ isFile: () => files.has(name), isDirectory: () => name === workDir, size: 100 }),
        readdirSync: dir => [...files].filter(name => path.dirname(name) === dir).map(name => path.basename(name)),
      },
      child_process: { spawn(executable, args) {
        calls.push({ executable, args });
        const child = new EventEmitter();
        child.stdout = new EventEmitter(); child.stderr = new EventEmitter();
        queueMicrotask(() => {
          const experiment = args[args.indexOf('--experiment') + 1];
          files.add(path.join(workDir, `prediction_${experiment}_test.csv`));
          child.stdout.emit('data', Buffer.from('{"type":"progress","percent":100}\n'));
          child.emit('close', 0);
        });
        return child;
      } },
    });
    const output = await service.runExtractorsOnSDF({
      sdfPath: sdf, workDir, taskId: 'test', pythonExe, env: {}, sse() {},
      selectedMethods: ['COSY', 'TOCSY'],
      extractorOptions: { cosy: { min_bonds: 2, max_bonds: 3 },
        tocsy: { min_bonds: 4, max_bonds: 6, include_diagonal: false } },
    });
    assert.equal(calls.length, 2);
    for (const [i, experiment] of ['COSY', 'TOCSY'].entries()) {
      const args = calls[i].args;
      assert.equal(args[0], pythonExe.endsWith('.exe') ? 'predict_COSY' : script);
      assert.equal(args[args.indexOf('--experiment') + 1], experiment);
      assert.equal(args[args.indexOf('--max_bonds') + 1], i ? '6' : '3');
      assert.equal(args[args.indexOf('--min_bonds') + 1], i ? '4' : '2');
      assert.ok(args.includes(i ? '--no_diagonal' : '--diagonal'));
      assert.ok(!args.includes('--mixing_time'));
      assert.match(output[`csv_${experiment}`], new RegExp(`prediction_${experiment}_test.csv$`));
    }
  });
}
