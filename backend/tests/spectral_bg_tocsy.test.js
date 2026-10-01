const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const { EventEmitter } = require('node:events');

function loadRoute(python, missingTocsy) {
  const handlers = {}, calls = [];
  const mocks = {
    express: { Router: () => ({ get() {}, post: (url, handler) => { handlers[url] = handler; } }) },
    fs: {
      existsSync: file => !(missingTocsy && file === path.join('/tocsy', '2rr')),
      statSync: () => ({ size: 100 }), readFileSync: () => '{}',
    },
    child_process: { spawn: (exe, args) => {
      calls.push({ exe, args });
      const child = new EventEmitter();
      child.stdout = new EventEmitter(); child.stderr = new EventEmitter();
      queueMicrotask(() => child.emit('close', 0));
      return child;
    } },
    uuid: { v4: () => 'test-export' },
    '../config/paths': { BG_SCRIPT: '/scripts/nmr_export_bg.py', getActiveOutputsDir: () => '/outputs' },
    '../utils/filesystem': { ensureDir() {} },
    '../utils/pathNormalizer': { normalizeHostPathToContainer: value => value },
    '../utils/logger': { logWithMemory() {} },
    '../utils/spawn': { getPythonExe: () => python, buildPythonEnv: () => ({}) },
    '../services/backgroundService': {},
    '../constants/defaults': require('../constants/defaults'),
  };
  vm.runInNewContext(fs.readFileSync(path.join(__dirname, '../routes/spectralBg.js'), 'utf8'), {
    require: name => mocks[name] || require(name), module: { exports: {} }, console,
  });
  return { handler: handlers['/api/spectral-bg/export'], calls };
}

for (const python of ['python', 'nmr_engine_production.exe']) {
  for (const missing of [false, true]) {
    test(`${python}: TOCSY export ${missing ? 'missing 2rr does not block COSY' : 'uses separate files and bounds'}`, async () => {
      const { handler, calls } = loadRoute(python, missing);
      let result, status = 200;
      await handler({ body: {
        pdata: { TOCSY: '/tocsy', COSY: '/cosy' },
        options: { tocsy: { f1_ppm_min: 2, f1_ppm_max: 4, f2_ppm_min: 1, f2_ppm_max: 5 } },
      } }, { status(code) { status = code; return this; }, json(data) { result = data; } });
      assert.equal(status, 200);
      assert.equal(result.outputs.COSY.ok, true);
      assert.equal(result.outputs.TOCSY.ok, !missing);
      if (missing) {
        assert.match(result.outputs.TOCSY.error, /missing 2rr/);
        assert.equal(calls.length, 1);
      } else {
        const args = calls[0].args;
        assert.equal(args[args.indexOf('--method') + 1], 'TOCSY');
        assert.equal(args[args.indexOf('--f1_ppm_min') + 1], '2');
        assert.equal(args[args.indexOf('--f2_ppm_max') + 1], '5');
        assert.match(args[args.indexOf('--out_png') + 1], /tocsy_bg\.png$/);
      }
    });
  }
}
