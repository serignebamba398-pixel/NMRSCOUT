const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

// Exercise both request handlers through params.json creation, stopping before R.
for (const handler of ['handleAnnotation', 'handleAnnotationFromPrediction']) {
  for (const method of ['GET', 'POST']) {
    for (const annotationType of ['COSY', 'TOCSY']) {
      test(`${handler} ${method} forwards ${annotationType} to R`, async () => {
        const writes = [];
        const mocks = {
          fs: {
            existsSync: () => true,
            statSync: () => ({ size: 10 }),
            copyFileSync: () => {},
            readdirSync: () => [`prediction_${annotationType}_prediction.csv`],
            writeFileSync: (file, content) => writes.push({ file, content }),
          },
          uuid: { v4: () => 'annotation-test' },
          'async-mutex': { Semaphore: class { async acquire() { return [1, () => {}]; } } },
          '../config': { DEFAULT_LIGHT_PARAMS: { nb_cores: 1 } },
          '../config/paths': {
            uploadsDir: '/uploads',
            getActiveOutputsDir: () => '/outputs',
          },
          '../utils/filesystem': { ensureDir: () => {}, safeResolve: path.join },
          '../utils/spawn': { getRscript: () => { throw new Error('R launch reached'); } },
          '../utils/logger': { logWithMemory: () => {} },
          '../middleware/sse': { sseInit: () => null, sseSend: () => {} },
        };
        const moduleObject = { exports: {} };
        const source = fs.readFileSync(path.join(__dirname, '../services/annotationService.js'), 'utf8');
        vm.runInNewContext(source, {
          module: moduleObject,
          require: name => Object.hasOwn(mocks, name) ? mocks[name]
            : name.startsWith('.') ? {} : require(name),
        });
        const values = {
          annotationType, predictionTaskId: 'prediction',
          tableauFile: 'tableau.csv', matrixFile: 'matrix.csv',
          // Request type must win over stale saved parameters.
          parameters: { annotation_type: annotationType === 'TOCSY' ? 'COSY' : 'TOCSY' },
        };
        if (method === 'GET') values.parameters = JSON.stringify(values.parameters);
        let response;
        const res = { status() { return this; }, json(value) { response = value; } };
        await moduleObject.exports[handler]({ body: values, query: values, headers: {} }, res, method);
        assert.equal(response?.error, 'R launch reached');
        const paramsFiles = writes.filter(({ file }) => path.basename(file) === 'params.json');
        assert.ok(paramsFiles.length > 0);
        for (const { content } of paramsFiles) {
          assert.equal(JSON.parse(content).annotation_type, annotationType);
        }
      });
    }
  }
}
