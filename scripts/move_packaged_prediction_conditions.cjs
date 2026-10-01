// Apply the existing source layout to an older packaged UI without rebuilding it.
const fs = require('fs');
const path = require('path');
const crypto = require('crypto');
const assert = require('assert');
const asar = require('@electron/asar');
const { Pickle } = require('../node_modules/@electron/asar/lib/pickle');
const parser = require('../frontend/node_modules/@babel/parser');

const root = path.resolve(__dirname, '..');
const archive = path.join(root, 'dist/win-unpacked/resources/app.asar');
const raw = asar.getRawHeader(archive);
const original = fs.readFileSync(archive);
const jsFiles = asar.listPackage(archive).filter(p => /frontend[\\/]build[\\/]static[\\/]js[\\/]main\.[^\\/]+\.js$/.test(p));
assert.equal(jsFiles.length, 1);
const jsPath = jsFiles[0].replace(/^[\\/]/, '');
const source = asar.extractFile(archive, jsPath).toString();
const ast = parser.parse(source);
function pathsForString(value) {
  const results = [];
  function walk(node, parents = []) {
    if (!node || typeof node !== 'object') return;
    if (node.type === 'StringLiteral' && node.value === value) results.push([...parents, node]);
    for (const [key, child] of Object.entries(node)) {
      if (['loc', 'start', 'end', 'extra'].includes(key)) continue;
      if (Array.isArray(child)) child.forEach(n => walk(n, [...parents, node]));
      else if (child && child.type) walk(child, [...parents, node]);
    }
  }
  walk(ast);
  assert.equal(results.length, 1, `Expected unique ${value}`);
  return results[0];
}
const prop = (obj, key) => obj?.properties?.find(p => p.key?.name === key)?.value;
const expPath = pathsForString('Experimental Conditions');
const section = [...expPath].reverse().find(n => n.type === 'CallExpression' && prop(n.arguments[1], 'className')?.value === 'settings-section');
assert(section);
const owner = [...expPath].reverse().find(n => n.type === 'ArrowFunctionExpression');
const returnNode = owner.body.body.find(n => n.type === 'ReturnStatement');
assert(returnNode && owner.body.type === 'BlockStatement');
const declarations = source.slice(owner.body.start + 1, returnNode.start);
assert(declarations.includes('freq_mhz') && declarations.includes('temperature_k'));
const moved = `(()=>{${declarations}return ${source.slice(section.start, section.end)}})()`;
const methodsPath = pathsForString('Prediction Methods');
const sourceTab = [...methodsPath].reverse().find(n => n.type === 'CallExpression' && prop(n.arguments[1], 'className')?.value === 'tab-content');
const children = prop(sourceTab.arguments[1], 'children');
assert.equal(children.type, 'ArrayExpression');
const summary = children.elements.find(n => n?.type === 'CallExpression' && prop(n.arguments[1], 'className')?.value === 'methods-summary');
assert(summary && summary.start < section.start);
const jsxCallee = section.callee.expressions?.at(-1);
assert.equal(jsxCallee?.type, 'MemberExpression');
const runtimeName = jsxCallee.object.name;
const hint = `(0,${runtimeName}.jsx)("p",{className:"hint",children:"Experimental conditions are configured in Source & Methods."})`;
let updated = source.slice(0, section.start) + hint + source.slice(section.end);
updated = updated.slice(0, summary.start) + moved + ',' + updated.slice(summary.start);
// The original source map describes the old layout, so do not advertise it.
updated = updated.replace(/\n?\/\/# sourceMappingURL=[^\r\n]*/g, '');
parser.parse(updated);

// Exercise the moved controls with the same captured state and event handlers.
const changes = [];
const jsx = (type, props) => ({ type, props });
const render = new Function('R', 'Z', 'h', runtimeName, `return ${moved}`);
const settings = { freq_mhz: 800, solvent: 'D2O', temperature_k: 303, ph: 6.5 };
const tree = render(settings, (...args) => changes.push(args), false, { jsx, jsxs: jsx });
function controls(node, out = []) {
  if (!node || typeof node !== 'object') return out;
  if (['input', 'select'].includes(node.type)) out.push(node);
  const child = node.props?.children;
  (Array.isArray(child) ? child : [child]).forEach(n => controls(n, out));
  return out;
}
const fields = controls(tree);
assert.deepEqual(fields.map(n => n.props.value), [800, 'D2O', 303, 6.5]);
fields.forEach((n, i) => n.props.onChange({ target: { value: ['600', 'CDCl3', '298', '7.4'][i] } }));
assert.deepEqual(changes, [['freq_mhz', 600], ['solvent', 'CDCl3'], ['temperature_k', 298], ['ph', 7.4]]);
assert(controls(render(settings, () => {}, true, { jsx, jsxs: jsx })).every(n => n.props.disabled));

const content = Buffer.from(updated);
const hash = b => crypto.createHash('sha256').update(b).digest('hex');
let entry = raw.header;
for (const segment of jsPath.split(/[\\/]/)) entry = entry.files[segment];
assert(!entry.unpacked);
const payload = original.subarray(8 + raw.headerSize);
entry.offset = String(payload.length);
entry.size = content.length;
const blockSize = 4 * 1024 * 1024;
const blocks = [];
for (let i = 0; i < content.length; i += blockSize) blocks.push(hash(content.subarray(i, i + blockSize)));
entry.integrity = { algorithm: 'SHA256', hash: hash(content), blockSize, blocks };
const header = Pickle.createEmpty();
header.writeString(JSON.stringify(raw.header));
const size = Pickle.createEmpty();
size.writeUInt32(header.toBuffer().length);
const patchedArchive = Buffer.concat([size.toBuffer(), header.toBuffer(), payload, content]);
const backupDir = path.join(root, 'backups', `prediction_conditions_${Date.now()}`);
fs.mkdirSync(backupDir, { recursive: true });
fs.writeFileSync(path.join(backupDir, 'app.asar'), original);
const staged = archive + '.conditions.tmp';
fs.writeFileSync(staged, patchedArchive);
assert.equal(asar.extractFile(staged, jsPath).toString(), updated);
// All packed files other than the changed JS retain their exact byte contents.
let verifiedFiles = 0;
function verify(node, prefix = '') {
  for (const [name, info] of Object.entries(node.files || {})) {
    const filename = path.join(prefix, name);
    if (info.files) verify(info, filename);
    else if (!info.unpacked && !info.link && filename !== jsPath) {
      assert(asar.extractFile(archive, filename).equals(asar.extractFile(staged, filename)), filename);
      verifiedFiles++;
    }
  }
}
verify(raw.header);
assert.equal(hash(fs.readFileSync(archive)), hash(original), 'Archive changed while preparing patch');
const report = { archive, jsPath, backupDir, verifiedUnchangedPackedFiles: verifiedFiles,
  controls: 'PASS: stored values, change handlers, execution lock', syntax: 'PASS', sha256: hash(patchedArchive),
  originalSha256: hash(original), status: 'ready' };
fs.writeFileSync(path.join(backupDir, 'report.json'), JSON.stringify(report, null, 2));
fs.renameSync(staged, archive);
asar.uncacheAll();
assert.equal(hash(asar.extractFile(archive, jsPath)), hash(content));
report.status = 'applied';
fs.writeFileSync(path.join(backupDir, 'report.json'), JSON.stringify(report, null, 2));
console.log(JSON.stringify(report, null, 2));
