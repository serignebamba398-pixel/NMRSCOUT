const fs = require('fs');
const path = require('path');
const assert = require('assert');
const asar = require('@electron/asar');
const parser = require('../frontend/node_modules/@babel/parser');
const traverse = require('../frontend/node_modules/@babel/traverse').default;
const babel = require('../frontend/node_modules/@babel/core');
const updatePackagedJs = require('./update_packaged_js.cjs');
const root = path.resolve(__dirname, '..');
const archive = path.join(root, 'dist/win-unpacked/resources/app.asar');
const filename = asar.listPackage(archive).find(p => /main\.[^\\/]+\.js$/.test(p)).replace(/^[\\/]/, '');
const original = asar.extractFile(archive, filename).toString();
assert(!original.includes('function MoleculeScoreSummary('), 'Molecule summary is already installed');
const ast = parser.parse(original);
const property = (obj, key) => obj.properties.find(p => (p.key.name || p.key.value) === key)?.value;
let cardPath;
traverse(ast, { ObjectProperty(p) {
  if (p.node.key.value === 'data-mol-source') {
    assert(!cardPath, 'Multiple molecule grids');
    cardPath = p;
  }
} });
assert(cardPath);
const cardProps = cardPath.parentPath.node;
const sourceName = cardPath.node.value.name;
const elements = property(cardProps, 'children').elements;
const oldLabel = elements.at(-1);
assert.equal(oldLabel.type, 'CallExpression');
const modal = cardPath.findParent(p => p.isFunctionDeclaration());
let scoresName;
let rowsCall;
modal.traverse({ ObjectProperty(p) {
  if (p.node.key.name === 'sourcesWithScore') scoresName = p.node.value.left?.name || p.node.value.name;
  if (p.node.key.name === 'rowsByType' && p.node.value.type === 'CallExpression'
      && p.node.value.arguments[0]?.name === sourceName) rowsCall = original.slice(p.node.value.start, p.node.value.end);
} });
assert(scoresName && rowsCall);
const runtime = oldLabel.callee.expressions.at(-1).object.name;
let component = fs.readFileSync(path.join(root, 'frontend/src/components/viewer/MoleculeScoreSummary.js'), 'utf8');
component = component.replace(/^import React from 'react';\s*/m, '').replace(/export default function /g, 'function ').replace(/export function /g, 'function ');
const compiled = babel.transformSync(component, { configFile: false, babelrc: false,
  presets: [[require.resolve('../frontend/node_modules/@babel/preset-react'), { runtime: 'classic', pragma: '__jsx' }]] }).code;
const replacement = `((__jsx)=>{${compiled};return MoleculeScoreSummary({source:${sourceName},rowsByType:${rowsCall},averageScore:(${scoresName}.find(item=>item.source===${sourceName})||{}).avgCoverage})})((tag,props,...children)=>(0,${runtime}.jsxs)(tag,{...props,children:children.length===1?children[0]:children}))`;
const updated = (original.slice(0, oldLabel.start) + replacement + original.slice(oldLabel.end))
  .replace(/\n?\/\/# sourceMappingURL=[^\r\n]*/g, '');
parser.parse(updated);
// Check the exact injected expression using realistic rows and ranking scores.
const evaluate = new Function(sourceName, scoresName, rowsCall.split('(')[0], runtime, `return ${replacement}`);
const jsx = (type, props) => ({ type, props });
const rendered = evaluate('NPA020091', [{ source: 'NPA020091', avgCoverage: 98 }],
  () => ({ '1H': { coverage_sub: 100, avg_variabilite: 0.0452 }, '13C': { coverage_sub: 80, avg_variabilite: 0.0824 } }), { jsxs: jsx });
const result = JSON.stringify(rendered);
for (const expected of ['NPA020091', '98.0', '90.00%', '0.0638']) assert(result.includes(expected), expected);
updatePackagedJs(archive, filename, updated, 'molecule_scores');
