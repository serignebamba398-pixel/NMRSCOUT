import { tocsyHighlights } from './tocsyHighlights';
import { findBondsInPath } from './molParser';

jest.mock('../api', () => ({ API_BASE_URL: '' }));

function molecule(symbols, edges) {
  const field = n => String(n).padStart(3);
  return ['test', '', '', field(symbols.length) + field(edges.length) + '  0  0  0  0            999 V2000',
    ...symbols.map(s => ' '.repeat(31) + s.padEnd(3)),
    ...edges.map(([a, b]) => field(a + 1) + field(b + 1) + '  1  0  0  0'), 'M  END'].join('\n');
}

const mol = molecule(['H', 'C', 'C', 'C', 'C', 'C', 'H', 'H'],
  [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [1, 7]]);

test('TOCSY highlights all six bonds between distant protons in both renderers', () => {
  const result = tocsyHighlights(mol, [], { matched_atom_indices: '0:6||6:0' });
  expect(result.bonds).toEqual([0, 1, 2, 3, 4, 5]);
  expect([...result.stickAtoms].sort()).toEqual([0, 1, 2, 3, 4, 5, 6]);
  expect([...result.crossCarbons]).toEqual([1, 2, 3, 4, 5]);
  expect(result.protons).toEqual(new Set([0, 6]));
  // Preserve the existing COSY display limit.
  expect(findBondsInPath(mol, [[0, 6]])).toEqual([]);
  expect(findBondsInPath(mol, [[0, 2]])).toEqual([0, 1]);
});

test('equivalent diagonal protons color each C-H bond without cross-peak carbon spheres', () => {
  const result = tocsyHighlights(mol, [], { matched_atom_indices: '0~7:7~0' });
  expect(result.bonds).toEqual([0, 6]);
  expect(result.stickAtoms).toEqual(new Set([0, 1, 7]));
  expect(result.crossCarbons.size).toBe(0);
});

test('equivalent cross-peak protons include every endpoint bond', () => {
  const result = tocsyHighlights(mol, [], { matched_atom_indices: '0~7:6' });
  expect(new Set(result.bonds)).toEqual(new Set([0, 1, 2, 3, 4, 5, 6]));
});

test('ppm fallback matches paired shifts and supports decimal commas', () => {
  const shifts = [1.1234, 0, 0, 0, 0, 0, 2.5, 1.1234];
  const result = tocsyHighlights(mol, shifts, { ppm_molecule_matched: '1,123|2,5' });
  expect(new Set(result.bonds)).toEqual(new Set([0, 1, 2, 3, 4, 5, 6]));
  expect(tocsyHighlights(mol, shifts, { ppm_molecule_matched: '1.123|9' }).bonds).toEqual([]);
});

test('missing rows and invalid indices produce no highlights', () => {
  expect(tocsyHighlights(mol).bonds).toEqual([]);
  expect(tocsyHighlights(mol, [], { matched_atom_indices: '0:999||-1:6||1:6||0:bad' }).protons.size).toBe(0);
});

test('a molecule with a null TOCSY row remains renderable without TOCSY highlights', () => {
  expect(tocsyHighlights(mol, [], null)).toEqual({
    protons: new Set(), bonds: [], stickAtoms: new Set(), crossCarbons: new Set(),
  });
});

test('path search handles cycles and disconnected fragments', () => {
  const ring = molecule(['H', 'C', 'C', 'C', 'H', 'H'],
    [[0, 1], [1, 2], [2, 3], [3, 1], [3, 4]]);
  expect(tocsyHighlights(ring, [], { matched_atom_indices: '0:4' }).bonds).toEqual([0, 3, 4]);
  expect(tocsyHighlights(ring, [], { matched_atom_indices: '0:5' }).bonds).toEqual([]);
});
