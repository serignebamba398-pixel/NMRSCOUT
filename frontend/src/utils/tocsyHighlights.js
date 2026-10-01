import { findBondsInPath, parseAtomSymbolsFromMolblock } from './molParser';

// Use the same molecular paths for the RDKit bonds and the 3D stick endpoints.
export function tocsyHighlights(molblock, shifts = [], row = {}) {
  // Missing experiment results are explicitly null in the combined viewer.
  if (row == null) {
    return { protons: new Set(), bonds: [], stickAtoms: new Set(), crossCarbons: new Set() };
  }
  const atoms = parseAtomSymbolsFromMolblock(molblock);
  const protons = new Set();
  const crossCarbons = new Set();
  const stickAtoms = new Set();
  const pairs = [];
  const crossPairs = [];
  const addGroups = (left, right) => {
    if (!left.length || !right.length) return;
    const diagonal = left.length === right.length && left.every(h => right.includes(h));
    [...left, ...right].forEach(h => protons.add(h));
    if (diagonal) {
      left.forEach(h => pairs.push([h, h]));
    } else {
      for (const h1 of left) for (const h2 of right) {
        pairs.push([h1, h2]);
        crossPairs.push([h1, h2]);
      }
    }
  };
  if (row.matched_atom_indices) {
    const group = value => [...new Set(String(value || '').split('~')
      .filter(s => /^\d+$/.test(s.trim())).map(Number)
      .filter(i => atoms[i] === 'H'))];
    for (const entry of row.matched_atom_indices.split('||')) {
      const parts = entry.split(':');
      if (parts.length === 2) addGroups(group(parts[0]), group(parts[1]));
    }
  } else {
    const round = v => Math.round(v * 1000) / 1000;
    const group = ppm => atoms.map((symbol, i) =>
      symbol === 'H' && Number.isFinite(shifts[i]) && round(shifts[i]) === round(ppm) ? i : -1
    ).filter(i => i >= 0);
    for (const entry of (row.ppm_molecule_matched || '').split('||')) {
      const parts = entry.split('|');
      if (parts.length !== 2 || parts.some(p => !p.trim())) continue;
      const values = parts.map(p => Number(p.replace(',', '.')));
      if (values.every(Number.isFinite)) addGroups(group(values[0]), group(values[1]));
    }
  }
  // TOCSY can span more than the five bonds used by the COSY display.
  const bonds = findBondsInPath(molblock, pairs, Infinity);
  const crossBonds = new Set(findBondsInPath(molblock, crossPairs, Infinity));
  const lines = (molblock || '').split(/\r?\n/);
  for (const bond of bonds) {
    const line = lines[4 + atoms.length + bond] || '';
    const endpoints = [Number(line.slice(0, 3)) - 1, Number(line.slice(3, 6)) - 1];
    for (const index of endpoints) {
      stickAtoms.add(index);
      if (crossBonds.has(bond) && atoms[index] === 'C') crossCarbons.add(index);
    }
  }
  return { protons, bonds, stickAtoms, crossCarbons };
}
