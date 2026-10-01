import React from 'react';
import { render, screen } from '@testing-library/react';
import MoleculeScoreSummary, { moleculeScoreSummary } from './MoleculeScoreSummary';

test('averages available values including zero and decimal commas', () => {
  expect(moleculeScoreSummary({
    '1H': { coverage_sub: '100', avg_variabilite: '0,0452' },
    '13C': { coverage_sub: 0, avg_variabilite: '0.0824' },
    COSY: { coverage_sub: '', avg_variabilite: null },
    HSQC: { coverage_sub: 'invalid', avg_variabilite: 'NaN' },
  }, 75)).toEqual({ score: 75, sub: 50, variance: 0.0638 });
});

test('does not turn missing scores into zero', () => {
  expect(moleculeScoreSummary({}, -1)).toEqual({ score: null, sub: null, variance: null });
  expect(moleculeScoreSummary({}, 0).score).toBe(0);
});

test('shows the source without the prediction suffix and formatted metrics', () => {
  render(<MoleculeScoreSummary source="1H_NPA020091_predicted" averageScore={100}
    rowsByType={{ '1H': { coverage_sub: 100, avg_variabilite: 0.0452 } }} />);
  ['1H_NPA020091', 'Score', '100.0', 'Sub', '100.00%', 'Var', '0.0452']
    .forEach(text => expect(screen.getByText(text)).toBeTruthy());
});
