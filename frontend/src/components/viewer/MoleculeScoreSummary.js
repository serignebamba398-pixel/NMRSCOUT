import React from 'react';

export function moleculeScoreSummary(rowsByType = {}, averageScore = null) {
  const number = value => {
    if (value == null || String(value).trim() === '') return null;
    const parsed = Number(String(value).trim().replace(/%$/, '').replace(',', '.'));
    return Number.isFinite(parsed) ? parsed : null;
  };
  const mean = values => {
    const available = values.map(number).filter(value => value !== null);
    return available.length ? available.reduce((sum, value) => sum + value, 0) / available.length : null;
  };
  const rows = Object.values(rowsByType || {}).filter(Boolean);
  const score = number(averageScore);
  return {
    score: score !== null && score >= 0 ? score : null,
    sub: mean(rows.map(row => row.coverage_sub)),
    variance: mean(rows.map(row => row.avg_variabilite)),
  };
}

export default function MoleculeScoreSummary({ source, rowsByType, averageScore }) {
  const displaySource = String(source ?? '').replace(/_predicted$/i, '');
  const summary = moleculeScoreSummary(rowsByType, averageScore);
  const format = (value, digits, suffix = '') => value === null ? 'N/A' : `${value.toFixed(digits)}${suffix}`;
  return (
    <div style={{ borderTop: '1px solid #e5e7eb', paddingTop: 8, fontSize: 11 }}>
      <div title={displaySource} style={{ fontWeight: 600, textAlign: 'center', overflowWrap: 'anywhere', marginBottom: 8 }}>
        {displaySource}
      </div>
      <dl style={{ display: 'grid', gridTemplateColumns: '1fr auto', gap: '4px 8px', margin: 0 }}>
        <dt title="Average score used to rank the molecules">Score</dt>
        <dd style={{ margin: 0, fontWeight: 700 }}>{format(summary.score, 1)}</dd>
        <dt title="Mean Sub score across spectra with an available value">Sub</dt>
        <dd style={{ margin: 0, fontWeight: 600 }}>{format(summary.sub, 2, '%')}</dd>
        <dt title="Mean variance across spectra with an available value">Var</dt>
        <dd style={{ margin: 0, fontWeight: 600 }}>{format(summary.variance, 4)}</dd>
      </dl>
    </div>
  );
}
