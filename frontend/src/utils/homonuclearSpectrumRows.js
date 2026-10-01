// Adapt atom-pair prediction CSVs to the existing annotation viewer shape.
// Already annotated rows keep their existing representation.
export function homonuclearSpectrumRows(rows, type) {
  if (!['COSY', 'TOCSY'].includes(type) || !rows.some(row => row.Filename && 'Delta1' in row && 'Delta2' in row)) return rows;
  const groups = new Map();
  for (const row of rows) {
    const source = row.Source || row.Filename;
    if (!source) continue;
    if (!groups.has(source)) groups.set(source, { ...row, Source: source, predictionPeaks: [] });
    const d1 = parseFloat(row.Delta1), d2 = parseFloat(row.Delta2);
    if (Number.isFinite(d1) && Number.isFinite(d2)) groups.get(source).predictionPeaks.push(row);
  }
  return [...groups.values()].map(row => ({
    ...row,
    ppm_molecule: row.predictionPeaks.map(peak => `${peak.Delta1}|${peak.Delta2}`).join('||'),
  }));
}
