export function tocsyPredictionOptions(options = {}) {
  const maxBonds = Number(options.max_bonds ?? 3);
  const minBonds = Number(options.min_bonds ?? 2);
  if (!Number.isInteger(maxBonds) || maxBonds < 1) {
    throw new Error('TOCSY: Maximum bonds must be an integer >= 1.');
  }
  if (!Number.isInteger(minBonds) || minBonds < 1 || minBonds > maxBonds) {
    throw new Error('TOCSY: Minimum bonds must be an integer >= 1 and <= Maximum bonds.');
  }
  const result = { min_bonds: minBonds, max_bonds: maxBonds, include_diagonal: options.include_diagonal ?? true };
  return result;
}
