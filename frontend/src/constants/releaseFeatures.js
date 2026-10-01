// Enable this when TOCSY is included in a published version.
export const TOCSY_ENABLED = false;

export const isMethodAvailable = (method) =>
  TOCSY_ENABLED || String(method).toUpperCase() !== 'TOCSY';

export const availableMethodEntries = (values = {}) =>
  Object.fromEntries(Object.entries(values).filter(([method]) => isMethodAvailable(method)));

export const isArtifactAvailable = (artifact) =>
  TOCSY_ENABLED || !/(^|[_\W])TOCSY([_\W]|$)/i.test(
    [artifact.type, artifact.file, artifact.url].filter(Boolean).join(' ')
  );

export const restoreAvailableSelections = (state) => ({
  ...state,
  ...(state.selectedPredMethods && {
    selectedPredMethods: state.selectedPredMethods.filter(isMethodAvailable),
  }),
  ...(state.selectedMethods && {
    selectedMethods: state.selectedMethods.filter(isMethodAvailable),
  }),
  ...(state.activeType && !isMethodAvailable(state.activeType) && { activeType: '1H' }),
  ...(state.toleranceOptimSelectedType && !isMethodAvailable(state.toleranceOptimSelectedType) && {
    toleranceOptimSelectedType: '1H',
  }),
});
