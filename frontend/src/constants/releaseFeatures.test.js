import {
  availableMethodEntries,
  isArtifactAvailable,
  isMethodAvailable,
  restoreAvailableSelections,
} from './releaseFeatures';

test('keeps TOCSY unavailable in the current release', () => {
  expect(isMethodAvailable('TOCSY')).toBe(false);
  expect(isMethodAvailable('COSY')).toBe(true);
  expect(availableMethodEntries({ COSY: '/cosy.csv', TOCSY: '/tocsy.csv' })).toEqual({
    COSY: '/cosy.csv',
  });
});

test('removes TOCSY from selections restored from older projects', () => {
  expect(restoreAvailableSelections({
    selectedPredMethods: ['1H', 'TOCSY'],
    selectedMethods: ['COSY', 'TOCSY'],
    activeType: 'TOCSY',
    toleranceOptimSelectedType: 'TOCSY',
  })).toEqual({
    selectedPredMethods: ['1H'],
    selectedMethods: ['COSY'],
    activeType: '1H',
    toleranceOptimSelectedType: '1H',
  });
});

test('hides TOCSY artifacts from results', () => {
  expect(isArtifactAvailable({ type: 'CSV_TOCSY', file: 'prediction_TOCSY_task.csv' })).toBe(false);
  expect(isArtifactAvailable({ type: 'CSV_COSY', file: 'prediction_COSY_task.csv' })).toBe(true);
});
