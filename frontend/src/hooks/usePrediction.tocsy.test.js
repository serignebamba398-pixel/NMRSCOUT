import React from 'react';
import { act, render, renderHook, fireEvent, screen } from '@testing-library/react';
import { usePrediction } from './usePrediction';
import { runPredictionFromSmiles, runPredictionFromCSV } from '../api';
import { TocsyPredictionFields } from '../components/prediction/TocsyPredictionFields';
import { homonuclearSpectrumRows } from '../utils/homonuclearSpectrumRows';
import { parseCSVSmart } from '../utils/csvParser';
import { tocsyPredictionOptions } from '../utils/tocsyPrediction';

jest.mock('../api', () => ({
  runPredictionFromSmiles: jest.fn().mockResolvedValue({}),
  runPredictionFromCSV: jest.fn().mockResolvedValue({}),
  API_BASE_URL: '',
}));

beforeEach(() => jest.clearAllMocks());

for (const mode of ['SMILES', 'CSV']) {
  test(`${mode} excludes unavailable TOCSY selections from the request`, async () => {
    let state = {
      mode, smiles: 'CCCC', molId: 'mol', csvUploadedName: 'molecules.csv',
      selectedPredMethods: ['COSY', 'TOCSY'],
      options: { cosy: { min_bonds: 2, max_bonds: 3 },
        tocsy: { min_bonds: '4', max_bonds: '6', mixing_time: '80', include_diagonal: false } },
    };
    const { result } = renderHook(() => usePrediction({ predictionState: state,
      setPredictionState: fn => { state = fn(state); }, addLog: jest.fn(), setShowLogs: jest.fn() }));
    await act(async () => { await result.current.runPrediction(); });
    const request = mode === 'SMILES' ? runPredictionFromSmiles : runPredictionFromCSV;
    const args = request.mock.calls[0];
    expect(args[1].methods).toEqual(['COSY']);
    expect(args[1].global.extractors.tocsy).toBeUndefined();
    expect(args[1].global.extractors.cosy.max_bonds).toBe(3);
    act(() => args[5]({ success: true, taskId: 'task', outputs: { csv_TOCSY: '/outputs/task/prediction_TOCSY_task.csv' } }));
    expect(state.artifacts).toEqual([expect.objectContaining({ type: 'CSV_TOCSY' })]);
  });
}

test('TOCSY exposes bounds and diagonal without a mixing time field', () => {
  const change = jest.fn();
  render(<TocsyPredictionFields options={{ min_bonds: 2, max_bonds: 6 }} onChange={change} />);
  expect(screen.queryByLabelText('Mixing time (ms)')).toBeNull();
  fireEvent.change(screen.getByLabelText('Maximum bonds'), { target: { value: '8' } });
  expect(change).toHaveBeenCalledWith('max_bonds', '8');
});

test('frontend ignores saved mixing time and rejects invalid bounds', () => {
  expect(tocsyPredictionOptions({ max_bonds: 20, mixing_time: '' })).toEqual({ min_bonds: 2, max_bonds: 20, include_diagonal: true });
  for (const options of [{ max_bonds: 0 }, { min_bonds: 7, max_bonds: 6 }]) {
    expect(() => tocsyPredictionOptions(options)).toThrow('TOCSY:');
  }
});

test('prediction CSV points and metadata survive adaptation to the shared viewer', () => {
  const { data } = parseCSVSmart('Filename,H1,Delta1,H2,Delta2,Bonds,SpinSystem,RelaySteps,MaxBonds\na.sdf,4,1.2,5,2.3,3,0,1,6\na.sdf,5,2.3,4,1.2,3,0,1,6');
  for (const type of ['COSY', 'TOCSY']) {
    const rows = homonuclearSpectrumRows(data, type);
    expect(rows).toHaveLength(1);
    expect(rows[0].Source).toBe('a.sdf');
    expect(rows[0].ppm_molecule).toBe('1.2|2.3||2.3|1.2');
    expect(rows[0].predictionPeaks[0].MaxBonds).toBe('6');
  }
  const annotation = [{ Source: 'a', ppm_molecule: '1|2' }];
  expect(homonuclearSpectrumRows(annotation, 'TOCSY')).toBe(annotation);
});
