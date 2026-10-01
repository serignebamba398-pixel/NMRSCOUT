import { useState } from 'react';
import { act, renderHook } from '@testing-library/react';
import { useSpectrumWindowSettings } from './useSpectrumWindowSettings';

const bounds = { f1_ppm_min: 1, f1_ppm_max: 4, f2_ppm_min: 2, f2_ppm_max: 7 };
function useWindow(options) {
  const [zoom, setZoom] = useState({});
  useSpectrumWindowSettings(options, setZoom);
  return { zoom, setZoom };
}

test('applies all four TOCSY limits without Bruker data and updates after edits', () => {
  const { result, rerender } = renderHook(useWindow, { initialProps: { TOCSY: bounds } });
  expect(result.current.zoom.TOCSY).toEqual({ f1: { min: 1, max: 4 }, f2: { min: 2, max: 7 } });
  rerender({ TOCSY: { ...bounds, f1_ppm_min: '0', f2_ppm_max: '8' } });
  expect(result.current.zoom.TOCSY).toEqual({ f1: { min: 0, max: 4 }, f2: { min: 2, max: 8 } });
});

test('ordinary rerenders and COSY changes preserve a manual TOCSY zoom', () => {
  const { result, rerender } = renderHook(useWindow, { initialProps: { TOCSY: bounds } });
  const manual = { f1: { min: 2, max: 3 }, f2: { min: 3, max: 4 } };
  act(() => result.current.setZoom({ TOCSY: manual }));
  rerender({ TOCSY: { ...bounds }, COSY: bounds });
  expect(result.current.zoom.TOCSY).toEqual(manual);
  expect(result.current.zoom.COSY).toEqual({ f1: { min: 1, max: 4 }, f2: { min: 2, max: 7 } });
});

test('normalizes reversed bounds and ignores incomplete or zero-width windows', () => {
  const { result, rerender } = renderHook(useWindow, { initialProps: { TOCSY: { ...bounds, f1_ppm_min: 5 } } });
  const valid = result.current.zoom.TOCSY;
  expect(valid.f1).toEqual({ min: 4, max: 5 });
  for (const value of ['', null, NaN, 4]) {
    rerender({ TOCSY: { ...bounds, f1_ppm_min: value } });
    expect(result.current.zoom.TOCSY).toEqual(valid);
  }
});
