import { COLOR_OPTIONS, normalizeColorLabel } from './colors';

test('palette labels are English and saved French labels normalize by color value', () => {
  expect(COLOR_OPTIONS.map(({ name }) => name)).toContain('Pure Green');
  expect(normalizeColorLabel({ name: 'Vert Pur', hex: '#00FF00' }).name).toBe('Pure Green');
  expect(normalizeColorLabel({ name: 'Vert Émeraude', hex: '#00A878' }).name).toBe('Emerald Green');
});