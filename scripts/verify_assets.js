#!/usr/bin/env node
/*
  Verify packaging/runtime assets are present before building.

  Usage:
    node scripts/verify_assets.js           # warnings only
    node scripts/verify_assets.js --strict  # exit 1 if missing
*/

const fs = require('fs');
const path = require('path');

const strict = process.argv.includes('--strict');

function exists(p) {
  try {
    fs.accessSync(p, fs.constants.F_OK);
    return true;
  } catch {
    return false;
  }
}

function firstExisting(paths) {
  for (const p of paths) if (exists(p)) return p;
  return null;
}

const root = path.resolve(__dirname, '..');

const checks = [];

// Electron-builder extraResources expect these to exist when packaging.
checks.push({
  label: 'frontend build',
  requiredFor: 'packaging',
  paths: [path.join(root, 'frontend', 'build', 'index.html')],
});

checks.push({
  label: 'bundled python dist (PyInstaller output)',
  requiredFor: 'packaging',
  paths: [path.join(root, 'backend', 'dist', 'nmr_engine_production')],
});

checks.push({
  label: 'backend scripts',
  requiredFor: 'packaging',
  paths: [path.join(root, 'backend', 'scripts')],
});

checks.push({
  label: 'backend models',
  requiredFor: 'packaging',
  paths: [path.join(root, 'backend', 'models')],
});

// R-Portable is optional for some workflows but required if you rely on R scripts.
checks.push({
  label: 'R-Portable (optional but recommended)',
  requiredFor: 'runtime (if using R scripts)',
  paths: [
    path.join(root, 'R-Portable', 'bin', 'Rscript.exe'),
    path.join(root, 'R-Portable', 'bin', 'Rscript'),
  ],
  optional: true,
});

let hasErrors = false;

console.log('[verify-assets] repo:', root);
console.log('[verify-assets] mode:', strict ? 'strict' : 'warn');

for (const c of checks) {
  const found = firstExisting(c.paths);
  if (found) {
    console.log(`[OK] ${c.label}: ${path.relative(root, found)}`);
    continue;
  }

  const msg = `[MISSING] ${c.label} (${c.requiredFor})`;
  if (c.optional) {
    console.log(`${msg} — optional`);
    continue;
  }

  console.log(msg);
  hasErrors = true;
}

if (hasErrors) {
  console.log('\nNext steps:');
  console.log('  - Build frontend:   cd frontend && npm run build');
  console.log('  - Build python exe: cd backend && powershell -ExecutionPolicy Bypass -File build_bundle_complete.ps1');
  console.log('  - Then package:     npm run dist');

  if (strict) process.exit(1);
}
