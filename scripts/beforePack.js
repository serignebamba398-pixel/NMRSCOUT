/**
 * beforePack hook for electron-builder
 * Recursively remove files with reserved Windows names
 * (nul, con, prn, aux, com1-com9, lpt1-lpt9) before packaging.
 * These names cause electron-builder to fail on Windows.
 */
const fs = require('fs');
const path = require('path');

// Reserved Windows basenames, matched case-insensitively.
const WINDOWS_RESERVED = new Set([
  'nul', 'con', 'prn', 'aux',
  'com1','com2','com3','com4','com5','com6','com7','com8','com9',
  'lpt1','lpt2','lpt3','lpt4','lpt5','lpt6','lpt7','lpt8','lpt9',
]);

function isReservedName(filename) {
  // Check the basename without its extension; "nul.txt" is also reserved.
  const base = path.basename(filename, path.extname(filename)).toLowerCase();
  return WINDOWS_RESERVED.has(base);
}

function deleteReservedFiles(dir) {
  let count = 0;
  let entries;
  try {
    entries = fs.readdirSync(dir);
  } catch {
    return count;
  }

  for (const entry of entries) {
    const fullPath = path.join(dir, entry);
    let stat;
    try {
      stat = fs.statSync(fullPath);
    } catch {
      continue;
    }

    if (stat.isDirectory()) {
      count += deleteReservedFiles(fullPath);
    } else if (isReservedName(entry)) {
      try {
        // Use the \\?\ prefix to access reserved Windows paths.
        const uncPath = '\\\\?\\' + path.resolve(fullPath);
        fs.unlinkSync(uncPath);
        console.log(`[beforePack] Supprimé fichier réservé : ${fullPath}`);
        count++;
      } catch (err) {
        // Fallback : cmd del
        try {
          const { execSync } = require('child_process');
          execSync(`del /f /q "\\\\?\\${path.resolve(fullPath)}"`, { stdio: 'pipe' });
          console.log(`[beforePack] Supprimé via cmd : ${fullPath}`);
          count++;
        } catch {
          console.warn(`[beforePack] Impossible de supprimer : ${fullPath}`);
        }
      }
    }
  }
  return count;
}

exports.default = async function(context) {
  const projectDir = context.packager.projectDir;
  console.log('[beforePack] Recherche de fichiers réservés Windows...');

  // Directories to scan.
  const dirsToScan = [
    path.join(projectDir, 'backend'),
    path.join(projectDir, 'frontend', 'build'),
  ];

  let total = 0;
  for (const dir of dirsToScan) {
    if (fs.existsSync(dir)) {
      total += deleteReservedFiles(dir);
    }
  }

  if (total > 0) {
    console.log(`[beforePack] ${total} fichier(s) réservé(s) supprimé(s)`);
  } else {
    console.log('[beforePack] Aucun fichier réservé trouvé');
  }
};
