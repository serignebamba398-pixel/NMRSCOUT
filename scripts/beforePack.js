/**
 * beforePack hook for electron-builder
 * Supprime récursivement tous les fichiers avec des noms réservés Windows
 * (nul, con, prn, aux, com1-com9, lpt1-lpt9) avant le packaging.
 * Ces fichiers font crasher electron-builder sur Windows.
 */
const fs = require('fs');
const path = require('path');

// Noms réservés Windows (insensible à la casse, sans extension)
const WINDOWS_RESERVED = new Set([
  'nul', 'con', 'prn', 'aux',
  'com1','com2','com3','com4','com5','com6','com7','com8','com9',
  'lpt1','lpt2','lpt3','lpt4','lpt5','lpt6','lpt7','lpt8','lpt9',
]);

function isReservedName(filename) {
  // Vérifier le nom sans extension (ex: "nul.txt" est aussi réservé)
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
        // Utiliser le préfixe \\?\ pour bypasser les restrictions Windows
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

  // Dossiers à scanner
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
