/**
 * Parser CSV intelligent avec détection automatique du délimiteur
 */

/**
 * Parse un texte CSV en objet { headers, data }
 * @param {string} text - Contenu du CSV
 * @returns {{ headers: string[], data: object[] }}
 */
export function parseCSVSmart(text) {
  if (!text || typeof text !== "string") {
    return { headers: [], data: [] };
  }

  // Détection du délimiteur en analysant les 5 premières lignes
  const firstLines = text.split(/\r?\n/).slice(0, 5).join("\n");
  const candidates = [";", ",", "\t"];
  let bestDelimiter = ";";
  let maxCount = -1;

  for (const delimiter of candidates) {
    const counts = firstLines
      .split(/\r?\n/)
      .map((line) => (line.match(new RegExp(`\\${delimiter}`, "g")) || []).length);
    
    const count = Math.max(...counts, 0);
    if (count > maxCount) {
      bestDelimiter = delimiter;
      maxCount = count;
    }
  }

  // Parsing avec gestion des guillemets doubles
  const rows = [];
  let currentRow = [];
  let currentValue = "";
  let inQuotes = false;

  const pushValue = () => {
    // Retirer BOM si présent au début
    if (currentValue && currentValue.charCodeAt(0) === 0xfeff) {
      currentValue = currentValue.slice(1);
    }
    currentRow.push(currentValue);
    currentValue = "";
  };

  const pushRow = () => {
    if (currentRow.length > 0) {
      rows.push(currentRow);
      currentRow = [];
    }
  };

  for (let i = 0; i < text.length; i++) {
    const char = text[i];

    if (inQuotes) {
      if (char === '"') {
        // Échappement de guillemet : "" → "
        if (text[i + 1] === '"') {
          currentValue += '"';
          i++;
        } else {
          inQuotes = false;
        }
      } else {
        currentValue += char;
      }
      continue;
    }

    if (char === '"') {
      inQuotes = true;
    } else if (char === bestDelimiter) {
      pushValue();
    } else if (char === "\n") {
      pushValue();
      pushRow();
    } else if (char !== "\r") {
      currentValue += char;
    }
  }

  // Dernière valeur/ligne
  if (currentValue.length > 0 || currentRow.length > 0) {
    pushValue();
    pushRow();
  }

  if (rows.length === 0) {
    return { headers: [], data: [] };
  }

  // Première ligne = headers
  const headers = rows[0].map((h) => (h || "").trim());

  // Lignes suivantes = data (filtrer lignes vides)
  const data = rows
    .slice(1)
    .filter((row) => row.some((cell) => (cell || "").trim() !== ""))
    .map((row) => {
      const obj = {};
      headers.forEach((header, index) => {
        obj[header] = (row[index] ?? "").trim();
      });
      return obj;
    });

  return { headers, data };
}

/**
 * Détermine si une ligne CSV est un header (première ligne)
 * @param {string} line - Ligne CSV
 * @returns {boolean}
 */
export const isHeaderLine = (line) => {
  if (!line) return false;
  
  const firstCell = line.split(/[;, \t]/)[0].trim().toLowerCase();
  if (!firstCell) return false;

  // Détection de headers communs
  const commonHeaders = ["ppm", "ppm_f1", "ppm_f2", "id", "smiles", "mol_id"];
  if (commonHeaders.includes(firstCell)) {
    return true;
  }

  // Si la première cellule n'est pas numérique, probablement un header
  const numericPattern = /^[-+]?\d*\.?\d+(e[-+]?\d+)?$/i;
  return !numericPattern.test(firstCell);
};