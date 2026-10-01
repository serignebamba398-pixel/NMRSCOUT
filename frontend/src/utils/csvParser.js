/**
 * CSV parser with delimiter detection.
 */

/**
 * Parse CSV text into { headers, data }.
 * @param {string} text - CSV text.
 * @returns {{ headers: string[], data: object[] }}
 */
export function parseCSVSmart(text) {
  if (!text || typeof text !== "string") {
    return { headers: [], data: [] };
  }

  // Detect the delimiter from the first five lines.
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

  // Parse quoted fields and escaped double quotes.
  const rows = [];
  let currentRow = [];
  let currentValue = "";
  let inQuotes = false;

  const pushValue = () => {
    // Remove a leading byte-order mark, if present.
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
        // Decode escaped double quotes: "" → ".
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

  // Flush the last field and row.
  if (currentValue.length > 0 || currentRow.length > 0) {
    pushValue();
    pushRow();
  }

  if (rows.length === 0) {
    return { headers: [], data: [] };
  }

  // Use the first row as column headers.
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
 * Check whether a CSV line looks like a header.
 * @param {string} line - CSV line.
 * @returns {boolean}
 */
export const isHeaderLine = (line) => {
  if (!line) return false;
  
  const firstCell = line.split(/[;, \t]/)[0].trim().toLowerCase();
  if (!firstCell) return false;

  // Recognize common column names.
  const commonHeaders = ["ppm", "ppm_f1", "ppm_f2", "id", "smiles", "mol_id"];
  if (commonHeaders.includes(firstCell)) {
    return true;
  }

  // Treat a nonnumeric first cell as a likely header.
  const numericPattern = /^[-+]?\d*\.?\d+(e[-+]?\d+)?$/i;
  return !numericPattern.test(firstCell);
};