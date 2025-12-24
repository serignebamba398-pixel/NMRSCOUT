/**
 * Formatting utilities for display
 */

/**
 * Formats a file size in bytes to a readable unit
 * @param {number} bytes - Size in bytes
 * @returns {string} Formatted size (ex: "1.5 MB")
 */
export const formatFileSize = (bytes) => {
  if (!Number.isFinite(bytes) || bytes < 0) return "—";
  
  const units = ["B", "KB", "MB", "GB", "TB"];
  let size = bytes;
  let unitIndex = 0;

  while (size >= 1024 && unitIndex < units.length - 1) {
    size /= 1024;
    unitIndex++;
  }

  return `${size.toFixed(unitIndex === 0 ? 0 : 1)} ${units[unitIndex]}`;
};

/**
 * Formats a duration in milliseconds to a readable string
 * @param {number} ms - Duration in milliseconds
 * @returns {string} Formatted duration (ex: "2m 30s")
 */
export const formatDuration = (ms) => {
  if (!Number.isFinite(ms) || ms < 0) return "—";

  const seconds = Math.floor(ms / 1000);
  const minutes = Math.floor(seconds / 60);
  const hours = Math.floor(minutes / 60);

  if (hours > 0) {
    return `${hours}h ${minutes % 60}m ${seconds % 60}s`;
  }
  if (minutes > 0) {
    return `${minutes}m ${seconds % 60}s`;
  }
  return `${seconds}s`;
};

/**
 * Shortens a CSV label for display
 * @param {string} filename - Full file name
 * @returns {string} Shortened label
 */
export const shortCsvLabel = (filename) => {
  if (!filename || typeof filename !== "string") return "—";
  
  const parts = filename.split("_");
  if (parts.length >= 3) {
    return `${parts[0]}_${parts[1]}_...`;
  }
  
  // If name is long, truncate
  if (filename.length > 40) {
    return filename.substring(0, 37) + "...";
  }
  
  return filename;
};

/**
 * Formats a number with thousands separators
 * @param {number} num - Number to format
 * @returns {string} Formatted number (ex: "1,234,567")
 */
export const formatNumber = (num) => {
  if (!Number.isFinite(num)) return "—";
  return num.toLocaleString("en-US");
};

/**
 * Formats a percentage with precision
 * @param {number} value - Value between 0 and 100
 * @param {number} decimals - Number of decimal places (default: 2)
 * @returns {string} Formatted percentage (ex: "45.67%")
 */
export const formatPercent = (value, decimals = 2) => {
  if (!Number.isFinite(value)) return "—";
  return `${value.toFixed(decimals)}%`;
};