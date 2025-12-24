const multer = require('multer');
const path = require('path');
const { uploadsDir } = require('../config/paths');
const { MAX_UPLOAD_MB } = require('../config');
const { logWithMemory } = require('../utils/logger');

function sanitizeFilename(name) {
  const base = (name || '')
    .toString()
    .replace(/[\u0000-\u001f<>:"/\\|?*]+/g, '_')
    .replace(/\.{2,}/g, '.');
  return base.substring(0, 100);
}

const storage = multer.diskStorage({
  destination: (req, file, cb) => cb(null, uploadsDir),
  filename: (req, file, cb) => {
    const ts = Date.now();
    cb(null, `${ts}_${sanitizeFilename(file.originalname)}`);
  },
});

const upload = multer({
  storage,
  limits: {
    fileSize: MAX_UPLOAD_MB * 1024 * 1024,
    files: 20,
    fields: 500,
    parts: 1000,
  },
  fileFilter: (req, file, cb) => {
    const name = (file.originalname || '').toLowerCase();
    const mimeType = (file.mimetype || '').toLowerCase();

    const validCSV =
      name.endsWith('.csv') &&
      (mimeType === 'text/csv' ||
        mimeType === 'application/csv' ||
        mimeType === 'text/plain' ||
        mimeType === 'application/vnd.ms-excel');

    const validSDF =
      name.endsWith('.sdf') &&
      (mimeType === 'chemical/x-mdl-sdfile' ||
        mimeType === 'text/plain' ||
        mimeType === 'application/octet-stream');

    if (!validCSV && !validSDF) {
      logWithMemory(`Upload refused: ${name} (type: ${mimeType})`, 'WARN');
      return cb(new Error('Invalid file: unsupported extension or MIME type'));
    }

    cb(null, true);
  },
});

module.exports = upload;