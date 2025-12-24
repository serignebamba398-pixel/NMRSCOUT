const { logWithMemory } = require('../utils/logger');

const notFoundHandler = (req, res) => {
  res.status(404).json({
    error: 'Route not found',
    method: req.method,
    url: req.url,
    timestamp: new Date().toISOString(),
  });
};

const errorHandler = (err, req, res, next) => {
  logWithMemory(`Global error: ${err.message}`, 'ERROR');

  if (
    err &&
    /Invalid file|Invalid extension|File too large|Extension invalide|Fichier trop volumineux|upload/i.test(
      err.message
    )
  ) {
    return res.status(400).json({ error: err.message });
  }

  res.status(500).json({
    error: 'Internal server error',
    details: process.env.NODE_ENV === 'development' ? err.message : undefined,
    timestamp: new Date().toISOString(),
  });
};

module.exports = { notFoundHandler, errorHandler };