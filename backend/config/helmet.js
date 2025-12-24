const helmet = require('helmet');

const isElectron = process.env.IS_ELECTRON === 'true';

// In Electron mode, disable security restrictions
// because we run locally and file:// needs access to localhost
const helmetConfig = isElectron
  ? (req, res, next) => next() // No-op middleware
  : helmet({
      crossOriginEmbedderPolicy: false,
      contentSecurityPolicy: {
        directives: {
          defaultSrc: ["'self'"],
          scriptSrc: ["'self'", "'unsafe-inline'", "'unsafe-eval'"],
          styleSrc: ["'self'", "'unsafe-inline'"],
          connectSrc: [
            "'self'",
            'http://localhost:4000',
            'http://localhost:3000',
            'http://localhost:5000',
            'http://127.0.0.1:4000',
            'http://127.0.0.1:3000',
            'http://127.0.0.1:5000',
          ],
          imgSrc: ["'self'", 'data:', 'blob:'],
          fontSrc: ["'self'", 'data:'],
          objectSrc: ["'none'"],
          mediaSrc: ["'self'"],
          frameSrc: ["'none'"],
        },
      },
    });

module.exports = helmetConfig;