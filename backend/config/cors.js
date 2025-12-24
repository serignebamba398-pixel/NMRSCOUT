const cors = require('cors');

const allowedOrigins = (process.env.CORS_ORIGINS || '')
  .split(',')
  .map(s => s.trim())
  .filter(Boolean);

const corsOrigins = [
  'http://localhost:3000',
  'http://127.0.0.1:3000',
  'http://localhost:5000',
  'http://127.0.0.1:5000',
  'null', // file:// origins
  process.env.FRONTEND_URL,
  ...allowedOrigins,
].filter(Boolean);

// In Electron mode, use a very permissive CORS config
const isElectron = process.env.IS_ELECTRON === 'true';

const corsMiddleware = isElectron
  ? cors({
      origin: true, // Accept any origin
      credentials: true,
      methods: ['GET', 'POST', 'PUT', 'DELETE', 'OPTIONS', 'HEAD'],
      allowedHeaders: ['Content-Type', 'Authorization', 'Cache-Control', 'Accept', 'X-Requested-With'],
      exposedHeaders: ['Content-Length', 'Content-Type'],
    })
  : cors({
      origin: function (origin, cb) {
        if (!origin) return cb(null, true);
        if (corsOrigins.indexOf(origin) !== -1) return cb(null, true);
        if (origin === 'null') return cb(null, true); // file:// sends 'null'
        return cb(null, false);
      },
      credentials: true,
      methods: ['GET', 'POST', 'PUT', 'DELETE', 'OPTIONS'],
      allowedHeaders: ['Content-Type', 'Authorization', 'Cache-Control'],
    });

module.exports = { corsMiddleware, corsOrigins };