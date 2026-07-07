const helmet = require('helmet');

const isElectron = process.env.IS_ELECTRON === 'true';
const isDev = process.env.NODE_ENV !== 'production';

// In Electron mode we still want security headers; we only relax CSP enough
// to allow local development (React dev server) and local API calls.
const connectSrc = [
  "'self'",
  'http://localhost:3000',
  'http://localhost:4000',
  'http://localhost:5000',
  'http://127.0.0.1:3000',
  'http://127.0.0.1:4000',
  'http://127.0.0.1:5000',
  // Add ws:// for WebSocket/SSE connections in Electron
  'ws://localhost:3000',
  'ws://localhost:5000',
  'ws://127.0.0.1:3000',
  'ws://127.0.0.1:5000',
];

const scriptSrc = ["'self'", "'unsafe-inline'"];
if (isDev) scriptSrc.push("'unsafe-eval'");

// Disable CSP entirely in Electron mode to avoid fetch blocking issues
const helmetConfig = isElectron
  ? helmet({
      crossOriginEmbedderPolicy: false,
      contentSecurityPolicy: false, // Disable CSP in Electron
      crossOriginResourcePolicy: false,
    })
  : helmet({
      crossOriginEmbedderPolicy: false,
      contentSecurityPolicy: {
        directives: {
          defaultSrc: ["'self'"],
          // CRA dev tooling may require eval in dev; keep it off in production.
          scriptSrc,
          styleSrc: ["'self'", "'unsafe-inline'"],
          connectSrc,
          imgSrc: ["'self'", 'data:', 'blob:'],
          fontSrc: ["'self'", 'data:'],
          objectSrc: ["'none'"],
          mediaSrc: ["'self'"],
          frameSrc: ["'none'"],
          baseUri: ["'self'"],
          formAction: ["'self'"],
          frameAncestors: ["'none'"],
        },
      },
      crossOriginResourcePolicy: undefined,
    });

module.exports = helmetConfig;