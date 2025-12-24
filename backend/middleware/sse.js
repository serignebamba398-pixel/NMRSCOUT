const { corsOrigins } = require('../config/cors');
const { getSystemMemoryUsage, logWithMemory } = require('../utils/logger');

function sseInit(res, origin) {
  const allowOrigin =
    origin && corsOrigins.includes(origin) ? origin : corsOrigins[0] || origin || '*';

  const headers = {
    'Content-Type': 'text/event-stream',
    'Cache-Control': 'no-cache, no-transform, no-store',
    Connection: 'keep-alive',
    'X-Accel-Buffering': 'no',
    'Keep-Alive': 'timeout=0, max=0',
    'Transfer-Encoding': 'chunked',
    Vary: 'Origin',
  };

  if (allowOrigin && allowOrigin !== '*') {
    headers['Access-Control-Allow-Origin'] = allowOrigin;
    headers['Access-Control-Allow-Credentials'] = 'true';
    headers['Access-Control-Expose-Headers'] = 'Content-Type';
    headers['Access-Control-Allow-Headers'] = 'Content-Type, Authorization, Cache-Control';
  } else {
    headers['Access-Control-Allow-Origin'] = '*';
  }

  res.writeHead(200, headers);
  res.setTimeout?.(0);
  res.flushHeaders?.();

  const hb = setInterval(() => {
    try {
      res.write(`: keepalive ${Date.now()}\n\n`);
    } catch {
      clearInterval(hb);
    }
  }, 10000);

  return hb;
}

function sseSend(res, data) {
  try {
    const json = JSON.stringify({
      ...data,
      timestamp: new Date().toISOString(),
      serverMemory: getSystemMemoryUsage(),
    });
    res.write(`data: ${json}\n\n`);
  } catch (e) {
    logWithMemory(`SSE send error: ${e.message}`, 'ERROR');
  }
}

module.exports = { sseInit, sseSend };