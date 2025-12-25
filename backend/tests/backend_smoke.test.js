const test = require('node:test');
const assert = require('node:assert/strict');
const { spawn } = require('node:child_process');
const net = require('node:net');
const os = require('node:os');
const path = require('node:path');
const fs = require('node:fs');

async function getFreePort() {
  return await new Promise((resolve, reject) => {
    const server = net.createServer();
    server.unref();
    server.on('error', reject);
    server.listen(0, '127.0.0.1', () => {
      const { port } = server.address();
      server.close(() => resolve(port));
    });
  });
}

async function waitForHealth(port, timeoutMs = 15000) {
  const start = Date.now();
  // Node 18+ has global fetch
  while (Date.now() - start < timeoutMs) {
    try {
      const res = await fetch(`http://127.0.0.1:${port}/health`);
      if (res.ok) return await res.json();
    } catch {}
    await new Promise((r) => setTimeout(r, 250));
  }
  throw new Error('Timed out waiting for /health');
}

function startBackend({ port }) {
  const logDir = fs.mkdtempSync(path.join(os.tmpdir(), 'annotation-nmr-backend-logs-'));
  const env = {
    ...process.env,
    PORT: String(port),
    NODE_ENV: 'production',
    LOG_DIR: logDir,
    // Ensure health endpoint doesn't depend on Electron env during tests
    IS_ELECTRON: 'false',
    IS_ELECTRON_BUNDLED: 'false',
  };

  const child = spawn(process.execPath, ['server.js'], {
    cwd: __dirname + '/..',
    env,
    stdio: ['ignore', 'pipe', 'pipe'],
  });

  let stdout = '';
  let stderr = '';
  child.stdout.on('data', (d) => (stdout += d.toString()));
  child.stderr.on('data', (d) => (stderr += d.toString()));

  child.getLogs = () => ({ stdout, stderr });
  return child;
}

async function stopProcess(child) {
  if (!child || child.killed) return;
  child.kill('SIGTERM');
  await new Promise((r) => setTimeout(r, 250));
  if (!child.killed) child.kill('SIGKILL');
}

test('backend /health responds OK', async () => {
  const port = await getFreePort();
  const child = startBackend({ port });
  try {
    const body = await waitForHealth(port);
    assert.equal(body.status, 'OK');
    assert.ok(body.timestamp);
    assert.ok(body.node);
  } catch (e) {
    const logs = child.getLogs?.() || {};
    e.message += `\n--- backend stdout ---\n${logs.stdout || ''}\n--- backend stderr ---\n${logs.stderr || ''}`;
    throw e;
  } finally {
    await stopProcess(child);
  }
});

test('backend /api/path/normalize normalizes in non-Electron mode', async () => {
  const port = await getFreePort();
  const child = startBackend({ port });
  try {
    await waitForHealth(port);

    const res = await fetch(`http://127.0.0.1:${port}/api/path/normalize`, {
      method: 'POST',
      headers: { 'content-type': 'application/json' },
      body: JSON.stringify({ rawPath: 'C:\\Users\\test\\outputs' }),
    });
    assert.ok(res.ok);
    const body = await res.json();
    assert.equal(body.ok, true);
    assert.equal(body.path, '/host/c/Users/test/outputs');
  } catch (e) {
    const logs = child.getLogs?.() || {};
    e.message += `\n--- backend stdout ---\n${logs.stdout || ''}\n--- backend stderr ---\n${logs.stderr || ''}`;
    throw e;
  } finally {
    await stopProcess(child);
  }
});
