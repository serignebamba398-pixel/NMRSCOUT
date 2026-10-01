// Replace one packed JS file while preserving all other archive contents.
const fs = require('fs');
const path = require('path');
const crypto = require('crypto');
const assert = require('assert');
const asar = require('@electron/asar');
const { Pickle } = require('../node_modules/@electron/asar/lib/pickle');

module.exports = function updatePackagedJs(archive, filename, updated, label) {
  const original = fs.readFileSync(archive);
  const raw = asar.getRawHeader(archive);
  const hash = value => crypto.createHash('sha256').update(value).digest('hex');
  const content = Buffer.from(updated);
  const payload = original.subarray(8 + raw.headerSize);
  let entry = raw.header;
  for (const name of filename.split(/[\\/]/)) entry = entry.files[name];
  assert(!entry.unpacked);
  entry.offset = String(payload.length);
  entry.size = content.length;
  const blockSize = 4 * 1024 * 1024;
  const blocks = [];
  for (let i = 0; i < content.length; i += blockSize) blocks.push(hash(content.subarray(i, i + blockSize)));
  entry.integrity = { algorithm: 'SHA256', hash: hash(content), blockSize, blocks };
  const header = Pickle.createEmpty();
  header.writeString(JSON.stringify(raw.header));
  const size = Pickle.createEmpty();
  size.writeUInt32(header.toBuffer().length);
  const output = Buffer.concat([size.toBuffer(), header.toBuffer(), payload, content]);
  const backupDir = path.resolve(__dirname, '..', 'backups', `${label}_${Date.now()}`);
  fs.mkdirSync(backupDir, { recursive: true });
  fs.writeFileSync(path.join(backupDir, 'app.asar'), original);
  const staged = archive + '.ui.tmp';
  fs.writeFileSync(staged, output);
  asar.uncacheAll();
  assert.equal(asar.extractFile(staged, filename).toString(), updated);
  let verified = 0;
  function check(node, prefix = '') {
    for (const [name, info] of Object.entries(node.files || {})) {
      const file = path.join(prefix, name);
      if (info.files) check(info, file);
      else if (!info.unpacked && !info.link && file !== filename) {
        assert(asar.extractFile(archive, file).equals(asar.extractFile(staged, file)), file);
        verified++;
      }
    }
  }
  check(raw.header);
  const report = { archive, filename, backupDir, verifiedUnchangedPackedFiles: verified,
    originalSha256: hash(original), sha256: hash(output), status: 'ready' };
  const reportPath = path.join(backupDir, 'report.json');
  fs.writeFileSync(reportPath, JSON.stringify(report, null, 2));
  assert.equal(hash(fs.readFileSync(archive)), hash(original));
  try {
    fs.renameSync(staged, archive);
  } catch (error) {
    if (error.code !== 'EPERM' && error.code !== 'EBUSY') throw error;
    console.log(JSON.stringify({ ...report, status: 'waiting_for_app_close', staged }, null, 2));
    process.exitCode = 2;
    return;
  }
  asar.uncacheAll();
  assert.equal(hash(asar.extractFile(archive, filename)), hash(content));
  report.status = 'applied';
  fs.writeFileSync(reportPath, JSON.stringify(report, null, 2));
  console.log(JSON.stringify(report, null, 2));
};
