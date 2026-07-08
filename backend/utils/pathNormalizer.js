function toForward(p) {
  return String(p || '')
    .trim()
    .replace(/\\/g, '/');
}

function lstripSlash(p) {
  return String(p).replace(/^\/+/, '');
}

/**
 * Normalize a host filesystem path to its container mount path.
 *
 * The install scripts (install.ps1 / install.sh) generate docker-compose.volumes.yml
 * which mounts host directories into the container under /host/:
 *
 *   Windows (all drives):
 *     C:/  -> /host/c/       D:/  -> /host/d/       E:/ (USB) -> /host/e/
 *
 *   macOS:
 *     /Users           -> /host/Users
 *     /Volumes/USB-NMR -> /host/Volumes/USB-NMR
 *
 *   Linux:
 *     /home              -> /host/home
 *     /media/<user>/usb  -> /host/media/<user>/usb
 *     /mnt/data          -> /host/mnt/data
 *     /run/media/...     -> /host/run/media/...
 *
 * This function converts paths the user pastes in the UI to the
 * corresponding path inside the container.
 */
function normalizeHostPathToContainer(raw) {
  if (!raw) return '';
  let p = toForward(raw);

  // Already normalized
  if (p.startsWith('/host/')) return p;

  // Windows: C:/anything -> /host/c/anything
  // Covers all drives (C:, D:, E: USB, etc.)
  if (/^[A-Za-z]:\//.test(p)) {
    return `/host/${p[0].toLowerCase()}/${lstripSlash(p.slice(3))}`;
  }

  // macOS & Linux absolute paths: prepend /host
  // Covers: /Users/..., /Volumes/..., /home/..., /media/..., /mnt/..., /run/media/...
  if (p.startsWith('/')) {
    return `/host${p}`;
  }

  // Relative path fallback
  return `/host/home/${lstripSlash(p)}`;
}

module.exports = { normalizeHostPathToContainer };