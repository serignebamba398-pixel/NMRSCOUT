function toForward(p) {
  return String(p || '')
    .trim()
    .replace(/\\/g, '/');
}

function lstripSlash(p) {
  return String(p).replace(/^\/+/, '');
}

const MAPS = [
  // Windows "C:/..." -> /host/c/...
  {
    test: p => /^[A-Za-z]:\//.test(p),
    map: p => `/host/${p[0].toLowerCase()}/${lstripSlash(p.slice(3))}`,
  },
  // Windows profils "*/Users/<login>/*" -> /host/winhome/...
  {
    test: p => /\/users\//i.test(p),
    map: p => p.replace(/^.*?\/users\//i, '/host/winhome/'),
  },
  // Linux/macOS $HOME -> /host/home/...
  {
    test: p => {
      const home = process.env.HOME || '/home';
      return p === home || p.startsWith(home + '/');
    },
    map: p => p.replace(process.env.HOME || '/home', '/host/home'),
  },
  // macOS absolu: /Users/<login>/... -> /host/home/...
  {
    test: p => /^\/Users\/[^/]+\//.test(p),
    map: p => p.replace(/^\/Users\/[^/]+/, '/host/home'),
  },
];

function normalizeHostPathToContainer(raw) {
  if (!raw) return '';
  let p = toForward(raw);

  if (p.startsWith('/host/')) return p;

  if (/^[A-Za-z]:\//.test(p)) {
    return `/host/${p[0].toLowerCase()}/${lstripSlash(p.slice(3))}`;
  }

  if (!p.startsWith('/')) {
    return `/host/home/${lstripSlash(p)}`;
  }

  for (const m of MAPS) {
    if (m.test(p)) return m.map(p);
  }

  return p;
}

module.exports = { normalizeHostPathToContainer };