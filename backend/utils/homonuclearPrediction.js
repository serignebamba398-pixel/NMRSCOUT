// Shared validation for HTTP requests and extractor invocations.
function normalizeHomonuclearOptions(experiment, options = {}) {
  const fail = message => { const error = new Error(`${experiment}: ${message}`); error.statusCode = 400; throw error; };
  const integer = (value, fallback, name) => {
    const n = value == null ? fallback : Number(value);
    if (typeof value === 'boolean' || !Number.isInteger(n) || n < 1) fail(`${name} must be an integer >= 1`);
    return n;
  };
  const result = { max_bonds: integer(options.max_bonds, 3, 'max_bonds') };
  result.min_bonds = integer(options.min_bonds, 2, 'min_bonds');
  if (result.min_bonds > result.max_bonds) fail('min_bonds must be <= max_bonds');
  const diagonal = options.include_diagonal ?? true;
  if (![true, false, 'true', 'false'].includes(diagonal)) fail('include_diagonal must be boolean');
  result.include_diagonal = diagonal === true || diagonal === 'true';
  return result;
}

function homonuclearArgs(experiment, options) {
  const config = normalizeHomonuclearOptions(experiment, options);
  const args = ['--experiment', experiment, '--max_bonds', String(config.max_bonds)];
  args.push('--min_bonds', String(config.min_bonds));
  args.push(config.include_diagonal ? '--diagonal' : '--no_diagonal');
  return args;
}

module.exports = { normalizeHomonuclearOptions, homonuclearArgs };
