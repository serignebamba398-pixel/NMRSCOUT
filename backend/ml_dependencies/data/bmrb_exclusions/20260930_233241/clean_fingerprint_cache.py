"""Subtract audited BMRB entries without rebuilding other database sources."""
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import pickle
import re
import shutil
import sys
import tempfile
from collections import Counter, defaultdict

AUDIT = Path(__file__).resolve().parent
ROOT = AUDIT.parents[4]
sys.path.insert(0, str(ROOT / 'backend' / 'scripts'))
import fingerprint_db_predictor as predictor


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def new_index():
    return {r: defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
            for r in predictor.SUPPORTED_RADII}


def main():
    cache = Path(predictor.DB_CACHE_FILE)
    backup = AUDIT / 'fingerprint_database.before_exclusions.pkl'
    report_path = AUDIT / 'fingerprint_cleanup_report.json'
    if report_path.exists() or backup.exists():
        raise RuntimeError('Previous cleanup artifacts exist; refusing to subtract twice')
    before_hash = digest(cache)
    with cache.open('rb') as stream:
        data = pickle.load(stream)
    assert data['version'] == '13.0'
    assert data.get('bmrb_ambiguity_revision', 0) >= predictor.REVISION
    assert predictor.BEYOND_AVAILABLE
    # Bypass __init__: it can automatically trigger a complete rebuild.
    engine = object.__new__(predictor.FingerprintPredictor)
    removed = new_index()
    composite = new_index()
    exact_candidates = defaultdict(set)
    blocks = re.findall(rb'.*?^\$\$\$\$[^\S\r\n]*(?:\r?\n|$)',
                        (AUDIT / 'excluded_entries.sdf').read_bytes(), re.M | re.S)
    assert len(blocks) == 191
    entry_counts = []
    with tempfile.TemporaryDirectory(prefix='bmrb_replay_', dir=AUDIT) as staging:
        predictor.BMRB_DIR = staging
        for i, block in enumerate(blocks, 1):
            (Path(staging) / 'entry.sdf').write_bytes(block)
            exact = {}
            with contextlib.redirect_stdout(io.StringIO()):
                count = engine._load_bmrb(removed, {}, exact, db_composite=composite)
            entry_counts.append(count)
            for key, solvents in exact.items():
                for solvent, nuclei in solvents.items():
                    for nucleus, atoms in nuclei.items():
                        for atom, value in atoms.items():
                            exact_candidates[key, solvent, nucleus, atom].add(value)
            if i % 25 == 0:
                print(f'Replayed {i}/{len(blocks)} excluded entries', flush=True)

    changes = []
    totals = {}
    # Verify every required multiplicity before changing any cache bucket.
    for name, contribution in [('indexes', removed), ('indexes_composite', composite)]:
        per_radius = Counter()
        for radius, fps in contribution.items():
            for fp, solvents in fps.items():
                for solvent, nuclei in solvents.items():
                    for nucleus, values in nuclei.items():
                        bucket = data.get(name, {}).get(radius, {}).get(fp, {}).get(solvent, {})
                        old = bucket.get(nucleus, [])
                        required = Counter(values)
                        if required - Counter(old):
                            raise RuntimeError(f'Missing contributions: {name}, R{radius}, {solvent}, {nucleus}; cache unchanged')
                        pending = required.copy()
                        kept = []
                        for value in old:
                            if pending[value]:
                                pending[value] -= 1
                            else:
                                kept.append(value)
                        assert Counter(old) == Counter(kept) + required
                        changes.append((bucket, nucleus, kept))
                        per_radius[radius] += len(values)
        totals[name] = dict(per_radius)
    print('All fingerprint contribution counts verified', totals, flush=True)

    exact_changes = []
    for (key, solvent, nucleus, atom), values in exact_candidates.items():
        bucket = data.get('exact_db', {}).get(key, {}).get(solvent, {}).get(nucleus, {})
        if atom in bucket and bucket[atom] in values:
            exact_changes.append((bucket, atom))
    for bucket, nucleus, kept in changes:
        bucket[nucleus] = kept
    for bucket, atom in exact_changes:
        del bucket[atom]
    for name in ['indexes', 'indexes_composite']:
        for fps in data[name].values():
            for fp in list(fps):
                for solvent in list(fps[fp]):
                    for nucleus in list(fps[fp][solvent]):
                        if not fps[fp][solvent][nucleus]:
                            del fps[fp][solvent][nucleus]
                    if not fps[fp][solvent]:
                        del fps[fp][solvent]
                if not fps[fp]:
                    del fps[fp]
    for key in list(data['exact_db']):
        solvents = data['exact_db'][key]
        for solvent in list(solvents):
            for nucleus in list(solvents[solvent]):
                if not solvents[solvent][nucleus]:
                    del solvents[solvent][nucleus]
            if not solvents[solvent]:
                del solvents[solvent]
        if not solvents:
            del data['exact_db'][key]
    active_fps = set().union(*(set(fps) for fps in data['indexes'].values()))
    data['fp_to_env'] = {fp: env for fp, env in data['fp_to_env'].items() if fp in active_fps}
    data['stats'] = {
        'fps_by_radius': {r: len(v) for r, v in data['indexes'].items()},
        'composite_by_radius': {r: len(v) for r, v in data['indexes_composite'].items()},
    }
    data['total_records'] = sum(data['stats']['fps_by_radius'].values())
    report = {
        'excluded_entries': len(blocks),
        'entries_with_no_replayed_shifts': sum(v == 0 for v in entry_counts),
        'removed_values_by_index_and_radius': totals,
        'removed_exact_atom_values': len(exact_changes),
        'exact_removal_policy': 'Remove only atom values matching excluded entries; provenance is not stored, so identical values from other sources cannot be distinguished.',
        'before_sha256': before_hash,
        'stats_after': data['stats'],
    }
    data['bmrb_exclusion_audit'] = report.copy()
    assert digest(cache) == before_hash, 'Cache changed during audit'
    shutil.copy2(cache, backup)
    assert digest(backup) == before_hash
    temporary = cache.with_suffix('.cleaning.tmp')
    with temporary.open('wb') as stream:
        pickle.dump(data, stream, protocol=pickle.HIGHEST_PROTOCOL)
    # Reload the serialized file and validate affected buckets before replacement.
    del changes, exact_changes, data
    with temporary.open('rb') as stream:
        verified = pickle.load(stream)
    assert verified['bmrb_exclusion_audit'] == report
    assert all(verified['indexes'].get(r) for r in predictor.SUPPORTED_RADII)
    for (key, solvent, nucleus, atom), values in exact_candidates.items():
        bucket = verified['exact_db'].get(key, {}).get(solvent, {}).get(nucleus, {})
        assert atom not in bucket or bucket[atom] not in values
    assert digest(cache) == before_hash, 'Cache changed before replacement'
    os.replace(temporary, cache)
    report['after_sha256'] = digest(cache)
    report['validation'] = 'PASS: multiplicities, backup, serialized reload, exact exclusions'
    report_path.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
