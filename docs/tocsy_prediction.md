# TOCSY prediction: structural systems and distance limits

The implementation keeps the user's revised `backend/scripts/predict_COSY_real.py`.
As confirmed on 2026-09-26, **both min_bonds and max_bonds filter TOCSY pairs**,
after structural spin systems have been detected. This supersedes the earlier
request to ignore min_bonds in TOCSY.

## Structural model

The existing `_group_parent_carbon`, `build_tocsy_spin_graph`, `find_spin_systems`
and `shortest_relay_steps` functions are retained. Groups on the same carbon or
on directly bonded proton-bearing carbons are connected. OH/NH/SH remain
excluded. Nonprotonated carbons and heteroatom bridges do not provide relay
nodes; alternative routes around cycles remain available. No additional rule
arbitrarily splits aromatic and aliphatic regions.

This is a topological approximation, not a calculation of J couplings,
transfer efficiency or peak intensities. Changing either distance limit does
not change components or their IDs. For a retained pair, `Bonds` is the
molecular shortest-path distance; `RelaySteps` is the distance in the proton
group graph. These quantities may differ, notably around rings.

## Parameters and compatibility

- COSY retains its direct distance-window algorithm.
- TOCSY retains a pair only if both groups share a spin system and
  `min_bonds <= Bonds <= max_bonds`.
- Bounds are integers, `1 <= min_bonds <= max_bonds`; defaults remain 2 and 3.
- Diagonals remain separately selectable and have `Bonds = RelaySteps = 0`.
- No mixing-time parameter is exposed or sent to Python. Users choose distance
  limits based on their experiment. Stale saved mixing-time options are ignored.
- TOCSY appends `MaxBonds` after the existing 12 CSV columns.
  COSY columns and both output filename patterns are preserved.

Example:

```sh
python backend/scripts/predict_COSY_real.py --sdf_path molecule.sdf --output_dir results --task_id demo --experiment TOCSY --min_bonds 2 --max_bonds 6 --diagonal
```

## Existing application path

`TocsyPredictionFields` is shared by both prediction forms. `usePrediction`
transmits its values as `options.extractors.tocsy` for SMILES and CSV requests.
The existing `/api/predict/start` endpoint validates them before starting SSE
or prediction. `extractorService` forwards the CLI options to the same Python
script. Electron continues to use its existing `predict_COSY` engine command,
with `--experiment TOCSY`; Docker invokes the script directly.

The response still exposes `outputs.csv_TOCSY`. Atom-pair CSVs with `Filename`,
`Delta1` and `Delta2` are grouped for the existing shared `Scatter2D` viewer;
annotation CSV rows keep their current format. Maximum bonds is displayed for TOCSY prediction results. Per-pair metadata remain in the CSV
and in the adapted rows; no new tooltip engine was introduced.

## Verification and deployment

```sh
python backend/tests/test_tocsy_prediction.py
node backend/tests/homonuclear_prediction.test.js
npm --prefix frontend test -- --watchAll=false --runInBand --runTestsByPath src/hooks/usePrediction.tocsy.test.js
npm --prefix frontend run build
```

Python tests need RDKit, pandas and NumPy. Backend tests check argument dispatch
for both standard Python and the Electron executable; they do not build a full
Docker image or Electron installer. Redistributing this change requires rebuilding
the Python bundle and frontend for Electron, or rebuilding the Docker images.
Existing CSVs are not modified: rerun prediction to generate the new metadata
and the selected TOCSY distance window.
