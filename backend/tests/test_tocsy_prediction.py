"""Run with an environment containing RDKit, NumPy and pandas."""
import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pandas as pd
from rdkit import Chem

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
sys.path.insert(0, str(SCRIPTS))
import predict_COSY_real as predictor


class TocsyPredictionTest(unittest.TestCase):
    def molecule(self, smiles):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        h2g, g2h = predictor.get_h_groups(mol)
        return mol, h2g, g2h

    def systems(self, mol, g2h):
        adjacency = predictor.build_tocsy_spin_graph(mol, g2h)
        return (adjacency, *predictor.find_spin_systems(adjacency))

    def pairs(self, mol, h2g, g2h, max_bonds=3, diagonal=True, min_bonds=2):
        with contextlib.redirect_stdout(io.StringIO()):
            return predictor.get_tocsy_pairs(mol, h2g, g2h, min_bonds, max_bonds, diagonal)

    def test_chain_depth_changes_pairs_not_systems(self):
        mol, h2g, g2h = self.molecule('CCCC')
        short = self.pairs(mol, h2g, g2h, 3)
        long = self.pairs(mol, h2g, g2h, 6)
        systems = lambda pairs: {h1: sid for h1, h2, _, _, sid in pairs if h1 == h2}
        self.assertEqual(systems(short), systems(long))
        self.assertEqual(len(set(systems(long).values())), 1)
        self.assertGreater(len(long), len(short))
        self.assertTrue(all(p[2] <= 3 for p in short))
        distant = self.pairs(mol, h2g, g2h, 6, min_bonds=4)
        self.assertEqual(systems(long), systems(distant))
        self.assertTrue(all(a == b or 4 <= bonds <= 6 for a, b, bonds, *_ in distant))
        self.assertLess(len(distant), len(long))

    def test_separate_networks_never_merge(self):
        for smiles in ['CC.CCC', 'CCOCC', 'CCNCC', 'CCSCC', 'CCc1ccccc1']:
            with self.subTest(smiles=smiles):
                mol, h2g, g2h = self.molecule(smiles)
                _, system, components = self.systems(mol, g2h)
                self.assertGreaterEqual(len(components), 2)
                pairs = self.pairs(mol, h2g, g2h, 20)
                self.assertTrue(all(system[h2g[a]] == system[h2g[b]] for a, b, *_ in pairs))
                self.assertEqual(predictor.build_tocsy_spin_graph(mol, g2h, 1, 20),
                                 predictor.build_tocsy_spin_graph(mol, g2h, 2, 3))

    def test_quaternary_carbon_blocks_relay(self):
        mol, h2g, g2h = self.molecule('CCC(C)(C)CC')
        _, system, components = self.systems(mol, g2h)
        self.assertEqual(len(components), 4)
        self.assertTrue(all(system[h2g[a]] == system[h2g[b]] for a, b, *_ in self.pairs(mol, h2g, g2h, 20)))

    def test_cycle_keeps_alternative_route(self):
        mol, h2g, g2h = self.molecule('CC1(C)CCCCC1')
        carbon_h = predictor.build_parent_to_hydrogens(mol)
        a, b = carbon_h[3][0], carbon_h[7][0]
        # Shortest molecular path crosses C1; relay travels the other way round.
        pairs = self.pairs(mol, h2g, g2h, 4)
        pair = next(p for p in pairs if p[:2] == (a, b))
        self.assertEqual(pair[2], 4)
        self.assertEqual(pair[3], 4)
        _, _, components = self.systems(mol, g2h)
        self.assertEqual(len(components), 3)  # ring network + two isolated methyls

    def test_cosy_matches_original_distance_rule(self):
        for smiles in ['CCCC', 'CCOCC', 'CC1(C)CCCCC1']:
            mol, h2g, _ = self.molecule(smiles)
            for lower, upper in [(2, 3), (2, 6), (3, 3)]:
                expected = {(h, h, 0, 0, None) for h in h2g}
                for a in h2g:
                    for b in h2g:
                        if a == b or h2g[a] == h2g[b]:
                            continue
                        route = Chem.GetShortestPath(mol, a, b)
                        if route and lower <= len(route) - 1 <= upper:
                            expected.add((a, b, len(route) - 1, 1, None))
                self.assertEqual(set(predictor.get_cosy_pairs(mol, h2g, lower, upper)), expected)

    def test_diagonal_symmetry_and_exchangeable_exclusion(self):
        mol, h2g, g2h = self.molecule('CCCO')
        self.assertTrue(all(mol.GetAtomWithIdx(h).GetNeighbors()[0].GetSymbol() == 'C' for h in h2g))
        for experiment in ['COSY', 'TOCSY']:
            for diagonal in [False, True]:
                pairs = (predictor.get_cosy_pairs(mol, h2g, include_diagonal=diagonal)
                         if experiment == 'COSY' else self.pairs(mol, h2g, g2h, diagonal=diagonal))
                self.assertEqual(sum(a == b for a, b, *_ in pairs), len(h2g) if diagonal else 0)
                pair_set = set(pairs)
                for a, b, bonds, relay, system in pairs:
                    self.assertIn((b, a, bonds, relay, system), pair_set)

    def test_1h_cli_single_file_and_folder(self):
        with tempfile.TemporaryDirectory() as folder:
            mol, _, _ = self.molecule('CCO')
            mol.SetProp('pred_shifts', ' '.join(str(atom.GetIdx() / 10) for atom in mol.GetAtoms()))
            sdf = Path(folder) / 'ethanol.sdf'
            with Chem.SDWriter(str(sdf)) as writer:
                writer.write(mol)
            for input_flag, input_path in [('--sdf_path', str(sdf)), ('--input_folder', folder)]:
                for keep_exchangeable in [False, True]:
                    with self.subTest(input_flag=input_flag, keep_exchangeable=keep_exchangeable):
                        output = Path(folder) / 'prediction_1H_test.csv'
                        args = [sys.executable, str(SCRIPTS / 'predict_1H_real.py'),
                                input_flag, input_path, '--output_csv', str(output)]
                        if keep_exchangeable:
                            args.append('--keep_exchangeable')
                        result = subprocess.run(args, capture_output=True, text=True)
                        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                        table = pd.read_csv(output)
                        self.assertEqual(len(table), 6 if keep_exchangeable else 5)
                        self.assertTrue(table.Delta_ppm.notna().all())

    def test_cli_csv_metadata_and_validation(self):
        with tempfile.TemporaryDirectory() as folder:
            mol, _, _ = self.molecule('CCCC')
            mol.SetProp('pred_shifts', ' '.join(str(atom.GetIdx() / 10) for atom in mol.GetAtoms()))
            sdf = Path(folder) / 'sample.sdf'
            with Chem.SDWriter(str(sdf)) as writer:
                writer.write(mol)
            base = [sys.executable, str(SCRIPTS / 'predict_COSY_real.py'), '--sdf_path', str(sdf), '--output_dir', folder, '--task_id', 'test']
            def run(extra):
                return subprocess.run(base + extra, capture_output=True, text=True)
            result = run(['--experiment', 'TOCSY', '--max_bonds', '6'])
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(json.loads(result.stdout.splitlines()[-1])['success'])
            output = Path(folder) / 'prediction_TOCSY_test.csv'
            df = pd.read_csv(output)
            self.assertEqual(df.columns.tolist(), ['Filename', 'H1', 'Delta1', 'H2', 'Delta2', 'Bonds', 'Group1', 'Group2', 'ECFP4', 'Method', 'RelaySteps', 'SpinSystem', 'MaxBonds'])
            self.assertTrue((df.MaxBonds == 6).all())
            self.assertEqual(run([]).returncode, 0)  # legacy COSY invocation
            self.assertTrue((Path(folder) / 'prediction_COSY_test.csv').exists())
            for extra in [['--max_bonds', '0'], ['--min_bonds', '7', '--max_bonds', '6']]:
                self.assertNotEqual(run(extra).returncode, 0)


if __name__ == '__main__':
    unittest.main()
