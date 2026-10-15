import json
import tempfile
import unittest
from pathlib import Path

from minisv.trio_finalize import classify, process_locus, run


def read(name, distance, parent_distance=100):
    return dict(name=name, length=1000, distances={'C0.0': distance, 'F0.0': parent_distance, 'M0.0': parent_distance})


def fixture():
    parents = [read(f'p{i}', 100, 0) for i in range(10)]
    return dict(child='child', trid='chr1_100_200_TR', models={'C0.0': 'A'*100, 'F0.0': 'A'*160, 'M0.0': 'A'*160},
                full_reads={'C': [read('c1', 30), read('c2', 40), read('c3', 60)], 'F': parents, 'M': parents},
                transmissions=[dict(total_edit_distance=60, paternal_child='C0.0', father='F0.0', maternal_child='C0.0', mother='M0.0')])


class FinalizeTests(unittest.TestCase):
    def test_threshold_boundary(self):
        data = fixture()
        self.assertEqual(classify(data, 'C0.0', -60, .05)['status'], 'REVIEW')
        self.assertEqual(classify(data, 'C0.0', -60, .06)['status'], 'RULE_PASS')

    def test_parent_rechecked_at_each_rate(self):
        data = fixture()
        data['full_reads']['F'][0] = read('parent_alt', 70)
        self.assertEqual(classify(data, 'C0.0', -60, .06)['status'], 'RULE_PASS')
        self.assertIn('PARENT_TARGET_SUPPORT', classify(data, 'C0.0', -60, .08)['reasons'])

    def test_low_parent_depth(self):
        data = fixture()
        data['full_reads']['M'] = data['full_reads']['M'][:9]
        self.assertIn('PARENT_SPANNING_DEPTH_LT_MIN', classify(data, 'C0.0', -60, .08)['reasons'])

    def test_duplicate_molecule_rejected(self):
        data = fixture()
        data['full_reads']['C'].append(data['full_reads']['C'][0])
        with self.assertRaises(ValueError):
            classify(data, 'C0.0', -60, .08)

    def test_local_only_event_and_parent_ties(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / 'local.json'
            p.write_text(json.dumps(fixture()))
            result = process_locus(('child', str(p), [.05, .06, .08]))
            self.assertEqual(len(result['events']), 1)
            event = result['events'][0]
            self.assertTrue(event['origin_ambiguous'])
            self.assertEqual(event['delta'], -60)
            self.assertEqual(event['evidence']['0.06']['status'], 'RULE_PASS')

    def test_end_to_end_without_c3_or_truth(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            local = root / 'local'
            local.mkdir()
            (local / 'x.json').write_text(json.dumps(fixture()))
            (root / 'c3.json').write_text('[]')
            manifest = root / 'manifest.json'
            manifest.write_text(json.dumps(dict(samples=[dict(sample='child', local_directory=str(local))], c3_audit=str(root / 'c3.json'))))
            summary = run(manifest, root / 'out', workers=1)
            self.assertEqual(summary['0.05']['rule_pass_events'], 0)
            self.assertEqual(summary['0.06']['rule_pass_events'], 1)
            self.assertEqual(json.loads((root / 'out/status.json').read_text())['stage'], 'COMPLETE')


if __name__ == '__main__':
    unittest.main()
