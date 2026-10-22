import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from read_breakdown_csv import bowtie_counts, build_table, percentage
from bacteria_with_decoys_csvConversion import summary_rows


class ReportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.output = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)

    def fixture(self, paired=False):
        (self.output / 'cutadapt.txt').write_text(
            'Total read pairs processed: 20\nPairs written (passing filters): 18\n' if paired else
            'Total reads processed: 20\nReads written (passing filters): 18\n')
        specs = [
            ('decoy', 'decoy', 'decoy', 2, None),
            ('bacteria_rna', 'bacteria', 'rna', 5, [1, 3, 1, 0]),
            ('host_rna', 'host', 'rna', 3, [0, 1, 2, 0]),
            ('bacteria_dna', 'bacteria', 'dna', 4, [0, 1, 1, 2]),
            ('host_dna', 'host', 'dna', 2, [1, 0, 0, 1]),
        ]
        remaining, stages = 18, []
        for name, organism, kind, aligned, categories in specs:
            if paired:
                log = f'{remaining} reads; of these:\n  {remaining} (100.00%) were paired; of these:\n    {remaining-aligned} (0.00%) aligned concordantly 0 times\n'
            else:
                log = f'{remaining} reads; of these:\n  {remaining} (100.00%) were unpaired; of these:\n    {remaining-aligned} (0.00%) aligned 0 times\n'
            (self.output / f'{name}.log').write_text(log)
            stage = dict(name=name, organism=organism, kind=kind, log=f'{name}.log')
            if categories is not None:
                stage['categories'] = f'{name}.json'
                data = dict(aligned=aligned, categories=dict(zip(('rRNA', 'CDS', 'other_RNA', 'intergenic_unannotated'), categories)),
                            cds_features_with_reads=['a'] if kind == 'rna' else ['a', 'b'])
                (self.output / f'{name}.json').write_text(json.dumps(data))
            if name == 'bacteria_dna':
                stage['coverage'] = 'coverage.txt'
            remaining -= aligned
            stages.append(stage)
        (self.output / 'coverage.txt').write_text('#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\nchr1\t1\t100\t1\t20\t20\t1\t40\t42\nchr2\t1\t200\t1\t60\t30\t2\t40\t42\n')
        manifest = {'samples': [dict(sample='sample', cutadapt_log='cutadapt.txt', stages=stages, leftover=[])]}
        self.save_manifest(manifest)
        return manifest

    def save_manifest(self, manifest):
        (self.output / 'run_manifest.json').write_text(json.dumps(manifest))

    def test_bowtie_empty_input_terse_report(self):
        path = self.output / 'empty.log'
        path.write_text('0 reads\n0.00% overall alignment rate\n')
        self.assertEqual(bowtie_counts(path), (0, 0))
        warnings = '[WARNING] Failed to launch x86-64-v3 version, staying with default\n' * 2
        path.write_text(warnings + '0 reads\n0.00% overall alignment rate\n')
        self.assertEqual(bowtie_counts(path), (0, 0))
        path.write_text('1 reads\n0.00% overall alignment rate\n')
        with self.assertRaises(ValueError):
            bowtie_counts(path)

    def test_single_read_categories_conserve_input(self):
        self.fixture()
        sample, counts, sources = build_table(self.output)[0]
        self.assertEqual(sample, 'sample')
        self.assertEqual(counts['Bacteria aligned'], 9)
        self.assertEqual([counts[f'Bacteria {category}'] for category in ('rRNA', 'CDS', 'other RNA', 'intergenic/unannotated')], [1, 4, 2, 2])
        self.assertEqual([counts[f'Host {category}'] for category in ('rRNA', 'CDS', 'other RNA', 'intergenic/unannotated')], [1, 1, 2, 1])
        self.assertEqual(counts['Bacteria non-rRNA'], 8)
        self.assertEqual(counts['Leftover reads'], 2)
        self.assertIn('bacteria_rna.json', sources['Bacteria CDS'])
        self.assertIn('bacteria_dna.json', sources['Bacteria CDS'])

    def test_pairs_counted_as_pairs_in_reports(self):
        self.fixture(paired=True)
        _, counts, _ = build_table(self.output)[0]
        self.assertEqual(counts['Bacteria aligned'], 9)
        self.assertEqual(counts['Total reads'], 20)
        self.assertEqual(percentage(counts['Bacteria CDS'], counts['Total reads']), '20.00%')

    def test_cli_writes_canonical_and_legacy_breakdown(self):
        self.fixture()
        subprocess.run([sys.executable, str(ROOT / 'read_breakdown_csv.py'), '--output-dir', str(self.output)], check=True)
        canonical = self.output / 'output_read_breakdown.csv'
        self.assertEqual(canonical.read_bytes(), (self.output / f'BACTERIA_{self.output.name}_read_breakdown.csv').read_bytes())
        with canonical.open() as handle:
            row = list(csv.DictReader(handle))[0]
        self.assertEqual(row['Host CDS (reads)'], '1')
        self.assertEqual(row['Host other RNA (% of total reads)'], '10.00%')

    def test_summary_aggregates_rna_dna_and_unique_cds_ids(self):
        self.fixture()
        self.assertEqual(summary_rows(self.output)[0], ('sample', 20, 18, 9, '26.67%', '1.6667', 4, 2))

    def test_inconsistent_stage_inputs_are_rejected(self):
        self.fixture()
        p = self.output / 'bacteria_dna.log'
        p.write_text(p.read_text().replace('8 reads;', '9 reads;'))
        with self.assertRaisesRegex(ValueError, 'Stage totals'):
            build_table(self.output)

    def test_mismatched_category_sum_is_rejected(self):
        self.fixture()
        p = self.output / 'bacteria_rna.json'
        data = json.loads(p.read_text())
        data['categories']['CDS'] += 1
        p.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, 'Category totals'):
            build_table(self.output)

    def test_missing_and_negative_categories_are_rejected(self):
        for invalid in ('missing', 'negative'):
            with self.subTest(invalid=invalid):
                self.fixture()
                p = self.output / 'host_rna.json'
                data = json.loads(p.read_text())
                if invalid == 'missing':
                    del data['categories']['CDS']
                else:
                    data['categories']['CDS'] = -1
                p.write_text(json.dumps(data))
                with self.assertRaises(ValueError):
                    build_table(self.output)

    def test_dna_before_rna_is_rejected(self):
        manifest = self.fixture()
        stages = manifest['samples'][0]['stages']
        stages[1]['kind'] = 'dna'
        self.save_manifest(manifest)
        with self.assertRaisesRegex(ValueError, 'RNA stage follows'):
            build_table(self.output)

    def test_host_disabled_has_zero_categories(self):
        manifest = self.fixture()
        stages = manifest['samples'][0]['stages']
        manifest['samples'][0]['stages'] = [stage for stage in stages if stage['organism'] != 'host']
        (self.output / 'bacteria_dna.log').write_text('11 reads; of these:\n  7 (63.64%) aligned 0 times\n')
        self.save_manifest(manifest)
        _, counts, sources = build_table(self.output)[0]
        self.assertEqual(counts['Host aligned'], 0)
        self.assertEqual(counts['Host CDS'], 0)
        self.assertEqual(counts['Leftover reads'], 7)
        self.assertEqual(sources['Host CDS'], 'Host alignment not configured')


if __name__ == '__main__':
    unittest.main()
