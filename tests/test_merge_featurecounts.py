import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from merge_featurecounts import ANNOTATION_COLUMNS, TABLE_PATHS, merge_featurecounts


class MergeFeatureCountsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.samples = [self.base / name for name in ('sample-a', 'sample-b')]
        for directory in self.samples:
            directory.mkdir()
        self.output = self.base / 'staging'

    def table(self, directory, relative=TABLE_PATHS[0], counts=(2, 0), sam=None,
              annotations=None, summary=True):
        sam = sam or str(directory / 'BACTERIA_sample.sam')
        path = directory / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        if annotations is None:
            annotations = [('gene1', 'chr1', '1', '100', '+', '100'),
                           ('gene2', 'chr1', '201', '300', '-', '100')]
        with path.open('w', newline='') as handle:
            handle.write('# Program:featureCounts v2.1.1; Command: ...\n')
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow((*ANNOTATION_COLUMNS, sam))
            for annotation, count in zip(annotations, counts):
                writer.writerow((*annotation, count))
        if summary:
            Path(str(path) + '.summary').write_text(f'Status\t{sam}\nAssigned\t{sum(counts)}\nUnassigned_Unmapped\t7\n')
        return path

    def read_rows(self, path, delimiter='\t'):
        with path.open(newline='') as handle:
            return list(csv.reader((line for line in handle if not line.startswith('#')), delimiter=delimiter))

    def test_combines_all_stage_matrices_and_summaries_in_sample_order(self):
        originals = {}
        for relative in TABLE_PATHS:
            for i, directory in enumerate(self.samples):
                path = self.table(directory, relative, counts=(i + 1, i * 3))
                originals[path] = path.read_bytes()
        generated = merge_featurecounts(self.samples, self.output)
        self.assertEqual(len(generated), 20)  # 8 matrices, 8 summaries, 4 CDS CSVs
        self.assertTrue(all(path.is_absolute() and path.is_file() for path in generated))
        for relative in TABLE_PATHS:
            rows = self.read_rows(self.output / relative)
            self.assertEqual(rows[0][:6], list(ANNOTATION_COLUMNS))
            self.assertEqual(rows[0][6:], [str(directory / 'BACTERIA_sample.sam') for directory in self.samples])
            self.assertEqual(rows[1][6:], ['1', '2'])
            self.assertEqual(rows[2][6:], ['0', '3'])
            summary = self.read_rows(Path(str(self.output / relative) + '.summary'))
            self.assertEqual(summary[0], ['Status', *rows[0][6:]])
            self.assertEqual(summary[1], ['Assigned', '1', '5'])
            if relative.name.endswith('_summary.txt'):
                self.assertEqual(self.read_rows((self.output / relative).with_suffix('.csv'), ','), rows)
        for path, content in originals.items():
            self.assertEqual(path.read_bytes(), content)

    def test_missing_stage_for_all_samples_is_allowed(self):
        self.assertEqual(merge_featurecounts(self.samples, self.output), [])
        for directory in self.samples:
            self.table(directory, summary=False)
        generated = merge_featurecounts(self.samples, self.output)
        self.assertEqual({path.name for path in generated}, {'featurecounts_BACTERIA_summary.txt', 'featurecounts_BACTERIA_summary.csv'})

    def test_partial_stage_or_summary_is_rejected_before_output(self):
        self.table(self.samples[0])
        with self.assertRaisesRegex(ValueError, 'Missing featureCounts tables'):
            merge_featurecounts(self.samples, self.output)
        self.assertFalse(self.output.exists())
        self.table(self.samples[1], summary=False)
        with self.assertRaisesRegex(ValueError, 'summary companions'):
            merge_featurecounts(self.samples, self.output)
        self.assertFalse(self.output.exists())

    def test_disagreed_annotation_order_preserves_existing_output(self):
        self.table(self.samples[0])
        self.table(self.samples[1], annotations=[('gene2', 'chr1', '201', '300', '-', '100'),
                                                ('gene1', 'chr1', '1', '100', '+', '100')])
        self.output.mkdir()
        target = self.output / TABLE_PATHS[0]
        target.write_text('existing output\n')
        with self.assertRaisesRegex(ValueError, 'Mismatched featureCounts annotation'):
            merge_featurecounts(self.samples, self.output)
        self.assertEqual(target.read_text(), 'existing output\n')
        self.assertEqual(sorted(path.name for path in self.output.iterdir()), [target.name])

    def test_gene_row_count_mismatch_is_rejected(self):
        self.table(self.samples[0])
        self.table(self.samples[1], counts=(4,))
        with self.assertRaisesRegex(ValueError, 'Different gene row counts'):
            merge_featurecounts(self.samples, self.output)

    def test_multiple_sample_columns_and_invalid_counts_are_rejected(self):
        for invalid in ('extra_column', 'negative_count'):
            with self.subTest(invalid=invalid):
                for directory in self.samples:
                    self.table(directory)
                path = self.samples[1] / TABLE_PATHS[0]
                text = path.read_text()
                if invalid == 'extra_column':
                    text = text.replace('Length\t', 'Length\tsecond_sample\t')
                else:
                    text = text.replace('\t2\n', '\t-2\n')
                path.write_text(text)
                with self.assertRaises(ValueError):
                    merge_featurecounts(self.samples, self.output)

    def test_duplicate_sample_columns_and_directories_are_rejected(self):
        for directory in self.samples:
            self.table(directory, sam='same.sam')
        with self.assertRaisesRegex(ValueError, 'Duplicate featureCounts sample columns'):
            merge_featurecounts(self.samples, self.output)
        with self.assertRaisesRegex(ValueError, 'Duplicate sample output directories'):
            merge_featurecounts([self.samples[0], self.samples[0]], self.output)

    def test_summary_header_and_status_must_match(self):
        for directory in self.samples:
            self.table(directory)
        path = Path(str(self.samples[1] / TABLE_PATHS[0]) + '.summary')
        original = path.read_text()
        path.write_text(original.replace('BACTERIA_sample.sam', 'OTHER.sam'))
        with self.assertRaisesRegex(ValueError, 'sample column disagrees'):
            merge_featurecounts(self.samples, self.output)
        path.write_text(original.replace('Assigned\t', 'WrongStatus\t'))
        with self.assertRaisesRegex(ValueError, 'Mismatched or duplicate featureCounts status'):
            merge_featurecounts(self.samples, self.output)

    def test_empty_gene_table_still_contains_all_sample_headers(self):
        for directory in self.samples:
            self.table(directory, counts=(), annotations=[])
        merge_featurecounts(self.samples, self.output)
        self.assertEqual(len(self.read_rows(self.output / TABLE_PATHS[0])), 1)
        self.assertEqual(len(self.read_rows(self.output / TABLE_PATHS[0])[0]), 8)

    def test_cli_includes_every_supplied_sample(self):
        for i, directory in enumerate(self.samples):
            self.table(directory, counts=(i, i + 1))
        subprocess.run([sys.executable, str(ROOT / 'merge_featurecounts.py'),
                        '--sample-dir', str(self.samples[0]), '--sample-dir', str(self.samples[1]),
                        '--output-dir', str(self.output)], check=True)
        self.assertEqual(self.read_rows(self.output / TABLE_PATHS[0])[1][6:], ['0', '1'])

    def test_many_samples_keep_every_count_column(self):
        directories = [self.base / f'worker-{i}' for i in range(40)]
        for i, directory in enumerate(directories):
            directory.mkdir()
            self.table(directory, counts=(i, 40 - i))
        merge_featurecounts(directories, self.output)
        rows = self.read_rows(self.output / TABLE_PATHS[0])
        self.assertEqual(rows[1][6:], [str(i) for i in range(40)])
        self.assertEqual(rows[2][6:], [str(40 - i) for i in range(40)])


if __name__ == '__main__':
    unittest.main()
