import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from merge_featurecounts import ANNOTATION_COLUMNS, TABLE_PATHS, merge_featurecounts
from csv_utils import CSV_FIELD_SIZE_LIMIT, LARGE_FIELD_REPORT, LARGE_FIELD_THRESHOLD


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
        self.assertEqual(sorted(path.name for path in self.output.iterdir()), sorted([target.name, LARGE_FIELD_REPORT]))
        self.assertIn('Status: FAILED', (self.output / LARGE_FIELD_REPORT).read_text())

    def test_large_human_annotations_merge_without_truncation_and_are_reported(self):
        long_chr = 'chr1;' * 30000 + 'LAST_CHROMOSOME'
        long_start = '101;' * 40000 + 'LAST_COORDINATE'
        annotations = [('human_large_gene', long_chr, long_start, '200', '+', '100')]
        originals = {}
        relative = Path('bowtie_alignments/host/featurecounts_HOST_summary.txt')
        for index, directory in enumerate(self.samples):
            path = self.table(directory, relative, annotations=annotations, counts=(index + 3,))
            originals[path] = path.read_bytes()
        previous = csv.field_size_limit(LARGE_FIELD_THRESHOLD)
        self.addCleanup(csv.field_size_limit, previous)
        with self.assertRaises(csv.Error):
            self.read_rows(path)
        merge_featurecounts(self.samples, self.output)
        rows = self.read_rows(self.output / relative)
        self.assertEqual(rows[1], [*annotations[0], '3', '4'])
        self.assertEqual(self.read_rows((self.output / relative).with_suffix('.csv'), ','), rows)
        self.assertEqual(csv.field_size_limit(), 64 * 1024 * 1024)
        report = (self.output / LARGE_FIELD_REPORT).read_text()
        for path, content in originals.items():
            self.assertIn(str(path), report)
            self.assertEqual(path.read_bytes(), content)
        for text in ('human_large_gene', long_chr, long_start, 'Oversized column: Chr (2)',
                     'Oversized column: Start (3)', 'Oversized cells: 4',
                     'Source records with oversized cells: 2', 'Status: COMPLETED'):
            self.assertIn(text, report)

    def test_parser_failure_keeps_complete_raw_record_in_explicit_report(self):
        value = 'chr1;' * 30000 + 'END_OF_PROBLEM_CELL'
        for directory in self.samples:
            self.table(directory, annotations=[('problem_gene', value, '1', '10', '+', '10')], counts=(2,))
        report = self.base / 'durable' / LARGE_FIELD_REPORT
        # Exercise the hard-cap error path without allocating a 64+ MiB fixture.
        with patch('csv_utils.CSV_FIELD_SIZE_LIMIT', LARGE_FIELD_THRESHOLD + 1):
            with self.assertRaisesRegex(ValueError, 'full raw record saved'):
                merge_featurecounts(self.samples, self.output, report_path=report)
        text = report.read_text()
        for expected in ('Status: FAILED', 'CSV parser errors: 1', 'problem_gene', value,
                         'Table record: 2', 'Physical lines: 3-3', 'BEGIN RAW RECORD'):
            self.assertIn(expected, text)
        self.assertFalse((self.output / TABLE_PATHS[0]).exists())

    def test_quoted_multiline_annotation_retains_comments_and_blank_lines(self):
        annotation = 'chr1;' * 30000 + '\n# annotation content, not a comment\n\nlast annotation line'
        for directory in self.samples:
            self.table(directory, annotations=[('multi_gene', annotation, '1', '10', '+', '10')], counts=(2,))
        merge_featurecounts(self.samples, self.output)
        with (self.output / TABLE_PATHS[0]).open(newline='') as handle:
            next(handle)  # merge comment
            rows = list(csv.reader(handle, delimiter='\t'))
        self.assertEqual(rows[1][1], annotation)
        self.assertIn(annotation, (self.output / LARGE_FIELD_REPORT).read_text())

    def test_parser_failure_does_not_cut_off_a_multiline_annotation(self):
        annotation = 'chr1;' * 30000 + '\nUNREAD_MULTILINE_END'
        for directory in self.samples:
            self.table(directory, annotations=[('multiline_gene', annotation, '1', '10', '+', '10')], counts=(2,))
        with patch('csv_utils.CSV_FIELD_SIZE_LIMIT', LARGE_FIELD_THRESHOLD + 1):
            with self.assertRaisesRegex(ValueError, 'full raw record saved'):
                merge_featurecounts(self.samples, self.output)
        report = (self.output / LARGE_FIELD_REPORT).read_text()
        self.assertIn(annotation, report)
        self.assertIn('UNREAD_MULTILINE_END', report)
        self.assertIn('Status: FAILED', report)

    def test_no_large_cells_report_is_refreshed_on_rerun(self):
        for directory in self.samples:
            self.table(directory)
        self.output.mkdir()
        report = self.output / LARGE_FIELD_REPORT
        report.write_text('stale oversized field report\n')
        merge_featurecounts(self.samples, self.output)
        self.assertIn('Oversized cells: 0', report.read_text())
        self.assertNotIn('stale oversized field', report.read_text())

    def test_report_cannot_overwrite_a_source_table(self):
        for directory in self.samples:
            self.table(directory)
        source = self.samples[0] / TABLE_PATHS[0]
        original = source.read_bytes()
        with self.assertRaisesRegex(ValueError, 'must not overwrite'):
            merge_featurecounts(self.samples, self.output, report_path=source)
        self.assertEqual(source.read_bytes(), original)

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
