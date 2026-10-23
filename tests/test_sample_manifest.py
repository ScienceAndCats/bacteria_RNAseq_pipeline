import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sample_manifest import build_manifest, discover_samples


class SampleManifestTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="sample manifest ")
        self.root = Path(self.temporary.name)
        self.addCleanup(self.temporary.cleanup)

    def touch(self, *names):
        for name in names:
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()

    def test_single_sorted_names_spaces_and_ignored_non_fastq_matches(self):
        self.touch("z.fastq.gz", "a quoted 'sample.fastq", "ignore.fastq.gz.log", "README")
        records = discover_samples(self.root, "*", "single")
        self.assertEqual([r["sample"] for r in records], ["a quoted 'sample", "z"])
        self.assertTrue(all(Path(r["read1"]).is_absolute() and r["read2"] is None for r in records))

    def test_glob_characters_in_directory_and_sample_are_literal(self):
        directory = self.root / "reads [2026]"
        directory.mkdir()
        (directory / "sample [A] $name.fastq").touch()
        records = discover_samples(directory, "*.fastq")
        self.assertEqual(records[0]["sample"], "sample [A] $name")

    def test_paired_suffixes_are_preserved_once_and_mixed_compression_allowed(self):
        self.touch("B_R2_001.fastq", "B_R1_001.fastq.gz", "A_R1.fastq", "A_R2.fastq")
        records = discover_samples(self.root, "*.fastq*", "paired")
        self.assertEqual([r["sample"] for r in records], ["A", "B_001"])
        self.assertEqual(Path(records[1]["read1"]).name, "B_R1_001.fastq.gz")
        self.assertEqual(Path(records[1]["read2"]).name, "B_R2_001.fastq")

    def test_r1_only_glob_discovers_existing_unselected_r2(self):
        self.touch("A_R1.fastq.gz", "A_R2.fastq.gz")
        self.assertEqual(len(discover_samples(self.root, "*_R1.fastq.gz", "paired")), 1)

    def test_r2_only_glob_discovers_existing_unselected_r1(self):
        self.touch("A_R1.fastq.gz", "A_R2.fastq.gz")
        self.assertEqual(len(discover_samples(self.root, "*_R2.fastq.gz", "paired")), 1)

    def test_orphan_and_wrong_suffix_are_rejected(self):
        self.touch("A_R1_001.fastq", "A_R2_002.fastq")
        with self.assertRaisesRegex(ValueError, "Missing R2"):
            discover_samples(self.root, "*", "paired")

    def test_orphan_r2_is_rejected(self):
        self.touch("A_R2.fastq")
        with self.assertRaisesRegex(ValueError, "Missing R1"):
            discover_samples(self.root, "*", "paired")

    def test_unmatched_selected_fastq_is_not_silently_ignored(self):
        self.touch("A_R1.fastq", "A_R2.fastq", "unpaired.fastq")
        with self.assertRaisesRegex(ValueError, "no valid _R1/_R2"):
            discover_samples(self.root, "*", "paired")

    def test_duplicate_single_plain_and_gzip_are_rejected(self):
        self.touch("A.fastq", "A.fastq.gz")
        with self.assertRaisesRegex(ValueError, "Duplicate sample"):
            discover_samples(self.root)

    def test_duplicate_paired_plain_and_gzip_are_rejected(self):
        self.touch("A_R1.fastq", "A_R1.fastq.gz", "A_R2.fastq")
        with self.assertRaisesRegex(ValueError, "Duplicate R1"):
            discover_samples(self.root, "*", "paired")

    def test_distinct_pairs_with_same_normalized_name_are_rejected(self):
        self.touch("A_R1_001.fastq", "A_R2_001.fastq", "A_001_R1.fastq", "A_001_R2.fastq")
        with self.assertRaisesRegex(ValueError, "Duplicate sample"):
            discover_samples(self.root, "*", "paired")

    def test_same_basename_in_different_folders_is_rejected(self):
        self.touch("one/A.fastq", "two/A.fastq")
        with self.assertRaisesRegex(ValueError, "Duplicate sample"):
            discover_samples(self.root, "**/*.fastq", "single")

    def test_control_character_and_empty_stem_rejected(self):
        self.touch("bad\nname.fastq")
        with self.assertRaisesRegex(ValueError, "control character"):
            discover_samples(self.root)
        (self.root / "bad\nname.fastq").unlink()
        self.touch(".fastq")
        with self.assertRaisesRegex(ValueError, "empty"):
            discover_samples(self.root, ".fastq")

    def test_empty_and_missing_directory_rejected(self):
        with self.assertRaisesRegex(ValueError, "No FASTQ"):
            discover_samples(self.root)
        with self.assertRaisesRegex(ValueError, "does not exist"):
            discover_samples(self.root / "missing")

    def test_bad_layout_and_escaping_globs_rejected(self):
        for layout in ("mixed", "", None):
            with self.assertRaisesRegex(ValueError, "READ_LAYOUT"):
                discover_samples(self.root, layout=layout)
        for pattern in ("../*.fastq", "/tmp/*.fastq"):
            with self.assertRaisesRegex(ValueError, "within FASTQ_DIR"):
                discover_samples(self.root, pattern=pattern)

    def test_maximum_limit_rejects_before_manifest_written(self):
        self.touch(*(f"sample_{i}.fastq" for i in range(41)))
        with self.assertRaisesRegex(ValueError, "41 samples.*maximum of 40"):
            build_manifest(self.root, max_samples=40)
        for value in (0, -1, True):
            with self.assertRaisesRegex(ValueError, "positive integer"):
                build_manifest(self.root, max_samples=value)

    def test_cli_writes_documented_schema(self):
        self.touch("test.fastq")
        destination = self.root / "run files/samples.json"
        script = Path(__file__).resolve().parents[1] / "sample_manifest.py"
        result = subprocess.run([sys.executable, str(script), "--fastq-dir", str(self.root), "--glob", "*.fastq", "--layout", "single", "--output", str(destination), "--max-samples", "40"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        manifest = json.loads(destination.read_text())
        self.assertEqual(manifest, {"read_layout": "single", "samples": [{"sample": "test", "read1": str(self.root / "test.fastq"), "read2": None}]})


if __name__ == "__main__":
    unittest.main()
