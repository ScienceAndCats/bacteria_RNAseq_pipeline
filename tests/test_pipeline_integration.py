"""Real-tool integration checks; run with the pipeline environment on PATH."""
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import json
import os
from pathlib import Path
import random
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, SimpleLocation
except ImportError:
    SeqIO = None

REPOSITORY = Path(__file__).resolve().parents[1]
TOOLS = ("cutadapt", "bowtie2", "bowtie2-build", "samtools", "featureCounts")
AVAILABLE = SeqIO is not None and all(shutil.which(tool) for tool in TOOLS)


def gzwrite(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def transcript(ident, sequence, kind, gene):
    record = SeqRecord(Seq(sequence), id=ident, name=ident, description="synthetic transcript")
    record.annotations["molecule_type"] = "RNA"
    mol_type = "mRNA" if kind == "CDS" else kind
    record.features = [SeqFeature(SimpleLocation(0, len(sequence), strand=1), type="source", qualifiers={"mol_type": [mol_type]})]
    record.features.append(SeqFeature(SimpleLocation(0, len(sequence), strand=1), type=kind, qualifiers={"locus_tag": [gene]}))
    return record


@unittest.skipUnless(AVAILABLE, "real bioinformatics tools and Biopython are required")
class PipelineIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="rnaseq integration ")
        cls.base = Path(cls.temporary.name)
        cls.references = cls.base / "reference folders"
        cls.cache = cls.base / "index cache"
        rng = random.Random(72835)
        randseq = lambda n: "".join(rng.choices("ACGT", k=n))
        cls.genomes = {role: randseq(8000) for role in ("bacteria", "host")}
        cls.decoy, cls.leftover = randseq(1000), randseq(700)
        cls.transcripts = {role: {"coding": g[3000:3300] + g[3800:4100], "rrna": g[1500:2100], "other": g[2400:2700]} for role, g in cls.genomes.items()}
        # Host junction sequence also occurs in bacterial genomic DNA. RNA must win.
        g = cls.genomes["bacteria"]
        cls.genomes["bacteria"] = g[:5000] + cls.transcripts["host"]["coding"] + g[5600:]
        cls.inputs = []
        for role in ("bacteria", "host"):
            genome, seqs = cls.genomes[role], cls.transcripts[role]
            cls.inputs.extend([(f"{role}_rna_CDS", seqs["coding"], 250), (f"{role}_rna_rRNA", seqs["rrna"], 100), (f"{role}_rna_other", seqs["other"], 10), (f"{role}_dna_CDS", genome, 200), (f"{role}_dna_intergenic", genome, 6500)])
            folder = cls.references / ("PAO1" if role == "bacteria" else "human")
            dna_folder = folder if role == "bacteria" else folder / "human_DNA"
            rna_folder = folder if role == "bacteria" else folder / "human_RNA"
            stem = "GCF_000000001.1_ASM1v1_genomic"
            gzwrite(dna_folder / f"{stem}.fna.gz", f">{role}\n{genome}\n")
            annotation = ["##gff-version 3"]
            for kind, start, stop, gene in (("CDS", 101, 1000, f"{role}_dna_gene"), ("rRNA", 1501, 2100, f"{role}_rrna"), ("tRNA", 2401, 2700, f"{role}_other"), ("CDS", 3001, 3300, f"{role}_splice"), ("CDS", 3801, 4100, f"{role}_splice")):
                annotation.append(f"{role}\tsynthetic\t{kind}\t{start}\t{stop}\t.\t+\t.\tlocus_tag={gene}")
            gzwrite(dna_folder / f"{stem}.gff.gz", "\n".join(annotation) + "\n")
            # Must not be included among genome sequences (duplicate ID by design).
            gzwrite(dna_folder / "ASM_cds_from_genomic.fna.gz", f">{role}\nAAAA\n")
            records = [transcript(f"{role}_coding_tx", seqs["coding"], "CDS", f"{role}_splice"), transcript(f"{role}_rrna_tx", seqs["rrna"], "rRNA", f"{role}_rrna"), transcript(f"{role}_other_tx", seqs["other"], "tRNA", f"{role}_other")]
            shards = [("ASM_rna_from_genomic", records)] if role == "bacteria" else [("human.1.rna", records[:2]), ("human.2.rna", records[2:])]
            for rna_stem, shard in shards:
                gzwrite(rna_folder / f"{rna_stem}.fna.gz", "".join(f">{r.id}\n{str(r.seq)}\n" for r in shard))
                with gzip.open(rna_folder / f"{rna_stem}.gbff.gz", "wt") as handle:
                    SeqIO.write(shard, handle, "genbank")
        cls.inputs.extend([("decoy", cls.decoy, 100), ("leftover", cls.leftover, 100)])
        folder = cls.references / "decoy"
        folder.mkdir(parents=True)
        fasta = folder / "Gina_Pangenome.fasta"
        fasta.write_text(f">decoy\n{cls.decoy}\n")
        with (cls.base / "decoy-build.log").open("w") as log:
            subprocess.run(["bowtie2-build", str(fasta), str(folder / "Gina_Pangenome")], stdout=log, stderr=subprocess.STDOUT, check=True)
        cls.original_files = cls.reference_snapshot()

    @classmethod
    def reference_snapshot(cls):
        return {p.relative_to(cls.references): (p.stat().st_size, p.stat().st_mtime_ns) for p in cls.references.rglob("*") if p.is_file()}

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def run_pipeline(self, name, layout="single", host=True, rerun=False, inputs=None):
        inputs = self.inputs if inputs is None else inputs
        directory = self.base / name
        reads = directory / "input reads"
        reads.mkdir(parents=True, exist_ok=True)
        for mate in ((1, 2) if layout == "paired" else (1,)):
            text = []
            for ident, sequence, start in inputs:
                offset = start if mate == 1 else start + 160
                bases = sequence[offset:offset + 100]
                self.assertEqual(len(bases), 100)
                if mate == 2:
                    bases = str(Seq(bases).reverse_complement())
                text.append(f"@{ident}\n{bases}\n+\n{'I' * 100}\n")
            filename = f"sample_R{mate}.fastq.gz" if layout == "paired" else "sample.fastq.gz"
            gzwrite(reads / filename, "".join(text))
        output = directory / "analysis output"
        config = directory / "run.env"
        settings = {"FASTQ_DIR": str(reads), "FASTQ_GLOB": "*.fastq.gz", "READ_LAYOUT": layout, "OUTPUT_DIR": str(output), "DECOY_REFERENCE_DIR": str(self.references / "decoy"), "BACTERIA_REFERENCE_DIR": str(self.references / "PAO1"), "HOST_REFERENCE_DIR": str(self.references / "human") if host else "", "REFERENCE_CACHE_DIR": str(self.cache), "THREADS": "2", "COUNT_STRAND": "forward", "GENE_POSITION_BINS": "5", "METAGENE_MIN_FEATURE_READS": "1", "SAMPLE_READS": "true", "SAMPLE_SIZE": str(len(inputs)), "SAMPLE_SEED": "17"}
        config.write_text("".join(f"{key}={shlex.quote(value)}\n" for key, value in settings.items()))
        log_path = directory / ("rerun.log" if rerun else "pipeline.log")
        env = dict(os.environ, MPLBACKEND="Agg", MPLCONFIGDIR=str(self.base / "matplotlib"))
        with log_path.open("w") as log:
            result = subprocess.run(["bash", str(REPOSITORY / "map_bacteria_with_decoys.sh"), str(config)], cwd=REPOSITORY, stdout=log, stderr=subprocess.STDOUT, env=env)
        if result.returncode:
            details="\n".join(f"{p.name}:\n{p.read_text()}" for p in output.rglob("*.bowtie*.txt"))
            self.fail(f"Pipeline exited {result.returncode}:\n{log_path.read_text()[-10000:]}\n{details}")
        with (output / "output_read_breakdown.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["Sample"], "sample")
        return output, rows[0]

    def assert_count(self, row, label, value):
        self.assertEqual(int(row[f"{label} (reads)"]), value, label)
        self.assertEqual(row[f"{label} (% of total reads)"], f"{100 * value / 12:.2f}%", label)

    def assert_outputs(self, output, row, host=True):
        for label, value in (("Total reads", 12), ("Reads passing cutadapt length filter", 12), ("Decoy aligned", 1), ("Bacteria aligned", 5 if host else 6), ("Bacteria CDS", 2), ("Bacteria rRNA", 1), ("Bacteria other RNA", 1), ("Bacteria intergenic/unannotated", 1 if host else 2), ("Host aligned", 5 if host else 0), ("Leftover reads", 1 if host else 5)):
            self.assert_count(row, label, value)
        for label, count in (("CDS", 2), ("rRNA", 1), ("other RNA", 1), ("intergenic/unannotated", 1)):
            self.assert_count(row, f"Host {label}", count if host else 0)
        m = json.loads((output / "run_manifest.json").read_text())
        stages = m["samples"][0]["stages"]
        names = ["decoy", "bacteria_rna", "host_rna", "bacteria_dna", "host_dna"] if host else ["decoy", "bacteria_rna", "bacteria_dna"]
        self.assertEqual([s["name"] for s in stages], names)
        for stage in stages:
            if stage["name"] == "decoy":
                continue
            data = json.loads(Path(stage["categories"]).read_text())
            self.assertEqual(sum(data["categories"].values()), data["aligned"])
            expected = 3 if stage["kind"] == "rna" else 2
            if not host and stage["name"] == "bacteria_dna":
                expected = 3
            self.assertEqual(data["aligned"], expected, stage["name"])
        for path in m["samples"][0]["leftover"]:
            with gzip.open(path, "rt") as handle:
                lines = handle.readlines()
            self.assertEqual(len(lines) // 4, 1 if host else 5)
        self.assertTrue((output / f"BACTERIA_{output.name}.csv").exists())
        self.assertEqual((output / "output_read_breakdown.csv").read_bytes(), (output / f"BACTERIA_{output.name}_read_breakdown.csv").read_bytes())
        self.assertTrue((output / "featurecounts_BACTERIA_summary.txt").exists())
        self.assertTrue(list(output.glob("bowtie_alignments/bacteria/rna/*.metagene_CDS_profile.png")))
        manifest = json.loads((output / "reference_manifest.json").read_text())
        self.assertEqual(manifest["decoy"]["index"], str(self.references / "decoy/Gina_Pangenome"))
        self.assertEqual(self.original_files, self.reference_snapshot(), "reference files must remain untouched")

    def exercise_array_workers(self, layout):
        # Execute the actual scheduler entry points locally; no Slurm service is
        # required to exercise worker isolation and the dependent final report.
        sys.path.insert(0, str(REPOSITORY))
        from slurm_workflow import create_plan
        directory = self.base / (f"array workers {layout}" + (" % parent" if layout == "paired" else ""))
        reads = directory / "input reads"
        reads.mkdir(parents=True)
        names = ["alpha 'first" + (" %" if layout == "paired" else ""), "zeta"]
        for name in names:
            for mate in ((1, 2) if layout == "paired" else (1,)):
                records = []
                for ident, sequence, start in self.inputs:
                    offset = start if mate == 1 else start + 160
                    bases = sequence[offset:offset + 100]
                    if mate == 2:
                        bases = str(Seq(bases).reverse_complement())
                    records.append(f"@{ident}\n{bases}\n+\n{'I' * 100}\n")
                filename = f"{name}_R{mate}.fastq.gz" if layout == "paired" else f"{name}.fastq.gz"
                gzwrite(reads / filename, "".join(records))
        output = directory / "project output"
        config = directory / "source config.env"
        settings = {
            "FASTQ_DIR": str(reads), "FASTQ_GLOB": "*.fastq.gz", "READ_LAYOUT": layout,
            "OUTPUT_DIR": str(output), "DECOY_REFERENCE_DIR": str(self.references / "decoy"),
            "BACTERIA_REFERENCE_DIR": str(self.references / "PAO1"),
            "HOST_REFERENCE_DIR": str(self.references / "human"),
            "REFERENCE_CACHE_DIR": str(self.cache), "THREADS": "2", "COUNT_STRAND": "forward",
            "GENE_POSITION_BINS": "5", "METAGENE_MIN_FEATURE_READS": "1",
            "SAMPLE_READS": "false", "CONDA_ENV": "",
        }
        config.write_text("".join(f"{k}={shlex.quote(v)}\n" for k, v in settings.items()))
        plan_path, plan = create_plan(config, REPOSITORY, settings)
        self.assertEqual(plan["array_concurrency"], 8)
        self.assertEqual(plan["max_samples"], 40)
        self.assertEqual([sample["sample"] for sample in plan["samples"]], names)
        # Workers must use frozen settings even if the user's source file changes.
        config.write_text("FASTQ_DIR=/missing/after/submission\nOUTPUT_DIR=/missing/after/submission\n")

        def stage(mode, index=None, expect_success=True):
            log_path = directory / f"{mode}_{index}.log"
            env = dict(os.environ, MPLBACKEND="Agg", MPLCONFIGDIR=str(directory / f"matplotlib-{index}"), SLURM_CPUS_PER_TASK="2")
            if index is not None:
                env["SLURM_ARRAY_TASK_ID"] = str(index)
            with log_path.open("w") as log:
                result = subprocess.run([sys.executable, str(REPOSITORY / "slurm_workflow.py"), "execute", mode, str(plan_path)], cwd=REPOSITORY, stdout=log, stderr=subprocess.STDOUT, env=env)
            if expect_success:
                self.assertEqual(result.returncode, 0, log_path.read_text()[-12000:])
            else:
                self.assertNotEqual(result.returncode, 0)
            return result

        stage("prepare")
        reference_manifest = Path(plan["reference_manifest"])
        self.assertTrue(reference_manifest.is_file())
        index_snapshot = {p: p.stat().st_mtime_ns for p in self.cache.rglob("*.bt2")}
        with ThreadPoolExecutor(max_workers=2) as executor:
            list(executor.map(lambda index: stage("sample", index), range(2)))
        self.assertEqual(index_snapshot, {p: p.stat().st_mtime_ns for p in self.cache.rglob("*.bt2")})
        self.assertFalse((output / "output_read_breakdown.csv").exists())
        for index, sample in enumerate(plan["samples"]):
            worker_output = Path(sample["output_dir"])
            manifest = json.loads((worker_output / "run_manifest.json").read_text())
            self.assertEqual(manifest["sample_index"], index)
            self.assertEqual(manifest["input_sample_manifest"], plan["sample_manifest"])
            self.assertEqual([entry["sample"] for entry in manifest["samples"]], [sample["sample"]])
            self.assertFalse((worker_output / "output_read_breakdown.csv").exists())
            self.assertFalse((worker_output / f"BACTERIA_{worker_output.name}.csv").exists())
        stage("finalize")
        with (output / "output_read_breakdown.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual([row["Sample"] for row in rows], names)
        for row in rows:
            for label, count in (("Total reads", 12), ("Reads passing cutadapt length filter", 12), ("Decoy aligned", 1), ("Bacteria aligned", 5), ("Host aligned", 5), ("Leftover reads", 1), ("Bacteria CDS", 2), ("Bacteria rRNA", 1), ("Bacteria other RNA", 1), ("Bacteria intergenic/unannotated", 1)):
                self.assert_count(row, label, count)
        with (output / f"BACTERIA_{output.name}.csv").open() as handle:
            self.assertEqual([row["Sample"] for row in csv.DictReader(handle)], names)
        combined = json.loads((output / "run_manifest.json").read_text())
        self.assertEqual([entry["sample"] for entry in combined["samples"]], names)
        # Every stage's gene matrix keeps both sample columns, including the
        # historical CSV companion and featureCounts assignment summary.
        merged_files = []
        for relative, organism in ((Path("."), "BACTERIA"), (Path("bowtie_alignments/bacteria/rna"), "BACTERIA"), (Path("bowtie_alignments/host"), "HOST"), (Path("bowtie_alignments/host/rna"), "HOST")):
            for feature in ("summary", "rRNA"):
                table = output / relative / f"featurecounts_{organism}_{feature}.txt"
                with table.open() as handle:
                    matrix = list(csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t"))
                self.assertEqual(len(matrix[0]), 8)
                self.assertTrue(all(len(row) == 8 for row in matrix))
                expected_sams = [f"{organism}_{sample}_22bp.sam" for sample in names]
                self.assertEqual([Path(column).name for column in matrix[0][6:]], expected_sams)
                counts = [sum(int(row[column]) for row in matrix[1:]) for column in (6, 7)]
                expected = 1 if feature == "summary" or relative.name == "rna" else 0
                self.assertEqual(counts, [expected, expected])
                summary = Path(str(table) + ".summary")
                with summary.open() as handle:
                    assignment_rows = list(csv.reader(handle, delimiter="\t"))
                self.assertEqual(assignment_rows[0], ["Status", *matrix[0][6:]])
                self.assertTrue(all(len(row) == 3 for row in assignment_rows))
                merged_files.extend([table, summary])
                if feature == "summary":
                    csv_table = table.with_suffix(".csv")
                    with csv_table.open() as handle:
                        csv_rows = list(csv.reader(line for line in handle if not line.startswith("#")))
                    self.assertEqual(csv_rows, matrix)
                    merged_files.append(csv_table)
        # Missing or stale worker output must not replace successful reports
        # or any previously published gene matrices/assignment summaries.
        successful_outputs = {path: path.read_bytes() for path in merged_files + [output / "output_read_breakdown.csv"]}

        def assert_outputs_preserved():
            self.assertEqual({path: path.read_bytes() for path in successful_outputs}, successful_outputs)
        worker_path = Path(plan["samples"][1]["output_dir"]) / "run_manifest.json"
        saved = worker_path.read_bytes()
        try:
            worker_path.unlink()
            stage("finalize", expect_success=False)
            assert_outputs_preserved()
            stale = json.loads(saved)
            stale["sample_index"] = 0
            worker_path.write_text(json.dumps(stale))
            stage("finalize", expect_success=False)
            assert_outputs_preserved()
        finally:
            worker_path.write_bytes(saved)
        missing_table = Path(plan["samples"][1]["output_dir"]) / "featurecounts_BACTERIA_summary.txt"
        original_table = missing_table.read_bytes()
        try:
            missing_table.unlink()
            stage("finalize", expect_success=False)
            assert_outputs_preserved()
        finally:
            missing_table.write_bytes(original_table)
        # A required stage table absent from every worker is still a failure;
        # it must not be mistaken for an optional annotation type.
        required_tables = [Path(sample["output_dir"]) / "featurecounts_BACTERIA_summary.txt" for sample in plan["samples"]]
        saved_tables = {path: path.read_bytes() for path in required_tables}
        try:
            for path in required_tables:
                path.unlink()
            stage("finalize", expect_success=False)
            assert_outputs_preserved()
        finally:
            for path, contents in saved_tables.items():
                path.write_bytes(contents)
        stage("sample", index=2, expect_success=False)
        self.assertEqual(self.original_files, self.reference_snapshot())

    def test_two_single_end_array_workers_and_combined_reports(self):
        self.exercise_array_workers("single")

    def test_two_paired_end_array_workers_and_combined_reports(self):
        self.exercise_array_workers("paired")

    def test_single_end_host_rna_priority_and_cache_reuse(self):
        output, row = self.run_pipeline("single host")
        self.assert_outputs(output, row)
        indices = {p: p.stat().st_mtime_ns for p in self.cache.rglob("*.bt2")}
        self.assertTrue(indices)
        output, row = self.run_pipeline("single host", rerun=True)
        self.assert_outputs(output, row)
        self.assertEqual(indices, {p: p.stat().st_mtime_ns for p in self.cache.rglob("*.bt2")})

    def test_paired_end_counts_fragments_once(self):
        output, row = self.run_pipeline("paired host", layout="paired")
        self.assert_outputs(output, row)

    def test_all_reads_consumed_by_rna_keep_empty_dna_stages_valid(self):
        inputs = [entry for entry in self.inputs if entry[0] in {"bacteria_rna_rRNA", "bacteria_rna_other"}]
        output, row = self.run_pipeline("empty residual stages", inputs=inputs)
        self.assertEqual(int(row["Total reads (reads)"]), 2)
        self.assertEqual(int(row["Bacteria aligned (reads)"]), 2)
        self.assertEqual(row["Bacteria aligned (% of total reads)"], "100.00%")
        self.assertEqual(int(row["Leftover reads (reads)"]), 0)
        manifest = json.loads((output / "run_manifest.json").read_text())
        for stage in manifest["samples"][0]["stages"]:
            if stage["organism"] == "decoy":
                continue
            expected = 2 if stage["name"] == "bacteria_rna" else 0
            self.assertEqual(json.loads(Path(stage["categories"]).read_text())["aligned"], expected)
        for path in manifest["samples"][0]["leftover"]:
            with gzip.open(path, "rt") as handle:
                self.assertEqual(handle.read(), "")

    def test_no_host_classifies_shared_dna_and_keeps_host_reads(self):
        output, row = self.run_pipeline("single no host", host=False)
        self.assert_outputs(output, row, host=False)


if __name__ == "__main__":
    unittest.main()
