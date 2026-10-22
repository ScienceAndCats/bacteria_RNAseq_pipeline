#!/usr/bin/env python3
"""Write the project summary, including bacterial RNA and genomic alignments."""

import argparse
import csv
import json
import re
from pathlib import Path

from read_breakdown_csv import build_table, read_manifest, report_path


HEADER = (
    "Sample", "Raw Reads", "Reads Written", "Overall Allignment", "Percent Covered",
    "Mean Coverage", "Reads Assigned to features", "Features with non-zero reads",
)


def coverage_metrics(path):
    """Return length-weighted genomic coverage from samtools coverage output."""
    text = path.read_text()
    if text.startswith("#rname\t"):
        rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
        lengths = [int(row["endpos"]) - int(row["startpos"]) + 1 for row in rows]
        length = sum(lengths)
        if not length:
            return "0.00%", "0.0000"
        covered = sum(int(row["covbases"]) for row in rows)
        depth = sum(float(row["meandepth"]) * size for row, size in zip(rows, lengths))
        return f"{covered * 100 / length:.2f}%", f"{depth / length:.4f}"
    # Retain compatibility with the earlier human-readable samtools coverage -m report.
    covered = re.search(r"Percent covered:\s*([\d.]+%)", text)
    depth = re.search(r"Mean coverage:\s*([\d.]+)x", text)
    if not covered or not depth:
        raise ValueError(f"Could not parse samtools coverage report: {path}")
    return covered.group(1), depth.group(1)


def legacy_feature_counts(output_dir):
    path = output_dir / "featurecounts_BACTERIA_summary.txt"
    with path.open() as handle:
        rows = list(csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t"))
    if not rows or len(rows[0]) < 7:
        raise ValueError(f"Unexpected featureCounts format: {path}")
    result = {}
    for column, sam in enumerate(rows[0][6:], 6):
        match = re.fullmatch(r"BACTERIA_(.+)_\d+bp\.sam", Path(sam).name)
        if not match:
            raise ValueError(f"Unexpected alignment name: {sam}")
        values = [int(row[column]) for row in rows[1:]]
        result[match.group(1)] = sum(values), sum(value > 0 for value in values)
    return result


def summary_rows(output_dir):
    table = build_table(output_dir)
    manifest = {item["sample"]: item for item in read_manifest(output_dir)} if (output_dir / "run_manifest.json").exists() else None
    legacy = legacy_feature_counts(output_dir) if manifest is None else None
    rows = []
    for sample, counts, _ in table:
        coverage = ("", "")
        if manifest is not None:
            bacterial = [stage for stage in manifest[sample]["stages"] if stage["organism"] == "bacteria"]
            genomic = [stage for stage in bacterial if stage["kind"] == "dna"]
            if len(genomic) > 1:
                raise ValueError(f"Multiple bacterial genomic coverage stages for {sample}")
            if genomic and genomic[0].get("coverage"):
                coverage = coverage_metrics(report_path(output_dir, genomic[0]["coverage"]))
            cds_reads = counts["Bacteria CDS"]
            features = set()
            have_features = True
            for stage in bacterial:
                data = json.loads(report_path(output_dir, stage["categories"]).read_text())
                identifiers = data.get("cds_features_with_reads")
                if identifiers is None:
                    have_features = False
                elif not isinstance(identifiers, list) or not all(isinstance(value, str) for value in identifiers):
                    raise ValueError(f"Invalid CDS feature identifiers for {sample}")
                else:
                    features.update(identifiers)
            genes = len(features) if have_features else ""
        else:
            cds_reads, genes = legacy[sample]
            logs = [path for path in output_dir.glob("*cutadapt_log.txt") if re.fullmatch(re.escape(sample) + r"_\d+bp\.cutadapt_log\.txt", path.name)]
            if len(logs) != 1:
                raise ValueError(f"Ambiguous cutadapt log for {sample}")
            length = re.search(r"_(\d+)bp\.cutadapt_log", logs[0].name).group(1)
            coverage = coverage_metrics(output_dir / f"bowtie_alignments/bacteria/BACTERIA_{sample}_{length}bp_coverage.txt")
        rows.append((sample, counts["Total reads"], counts["Reads passing cutadapt length filter"],
                     counts["Bacteria aligned"], *coverage, cds_reads, genes))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    rows = summary_rows(output_dir)
    destination = args.output or output_dir / f"BACTERIA_{output_dir.name}.csv"
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(HEADER)
        writer.writerows(rows)


if __name__ == "__main__":
    main()
