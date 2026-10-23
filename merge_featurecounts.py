#!/usr/bin/env python3
"""Merge isolated sample featureCounts outputs without loading gene matrices.

All source tables must contain the same ordered annotation rows and exactly one
sample column. Matrix and summary columns follow the supplied directory order.
The destination should be a finalizer staging directory, published only after
all project reports have been validated.
"""

import argparse
import csv
from contextlib import ExitStack
from itertools import zip_longest
from pathlib import Path
import re
import tempfile


ANNOTATION_COLUMNS = ("Geneid", "Chr", "Start", "End", "Strand", "Length")
STAGES = (
    (Path("."), "BACTERIA"),
    (Path("bowtie_alignments/bacteria/rna"), "BACTERIA"),
    (Path("bowtie_alignments/host"), "HOST"),
    (Path("bowtie_alignments/host/rna"), "HOST"),
)
TABLE_PATHS = tuple(
    directory / f"featurecounts_{organism}_{feature}.txt"
    for directory, organism in STAGES
    for feature in ("summary", "rRNA")
)
MERGE_COMMENT = "# Merged featureCounts tables; sample columns follow input directory order.\n"


def _rows(handle):
    return csv.reader((line for line in handle if not line.startswith("#") and line.strip()), delimiter="\t")


def _count(value, path, row_number):
    if re.fullmatch(r"\d+", value) is None:
        raise ValueError(f"Invalid featureCounts count in {path}, row {row_number}: {value!r}")
    return value


def _merge_matrix(sources, destination):
    """Stream corresponding rows from one open file per sample."""
    with ExitStack() as stack:
        readers = [_rows(stack.enter_context(path.open(newline=""))) for path in sources]
        headers = [next(reader, None) for reader in readers]
        for path, header in zip(sources, headers):
            if header is None or len(header) != 7 or tuple(header[:6]) != ANNOTATION_COLUMNS or not header[6]:
                raise ValueError(f"Expected six annotation columns and one sample column in {path}")
        sample_columns = [header[6] for header in headers]
        if len(sample_columns) != len(set(sample_columns)):
            raise ValueError(f"Duplicate featureCounts sample columns in {sources}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = stack.enter_context(tempfile.NamedTemporaryFile(mode="w", newline="", prefix=f".{destination.name}.", dir=destination.parent, delete=False))
        temporary_path = Path(temporary.name)
        try:
            temporary.write(MERGE_COMMENT)
            writer = csv.writer(temporary, delimiter="\t", lineterminator="\n")
            writer.writerow((*ANNOTATION_COLUMNS, *sample_columns))
            for number, rows in enumerate(zip_longest(*readers), 2):
                if any(row is None for row in rows):
                    raise ValueError(f"Different gene row counts in featureCounts tables: {sources}")
                for path, row in zip(sources, rows):
                    if len(row) != 7:
                        raise ValueError(f"Expected one sample count in {path}, row {number}")
                    _count(row[6], path, number)
                annotation = rows[0][:6]
                if any(row[:6] != annotation for row in rows[1:]):
                    raise ValueError(f"Mismatched featureCounts annotation at row {number}: {sources}")
                writer.writerow((*annotation, *(row[6] for row in rows)))
            temporary.flush()
            temporary.close()
            temporary_path.replace(destination)
        except BaseException:
            temporary_path.unlink(missing_ok=True)
            raise
    return sample_columns


def _merge_summary(sources, destination, expected_columns):
    with ExitStack() as stack:
        readers = [_rows(stack.enter_context(path.open(newline=""))) for path in sources]
        for path, reader, sample in zip(sources, readers, expected_columns):
            header = next(reader, None)
            if header != ["Status", sample]:
                raise ValueError(f"featureCounts summary sample column disagrees with its table: {path}")
        temporary = stack.enter_context(tempfile.NamedTemporaryFile(mode="w", newline="", prefix=f".{destination.name}.", dir=destination.parent, delete=False))
        temporary_path = Path(temporary.name)
        try:
            writer = csv.writer(temporary, delimiter="\t", lineterminator="\n")
            writer.writerow(("Status", *expected_columns))
            seen_statuses = set()
            for number, rows in enumerate(zip_longest(*readers), 2):
                if any(row is None for row in rows):
                    raise ValueError(f"Different status row counts in featureCounts summaries: {sources}")
                for path, row in zip(sources, rows):
                    if len(row) != 2 or not row[0]:
                        raise ValueError(f"Invalid featureCounts status row in {path}, row {number}")
                    _count(row[1], path, number)
                status = rows[0][0]
                if any(row[0] != status for row in rows[1:]) or status in seen_statuses:
                    raise ValueError(f"Mismatched or duplicate featureCounts status at row {number}: {sources}")
                seen_statuses.add(status)
                writer.writerow((status, *(row[1] for row in rows)))
            if "Assigned" not in seen_statuses:
                raise ValueError(f"No Assigned status in featureCounts summaries: {sources}")
            temporary.flush()
            temporary.close()
            temporary_path.replace(destination)
        except BaseException:
            temporary_path.unlink(missing_ok=True)
            raise


def _write_csv(source, destination):
    with source.open(newline="") as handle, destination.open("w", newline="") as output:
        output.write(MERGE_COMMENT)
        writer = csv.writer(output, lineterminator="\n")
        writer.writerows(_rows(handle))


def merge_featurecounts(sample_directories, destination):
    """Return generated files after merging all available stage/type matrices.

    Stage/type tables and their optional .summary files must each be present for
    every sample or absent for every sample. Sources are never modified. CDS
    matrices also receive the historical .csv companion; rRNA tables retain TSV.
    """
    directories = [Path(directory).resolve() for directory in sample_directories]
    destination = Path(destination).resolve()
    if not directories:
        raise ValueError("At least one sample output directory is required")
    if len(directories) != len(set(directories)):
        raise ValueError("Duplicate sample output directories")
    for directory in directories:
        if not directory.is_dir():
            raise ValueError(f"Sample output directory does not exist: {directory}")
        if destination == directory:
            raise ValueError(f"Destination must differ from sample output directory: {directory}")

    plans = []
    for relative in TABLE_PATHS:
        sources = [directory / relative for directory in directories]
        present = [path.is_file() for path in sources]
        summary_sources = [Path(str(path) + ".summary") for path in sources]
        summary_present = [path.is_file() for path in summary_sources]
        if any(present) and not all(present):
            missing = [str(path) for path, exists in zip(sources, present) if not exists]
            raise ValueError(f"Missing featureCounts tables for some samples: {', '.join(missing)}")
        if any(summary_present) and (not all(summary_present) or not all(present)):
            raise ValueError(f"Missing featureCounts tables or summary companions for stage: {relative}")
        if all(present):
            plans.append((relative, sources, summary_sources if all(summary_present) else None))

    generated = []
    for relative, sources, summary_sources in plans:
        target = destination / relative
        columns = _merge_matrix(sources, target)
        generated.append(target)
        if summary_sources:
            summary_target = Path(str(target) + ".summary")
            _merge_summary(summary_sources, summary_target, columns)
            generated.append(summary_target)
        if relative.name.endswith("_summary.txt"):
            csv_target = target.with_suffix(".csv")
            _write_csv(target, csv_target)
            generated.append(csv_target)
    return generated


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-dir", action="append", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        merge_featurecounts(args.sample_dir, args.output_dir)
    except (ValueError, OSError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
