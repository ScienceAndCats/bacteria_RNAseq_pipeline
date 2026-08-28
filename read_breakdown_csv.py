#!/usr/bin/env python3
"""Create a per-sample read disposition table from pipeline reports."""

import argparse
import csv
import re
from pathlib import Path


ATTRIBUTES = (
    "Total reads",
    "Reads passing cutadapt length filter",
    "Decoy aligned",
    "Bacteria aligned",
    "Bacteria non-rRNA",
    "Bacteria rRNA",
    "Host aligned",
    "Host non-rRNA",
    "Host rRNA",
    "Leftover reads",
)


def report_number(report, pattern, path):
    match = re.search(pattern, report, re.MULTILINE)
    if not match:
        raise ValueError(f"Could not find the required count in {path}")
    return int(match.group(1).replace(",", ""))


def cutadapt_counts(path):
    report = path.read_text()
    total = report_number(
        report, r"Total (?:read pairs|reads) processed:\s*([\d,]+)", path
    )
    passing = report_number(
        report, r"(?:Pairs|Reads) written \(passing filters\):\s*([\d,]+)", path
    )
    return total, passing


def bowtie_counts(path):
    """Return input and aligned read/read-pair counts from a Bowtie2 report."""
    report = path.read_text()
    total = report_number(report, r"^\s*([\d,]+) reads; of these:", path)
    if "were paired; of these:" in report:
        unaligned = report_number(
            report, r"^\s*([\d,]+) \([^\n]+\) aligned concordantly 0 times$", path
        )
    else:
        unaligned = report_number(
            report, r"^\s*([\d,]+) \([^\n]+\) aligned 0 times$", path
        )
    return total, total - unaligned


def sample_from_log(path):
    match = re.fullmatch(r"(.+)_\d+bp\.cutadapt_log\.txt", path.name)
    if not match:
        raise ValueError(f"Unexpected cutadapt log name: {path}")
    return match.group(1)


def assigned_by_sample(path, prefix):
    """Read the Assigned row of a featureCounts .summary file."""
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or rows[0][0] != "Status":
        raise ValueError(f"Unexpected featureCounts summary format: {path}")
    assigned = next((row for row in rows[1:] if row[0] == "Assigned"), None)
    if assigned is None:
        raise ValueError(f"No Assigned row in {path}")
    result = {}
    for sam, count in zip(rows[0][1:], assigned[1:]):
        name = Path(sam).name
        match = re.fullmatch(rf"{prefix}_(.+)_\d+bp\.sam", name)
        if not match:
            raise ValueError(f"Unexpected alignment name in {path}: {sam}")
        result[match.group(1)] = int(count)
    return result


def percentage(count, total):
    return f"{count * 100 / total:.2f}%" if total else "0.00%"


def build_table(output_dir):
    cutadapt_logs = sorted(output_dir.glob("*cutadapt_log.txt"), key=lambda p: p.name.casefold())
    if not cutadapt_logs:
        raise ValueError(f"No cutadapt logs found in {output_dir}")

    bacteria_rrna_path = output_dir / "featurecounts_BACTERIA_rRNA.txt.summary"
    bacteria_rrna = assigned_by_sample(bacteria_rrna_path, "BACTERIA")
    host_rrna_path = output_dir / "bowtie_alignments/host/featurecounts_HOST_rRNA.txt.summary"
    host_enabled = host_rrna_path.exists()
    host_rrna = assigned_by_sample(host_rrna_path, "HOST") if host_enabled else {}

    rows = []
    for cutadapt_path in cutadapt_logs:
        sample = sample_from_log(cutadapt_path)
        total, passing = cutadapt_counts(cutadapt_path)
        # Allow a configured minimum length other than the historical 22 bp default.
        length = re.search(r"_(\d+)bp\.cutadapt_log", cutadapt_path.name).group(1)
        decoy_path = output_dir / f"bowtie_alignments/decoy/{sample}_{length}bp.mapped_to_other_bugs.bowtie2.txt"
        bacteria_path = output_dir / f"bowtie_alignments/bacteria/BACTERIA_{sample}_{length}bp.bowtie_output.txt"
        decoy_input, decoy_aligned = bowtie_counts(decoy_path)
        bacteria_input, bacteria_aligned = bowtie_counts(bacteria_path)
        if decoy_input != passing or bacteria_input != passing - decoy_aligned:
            raise ValueError(f"Stage totals are inconsistent for sample {sample}")

        host_path = output_dir / f"bowtie_alignments/host/HOST_{sample}_{length}bp.bowtie_output.txt"
        if host_enabled:
            host_input, host_aligned = bowtie_counts(host_path)
            if host_input != bacteria_input - bacteria_aligned:
                raise ValueError(f"Host stage total is inconsistent for sample {sample}")
        else:
            host_aligned = 0

        bacterial_rrna = bacteria_rrna[sample]
        host_rrna_count = host_rrna.get(sample, 0)
        if bacterial_rrna > bacteria_aligned or host_rrna_count > host_aligned:
            raise ValueError(f"rRNA count exceeds aligned count for sample {sample}")
        counts = {
            "Total reads": total,
            "Reads passing cutadapt length filter": passing,
            "Decoy aligned": decoy_aligned,
            "Bacteria aligned": bacteria_aligned,
            "Bacteria non-rRNA": bacteria_aligned - bacterial_rrna,
            "Bacteria rRNA": bacterial_rrna,
            "Host aligned": host_aligned,
            "Host non-rRNA": host_aligned - host_rrna_count,
            "Host rRNA": host_rrna_count,
            "Leftover reads": passing - decoy_aligned - bacteria_aligned - host_aligned,
        }
        sources = {
            "Total reads": f"{cutadapt_path.name}: Total reads/read pairs processed",
            "Reads passing cutadapt length filter": f"{cutadapt_path.name}: Reads/pairs written (passing filters)",
            "Decoy aligned": f"{decoy_path.relative_to(output_dir)}: input minus aligned 0 times",
            "Bacteria aligned": f"{bacteria_path.relative_to(output_dir)}: input minus aligned 0 times",
            "Bacteria rRNA": f"{bacteria_rrna_path.relative_to(output_dir)}: Assigned",
            "Host aligned": (f"{host_path.relative_to(output_dir)}: input minus aligned 0 times" if host_enabled else "Host alignment not configured"),
            "Host rRNA": (f"{host_rrna_path.relative_to(output_dir)}: Assigned" if host_enabled else "Host alignment not configured"),
        }
        sources["Bacteria non-rRNA"] = "Calculated: Bacteria aligned minus Bacteria rRNA"
        sources["Host non-rRNA"] = "Calculated: Host aligned minus Host rRNA"
        sources["Leftover reads"] = "Calculated: passing filter minus decoy, bacteria, and host aligned"
        rows.append((sample, counts, sources))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    destination = args.output or output_dir / f"BACTERIA_{output_dir.name}_read_breakdown.csv"
    header = ["Sample"]
    for attribute in ATTRIBUTES:
        header.extend((f"{attribute} (reads)", f"{attribute} (% of total reads)", f"{attribute} source"))
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for sample, counts, sources in build_table(output_dir):
            row = [sample]
            for attribute in ATTRIBUTES:
                row.extend((counts[attribute], percentage(counts[attribute], counts["Total reads"]), sources[attribute]))
            writer.writerow(row)


if __name__ == "__main__":
    main()
