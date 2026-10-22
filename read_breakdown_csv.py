#!/usr/bin/env python3
"""Create a per-sample read disposition table from pipeline reports."""

import argparse
import csv
import json
import re
from pathlib import Path


ATTRIBUTES = (
    "Total reads",
    "Reads passing cutadapt length filter",
    "Decoy aligned",
    "Bacteria aligned",
    "Bacteria non-rRNA",
    "Bacteria rRNA",
    "Bacteria CDS",
    "Bacteria other RNA",
    "Bacteria intergenic/unannotated",
    "Host aligned",
    "Host non-rRNA",
    "Host rRNA",
    "Host CDS",
    "Host other RNA",
    "Host intergenic/unannotated",
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
    # Bowtie2 omits the detailed hierarchy when its input FASTQ is empty.
    # Its CPU-dispatch wrapper may prefix warning lines even for valid runs.
    count_report = "\n".join(line for line in report.splitlines() if not line.lstrip().startswith("[WARNING]"))
    if re.fullmatch(r"\s*0 reads\s*\n0\.00% overall alignment rate\s*", count_report):
        return 0, 0
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


def build_legacy_table(output_dir):
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
        for organism in ("Bacteria", "Host"):
            for category in ("CDS", "other RNA", "intergenic/unannotated"):
                attribute = f"{organism} {category}"
                counts[attribute] = None
                sources[attribute] = "Unavailable in legacy reports; rerun the folder-reference pipeline"
        rows.append((sample, counts, sources))
    return rows


CATEGORY_LABELS = {
    "rRNA": "rRNA",
    "CDS": "CDS",
    "other_RNA": "other RNA",
    "intergenic_unannotated": "intergenic/unannotated",
}


def report_path(output_dir, value):
    path = Path(value)
    return path if path.is_absolute() else output_dir / path


def source_path(output_dir, path):
    try:
        return str(path.relative_to(output_dir))
    except ValueError:
        return str(path)


def read_manifest(output_dir):
    with (output_dir / "run_manifest.json").open() as handle:
        manifest = json.load(handle)
    samples = manifest.get("samples")
    if not isinstance(samples, list) or not samples:
        raise ValueError("run_manifest.json has no samples")
    names = [sample["sample"] for sample in samples]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate sample names in run_manifest.json")
    return samples


def nonnegative_count(value, context):
    if type(value) is not int or value < 0:
        raise ValueError(f"Invalid nonnegative read count for {context}: {value!r}")
    return value


def build_manifest_table(output_dir):
    rows = []
    for sample_info in read_manifest(output_dir):
        sample = sample_info["sample"]
        cutadapt_path = report_path(output_dir, sample_info["cutadapt_log"])
        total, passing = cutadapt_counts(cutadapt_path)
        if passing > total:
            raise ValueError(f"Passing reads exceed input for {sample}")
        counts = {attribute: 0 for attribute in ATTRIBUTES}
        sources = {attribute: "Host alignment not configured" for attribute in ATTRIBUTES}
        counts["Total reads"], counts["Reads passing cutadapt length filter"] = total, passing
        cutadapt_source = source_path(output_dir, cutadapt_path)
        sources["Total reads"] = f"{cutadapt_source}: Total reads/read pairs processed"
        sources["Reads passing cutadapt length filter"] = f"{cutadapt_source}: Reads/pairs written (passing filters)"
        remaining = passing
        logs = {"bacteria": [], "host": []}
        categories = {"bacteria": [], "host": []}
        stages = sample_info.get("stages", [])
        if not stages or stages[0].get("organism") != "decoy":
            raise ValueError(f"Missing initial decoy stage for {sample}")
        stage_names = set()
        saw_dna = False
        for stage in stages:
            name, organism, kind = stage["name"], stage["organism"], stage["kind"]
            if name in stage_names or organism not in ("decoy", "bacteria", "host"):
                raise ValueError(f"Invalid or duplicate stage for {sample}: {name}")
            stage_names.add(name)
            if kind == "rna" and saw_dna:
                raise ValueError(f"RNA stage follows a DNA stage for {sample}: {name}")
            if kind == "dna":
                saw_dna = True
            log_path = report_path(output_dir, stage["log"])
            stage_input, aligned = bowtie_counts(log_path)
            if stage_input != remaining or not 0 <= aligned <= stage_input:
                raise ValueError(f"Stage totals are inconsistent for sample {sample}, stage {name}")
            remaining -= aligned
            log_source = source_path(output_dir, log_path)
            if organism == "decoy":
                if kind != "decoy" or stage is not stages[0]:
                    raise ValueError(f"Unexpected decoy stage for {sample}: {name}")
                counts["Decoy aligned"] = aligned
                sources["Decoy aligned"] = f"{log_source}: input minus aligned 0 times"
                continue
            if kind not in ("rna", "dna"):
                raise ValueError(f"Unknown reference kind in stage {name}: {kind}")
            category_path = report_path(output_dir, stage["categories"])
            data = json.loads(category_path.read_text())
            category_counts = data.get("categories", {})
            if set(category_counts) != set(CATEGORY_LABELS):
                raise ValueError(f"Missing or unknown categories in {category_path}")
            classified_aligned = nonnegative_count(data.get("aligned"), category_path)
            for category, value in category_counts.items():
                nonnegative_count(value, f"{category_path}: {category}")
            if classified_aligned != aligned or sum(category_counts.values()) != aligned:
                raise ValueError(f"Category totals do not match aligned reads in {category_path}")
            label = organism.title()
            counts[f"{label} aligned"] += aligned
            for category, category_label in CATEGORY_LABELS.items():
                counts[f"{label} {category_label}"] += category_counts[category]
            logs[organism].append(log_source)
            categories[organism].append(source_path(output_dir, category_path))
        if not logs["bacteria"]:
            raise ValueError(f"No bacterial alignment stages for {sample}")
        for organism in ("bacteria", "host"):
            label = organism.title()
            counts[f"{label} non-rRNA"] = counts[f"{label} aligned"] - counts[f"{label} rRNA"]
            if logs[organism]:
                sources[f"{label} aligned"] = "; ".join(logs[organism]) + ": RNA + DNA aligned reads/read pairs"
                for category_label in CATEGORY_LABELS.values():
                    sources[f"{label} {category_label}"] = "; ".join(categories[organism]) + ": mutually exclusive annotation category"
                sources[f"{label} non-rRNA"] = f"Calculated: {label} aligned minus {label} rRNA"
        counts["Leftover reads"] = remaining
        sources["Leftover reads"] = "Calculated: passing filter minus decoy, bacteria, and host aligned"
        rows.append((sample, counts, sources))
    return rows


def build_table(output_dir):
    output_dir = Path(output_dir).resolve()
    if (output_dir / "run_manifest.json").exists():
        return build_manifest_table(output_dir)
    return build_legacy_table(output_dir)


def write_table(destination, rows):
    header = ["Sample"]
    for attribute in ATTRIBUTES:
        header.extend((f"{attribute} (reads)", f"{attribute} (% of total reads)", f"{attribute} source"))
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for sample, counts, sources in rows:
            row = [sample]
            for attribute in ATTRIBUTES:
                count = counts[attribute]
                row.extend((count, percentage(count, counts["Total reads"]) if count is not None else "", sources[attribute]))
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    destinations = [args.output] if args.output else [
        output_dir / "output_read_breakdown.csv",
        output_dir / f"BACTERIA_{output_dir.name}_read_breakdown.csv",
    ]
    rows = build_table(output_dir)
    for destination in destinations:
        write_table(destination, rows)


if __name__ == "__main__":
    main()
