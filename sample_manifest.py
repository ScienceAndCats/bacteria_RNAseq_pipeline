#!/usr/bin/env python3
"""Discover FASTQ samples once for serial runs or isolated array workers."""

import argparse
import glob
import json
from pathlib import Path
import re
import unicodedata

PAIRED_NAME = re.compile(r"^(?P<prefix>.+)_R(?P<mate>[12])(?P<suffix>_.*)?\.fastq(?:\.gz)?$")


def validate_text(value, label):
    """Reject control characters before a filename enters manifests or scripts."""
    if not value or any(unicodedata.category(character) == "Cc" for character in value):
        raise ValueError(f"{label} is empty or contains a control character")
    return value


def validate_sample_name(value):
    validate_text(value, "Sample name")
    if value in (".", "..") or "/" in value:
        raise ValueError(f"Invalid sample name: {value!r}")
    return value


def is_fastq(path):
    return path.name.endswith((".fastq", ".fastq.gz"))


def discover_samples(directory, pattern="*.fastq*", layout="single"):
    """Return deterministic sample IDs and absolute mate paths.

    Non-FASTQ glob matches are ignored. In paired mode selected FASTQs must
    have an R1/R2 marker, and their matching mate must exist even if the glob
    selected only R1 or R2. Mates can independently be plain or gzip files.
    Each suffix after the read marker must agree (for example ``_R1_001``).
    """
    if layout not in ("single", "paired"):
        raise ValueError("READ_LAYOUT must be single or paired")
    root = Path(directory).expanduser().resolve()
    validate_text(str(root), "FASTQ directory")
    if not root.is_dir():
        raise ValueError(f"FASTQ directory does not exist: {root}")
    validate_text(pattern, "FASTQ glob")
    if Path(pattern).is_absolute() or ".." in Path(pattern).parts:
        raise ValueError("FASTQ glob must stay within FASTQ_DIR")
    selected = [Path(value) for value in sorted(glob.glob(glob.escape(str(root)) + "/" + pattern, recursive=True))
                if Path(value).is_file() and is_fastq(Path(value))]
    if not selected:
        raise ValueError(f"No FASTQ files found in {root} matching {pattern}")
    samples = {}
    selected_paths = set()
    for path in selected:
        validate_text(str(path), "FASTQ path")
        if layout == "single":
            sample = path.name.removesuffix(".gz").removesuffix(".fastq")
            record = {"sample": validate_sample_name(sample), "read1": str(path.resolve()), "read2": None}
        else:
            match = PAIRED_NAME.fullmatch(path.name)
            if match is None:
                raise ValueError(f"Paired FASTQ has no valid _R1/_R2 marker: {path.name}")
            prefix, suffix = match["prefix"], match["suffix"] or ""
            sample = validate_sample_name(prefix + suffix)
            mates = []
            for number in (1, 2):
                plain = path.parent / f"{prefix}_R{number}{suffix}.fastq"
                candidates = [p for p in (plain, Path(str(plain) + ".gz")) if p.is_file()]
                if not candidates:
                    raise ValueError(f"Missing R{number} mate for {path.name}; expected {plain.name}[.gz]")
                if len(candidates) != 1:
                    raise ValueError(f"Duplicate R{number} FASTQ files for sample {sample}: plain and gzip copies")
                validate_text(str(candidates[0]), "FASTQ path")
                mates.append(str(candidates[0].resolve()))
            record = {"sample": sample, "read1": mates[0], "read2": mates[1]}
        prior = samples.get(record["sample"])
        if prior is not None:
            # A paired sample is normally encountered once for each selected mate.
            if layout == "paired" and prior == record:
                continue
            raise ValueError(f"Duplicate sample name after FASTQ normalization: {record['sample']}")
        paths = {record["read1"]} | ({record["read2"]} if record["read2"] else set())
        if paths & selected_paths or len(paths) != (2 if layout == "paired" else 1):
            raise ValueError(f"FASTQ path is reused by multiple samples: {record['sample']}")
        selected_paths.update(paths)
        samples[record["sample"]] = record
    return [samples[name] for name in sorted(samples)]


def build_manifest(directory, pattern="*.fastq*", layout="single", max_samples=None):
    if max_samples is not None and (type(max_samples) is not int or max_samples < 1):
        raise ValueError("Maximum sample count must be a positive integer")
    samples = discover_samples(directory, pattern, layout)
    if max_samples is not None and len(samples) > max_samples:
        raise ValueError(f"Discovered {len(samples)} samples, exceeding the configured maximum of {max_samples}")
    return {"read_layout": layout, "samples": samples}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq-dir", type=Path, required=True)
    parser.add_argument("--glob", dest="pattern", default="*.fastq*")
    parser.add_argument("--layout", choices=("single", "paired"), default="single")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-samples", type=int)
    args = parser.parse_args()
    try:
        manifest = build_manifest(args.fastq_dir, args.pattern, args.layout, args.max_samples)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_name(args.output.name + ".tmp")
        temporary.write_text(json.dumps(manifest, indent=2) + "\n")
        temporary.replace(args.output)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
