#!/usr/bin/env python3
"""Count mutually exclusive annotation categories for mapped reads or proper pairs.

Inputs must group paired primary alignments by query name (Bowtie2 --reorder).
Coordinates come from aligned CIGAR blocks, never the unsequenced fragment span.
"""

import argparse
from bisect import bisect_left
from collections import defaultdict
import gzip
from itertools import groupby
import json
from pathlib import Path
from urllib.parse import unquote

import pysam


CATEGORIES = ("rRNA", "CDS", "other_RNA", "intergenic_unannotated")
BIN_SIZE = 16384
TRANSCRIPTS = {"mrna", "transcript", "primary_transcript"}
UTR_TYPES = {"utr", "five_prime_utr", "three_prime_utr", "5utr", "3utr", "5'utr", "3'utr"}
RNA_TYPES = {
    "trna", "tmrna", "ncrna", "nc_rna", "srna", "misc_rna", "snrna", "scrna",
    "snorna", "mirna", "pirna", "antisense_rna", "rnase_p_rna", "rnase_mrp_rna",
    "telomerase_rna", "guide_rna", "ribozyme", "lnc_rna", "lncrna",
}


def attributes(text):
    result = {}
    for item in text.strip().strip(";").split(";"):
        if "=" in item:
            key, value = item.strip().split("=", 1)
            result[key] = unquote(value)
    return result


def gff_rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Expected 9 GFF columns at {path}:{number}")
            start, end = int(fields[3]) - 1, int(fields[4])
            if start < 0 or end <= start:
                raise ValueError(f"Invalid GFF coordinates at {path}:{number}")
            yield fields[0], start, end, fields[6], fields[2].lower(), attributes(fields[8])


def type_category(feature_type):
    if feature_type == "cds":
        return "CDS"
    if feature_type == "rrna":
        return "rRNA"
    if feature_type in TRANSCRIPTS or feature_type in UTR_TYPES or feature_type in RNA_TYPES or "rna" in feature_type:
        return "other_RNA"
    return None


class AnnotationIndex:
    """Merged category intervals plus a binned CDS-to-feature lookup."""

    def __init__(self, path, reference_kind="dna"):
        # A genomic transcript interval includes introns. Where exon records are
        # available, use those exons instead of the encompassing RNA interval.
        exon_parents = set()
        parent_categories = {}
        parent_links = {}
        for _, _, _, _, typ, attrs in gff_rows(path):
            if typ == "exon":
                exon_parents.update(attrs.get("Parent", "").split(","))
            if attrs.get("ID"):
                parent_categories[attrs["ID"]] = type_category(typ)
                parent_links[attrs["ID"]] = attrs.get("Parent", "").split(",")

        def ancestor_category(ident, seen=None):
            seen = set() if seen is None else seen
            if ident in seen:
                return None
            seen.add(ident)
            own = parent_categories.get(ident)
            found = [own] if own else []
            for parent in parent_links.get(ident, ()):
                parent_category = ancestor_category(parent, seen)
                if parent_category:
                    found.append(parent_category)
            return next((category for category in CATEGORIES if category in found), None)

        intervals = defaultdict(list)
        self.cds_bins = defaultdict(list)
        self.references = set()
        self.reference_ends = {}
        for ref, start, end, strand, typ, attrs in gff_rows(path):
            self.references.add(ref)
            category = type_category(typ)
            if typ == "exon":
                categories = [ancestor_category(parent) for parent in attrs.get("Parent", "").split(",")]
                category = "rRNA" if "rRNA" in categories else "other_RNA"
            if category is None:
                continue
            if reference_kind == "dna":
                if attrs.get("ID") in exon_parents:
                    continue
                # mRNA without exons cannot define exonic genomic positions.
                if typ in TRANSCRIPTS:
                    continue
            self.reference_ends[ref] = max(self.reference_ends.get(ref, 0), end)
            intervals[(ref, strand, category)].append((start, end))
            if category == "CDS":
                ident = next((attrs[key] for key in ("locus_tag", "locus", "gene", "ID", "Parent") if attrs.get(key)), f"{ref}:{start + 1}-{end}:{strand}")
                for bin_number in range(start // BIN_SIZE, (end - 1) // BIN_SIZE + 1):
                    self.cds_bins[(ref, strand, bin_number)].append((start, end, ident))
        if not intervals:
            raise ValueError(f"No relevant CDS, rRNA, exon/UTR or other RNA features in {path}")
        self.intervals = {}
        for key, values in intervals.items():
            merged = []
            for start, end in sorted(values):
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
                else:
                    merged.append((start, end))
            self.intervals[key] = ([item[0] for item in merged], [item[1] for item in merged])

    def classify(self, ref, blocks, strand):
        strands = ("+", "-", ".", "?") if strand is None else (strand, ".", "?")
        for category in CATEGORIES[:-1]:
            for annotation_strand in strands:
                starts, ends = self.intervals.get((ref, annotation_strand, category), ((), ()))
                for start, end in blocks:
                    index = bisect_left(starts, end) - 1
                    if index >= 0 and ends[index] > start:
                        features = set()
                        if category == "CDS":
                            for block_start, block_end in blocks:
                                for selected_strand in strands:
                                    for bin_number in range(block_start // BIN_SIZE, (block_end - 1) // BIN_SIZE + 1):
                                        for feature_start, feature_end, ident in self.cds_bins.get((ref, selected_strand, bin_number), ()):
                                            if feature_start < block_end and feature_end > block_start:
                                                features.add(ident)
                        return category, features
        return "intergenic_unannotated", set()


def classify_alignments(sam, gff, paired=False, strand="forward", reference_kind="dna"):
    annotation = AnnotationIndex(gff, reference_kind)
    counts = dict.fromkeys(CATEGORIES, 0)
    features = set()
    aligned = 0

    def add(reads):
        nonlocal aligned
        first = reads[0]
        reference = first.reference_name
        if any(read.reference_name != reference for read in reads):
            raise ValueError(f"Proper pair has different reference names: {first.query_name}")
        selected_strand = None if strand == "unstranded" else ("-" if first.is_reverse else "+")
        if strand == "reverse":
            selected_strand = "+" if selected_strand == "-" else "-"
        blocks = [block for read in reads for block in read.get_blocks() if block[0] < block[1]]
        category, matched_features = annotation.classify(reference, blocks, selected_strand)
        counts[category] += 1
        features.update(matched_features)
        aligned += 1

    mode = "rb" if str(sam).endswith((".bam", ".cram")) else "r"
    with pysam.AlignmentFile(str(sam), mode) as handle:
        reference_lengths = dict(zip(handle.references, handle.lengths))
        missing = annotation.references - reference_lengths.keys()
        if missing:
            raise ValueError(f"GFF reference names are absent from the SAM header: {', '.join(sorted(missing)[:10])}")
        too_short = [ref for ref, end in annotation.reference_ends.items() if end > reference_lengths[ref]]
        if too_short:
            raise ValueError(f"GFF coordinates exceed SAM reference lengths: {', '.join(sorted(too_short)[:10])}")
        for query_name, group in groupby(handle, key=lambda read: read.query_name):
            primary = [read for read in group if not (read.is_unmapped or read.is_secondary or read.is_supplementary)]
            if not paired:
                for read in primary:
                    add([read])
                continue
            proper = [read for read in primary if read.is_proper_pair and not read.mate_is_unmapped]
            if not proper:
                continue
            if len(proper) != 2 or sum(read.is_read1 for read in proper) != 1 or sum(read.is_read2 for read in proper) != 1:
                raise ValueError(f"Expected adjacent primary mates for {query_name}; use Bowtie2 --reorder")
            proper.sort(key=lambda read: not read.is_read1)
            add(proper)
    return {"aligned": aligned, "categories": counts, "cds_features_with_reads": sorted(features)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sam", required=True, type=Path)
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--paired", action="store_true")
    parser.add_argument("--strand", choices=("forward", "reverse", "unstranded"), default="forward")
    parser.add_argument("--reference-kind", choices=("rna", "dna"), default="dna")
    args = parser.parse_args()
    try:
        result = classify_alignments(args.sam, args.gff, args.paired, args.strand, args.reference_kind)
    except (ValueError, OSError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
