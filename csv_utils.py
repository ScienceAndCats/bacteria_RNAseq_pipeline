"""CSV limits and streaming diagnostics for long featureCounts annotations."""

import csv
from pathlib import Path
import sys


# Match longread_bacteria_RNAseq_pipeline's large-annotation fix (04533fb).
CSV_FIELD_SIZE_LIMIT = 64 * 1024 * 1024
LARGE_FIELD_THRESHOLD = 131072
LARGE_FIELD_REPORT = "large_csv_fields.txt"
ANNOTATION_COLUMNS = ("Geneid", "Chr", "Start", "End", "Strand", "Length")


def configure_csv_field_size_limit():
    """Apply the limit in each process that actually reads CSV/TSV data."""
    csv.field_size_limit(CSV_FIELD_SIZE_LIMIT)


class LargeFieldReport:
    """Write complete offending fields incrementally, including on failures."""

    def __init__(self, path):
        self.path = Path(path).resolve()
        self.cells = 0
        self.records = 0
        self.parse_errors = 0

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = self.path.open("w", encoding="utf-8", newline="")
        self.handle.write(
            "Large CSV/TSV field diagnostics\n"
            f"Reporting threshold: {LARGE_FIELD_THRESHOLD} characters (old CSV limit)\n"
            f"CSV field limit: {CSV_FIELD_SIZE_LIMIT} characters (64 Mi)\n"
            "Each source record with an oversized cell is recorded once.\n"
            "Annotation values and oversized cells are included IN FULL, without truncation.\n"
        )
        self.handle.flush()
        print(f"CSV field diagnostics: {self.path}", file=sys.stderr, flush=True)
        return self

    def record(self, path, record_number, line_number, header, row):
        oversized = [(index, value) for index, value in enumerate(row)
                     if len(value) > LARGE_FIELD_THRESHOLD]
        if not oversized:
            return
        self.records += 1
        self.cells += len(oversized)
        self.handle.write(
            f"\n--- Source record {self.records} ---\nSource file: {path}\n"
            f"Table record: {record_number} (header = 1)\n"
            f"Physical line ending record: {line_number}\n"
        )
        for index, value in oversized:
            name = header[index] if index < len(header) else f"column_{index + 1}"
            self.handle.write(f"Oversized column: {name} ({index + 1}); characters: {len(value)}\n")
        # Include all annotation context, even when e.g. only Start is long.
        included = set()
        for index, name in enumerate(header):
            if index < len(row) and name in ANNOTATION_COLUMNS:
                self._value(name, row[index])
                included.add(index)
        for index, value in oversized:
            if index not in included:
                name = header[index] if index < len(header) else f"column_{index + 1}"
                self._value(name, value)
        self.handle.flush()

    def _value(self, name, value):
        self.handle.write(f"\n{name} (full value; {len(value)} characters):\n<<<BEGIN VALUE>>>\n")
        self.handle.write(value)
        self.handle.write("\n<<<END VALUE>>>\n")

    def parser_error(self, path, record_number, lines, remaining_input, error):
        self.parse_errors += 1
        self.handle.write(
            f"\n--- CSV parser error ---\nSource file: {path}\n"
            f"Table record: {record_number} (header = 1)\nError: {error}\n"
            "The parser could not finish this record. Raw text starting at this record\n"
            "follows, including the remaining table so a multiline field is not cut off.\n"
        )
        if lines:
            self.handle.write(f"Physical lines: {lines[0][0]}-{lines[-1][0]}\n")
        self.handle.write("<<<BEGIN RAW RECORD>>>\n")
        for _, line in lines:
            self.handle.write(line)
        # The failure may occur before the end of a quoted multiline field.
        # Stream the remainder rather than truncate it or raise the parsing cap.
        for line in remaining_input:
            self.handle.write(line)
        self.handle.write("\n<<<END RAW RECORD>>>\n")
        self.handle.flush()

    def __exit__(self, exc_type, exc, traceback):
        try:
            self.handle.write(
                f"\nStatus: {'FAILED' if exc_type else 'COMPLETED'}\n"
                f"Oversized cells: {self.cells}\nSource records with oversized cells: {self.records}\n"
                f"CSV parser errors: {self.parse_errors}\n"
            )
            if not self.cells and not self.parse_errors:
                self.handle.write("No cells exceeded the old CSV field limit in the records read.\n")
            if exc_type:
                self.handle.write(f"Failure: {exc_type.__name__}: {exc}\n")
        finally:
            self.handle.close()


def featurecounts_rows(handle, report=None):
    """Read TSV rows with source context and bounded, per-record diagnostics."""
    configure_csv_field_size_limit()
    pending = []

    def lines():
        for number, line in enumerate(handle, 1):
            # Comments/blanks inside a quoted multiline cell are data.
            if not pending and (line.startswith("#") or not line.strip()):
                continue
            pending.append((number, line))
            yield line

    reader = csv.reader(lines(), delimiter="\t")
    header = None
    record_number = 0
    source = str(Path(handle.name).resolve())
    while True:
        pending.clear()
        try:
            row = next(reader)
        except StopIteration:
            return
        except csv.Error as exc:
            if report is not None:
                report.parser_error(source, record_number + 1, pending, handle, exc)
            detail = f"; full raw record saved in {report.path}" if report else ""
            raise ValueError(f"Cannot parse {source}, record {record_number + 1}: {exc}{detail}") from exc
        record_number += 1
        if header is None:
            header = row
        if report is not None:
            report.record(source, record_number, pending[-1][0], header, row)
        pending.clear()
        yield row
