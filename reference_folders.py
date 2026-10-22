#!/usr/bin/env python3
"""Discover NCBI reference folders and prepare read-only-source Bowtie2 caches."""

import argparse
import fcntl
import gzip
import hashlib
import json
import re
import subprocess
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote, unquote

FASTA_SUFFIXES = ('.fa', '.fasta', '.fna')
RNA_TYPES = {'rrna', 'trna', 'tmrna', 'ncrna', 'srna', 'misc_rna', 'snrna',
             'snorna', 'scrna', 'mirna', 'pirna', 'antisense_rna', 'ribozyme'}
TRANSCRIPTS = {'mrna', 'transcript', 'primary_transcript'}
EXON_TYPES = {'exon', 'five_prime_utr', 'three_prime_utr', 'utr', '5utr', '3utr'}
CACHE_VERSION = 5


def text_open(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def stem(path):
    name = path.name.removesuffix('.gz')
    return name.rsplit('.', 1)[0]


def is_fasta(path):
    return path.name.removesuffix('.gz').lower().endswith(FASTA_SUFFIXES)


def is_genbank(path):
    return path.name.removesuffix('.gz').lower().endswith(('.gbff', '.gb', '.gbk'))


def is_rna(path):
    name = path.name.lower()
    return 'rna_from_genomic' in name or re.search(r'(?:^|[_.-])rna(?:[_.-]|$)', name) is not None or any(p.name.lower().endswith('_rna') for p in path.parents)


def files_under(folder):
    folder = Path(folder).resolve()
    if not folder.is_dir():
        raise ValueError(f'Reference directory does not exist: {folder}')
    return sorted(p for p in folder.rglob('*') if p.is_file())


def complete_index(prefix):
    return any(all(Path(f'{prefix}.{part}.{ext}').is_file() and Path(f'{prefix}.{part}.{ext}').stat().st_size > 0
                   for part in ('1', '2', '3', '4', 'rev.1', 'rev.2')) for ext in ('bt2', 'bt2l'))


def index_prefixes(files):
    return sorted({Path(str(p).rsplit('.1.', 1)[0]) for p in files
                   if re.search(r'(?<!rev)\.1\.bt2l?$', p.name) and complete_index(str(p).rsplit('.1.', 1)[0])})


def signature(files):
    data = [(str(p), p.stat().st_size, p.stat().st_mtime_ns) for p in sorted(files)]
    return hashlib.sha256(json.dumps([CACHE_VERSION, data]).encode()).hexdigest()[:20]


def combine_fasta(paths, destination):
    """Copy sequences one line at a time; retain only identifiers and lengths."""
    lengths = {}
    descriptions = {}
    with destination.open('w') as output:
        for path in paths:
            ident = None
            with text_open(path) as source:
                for line in source:
                    if line.startswith('>'):
                        description = line[1:].strip()
                        ident = description.split()[0] if description else ''
                        if not ident or ident in lengths:
                            raise ValueError(f'Empty or duplicate FASTA identifier {ident!r} in {path}')
                        lengths[ident] = 0
                        descriptions[ident] = description
                        output.write(f'>{ident}\n')
                    elif line.strip():
                        if ident is None:
                            raise ValueError(f'Sequence before FASTA header in {path}')
                        sequence = ''.join(line.split())
                        lengths[ident] += len(sequence)
                        output.write(sequence + '\n')
    if not lengths or any(n == 0 for n in lengths.values()):
        raise ValueError(f'Empty reference sequence in {paths}')
    return lengths, descriptions


def combine_reference(paths, destination):
    """Use FASTA sequences when present, or the ORIGIN sequence in GBFF."""
    if all(is_fasta(p) for p in paths):
        return combine_fasta(paths, destination)
    from Bio import SeqIO
    lengths, descriptions, fasta_digests, genbank_seen = {}, {}, {}, set()
    # Put FASTA first so differently named GBFF annotation files can carry
    # their same sequences without duplicating those reference records.
    paths = sorted(paths, key=lambda path: not is_fasta(path))
    with destination.open('w') as output:
        for path in paths:
            with text_open(path) as source:
                for record in SeqIO.parse(source, 'fasta' if is_fasta(path) else 'genbank'):
                    if is_genbank(path):
                        if record.id in genbank_seen:
                            raise ValueError(f'Duplicate GenBank identifier {record.id} in {path}')
                        genbank_seen.add(record.id)
                    if record.id in lengths:
                        if not is_genbank(path) or record.id not in fasta_digests:
                            raise ValueError(f'Duplicate reference identifier {record.id} in {path}')
                        if len(record) != lengths[record.id]:
                            raise ValueError(f'GenBank/FASTA lengths disagree for {record.id}')
                        digest = hashlib.sha256()
                        try:
                            for start in range(0, len(record), 65536):
                                digest.update(str(record.seq[start:start + 65536]).upper().encode())
                        except ValueError as exc:
                            raise ValueError(f'GBFF has no usable ORIGIN sequence for {record.id}: {path}') from exc
                        if digest.digest() != fasta_digests[record.id]:
                            raise ValueError(f'GenBank/FASTA sequences disagree for {record.id}')
                        continue
                    if not len(record):
                        raise ValueError(f'Empty reference sequence {record.id} in {path}')
                    lengths[record.id] = len(record)
                    descriptions[record.id] = f'{record.id} {record.description}'
                    output.write(f'>{record.id}\n')
                    digest = hashlib.sha256()
                    try:
                        for start in range(0, len(record), 65536):
                            sequence = str(record.seq[start:start + 65536])
                            output.write(sequence + '\n')
                            digest.update(sequence.upper().encode())
                    except ValueError as exc:
                        raise ValueError(f'GBFF has no usable ORIGIN sequence for {record.id}: {path}') from exc
                    if is_fasta(path):
                        fasta_digests[record.id] = digest.digest()
    if not lengths:
        raise ValueError(f'No reference sequences found in {paths}')
    return lengths, descriptions


def attrs(value):
    result = {}
    for entry in value.strip().split(';'):
        if '=' in entry:
            key, val = entry.split('=', 1)
            result[key] = unquote(val)
    return result


def biological_type(value):
    low = value.lower()
    return low == 'cds' or low in TRANSCRIPTS or low in EXON_TYPES or low in RNA_TYPES or ('rna' in low and low not in {'mrna_region'})


def write_feature(handle, seqid, kind, start, end, strand, ident, extra=None):
    attributes = {'ID': ident, 'locus_tag': ident}
    attributes.update(extra or {})
    encoded = ';'.join(f'{key}={quote(str(val), safe="_:.-|,")}' for key, val in attributes.items())
    handle.write(f'{seqid}\treference_folders\t{kind}\t{start}\t{end}\t.\t{strand}\t.\t{encoded}\n')


def normalized_gff(paths, destination, lengths, require_all=False):
    """Normalize identifier attributes while preserving explicit annotation types."""
    metadata = defaultdict(list)
    annotated = set()
    count = 0
    with destination.open('w') as output:
        output.write('##gff-version 3\n')
        for path in paths:
            with text_open(path) as source:
                for number, line in enumerate(source, 1):
                    if line.startswith('##FASTA'):
                        break
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) != 9:
                        raise ValueError(f'Malformed GFF record at {path}:{number}')
                    seqid, _, kind, start, end, _, strand, _, attr = parts
                    if not biological_type(kind) and kind.lower() != 'gene':
                        continue
                    start, end = int(start), int(end)
                    if seqid not in lengths or start < 1 or end < start or end > lengths[seqid]:
                        raise ValueError(f'Annotation and FASTA coordinates disagree at {path}:{number} ({seqid}:{start}-{end})')
                    a = attrs(attr)
                    ident = next((a[key] for key in ('locus_tag', 'locus', 'gene', 'ID', 'Name', 'Parent') if a.get(key)), f'feature_{count + 1}')
                    # Preserve hierarchy attributes for existing positional profiles.
                    write_feature(output, seqid, kind, start, end, strand, ident,
                                  {key: a[key] for key in ('ID', 'Parent', 'gene', 'Name') if a.get(key)})
                    for key in ('locus_tag', 'locus', 'gene', 'ID', 'Name'):
                        if a.get(key):
                            metadata[a[key]].append(kind)
                    count += biological_type(kind)
                    if biological_type(kind):
                        annotated.add(seqid)
    if not count:
        raise ValueError(f'No CDS or RNA annotation features found in {paths}')
    if require_all and annotated != set(lengths):
        raise ValueError(f'Missing RNA annotation for FASTA sequences: {", ".join(sorted(set(lengths) - annotated)[:10])}')
    return metadata


def genbank_gff(paths, destination, lengths):
    from Bio import SeqIO
    metadata = defaultdict(list)
    found = set()
    count = 0
    with destination.open('w') as output:
        output.write('##gff-version 3\n')
        for path in paths:
            with text_open(path) as source:
                for record in SeqIO.parse(source, 'genbank'):
                    if record.id not in lengths:
                        raise ValueError(f'GenBank sequence {record.id} has no matching FASTA sequence in {path}')
                    if len(record) != lengths[record.id]:
                        raise ValueError(f'GenBank/FASTA lengths disagree for {record.id}')
                    found.add(record.id)
                    transcript = any(f.qualifiers.get('mol_type', [''])[0].lower() in {'mrna', 'transcribed rna', 'ncrna', 'rrna', 'trna'} for f in record.features if f.type == 'source')
                    if transcript and not any(f.type.lower() in TRANSCRIPTS or f.type.lower() in RNA_TYPES for f in record.features):
                        kind = next((f.qualifiers['mol_type'][0] for f in record.features if f.type == 'source' and f.qualifiers.get('mol_type')), 'transcript')
                        kind = kind if kind.lower() != 'transcribed rna' else 'transcript'
                        write_feature(output, record.id, kind, 1, len(record), '+', record.id)
                        count += 1
                    for feature in record.features:
                        if not biological_type(feature.type) and feature.type != 'gene':
                            continue
                        if feature.location is None:
                            raise ValueError(f'Feature without location in {path}: {record.id}')
                        ident = next((feature.qualifiers[key][0] for key in ('locus_tag', 'locus', 'gene', 'transcript_id', 'protein_id') if feature.qualifiers.get(key)), f'{record.id}_{feature.type}_{count + 1}')
                        for part in feature.location.parts:
                            start, end = int(part.start) + 1, int(part.end)
                            if start < 1 or end > lengths[record.id] or end < start:
                                raise ValueError(f'GenBank feature outside FASTA sequence: {record.id}')
                            strand = '-' if part.strand == -1 else '+'
                            write_feature(output, record.id, feature.type, start, end, strand, ident)
                            if feature.type.lower() in TRANSCRIPTS:
                                # GenBank joins describe transcript exon pieces. Preserve
                                # these separately so genomic introns are not called RNA.
                                write_feature(output, record.id, 'exon', start, end, strand, ident)
                        for key in ('locus_tag', 'locus', 'gene', 'transcript_id', 'protein_id'):
                            for value in feature.qualifiers.get(key, []):
                                metadata[value].append(feature.type)
                        count += biological_type(feature.type)
    if found != set(lengths):
        raise ValueError(f'Missing GenBank annotation for FASTA sequences: {", ".join(sorted(set(lengths) - found)[:10])}')
    if not count:
        raise ValueError(f'No CDS or RNA annotation features found in {paths}')
    return metadata


def header_rna_gff(descriptions, lengths, metadata, destination):
    with destination.open('w') as output:
        output.write('##gff-version 3\n')
        for seqid, description in descriptions.items():
            a = dict(re.findall(r'\[([^=\]]+)=([^\]]*)\]', description))
            kind = a.get('gbkey', '')
            ident = next((a[k] for k in ('locus_tag', 'gene', 'transcript_id') if a.get(k)), seqid)
            matches = [t for k in ('locus_tag', 'gene', 'transcript_id') for t in metadata.get(a.get(k, ''), [])]
            if kind.lower() in {'rna', 'gene', ''}:
                kinds = sorted({t for t in matches if biological_type(t) and t.lower() != 'cds'})
                if len(kinds) == 1:
                    kind = kinds[0]
            if not biological_type(kind) or kind.lower() in TRANSCRIPTS or kind.lower() == 'cds':
                raise ValueError(f'Cannot determine RNA/CDS coordinates for {seqid}; provide matching RNA GBFF annotation (or noncoding RNA gbkey/locus_tag FASTA headers)')
            write_feature(output, seqid, kind, 1, lengths[seqid], '+', ident)


def annotation_seqids(path):
    result = set()
    with text_open(path) as source:
        if is_genbank(path):
            accession, version = None, None
            for line in source:
                if line.startswith('ACCESSION'):
                    accession = line.split()[1]
                elif line.startswith('VERSION'):
                    version = line.split()[1]
                elif line.startswith('//'):
                    if version or accession:
                        result.add(version or accession)
                    accession, version = None, None
        else:
            for line in source:
                if line.startswith('##FASTA'):
                    break
                if line.strip() and not line.startswith('#'):
                    result.add(line.split('\t', 1)[0])
    return result


def matching_annotations(fastas, files, extension, lengths):
    wanted = {stem(p) for p in fastas}
    candidates = [p for p in files if p.name.removesuffix('.gz').lower().endswith(extension)]
    exact = [p for p in candidates if stem(p) in wanted]
    if exact:
        return exact
    # Renamed annotations are selected by reference sequence IDs rather than
    # requiring the user to restore a shared basename.
    selected, seen = [], set()
    for path in candidates:
        seqids = annotation_seqids(path)
        if not seqids.intersection(lengths):
            continue
        if seen.intersection(seqids):
            raise ValueError(f'Ambiguous annotation files contain the same reference sequences: {path}')
        selected.append(path)
        seen.update(seqids)
    return selected


def validate_cache_location(folder, cache):
    folder, cache = Path(folder).resolve(), Path(cache).resolve()
    if cache == folder or folder in cache.parents:
        raise ValueError(f'Reference cache must be outside the input reference folder: {folder}')


def fasta_identifiers(paths):
    identifiers = set()
    for path in paths:
        with text_open(path) as source:
            for line in source:
                if line.startswith('>') and line[1:].strip():
                    identifiers.add(line[1:].split()[0])
    return identifiers


def build_index(fasta, prefix, threads):
    if complete_index(prefix):
        return
    with Path(f'{prefix}.build.log').open('w') as log:
        result = subprocess.run(['bowtie2-build', '--threads', str(threads), str(fasta), str(prefix)], stdout=log, stderr=subprocess.STDOUT)
    if result.returncode or not complete_index(prefix):
        raise ValueError(f'bowtie2-build failed; see {prefix}.build.log')


def extract_genomic_rna(fasta, gff, destination):
    """Extract annotated noncoding RNA when no separate RNA FASTA was supplied."""
    import pysam
    if not Path(str(fasta) + '.fai').exists():
        pysam.faidx(str(fasta))
    descriptions = {}
    with pysam.FastaFile(str(fasta)) as reference, destination.open('w') as output, gff.open() as annotations:
        for i, line in enumerate(annotations):
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if fields[2].lower() not in RNA_TYPES:
                continue
            seq = reference.fetch(fields[0], int(fields[3]) - 1, int(fields[4]))
            if fields[6] == '-':
                seq = seq.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]
            a = attrs(fields[8])
            seqid = f'annotated_rna_{i}'
            description = f'{seqid} [gbkey={fields[2]}] [locus_tag={a["locus_tag"]}]'
            descriptions[seqid] = description
            output.write(f'>{description}\n{seq}\n')
    return descriptions


def prepare_organism(folder, cache, role, threads):
    validate_cache_location(folder, cache)
    files = files_under(folder)
    fastas = [p for p in files if is_fasta(p) and 'cds_from_genomic' not in p.name.lower()]
    # A GBFF can supply both sequence and annotation. Avoid reading the same
    # sequence twice when its preferred FASTA counterpart also exists.
    fasta_stems = {stem(p) for p in fastas}
    unmatched_gbffs = [p for p in files if is_genbank(p) and stem(p) not in fasta_stems
                      and 'cds_from_genomic' not in p.name.lower()]
    rna = [p for p in fastas if is_rna(p)]
    dna = [p for p in fastas if p not in rna]
    if unmatched_gbffs:
        # A renamed annotation may not have "RNA" in its filename. Its
        # sequence IDs determine which matching FASTA stage it belongs to.
        rna_ids, dna_ids = fasta_identifiers(rna), fasta_identifiers(dna)
        for path in unmatched_gbffs:
            seqids = annotation_seqids(path)
            matches_rna, matches_dna = bool(seqids & rna_ids), bool(seqids & dna_ids)
            if matches_rna and matches_dna:
                raise ValueError(f'GBFF mixes genomic and RNA reference sequences: {path}')
            (rna if matches_rna or not matches_dna and is_rna(path) else dna).append(path)
    if not dna:
        raise ValueError(f'No genomic FASTA or sequence-bearing GBFF found in {folder}')
    # Different assembled genomes in one directory cannot be silently mixed.
    assemblies = {re.search(r'GC[AF]_\d+\.\d+', p.name).group(0) for p in dna if re.search(r'GC[AF]_\d+\.\d+', p.name)}
    if len(assemblies) > 1:
        raise ValueError(f'Multiple genome assemblies in {folder}: {sorted(assemblies)}')
    relevant = [p for p in files if is_fasta(p) or p.name.removesuffix('.gz').lower().endswith(('.gff', '.gff3', '.gbff', '.gb', '.gbk'))]
    directory = cache / f'{role}-{signature(relevant)}'
    directory.mkdir(parents=True, exist_ok=True)
    ready = directory / 'ready.json'
    with (directory / '.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if ready.exists():
            result = json.loads(ready.read_text())
            if all(stage is None or complete_index(stage['index']) and Path(stage['gff']).is_file() for stage in result.values()):
                return result
        dna_fasta = directory / 'dna.fa'
        lengths, _ = combine_reference(dna, dna_fasta)
        gff = matching_annotations(dna, files, ('.gff', '.gff3'), lengths)
        gbff = matching_annotations(dna, files, ('.gbff', '.gb', '.gbk'), lengths) if not gff else []
        dna_gff = directory / 'dna.gff3'
        if gff:
            metadata = normalized_gff(gff, dna_gff, lengths)
        elif gbff:
            metadata = genbank_gff(gbff, dna_gff, lengths)
        else:
            raise ValueError(f'No matching genomic GFF/GBFF annotations found in {folder}')
        result = {'dna': {'index': str(directory / 'dna'), 'gff': str(dna_gff)}, 'rna': None}
        rna_fasta, rna_gff = directory / 'rna.fa', directory / 'rna.gff3'
        if rna:
            rna_lengths, descriptions = combine_reference(rna, rna_fasta)
            rna_gffs = matching_annotations(rna, files, ('.gff', '.gff3'), rna_lengths)
            rna_gbffs = matching_annotations(rna, files, ('.gbff', '.gb', '.gbk'), rna_lengths) if not rna_gffs else []
            if rna_gffs:
                normalized_gff(rna_gffs, rna_gff, rna_lengths, require_all=True)
            elif rna_gbffs:
                genbank_gff(rna_gbffs, rna_gff, rna_lengths)
            else:
                header_rna_gff(descriptions, rna_lengths, metadata, rna_gff)
        else:
            descriptions = extract_genomic_rna(dna_fasta, dna_gff, rna_fasta)
            if descriptions:
                rna_lengths = {ident: len_ for ident, len_ in fasta_lengths(rna_fasta).items()}
                header_rna_gff(descriptions, rna_lengths, metadata, rna_gff)
        if rna or descriptions:
            result['rna'] = {'index': str(directory / 'rna'), 'gff': str(rna_gff)}
            build_index(rna_fasta, directory / 'rna', threads)
        build_index(dna_fasta, directory / 'dna', threads)
        ready.write_text(json.dumps(result, indent=2) + '\n')
        return result


def fasta_lengths(path):
    lengths = {}
    with path.open() as source:
        for line in source:
            if line.startswith('>'):
                ident = line[1:].split()[0]
                lengths[ident] = 0
            elif line.strip():
                lengths[ident] += len(line.strip())
    return lengths


def prepare_decoy(folder, cache, threads):
    validate_cache_location(folder, cache)
    files = files_under(folder)
    indices = index_prefixes(files)
    if len(indices) > 1:
        raise ValueError(f'Multiple complete decoy Bowtie2 indices in {folder}: {indices}')
    if indices:
        return {'index': str(indices[0])}
    fastas = [p for p in files if is_fasta(p)]
    fasta_stems = {stem(p) for p in fastas}
    sequences = fastas + [p for p in files if is_genbank(p) and stem(p) not in fasta_stems]
    if not sequences:
        raise ValueError(f'No complete Bowtie2 index, FASTA or sequence-bearing GBFF in decoy folder: {folder}')
    directory = cache / f'decoy-{signature(sequences)}'
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / '.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not complete_index(directory / 'decoy'):
            combine_reference(sequences, directory / 'decoy.fa')
            build_index(directory / 'decoy.fa', directory / 'decoy', threads)
    return {'index': str(directory / 'decoy')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--decoy-dir', required=True, type=Path)
    parser.add_argument('--bacteria-dir', required=True, type=Path)
    parser.add_argument('--host-dir', type=Path)
    parser.add_argument('--cache-dir', required=True, type=Path)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('--threads must be positive')
    args.cache_dir = args.cache_dir.resolve()
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    try:
        manifest = {'decoy': prepare_decoy(args.decoy_dir, args.cache_dir, args.threads),
                    'bacteria': prepare_organism(args.bacteria_dir, args.cache_dir, 'bacteria', args.threads),
                    'host': prepare_organism(args.host_dir, args.cache_dir, 'host', args.threads) if args.host_dir else None}
    except (ValueError, OSError) as exc:
        parser.exit(1, f'ERROR: {exc}\n')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + '.tmp')
    temporary.write_text(json.dumps(manifest, indent=2) + '\n')
    temporary.replace(args.output)


if __name__ == '__main__':
    main()
