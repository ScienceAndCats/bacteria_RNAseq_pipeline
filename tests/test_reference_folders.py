"""Reference discovery tests use tiny NCBI-style layouts and mocked index builds."""
import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import reference_folders as refs


class ReferenceFoldersTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cache = self.root / 'cache'
        self.cache.mkdir()

    @staticmethod
    def build(fasta, prefix, threads):
        for part in ('1', '2', '3', '4', 'rev.1', 'rev.2'):
            Path(f'{prefix}.{part}.bt2').write_text('index')

    def write_gz(self, path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(path, 'wt') as output:
            output.write(value)

    def bacteria(self):
        folder = self.root / 'PAO1'
        self.write_gz(folder / 'GCF_000006765.1_genomic.fna.gz', '>chromosome\n' + 'ACGT' * 100 + '\n')
        self.write_gz(folder / 'GCF_000006765.1_genomic.gff.gz',
                      '##gff-version 3\nchromosome\ttest\tCDS\t1\t100\t.\t+\t0\tID=cds;gene=gene1\n'
                      'chromosome\ttest\trRNA\t201\t300\t.\t-\t.\tID=rna;locus_tag=rna1\n')
        self.write_gz(folder / 'GCF_000006765.1_cds_from_genomic.fna.gz', '>excluded_CDS\nACGT\n')
        return folder

    def test_bacterial_ncbi_layout_and_cache(self):
        folder = self.bacteria()
        self.write_gz(folder / 'GCF_000006765.1_rna_from_genomic.fna.gz',
                      '>rna [locus_tag=rna1] [gbkey=rRNA]\n' + 'ACGT' * 25 + '\n')
        with patch.object(refs, 'build_index', side_effect=self.build) as builder:
            manifest = refs.prepare_organism(folder, self.cache, 'bacteria', 2)
            self.assertEqual(builder.call_count, 2)
            self.assertNotIn('excluded_CDS', Path(manifest['dna']['index'] + '.fa').read_text())
            self.assertIn('\trRNA\t1\t100\t', Path(manifest['rna']['gff']).read_text())
            self.assertIn('locus_tag=gene1', Path(manifest['dna']['gff']).read_text())
            self.assertEqual(refs.prepare_organism(folder, self.cache, 'bacteria', 2), manifest)
            self.assertEqual(builder.call_count, 2)
        self.assertEqual(len(list(folder.iterdir())), 4)

    def test_extract_noncoding_rna_when_no_rna_fasta(self):
        folder = self.bacteria()
        with patch.object(refs, 'build_index', side_effect=self.build):
            manifest = refs.prepare_organism(folder, self.cache, 'bacteria', 2)
        self.assertIsNotNone(manifest['rna'])
        self.assertIn('rRNA', Path(manifest['rna']['gff']).read_text())

    def test_duplicate_sequence_identifiers_fail(self):
        first, second = self.root / 'first.fa', self.root / 'second.fa'
        first.write_text('>duplicate\nACGT\n')
        second.write_text('>duplicate\nTGCA\n')
        with self.assertRaisesRegex(ValueError, 'duplicate FASTA'):
            refs.combine_fasta([first, second], self.root / 'combined.fa')

    def test_mismatched_gff_fails(self):
        folder = self.bacteria()
        self.write_gz(folder / 'GCF_000006765.1_genomic.gff.gz',
                      'other\ttest\tCDS\t1\t100\t.\t+\t0\tlocus_tag=gene\n')
        with self.assertRaisesRegex(ValueError, 'coordinates disagree'):
            refs.prepare_organism(folder, self.cache, 'bacteria', 2)

    def test_multiple_assemblies_fail(self):
        folder = self.bacteria()
        self.write_gz(folder / 'GCF_000006766.1_genomic.fna.gz', '>other\nACGT\n')
        with self.assertRaisesRegex(ValueError, 'Multiple genome assemblies'):
            refs.prepare_organism(folder, self.cache, 'bacteria', 2)

    def test_prebuilt_large_decoy_index_reused(self):
        folder = self.root / 'decoy'
        folder.mkdir()
        for part in ('1', '2', '3', '4', 'rev.1', 'rev.2'):
            (folder / f'Gina_Pangenome.{part}.bt2l').write_text('index')
        (folder / 'Gina_Pangenome.fasta').write_text('>decoy\nACGT\n')
        with patch.object(refs, 'build_index') as builder:
            manifest = refs.prepare_decoy(folder, self.cache, 2)
            builder.assert_not_called()
        self.assertEqual(manifest['index'], str(folder / 'Gina_Pangenome'))

    def test_host_split_rna_genbank_and_genome_gff(self):
        try:
            from Bio import SeqIO
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            from Bio.SeqFeature import SeqFeature, SimpleLocation
        except ImportError:
            self.skipTest('Biopython unavailable')
        folder = self.root / 'human'
        self.write_gz(folder / 'human_DNA' / 'GCF_000001405.40_genomic.fna.gz', '>chr1\n' + 'ACGT' * 100 + '\n')
        self.write_gz(folder / 'human_DNA' / 'GCF_000001405.40_genomic.gff.gz',
                      'chr1\ttest\tCDS\t1\t100\t.\t+\t0\tlocus_tag=gene1\n')
        for i in (1, 2):
            self.write_gz(folder / 'human_RNA' / f'human.{i}.rna.fna.gz', f'>NM_00{i}.1\n' + 'ACGT' * 50 + '\n')
            record = SeqRecord(Seq('ACGT' * 50), id=f'NM_00{i}.1', name=f'NM_00{i}', description='transcript')
            record.annotations['molecule_type'] = 'mRNA'
            record.features = [SeqFeature(SimpleLocation(0, 200), type='source', qualifiers={'mol_type': ['mRNA']}),
                               SeqFeature(SimpleLocation(20, 150), type='CDS', qualifiers={'gene': [f'gene{i}']})]
            with gzip.open(folder / 'human_RNA' / f'human.{i}.rna.gbff.gz', 'wt') as output:
                SeqIO.write(record, output, 'genbank')
        with patch.object(refs, 'build_index', side_effect=self.build):
            manifest = refs.prepare_organism(folder, self.cache, 'host', 2)
        annotation = Path(manifest['rna']['gff']).read_text()
        self.assertIn('NM_001.1\treference_folders\tCDS\t21\t150', annotation)
        self.assertIn('NM_002.1\treference_folders\tmRNA\t1\t200', annotation)

    def test_supplied_unannotated_transcript_fails(self):
        folder = self.bacteria()
        self.write_gz(folder / 'mystery.rna.fna.gz', '>unknown\nACGT\n')
        with self.assertRaisesRegex(ValueError, 'matching RNA GBFF'):
            refs.prepare_organism(folder, self.cache, 'bacteria', 2)

    def test_renamed_annotations_need_no_shared_basename(self):
        folder = self.bacteria()
        (folder / 'GCF_000006765.1_genomic.gff.gz').rename(folder / 'genes.gff.gz')
        self.write_gz(folder / 'arbitrary.rna.fna.gz', '>rna1\nACGTACGT\n')
        self.write_gz(folder / 'rna_annotations.gff.gz',
                      'rna1\ttest\trRNA\t1\t8\t.\t+\t.\tlocus_tag=rna1\n')
        with patch.object(refs, 'build_index', side_effect=self.build):
            manifest = refs.prepare_organism(folder, self.cache, 'bacteria', 2)
        self.assertIn('\trRNA\t1\t8\t', Path(manifest['rna']['gff']).read_text())

    def test_cache_cannot_be_inside_input_reference_folder(self):
        folder = self.bacteria()
        with self.assertRaisesRegex(ValueError, 'outside the input'):
            refs.prepare_organism(folder, folder / '.cache', 'bacteria', 2)
        with self.assertRaisesRegex(ValueError, 'outside the input'):
            refs.prepare_decoy(folder, folder / '.cache', 2)

    def test_partial_rna_annotations_fail_instead_of_unannotated_counts(self):
        folder = self.bacteria()
        for i in (1, 2):
            self.write_gz(folder / f'RNA.{i}.rna.fna.gz', f'>rna{i}\nACGTACGT\n')
        self.write_gz(folder / 'RNA.1.rna.gff.gz',
                      'rna1\ttest\trRNA\t1\t8\t.\t+\t.\tlocus_tag=rna1\n')
        with self.assertRaisesRegex(ValueError, 'Missing RNA annotation.*rna2'):
            refs.prepare_organism(folder, self.cache, 'bacteria', 2)

    def test_ambiguous_renamed_annotations_fail(self):
        folder = self.bacteria()
        original = folder / 'GCF_000006765.1_genomic.gff.gz'
        (folder / 'genes1.gff.gz').write_bytes(original.read_bytes())
        original.rename(folder / 'genes2.gff.gz')
        with self.assertRaisesRegex(ValueError, 'Ambiguous annotation files'):
            refs.prepare_organism(folder, self.cache, 'bacteria', 2)

    def test_gff_retains_exons_and_utrs(self):
        folder = self.bacteria()
        self.write_gz(folder / 'GCF_000006765.1_genomic.gff.gz',
                      'chromosome\ttest\tgene\t1\t300\t.\t+\t.\tID=gene;gene=gene1\n'
                      'chromosome\ttest\tmRNA\t1\t300\t.\t+\t.\tID=transcript;Parent=gene\n'
                      'chromosome\ttest\texon\t1\t100\t.\t+\t.\tID=exon1;Parent=transcript\n'
                      'chromosome\ttest\tfive_prime_UTR\t1\t20\t.\t+\t.\tID=utr;Parent=transcript\n'
                      'chromosome\ttest\tCDS\t21\t100\t.\t+\t0\tID=cds;Parent=transcript;gene=gene1\n')
        with patch.object(refs, 'build_index', side_effect=self.build):
            manifest = refs.prepare_organism(folder, self.cache, 'bacteria', 2)
        annotation = Path(manifest['dna']['gff']).read_text()
        self.assertIn('\texon\t1\t100\t', annotation)
        self.assertIn('\tfive_prime_UTR\t1\t20\t', annotation)
        self.assertIn('Parent=transcript', annotation)

    def test_genbank_only_genomic_and_rna_references(self):
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
        folder = self.root / 'host'
        dna = SeqRecord(Seq('ACGT' * 100), id='chr1', name='chr1', description='genome')
        dna.annotations['molecule_type'] = 'DNA'
        dna.features = [SeqFeature(CompoundLocation([SimpleLocation(0, 100), SimpleLocation(200, 300)]),
                                   type='mRNA', qualifiers={'gene': ['gene1']}),
                        SeqFeature(SimpleLocation(20, 90), type='CDS', qualifiers={'gene': ['gene1']})]
        rna = SeqRecord(Seq('ACGT' * 50), id='NM_001.1', name='NM_001', description='transcript')
        rna.annotations['molecule_type'] = 'mRNA'
        rna.features = [SeqFeature(SimpleLocation(0, 200), type='source', qualifiers={'mol_type': ['mRNA']}),
                        SeqFeature(SimpleLocation(20, 150), type='CDS', qualifiers={'gene': ['gene1']})]
        for name, record in [('host_DNA/genomic.gbff.gz', dna), ('host_RNA/host.1.rna.gbff.gz', rna)]:
            path = folder / name
            path.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(path, 'wt') as output:
                SeqIO.write(record, output, 'genbank')
        with patch.object(refs, 'build_index', side_effect=self.build):
            manifest = refs.prepare_organism(folder, self.cache, 'host', 2)
        genomic = Path(manifest['dna']['gff']).read_text()
        self.assertIn('\texon\t1\t100\t', genomic)
        self.assertIn('\texon\t201\t300\t', genomic)
        self.assertNotIn('\texon\t1\t300\t', genomic)
        self.assertIn('NM_001.1', Path(manifest['rna']['gff']).read_text())

    def test_renamed_genbank_annotation_does_not_duplicate_fasta_sequence(self):
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, SimpleLocation
        folder = self.root / 'organism'
        folder.mkdir()
        (folder / 'genome.fna').write_text('>chr1\n' + 'ACGT' * 100 + '\n')
        (folder / 'transcripts.rna.fna').write_text('>NM_001.1\n' + 'ACGT' * 50 + '\n')
        dna = SeqRecord(Seq('ACGT' * 100), id='chr1', name='chr1', description='genome')
        dna.annotations['molecule_type'] = 'DNA'
        dna.features = [SeqFeature(SimpleLocation(20, 90), type='CDS', qualifiers={'gene': ['gene1']})]
        rna = SeqRecord(Seq('ACGT' * 50), id='NM_001.1', name='NM_001', description='transcript')
        rna.annotations['molecule_type'] = 'mRNA'
        rna.features = [SeqFeature(SimpleLocation(0, 200), type='source', qualifiers={'mol_type': ['mRNA']}),
                        SeqFeature(SimpleLocation(20, 150), type='CDS', qualifiers={'gene': ['gene1']})]
        for filename, record in [('annotations.gbff', dna), ('transcript_annotations.gbff', rna)]:
            with (folder / filename).open('w') as output:
                SeqIO.write(record, output, 'genbank')
        with patch.object(refs, 'build_index', side_effect=self.build):
            manifest = refs.prepare_organism(folder, self.cache, 'bacteria', 2)
        self.assertEqual(Path(manifest['dna']['index'] + '.fa').read_text().count('>'), 1)
        self.assertEqual(Path(manifest['rna']['index'] + '.fa').read_text().count('>'), 1)
        self.assertIn('\tCDS\t21\t150\t', Path(manifest['rna']['gff']).read_text())

    def test_fasta_genbank_sequence_mismatch_fails_and_gbff_only_shards_retained(self):
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        fasta = self.root / 'sequence.fa'
        fasta.write_text('>id1\nACGT\n')
        records = [SeqRecord(Seq('ACGT'), id='id1'), SeqRecord(Seq('AAAA'), id='id2')]
        for record in records:
            record.annotations['molecule_type'] = 'DNA'
        gbff = self.root / 'othername.gbff'
        SeqIO.write(records, gbff, 'genbank')
        lengths, _ = refs.combine_reference([fasta, gbff], self.root / 'out.fa')
        self.assertEqual(lengths, {'id1': 4, 'id2': 4})
        records[0].seq = Seq('GGGG')
        SeqIO.write(records, gbff, 'genbank')
        with self.assertRaisesRegex(ValueError, 'sequences disagree'):
            refs.combine_reference([fasta, gbff], self.root / 'out.fa')

    def test_duplicate_genbank_ids_not_silently_deduplicated(self):
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        fasta = self.root / 'sequence.fa'
        fasta.write_text('>id1\nACGT\n')
        record = SeqRecord(Seq('ACGT'), id='id1')
        record.annotations['molecule_type'] = 'DNA'
        first, second = self.root / 'first.gbff', self.root / 'second.gbff'
        SeqIO.write(record, first, 'genbank')
        SeqIO.write(record, second, 'genbank')
        with self.assertRaisesRegex(ValueError, 'Duplicate GenBank identifier'):
            refs.combine_reference([fasta, first, second], self.root / 'out.fa')


if __name__ == '__main__':
    unittest.main()
