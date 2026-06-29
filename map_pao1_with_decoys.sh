#load modules
module load anaconda3
conda activate MappingTools

#cutadapt (Add concatenate later)
fastq_filenames=($(ls *.fastq))
for i in "${fastq_filenames[@]}"
do
  i_basename=$(basename $i .fastq)
  echo "trimming adapters from $i_basename ..."
  cutadapt -m 22 -j 40 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o "$i_basename"_22bp.trim.fastq "$i_basename".fastq > "$i_basename"_22bp.cutadapt_log.txt
done

#mapping to PAO1 decoys
trim_filenames=($(ls *.trim.fastq))
for i in "${trim_filenames[@]}"
do
  i_basename=$(basename $i .trim.fastq)
  echo "map to PA decoys $i_basename ..."
  bowtie2 -p 24 -x ~/r-mwhiteley3-0/referencegenomes/Gina_Pangenome -U "$i_basename".trim.fastq -S "$i_basename".mapped_to_other_bugs.sam --un-gz "$i_basename"_unmapped_to_other_bugs.fastq > "$i_basename".mapped_to_other_bugs.bowtie2.txt
done

#mapping with bowtie2 PAO1
PA_trim_filenames=($(ls *unmapped_to_other_bugs.fastq))
for i in "${PA_trim_filenames[@]}"
do
  i_basename=$(basename $i .fastq)
  echo "map to pseudomonas $i_basename ...."
  bowtie2 --end-to-end -p 16 -x ~/r-mwhiteley3-0/referencegenomes/Pseudomonas_aeruginosa_PAO1_107 -q -U "$i_basename".fastq -S PAO1_"$i_basename".sam 2>> PAO1_"$i_basename".bowtie_output.txt
done

PAO1_sam_filenames=($(ls PAO1*.sam))


#Create featurecounts summary file
#Frome Gina's paper: featureCounts v2.0.1 was used to assign mapped reads to PAO1 genes with the flags -s 1 (stranded) and -O (allowMultiOverlap) so that each read was assigned to a single locus or to neighboring genes
featureCounts -T 24 -a ~/r-mwhiteley3-0/referencegenomes/Pseudomonas_aeruginosa_PAO1_107.gff3 -O -s 1 -g locus -t CDS -o featurecounts_PAO1_summary.txt PAO1*.sam
sed 's/\t/,/g' featurecounts_PAO1_summary.txt > featurecounts_PAO1_summary.csv
#mv featurecounts_PAO1_summary.txt.summary featurecounts_PAO1_summary_FILE_MAY_DIFFER_FROM_ACTUAL_COUNTS_SUM.summary

### Calculate coverage with respect to PAO1 with samtools ###
for i in "${PAO1_sam_filenames[@]}"
do
    i_basename=$(basename $i .sam)
    echo "analyzing $i_basename ..."
    samtools view -@ 16 -bS "$i_basename".sam > "$i_basename".bam
    samtools sort -@ 16 -o "$i_basename".sorted.bam "$i_basename".bam
    samtools coverage -mA "$i_basename".sorted.bam
    samtools coverage -m -o "$i_basename"_coverage.txt "$i_basename".sorted.bam
done

### Runs python script to convert results to csv file ###
python3 ~/r-mwhiteley3-0/scripts/pao1_with_decoys_csvConversion.py
