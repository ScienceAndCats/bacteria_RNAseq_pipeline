#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE="${1:-config.env}"
[[ -f "$CONFIG_FILE" ]] || { echo "Config file not found: $CONFIG_FILE" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG_FILE"
: "${FASTQ_DIR:=.}"
: "${FASTQ_GLOB:=*.fastq}"
: "${READ_LAYOUT:=single}"
: "${OUTPUT_DIR:=.}"
: "${SAMPLE_READS:=false}"
: "${SAMPLE_SIZE:=5000}"
: "${SAMPLE_SEED:=1}"
: "${DECOY_REFERENCE_DIR:?Set DECOY_REFERENCE_DIR to the decoy folder}"
: "${BACTERIA_REFERENCE_DIR:?Set BACTERIA_REFERENCE_DIR to the bacterial folder}"
: "${HOST_REFERENCE_DIR:=}"
: "${ADAPTER_SINGLE:=${ADAPTER_SEQUENCE:-AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC}}"
: "${ADAPTER_R1:=CTGTCTCTTATACACATCT}"
: "${ADAPTER_R2:=CTGTCTCTTATACACATCT}"
: "${MIN_READ_LENGTH:=22}"
: "${THREADS:=16}"
: "${GENE_POSITION_BINS:=100}"
: "${METAGENE_MIN_FEATURE_READS:=10}"
: "${COUNT_STRAND:=forward}"
: "${CSV_CONVERSION_SCRIPT:=bacteria_with_decoys_csvConversion.py}"
for var in THREADS MIN_READ_LENGTH GENE_POSITION_BINS METAGENE_MIN_FEATURE_READS; do
  [[ "${!var}" =~ ^[1-9][0-9]*$ ]] || { echo "$var must be a positive integer" >&2; exit 1; }
done
[[ "$READ_LAYOUT" == single || "$READ_LAYOUT" == paired ]] || { echo 'READ_LAYOUT must be single or paired' >&2; exit 1; }
case "$COUNT_STRAND" in forward) count_strand=1;; reverse) count_strand=2;; unstranded) count_strand=0;; *) echo 'COUNT_STRAND must be forward, reverse, or unstranded' >&2; exit 1;; esac
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ "$CSV_CONVERSION_SCRIPT" == /* ]] || CSV_CONVERSION_SCRIPT="$SCRIPT_DIR/$CSV_CONVERSION_SCRIPT"
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
FASTQ_DIR="$(cd "$FASTQ_DIR" && pwd)"
: "${REFERENCE_CACHE_DIR:=$OUTPUT_DIR/.reference_cache}"
reference_args=(--decoy-dir "$DECOY_REFERENCE_DIR" --bacteria-dir "$BACTERIA_REFERENCE_DIR" --cache-dir "$REFERENCE_CACHE_DIR" --threads "$THREADS" --output "$OUTPUT_DIR/reference_manifest.json")
[[ -z "$HOST_REFERENCE_DIR" ]] || reference_args+=(--host-dir "$HOST_REFERENCE_DIR")
python3 "$SCRIPT_DIR/reference_folders.py" "${reference_args[@]}"
# Values originate in the JSON manifest and are shell-quoted by shlex.quote.
reference_shell=$(python3 - "$OUTPUT_DIR/reference_manifest.json" <<'PYREF'
import json, shlex, sys
m=json.load(open(sys.argv[1])); stages=[('decoy','decoy','decoy',m['decoy'])]
for kind in ('rna','dna'):
    for organism in ('bacteria','host'):
        if m.get(organism) and m[organism].get(kind):
            stages.append((f'{organism}_{kind}',organism,kind,m[organism][kind]))
for name,values in [('stage_names',[s[0] for s in stages]),('stage_orgs',[s[1] for s in stages]),('stage_kinds',[s[2] for s in stages]),('stage_indexes',[s[3]['index'] for s in stages]),('stage_gffs',[s[3].get('gff','') for s in stages])]:
    print(name+'=('+ ' '.join(shlex.quote(str(v)) for v in values)+')')
PYREF
)
eval "$reference_shell"
shopt -s nullglob
matched_fastqs=("$FASTQ_DIR"/$FASTQ_GLOB)
if (( ${#matched_fastqs[@]} == 0 )); then
  echo "No FASTQ files found in $FASTQ_DIR matching $FASTQ_GLOB" >&2
  exit 1
fi

read1_filenames=()
read2_filenames=()
if [[ "$READ_LAYOUT" == "paired" ]]; then
  for i in "${matched_fastqs[@]}"; do
    if [[ "$i" =~ _R1(_[^/]*)?\.fastq(\.gz)?$ ]]; then
      mate="${i%_R1*}_R2${i##*_R1}"
      [[ -f "$mate" ]] || { echo "ERROR: Missing R2 mate for $i (expected $mate)" >&2; exit 1; }
      read1_filenames+=("$i")
      read2_filenames+=("$mate")
    elif [[ "$i" =~ _R2(_[^/]*)?\.fastq(\.gz)?$ ]]; then
      mate="${i%_R2*}_R1${i##*_R2}"
      [[ -f "$mate" ]] || { echo "ERROR: Missing R1 mate for $i (expected $mate)" >&2; exit 1; }
    fi
  done
  (( ${#read1_filenames[@]} > 0 )) || { echo "No paired FASTQs named *_R1[_suffix].fastq[.gz] and *_R2[_suffix].fastq[.gz] found" >&2; exit 1; }
else
  read1_filenames=("${matched_fastqs[@]}")
fi

# Randomly subsample complete FASTQ records with reservoir sampling. This keeps
# memory use bounded by SAMPLE_SIZE even for very large input files.
sample_fastq() {
  local input_file="$1"
  local output_file="$2"
  python3 - "$input_file" "$output_file" "$SAMPLE_SIZE" "$SAMPLE_SEED" <<'PY'
import gzip
import hashlib
import random
import sys

input_path, output_path, sample_size, seed = sys.argv[1:]
sample_size = int(sample_size)
rng = random.Random(f"{seed}:{hashlib.sha256(input_path.encode()).hexdigest()}")
opener = gzip.open if input_path.endswith(".gz") else open
reservoir = []

with opener(input_path, "rt") as source:
    read_count = 0
    while True:
        record = [source.readline() for _ in range(4)]
        if not record[0]:
            break
        if any(line == "" for line in record):
            raise SystemExit(f"Incomplete FASTQ record in {input_path}")
        if read_count < sample_size:
            reservoir.append(record)
        else:
            replacement = rng.randrange(read_count + 1)
            if replacement < sample_size:
                reservoir[replacement] = record
        read_count += 1

with open(output_path, "w") as destination:
    for record in reservoir:
        destination.writelines(record)
PY
}

# Sample the same record numbers from both mates so pairing is preserved.
sample_fastq_pair() {
  python3 - "$1" "$2" "$3" "$4" "$SAMPLE_SIZE" "$SAMPLE_SEED" <<'PY'
import gzip, hashlib, random, sys
r1, r2, out1, out2, size, seed = sys.argv[1:]
size = int(size)
opener = lambda path: gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")
rng = random.Random(f"{seed}:{hashlib.sha256((r1 + ':' + r2).encode()).hexdigest()}")
reservoir = []
with opener(r1) as first, opener(r2) as second:
    count = 0
    while True:
        a, b = [first.readline() for _ in range(4)], [second.readline() for _ in range(4)]
        if not a[0] and not b[0]: break
        if not a[0] or not b[0] or any(x == "" for x in a + b):
            raise SystemExit(f"Mates have unequal or incomplete FASTQ records: {r1}, {r2}")
        if count < size: reservoir.append((a, b))
        else:
            replacement = rng.randrange(count + 1)
            if replacement < size: reservoir[replacement] = (a, b)
        count += 1
with open(out1, "w") as first, open(out2, "w") as second:
    for a, b in reservoir: first.writelines(a); second.writelines(b)
PY
}

case "${SAMPLE_READS,,}" in
  true)
    [[ "$SAMPLE_SIZE" =~ ^[1-9][0-9]*$ ]] || { echo "SAMPLE_SIZE must be a positive integer" >&2; exit 1; }
    sample_dir="$OUTPUT_DIR/.bacteria_sampled_fastq"
    rm -rf "$sample_dir"
    mkdir -p "$sample_dir"
    sampled_r1=(); sampled_r2=()
    for idx in "${!read1_filenames[@]}"; do
      out1="$sample_dir/$(basename "${read1_filenames[$idx]%.gz}")"
      if [[ "$READ_LAYOUT" == "paired" ]]; then
        out2="$sample_dir/$(basename "${read2_filenames[$idx]%.gz}")"
        echo "sampling up to $SAMPLE_SIZE read pairs from $(basename "${read1_filenames[$idx]}") and $(basename "${read2_filenames[$idx]}") ..."
        sample_fastq_pair "${read1_filenames[$idx]}" "${read2_filenames[$idx]}" "$out1" "$out2"
        sampled_r2+=("$out2")
      else
        echo "sampling up to $SAMPLE_SIZE reads from $(basename "${read1_filenames[$idx]}") ..."
        sample_fastq "${read1_filenames[$idx]}" "$out1"
      fi
      sampled_r1+=("$out1")
    done
    read1_filenames=("${sampled_r1[@]}"); read2_filenames=("${sampled_r2[@]}")
    ;;
  false) ;;
  *) echo "SAMPLE_READS must be true or false" >&2; exit 1 ;;
esac

# Trim adapters while retaining paired mates together.
trim_r1=(); trim_r2=(); sample_names=()
for idx in "${!read1_filenames[@]}"; do
  filename=$(basename "${read1_filenames[$idx]}"); filename=${filename%.gz}; filename=${filename%.fastq}
  sample="$filename"
  if [[ "$READ_LAYOUT" == "paired" ]]; then sample="${filename%_R1*}${filename##*_R1}"; fi
  out1="$OUTPUT_DIR/${sample}_${MIN_READ_LENGTH}bp.trim.fastq"
  if [[ "$READ_LAYOUT" == "paired" ]]; then
    out1="$OUTPUT_DIR/${sample}_R1_${MIN_READ_LENGTH}bp.trim.fastq"
    out2="$OUTPUT_DIR/${sample}_R2_${MIN_READ_LENGTH}bp.trim.fastq"
    cutadapt -m "$MIN_READ_LENGTH" -j "$THREADS" -a "$ADAPTER_R1" -A "$ADAPTER_R2" \
      -o "$out1" -p "$out2" "${read1_filenames[$idx]}" "${read2_filenames[$idx]}" \
      > "$OUTPUT_DIR/${sample}_${MIN_READ_LENGTH}bp.cutadapt_log.txt"
    trim_r2+=("$out2")
  else
    cutadapt -m "$MIN_READ_LENGTH" -j "$THREADS" -a "$ADAPTER_SINGLE" -o "$out1" "${read1_filenames[$idx]}" \
      > "$OUTPUT_DIR/${sample}_${MIN_READ_LENGTH}bp.cutadapt_log.txt"
  fi
  trim_r1+=("$out1"); sample_names+=("$sample")
done

# Only concordant pairs are assigned at a paired stage. This makes the
# --un-conc files and per-pair counts agree, without mixed/discordant reuse.
events_file=$(mktemp "$OUTPUT_DIR/.run_events.XXXXXX")
trap 'rm -f "$events_file"' EXIT
for idx in "${!trim_r1[@]}"; do
  sample="${sample_names[$idx]}"; basename="${sample}_${MIN_READ_LENGTH}bp"
  current1="${trim_r1[$idx]}"; current2=""
  [[ "$READ_LAYOUT" != paired ]] || current2="${trim_r2[$idx]}"
  for stage in "${!stage_names[@]}"; do
    organism="${stage_orgs[$stage]}"; kind="${stage_kinds[$stage]}"
    directory="$OUTPUT_DIR/bowtie_alignments/$organism"
    [[ "$kind" != rna ]] || directory="$directory/rna"
    mkdir -p "$directory"
    prefix="$directory/${organism^^}_${basename}"
    if [[ "$organism" == decoy ]]; then
      prefix="$directory/${basename}.mapped_to_other_bugs"
      log="$prefix.bowtie2.txt"
    else
      log="$prefix.bowtie_output.txt"
    fi
    echo "mapping $sample to ${stage_names[$stage]} ..."
    if [[ "$READ_LAYOUT" == paired ]]; then
      unmapped="$prefix.unmapped_%.fastq"
      input_args=(-1 "$current1" -2 "$current2" --no-mixed --no-discordant)
      output_args=(--un-conc "$unmapped")
    else
      unmapped="$prefix.unmapped.fastq"
      input_args=(-U "$current1"); output_args=(--un "$unmapped")
    fi
    bowtie2 --end-to-end --reorder -p "$THREADS" -x "${stage_indexes[$stage]}" \
      "${input_args[@]}" -S "$prefix.sam" "${output_args[@]}" 2> "$log"
    # Bowtie2's gzip wrapper does not quote paths containing spaces. Compress
    # its plain FASTQ outputs ourselves with safe argument boundaries instead.
    if [[ "$READ_LAYOUT" == paired ]]; then
      gzip -f -- "${unmapped/\%/1}" "${unmapped/\%/2}"
    else
      gzip -f -- "$unmapped"
    fi
    unmapped="$unmapped.gz"
    categories=""; coverage=""
    if [[ "$organism" != decoy ]]; then
      categories="$prefix.categories.json"
      category_args=(--sam "$prefix.sam" --gff "${stage_gffs[$stage]}" --reference-kind "$kind" --strand "$COUNT_STRAND" --output "$categories")
      [[ "$READ_LAYOUT" != paired ]] || category_args+=(--paired)
      python3 "$SCRIPT_DIR/alignment_categories.py" "${category_args[@]}"
      samtools sort -@ "$THREADS" -o "$prefix.sorted.bam" "$prefix.sam"
      samtools index -@ "$THREADS" "$prefix.sorted.bam"
      coverage="$prefix.coverage.tsv"
      samtools coverage -o "$coverage" "$prefix.sorted.bam"
      # RNA and genomic coordinates stay separate; profile each bacterial stage.
      profile_supported=$(python3 - "${stage_gffs[$stage]}" <<'PYPROFILE'
import sys
with open(sys.argv[1]) as handle:
    types={line.split('\t')[2].lower() for line in handle if line.strip() and not line.startswith('#')}
print('yes' if any(t=='cds' or ('rna' in t and t not in {'mrna','transcript','primary_transcript'}) for t in types) else 'no')
PYPROFILE
)
      if [[ "$organism" == bacteria && "$profile_supported" == yes ]]; then
        profile_args=(--bam "$prefix.sorted.bam" --gff "${stage_gffs[$stage]}" --gene-id locus_tag --bins "$GENE_POSITION_BINS" --min-feature-reads "$METAGENE_MIN_FEATURE_READS" --output-prefix "$prefix")
        [[ "$READ_LAYOUT" != paired ]] || profile_args+=(--paired)
        python3 "$SCRIPT_DIR/gene_position_profile.py" "${profile_args[@]}"
      fi
    fi
    # Persist explicit current-run paths, rather than discovering stale outputs.
    python3 - "$events_file" "$sample" "$OUTPUT_DIR/${basename}.cutadapt_log.txt" "${stage_names[$stage]}" "$organism" "$kind" "$log" "$categories" "$coverage" <<'PYEVENT'
import json, sys
path,sample,cutadapt,name,organism,kind,log,categories,coverage=sys.argv[1:]
with open(path,'a') as out:
    out.write(json.dumps({'sample':sample,'cutadapt_log':cutadapt,'stage':{'name':name,'organism':organism,'kind':kind,'log':log,'categories':categories or None,'coverage':coverage or None}})+'\n')
PYEVENT
    if [[ "$READ_LAYOUT" == paired ]]; then current1="${unmapped/\%/1}"; current2="${unmapped/\%/2}"; else current1="$unmapped"; fi
  done
  leftover_dir="$OUTPUT_DIR/bowtie_alignments/leftover"; mkdir -p "$leftover_dir"
  if [[ "$READ_LAYOUT" == paired ]]; then
    cp -- "$current1" "$leftover_dir/${basename}_leftover_1.fastq.gz"
    cp -- "$current2" "$leftover_dir/${basename}_leftover_2.fastq.gz"
  else
    cp -- "$current1" "$leftover_dir/${basename}_leftover.fastq.gz"
  fi
done

# Retain gene-level featureCounts tables for both coordinate systems. The
# mutually exclusive read disposition comes from alignment_categories.py.
featurecounts_pair_args=()
[[ "$READ_LAYOUT" != paired ]] || featurecounts_pair_args=(-p --countReadPairs -B -C)
for stage in "${!stage_names[@]}"; do
  organism="${stage_orgs[$stage]}"; kind="${stage_kinds[$stage]}"
  [[ "$organism" != decoy ]] || continue
  directory="$OUTPUT_DIR/bowtie_alignments/$organism"
  [[ "$kind" != rna ]] || directory="$directory/rna"
  sam_files=()
  for sample in "${sample_names[@]}"; do sam_files+=("$directory/${organism^^}_${sample}_${MIN_READ_LENGTH}bp.sam"); done
  count_dir="$directory"
  [[ "$organism" != bacteria || "$kind" != dna ]] || count_dir="$OUTPUT_DIR"
  for feature in CDS rRNA; do
    # featureCounts errors when the requested type is absent; no table is needed
    # for a type absent from this reference's normalized annotation.
    if awk -F '\t' -v type="$feature" '!/^#/ && $3==type {found=1; exit} END {exit !found}' "${stage_gffs[$stage]}"; then
      count_output="$count_dir/featurecounts_${organism^^}_${feature}.txt"
      [[ "$feature" != CDS ]] || count_output="$count_dir/featurecounts_${organism^^}_summary.txt"
      featureCounts -T "$THREADS" -a "${stage_gffs[$stage]}" -s "$count_strand" "${featurecounts_pair_args[@]}" -g locus_tag -t "$feature" \
        -o "$count_output" "${sam_files[@]}"
      if [[ "$feature" == CDS ]]; then sed 's/\t/,/g' "$count_output" > "${count_output%.txt}.csv"; fi
    fi
  done
done
python3 - "$events_file" "$OUTPUT_DIR" "$READ_LAYOUT" "${sample_names[@]}" <<'PYMANIFEST'
import json, sys
from pathlib import Path
path,directory,layout,*names=sys.argv[1:]; out=Path(directory)
samples={}
for line in Path(path).read_text().splitlines():
    e=json.loads(line); sample=samples.setdefault(e['sample'],{'sample':e['sample'],'cutadapt_log':e['cutadapt_log'],'stages':[]})
    sample['stages'].append(e['stage'])
for name,sample in samples.items():
    stem=Path(sample['cutadapt_log']).name.removesuffix('.cutadapt_log.txt')
    suffixes=['_1','_2'] if layout=='paired' else ['']
    sample['leftover']=[str(out/'bowtie_alignments/leftover'/f'{stem}_leftover{s}.fastq.gz') for s in suffixes]
m={'read_layout':layout,'samples':[samples[n] for n in names]}
tmp=out/'run_manifest.json.tmp'; tmp.write_text(json.dumps(m,indent=2)+'\n'); tmp.replace(out/'run_manifest.json')
PYMANIFEST
(
  cd "$OUTPUT_DIR"
  python3 "$CSV_CONVERSION_SCRIPT"
  python3 "$SCRIPT_DIR/read_breakdown_csv.py"
)
