#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_NAME=$(basename "$0")

usage() {
  cat <<EOF
Usage: $SCRIPT_NAME [RAWDIR] [GENOME_FASTA] [WORKDIR] [prefix] [suffix]

RAWDIR  - directory containing FASTQ files; searches for pairs matching:
         **/*_R1.f*q*  and  **/*_1.f*q*

GENOME_FASTA - path to a FASTA file for bwa mem
               (an index is required: bwa index genome.fa)

WORKDIR - where to store results/temporary files (default: \$PWD)

prefix - prefix (in the filename without extension) to remove if it matches
suffix - suffix (in the filename without extension) to remove if it matches

Options:
  -h, --help    Show this help message and exit.
EOF
}

if [[ "${1-}" == "-h" || "${1-}" == "--help" ]]; then
  usage
  exit 0
fi

if (( $# < 2 )); then
  echo "ERROR: at least 2 arguments are required: RAWDIR and GENOME_FASTA" >&2
  usage >&2
  exit 1
fi

if (( $# > 5 )); then
  echo "ERROR: too many arguments (max 5)" >&2
  usage >&2
  exit 1
fi

RAWDIR="$(realpath "$1")"
GENOME="$(realpath "$2")"

GNAME=$(basename -s .fa $GENOME)

WORKDIR="${3:-$PWD}"
mkdir -p -- "$WORKDIR"
WORKDIR="$(realpath "$WORKDIR")"

prefix="${4-}"
suffix="${5-}"

cd "$WORKDIR" || { echo "Failed to change directory: $WORKDIR" >&2; exit 1; }

mkdir -p -- "$WORKDIR/BAMs" "$WORKDIR/QC" "$WORKDIR/logs"

# Original stdout, sent to the terminal
exec 3>&1

RUN_TS="$(date +%F_%H-%M-%S)"
LOGDIR="$WORKDIR/logs/$RUN_TS"
QCDIR="$WORKDIR/QC/$RUN_TS"
mkdir -p -- "$LOGDIR" "$QCDIR"

MASTER_LOG="$LOGDIR/${SCRIPT_NAME}.log"

# Send all output (stdout + stderr) to the master run log
exec >"$MASTER_LOG" 2>&1

echo "Working directory: $WORKDIR" >&3
echo "LOGDIR: $LOGDIR" >&3
echo "MASTER LOG: $MASTER_LOG" >&3
echo "QCDIR: $QCDIR" >&3

TEMPDIR="$(mktemp -d -p "$WORKDIR" temp1.XXXXXX)"
trap 'rm -rf -- "$TEMPDIR"' EXIT
cd "$TEMPDIR"

shopt -s nullglob globstar

# Avoid processing the same file twice if the glob patterns overlap
declare -A seen_R1=()

n_found=0
n_processed=0
n_skipped=0

detect_ext() {
  local base="$1"
  case "$base" in
    *.fastq.gz) echo ".fastq.gz" ;;
    *.fq.gz)    echo ".fq.gz" ;;
    *.fastq)    echo ".fastq" ;;
    *.fq)       echo ".fq" ;;
    *)          echo "" ;;
  esac
}

run_sample() {
  local sample_id="$1"
  local R1="$2"
  local R2="$3"

  echo "========== Processing sample: $sample_id =========="
  echo " R1: $R1"
  echo " R2: $R2"

  echo "[$sample_id] Running fastp..."
  fastp -i "$R1" -I "$R2" \
    -o "${sample_id}.out.R1.fq.gz" -O "${sample_id}.out.R2.fq.gz" \
    --html "$QCDIR/${sample_id}.fastp_report.html" \
    --json "$QCDIR/${sample_id}.fastp_report.json"

  echo "[$sample_id] Running BWA-MEM..."
  local bam="${sample_id}_to_${GNAME}.aligned.sorted.bam"
  bwa mem -t 8 "$GENOME" \
    "${sample_id}.out.R1.fq.gz" "${sample_id}.out.R2.fq.gz" \
    | samtools sort -@ 8 -o "$bam" -

  echo "[$sample_id] Indexing BAM..."
  samtools index "$bam"

  local bai=""
  if [[ -f "$bam.bai" ]]; then
    bai="$bam.bai"
  elif [[ -f "${bam%.bam}.bai" ]]; then
    bai="${bam%.bam}.bai"
  else
    echo "ERROR: no .bai index found for $bam" >&2
    exit 1
  fi

  echo "[$sample_id] Cleaning up..."
  rm -f -- "${sample_id}.out.R1.fq.gz" "${sample_id}.out.R2.fq.gz"

  mv -- "$bam" "$bai" "$WORKDIR/BAMs/"

  echo "Finished sample: $sample_id"
}

for R1 in "$RAWDIR"/**/*_R1*.f*q* "$RAWDIR"/**/*_1*.f*q*; do
  # Ensure uniqueness by realpath
  R1_real="$(realpath "$R1" 2>/dev/null || true)"
  [[ -n "$R1_real" ]] || continue
  [[ -n "${seen_R1[$R1_real]-}" ]] && continue
  seen_R1["$R1_real"]=1

  ((++n_found))

  dir="$(dirname "$R1")"
  base="$(basename "$R1")"

  ext="$(detect_ext "$base")"
  if [[ -z "$ext" ]]; then
    echo "SKIP: unrecognized file extension: $base" >&2
    ((++n_skipped))
    continue
  fi

  stem="${base%$ext}"

  # Find R2 and the sample stem
  if [[ "$stem" == *"_R1"* ]]; then
    R2_stem="${stem/_R1/_R2}"
    sample_stem="${stem%%_R1*}"
  elif [[ "$stem" == *"_1" ]]; then
    R2_stem="${stem%_1}_2"
    sample_stem="${stem%_1}"
  elif [[ "$stem" == *"_1_"* ]]; then
    R2_stem="${stem/_1_/_2_}"
    sample_stem="${stem%%_1_*}"
  else
    R2_stem="${stem/_1/_2}"
    sample_stem="${stem%%_1*}"
  fi

  R2="$dir/$R2_stem$ext"
  if [[ ! -f "$R2" ]]; then
    echo "ERROR: no mate found for $R1 (expected $R2)" >&2
    exit 1
  fi

  # Apply prefix/suffix to the filename without extension
  sample_id="$sample_stem"

  if [[ -n "$prefix" && "$sample_id" != "$prefix"* ]]; then
    echo "SKIP: prefix does not match for '$sample_id' (expected '$prefix...')" >&2
    ((++n_skipped))
    continue
  fi
  if [[ -n "$prefix" ]]; then
    sample_id="${sample_id#"$prefix"}"
  fi

  if [[ -n "$suffix" && "$sample_id" != *"$suffix" ]]; then
    echo "SKIP: suffix does not match for '$sample_id' (expected '...$suffix')" >&2
    ((++n_skipped))
    continue
  fi
  if [[ -n "$suffix" ]]; then
    sample_id="${sample_id%"$suffix"}"
  fi

  # Per-sample log: output is duplicated both in the master log (stdout)
  # and in the sample log
  sample_log="$LOGDIR/${sample_id}.log"

  echo "Processing sample: $sample_id" >&3
  echo "Sample log: $sample_log" >&3

  # IMPORTANT: tee writes to the sample log, then the stream continues
  # into the master log (because stdout is already redirected there)
  {
    run_sample "$sample_id" "$R1" "$R2"
  } 2>&1 | tee -a "$sample_log"

  ((++n_processed))
done

# Summary
echo "----- SUMMARY -----" >&3
echo "Found R1: $n_found" >&3
echo "Processed: $n_processed" >&3
echo "Skipped: $n_skipped" >&3
echo "Master log: $MASTER_LOG" >&3

if (( n_found == 0 )); then
  echo "ERROR: no R1 files found in directory '$RAWDIR'" >&2
  exit 1
fi

if (( n_processed == 0 )); then
  echo "ERROR: no samples were processed (all were skipped due to prefix/suffix filters or format errors)" >&2
  exit 1
fi
