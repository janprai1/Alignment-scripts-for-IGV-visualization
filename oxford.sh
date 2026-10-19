#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_NAME=$(basename "$0")

on_error() {
  local exit_code=$?
  local line_no=${1:-unknown}
  local cmd=${2:-unknown}
  local src=${3:-$SCRIPT_NAME}

  {
    echo "ERROR: команда завершилась с кодом $exit_code"
    echo "  файл:   $src"
    echo "  строка: $line_no"
    echo "  команда: $cmd"
  } >&2

  exit "$exit_code"
}

trap 'on_error "${LINENO}" "${BASH_COMMAND}" "${BASH_SOURCE[0]}"' ERR

usage() {
  cat <<EOF
Usage: $SCRIPT_NAME [RAWDIR] [GENOME_FASTA] [WORKDIR] [prefix] [suffix]

RAWDIR  - папка с fastq; ищем fastq файлы

GENOME_FASTA - путь к FASTA для minimap2 (индекс не нужен, : minimap2 index genome.fa)

WORKDIR - куда складывать результаты/временные файлы (по умолчанию \$PWD)

prefix - префикс (в имени без расширения), который надо отрезать, если совпал
suffix - суффикс (в имени без расширения), который надо отрезать, если совпал

Опции:
  -h, --help    Показать справку и выйти.
EOF
}

if [[ "${1-}" == "-h" || "${1-}" == "--help" ]]; then
  usage
  exit 0
fi

if (( $# < 2 )); then
  echo "ERROR: нужно минимум 2 аргумента: RAWDIR и GENOME_FASTA" >&2
  usage >&2
  exit 1
fi

if (( $# > 5 )); then
  echo "ERROR: слишком много аргументов (max 5)" >&2
  usage >&2
  exit 1
fi

command -v minimap2 >/dev/null || { echo "ERROR: minimap2 не найден" >&2; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools не найден" >&2; exit 1; }
command -v realpath >/dev/null || { echo "ERROR: realpath не найден" >&2; exit 1; }
[[ -d "$1" ]] || { echo "ERROR: RAWDIR не существует: $1" >&2; exit 1; }
[[ -e "$2" ]] || { echo "ERROR: GENOME_FASTA не существует: $2" >&2; exit 1; }

RAWDIR="$(realpath "$1")"
GENOME="$(realpath "$2")"

WORKDIR="${3:-$PWD}"
mkdir -p -- "$WORKDIR"
WORKDIR="$(realpath "$WORKDIR")"

prefix="${4-}"
suffix="${5-}"

cd "$WORKDIR" || { echo "Не удалось перейти в директорию: $WORKDIR" >&2; exit 1; }

mkdir -p -- "$WORKDIR/BAMs" "$WORKDIR/QC" "$WORKDIR/logs"

# Оригинальный stdout в терминал
exec 3>&1

RUN_TS="$(date +%F_%H-%M-%S)"
LOGDIR="$WORKDIR/logs/$RUN_TS"
QCDIR="$WORKDIR/QC/$RUN_TS"
mkdir -p -- "$LOGDIR" "$QCDIR"

MASTER_LOG="$LOGDIR/${SCRIPT_NAME}.log"

# Всё (stdout+stderr) в общий лог запуска
exec >"$MASTER_LOG" 2>&1

echo "Рабочая папка: $WORKDIR" >&3
echo "LOGDIR: $LOGDIR" >&3
echo "MASTER LOG: $MASTER_LOG" >&3
echo "QCDIR: $QCDIR" >&3

TEMPDIR="$(mktemp -d -p "$WORKDIR" temp1.XXXXXX)"
trap 'rm -rf -- "$TEMPDIR"' EXIT
cd "$TEMPDIR"

shopt -s nullglob globstar

# чтобы не обработать один и тот же файл дважды, если маски пересекутся
declare -A seen_FASTQ=()

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
  local FASTQ="$2"

  echo "========== Processing sample: $sample_id =========="
  echo " FASTQ: $FASTQ"

  echo "[$sample_id] Running minimap2..."
  
  local bam="${sample_id}.sorted.bam"
  minimap2 -t 8 -a -x map-ont -Y --MD --cs "$GENOME" "$FASTQ" \
  | samtools sort -O BAM -o "$bam" -
  
  echo "[$sample_id] Indexing BAM..."
  samtools index "$bam"

  local bai=""
  if [[ -f "$bam.bai" ]]; then
    bai="$bam.bai"
  elif [[ -f "${bam%.bam}.bai" ]]; then
    bai="${bam%.bam}.bai"
  else
    echo "ERROR: не найден bai для $bam" >&2
    exit 1
  fi

  echo "[$sample_id] Cleaning up..."

  mv -- "$bam" "$bai" "$WORKDIR/BAMs/"

  echo "Finished sample: $sample_id"
}

for FASTQ in \
  "$RAWDIR"/**/*.fastq.gz \
  "$RAWDIR"/**/*.fq.gz \
  "$RAWDIR"/**/*.fastq \
  "$RAWDIR"/**/*.fq
do
  echo "обнаружен файл" $FASTQ >&3
  # уникальность по realpath
  FASTQ_real="$(realpath "$FASTQ" 2>/dev/null || true)"
  [[ -n "$FASTQ_real" ]] || continue
  [[ -n "${seen_FASTQ[$FASTQ_real]-}" ]] && continue
  seen_FASTQ["$FASTQ_real"]=1

  ((++n_found))

  dir="$(dirname "$FASTQ")"

  # Применяем prefix/suffix к имени без расширения
  base="$(basename "$FASTQ")"
  
  ext="$(detect_ext "$base")"
  if [[ -z "$ext" ]]; then
    echo "SKIP: не понимаю расширение: $base" >&2
    ((++n_skipped))
    continue
  fi

  sample_id="${base%$ext}"
  
  if [[ -n "$prefix" && "$sample_id" != "$prefix"* ]]; then
    echo "SKIP: prefix не совпал для '$sample_id' (ожидался '$prefix...')" >&2
    ((++n_skipped))
    continue
  fi
  if [[ -n "$prefix" ]]; then
    sample_id="${sample_id#"$prefix"}"
  fi

  if [[ -n "$suffix" && "$sample_id" != *"$suffix" ]]; then
    echo "SKIP: suffix не совпал для '$sample_id' (ожидался '...$suffix')" >&2
    ((++n_skipped))
    continue
  fi
  if [[ -n "$suffix" ]]; then
    sample_id="${sample_id%"$suffix"}"
  fi

  # Лог на сэмпл: вывод дублируется и в общий лог (stdout), и в sample log
  sample_log="$LOGDIR/${sample_id}.log"

  echo "Processing sample: $sample_id" >&3
  echo "Sample log: $sample_log" >&3

  # ВАЖНО: tee пишет в sample log, а дальше поток идёт в общий лог (т.к. stdout уже туда)
  {
    run_sample "$sample_id" "$FASTQ"
  } 2>&1 | tee -a "$sample_log"

  ((++n_processed))
done

# Итоги
echo "----- SUMMARY -----" >&3
echo "Found FASTQ: $n_found" >&3
echo "Processed: $n_processed" >&3
echo "Skipped: $n_skipped" >&3
echo "Master log: $MASTER_LOG" >&3

if (( n_found == 0 )); then
  echo "ERROR: в каталоге '$RAWDIR' не найдено ни одного FASTQ файла" >&2
  exit 1
fi

if (( n_processed == 0 )); then
  echo "ERROR: не обработано ни одного sample (все пропущены из-за prefix/suffix или ошибок формата)" >&2
  exit 1
fi
