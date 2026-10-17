#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_NAME=$(basename "$0")

usage() {
  cat <<EOF
Usage: $SCRIPT_NAME [RAWDIR] [WORKDIR] [prefix] [suffix]

RAWDIR  - папка с fastq; ищем пары по:
         **/*_R1.f*q*  и  **/*_1.f*q*

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

if (( $# < 1 )); then
  echo "ERROR: нужно минимум 1 аргумент: RAWDIR " >&2
  usage >&2
  exit 1
fi

if (( $# > 4 )); then
  echo "ERROR: слишком много аргументов (max 4)" >&2
  usage >&2
  exit 1
fi

RAWDIR="$(realpath "$1")"

WORKDIR="${2:-$PWD}"
mkdir -p -- "$WORKDIR"
WORKDIR="$(realpath "$WORKDIR")"

prefix="${3-}"
suffix="${4-}"

cd "$WORKDIR" || { echo "Не удалось перейти в директорию: $WORKDIR" >&2; exit 1; }

mkdir -p -- "$WORKDIR/scaffolds" "$WORKDIR/contigs" "$WORKDIR/fastg" "$WORKDIR/QC" "$WORKDIR/logs"

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
# trap 'rm -rf -- "$TEMPDIR"' EXIT
cd "$TEMPDIR"

shopt -s nullglob globstar

# чтобы не обработать один и тот же файл дважды, если маски пересекутся
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

  echo "[$sample_id] Running SPAdes..."
  
  spades.py --plasmid -1 "${sample_id}.out.R1.fq.gz" -2 "${sample_id}.out.R2.fq.gz" -o "${sample_id}"

  if [ -f ./${sample_id}/assembly_graph.fastg ]; then
    mv ./${sample_id}/assembly_graph.fastg $WORKDIR/fastg/${sample_id}_graph.fastg
  else
    echo "WARN: нет файла ./${sample_id}/assembly_graph.fastg — пропускаю" >&2
  fi

  if [ -f ./${sample_id}/contigs.fasta ]; then
    mv ./${sample_id}/contigs.fasta $WORKDIR/contigs/${sample_id}_contigs.fasta
  else
    echo "WARN: нет файла ./${sample_id}/contigs.fasta — пропускаю" >&2
  fi

  if [ -f ./${sample_id}/scaffolds.fasta ]; then
    mv ./${sample_id}/scaffolds.fasta $WORKDIR/scaffolds/${sample_id}_scaffolds.fasta
  else
    echo "WARN: нет файла ./${sample_id}/scaffolds.fasta — пропускаю" >&2
  fi
  
  }

for R1 in "$RAWDIR"/**/*_R1*.f*q* "$RAWDIR"/**/*_1*.f*q*; do
  # уникальность по realpath
  R1_real="$(realpath "$R1" 2>/dev/null || true)"
  [[ -n "$R1_real" ]] || continue
  [[ -n "${seen_R1[$R1_real]-}" ]] && continue
  seen_R1["$R1_real"]=1

  ((++n_found))

  dir="$(dirname "$R1")"
  base="$(basename "$R1")"

  ext="$(detect_ext "$base")"
  if [[ -z "$ext" ]]; then
    echo "SKIP: не понимаю расширение: $base" >&2
    ((++n_skipped))
    continue
  fi

  stem="${base%$ext}"

  # Находим R2 и sample_stem
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
    echo "ERROR: нет пары для $R1 (ожидался $R2)" >&2
    exit 1
  fi

  # Применяем prefix/suffix к имени без расширения
  sample_id="$sample_stem"

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
    run_sample "$sample_id" "$R1" "$R2"
  } 2>&1 | tee -a "$sample_log"

  ((++n_processed))
done

# Итоги
echo "----- SUMMARY -----" >&3
echo "Found R1: $n_found" >&3
echo "Processed: $n_processed" >&3
echo "Skipped: $n_skipped" >&3
echo "Master log: $MASTER_LOG" >&3

if (( n_found == 0 )); then
  echo "ERROR: в каталоге '$RAWDIR' не найдено ни одного R1 файла" >&2
  exit 1
fi

if (( n_processed == 0 )); then
  echo "ERROR: не обработано ни одного sample (все пропущены из-за prefix/suffix или ошибок формата)" >&2
  exit 1
fi
