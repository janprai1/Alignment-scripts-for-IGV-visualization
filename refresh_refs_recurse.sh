#!/usr/bin/env bash
set -euo pipefail

DIR="${1:-}"

if [[ -z "$DIR" ]]; then
    echo "Usage: $0 /path/to/folder"
    exit 1
fi

if [[ ! -d "$DIR" ]]; then
    echo "Error: directory does not exist: $DIR"
    exit 1
fi

found=false

while IFS= read -r -d '' fa; do
    found=true
    echo "Checking: $fa"

    if [[ ! -f "${fa}.fai" ]]; then
        echo "  Missing ${fa}.fai — running samtools faidx"
        samtools faidx "$fa"
    else
        echo "  ${fa}.fai exists"
    fi

    if [[ ! -f "${fa}.bwt" ]]; then
        echo "  Missing ${fa}.bwt — running bwa index"
        bwa index "$fa"
    else
        echo "  ${fa}.bwt exists"
    fi
done < <(
    find "$DIR" -type f \( -name '*.fa' -o -name '*.fasta' \) -print0
)

if [[ "$found" == false ]]; then
    echo "No .fa or .fasta files found in $DIR or its subfolders"
    exit 0
fi

echo "Done."