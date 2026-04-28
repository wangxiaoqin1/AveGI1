#!/bin/bash
# step1_screen_positive.sh - 筛选AveGI1阳性菌株
# Usage: ./step1_screen_positive.sh <int.fa> <genomes_dir> <output_list.txt> [identity]

INT_FASTA=$1
GENOMES_DIR=$2
OUT_LIST=$3
IDENTITY=${4:-90}

if [ $# -lt 3 ]; then
    echo "Usage: $0 <int.fa> <genomes_dir> <output_list.txt> [identity]"
    exit 1
fi

> "$OUT_LIST"

for genome in "$GENOMES_DIR"/*.{fa,fasta,fna,fsa} 2>/dev/null; do
    [ -f "$genome" ] || continue
    base=$(basename "$genome" | sed 's/\.[^.]*$//')
    if blastn -query "$INT_FASTA" -db "$genome" -perc_identity "$IDENTITY" -evalue 1e-5 -outfmt 6 | head -1 | grep -q .; then
        echo "$genome" >> "$OUT_LIST"
        echo "POSITIVE: $base"
    else
        echo "negative: $base"
    fi
done

count=$(wc -l < "$OUT_LIST")
echo "Positive strains: $count" >&2
