#!/bin/bash
# step2_extract_island.sh - 提取基因岛序列
# Usage: ./step2_extract_island.sh <positive_list.txt> <int.fa> <fhub.fa> <islands_dir> [identity] [expand]

POSITIVE_LIST=$1
INT_FASTA=$2
FHUB_FASTA=$3
OUT_DIR=$4
IDENTITY=${5:-90}
EXPAND=${6:-100}

if [ $# -lt 4 ]; then
    echo "Usage: $0 <positive_list.txt> <int.fa> <fhub.fa> <islands_dir> [identity] [expand]"
    exit 1
fi

mkdir -p "$OUT_DIR"

while IFS= read -r genome; do
    [ -f "$genome" ] || continue
    base=$(basename "$genome" | sed 's/\.[^.]*$//')
    echo ">> $base"

    # Check single contig (finished genome)
    if [ $(grep -c "^>" "$genome") -ne 1 ]; then
        echo "  Not a finished genome, skip"
        continue
    fi

    # Index if needed
    if [ ! -f "${genome}.fai" ]; then
        samtools faidx "$genome"
    fi

    # Locate int
    int_hit=$(blastn -query "$INT_FASTA" -db "$genome" -perc_identity "$IDENTITY" -evalue 1e-5 -outfmt "6 sseqid sstart send" | head -1)
    [ -z "$int_hit" ] && { echo "  int not found"; continue; }
    int_contig=$(echo "$int_hit" | cut -f1)
    int_start=$(echo "$int_hit" | cut -f2)
    int_end=$(echo "$int_hit" | cut -f3)
    if [ $int_start -gt $int_end ]; then
        tmp=$int_start; int_start=$int_end; int_end=$tmp
    fi

    # Locate fhub
    fhub_hit=$(blastn -query "$FHUB_FASTA" -db "$genome" -perc_identity "$IDENTITY" -evalue 1e-5 -outfmt "6 sseqid sstart send" | head -1)
    [ -z "$fhub_hit" ] && { echo "  fhub not found"; continue; }
    fhub_contig=$(echo "$fhub_hit" | cut -f1)
    fhub_start=$(echo "$fhub_hit" | cut -f2)
    fhub_end=$(echo "$fhub_hit" | cut -f3)
    if [ $fhub_start -gt $fhub_end ]; then
        tmp=$fhub_start; fhub_start=$fhub_end; fhub_end=$tmp
    fi

    # Must be same contig
    [ "$int_contig" != "$fhub_contig" ] && { echo "  int and fhub on different contigs"; continue; }

    # Compute region
    start=$(( (int_start < fhub_start ? int_start : fhub_start) - EXPAND ))
    end=$(( (int_end > fhub_end ? int_end : fhub_end) + EXPAND ))
    [ $start -lt 1 ] && start=1

    echo "  Extracting $int_contig:$start-$end"
    samtools faidx "$genome" "${int_contig}:${start}-${end}" > "$OUT_DIR/${base}_island.fasta" 2>/dev/null
    if [ -s "$OUT_DIR/${base}_island.fasta" ]; then
        echo "  Saved: $OUT_DIR/${base}_island.fasta"
    else
        echo "  Extraction failed"
    fi
done < "$POSITIVE_LIST"
