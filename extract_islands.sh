#!/bin/bash
set -e

INT_FASTA="$1"
FHUB_FASTA="$2"
GENOMES_DIR="$3"
OUT_DIR="$4"
THREADS="${5:-1}"
IDENTITY="${6:-90}"
EXPAND="${7:-100}"

if [ $# -lt 4 ]; then
    echo "Usage: $0 <int.fa> <fhub.fa> <genomes_dir> <out_dir> [threads] [identity] [expand]"
    exit 1
fi

mkdir -p "$OUT_DIR"/{positive_strains,islands,prokka,resfinder}

positive_list="$OUT_DIR/positive_strains/positive_genomes.list"
> "$positive_list"

for genome in "$GENOMES_DIR"/*.{fa,fasta,fna,fsa} 2>/dev/null; do
    [ -f "$genome" ] || continue
    base=$(basename "$genome" | sed 's/\.[^.]*$//')
    if blastn -query "$INT_FASTA" -db "$genome" -perc_identity "$IDENTITY" -evalue 1e-5 -outfmt 6 | head -1 | grep -q .; then
        echo "$genome" >> "$positive_list"
    fi
done

pos_count=$(wc -l < "$positive_list")
if [ $pos_count -eq 0 ]; then
    exit 0
fi

while IFS= read -r genome; do
    base=$(basename "$genome" | sed 's/\.[^.]*$//')
    if [ $(grep -c "^>" "$genome") -ne 1 ]; then
        continue
    fi
    makeblastdb -in "$genome" -dbtype nucl -parse_seqids >/dev/null 2>&1
    int_hit=$(blastn -query "$INT_FASTA" -db "$genome" -perc_identity "$IDENTITY" -evalue 1e-5 -outfmt "6 sseqid sstart send" | head -1)
    [ -z "$int_hit" ] && continue
    int_contig=$(echo "$int_hit" | cut -f1)
    int_start=$(echo "$int_hit" | cut -f2)
    int_end=$(echo "$int_hit" | cut -f3)
    [ $int_start -gt $int_end ] && { tmp=$int_start; int_start=$int_end; int_end=$tmp; }
    fhub_hit=$(blastn -query "$FHUB_FASTA" -db "$genome" -perc_identity "$IDENTITY" -evalue 1e-5 -outfmt "6 sseqid sstart send" | head -1)
    [ -z "$fhub_hit" ] && continue
    fhub_contig=$(echo "$fhub_hit" | cut -f1)
    fhub_start=$(echo "$fhub_hit" | cut -f2)
    fhub_end=$(echo "$fhub_hit" | cut -f3)
    [ $fhub_start -gt $fhub_end ] && { tmp=$fhub_start; fhub_start=$fhub_end; fhub_end=$tmp; }
    [ "$int_contig" != "$fhub_contig" ] && continue
    start=$(( (int_start < fhub_start ? int_start : fhub_start) - EXPAND ))
    end=$(( (int_end > fhub_end ? int_end : fhub_end) + EXPAND ))
    [ $start -lt 1 ] && start=1
    samtools faidx "$genome" "${int_contig}:${start}-${end}" > "$OUT_DIR/islands/${base}_island.fasta" 2>/dev/null
    [ ! -s "$OUT_DIR/islands/${base}_island.fasta" ] && continue
    prokka "$OUT_DIR/islands/${base}_island.fasta" --outdir "$OUT_DIR/prokka/${base}" --prefix "$base" --cpus "$THREADS" --force >/dev/null 2>&1
    run_resfinder.py -ifa "$OUT_DIR/islands/${base}_island.fasta" -o "$OUT_DIR/resfinder/${base}" -t "$THREADS" >/dev/null 2>&1
done < "$positive_list"
