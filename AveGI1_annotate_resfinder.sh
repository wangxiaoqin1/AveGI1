#!/bin/bash
# step3_annotate_resfinder.sh - Prokka and ResFinder on extracted islands
# Usage: ./step3_annotate_resfinder.sh <islands_dir> <prokka_outdir> <resfinder_outdir> [threads]

ISLANDS_DIR=$1
PROKKA_OUT=$2
RESFINDER_OUT=$3
THREADS=${4:-1}

if [ $# -lt 3 ]; then
    echo "Usage: $0 <islands_dir> <prokka_outdir> <resfinder_outdir> [threads]"
    exit 1
fi

mkdir -p "$PROKKA_OUT" "$RESFINDER_OUT"

for island_fasta in "$ISLANDS_DIR"/*_island.fasta; do
    [ -f "$island_fasta" ] || continue
    base=$(basename "$island_fasta" | sed 's/_island\.fasta$//')
    echo "Processing: $base"

    # Prokka
    prokka "$island_fasta" --outdir "$PROKKA_OUT/$base" --prefix "$base" --cpus "$THREADS" --force >/dev/null 2>&1
    if [ -d "$PROKKA_OUT/$base" ]; then
        echo "  Prokka done"
    else
        echo "  Prokka failed"
    fi

    # ResFinder
    run_resfinder.py -ifa "$island_fasta" -o "$RESFINDER_OUT/$base" -t "$THREADS" >/dev/null 2>&1
    if [ -f "$RESFINDER_OUT/$base/ResFinder_results_tab.txt" ]; then
        echo "  ResFinder results:"
        cat "$RESFINDER_OUT/$base/ResFinder_results_tab.txt" | sed 's/^/    /'
    else
        echo "  ResFinder failed or no resistance genes"
    fi
done
