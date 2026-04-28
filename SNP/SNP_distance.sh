#!/bin/bash
set -euo pipefail

GENOMES_DIR=${1:-./genomes}
ANNOT_DIR=${2:-./annotations}
PANAROO_DIR=${3:-./panaroo_results}

THREADS_PROKKA=8
THREADS_PANAROO=32
CORE_THRESHOLD=0.95
BOOTSTRAP=1000
NMAX=5000

mkdir -p "$ANNOT_DIR" "$PANAROO_DIR"

if command -v parallel >/dev/null 2>&1; then
    parallel -j "$THREADS_PROKKA" 'base=$(basename {} .fasta); prokka --outdir "'"$ANNOT_DIR"'/$base" --prefix "$base" --force {}' ::: "$GENOMES_DIR"/*.fasta
else
    for genome in "$GENOMES_DIR"/*.fasta; do
        base=$(basename "$genome" .fasta)
        prokka --outdir "$ANNOT_DIR/$base" --prefix "$base" --force "$genome"
    done
fi

gff_files=()
for d in "$ANNOT_DIR"/*; do
    [ -d "$d" ] && gff_files+=("$d/$(basename "$d").gff")
done
panaroo -i "${gff_files[@]}" -o "$PANAROO_DIR" --clean-mode moderate --alignment core --core_threshold "$CORE_THRESHOLD" --threads "$THREADS_PANAROO"

snp-sites -c -o core_snps.aln "$PANAROO_DIR/core_gene_alignment.aln"
snp-dists core_snps.aln > snp_distance_matrix.tsv
iqtree2 --safe -s core_snps.aln -m MFP -B "$BOOTSTRAP" -nm "$NMAX" -T AUTO --prefix my_core_genome_tree
