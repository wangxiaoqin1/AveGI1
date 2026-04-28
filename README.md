# Genomic Island Extraction Pipeline

## Overview

This pipeline identifies AveGI1‑positive strains based on the AveGI1 integrase gene (`int`). For each positive strain, it extracts the complete genomic island located between the integrase gene (`int`) and the `fhub` gene using the final complete genome (single‑contig finished assembly). The extracted island is then annotated with Prokka and screened for antimicrobial resistance genes using ResFinder.

## Dependencies

- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [SAMtools](http://www.htslib.org/download/)
- [Prokka](https://github.com/tseemann/prokka)
- [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder/) (local version)

## Installation

```bash
# Create a conda environment (recommended)
conda create -n island_pipeline blast samtools prokka
conda activate island_pipeline

# Install ResFinder locally
git clone https://bitbucket.org/genomicepidemiology/resfinder.git
cd resfinder
python setup.py install
