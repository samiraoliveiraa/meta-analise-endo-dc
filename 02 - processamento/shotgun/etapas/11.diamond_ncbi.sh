#!/bin/bash
#SBATCH --job-name=diamond_anotacao
#SBATCH --output=diamond_anotacao.log
#SBATCH --error=diamond_anotacao.err
#SBATCH --ntasks=16

gene_catalog="caminho do catálogo de genes preditos"
output_dir="caminho do output após alinhamento do catálogo contra o banco de dados"
mkdir -p "$output_dir"

ncbi_nr_db="/caminho/database/reference_nr.dmnd"

# anotação taxonômica (ncbi-nr)
diamond blastx \
  --query "$gene_catalog" \
  --db "$ncbi_nr_db" \
  --out "$output_dir/gene_ncbi.tsv" \
  --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
  --max-target-seqs 1 \
  --evalue 1e-5 \
  --threads 16

