#!/bin/bash
#SBATCH --job-name=diamond_anotacao
#SBATCH --output=diamond_anotacao.log
#SBATCH --error=diamond_anotacao.err
#SBATCH --ntasks=16

gene_catalog="caminho do catálogo de genes predito"
output_dir="caminho de output para a anotação funcional"
mkdir -p "$output_dir"

eggnog_db="/caminho/database/eggnog.dmnd"

# anotação funcional (eggnog)
diamond blastx \
  --query "$gene_catalog" \
  --db "$eggnog_db" \
  --out "$output_dir/gene_eggnog.tsv" \
  --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
  --max-target-seqs 1 \
  --evalue 1e-5 \
  --threads 16

