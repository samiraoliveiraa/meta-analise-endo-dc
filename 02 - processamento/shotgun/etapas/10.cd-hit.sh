#!/bin/bash
#SBATCH --job-name=cdhit
#SBATCH --output=cdhit.log
#SBATCH --error=cdhit.err
#SBATCH --ntasks=16

# diretório onde estão os .ffn preditos
ffn_dir="CAMINHO ONDE ESTÃO OS ARQUIVOS .ffn"
out_dir="CAMINHO DE OUTPUT PARA A GERAÇÃO DO CATÁLOGO DE GENES"
mkdir -p "$out_dir"

# concatena todos os genes preditos em um único arquivo
cat "$ffn_dir"/*_genes_ge100.ffn > "$out_dir/genes_concatenados.ffn"

# roda cd-hit-est (para sequências de DNA)
cd-hit-est -i "$out_dir/genes_concatenados.ffn" \
           -o "$out_dir/gene_catalog.ffn" \
           -c 0.95 -aS 0.9 -G 0 -g 1 -d 0 \
           -T 16 -M 8000

