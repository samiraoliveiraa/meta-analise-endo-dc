#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc.log
#SBATCH --error=fastqc.err
#SBATCH --ntasks=16

artigo="NOME DO ARTIGO"
input_dir="CAMINHO DOS DADOS ORIGINAIS"
output_dir="CAMINHO DO OUTPUT APÓS O FASTQC"
mkdir -p "$output_dir"

for file in "$input_dir"/.fastq; do
            fastqc "$file" -o "$output_dir"
        done
