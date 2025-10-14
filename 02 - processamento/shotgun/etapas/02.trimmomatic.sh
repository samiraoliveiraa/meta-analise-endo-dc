#!/bin/bash
#SBATCH -n 16
#SBATCH --job-name=trimmomatic
#SBATCH --output=trimmomatic.log
#SBATCH --error=trimmomatic.err
#SBATCH --ntasks=16

artigo="NOME DO ARTIGO"
input_dir="CAMINHO DOS DADOS ORIGINAIS"
output_dir="CAMINHO PARA O OUTPUT APÓS A TRIMMAGEM"
mkdir -p $output_dir

for file in $input_dir/*.fastq; do
    base_name=$(basename "$file" _1.fastq)

    forward="${input_dir}/${base_name}_1.fastq"
    reverse="${input_dir}/${base_name}_2.fastq"

    base_out="${output_dir}/${base_name}.fastq"

    trimmomatic PE -threads 16 \
        "$forward" "$reverse" -baseout "$base_out" \
        ILLUMINACLIP:/caminho/dos/adaptadores/TruSeq3-PE.fa:2:30:10 \
        LEADING:25 TRAILING:25 SLIDINGWINDOW:4:20 MINLEN:75
done
