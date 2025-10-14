#!/bin/bash
#SBATCH --job-name=soapaligner
#SBATCH --output=soapaligner.log
#SBATCH --error=soapaligner.err
#SBATCH --ntasks=16

input_dir="CAMINHO DOS DADOS TRIMMADOS"
output_dir="CAMINHO PARA O OUTPUT DO ALINHAMENTO COM O GENOMA HUMANO"
mkdir -p "$output_dir"

# Aqui, apenas os arquivos pareados serão utilizados
for r1 in "$input_dir"/*_1P.fastq; do
    r2="${r1/_1P.fastq/_2P.fastq}"
    base_name=$(basename "$r1" "_1P.fastq")

    soap -D /caminho/do/genoma_humano/indexado/Homo_sapiens.GRCh38.dna.primary_assembly.fa.index \
         -a "$r1" \
         -b "$r2" \
         -o "$output_dir/${base_name}_aligned.soap" \
	 -u "$output_dir/${base_name}_unmapped.fastq" \
	 -2 "$output_dir/${base_name}_unmapped_2.fastq" \
         -s 50 -l 30 -v 5 -m 150 -x 400 -p 16 
	 
done

