#!/bin/bash
#SBATCH --job-name=soapdenovo
#SBATCH --output=soapdenovo.log
#SBATCH --error=soapdenovo.err
#SBATCH --ntasks=16

input_dir="CAMINHO DOS ARQUIVOS DE CONFIGURAÇÃO"
output_dir="CAMINHO DO OUTPUT DOS ARQUIVOS APÓS A MONTAGEM"
mkdir -p "$output_dir"

for cfg in "${input_dir}"/config_*.txt; do
    sample=$(basename "$cfg" .txt)

    SOAPdenovo-127mer all \
        -s "$cfg" \
        -K 77 \
        -o "${output_dir}/${sample}_denovo" \
        -d 1 -M 3 -R -u -F \
        -p 16
done

