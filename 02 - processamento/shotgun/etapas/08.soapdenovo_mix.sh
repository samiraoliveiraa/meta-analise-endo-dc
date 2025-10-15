#!/bin/bash
#SBATCH --job-name=mix_soapdenovo
#SBATCH --output=mix_soapdenovo.log
#SBATCH --error=mix_soapdenovo.err
#SBATCH --ntasks=16

input_dir="CAMINHO DOS ARQUIVOS DE CONFIGURAÇÃO DAS READS NÃO MAPEADAS CONTRA AS SEQUÊNCIAS MONTADAS"
output_dir="CAMINHO DO OUTPUT DAS MONTAGENS DESTAS READS"
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

