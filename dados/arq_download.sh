#!/bin/bash
#SBATCH --job-name=download-data
#SBATCH --output=download.log
#SBATCH --error=download.err
#SBATCH --time=50:00:00
#SBATCH --ntasks=1

bash
conda activate bioinfo # ativar o ambiente que contém o sra_tools instalado

artigo=$(basename "$PWD")

base_dir="/home/public/microbiota_jubs_sams/01-dados/intestinal/amplicon/$artigo"

mkdir -p "$base_dir"

for arq_grupo in ${artigo}_*.txt; do

    [ -e "$arq_grupo" ] || continue
    grupo=$(echo "$arq_grupo" | sed -E "s/${artigo}_(.*)\.txt/\1/")

    grupo_dir="${base_dir}/${grupo}"
    mkdir -p "$grupo_dir"

    while IFS= read -r id || [[ -n "$id" ]]; do
        if [[ -n "$id" ]]; then
            echo "Baixando $id para $grupo_dir"
            fasterq-dump "$id" -O "$grupo_dir" --split-files --threads 4
        fi
    done < "$arq_grupo"

done