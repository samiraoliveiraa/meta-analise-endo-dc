#!/bin/bash
#SBATCH -J qiime2                 # Job name
#SBATCH -p batch-AMD              # Partition/queue
#SBATCH -n 1                      # Number of nodes
#SBATCH --ntasks-per-node=1       # 10 CPU cores per node

source ~/.bashrc
conda activate qiime2-2024.5

DIR_NAME=$(basename "$PWD")

mkdir -p temp

# Importando dados Single-End
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "${DIR_NAME}.tsv" \
  --output-path "temp/demux_${DIR_NAME}.qza" \
  --input-format SingleEndFastqManifestPhred33V2

# Remoção de primers
#qiime cutadapt trim-single \
#  --i-demultiplexed-sequences "temp/demux_${DIR_NAME}.qza" \
#  --p-adapter CCTACGGGRSGCAGCAG \
#  --o-trimmed-sequences "temp/demux_${DIR_NAME}-trimmed.qza"

# Visualização
qiime demux summarize \
  --i-data "temp/demux_${DIR_NAME}.qza" \
  --o-visualization "demux_${DIR_NAME}.qzv"
