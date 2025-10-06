#!/bin/bash
#SBATCH -J qiime2                 # Job name
#SBATCH -p batch-AMD              # Partition/queue
#SBATCH -n 1                      # Number of nodes
#SBATCH --ntasks-per-node=1       # 10 CPU cores per node
#SBATCH --cpus-per-task=12        # CPUs per task

source ~/.bashrc
conda activate qiime2-2024.5

qiime rescript get-silva-data \
   --p-version 138.1 \
   --p-target SSURef_NR99 \
   --o-silva-sequences silva-138.1-ssu-nr99-seqs.qza \
   --o-silva-taxonomy silva-138.1-ssu-nr99-tax.qza

qiime rescript reverse-transcribe \
  --i-rna-sequences silva-138.1-ssu-nr99-seqs.qza \
  --o-dna-sequences silva-138.1-ssu-nr99-seqs-dna.qza
