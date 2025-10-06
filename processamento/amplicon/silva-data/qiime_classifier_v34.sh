#!/bin/bash
#SBATCH -J qiime2                 # Job name
#SBATCH -p batch-AMD              # Partition/queue
#SBATCH -n 1                      # Number of nodes
#SBATCH --ntasks-per-node=1       # 10 CPU cores per node
#SBATCH --cpus-per-task=12        # CPUs per task

source ~/.bashrc
conda activate qiime2-2024.5

qiime feature-classifier extract-reads \
  --i-sequences silva-138.1-ssu-nr99-seqs-dna.qza \
  --p-f-primer CCTACGGGRSGCAGCAG \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-n-jobs 12 \
  --o-reads silva-138.1-v34-seqs.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138.1-v34-seqs.qza \
  --i-reference-taxonomy silva-138.1-ssu-nr99-tax.qza \
  --o-classifier silva-138.1-v34-classifier.qza
