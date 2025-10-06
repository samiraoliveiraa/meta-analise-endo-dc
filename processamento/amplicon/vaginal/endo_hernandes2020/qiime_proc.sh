#!/bin/bash
#SBATCH -J qiime2                 # Job name
#SBATCH -p batch-AMD              # Partition/queue
#SBATCH -n 1                      # Number of nodes
#SBATCH --ntasks-per-node=1       # 10 CPU cores per node
#SBATCH --cpus-per-task=12         # CPUs per task

source ~/.bashrc
conda activate qiime2-2024.5

DIR_NAME=$(basename "$PWD")

qiime dada2 denoise-single \
  --i-demultiplexed-seqs "temp/demux_${DIR_NAME}.qza" \
  --p-trim-left 20 \
  --p-trunc-len 200 \
  --p-n-threads 4 \
  --o-table "temp/table_${DIR_NAME}.qza" \
  --o-representative-sequences "temp/rep-seqs_${DIR_NAME}.qza" \
  --o-denoising-stats "temp/stats_${DIR_NAME}.qza"

qiime feature-classifier classify-sklearn \
  --i-classifier /home/public/microbiota_jubs_sams/03-qiime/silva-data/silva-138.1-v34-classifier.qza \
  --i-reads "temp/rep-seqs_${DIR_NAME}.qza" \
  --o-classification "temp/taxonomy_${DIR_NAME}.qza" \
  --p-n-jobs 4

qiime taxa filter-table \
  --i-table "temp/table_${DIR_NAME}.qza" \
  --i-taxonomy "temp/taxonomy_${DIR_NAME}.qza" \
  --p-include d__Bacteria \
  --p-exclude mitochondria \
  --o-filtered-table "temp/table_${DIR_NAME}_filtered.qza"

qiime taxa filter-seqs \
  --i-sequences "temp/rep-seqs_${DIR_NAME}.qza" \
  --i-taxonomy "temp/taxonomy_${DIR_NAME}.qza" \
  --p-include d__Bacteria \
  --p-exclude mitochondria \
  --o-filtered-sequences "temp/sequences_${DIR_NAME}_filtered.qza"

qiime taxa barplot \
  --i-table "temp/table_${DIR_NAME}_filtered.qza" \
  --i-taxonomy "temp/taxonomy_${DIR_NAME}.qza" \
  --m-metadata-file "${DIR_NAME}_metadata.tsv" \
  --o-visualization "bar-plots_${DIR_NAME}.qzv"

qiime diversity core-metrics \
  --i-table "temp/table_${DIR_NAME}_filtered.qza" \
  --p-sampling-depth 10000 \
  --m-metadata-file "${DIR_NAME}_metadata.tsv" \
  --output-dir diversity

cd diversity

qiime tools export \
  --input-path "bray_curtis_distance_matrix.qza" \
  --output-path pasta_temp

mv pasta_temp/distance-matrix.tsv ../bc_${DIR_NAME}.tsv

rm -r pasta_temp

qiime diversity alpha-group-significance \
  --i-alpha-diversity "shannon_vector.qza" \
  --m-metadata-file "../${DIR_NAME}_metadata.tsv" \
  --o-visualization "../shannon_${DIR_NAME}.qzv"

cd ..

qiime picrust2 full-pipeline \
   --i-table "temp/table_${DIR_NAME}_filtered.qza" \
   --i-seq "temp/sequences_${DIR_NAME}_filtered.qza" \
   --output-dir picrust2 \
   --p-placement-tool sepp \
   --p-threads 12 \
   --p-hsp-method pic \
   --p-max-nsti 2 \

cd picrust2

qiime taxa barplot \
  --i-table "ko_metagenome.qza" \
  --m-metadata-file "../${DIR_NAME}_metadata.tsv" \
  --o-visualization "ko_barplot_${DIR_NAME}.qzv"

qiime taxa barplot \
  --i-table "ec_metagenome.qza" \
  --m-metadata-file "../${DIR_NAME}_metadata.tsv" \
  --o-visualization "ec_barplot_${DIR_NAME}.qzv"

qiime tools export \
  --input-path "pathway_abundance.qza" \
  --output-path pasta_temp

biom convert \
  -i pasta_temp/feature-table.biom \
  -o "../pathway_${DIR_NAME}.tsv" \
  --to-tsv

rm -r pasta_temp

cd ..

qiime feature-table transpose \
  --i-table "temp/table_${DIR_NAME}_filtered.qza" \
  --o-transposed-feature-table "temp/table_${DIR_NAME}_trans.qza"

qiime metnet generateFeatures \
  --i-frequency "temp/table_${DIR_NAME}_trans.qza" \
  --i-taxa "temp/taxonomy_${DIR_NAME}.qza" \
  --p-selection AGREDA \
  --p-level g \
  --o-reactions "temp/output_reactions.qza" \
  --o-subsystems "temp/output_subsystems.qza" \
  --o-xmatrix "temp/output_X_matrix.qza"

qiime metnet differentialReactions \
  --i-reactions "temp/output_reactions.qza" \
  --m-metadata-file "${DIR_NAME}_metadata.tsv" \
  --m-metadata-column group \
  --p-condition-name ${DIR_NAME%_*} \
  --p-control-name controle \
  --p-selection-model AGREDA \
  --o-differential-analysis "temp/reactions_differential.qza"

qiime metnet differentialSubSystems \
  --i-subsystems "temp/output_subsystems.qza" \
  --m-metadata-file "${DIR_NAME}_metadata.tsv" \
  --m-metadata-column group \
  --p-condition-name ${DIR_NAME%_*} \
  --p-control-name controle \
  --o-differential-analysis "temp/subsystems_differential.qza"

qiime metadata tabulate \
  --m-input-file "temp/subsystems_differential.qza" \
  --o-visualization "subsystems_${DIR_NAME}.qzv"