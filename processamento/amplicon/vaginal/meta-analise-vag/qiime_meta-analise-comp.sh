#!/bin/bash
#SBATCH -J qiime2_meta           # Nome do Job
#SBATCH -p batch-AMD              # Partição/Fila
#SBATCH -n 1                      # Número de nós
#SBATCH --ntasks-per-node=1       # Tarefas por nó
#SBATCH --cpus-per-task=12         # CPUs por tarefa para paralelizar

source ~/.bashrc
conda activate qiime2-2024.5

mkdir -p temp

base_dir="/home/public/microbiota_jubs_sams/03-qiime/vaginal/"
output_path="/home/public/microbiota_jubs_sams/03-qiime/vaginal/meta-analise-vag/temp"

STUDIES="endo_hernandes2020 endo_jimenez2024"


for DIR_NAME in $STUDIES; do

    estudo_path="${base_dir}/${DIR_NAME}"

    echo "Processando o estudo: ${DIR_NAME}"

    qiime feature-classifier classify-sklearn \
      --i-classifier "/home/public/microbiota_jubs_sams/03-qiime/silva-data/silva-138.1-classifier.qza" \
      --i-reads "${estudo_path}/temp/rep-seqs_${DIR_NAME}.qza" \
      --o-classification "${output_path}/taxonomy_${DIR_NAME}.qza" \
      --p-n-jobs 4

     qiime taxa filter-table \
      --i-table "${estudo_path}/temp/table_${DIR_NAME}.qza" \
      --i-taxonomy "${output_path}/taxonomy_${DIR_NAME}.qza" \
      --p-include d__Bacteria \
      --p-exclude mitochondria \
      --o-filtered-table "${output_path}/table_${DIR_NAME}_filtered.qza"

    qiime taxa collapse \
      --i-table "${output_path}/table_${DIR_NAME}_filtered.qza" \
      --i-taxonomy "${output_path}/taxonomy_${DIR_NAME}.qza" \
      --p-level 6 \
      --o-collapsed-table "${output_path}/table_genero_${DIR_NAME}.qza"

    cp ${estudo_path}/temp/rep-seqs* ${output_path}

    echo "Estudo ${DIR_NAME} finalizado com sucesso."

done

qiime feature-table merge \
  --i-tables temp/table_genero*.qza \
  --o-merged-table temp/abundancia_genero_unificada.qza

qiime feature-table merge \
  --i-tables temp/table_*filtered.qza \
  --o-merged-table temp/abundancia_unificada.qza

qiime feature-table merge-seqs \
  --i-data temp/rep-seqs* \
  --o-merged-data temp/sequencias_unificadas.qza

qiime feature-table merge-taxa \
  --i-data temp/taxonomy* \
  --o-merged-data temp/taxonomia_unificada.qza

qiime diversity core-metrics \
  --i-table temp/abundancia_genero_unificada.qza \
  --p-sampling-depth 10000 \
  --m-metadata-file meta-analise-vag_metadata.tsv \
  --output-dir diversity

cd diversity

mkdir -p pasta_temp

qiime tools export \
  --input-path bray_curtis_distance_matrix.qza \
  --output-path pasta_temp

mv pasta_temp/distance-matrix.tsv ../bc_meta-analise-vag.tsv

rm -r pasta_temp

qiime diversity beta-group-significance \
  --i-distance-matrix bray_curtis_distance_matrix.qza \
  --m-metadata-file ../meta-analise-vag_metadata.tsv \
  --m-metadata-column group \
  --p-method permanova \
  --o-visualization ../beta_meta-analise-vag.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity shannon_vector.qza \
  --m-metadata-file ../meta-analise-vag_metadata.tsv \
  --o-visualization ../shannon_meta-analise-vag.qzv

cd ..

qiime picrust2 full-pipeline \
   --i-table "temp/abundancia_unificada.qza" \
   --i-seq "temp/sequencias_unificadas.qza" \
   --output-dir picrust2 \
   --p-placement-tool sepp \
   --p-threads 12 \
   --p-hsp-method pic \
   --p-max-nsti 2 

cd picrust2

qiime taxa barplot \
  --i-table "ko_metagenome.qza" \
  --m-metadata-file "../meta-analise-vag_metadata.tsv" \
  --o-visualization "ko_barplot_meta-analise-vag.qzv"

qiime taxa barplot \
  --i-table "ec_metagenome.qza" \
  --m-metadata-file "../meta-analise-vag_metadata.tsv" \
  --o-visualization "ec_barplot_meta-analise-vag.qzv"

qiime tools export \
  --input-path "pathway_abundance.qza" \
  --output-path pasta_temp

biom convert \
  -i pasta_temp/feature-table.biom \
  -o "../pathway_meta-analise-vag.tsv" \
  --to-tsv

rm -r pasta_temp

cd ..

qiime feature-table transpose \
  --i-table "temp/abundancia_unificada.qza" \
  --o-transposed-feature-table "temp/table_meta-analise-vag_trans.qza"

qiime metnet generateFeatures \
  --i-frequency "temp/table_meta-analise-vag_trans.qza" \
  --i-taxa "temp/taxonomia_unificada.qza" \
  --p-selection AGREDA \
  --p-level g \
  --o-reactions "temp/output_reactions.qza" \
  --o-subsystems "temp/output_subsystems.qza" \
  --o-xmatrix "temp/output_X_matrix.qza"

qiime metnet differentialReactions \
  --i-reactions "temp/output_reactions.qza" \
  --m-metadata-file "meta-analise-vag_metadata.tsv" \
  --m-metadata-column group \
  --p-condition-name "Endometriose" \
  --p-control-name "Controle" \
  --p-selection-model AGREDA \
  --o-differential-analysis "temp/reactions_differential_endo.qza"

qiime metnet differentialReactions \
  --i-reactions "temp/output_reactions.qza" \
  --m-metadata-file "meta-analise-vag_metadata.tsv" \
  --m-metadata-column group \
  --p-condition-name "Dor Crônica" \
  --p-control-name "Controle" \
  --p-selection-model AGREDA \
  --o-differential-analysis "temp/reactions_differential_dor.qza"

qiime metnet differentialSubSystems \
  --i-subsystems "temp/output_subsystems.qza" \
  --m-metadata-file "meta-analise-vag_metadata.tsv" \
  --m-metadata-column group \
  --p-condition-name "Endometriose" \
  --p-control-name "Controle" \
  --o-differential-analysis "temp/subsystems_differential_endo.qza"

qiime metnet differentialSubSystems \
  --i-subsystems "temp/output_subsystems.qza" \
  --m-metadata-file "meta-analise-vag_metadata.tsv" \
  --m-metadata-column group \
  --p-condition-name "Dor Crônica" \
  --p-control-name "Controle" \
  --o-differential-analysis "temp/subsystems_differential_dor.qza"

qiime metadata tabulate \
  --m-input-file "temp/subsystems_differential_endo.qza" \
  --o-visualization "subsystems_meta-analise-vag_endo.qzv"

qiime metadata tabulate \
  --m-input-file "temp/subsystems_differential_cppendo.qza" \
  --o-visualization "subsystems_meta-analise-vag_dor.qzv"

