#!/bin/bash
#SBATCH --job-name=metagene
#SBATCH --output=metagene.log
#SBATCH --error=metagene.err
#SBATCH --ntasks=16

# caminhos arquivos
scaf_single="CAMINHO DAS READS MONTADAS APÓS SOAPDENOVO"
scaf_mixed="CAMINHO DAS READS MONTADAS NOVAMENTE APÓS SOAPDENOVO"
out_dir="CAMINHO DOS OUTPUTS A SEREM GERADOS"
mkdir -p "$out_dir"

# caminho metagene
mgm_dir="/caminho/do/software/MetaGeneMark_linux_64"
mgm_bin="${mgm_dir}/mgm/gmhmmp"
mgm_mod="${mgm_dir}/mgm/MetaGeneMark_v1.mod"
export GM_KEY="${mgm_dir}/gm_key"

for scaf in "$scaf_single"/*.scafSeq; do
    sample=$(basename "$scaf" .scafSeq)

    mixed_scaf="${scaf_mixed}/${sample}_mixed_denovo.scafSeq"
    combined_fasta="${out_dir}/${sample}_scaffolds_ge500.fasta"
    > "$combined_fasta"

    # filtrando as amostras normais (sem ser o mixed)
    if [[ -f "$scaf" ]]; then
        seqkit seq -m 500 "$scaf" | sed "s/^>/>${sample}_/" >> "$combined_fasta"
    fi

    # filtrando (>500 bp) para as mesmas amostras, mas o mixed
    if [[ -f "$mixed_scaf" ]]; then
        seqkit seq -m 500 "$mixed_scaf" | sed "s/^>/>${sample}_mixed_/" >> "$combined_fasta"
    fi

    genes_raw="${out_dir}/${sample}_genes.gff"
    genes_faa="${out_dir}/${sample}_proteins.faa"
    genes_ffn="${out_dir}/${sample}_genes.ffn"

    $mgm_bin -a -d -f G \
    -m "$mgm_mod" \
    "$combined_fasta" \
    -o "$genes_raw" \
    -A "$genes_faa" \
    -D "$genes_ffn"

    genes_filtered="${out_dir}/${sample}_genes_ge100.ffn"
    seqkit seq -m 100 "$genes_ffn" > "$genes_filtered"

done


