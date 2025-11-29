#!/bin/bash
#SBATCH --job-name=soapaligner_mix
#SBATCH --output=soapaligner_mix.log
#SBATCH --error=soapaligner_mix.err
#SBATCH --ntasks=16

soapdenovo_dir="CAMINHO DOS MONTADOS APÓS SOAPDENOVO"
reads_dir="CAMINHO DAS READS APÓS O SOAPALIGNER (REMOÇÃO DO GENOMA HUMANO)"
output_dir="OUTPUT DO ALINHAMENTO ENTRE AS READS ORIGINAIS E AS MONTAGENS"
mkdir -p "$output_dir"

for scaf in "$soapdenovo_dir"/*.scafSeq; do
    sample=$(basename "$scaf" ".scafSeq")

    2bwt-builder "$scaf"

    sra_id=$(echo "$sample" | grep -o 'SRR[0-9]\+')
    r1="${reads_dir}/${sra_id}_unmapped_converted.fastq"

    if [[ -f "$r1" ]]; then
        soap -D "${scaf}.index" \
             -a "$r1" \
             -o "${output_dir}/${sra_id}_vs_scaftigs.soap" \
             -u "${output_dir}/${sra_id}_unmapped_again.fastq" \
             -s 50 -l 30 -v 5 -p 16
    fi
done


