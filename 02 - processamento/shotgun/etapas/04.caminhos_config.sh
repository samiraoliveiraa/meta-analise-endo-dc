#!/bin/bash
#SBATCH --job-name=config_txt_soap
#SBATCH --output=config_soap.log
#SBATCH --error=config_soap.err
#SBATCH --ntasks=16

# nesse script, há a conversão dos arquivos fasta não mapeados do soapaligner para fastq
# e a criação dos arquivos configs para input no soapdenovo

input_dir="CAMINHO DOS ARQUIVOS APÓS SOAPALIGNER CONTRA O GENOMA HUMANO"
output_dir_conversao="CAMINHO DE OUTPUT APÓS A CONVERSÃO EM FASTQ"
output_dir_config="CAMINHO DE OUTPUT PARA OS ARQUIVOS DE CONFIGURAÇÃO"

mkdir -p $output_dir_conversao
mkdir -p $output_dir_config

for fa in "${input_dir}"/*_unmapped.fastq; do
    [ -e "$fa" ] || continue
    sample=$(basename "$fa" _unmapped.fastq)

    fq="${output_dir_conversao}/${sample}_unmapped_converted.fastq"

    seqtk seq -F fastq "$fa" > "$fq"
    cat > "${output_dir_config}"/config_${sample}.txt <<EOF
max_rd_len=150

[LIB]
avg_ins=350
reverse_seq=0
asm_flags=3
rank=1
q=${fq}
EOF
done


