#!/bin/bash

# # Запуск в homology_env
source  /opt/miniconda3/bin/activate homology_env
echo "Running in homology env"
bash ./scripts/nuc_mmseq.sh
conda deactivate

cd ~
# Запуск в python3
rm -rf 2_filtrated
mkdir 2_filtrated
cd ./2_filtrated
source  /opt/miniconda3/bin/activate python3
echo "Running in python3 env"
find /home/v_subrakova/files/new_split_faa -name '*.faa' -exec cat {} + > concat.faa
cp /home/s_borovikova/2_nuc_mmseq2_results/4_create_repBD/headers.txt ./
python3 ../scripts/filtrate_faa.py
conda deactivate
echo "Количесвто белков в исходном файле $(grep '>' concat.faa| wc -l)"
echo "Количесвто белков в headers $(wc -l headers.txt)"
echo "Количесвто белков после фильтрации $(grep '>' concat_filtrated.faa| wc -l)"

cd ~
# Запуск в homology_env
source  /opt/miniconda3/bin/activate homology_env
echo "Running in homology env"
bash ./scripts/protein_mmseq.sh
echo "Number of unique clusters $(awk '{print $1}' ./2_mmseq2_results/3_output_tsv/DB_clu.tsvDB_clu.tsv | sort | uniq | wc -l)"
echo "Number of unique proteins $(awk '{print $2}' ./2_mmseq2_results/3_output_tsv/DB_clu.tsvDB_clu.tsv | sort | uniq | wc -l)"
conda deactivate