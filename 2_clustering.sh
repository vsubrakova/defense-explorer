#!/bin/bash

# get representative sequences
bash ./scripts/nuc_mmseq.sh

# create folder for protein selection results
mkdir ./data/filtrated_faa
cd ./data/filtrated_faa
find ../splitted_faa -name '*.fasta' -exec cat {} + > ./concat.faa
cp ../nuc_mmseq2_results/4_create_repBD/headers.txt .

# select representative sequences
python3 ../scripts/filtrate_faa.py
echo "Number of proteins in start file $(grep '>' concat.faa| wc -l)"
echo "Number of headers in headers $(wc -l headers.txt)"
echo "Number of proteins after filtration $(grep '>' concat_filtrated.faa| wc -l)"
cd ../../

# get protein clusters
bash ./scripts/protein_mmseq.sh
echo "Number of unique clusters $(awk '{print $1}' ./data/mmseq2_results/3_output_tsv/DB_clu.tsv | sort | uniq | wc -l)"
echo "Number of unique proteins $(awk '{print $2}' ./data/mmseq2_results/3_output_tsv/DB_clu.tsv | sort | uniq | wc -l)"
