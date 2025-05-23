# Protein clustering

# create folder for Protein clustering results
mkdir mmseq2_results
# copy filtreated protein seqs into directory
cp ./filtrated_faa/concat_filtrated.faa ./mmseq2_results

cd ./mmseq2_results
# create folder for logs
mkdir 0_logs
# create folder for database
mkdir 1_createDB
# create database
mmseqs createdb ./concat_filtrated.faa  ./1_createDB/DB -v 3 &> ./0_logs/creatdb_log.txt
# create folder for clustering results and temporary files
mkdir -p 2_cluster_results/tmp
# clustering protein seqs with params --min-seq-id 0.4 -c 0.8
mmseqs cluster ./1_createDB/DB ./2_cluster_results/DB_clu ./2_cluster_results/tmp --min-seq-id 0.4 -c 0.8 --cluster-reassign --threads 30 1>./0_logs/cluster_log.txt
# create folder for final table
mkdir 3_output_tsv
# create table.tsv with 'cluster_id - protein_id' data
mmseqs createtsv ./1_createDB/DB ./1_createDB/DB ./2_cluster_results/DB_clu ./3_output_tsv/DB_clu.tsv --threads 25 1>./0_logs/output_tsv_log.txt
