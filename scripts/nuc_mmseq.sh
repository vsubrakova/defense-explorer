# Representative Sequence Selection

# create folder for nucleotide seqs 
echo "Nucleic acid clustering"

mkdir nuc_mmseq2_results
cd ./nuc_mmseq2_results
# create folder for logs
mkdir 0_logs

# create folder for database
mkdir 1_createDB
# merge all parsed files
find ../test_data/nuc_files -name '*.fna' -exec cat {} + > concat.fna
echo "Number of nucleotide seqs in files: $(grep '>' concat.fna | wc -l)"

# create database
mmseqs createdb ./concat.fna  ./1_createDB/DB_nuc -v 3 &> ./0_logs/creatdb_log.txt
# create folder for clustering results and temporary files
mkdir -p 2_cluster_results/tmp
# clustering nucleotide seqs with params (--min-seq-id 0.95 -c 0.35 --cov-mode 1 --cluster-mode 2)
mmseqs cluster ./1_createDB/DB_nuc ./2_cluster_results/DB_clu_nuc ./2_cluster_results/tmp --min-seq-id 0.95 -c 0.35 --cov-mode 1 --cluster-mode 2 --cluster-reassign --threads 25 1>./0_logs/cluster_log.txt

# create folder for subDB 
mkdir 3_create_subDB
mmseqs createsubdb ./2_cluster_results/DB_clu_nuc ./1_createDB/DB_nuc ./3_create_subDB/DB_rep_clu_nuc 1>./0_logs/subdb_log.txt

# create folder for representative seqs
mkdir 4_create_repBD
mmseqs convert2fasta ./3_create_subDB/DB_rep_clu_nuc ./4_create_repBD/DB_rep_clu_nuc.fasta 1>./0_logs/fasta_log.txt
# calculate final number of representative seqs
echo "Number of representative seqs: $(grep '>' ./4_create_repBD/DB_rep_clu_nuc.fasta| wc -l)"
grep '>' ./4_create_repBD/DB_rep_clu_nuc.fasta > ./4_create_repBD/headers.txt
