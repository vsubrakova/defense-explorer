# project_KOT1

## Preparations

Download the test dataset by running:
```(bash)
git clone git@github.com:vsubrakova/project_KOT1.git
cd project_KOT1
unzip loci.zip
```
Also we would need to prepare environment and install Spacedust:
```(bash)
mamba create -n spacedust
mamba activate spacedust
# Install MMseqs2
wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
tar xvf mmseqs-linux-avx2.tar.gz
export PATH=$(pwd)/mmseqs/bin:$PATH
# I have Linux AVX2 build (check using: cat /proc/cpuinfo | grep avx2). For other builds please refer to spacedust manual
wget https://mmseqs.com/spacedust/spacedust-linux-avx2.tar.gz 
tar xvzf spacedust-linux-avx2.tar.gz; export PATH=$(pwd)/spacedust/bin/:$PATH
sudo apt install prodigal

wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz
tar xvzf foldseek-linux-avx2.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH

#insert your path to foldseek
/home/vera/Desktop/project_KOT1/foldseek/bin/foldseek databases Alphafold/UniProt refFoldseekDB tmpFolder
databases Alphafold/UniProt refFoldseekDB tmpFolder 

```
Run Prodigal annotation
```(bash)
mkdir -p prodigal_output

for fna in split_fna/*.fna; do
    base=$(basename "$fna" .fna)
    prodigal -i "$fna" \
             -a "prodigal_output/${base}.faa" \
             -d "prodigal_output/${base}_genes.fna" \
             -o "prodigal_output/${base}.gff" \
             -f gff
done
```
## Running Spacedust

Create annotated protein database for Spacedust
```(bash)
spacedust createsetdb prodigal_output/*.faa setDB tmpFolder --gff-dir prodigal_output/ --gff-type CDS
```

Homology search:
```(bash)
spacedust clustersearch setDB setDB result.tsv tmpFolder --threads 1

# With foldseek
spacedust clustersearch setDB setDB result_foldseek.tsv tmpFolder --search-mode 1

spacedust --remove-tmp-files
```

