# Search of novel immunity systems in the metagenomes of microbial communities

Bacteria are constantly exposed to the threat of infection by bacteriophages, and in response they have developed a variety of defense systems. The best known systems include the restriction-modification (R-M) system, which recognizes and cuts foreign DNA, and the CRISPR-Cas system, which provides acquired immunity by “remembering” previously attacked phages. In the genomes of microorganisms, immune systems often colocalize with each other, forming protective islands. These regions are an excellent target for finding new immune systems, since previously unexplored systems can colocalize with already known ones.

Thus, the goal of the project is to **find new immune systems in bacteria** by analyzing genes that colocalize with already known systems on a new metagenome dataset.

## Preparations

Download the test dataset by running:
```(bash)
git clone git@github.com:vsubrakova/project_KOT1.git
cd project_KOT1
unzip loci.zip
```

## Pipeline

### MMseqs2 part 0 - protein clustering

Run protein clustering with MMseqs2:
```(bash)
mkdir ./result_mmseq_270
mmseqs createdb ./split_faa/*.faa ./result_mmseq_270/DB_270
mkdir ./result_mmseq_270/tmp
mmseqs cluster --min-seq-id 0.3 -c 0.8 ./result_mmseq_270/DB_270 ./result_mmseq_270/DB_270_clu tmp
mmseqs createtsv ./result_mmseq_270/DB_270 ./result_mmseq_270/DB_270 ./result_mmseq_270/DB_270_clu ./result_mmseq_270/DB_270_clu.tsv
```

### Python part 1 - similarity matrix creation

### R part 2 - co-localisation of genes
