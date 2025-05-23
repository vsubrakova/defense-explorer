mamba activate spacedust

mkdir splitted_faa splitted_fna


find ./data -name "*.gbk" | \
xargs -P 2 -I {} python ./scripts/split_gbk_by_protocluster_faa.py {} ./splitted_faa


find ./data -name "*.gbk" | \
xargs -P 2 -I {} python ./scripts/split_gbk_by_protocluster_fna.py {} ./splitted_fna

