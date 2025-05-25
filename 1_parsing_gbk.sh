mkdir ./data/splitted_faa ./data/splitted_fna


find ./data/gbk -name "*.gbk" | \
xargs -P 2 -I {} python ./scripts/split_gbk_by_protocluster_faa.py {} ./data/splitted_faa


find ./data/gbk -name "*.gbk" | \
xargs -P 2 -I {} python ./scripts/split_gbk_by_protocluster_fna.py {} ./data/splitted_fna

