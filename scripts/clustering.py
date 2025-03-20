import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import csr_matrix


# читаем TSV результат таблицы mmseq2
DB_270_clu = pd.read_csv('/home/vera/Desktop/project_KOT1/filtered_DB_270.csv', sep='\t')
# группируем данные по cluster_id и вычисляем количество белков в кластере mmseq2
DB_270_sizes = DB_270_clu.groupby('cluster_id')['protein_id'].count().reset_index(name='cluster_size')
# отбираем кластеры с одним белком
clusters_to_remove = DB_270_sizes[DB_270_sizes['cluster_size'] == 1]['cluster_id'].tolist()
# создаем отфильтрованный датасет
filtered_DB_270 = DB_270_clu[~DB_270_clu['cluster_id'].isin(clusters_to_remove)]

cluster_names = filtered_DB_270["cluster_id"].unique()

# Создаем датафрейм для заполнения 0 и 1
# Считаем количетво файлов= количество геномов
genome_names = list(set(i.split('_')[0] for i in cluster_names))
genomes_number = len(genome_names)

# Считаем количетво кластеров
cluster_number = len(filtered_DB_270['cluster_id'].unique())
cluster_names = list(filtered_DB_270['cluster_id'].unique())
columns_names =  cluster_names

# Создаем датафрейм заполненный 0, количество колонок = cluster_number+1 (сами колнки называются по белкам), количество строк = genomes_number
genome_cluster = pd.DataFrame(0, index=np.arange(genomes_number), columns=columns_names)
genome_cluster

genome_name = pd.Series(genome_names, name="genome_name")
genome_cluster = pd.concat([genome_name, genome_cluster], axis=1)
genome_cluster


for cluster_id, protein_id in zip(filtered_DB_270.cluster_id, filtered_DB_270.protein_id):
    for cluster in genome_cluster.columns:
        if cluster == cluster_id:
            if protein_id.split('_')[0] in list(genome_cluster['genome_name']):
                ind = genome_cluster[genome_cluster.genome_name == protein_id.split('_')[0]].index
                genome_cluster.loc[ind, cluster] = 1
        else:
            continue

# Convert the matrix to a sparse format
A_sparse = csr_matrix(genome_cluster.drop('genome_name', axis=1).to_numpy())

intersection_matrix = A_sparse.T @ A_sparse

set_sizes = np.array(A_sparse.sum(axis=0)).flatten()

union_matrix = np.add.outer(set_sizes, set_sizes) - intersection_matrix.toarray()

jaccard_matrix = intersection_matrix.toarray() / union_matrix
jaccard_matrix[np.isnan(jaccard_matrix)] = 0

sns.clustermap(jaccard_matrix,
            cmap="YlGnBu",
            #    xticklabels=names,
            #    yticklabels=names,
            figsize=(12, 10))


jaccard_matrix = pd.DataFrame(jaccard_matrix)
jaccard_matrix.index = genome_cluster.columns[1:]
jaccard_matrix.colnames = genome_cluster.columns[1:]


jaccard_matrix.to_csv("/home/vera/Desktop/project_KOT1/scripts/jaccard_genome_cluster.csv", sep="\t")
