import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re
import ast

def draw_heatmap_for_modules(filtered_db_path, modules_annotation_path, modules_df):
    """
    Analyze immune content in modules and plot summary graphs.

    Parameters:
    -----------
    filtered_db_path : str
        Path to filtered_DB_2.csv containing cluster_id and protein_id.
    modules_annotation_path : str
        Path to modules_annotation.csv containing immune gene annotations.
    modules_df : pd.DataFrame
        DataFrame with modules info, must contain 'module_id', 'cluster_ids', and 'mmseq_protein_ids' columns.
        'cluster_ids' can be lists or string representations of lists.

    Returns:
    --------
    None
        Plots heatmap and histogram summarizing immune fraction per module.
    """

    # Load filtered DB with protein to cluster mapping
    filtered_DB = pd.read_csv(filtered_db_path, sep=',')
    
    # Group proteins by cluster
    clusters = filtered_DB.groupby('cluster_id')['protein_id'].apply(list).reset_index()
    clusters['size'] = clusters['protein_id'].apply(len)
    
    # Function to safely parse cluster_ids strings
    def safe_literal_eval(val):
        try:
            return ast.literal_eval(val)
        except (ValueError, SyntaxError):
            return val

    # Normalize cluster names by stripping region suffixes
    def normalize_cluster_name(cluster_name):
        cluster_name = re.sub(r'_protocluster_\d+_upstream', '', cluster_name)
        cluster_name = re.sub(r'_protocluster_\d+_downstream', '', cluster_name)
        return cluster_name

    # Parse and normalize cluster_ids in modules_df
    modules_df['cluster_ids'] = modules_df['cluster_ids'].apply(safe_literal_eval)
    modules_df['normalized_cluster_id'] = modules_df['cluster_ids'].apply(
        lambda x: [normalize_cluster_name(c) for c in x]
    )
    
    # Map protein lists to modules
    module_protein_ids = {}
    for module in modules_df['module_id'].unique():
        pcs_list = modules_df.loc[modules_df['module_id'] == module, 'cluster_ids'].values[0]
        cluster_proteins = []
        for pcs_name in pcs_list:
            proteins = clusters.loc[clusters['cluster_id'] == pcs_name, 'protein_id']
            if not proteins.empty:
                # Explode protein list into flat list
                protein_list = proteins.iloc[0]
            else:
                protein_list = []
            cluster_proteins.append(protein_list)
        module_protein_ids[module] = cluster_proteins
    modules_df['mmseq_protein_ids'] = modules_df['module_id'].map(module_protein_ids)

    # Calculate immune island protein counts and fractions
    modules_df['number_mmseq_proteins'] = modules_df['mmseq_protein_ids'].apply(len)
    modules_df['immuneisland_count'] = modules_df['mmseq_protein_ids'].apply(
        lambda x: sum('immuneisland' in str(protein) for protein in x)
    )
    modules_df['immuneisland_fraction'] = modules_df['immuneisland_count'] / modules_df['number_mmseq_proteins']

    # Load immune gene annotations
    immune_genes = pd.read_csv(modules_annotation_path)
    immune_genes['mmseq_protein_ids'] = None

    # Assign protein lists to immune_genes
    for _, module_row in modules_df.iterrows():
        normalized_pcs = module_row['normalized_cluster_id']
        protein_lists = module_row['mmseq_protein_ids']
        if len(normalized_pcs) != len(protein_lists):
            continue
        for cluster, proteins in zip(normalized_pcs, protein_lists):
            match_indices = immune_genes.index[immune_genes['cluster_id'] == cluster]
            for idx in match_indices:
                immune_genes.at[idx, 'mmseq_protein_ids'] = proteins

    # Aggregate immune data by module
    aggregated_data = immune_genes.groupby('module_id').agg(
        total_cluster_size=('cluster_size', 'sum'),
        total_immune=('immune', 'sum'),
        protein_ids=('mmseq_protein_ids', lambda x: sum([p for p in x if isinstance(p, list)], []))
    ).reset_index()

    aggregated_data['immune_to_total_ratio'] = aggregated_data['total_immune'] / aggregated_data['total_cluster_size']

    # Create bins for module size
    aggregated_data['size_bin'] = pd.cut(
        aggregated_data['total_cluster_size'],
        bins=[0, 5, 10, 15, 500000],
        labels=['<5', '5-9', '10-14', '15+'],
        right=False
    )

    # Heatmap data pivot table
    heatmap_data = aggregated_data.pivot_table(
        index='size_bin',
        columns=pd.cut(aggregated_data['immune_to_total_ratio'], bins=5),
        values='total_cluster_size',
        aggfunc='sum',
        fill_value=0
    )

    # Plot heatmap
    plt.figure(figsize=(6, 4))
    sns.heatmap(
        heatmap_data,
        annot=True,
        fmt='d',
        cmap='Blues',
        cbar_kws={'label': 'Number of proteins'},
        linewidths=0.5,
        annot_kws={"size": 11}
    )
    plt.title('Number of proteins in modules vs Immune island fraction')
    plt.xlabel('Immune island fraction')
    plt.ylabel('Number of proteins in a module')
    plt.xticks(rotation=45, fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.show()

    # Plot histogram of immune island fraction
    plt.figure(figsize=(12, 6))
    ax = sns.histplot(
        data=modules_df,
        x='immuneisland_fraction',
        bins=20,
        edgecolor='black'
    )
    for patch in ax.patches:
        height = patch.get_height()
        if height > 0:
            ax.text(
                x=patch.get_x() + patch.get_width() / 2,
                y=height + 0.5,
                s=f'{int(height)}',
                ha='center',
                va='bottom',
                fontsize=10
            )
    plt.title('Distribution of Immune Island fraction in modules')
    plt.xlabel('Immune Island fraction')
    plt.ylabel('Number of modules')
    plt.grid(axis='y', alpha=0.2)
    plt.tight_layout()
    plt.show()