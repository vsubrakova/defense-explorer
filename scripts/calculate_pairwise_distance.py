import pandas as pd
import re
import ast

def calculate_pairwise_distance(result: pd.DataFrame, filtered_DB: pd.DataFrame):
    """
    Analyze protein distances between clusters within modules and plot the average distance histogram.
    Adapted to modules DataFrame where cluster_ids are string representations of lists, and cluster_inds are tuples.

    Parameters:
    -----------
    result : pd.DataFrame
        DataFrame containing module information including 'module_id', 'cluster_ids' (list-like string).
    filtered_DB : pd.DataFrame
        DataFrame containing protein data with 'cluster_id' and 'protein_id' columns.

    Returns:
    --------
    pairwise_df : pd.DataFrame
        DataFrame with pairwise protein prefix distances between clusters per module.
    """

    # # Convert cluster_ids string to actual list for each row
    modules = result.copy()

    # Group proteins by cluster_id and aggregate protein IDs into lists
    clusters = filtered_DB.groupby('cluster_id')['protein_id'].apply(list).reset_index()
    clusters['size'] = clusters['protein_id'].apply(len)

    # Map module IDs to protein lists from their respective clusters
    module_protein_ids = {}
    for module in modules['module_id'].unique():
        pcs_list = modules.loc[modules['module_id'] == module, 'cluster_ids'].values[0]
        
        cluster_proteins = []
        for pcs_name in pcs_list:
            proteins = clusters.loc[clusters['cluster_id'] == pcs_name, 'protein_id']
            if not proteins.empty:
                proteins_list = proteins.values[0]
            else:
                proteins_list = []
            cluster_proteins.append(proteins_list)
        
        module_protein_ids[module] = cluster_proteins

    modules['mmseq_protein_ids'] = modules['module_id'].map(module_protein_ids)
    modules['number_mmseq_proteins'] = modules['mmseq_protein_ids'].apply(len)

    # Helper functions to parse protein IDs
    def extract_prefix(protein_id):
        match = re.match(r'(.+)_\d+$', protein_id)
        return match.group(1) if match else None

    def extract_number(protein_id):
        match = re.search(r'_(\d+)$', protein_id)
        return int(match.group(1)) if match else None

    def compute_distance(numbers1, numbers2):
        total_distance = 0
        num_pairs = 0
        for n1 in numbers1:
            for n2 in numbers2:
                total_distance += abs(n1 - n2)
                num_pairs += 1
        return (total_distance / num_pairs) - 1 if num_pairs > 0 else None

    cluster_pair_rows = []

    # Calculate distances between proteins with the same prefix in different clusters for each module
    for idx, row in modules.iterrows():
        module_id = row['module_id']
        pcs_list = row['cluster_ids']
        protein_lists = row['mmseq_protein_ids']
        
        prefix_dict = {}
        for cluster_number, cluster in enumerate(protein_lists):
            for protein in cluster:
                prefix = extract_prefix(protein)
                number = extract_number(protein)
                if prefix and number is not None:
                    prefix_dict.setdefault(prefix, []).append((cluster_number, number))
        
        for prefix, proteins in prefix_dict.items():
            cluster_dict = {}
            for cluster_number, number in proteins:
                cluster_dict.setdefault(cluster_number, []).append(number)
            
            for cluster_id_1, numbers_1 in cluster_dict.items():
                for cluster_id_2, numbers_2 in cluster_dict.items():
                    if cluster_id_1 < cluster_id_2:
                        distance = compute_distance(numbers_1, numbers_2)
                        if distance is not None:
                            cluster_pair_rows.append({
                                'module_id': module_id,
                                'cluster_ids': ', '.join(pcs_list),
                                'mmseq_protein_ids': protein_lists,
                                'prefix': prefix,
                                'numbers_1': numbers_1,
                                'numbers_2': numbers_2,
                                'distance': distance
                            })

    pairwise_df = pd.DataFrame(cluster_pair_rows)
    
    return pairwise_df
