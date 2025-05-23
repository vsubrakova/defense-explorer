import pandas as pd


def add_distance(result: pd.DataFrame, distance: pd.DataFrame) -> pd.DataFrame:
    """Add average pairwise distance between clusters in one module"""
    new_distance = (
        distance[["module_id", "distance"]]
        .groupby("module_id")
        .agg({"distance": list})
        .reset_index()
    )
    new_distance["average_distance"] = new_distance.distance.apply(
        _calculate_average_distance
    )
    new_distance = new_distance.drop("distance", axis=1)
    result_withdistance = result.merge(new_distance, on="module_id", how="left")
    return result_withdistance


def _calculate_average_distance(distanceces: list):
    """Add average pairwise distance between clusters in one module"""
    counter = 0
    for distance in distanceces:
        counter += distance
    return counter / (len(distanceces))


def add_coocurance(filtered_DB: pd.DataFrame, result: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate protein co-occurrence across clusters within modules.

    Parameters:
    -----------
        filtered_DB: DataFrame with protein-cluster assignments containing:
            - protein_id: Protein identifiers (may have suffixes like '_1')
            - cluster_id: Cluster identifiers
        result: DataFrame with module-cluster relationships containing:
            - module_id: Module identifiers
            - cluster_ids: List of cluster IDs for each module

    Returns:
    -----------
        Original result DataFrame with added 'coocurance' column showing:
        - Number of proteins shared by all clusters in each module
    """
    # Group proteins of cluster
    filtered_DB_copy = filtered_DB.copy()
    # Preprocess protein IDs
    filtered_DB_copy.protein_id = filtered_DB_copy.protein_id.apply(
        lambda x: x.rsplit("_", 1)[0]
    )
    # Group proteins by cluster - create sets of unique proteins per cluster
    cluster_protein = filtered_DB_copy.groupby("cluster_id").agg({"protein_id": set})

    # Prepare module-cluster mapping by exploding lists of clusters
    module_cluster = (
        result[["module_id", "cluster_ids"]]
        .set_index(["module_id"])
        .apply(pd.Series.explode)
        .reset_index()
        .rename(columns={"cluster_ids": "cluster_id"})
    )
    # Merge cluster proteins with module-cluster mapping
    # Results in module → cluster → proteins structure
    module_protein = (
        cluster_protein.merge(module_cluster, on="cluster_id", how="right")[
            ["module_id", "protein_id"]
        ]
        .groupby("module_id")
        .agg({"protein_id": list})
        .reset_index()
    )
    module_protein["coocurance"] = module_protein.protein_id.apply(
        lambda x: len(set.intersection(*x))
    )  # Intersection of all protein sets
    module_coocurance = module_protein.drop("protein_id", axis=1)

    result_with_coocurance = result.merge(module_coocurance, on="module_id")
    return result_with_coocurance


def _cluster_immuneproteins(
    DB_clu: pd.DataFrame, raw_annotation: pd.DataFrame
) -> pd.DataFrame:
    """
    Process cluster data to count immune proteins in each cluster.

    Parameters:
    -----------
    DB_clu : pd.DataFrame
        Input cluster data containing protein clusters. Expected columns:
        - cluster_id: cluster identifiers
        - protein_id: protein identifiers

    raw_annotation : pd.DataFrame
        Reference annotation data containing known sequences. Expected column:
        - seq_id: sequence identifiers to compare against

    Returns:
    --------
    pd.DataFrame
        Processed data with columns:
        - cluster_id: transformed cluster identifiers
        - immune: count of immune proteins in each cluster
    """
    raw_data = DB_clu.copy()
    raw_data.loc[:, "cluster_id"] = raw_data["cluster_id"].apply(_name_transformation)
    raw_data.loc[:, "protein_id"] = raw_data["protein_id"].apply(_name_transformation)
    raw_data["immune"] = ~raw_data["protein_id"].isin(raw_annotation["seq_id"])
    raw_data = raw_data.drop("protein_id", axis=1)
    raw_data = raw_data.groupby("cluster_id").agg({"immune": "sum"}).reset_index()
    return raw_data


def _name_transformation(name: str) -> str:
    """
    Delete from cluster and protein names parts protocluster and region
    """
    split_parts = name.split("_protocluster_", 1)
    last_number = split_parts[1].rsplit("_", 1)[1]
    name_transformed = split_parts[0] + "_" + last_number
    return name_transformed


def add_immune_numbers(
    result: pd.DataFrame, DB_clu: pd.DataFrame, raw_annotation: pd.DataFrame
) -> pd.DataFrame:
    """
    Adds immune protein counts to each module in the result DataFrame.

    Parameters:
    -----------
    result : pd.DataFrame
        DataFrame containing module information with columns:
        - cluster_ids: list-like of cluster identifiers per module
        - module_id: module identifiers

    DB_clu : pd.DataFrame
        Cluster data containing protein information with columns:
        - cluster_id: cluster identifiers
        - protein_id: protein identifiers

    raw_annotation : pd.DataFrame
        Reference annotation data containing known sequences with column:
        - seq_id: sequence identifiers to compare against

    Returns:
    --------
    pd.DataFrame
        Original result DataFrame with added column:
        - immune: count of immune proteins aggregated per module
    """
    cluster_module = (
        result[["cluster_ids", "module_id"]]
        .set_index(["module_id"])
        .apply(pd.Series.explode)
        .reset_index()
        .rename(columns={"cluster_ids": "cluster_id"})
    )
    cluster_module.cluster_id = cluster_module.cluster_id.apply(_name_transformation)
    immune_proteins = _cluster_immuneproteins(DB_clu, raw_annotation)
    module_immune = (
        cluster_module.merge(immune_proteins, on="cluster_id")
        .groupby("module_id")
        .agg({"immune": "sum"})
        .reset_index()
    )
    result_with_immune = result.merge(module_immune, on="module_id")
    return result_with_immune
