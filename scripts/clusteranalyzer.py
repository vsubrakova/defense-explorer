import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from statsmodels.stats.multitest import multipletests
from scipy.stats import hypergeom
import markov_clustering as mc

pd.set_option("display.max_colwidth", 200)


def filtrate_low_clusters(data: pd.DataFrame, min_cluster_size: int = 0):
    """
    Filtrate mmseqs2 results from low abudant clusters

    Parameters
    ----------
    data : pd.DataFrame
        Input DataFrame with at least two columns:
        - 'protein_id': Unique protein identifiers (e.g., "genomeA_123").
        - 'cluster_id': Cluster identifiers for each protein.

    min_cluster_size :int
        Minimal number of proteins in cluster (default: 0)

    Returns
    -------
    pd.DataFrame: filtrated dataframe

    """
    data = data.drop_duplicates()
    print(
        f'Number of unique regions : {len(data.protein_id.apply(lambda x: x.rsplit("_",1)[0]).unique())}'
    )
    print(f"Number of unique proteins: {len(data.protein_id.unique())}")
    print(f"Number of unique clusters : {len(data.cluster_id.unique())}")

    # Create dataframe with cluster names and sizes
    cluster_sizes = (
        data.groupby("cluster_id")["protein_id"]
        .count()
        .reset_index(name="cluster_size")
    )

    # Remove clusters with n sizes
    clusters_to_remove = cluster_sizes[
        cluster_sizes["cluster_size"] < min_cluster_size
    ]["cluster_id"].tolist()
    filtered_data = data[~data["cluster_id"].isin(clusters_to_remove)]

    print(
        f'Number of unique regions (after filtration) : {len(filtered_data.protein_id.apply(lambda x: x.rsplit("_",1)[0]).unique())}'
    )
    print(
        f"Number of unique proteins (after filtration) : {len(filtered_data.protein_id.unique())}"
    )
    print(
        f"Number of unique clusters (after filtration) : {len(filtered_data.cluster_id.unique())}"
    )
    print(f" Delete clusters containing less than {min_cluster_size} proteins")

    return filtered_data, cluster_sizes


def create_one_zero_matrix(
    data: pd.DataFrame,
) -> tuple[pd.Series, pd.Series, csr_matrix]:
    """
    Creates a binary (1/0) sparse matrix representing the presence of protein clusters in genomes.

    The function processes a DataFrame containing protein IDs and their associated cluster IDs,
    then constructs a sparse matrix where:
    - Rows correspond to genomes (derived from protein IDs).
    - Columns correspond to protein clusters.
    - A value of 1 indicates that a genome contains at least one protein from the cluster.

    Parameters
    ----------
    data : pd.DataFrame
        Input DataFrame with at least two columns:
        - 'protein_id': Unique protein identifiers (e.g., "genomeA_123").
        - 'cluster_id': Cluster identifiers for each protein.

    Returns
    -------
    tuple[pd.Series, pd.Series, csr_matrix]
        A tuple containing:
        - genome_names: pd.Series of unique genome IDs (extracted from protein_id).
        - cluster_names: pd.Series of unique cluster IDs.
        - sparse_matrix: Binary sparse matrix in CSR format (genomes × clusters).

    """
    # Create table protein_genome (with columns protein_id and genome_id)
    genome_names = data.protein_id.apply(lambda x: x.rsplit("_", 1)[0]).rename(
        "genome_id"
    )  # extract genome names from protein_id
    protein_names = data.protein_id
    protein_genome = pd.concat([protein_names, genome_names], axis=1)

    #  Merge protein_genome with cluster IDs and remove duplicates (1 genome/cluster pair = 1 entry)
    protein_genome_cluster = protein_genome.merge(data, on="protein_id", how="outer")
    genome_cluster = protein_genome_cluster[
        ["genome_id", "cluster_id"]
    ].drop_duplicates()

    # Convert to categorical for efficient encoding
    genome_cluster = genome_cluster.astype(
        {"genome_id": "category", "cluster_id": "category"}
    )
    row = genome_cluster.genome_id.cat.codes.to_numpy()  # Genome indices
    col = genome_cluster.cluster_id.cat.codes.to_numpy()  # Cluster indices
    genome_names = genome_cluster.genome_id.cat.categories
    cluster_names = genome_cluster.cluster_id.cat.categories

    # Create sparse matrix (1 if genome has at least 1 protein in the cluster)
    sparse_matrix = csr_matrix(
        (np.ones(len(row), dtype=np.int8), (row, col)),
        shape=(len(genome_names), len(cluster_names)),
    )

    print(
        f"Sparse matrix with {len(genome_names)} rows and {len(cluster_names)} columns was created"
    )

    return genome_names, cluster_names, sparse_matrix


def create_coocurance_matrix(sparse_matrix: csr_matrix) -> csr_matrix:
    """
    Creates a co-occurrence matrix of protein clusters from a genome-cluster binary matrix.

    The co-occurrence matrix represents how often pairs of protein clusters appear together
    in the same genomes. Diagonal elements (self-co-occurrence) are set to 0.

    Parameters
    ----------
    sparse_matrix : csr_matrix
        Binary sparse matrix of shape (n_genomes, n_clusters) where:
        - 1 indicates the presence of a cluster in a genome
        - 0 indicates absence

    Returns
    -------
    csr_matrix
        Symmetric co-occurrence matrix of shape (n_clusters, n_clusters) where:
        - Element [i, j] counts genomes where both clusters i and j are present
        - Diagonal elements are forced to 0

    """
    cooccurrence = sparse_matrix.T.dot(sparse_matrix)
    cooccurrence.setdiag(0)  # Exclude self-pairs
    print(
        f"Coocurance matrix {cooccurrence.shape[0]} rows and {cooccurrence.shape[1]} columns was created"
    )
    return cooccurrence


def calculate_pvals(
    genome_names: pd.Series, cooccurrence_matrix: csr_matrix, sparse_matrix: csr_matrix
) -> np.ndarray:
    """
    Calculates p-values for co-occurrence of protein clusters using the hypergeometric test.

    For each cluster pair, tests whether their observed co-occurrence is significantly
    higher than expected by chance, given:
    - The total number of genomes (N)
    - The prevalence of each individual cluster

    Parameters
    ----------
    genome_names : pd.Index
        Index of genome identifiers (length = N)
    cooccurrence_matrix : csr_matrix
        Co-occurrence matrix from create_cooccurrence_matrix()
    sparse_matrix : csr_matrix
        Original binary matrix of shape (n_genomes, n_clusters)

    Returns
    -------
    np.ndarray
        Square matrix of p-values (n_clusters × n_clusters) where:
        - p-values < 0.05 indicate significant co-occurrence
        - Diagonal elements are set to 1.0 (no self-comparison)
        - NaN values are replaced with 1.0

    Notes
    -----
    Uses scipy.stats.hypergeom.sf for survival function (1 - CDF) calculation.
    Multiple testing correction should be applied separately (e.g., FDR).
    """
    n_total = len(genome_names)  # Total number of genomes

    # Cluster prevalences (number of genomes containing each cluster)
    pc_counts = np.array(sparse_matrix.sum(axis=0)).flatten()

    # Convert co-occurrence to dense array for vectorized operations
    coocurance = cooccurrence_matrix.toarray()
    cluster_1 = pc_counts
    cluster_2 = pc_counts[:, None]  # Transposed view

    # Vectorized hypergeometric test (sf = survival function = 1 - CDF)
    pvals = hypergeom.sf(
        coocurance - 1, n_total, cluster_1, cluster_2
    )  # P(X ≥ observed_cooccurrence)
    pvals = np.nan_to_num(pvals, nan=1.0)  # Handle NaN
    np.fill_diagonal(pvals, 1.0)  # Exclude self-comparisons

    return pvals


def correct_pvals(
    pvals: np.ndarray, method: str = "fdr_bh", alpha: float = 0.05
) -> np.ndarray:
    """
    Corrects p-values for multiple testing

    Applies multiple testing correction to the upper triangle of a symmetric
    p-value matrix and mirrors the results to maintain symmetry. The diagonal
    remains unchanged (set to 1.0).

    Parameters
    ----------
    pvals : np.ndarray
        Square symmetric matrix of raw p-values. Must have shape (n, n).
    method : str, optional
        Multiple testing correction method. Supported methods:
        - 'fdr_bh': Benjamini-Hochberg FDR (default)
        - 'bonferroni': Bonferroni correction
    alpha : float, optional
        Significance level (default: 0.05). Used for method calculations.

    Returns
    -------
    np.ndarray
        Symmetric matrix of corrected p-values with same shape as input.
        Diagonal elements are always 1.0.

    """
    # Extract upper triangle (excluding diagonal)
    rows, cols = np.triu_indices(pvals.shape[0], k=1)
    p_values = pvals[rows, cols]

    # Apply multiple testing correction
    _, corrected, _, _ = multipletests(p_values, alpha=alpha, method=method)

    # Reconstruct full symmetric matrix
    corrected_pvals = np.ones_like(pvals)
    corrected_pvals[rows, cols] = corrected
    corrected_pvals[cols, rows] = corrected

    return corrected_pvals


def calculate_similarity_matrix(
    corrected_pvals: np.ndarray,
    threshold: float = 0,
    min_pval: float = 1e-10,
) -> np.ndarray:
    """
    Converts corrected p-values into a similarity matrix using negative log10 transformation.

    The similarity score between two clusters is computed as:
    similarity = -log10(p-value) if -log10(p-value) > threshold, else 0

    Parameters
    ----------
    corrected_pvals : np.ndarray
        A square matrix of corrected p-values (shape [n, n]) with values in range [0, 1].
    threshold : float, optional
        The significance threshold in -log10 space (default: 1.0, corresponding to p=0.1).
        Similarities below this value will be set to 0.
    min_pval : float, optional
        Minimum p-value to substitute when p-value = 0 to avoid infinite values (default: 1e-10).

    Returns
    -------
    np.ndarray
        The similarity matrix (shape [n, n]) with:
        - Values representing -log10 transformed p-values
        - Diagonal elements set to 0
        - All values ≥ 0 and clipped to [0, -log10(min_pval)]
    """
    # Protect against log(0)
    protected_pvals = np.where(corrected_pvals == 0, min_pval, corrected_pvals)

    # Calculate similarity scores
    similarity = -np.log10(protected_pvals)

    # Apply threshold and clean up
    similarity[similarity <= threshold] = 0
    similarity = np.clip(similarity, 0, None)
    np.fill_diagonal(similarity, 0)

    return similarity


def adjust_similarity_matrix(similarity: np.ndarray) -> np.ndarray:
    """
    Normalizes a similarity matrix by row sums to create a probability transition matrix.

    Each element is divided by its row sum, converting similarities to transition probabilities.
    Handles zero-row sums gracefully by replacing NaN/Inf values with zeros.

    Parameters
    ----------
    similarity : np.ndarray
        Square non-negative similarity matrix (shape [n, n]).

    Returns
    -------
    np.ndarray
        Normalized transition probability matrix (shape [n, n]) where:
        - Each row sums to 1 (except zero rows which remain zero)
        - All values are in range [0, 1]
        - Diagonal elements may be non-zero

    """

    # Calculate row sums with protection against division by zero
    with np.errstate(divide="ignore", invalid="ignore"):
        row_sums = similarity.sum(axis=1)
        similarity_norm = similarity / row_sums[:, np.newaxis]

    # Handle NaN/Inf from zero rows
    similarity_norm = np.nan_to_num(similarity_norm, nan=0.0)

    return similarity_norm


def clusterize(similarity: np.ndarray, inflation: float) -> list[list[int]]:
    """
    Performs clustering using the Markov Cluster Algorithm (MCL).

    Parameters:
    -----------
    similarity : np.ndarray
        A similarity matrix representing pairwise similarities between elements.

    inflation : float
        The inflation parameter for MCL algorithm (typically 1.5-3.0).
        Higher values lead to more but smaller clusters.
        Lower values produce fewer but larger clusters.

    Returns:
    --------
    list[list[int]]
        A list of clusters, where each cluster is represented as a list of node indices.
        Empty clusters are automatically filtered out.

    Notes:
    ------
    - The function automatically converts input to sparse CSR format for efficiency.
    - Requires the markov_clustering package to be installed.

    """
    similarity_sparse = csr_matrix(similarity)

    # Run clusterization
    result = mc.run_mcl(similarity_sparse, inflation=inflation)
    modules = mc.get_clusters(result)

    print(f"Get {len(modules)} modules!")

    return modules


def assemble_results(
    modules: list[list[int]],
    cluster_sizes: pd.Series,
    cluster_names: pd.Index,  # Изменили аннотацию на pd.Index
    min_module_size: int = 1,
) -> pd.DataFrame:
    """
    Assembles clustering results into a structured DataFrame.

    Parameters:
    -----------
    modules : list[list[int]]
        List of modules where each module contains cluster indices
    cluster_sizes : pd.Series
        Series with cluster_id as index
    cluster_names : pd.Index
        Index containing cluster names (mapping from positions to cluster IDs)
    min_module_size : int, optional
        Minimum module size to include (default: 1)

    Returns:
    --------
    pd.DataFrame
        Results DataFrame with module information
    """
    # Preprocess data for faster access
    cluster_size_map = cluster_sizes.set_index("cluster_id")["cluster_size"].to_dict()

    results = []
    for module_id, cluster_indices in enumerate(modules):
        if len(cluster_indices) <= min_module_size:
            continue

        # Vectorized operations
        module_clusters = [cluster_names[idx] for idx in cluster_indices]
        nodes_sizes = [cluster_size_map[name] for name in module_clusters]

        results.append(
            {
                "module_id": module_id + 1,
                "module_size": len(cluster_indices),
                "cluster_ids": module_clusters,
                "cluster_inds": cluster_indices,
                "cluster_sizes": nodes_sizes,
            }
        )

    return pd.DataFrame(results)
