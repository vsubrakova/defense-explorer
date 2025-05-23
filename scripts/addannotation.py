import pandas as pd

# take GO annotation
import requests
from concurrent.futures import ThreadPoolExecutor

# Set annotation
def annotation_transformation(annotation: pd.DataFrame) -> pd.DataFrame:
    """
    Processes and transforms protein annotation DataFrame by:
     1. Renaming columns
     2. Cleaning GO terms
     3. Filtering by e-value
     4. Splitting GO terms into lists

    Parameters
    ----------
    annotation : pd.DataFrame
        Input annotation DataFrame containing protein sequence annotations.
        Expected columns:
        - 'seq_id': protein sequence identifier (will be renamed to 'protein_id')
        - 'GO': Gene Ontology terms (may contain '(InterPro)' tags and pipe-separated values)
        - 'e_value': significance value for filtering

    Returns
    -------
    pd.DataFrame
        Transformed annotation DataFrame with:
        - Renamed 'seq_id' column to 'protein_id'
        - Cleaned GO terms (removed '(InterPro)' tags and split by '|')
        - Only records with e_value < 1e-10
        - GO terms as lists of strings
    """
    annotation = annotation.rename(columns={"seq_id": "protein_id"})  # rename column

    # Point GO annotation type
    annotation["GO"] = annotation["GO"].astype(str)
    annotation.loc[(annotation.GO == "0") | (annotation.GO == "-"), "GO"] = "-"

    # Filtrate only significant records
    annotation = annotation[annotation.e_value < 1e-10]

    # Process GO terms
    annotation.GO = annotation.GO.apply(
        lambda x: x.replace("(InterPro)", "").split("|")
    )
    return annotation


def process_go_annotations(
    annotation: pd.DataFrame, max_workers: int = 10
) -> pd.DataFrame:
    """
    Processes GO annotations by fetching term names from QuickGO API in parallel.

    For each GO ID in the input DataFrame, retrieves the corresponding term name
    from EBI's QuickGO service, using caching to avoid duplicate requests and
    parallel processing for performance.

    Parameters
    ----------
    annotation : pd.DataFrame
        Input DataFrame containing protein annotations with GO terms.
        Must contain columns:
        - 'protein_id': Identifier for proteins
        - 'GO': List of GO terms (either as list or pipe-separated string)

    max_workers : int, optional
        Maximum number of threads to use for parallel API requests (default: 10)

    Returns
    -------
    pd.DataFrame
        Processed DataFrame with columns:
        - 'protein_id': Original protein identifiers
        - 'GO_terms': Corresponding GO term names (or '-' for missing terms)

    Notes
    -----
    - Requires internet connection to access EBI QuickGO API
    - Uses in-memory caching to avoid duplicate API requests
    - Empty/missing GO terms are represented as '-'
    - The function will explode list-type GO columns automatically
    """
    # Initialize cache for storing fetched GO terms
    go_term_cache = {}

    def get_go_term_cached(go_id: str) -> tuple:
        """
        Helper function to fetch GO term with caching.

        Args:
            go_id: GO identifier (e.g., 'GO:0003677')

        Returns:
            Tuple of (original_go_id, term_name) or (go_id, '-') for invalid terms
        """
        if go_id in go_term_cache:
            return (go_id, go_term_cache[go_id])
        if go_id == "-":
            return (go_id, "-")
        url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go_id}"
        response = requests.get(url).json()
        term = response["results"][0]["name"]
        go_term_cache[go_id] = term
        return (go_id, term)

    # Prepare data - keep only needed columns and explode GO lists
    annotation = annotation[["protein_id", "GO"]].copy()
    exploded = annotation.explode("GO")
    unique_go_ids = [go_id for go_id in exploded["GO"].unique()]

    go_terms = {}
    # Parallel fetching of GO terms using ThreadPool
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(get_go_term_cached, go_id) for go_id in unique_go_ids
        ]

        for future in futures:
            go_id, term = future.result()
            if term is not None:
                go_terms[go_id] = term

    # Map fetched terms back to exploded DataFrame
    exploded["GO_terms"] = exploded["GO"].map(go_terms).fillna(exploded["GO"])
    exploded = exploded.drop_duplicates()
    exploded = exploded.drop("GO", axis=1)
    return exploded


def process_pfam_annotations(annotation: pd.DataFrame) -> pd.DataFrame:
    """
    Processes and filters PFAM domain annotations from InterProScan results.

    Parameters
    ----------
    annotation : pd.DataFrame
        DataFrame containing InterProScan annotation results. Must include:
        - 'protein_id': Protein identifier column
        - 'analysis_type': Type of domain annotation (must contain 'Pfam')
        - 'signature_name': PFAM domain identifier/name

    Returns
    -------
    pd.DataFrame
        Processed DataFrame containing:
        - 'protein_id': Original protein identifiers
        - 'pfam_term': PFAM domain names/identifiers

    """
    pfam_annotation = annotation[annotation.analysis_type == "Pfam"].copy()
    pfam_annotation = pfam_annotation[["protein_id", "signature_name"]]
    pfam_annotation = pfam_annotation.drop_duplicates()
    pfam_annotation = pfam_annotation.rename(columns={"signature_name": "pfam_term"})
    return pfam_annotation


def transform_cluster_protein_data(data: pd.DataFrame) -> pd.DataFrame:
    """
    Transforms cluster-protein association data by:
    1. Extracting region information from protein IDs
    2. Creating one-hot encoded columns for each region type
    3. Applying standardized naming transformation to protein and cluster IDs

    Parameters
    ----------
    data : pd.DataFrame
        Input DataFrame containing protein-cluster associations with columns:
        - 'protein_id': Protein identifier (expected to contain region info in format '*_region_*')
        - 'cluster_id': Cluster identifier
        - Other columns will be preserved

    Returns
    -------
    pd.DataFrame
        Transformed DataFrame with additional columns:
        - 'region_[type]': One-hot encoded columns for each region type (upstream/downstream/immuneisland)
        - 'protein_id': Transformed protein identifiers
        - 'cluster_id': Transformed cluster identifiers
        - All original columns (except temporary ones)
    """
    # Excrete and calculate region name (upstream, downstream, immuneisland) from protein_id and cluster_id
    data["region"] = data["protein_id"].apply(lambda x: x.split("_")[-2])
    data = pd.get_dummies(data, columns=["region"])
    # Rename protein and clusters names
    data.loc[:, "protein_id"] = data["protein_id"].apply(name_transformation)
    data.loc[:, "cluster_id"] = data["cluster_id"].apply(name_transformation)
    return data


def name_transformation(name: str) -> str:
    """
    Delete from cluster and protein names parts protocluster and region
    """
    split_parts = name.split("_protocluster_", 1)
    last_number = split_parts[1].rsplit("_", 1)[1]
    name_transformed = split_parts[0] + "_" + last_number
    return name_transformed


def combine_annotations_data(
    data: pd.DataFrame, GO_annotation: pd.DataFrame, pfam_annotation: pd.DataFrame
) -> pd.DataFrame:
    """
    Combines and aggregates protein cluster data with GO and PFAM annotations.

    Processes and merges three data sources to create comprehensive cluster summaries:
    1. Calculates region composition for each cluster
    2. Merges and aggregates GO annotations
    3. Merges and aggregates PFAM annotations
    4. Computes final cluster statistics

    Parameters
    ----------
    data : pd.DataFrame
        Cluster-protein association data containing:
        - 'cluster_id': Cluster identifiers
        - 'protein_id': Protein identifiers
        - 'region_*' columns: One-hot encoded region types (downstream/immuneisland/upstream)

    GO_annotation : pd.DataFrame
        GO term annotations with columns:
        - 'protein_id': Protein identifiers
        - 'GO_terms': Associated GO terms

    pfam_annotation : pd.DataFrame
        PFAM domain annotations with columns:
        - 'protein_id': Protein identifiers
        - 'pfam_term': Associated PFAM domains

    Returns
    -------
    pd.DataFrame
        Aggregated cluster data with columns:
        - 'cluster_id': Cluster identifier
        - 'cluster_size': Number of unique proteins in cluster
        - 'GO_terms': Processed GO term annotations (from count_annotat)
        - 'pfam_term': Processed PFAM annotations (from count_annotat)
        - 'region': Dominant region type (from process_row)

    Notes
    -----
    - Uses helper functions:
      * process_row(): Determines dominant region type
      * count_annotat(): Processes annotation lists into summary format
    - Missing annotations are filled with '-'
    - Final output contains one row per unique cluster_id
    """
    region_sums = (
        data.drop_duplicates(["cluster_id", "protein_id"])
        .groupby("cluster_id")
        .agg(
            {
                "region_downstream": "sum",
                "region_immuneisland": "sum",
                "region_upstream": "sum",
            }
        )
        .reset_index()
    )
    region_sums["region"] = region_sums.apply(process_row, axis=1)
    # Delete unnecissary columns before merge
    region_sums = region_sums.drop(
        ["region_downstream", "region_immuneisland", "region_upstream"], axis=1
    )
    # Merge and delete columns from data
    data = data.drop(
        ["region_downstream", "region_immuneisland", "region_upstream"], axis=1
    )
    data = data.merge(region_sums, on="cluster_id", how="left")

    data_GO = data.merge(GO_annotation, on="protein_id", how="left")
    data_GO = data_GO.fillna("-")
    go_sums = (
        data_GO[["cluster_id", "protein_id", "GO_terms"]]
        .drop_duplicates()
        .drop("protein_id", axis=1)
    )
    go_sums = go_sums.groupby("cluster_id").agg({"GO_terms": list}).reset_index()
    go_sums.GO_terms = go_sums.GO_terms.apply(count_annotat)
    data_GO = data_GO.drop("GO_terms", axis=1)
    data_GO = data_GO.merge(go_sums, on="cluster_id", how="left")
    data_GO = data_GO.fillna("-")

    data_GO_pfam = data_GO.merge(pfam_annotation, on="protein_id", how="left")
    data_GO_pfam = data_GO_pfam.fillna("-")
    pfam_sums = (
        data_GO_pfam[["cluster_id", "protein_id", "pfam_term"]]
        .drop_duplicates()
        .drop("protein_id", axis=1)
    )
    pfam_sums = pfam_sums.groupby("cluster_id").agg({"pfam_term": list}).reset_index()
    pfam_sums.pfam_term = pfam_sums.pfam_term.apply(count_annotat)
    data_GO_pfam = data_GO_pfam.drop("pfam_term", axis=1)
    data_GO_pfam = data_GO_pfam.merge(pfam_sums, on="cluster_id", how="left")
    data_GO_pfam = data_GO_pfam.fillna("-")
    data_GO_pfam = (
        data_GO_pfam.groupby("cluster_id")
        .agg(
            {
                "protein_id": "nunique",
                "GO_terms": "first",
                "pfam_term": "first",
                "region": "first",
            }
        )
        .reset_index()
    )
    data_GO_pfam = data_GO_pfam.rename(columns={"protein_id": "cluster_size"})

    return data_GO_pfam


def process_row(row: pd.Series) -> str:
    """Formats region counts from a DataFrame row into a standardized string representation."""
    return f'downstream: {row["region_downstream"]}, immuneisland: {row["region_immuneisland"]}, upstream: {row["region_upstream"]}'


def count_annotat(annots: list):
    """
    Counts and formats biological annotations from a list into a summarized string."""
    count_dic = {}
    sentence = ""
    # Count occurrences of each annotation
    for annot in annots:
        if count_dic.get(annot):
            count_dic[annot] += 1
        else:
            count_dic[annot] = 1
    # Format non-empty annotations
    for annot, counters in count_dic.items():
        if annot != "-":
            sentence = sentence + f"{annot} : {counters}; "

    if sentence == "":
        sentence = "-"
    return sentence
