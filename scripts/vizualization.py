import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def vizualize_region_distribution(clusters: pd.DataFrame) -> None:
    """
    Visualize the distribution of protein regions as a pie chart.
    
    The pie chart shows the proportion of proteins belonging to upstream, 
    immuneisland, and downstream regions based on their protein_id values.
    
    Parameters:
    -----------
    clusters : pd.DataFrame
        DataFrame containing protein cluster information. Must contain a 'protein_id' column
        with strings that can be matched against the region names.
        
    Returns:
    --------
    None
        Displays a matplotlib pie chart.
    """
    regions = ["upstream", "immuneisland", "downstream"]
    counts = {region: clusters.protein_id.str.contains(region, case=False, regex=True).sum() 
               for region in regions}
    plt.figure(figsize=(4, 4))
    colors = ['#9ecae1', '#6baed6', '#3182bd']  # Blue colours
    plt.pie(counts.values(),
        labels=counts.keys(),
        colors=colors,
        autopct='%1.1f%%',
        textprops={'fontsize': 15})
    plt.title("Distribution of protein's regions")
    plt.tight_layout()
    plt.show()

def vizualize_size_distribution(clusters: pd.DataFrame) -> None:
    """
    Visualize the distribution of cluster sizes as a histogram with log scale y-axis.
    
    The histogram shows how many clusters exist for each size category (number of proteins
    per cluster). The y-axis is log-scaled to better visualize the distribution.
    
    Parameters:
    -----------
    clusters : pd.DataFrame
        DataFrame containing protein cluster information. Must contain a 'cluster_id' column
        used to count cluster sizes.
        
    Returns:
    --------
    None
        Displays a matplotlib histogram.
    """
    plt.figure(figsize=(4, 4))
    sns.histplot(clusters.cluster_id.value_counts(), bins=50, kde=False)
    plt.xlabel("Сluster sizes (number of proteins)", fontsize=10)
    plt.ylabel("Number of clusters (log scale)", fontsize=10)
    plt.yscale('log')
    plt.title("Distribution of PC sizes")
    plt.show()

