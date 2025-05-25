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

    """
    plt.figure(figsize=(4, 4))
    sns.histplot(clusters.cluster_id.value_counts(), bins=50, kde=False)
    plt.xlabel("Сluster sizes (number of proteins)", fontsize=10)
    plt.ylabel("Number of clusters (log scale)", fontsize=10)
    plt.yscale('log')
    plt.title("Distribution of PC sizes")
    plt.show()


def vizualize_heatmap(modules_annotation:pd.DataFrame)-> None:
   """
    Analyze immune gene distribution across protein modules and generate heatmap.
    
    Parameters:
    -----------
        modules_annotation: DataFrame with columns:
            - module_id: str/int - module identifiers
            - cluster_sizes: int - cluster sizes
            - immune: int - immune gene counts

    """
    # Aggregate immune data by module
    histogram_data = modules_annotation.groupby('module_id').agg(
        total_cluster_size=('cluster_sizes', 'sum'),
        total_immune=('immune', 'first'),
        #protein_ids=('mmseq_protein_ids', lambda x: sum([p for p in x if isinstance(p, list)], []))
    ).reset_index()

    histogram_data['immune_to_total_ratio'] = histogram_data['total_immune'] / histogram_data['total_cluster_size']

    # Create bins for module size
    histogram_data['size_bin'] = pd.cut(
        histogram_data['total_cluster_size'],
        bins=[0, 50, 100, 150, 500000],
        labels=['0-49','50-99', '100-149', '150+'],
        right=False
    )

    # Heatmap data pivot table
    heatmap_data = histogram_data.pivot_table(
        index='size_bin',
        columns=pd.cut(histogram_data['immune_to_total_ratio'], bins=5),
        values='total_cluster_size',
        aggfunc='sum',
        fill_value=0,
        observed=False
    )

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
    

def vizualize_module_distribution(modules: pd.DataFrame):
    """
    Plot histogram showing distribution of protein cluster sizes across modules.
    
    Parameters:
    -----------
        modules: DataFrame containing module data with:
            - module_size: int - size of each module (number of protein clusters)
            
    """
    plt.figure(figsize=(3, 3))
    bin_edges = np.arange(modules['module_size'].min() - 0.5, modules['module_size'].max() + 1.5, 1)
    counts, bins, patches = plt.hist(modules['module_size'], bins=bin_edges, edgecolor='black')
    
    plt.xticks(np.arange(modules['module_size'].min(), modules['module_size'].max() + 1))
    
    plt.title('Distribution of Module Sizes', fontsize=14)
    plt.xlabel('Module size (# of protein clusters)', fontsize=12)
    plt.ylabel('Number of modules', fontsize=12)
    
    plt.tight_layout()
    plt.show()
    