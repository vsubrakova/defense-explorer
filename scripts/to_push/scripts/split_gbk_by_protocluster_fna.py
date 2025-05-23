import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extract_nt_sequence(GenBank_file, feature):
    """
    Extracts nucleotide sequence from a GenBank record using feature coordinates.

    Parameters:
    -----------
    GenBank_file : GenBank file with the sequence data.
    feature : The feature within the record from which the sequence will be extracted. Typically a CDS feature with coordinates specifying the start and end of the sequence.

    Returns:
    --------
    tuple : A tuple containing the nucleotide sequence as a string, start position, and end position. 
    """
    try:
        start = int(feature.location.start) + 1  # Convert to 1-based indexing
        end = int(feature.location.end)
        return str(GenBank_file.seq[start-1:end]), start, end  # Python uses 0-based indexing
    except Exception as e:
        print(f"Error extracting sequence: {e}")
        return None, None, None

def get_protein_id(feature):
    """
    Extract the protein ID from a GenBank feature.
    
    Parameters:
    -----------
    feature : SeqFeature
        The feature from which the protein ID will be extracted. Typically a CDS feature.

    Returns:
    --------
    str
        The protein ID, or "unknown" if it cannot be found.
    """
    if "protein_id" in feature.qualifiers:
        return feature.qualifiers["protein_id"][0]
    
    cds_id = feature.qualifiers.get("ID", ["unknown"])[0]
    parts = cds_id.split('_')
    return parts[-1] if len(parts) >= 3 else "unknown"

def split_gbk_by_protocluster_fna(filepath, output_dir):
    """
    Splits a GenBank file by protoclusters and extracts nucleotide sequences of CDS features 
    for each protocluster, organizing them into functional regions (upstream, immune island, downstream).

    The resulting sequences are saved in separate FASTA files, one for each protocluster region.

    Parameters:
    -----------
    filepath : str
        Path to the GenBank file to be processed.
    output_dir : str
        Directory where the resulting FASTA files will be saved. This directory is created if it doesn't exist.

    Returns:
    --------
    None
        This function doesn't return anything but saves the sequences into FASTA files for each protocluster.
    """
    os.makedirs(output_dir, exist_ok=True)
    basename = os.path.splitext(os.path.basename(filepath))[0]
    records = list(SeqIO.parse(filepath, "genbank"))

    # Group CDS features by their protocluster ID
    clusters = {}
    for record in records:
        for feature in record.features:
            if feature.type == "CDS" and "protocluster_id" in feature.qualifiers:
                cluster_id = feature.qualifiers["protocluster_id"][0]
                clusters.setdefault(cluster_id, []).append((record, feature))

    # Process each protocluster
    for cluster_id, record_features in clusters.items():
        # Find functional CDS features by checking for the presence of functional annotations
        functional_indices = [
            i for i, (_, feat) in enumerate(record_features)
            if "function" in feat.qualifiers or "gene_functions" in feat.qualifiers
        ]

        if not functional_indices:
            # Skip protoclusters with no functional CDS features
            print(f"Skipping protocluster {cluster_id} - no functional CDS features found")
            continue

        # Define regions (upstream, immune island, downstream)
        first_idx = functional_indices[0]
        last_idx = functional_indices[-1]
        regions = {
            "upstream": record_features[:first_idx],
            "immuneisland": record_features[first_idx:last_idx+1],
            "downstream": record_features[last_idx+1:]
        }

        # Process each region and write to FASTA files
        for region_name, features in regions.items():
            if not features:
                continue

            out_file = os.path.join(output_dir, f"{basename}_protocluster_{cluster_id}_{region_name}.fna")

            seq_records = []
            for record, feat in features:
                nt_seq, start, end = extract_nt_sequence(record, feat)
                if nt_seq:
                    protein_id = get_protein_id(feat)
                    header = f"{basename}_protocluster_{cluster_id}_{region_name}_{protein_id}"
                    seq_records.append(SeqRecord(
                        Seq(nt_seq),
                        id=header,
                        description=f""
                    ))

            if seq_records:
                SeqIO.write(seq_records, out_file, "fasta")
                #print(f"Saved: {out_file}")

if __name__ == "__main__":
    """
    Main function that takes command-line arguments for input and output paths,
    and calls split_gbk_by_protocluster_fna to process the GenBank file.
    """
    if len(sys.argv) != 3:
        print("Usage: python split_nt.py input.gbk output_folder/")
        sys.exit(1)
    
    split_gbk_by_protocluster_fna(sys.argv[1], sys.argv[2])

