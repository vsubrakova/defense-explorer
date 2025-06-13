import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def split_gbk_by_protocluster_faa(filepath, output_dir):
    """
    Splits a GenBank file into multiple FASTA files based on protocluster IDs and functional CDS features, splitting them into upstream, immune island, and downstream regions.

    Parameters:
    -----------
    filepath : str
        Path to the GenBank file to be processed.
    output_dir : str
        Directory where the resulting FASTA files will be saved.

    Returns:
    --------
    None
    """

    os.makedirs(output_dir, exist_ok=True)

    basename = os.path.splitext(os.path.basename(filepath))[0]
    sample_name = basename.split("_")[0]

    records = list(SeqIO.parse(filepath, "genbank"))

    clusters = {}

    # Iterate over each record and its features to identify protoclusters
    for record in records:
        for feature in record.features:
            if feature.type == "CDS" and "protocluster_id" in feature.qualifiers:
                cluster_id = feature.qualifiers["protocluster_id"][0]
                if cluster_id not in clusters:
                    clusters[cluster_id] = []
                clusters[cluster_id].append(feature)

    # Process each cluster found in the GenBank records
    for cluster_id, cds_list in clusters.items():
        # Find immune island boundaries based on functional annotations
        indices_with_function = [
            i
            for i, feat in enumerate(cds_list)
            if "function" in feat.qualifiers or "gene_functions" in feat.qualifiers
        ]

        # Skip clusters with no functional CDS features
        if not indices_with_function:
            print(f"Protocluster {cluster_id}: no functional CDS found - skipping.")
            continue

        # Identify the indices of the first and last functional CDS
        first_idx = indices_with_function[0]
        last_idx = indices_with_function[-1]

        # Divide the CDS list into upstream, immune island, and downstream regions
        regions = {
            "upstream": cds_list[:first_idx],
            "immuneisland": cds_list[first_idx : last_idx + 1],
            "downstream": cds_list[last_idx + 1 :],
        }

        # Process each region and save as a separate FASTA file
        for part, feats in regions.items():
            if not feats:
                continue

            out_file = os.path.join(
                output_dir, f"{basename}_protocluster_{cluster_id}_{part}.fasta"
            )

            seq_records = []
            for feat in feats:
                try:
                    aa_seq = feat.qualifiers["translation"][0]

                    # Get the contig and protein ID for the header
                    contig_id = feat.qualifiers.get("ID", ["unknown"])[0]
                    protein_id = contig_id.split("_")[-1]
                    header = f"{basename}_protocluster_{cluster_id}_{part}_{protein_id}"

                    # Create a SeqRecord and append to the list
                    seq_record = SeqRecord(Seq(aa_seq), id=header, description="")
                    seq_records.append(seq_record)
                except IndexError:
                    print(f"Empty 'translation' in protocluster {cluster_id}, region {part} — skipping.")

            if not seq_records:
                continue

            with open(out_file, "w") as f:
                SeqIO.write(seq_records, f, "fasta")
            # print(f"Saved: {out_file}")


if __name__ == "__main__":
    """
    Main function that takes command-line arguments for input and output paths,
    and calls split_gbk_by_protocluster to process the GenBank file.
    """
    if len(sys.argv) != 3:
        print("Usage: python split.py input.gbk output_folder/")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    split_gbk_by_protocluster_faa(input_path, output_path)
