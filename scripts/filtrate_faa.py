from Bio import SeqIO

# open file with representative headers
with open("./headers.txt", "r") as f:
    target_headers = {line.strip().replace(">", "") for line in f}

# open output file for writing (proteins presented in headers)
with open("concat_filtrated.faa", "w") as out_file:
    # open file with proteins
    for record in SeqIO.parse("./concat.faa", "fasta"):
        if record.description in target_headers:
            SeqIO.write(record, out_file, "fasta")

print(f"Ready! Filtered protein seqs in concat_filtrated.faa")
