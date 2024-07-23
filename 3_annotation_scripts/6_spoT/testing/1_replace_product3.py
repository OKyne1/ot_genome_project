import os
import sys
from Bio import SeqIO

def update_product_names(file_path):
    # First, remove lines containing /gene="spoT"
    with open(file_path, "r") as file:
        lines = file.readlines()

    with open(file_path, "w") as file:
        for line in lines:
            if '/gene="spoT"' not in line:
                file.write(line)
            else:
                print(f'Removed line: {line.strip()}')

    # Parse the file again to update product names as needed
    records = []
    for record in SeqIO.parse(file_path, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if "product" in feature.qualifiers:
                    product = feature.qualifiers["product"][0]
                    translation = feature.qualifiers.get("translation", [""])[0]
                    gene_length_aa = len(translation)
                    gene_length_nt = (feature.location.end - feature.location.start)
                    
                    if "ppGpp" in product or "SpoT" in product:
                        if gene_length_aa > 1200 or gene_length_nt > 3600:  # 1200 amino acids or 3600 nucleotides
                            feature.qualifiers["product"] = ["Orientia SpoT-synthetase"]
                            new_product = "Orientia SpoT-synthetase"
                        else:
                            feature.qualifiers["product"] = ["truncated RelA/SpoT homolog (RSH) protein"]
                            new_product = "truncated RelA/SpoT homolog (RSH) protein"
                        
                        print(f"Gene length (aa): {gene_length_aa}, Gene length (nt): {gene_length_nt}, New product: {new_product}")
        records.append(record)
    
    with open(file_path, "w") as output_handle:
        SeqIO.write(records, output_handle, "genbank")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <gbff_file1> <gbff_file2> ...")
        sys.exit(1)
    
    for file_path in sys.argv[1:]:
        if not os.path.isfile(file_path):
            print(f"File not found: {file_path}")
            continue
        update_product_names(file_path)
        print(f"Updated file: {file_path}")
