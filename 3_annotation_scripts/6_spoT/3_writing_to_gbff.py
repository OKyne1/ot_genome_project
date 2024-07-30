from Bio import SeqIO
import os

def process_gbff(input_file, locus_tag_file):
    # Read locus tags and values from the text file
    locus_tags_to_modify = {}
    with open(locus_tag_file, "r") as f:
        for line in f:
            locus_tag, percent_aligned, completeness, name = line.strip().split("|")
            if completeness == "complete":
                name_suffix = name.split("_")[-1]
                locus_tags_to_modify[locus_tag] = name_suffix

    with open(input_file, "r") as f:
        records = list(SeqIO.parse(f, "genbank"))

    for record in records:
        for feature in record.features:
            if "locus_tag" in feature.qualifiers:
                locus_tag = feature.qualifiers["locus_tag"][0]
                if locus_tag in locus_tags_to_modify:
                    if feature.type == "CDS":
                        if "product" in feature.qualifiers:
                            name_suffix = locus_tags_to_modify[locus_tag]
                            feature.qualifiers["product"] = [f"full length bifunctional spoT"]
                        feature.qualifiers["gene"] = ["spoT"]
                    if feature.type == "gene":
                        feature.qualifiers["gene"] = ["spoT"]

    output_filename = os.path.splitext(input_file)[0] + "_spotted.gbff"
    with open(output_filename, "w") as f:
        SeqIO.write(records, f, "genbank")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Overwrite product description and add /gene='spoT' to genes in a GenBank file.")
    parser.add_argument("input_file", help="Input GenBank file")
    parser.add_argument("locus_tag_file", help="Text file containing locus tags")
    args = parser.parse_args()

    process_gbff(args.input_file, args.locus_tag_file)