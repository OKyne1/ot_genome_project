import argparse
from Bio import SeqIO

def update_product_from_txt(genbank_file, txt_file, domain):
    # Read locus tags and product descriptions from the txt file into a dictionary
    locus_to_product = {}
    with open(txt_file, 'r') as txt_handle:
        for line in txt_handle:
            locus_tag, product = line.strip().split('\t')
            locus_to_product[locus_tag] = product

    # Update product descriptions and remove gene names in the GenBank file
    with open(genbank_file, 'r') as gbff_handle:
        records = list(SeqIO.parse(gbff_handle, 'genbank'))

    for record in records:
        for feature in record.features:
            if feature.type == 'gene' and 'gene' in feature.qualifiers:
                gene_names = feature.qualifiers['gene']
                filtered_gene_names = [name for name in gene_names if name.lower() in ('trag', 'bamd')]
                if filtered_gene_names:
                    feature.qualifiers['gene'] = filtered_gene_names
                else:
                    del feature.qualifiers['gene']
            elif feature.type == 'CDS' and 'locus_tag' in feature.qualifiers and feature.qualifiers['locus_tag'][0] in locus_to_product:
                gene_names = feature.qualifiers.get('gene', [])
                gene_names_lower = [name.lower() for name in gene_names] if gene_names else []
                if domain == 'tpr' and ('trag' in gene_names_lower or 'bamd' in gene_names_lower):
                    continue  # Skip modification for traG or bamD in the 'CDS' section
                filtered_gene_names = [name for name in gene_names if name.lower() in ('trag', 'bamd')]
                if filtered_gene_names:
                    feature.qualifiers['gene'] = filtered_gene_names
                else:
                    if 'gene' in feature.qualifiers:
                        del feature.qualifiers['gene']  # Remove gene name if not 'trag' or 'bamd'

                locus_tag = feature.qualifiers['locus_tag'][0]
                if locus_tag in locus_to_product:
                    if 'product' in feature.qualifiers:  # Check if product is already set
                        # Replace the existing product with the new one
                        feature.qualifiers['product'] = [locus_to_product[locus_tag]]
                    else:
                        # Set the product if it's not already present
                        feature.qualifiers['product'] = [locus_to_product[locus_tag]]

    # Write the updated records back to the GenBank file
    with open(genbank_file, 'w') as gbff_handle:
        SeqIO.write(records, gbff_handle, 'genbank')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update product descriptions in a GenBank file from a text file.")
    parser.add_argument("genbank_file", help="Input GenBank file path")
    parser.add_argument("txt_file", help="Input text file with locus tags and product descriptions")
    parser.add_argument("-dom", "--domain", choices=['ank', 'tpr'], default='ank', help="Domain type")
    args = parser.parse_args()

    update_product_from_txt(args.genbank_file, args.txt_file, args.domain)
    print("Product descriptions and gene names updated successfully.")
