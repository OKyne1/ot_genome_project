import argparse
import os

def adjust_bed_file(input_bed):
    # Define the output file name
    output_bed = input_bed.replace('.bed', '_contiguous.bed')

    # Read all regions from the BED file
    regions = []
    with open(input_bed, 'r') as infile:
        for line in infile:
            chrom, start, end = line.strip().split()[:3]
            start, end = int(start), int(end)
            regions.append((chrom, start, end))

    # Function to generate new contiguous regions
    def generate_contiguous_regions(regions):
        if not regions:
            return []
        
        new_regions = []
        current_start = 0
        
        for chrom, start, end in regions:
            new_end = current_start + (end - start)
            new_regions.append((chrom, current_start, new_end))
            current_start = new_end
        
        return new_regions

    # Generate new contiguous regions
    new_regions = generate_contiguous_regions(regions)

    # Write the output BED file with the new format
    with open(output_bed, 'w') as outfile:
        for chrom, start, end in new_regions:
            outfile.write(f"{chrom}\t{start}\t{end}\n")
    
    print(f"Contiguous BED file saved to {output_bed}")

def main():
    parser = argparse.ArgumentParser(description="Adjust BED files to be contiguous while preserving original order.")
    parser.add_argument("bed_files", nargs='+', help="Paths to the input BED files")
    
    args = parser.parse_args()
    
    for bed_file in args.bed_files:
        if os.path.isfile(bed_file):
            adjust_bed_file(bed_file)
        else:
            print(f"File not found: {bed_file}")

if __name__ == "__main__":
    main()
