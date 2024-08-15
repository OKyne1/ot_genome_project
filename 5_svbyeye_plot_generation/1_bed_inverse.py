import os
import tempfile
import subprocess
import argparse
from Bio import SeqIO

def generate_genome_file(fasta_file):
    """
    Generate a temporary genome file from a FASTA genome file.

    :param fasta_file: Path to the input FASTA genome file.
    :return: Path to the temporary genome file.
    """
    # Create a temporary file
    tmp_genome_file = tempfile.NamedTemporaryFile(delete=False, mode='w')
    
    try:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            # Write the contig/chromosome name and its length
            tmp_genome_file.write(f"{record.id}\t{len(record.seq)}\n")
        
        tmp_genome_file.close()
        return tmp_genome_file.name
    
    except Exception as e:
        print(f"An error occurred while generating genome file: {e}")
        tmp_genome_file.close()
        os.remove(tmp_genome_file.name)
        return None

def inverse_bed(fasta_file, input_bed, output_bed):
    """
    Inverse a BED file using bedtools complement with a temporary genome file.

    :param fasta_file: Path to the input FASTA genome file.
    :param input_bed: Path to the input BED file.
    :param output_bed: Path to the output (inverted) BED file.
    """
    # Generate the genome file
    genome_file = generate_genome_file(fasta_file)
    
    if genome_file is None:
        print("Failed to generate genome file.")
        return
    
    try:
        # Construct the bedtools complement command
        command = f"bedtools complement -i {input_bed} -g {genome_file}"
        
        # Execute the command and capture the output
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        
        # Check for errors
        if result.returncode != 0:
            print(f"Error running bedtools complement: {result.stderr}")
            return
        
        # Write the output to the output_bed file
        with open(output_bed, 'w') as out_file:
            out_file.write(result.stdout)
        
        print(f"Inverted BED file saved to {output_bed}")
    
    except Exception as e:
        print(f"An error occurred: {e}")
    
    finally:
        # Clean up the temporary genome file
        os.remove(genome_file)

def main():
    parser = argparse.ArgumentParser(description="Generate an inverse BED file using bedtools complement.")
    parser.add_argument("genome", help="Path to the input FASTA genome file")
    parser.add_argument("bed", help="Path to the input BED file")
    parser.add_argument("output", help="Path to the output (inverted) BED file")
    
    args = parser.parse_args()
    
    inverse_bed(args.genome, args.bed, args.output)

if __name__ == "__main__":
    main()
