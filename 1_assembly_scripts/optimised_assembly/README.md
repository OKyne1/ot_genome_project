# Pipeline for Assembly of Highly Repetitive Intracellular Bacteria Using Nanopore Reads
## Read Filtering
## Initial Optimisation Tests
### Different Assemblers
Canu and flye were both tested. In this analysis flye performed significantly better and required ~60x fewer CPU hours. Consequently, in the future all assemblies will be performed using flye.

### Flye Assembly Optimisation
#### Testing --min-overlap Parameter
This parameter was tested between 1000 and 10,000 at 1000 bp intervals for data from barcode 05 and the 3 combined barcodes (from the 2024 karp_2 sequencing run). Generally, the number of contigs reduces as the min-overlap increases and the longest contig seems to be optimal for around 7000 or 8000 bp. However, the performance of this was worse than other assembly attempts tested (below) and so this approach was dropped.

## 1. Read filtering

## 2. Pre-assembly filtering and flye assembly

## 3. Medaka polishing with ONT reads

## 4. Illumina polishing (Optional)

