# *O. tsutsugamushi* Genome Project

This is the repo for code relating to assembly and annotation of O. tsutsugamushi genomes from nanopore sequencing.

## Contents
- Assembly
- Annotation
- Analysis

## Assembly
Detailed information of assembly can be found [here](https://github.com/OKyne1/ot_genome_project/blob/main/1_assembly_scripts/optimised_assembly/README.md). The current approach involves:

1. Initial minimap2 filtering to remove host reads
2. Filtering out of shorter reads (<5000)
3. Flye assembly
4. Medaka polishing with all reads
5. Possibly Illumina polishing

### Assembly QC
The 3 C's of genome assembly:

1. Completeness
2. Contiguity
3. Correctness

#### Correctness
Traditionally BUSCO is the main method to test completeness. This is something we use in our analysis. However, unfortunately the regions that are hardest to assemble and polish are RAGEs. As these don't contain concerved single copy genes an assembly can score well on BUSCO but have low completeness. Additionally, RAGE regions are the hardest parts to polish and in our tests show the highest error rates (assemblies vs reference).

Thus, whilst this method is important it has its limitations. That said, it does correctly show the [Karp 2022 genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_022936085.1/) has lower levels of correctness due to many homopolymer errors.

#### Completeness
Assessment methods:

1. Is the genome closed?
2. Alignment to reference (for test assemblies)
3. Number and length of contigs

## Annotation
Obligate intracellular bacteria are often hard to annotate due to duplication of genes followed by trucation, pseudogenisation or degredation of genes. This makes it hard to identify the functional importance and identity of genes in *O. tsutsugamushi*. Additionally, as there is little or no recombination between strains there can be lots of sequence divergence, yet another problem in annotation of this bacterium.

### Gene Annotation
The annotation approach used for *O. tsutsugamushi* uses Bakta with user provided genes. The user provided genes are comprised of; heavily cleaned names from the Giengkam 2023 paper and additional annotations from NCBI for certain genes which annotated poorly and have low sequence homology between different copies of the gene.

### RAGE Annotation
