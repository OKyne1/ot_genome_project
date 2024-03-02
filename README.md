# *O. tsutsugamushi* Genome Project

This is the repo for code relating to assembly and annotation of O. tsutsugamushi genomes from nanopore sequencing.

## Contents
- Assembly
- Annotation
- Analysis

## Assembly
Detailed information of assembly can be found [here](https://github.com/OKyne1/ot_genome_project/blob/main/1_assembly_scripts/optimised_assembly/README.md). The current approach involves;
1. Initial minimap2 filtering to remove host reads
2. Filtering out of shorter reads (<5000)
3. Flye assembly
4. Medaka polishing with all reads
5. Possibly Illumina polishing

## Annotation
Obligate intracellular bacteria are often hard to annotate due to duplication of genes followed by trucation, pseudogenisation or degredation of genes. This makes it hard to identify the functional importance and identity of genes in *O. tsutsugamushi*. Additionally, as there is little or no recombination between strains there can be lots of sequence divergence, yet another problem in annotation of this bacterium.

### Gene Annotation
The annotation approach used for *O. tsutsugamushi* uses Bakta with user provided genes. The user provided genes are comprised of; heavily cleaned names from the Giengkam 2023 paper and additional annotations from NCBI for certain genes which annotated poorly and have low sequence homology between different copies of the gene.

### RAGE Annotation
