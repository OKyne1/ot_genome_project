# *O. tsutsugamushi* Genome Project

This is the repo for code relating to assembly and annotation of O. tsutsugamushi genomes from nanopore sequencing.

## Contents
- Assembly
- Annotation
- Analysis

## Assembly
This section contains;
1. *O. tsutsugamushi* optimised assembly scripts
2. Assembly quality control steps
3. Other, a collection of scripts that are useful, but not essential to the assembly pipeline (or were used in optimisation)

### Optimised assembly scripts
Tests performed to optimise:
1. Different assemblers (canu and flye)
2. flye with different assembly coverage
3. Prefiltering of reads for different coverage or lengths
4. Nanosim simulation of data and then testing the above prefiltering
### Assembly QC

### Other stuff

## Annotation
Obligate intracellular bacteria are often hard to annotate due to duplication of genes followed by trucation, pseudogenisation or degredation of genes. This makes it hard to identify the functional importance and identity of genes in *O. tsutsugamushi*. Additionally, as there is little or no recombination between strains there can be lots of sequence divergence, yet another problem in annotation of this bacterium.

### Gene Annotation
The annotation approach used for *O. tsutsugamushi* uses Bakta with user provided genes. The user provided genes are comprised of; heavily cleaned names from the Giengkam 2023 paper and additional annotations from NCBI for certain genes which annotated poorly and have low sequence homology between different copies of the gene.

### RAGE Annotation
