# Annotation of RAGE regions
## Code structure 
1. Defining RAGE derived regions
2. Defining boundaries for complete RAGEs
3. Validating boundaries within RAGE derived regions
4. Identifying complete RAGEs based on required genes

![workflow](./diagrams/code_workflow.png)

These stages are tied together in the script `main.sh` so that processing of genbank files can be performed in one stage to define the RAGEs derived regions and complete RAGEs. Output files .bed files.

### Code usage
1. Create and activate the rage conda environment. This can be created from `rage_environment.yml`.
2. The script can be run like this (locally, not adjusted for the cluster):  `bash ../main.sh <gbff_file> <gbff_file> ... <gbff_file>`. 

## RAGE derived regions
These are currently defined by a list of different rage genes. These were determined by identifying genes in the bakta annotation which were considered in RAGE regions by Jeanne's paper. 

This list of RAGE genes was generated from the Karp, Gilliam and Boryong genomes, then tested on other genomes to validate the success.

### Visualising success
![kato rage derived regions](./diagrams/kato_rage_derived.png)
This IGV image of the different boundaries shows that my script closely resembles the manual annotations (note manual annotations have additional breaks for boundaries between RAGEs). The 2 lines below the manual annotated boundaries shows the effects of allowing one or 2 mis-matched genes within a RAGE region.

The effect here is not that large and so the current script only allows 1 non-match. However it can easily be modified to permit 2. The code was developed to permit this, and also not to allow these non-matches to be adjacent to eachother. 

### Rules for the RAGE derived regions
1. Needs to be more than 2 genes next to eachother
2. A region can contain a single non-matched gene (but not at the edges)
3. Genes must match the `RAGE_gene_list` but they cannot match to the `exclusion_list` (this was generated so that genes which match the list but contain more information and should be excluded)


## Complete RAGE boundaries
### Rules for boundaries are:
<img src="https://github.com/OKyne1/ot_genome_project/blob/main/2_annotation_scripts/3_rage_classification/diagrams/rage_boundaries_conditions.png" width="500">
The black lines mark the boundary, if a gene is not included in the boundary then the script will restart on this line to identify additional boundaries using this gene.

### Edge case:
<img src="https://github.com/OKyne1/ot_genome_project/blob/main/2_annotation_scripts/4_rage_classification/diagrams/edge_case.png" width="800">

The script currently allows 2 RAGE regions to share a dnaA gene. This is unlikely to cause any issues as the chance of 2 complete RAGEs next to eachother is extremely unlikely and not problematic.


## Validating boundaries
Boundaries identified are validated against the RAGE derived regions. Only those contained within this .bed file are kept for required gene testing. This step uses the package bed tools.


## Testing for required genes
To be a complete RAGE, it must have the boundries (already identified), all tra genes, at least 1 transposase and at least 1 cargo gene (currently using those defined in the paper). The presence of these genes is tested on the validated boundaries and the regions which meet these criteria are outputted in a new .bed file.


## Initial results
These scripts identify **100%** of complete RAGEs (3/3) but it also identifies **17/52** of the complete RAGEs with truncated genes. This is probably better than underclassification, however, we need to find a way to identify which of these are not truncated.


## RAGE gene completeness
Overclassification of RAGEs, is a result of truncated or degraded genes not being identified. Using a Blast+ script 

Method:
1. Blast+ search of the genomes genes against a database of complete tra genes and complete integrases
2. Testing the percentage of the genome aligned

This was successful in identifying >95% of all cases (when testing 7 genomes using a db of the 8th), consequently the coverage will be increased when the all 8 are used in the db.

## Current limitations
- Need to make it work for contigs
- Handle the identification of complete RAGEs with truncated genes?
- Bed files are currently indexed wrong. Need to -1 from all start entries (not done yet as it makes comparisons easier)
- Need to add in the information from the completeness checks (and trak1 and trak2)
- Add complete, as a requirement for the complete RAGEs
