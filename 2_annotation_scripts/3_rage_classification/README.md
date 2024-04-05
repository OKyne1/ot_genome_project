# Annotation of RAGE regions
## Code structure 
1. Defining RAGE derived regions
2. Defining boundaries for complete RAGEs
3. Validating boundaries within RAGE derived regions
4. Identifying complete RAGEs based on required genes

These stages are tied together in the script `main.sh` so that processing of genbank files can be performed in one stage to define the RAGEs derived regions and complete RAGEs.

## RAGE derived regions
These are currently defined by a list of different rage genes. These were determined by identifying genes in the bakta annotation which were considered in RAGE regions by Jeanne's paper. 

This list of RAGE genes was generated from the Karp, Gilliam and Boryong genomes, then tested on some other genomes to validate the success.

### Visualising success
![kato rage derived regions](https://github.com/OKyne1/ot_genome_project/blob/main/2_annotation_scripts/3_rage_classification/diagrams/kato_rage_derived.png)
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

### Edge case currently not covered:
![Edge case](https://github.com/OKyne1/ot_genome_project/blob/main/2_annotation_scripts/3_rage_classification/diagrams/edge_case.png)
Currently, the script doesn't cover this rare edge case where there are 2 integrases 'facing' eachother with a dnaA gene sandwiched between them. The script will only accept the first RAGE region and not the second.

This likely doesn't will only have a very minor influence. However, it may be worth considering in the future.

### Modifications to make
1. Make the code handle contigs
2. Do I want to consider the edge case
3. Break the code into more readable chunks
4. It will also currently recognise a single integrase when it is next to another one (-->-->)

## Validating boundaries
Boundaries identified are validated against the RAGE derived regions. Only those contained within this .bed file are kept for required gene testing. This step uses the package bed tools.

## Testing for required genes
To be a complete RAGE, it must have the boundries (already identified), all tra genes, at least 1 transposase and at least 1 cargo gene (currently using those defined in the paper).

## Initial results
These scripts identify 100% of complete RAGEs (3/3) but it also identifies 17/52 of the complete RAGEs with truncated genes.

**Next steps**:
Investigate whether these complete RAGEs with truncated genes contain pseudo or truncated descriptions. If they do this may be a method to exclude them.

## Current limitations
- Edge case mentioned above
- Need to make it work for contigs
- Need to tie the scripts together
- May want to make the scripts more readable/shorter by breaking them into smaller functions/chunks
- Handle the identification of complete RAGEs with truncated genes
