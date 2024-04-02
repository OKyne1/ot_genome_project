# Annotation of RAGE regions
Structure:
1. Defining RAGE derived regions
2. Defining boundaries for complete RAGEs
3. Validating boundaries within RAGE derived regions
4. Identifying complete RAGEs based on required genes

## RAGE derived regions
These are currently defined by a list of different 

## Complete RAGE boundaries
### Rules for boundaries are:
![RAGE boundaries](https://github.com/OKyne1/ot_genome_project/blob/main/2_annotation_scripts/3_rage_classification/rage_boundaries_conditions.png)
The black lines mark the boundary, if a gene is not included in the boundary then the script will restart on this line to identify additional boundaries using this gene.

### Edge case currently not covered:
![Edge case](https://github.com/OKyne1/ot_genome_project/blob/main/2_annotation_scripts/3_rage_classification/edge_case.png)
Currently, the script doesn't cover this rare edge case where there are 2 integrases 'facing' eachother with a dnaA gene sandwiched between them. The script will only accept the first RAGE region and not the second.

## Validating boundaries

## Testing for required genes
