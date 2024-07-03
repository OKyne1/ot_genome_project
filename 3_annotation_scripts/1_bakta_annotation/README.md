# Bakta Annotation with Curated Proteins
This approach uses bakta, a rapid genome annotation tool (which doesn't require species specific databases) and enables user provided genes.

## Curated Proteins
1. Gene annotations from Giengkam et al. 2023
2. Genes from NCBI for increased coverage.

### Manual Annotations
Substantial cleaning of the dataset from Giengkam et al. 2023 was performed to create this database. This involved the removal of any truncated, pseudo or degraded sequences as well as removal of domain based annotations (e.g. ank and tprs) as these could result in less specific annotations.

### NCBI Proteins
*Orientia tsutsugamushi* protein sequences from NCBI: cinA, scaA, scaB, scaD, scaE, secA, traA, traG, tsa22, tsa47 and tsa56.

These were selected based on poor annotation with bakta with the Giengkam annotations (in numerous tests). Truncated or incomplete genes were excluded in the search.

The minimum identity was reduced from 90% to 80% for some genes; scaA, scaC, tsa22 and tsa56. This was deemed reasonable due to the high sequence varaibility for these genes likely due to selection pressure for diversity.

### New Names
As bakta requires a gene name and the gene transfer agents identified in Orientia, do not have names, I called these otgtaA-G in the file supplied to bakta.
