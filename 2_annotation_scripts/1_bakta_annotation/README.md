# Bakta Annotation with Curated Proteins
This approach uses bakta, a rapid genome annotation tool (which doesn't require species specific databases) and enables user provided genes.

## Curated Proteins
1. Gene annotations from Giengkam et al. 2023, with truncated, pseudo or degraded proteins removed. Also all domain based annotations were removed (e.g. Ank and tpr) to priorise more specific naming
2. Genes from NCBI for increased coverage.

### NCBI Proteins
*Orientia tsutsugamushi* protein sequences from NCBI: cinA, scaA, scaB, scaD, scaE, secA, traA, traG, tsa22, tsa47 and tsa56.

These were selected based on poor annotation with bakta with the Giengkam annotations (in numerous tests). Truncated or incomplete genes were excluded in the search.

The minimum identity was reduced from 90% to 80% for some genes; scaA, scaC, tsa22 and tsa56. This was deemed reasonable due to the high sequence varaibility for these genes likely due to selection pressure for diversity.
