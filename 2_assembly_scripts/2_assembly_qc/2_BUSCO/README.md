# BUSCO

BUSCO (Benchmarking Universal Single-Copy Orthologs) is a bioinformatics tool used to assess the completeness of genome assemblies by evaluating the presence of conserved single-copy genes. It provides a quantitative measure of genome quality by reporting the proportion of these orthologs that are complete, fragmented, or missing.

Unfortunately for _Orientia_ genomes, even genomes present in many contigs they will normally still look complete in BUSCO. Additionally for the database used, _Orientia_ doesn't contain 100% of the genes in the rickettsiales database.

Karp results in a BUSCO output with 92.9% of the genes complete (not the 100%) that would normally be expected (as a result of the genome reduction observed in intracellular bacteria). Here are the full details for the complete Karp genome:

| Input_file                        | Dataset             | Complete | Single | Duplicated | Fragmented | Missing | n_markers | Scaffold N50 | Contigs N50 | Percent gaps | Number of scaffolds |
|-----------------------------------|---------------------|----------|--------|------------|------------|---------|-----------|--------------|-------------|--------------|--------------------|
| Ot_karp_SAMEA104570320.fna        | rickettsiales_odb10 | 92.9     | 92.6   | 0.3        | 0.8        | 6.3     | 364       | 2469803      | 2469803     | 0.000%       | 1                  |


Usage:
`busco -i genomes/ -l rickettsiales_odb10 -o busco_output -m genome`

The outputs can then be plotted using the `generate_plot.py` which is part of the BUSCO package.

BUSCO will need to be present in the active environment to run this code. It also requires internet signal and so cannot be submitted as a job on the computer cluster.