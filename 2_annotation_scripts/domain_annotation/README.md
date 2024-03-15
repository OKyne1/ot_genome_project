# Domain Annotation
Some important domains (ankyrin and tetratricopeptide repeats) are not annotated as thougherly or completely by bakta as desired. Consequently we are implementing an alternative approach for the annotation of these features.

## Structure of Domain Annotation Code
1. Parsing gbff files from bakta to faa files containing locus_tag and sequence information
2. hmmersearch to identify the presence of these domains
3. Parsing hmmer outputs to give the domain type and counts
4. Adding this information back into the gbff files
