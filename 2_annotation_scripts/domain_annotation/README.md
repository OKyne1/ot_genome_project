# Domain Annotation
Some important domains (ankyrin and tetratricopeptide repeats) are not annotated as thougherly or completely by bakta as desired. Consequently we are implementing an alternative approach for the annotation of these features.

## Structure of Domain Annotation Code
1. Parsing gbff files from bakta to faa files containing locus_tag and sequence information
2. hmmersearch to identify the presence of these domains
3. Adding this information back into the gbff files
I am planning to create a main script which links these different scripts sequentially (essentially making a package)

Initially I used --domtblout, but by changing this to --tblout additional parsing to extract counts is no longer required.

## To Do
1. Decide on hmmersearch thresholds (coverage and identity type things?)
2. Work out why the hmmer database is currently insufficient for tpr annotation, find alternatives?
3. Find a way (and location) to add the data back into the gbff file
