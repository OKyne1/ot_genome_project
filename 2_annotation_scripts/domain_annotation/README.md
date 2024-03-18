# Domain Annotation
Some important domains (ankyrin and tetratricopeptide repeats) are not annotated as thougherly or completely by bakta as desired. Consequently we are implementing an alternative approach for the annotation of these features.

## Structure of Domain Annotation Code
1. Parsing gbff files from bakta to faa files containing locus_tag and sequence information
2. hmmersearch to identify the presence of these domains
3. Adding this information back into the gbff files
I am planning to create a main script which links these different scripts sequentially (essentially making a package)

Initially I used --domtblout, but by changing this to --tblout additional parsing to extract counts is no longer required.

## hmmersearch parameters

## HMM file origins and choices
### Anks
Choice of pfam .hmm file for anks was simple as there was only a single model (PF00023). 

### TPRs
TPR domain .hmm files were more challenging due to an abundance of different files. This is the method used to identify the .hmm files to use for tpr identification:
1. Interpro search for "Tetratricopeptide repeat"
2. Export the information as a .tsv file
3. Filter out all non-PFAM entries
4. Remove all enties with names other than exactly "Tetratricopeptide repeat"
5. Identify the short name for each tpr type.
This resulted in 19 different tpr .hmm models which could be concatenated into a single file.

### TPR .hmm files used:
| TPR   | PF ID   |
|-------|---------|
|	TPR_1	|	PF00515	|
|	TPR_2	|	PF07719	|
|	TPR_3	|	PF07720	|
|	TPR_4	|	PF07721	|
|	TPR_5	|	PF12688	|
|	TPR_6	|	PF13174	|
|	TPR_7	|	PF13176	|
|	TPR_8	|	PF13181	|
|	TPR_9	|	PF13371	|
|	TPR_10	|	PF13374	|
|	TPR_12	|	PF13424	|
|	TPR_14	|	PF13428	|
|	TPR_15	|	PF13429	|
|	TPR_17	|	PF13431	|
|	TPR_16	|	PF13432	|
|	TPR_18	|	PF13512	|
|	TPR_19	|	PF14559	|
|	TPR_20	|	PF14561	|
|	TPR_22	|	PF18833	|




## To Do
1. Decide on hmmersearch thresholds (coverage and identity type things?)
2. Work out why the hmmer database is currently insufficient for tpr annotation, find alternatives?
3. Find a way (and location) to add the data back into the gbff file
