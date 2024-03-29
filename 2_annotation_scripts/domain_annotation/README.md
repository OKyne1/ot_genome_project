# Domain Annotation
Some important domains (ankyrin and tetratricopeptide repeats) are not annotated as thoroughly or completely by bakta as desired. Consequently we are implementing an alternative approach for the annotation of these domains.

## Structure of Domain Annotation Code
1. Parsing gbff files from bakta to faa files containing locus_tag and sequence information
2. hmmersearch to identify the presence of tprs or anks domains
3. Parsing hmmer output
4. Adding this information back into the gbff files

These 4 scripts are linked by [domain_annotation.sh](https://github.com/OKyne1/ot_genome_project/blob/main/2_annotation_scripts/domain_annotation/domain_annotation_package/scripts/domain_annotation.sh) which requires all files from the domain_annotation_package directory (excluding test directory) in the relative positions they are currently found in, and correct packages installed (see [environment.yml](https://github.com/OKyne1/ot_genome_project/blob/main/2_annotation_scripts/domain_annotation/domain_annotation_package/environment.yml)). The command line usage is then: `bash domain_annotaiton.sh <gbff files>`

## Parsing gbff files to faa files
hmmersearch requires a .faa file input. To produce this, gbff files from bakta were converted to .faa files with the locus tag as the header.
The script [gen_2_faa.py](https://github.com/OKyne1/ot_genome_project/blob/main/2_annotation_scripts/domain_annotation/1_gbff_2_faa/gen_2_faa.py)

## hmmersearch
The output .faa files are used for `hmmersearch`.
### Thresholds
By default the -E theshold is 10, however this is way higher than we would like.
The threshold used by bakta is `-E 1E-10`, we used this. This means that per hit there will be on average 10^-10 false positivies. Later on we may want to consider reducing the stringency of this threshold.

### Output Format
Initially I used `--domtblout`, but by changing this to `--tblout`. This produces a simplified table which means additional parsing for things like anks is not required.

### HMM file origins and choices
#### Anks
Choice of pfam .hmm file for anks was simple as there was only a single model (PF00023). 

My initial analysis showed that in both boryong and karp, ank proteins were annotated as ankyrin or unnamed and the product had some information (though not always the repeat number). 

This suggests that it will be worth overwritting the product with the number of ankyrin repeats.

#### TPRs
TPR domain .hmm files were more challenging due to an abundance of different files. This is the method used to identify the .hmm files to use for tpr identification:
1. Interpro search for "Tetratricopeptide repeat"
2. Export the information as a .tsv file
3. Filter out all non-PFAM entries
4. Remove all enties with names other than exactly "Tetratricopeptide repeat"
5. Identify the short name for each tpr type.
6. Download and concatenate the files

This resulted in 19 different tpr .hmm models which could be concatenated into a single file.

##### TPR .hmm files used:
| TPR   | PFAM ID   |
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

Analysis of the product names (from bakta) for the tprs annotated through hmmer show:
1. Most proteins are correctly annotated with tpr product details, and with a gene name of tpr or no name
2. There are some proteins which are called pilW but shouldnt be
3. There are a couple of cases where there are other details/product descriptions (not pilW or tpr)

## Parsing the Hmmer Output
Different approaches are required for both anks and tprs. For anks we give details on both the repeat and the number of repeats, where as for tprs, this isn't possible as there is no standardised database.

For parsing, the code extracts the locus tag from the hmmer file. Then if -dom ank is specified it extracts the number of ank repeats. But if -dom tpr is specified it outputs "Tetratricopeptide repeat protein" and doesn't specify the repeat number (and excluded repeated values)

## Overwriting Product
### Location
**Gene name**: I don't think a protein should be named after a domain. So, I plan to remove any of these if present.

**Product**: This is where I plant to store the TPR and ANK information. From what I've pulled out of the gbff files, it looks like I can just overwrite this information for (nearly) all cases (except bamD and traG).

### Script
The script [4_overwriting_gbff.py](https://github.com/OKyne1/ot_genome_project/blob/main/2_annotation_scripts/domain_annotation/domain_annotation_package/scripts/4_overwriting_gbff.py) can takes a gbff file and txt file. This overwrites the product in the gbff file based on the data in the txt file (col1=locus_tag and col2=product) and removes any gene names present.

Investigation of the products for existing entries (8 genomes) showed that overwriting anks would cause no errors. However, overwriting tprs would cause loss of 9 legitimate **traG** proteins and also 8 legitimate **bamD**. Consequently, when overwriting tpr information, entries with gene names of traG or bamD are excluded. Those with product, but not gene name are just overwritten and they were not legitimate cases (of traG and there were none for bamD).

Other things overwritten in tprs; "**SycD/LcrH family type III secretion system chaperone**" and **pilW** proteins. There seems to be extremely weak evidence for these entries and thus it is reasonable to overwrite.

## Current code limitations
Only processes tprs or anks, if additional domains are desired the hmmer file will need to be added to the hmm_files directory and script 3 and 4 will need to be modified to take/process additional domains. This needs to consider what to overwrite and whether to include the number of repeats.

We have also not investigated things annotated as tprs or anks that aren't detected by bakta. This may be worth investigating in the furture, but is unlikely to cause any major problems.
