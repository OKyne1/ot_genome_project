# Domain Annotation
Some important domains (ankyrin and tetratricopeptide repeats) are not annotated as thoroughly or completely by bakta as desired. Consequently we are implementing an alternative approach for the annotation of these features.

## Structure of Domain Annotation Code
1. Parsing gbff files from bakta to faa files containing locus_tag and sequence information
2. hmmersearch to identify the presence of these domains
3. Adding this information back into the gbff files
I am planning to create a main script which links these different scripts sequentially (essentially making a package)

## Parsing gbff files to faa files
hmmersearch requires a .faa file input. To produce this, gbff files from bakta were converted to .faa files with the locus tag as the header.
The script [gen_2_faa.py](https://github.com/OKyne1/ot_genome_project/blob/main/2_annotation_scripts/domain_annotation/1_gbff_2_faa/gen_2_faa.py)

## hmmersearch parameters
### Thresholds
By default the -E theshold is 10, however this is way higher than we would like.
The threshold used by bakta is -E 1E-10, we used this. This means that per hit there will be on average 10^-10 false positivies. Later on we may want to consider reducing the stringency of this threshold.

### Output Format
Initially I used **--domtblout**, but by changing this to **--tblout**. This produces a simplified table which means additional parsing for things like anks is not required.

## HMM file origins and choices
### Anks
Choice of pfam .hmm file for anks was simple as there was only a single model (PF00023). 

My initial analysis showed that in both boryong and karp, ank proteins were annotated as ankyrin or unnamed and the product had some information (though not always the repeat number). 

This suggests that it will be worth overwritting the product with the number of ankyrin repeats.

### TPRs
TPR domain .hmm files were more challenging due to an abundance of different files. This is the method used to identify the .hmm files to use for tpr identification:
1. Interpro search for "Tetratricopeptide repeat"
2. Export the information as a .tsv file
3. Filter out all non-PFAM entries
4. Remove all enties with names other than exactly "Tetratricopeptide repeat"
5. Identify the short name for each tpr type.
6. Download and concatenate the files

This resulted in 19 different tpr .hmm models which could be concatenated into a single file.

#### TPR .hmm files used:
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

## Adding back to gbff
### Location
**Gene name**: I don't think a protein should be named after a domain. So, I plan to remove any of these present.

**Product**: This is where I plant to store the TPR and ANK information. From what I've pulled out of the gbff files, it looks like I can just overwrite this information for (nearly) all cases.

### Parsing the Hmmer Output
Different approaches are required for both anks and tprs. For anks we give details on both the repeat and the number of repeats, where as for tprs, this isn't possible as there is no standardised database.

#### Parsing Ank Hmmer Output
Need to convert to a txt file with this format:

locus_tag    product

Product should be "Ankyrin repeat protein with "$number" ankyrin repeats"

#### Parsing tpr Hmmer Output

### Overwriting Product
The script [overwriting_gbff.py](https://github.com/OKyne1/ot_genome_project/blob/main/2_annotation_scripts/domain_annotation/3_adding_2_gbff/overwriting_gbff.py) can takes a gbff file and txt file. This overwrites the product in the gbff file based on the data in the txt file (col1=locus_tag and col2=product).

### Other annotations
Currently this work hasn't tried to detect annotations which aren't detected by hmmer. This may be a limitation and could be worth investigating.

## To Do
Currently the tprs are producing multiple tpr hits of different types on the same locus tag. Are these unique or is it that these are duplicated hits? Yes these hits are being duplicated. For a location it can produce multiple different hits for the different models.
This means that the number of tprs is hard to classify using this data, however the it is still possible to identify the presence of a tpr. What is important in classification of these proteins?
