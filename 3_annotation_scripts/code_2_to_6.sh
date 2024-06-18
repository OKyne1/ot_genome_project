# post bakta code

# 2. COG: repeated per gbff (RAGE env)
# (rage) yes2coding@MIC-NB39:/mnt/c/Users/oakem/Documents/230904_moru/240119_sequencing/240530_ot_genomes/clinical_strain_analysis/gbff_file_processing/processing$ 
python ../../../../../../ot_genome_project/2_annotation_scripts/2_cog_removal/scripts/1_cog_removal.py wg-ot-23-0*.gbff


# 3. domain annotation (hmmer env)
# (hmmer) yes2coding@MIC-NB39:/mnt/c/Users/oakem/Documents/230904_moru/240119_sequencing/240530_ot_genomes/clinical_strain_analysis/gbff_file_processing/processing$ 
bash ../../../../../../ot_genome_project/2_annotation_scripts/3_domain_annotation/domain_annotation_package/scripts/domain_annotation.sh modif*.gbff
# probably want to rm -r hmmer_outputs faa_files


# 4. ankyrin annotation (rage env)
# (rage) yes2coding@MIC-NB39:/mnt/c/Users/oakem/Orientia Biology Dropbox/Oakem Kyne/PC/Documents/230904_moru/240119_sequencing/240530_ot_genomes/clinical_strain_analysis/gbff_file_processing/processing$ 
bash ../../../../../../ot_genome_project/2_annotation_scripts/4_ank_annotation/scripts/gene_completeness_master_script.sh modified_wg-ot-23-0*.gbff
# This will result in an error if there are extra gbff files in the directory which arent specified - but the code probably still will have worked. Output should be the name of the gbff file with _anked.gbff appended

# 5. completeness checks (rage env)
# (rage) yes2coding@MIC-NB39:/mnt/c/Users/oakem/Orientia Biology Dropbox/Oakem Kyne/PC/Documents/230904_moru/240119_sequencing/240530_ot_genomes/clinical_strain_analysis/gbff_file_processing/processing$ 
bash ../../../../../../ot_genome_project/2_annotation_scripts/5_completeness_checks/scripts/gene_completeness_master_script.sh modified_wg-ot-23-0*anked.gbff
# need to clean up some of the files here, remove *.txt *.faa

# 6. rage classification (rage env) - issues are occuring...
