# filtering vcf file. Requirements:
### allow 1/0, 0/1 and 1/1
### only include snp

# Original
# bcftools_viewCommand=view --include 'FMT/GT="1/1" && QUAL>=100 && FMT/DP>=10 && (FMT/AO)/(FMT/DP)>=0.05' core.raw.vcf

# bcftools is a real pain! I needed to install gsl as gsl=2.5 otherwise it wouldnt work and conda doesn't install it correctly

bcftools view --include 'QUAL>=100 && FMT/DP>=10 && (FMT/AO)/(FMT/DP)>=0.05' core.raw.vcf > core.processing.vcf
bcftools annotate --remove ^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL core.processing.vcf > core.processing2.vcf
awk 'NR<=27 || /TYPE=snp/' core.processing2.vcf > core.filt.vcf


# Takes a bed file and filters the vcf file to remove positions in the bed file
echo -e "LS398551.1\t5000000" | bedtools complement -i gilliam_rage_derived.bed -g <(echo -e "LS398551.1\t5000000") > inverse_regions.bed && bgzip -c core.filt.vcf > core.filt.vcf.gz && tabix -p vcf core.filt.vcf.gz && bcftools view -R inverse_regions.bed core.filt.vcf.gz -o core_regions.filt.vcf


# Reformating vcf file
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/AB\t%INFO/AO\t%INFO/DP\t%INFO/QA\t%INFO/QR\t%INFO/RO\t%INFO/TYPE\n' core_regions.filt.vcf > snp_table.txt
