# Plot generation method
# If using IR regions only:
# 1. Use the rage derived bed file to mask the genome (bedtools maskfasta -fi input.fasta -bed rage_derived_regions.bed -fo output_masked.fasta)
# 2. Remove masked regions from the fasta files (sed -i 's/N//g' filename)
# 3. Combine fasta files into a single file
# 4. Rename genomes, these are the ones which will be on the plot
#
# If using whole genomes:
# 1. Combine genomes (use cat)
# 2. Rename genomes, these are the ones which will be on the plot
#
# Multiple sequence alignment (both):
# 1. multiple sequence alignment (minimap2 -x asm20 -c --eqx -D -P -dual=no 31.fasta_masked.fasta GCA_900327255.1_UT76_genomic_masked.fasta GCA_900327275.1_Karp_genomic_masked.fna > output.paf)
# 2. Use the paf file in the script below
#
# Generation of the files for rage annotation
# - bed files contain the information
# - If plotting IR only these will need preprocessing
#
# IR bed file preprocessing
# - Need to create an inverse of the bed files (using the genome length)
# - Then need to remove the spaces the RAGE regions to give the boundaries
# Done using 1_bed_inverse.py and 2_IR_only_bed_conversion.py to convert the bed files --> this isn't very well annotated, give chatgpt the code and ask how it works.
#
# For pairwise use whole genome use `minimap2 -x asm20 -c --eqx -secondary=no {target.fasta} {query.fasta} > {output.alignment}`
-------------------------------------------------------------------------------
# Plotting the paf alignments to visualise the synteny
library(SVbyEye)
library(GenomicRanges)

#-------------------------------------------------------------------------------
setwd("C:/Users/oakem/Documents/230904_moru/240726_svbyeye_all_genomes/paf_alignment")
# load paf file
paf.file <- "all_genomes_wuj_rotated_ir_only.paf"
annot_file <- "all_17_r14_IR_contiguous.bed"  # Replace with your actual BED file path

#-------------------------------------------------------------------------------

# Read in PAF
paf.table <- readPaf(
  paf.file = paf.file,
  include.paf.tags = TRUE, restrict.paf.tags = "cg"
)

# Read the BED file, assuming it has columns: chrom, start, end
annot.gr <- read.table(annot_file, header = FALSE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
colnames(annot.gr) <- c("seqnames", "start", "end")  # Assign column names
annot.gr$strand <- "+"  # I think Granges objects need the direction - this is just a place holder value


# Convert data.frames to GRanges objects
annot_gr_ranges <- GRanges(seqnames = annot.gr$seqnames,
                           ranges = IRanges(start = annot.gr$start,
                                            end = annot.gr$end),
                           strand = annot.gr$strand)

seqnames.order <- c("wgot004", "wgot005", "wgot019", "UT76", "TW1", "wuj2014", "UT176", "Karp", "Boryong", "Gilliam", "Ikeda", "TA763", "TA686", "TW22", "wgot003", "Kato", "wgot013")
# Assuming plt is created with plotAVA() as before
plt <- plotAVA(paf.table = paf.table, color.by = "direction", seqnames.order = seqnames.order)


# Add annotations to the plot
plt_with_annotation <- addAnnotation(ggplot.obj = plt, annot.gr = annot_gr_ranges,
                                     coordinate.space = 'self', y.label.id = 'seqnames',
                                     annotation.level = 0, shape = "rectangle")

plt_with_annotation


?plotAVA

#-------------------------------------------------------------------------------
# colouring by identity
plt2 <- plotAVA(paf.table = paf.table, color.by = "identity", seqnames.order = seqnames.order, perc.identity.breaks = c(90, 95, 99, 99.5, 99.6, 99.7, 99.8, 99.9))


# Add annotations to the plot
plt_with_annotation2 <- addAnnotation(ggplot.obj = plt2, annot.gr = annot_gr_ranges,
                                     coordinate.space = 'self', y.label.id = 'seqnames',
                                     annotation.level = 0, shape = "rectangle")
plt_with_annotation2
#-------------------------------------------------------------------------------
# adding in deletions/insertions
# colouring by identity
plt3 <- plotAVA(paf.table = paf.table, color.by = "identity", seqnames.order = seqnames.order, perc.identity.breaks = c(90, 95, 99, 99.5, 99.6, 99.7, 99.8, 99.9), min.deletion.size = 5000, 
                min.insertion.size = 5000, highlight.sv = 'outline')


# Add annotations to the plot
plt_with_annotation3 <- addAnnotation(ggplot.obj = plt3, annot.gr = annot_gr_ranges,
                                      coordinate.space = 'self', y.label.id = 'seqnames',
                                      annotation.level = 0, shape = "rectangle")
plt_with_annotation3
