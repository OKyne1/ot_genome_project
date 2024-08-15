# Plot generation method for pairwise whole genome analysis
#
# 1. Rename fasta files (after ">" if desired)
# 2. Perform sequence alignment (minimap2 -x asm20 -c –eqx –secondary=no {target.fasta} {query.fasta} > {output.alignment})
# 3. Use the paf file in the code below
# 4. Load the bed files for complete and rage derived regions in the code below
# Note: this will not work if a complete RAGE bed file is empty and so the appropriate region should be commented/removed
#-------------------------------------------------------------------------------
library("SVbyEye")
library(GenomicRanges)
setwd("C:/Users/oakem/Documents/230904_moru/240726_svbyeye_all_genomes/paf_alignment_full_genome/wgot004_005_019_whole_genome_synteny/19vs5/")


#-------------------------------------------------------------------------------
# wuj2014-tw1 (post-wuj2014 rotation)
paf.file2 <- "19_5.paf"

# Read in PAF
paf.table2 <- readPaf(
  paf.file = paf.file2,
  include.paf.tags = TRUE, restrict.paf.tags = "cg"
)
plt <- plotMiro(paf.table = paf.table2, color.by = "identity", binsize = 1000, 
                perc.identity.breaks = c(99.5, 99.6, 99.7, 99.8, 99.9), 
                min.deletion.size = 1000, min.insertion.size = 1000,  highlight.sv = "outline")
plt
#-------------------------------------------------------------------------------
# Adding in RAGE boundaries

# Annotations for 4
annot_file <- "modified_wgot005_anked_completeness_checked_spotted_rage_derived.bed"  # Replace with your actual file path if different
annot.gr <- read.table(annot_file, header = FALSE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
colnames(annot.gr) <- c("seqnames", "start", "end")  # Assign column names
annot.gr$strand <- "+"  # I think Granges objects need the direction - this is just a place holder value
annot_gr_ranges <- GRanges(seqnames = annot.gr$seqnames,
                           ranges = IRanges(start = annot.gr$start,
                                            end = annot.gr$end),
                           strand = annot.gr$strand)

# Annotations for 19
annot_file2 <- "modified_wgot019_anked_completeness_checked_spotted_rage_derived.bed"  # Replace with your actual file path if different
annot.gr2 <- read.table(annot_file2, header = FALSE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
colnames(annot.gr2) <- c("seqnames", "start", "end")  # Assign column names
annot.gr2$strand <- "+"  # I think Granges objects need the direction - this is just a place holder value
annot_gr_ranges2 <- GRanges(seqnames = annot.gr2$seqnames,
                           ranges = IRanges(start = annot.gr2$start,
                                            end = annot.gr2$end),
                           strand = annot.gr2$strand)
######## No RAGEs so these bits were commented for this comparison. Change if there are complete RAGEs in the genomes being compared
# # Annotation for 5 complete RAGE
# annot_file3 <- "modified_wgot004_anked_completeness_checked_spotted_complete_RAGEs.bed"  # Replace with your actual file path if different
# annot.gr3 <- read.table(annot_file3, header = FALSE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
# colnames(annot.gr3) <- c("seqnames", "start", "end")  # Assign column names
# annot.gr3$strand <- "+"  # I think Granges objects need the direction - this is just a place holder value
# annot_gr_ranges3 <- GRanges(seqnames = annot.gr3$seqnames,
#                             ranges = IRanges(start = annot.gr3$start,
#                                              end = annot.gr3$end),
#                             strand = annot.gr3$strand)
# 
# # Annotation for 19 complete RAGE
# #annot_file4 <- "modified_wgot019_anked_completeness_checked_spotted_complete_RAGEs.bed"  # Replace with your actual file path if different
# #annot.gr4 <- read.table(annot_file4, header = FALSE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
# #colnames(annot.gr4) <- c("seqnames", "start", "end")  # Assign column names
# #annot.gr4$strand <- "+"  # I think Granges objects need the direction - this is just a place holder value
# #annot_gr_ranges4 <- GRanges(seqnames = annot.gr4$seqnames,
#  #                           ranges = IRanges(start = annot.gr4$start,
#   #                                           end = annot.gr4$end),
#    #                         strand = annot.gr4$strand)


plt1 <- addAnnotation(ggplot.obj = plt, 
                      annot.gr = annot_gr_ranges, coordinate.space = "target", 
                      annotation.label = 'RAGE Dervied', shape = "rectangle")
plt2 <- addAnnotation(ggplot.obj = plt1, 
                      annot.gr = annot_gr_ranges2, coordinate.space = "query", 
                      annotation.label = 'RAGE Dervied', shape = "rectangle")
# plt3 <- addAnnotation(ggplot.obj = plt2, 
#                       annot.gr = annot_gr_ranges3, coordinate.space = "target", 
#                       annotation.label = 'Complete RAGE', shape = "rectangle")
#plt4 <- addAnnotation(ggplot.obj = plt3, 
 #                     annot.gr = annot_gr_ranges3, coordinate.space = "query", 
  #                    annotation.label = 'Complete RAGE', shape = "rectangle")
plt2

