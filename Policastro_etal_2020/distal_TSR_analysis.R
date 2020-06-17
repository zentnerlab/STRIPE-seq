library(GenomicFeatures)
library(tsrexplorer)
library(tidyverse)
library(rtracklayer)
library(ChIPseeker)
library(org.Hs.eg.db)

# This script generates a set of consensus TSRs (in this case, TSRs detected in all three replicates of K562 STRIPE-seq)
# and annotates them relative to known TSSs using annotations derived from TxDb object. In this case, we consider TSRs < 1 kb
# from an annotated TSS to be 'proximal' and TSRs >= 1 kb from an annotated TSS to be 'distal.' After binning TSRs in this way,
# proximal and distal TSRs are overlapped with a set of K562 enhancers from EnhancerAtlas 2.0 (hg19 coordinates that were 
# lifted over to hg38). The numbers of proximal and distal TSRs overlapping and not overlapping enhancer annotations are then
# used to perform a chi-square test to see if there is a significant difference in the degree of overlap of proximal and distal
# TSRs with enhancers.

distal_tsr_dir <- file.path("human_work", "distal_TSRs")

if (!dir.exists(distal_tsr_dir)) {
    message("Creating directory 'human_work/distal_TSRs'...")
    dir.create(distal_tsr_dir)
} else {
    message("Directory 'human_work/distal_TSRs' already exists.")
}

# STRIPE-seq TSSs
STRIPE_TSSs <- map(list.files("human_data/STRIPE_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                       makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                                start.field = "TSS", end.field = "TSS")) %>%
    set_names(c("K562_100ng_1","K562_100ng_2","K562_100ng_3"))

# STRIPE-seq TSRs
STRIPE_TSRs <- map(list.files("human_data/STRIPE_TSRs/", full.names = TRUE), ~ read.delim(.x) %>%
                       makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                                start.field = "start", end.field = "end")) %>%
    set_names(c("K562_100ng_1","K562_100ng_2","K562_100ng_3"))

exp <- tsr_explorer(STRIPE_TSSs,STRIPE_TSRs)

# Generate TxDb from GTF
txdb <- makeTxDbFromGFF(file.path("human_data", "Homo_sapiens.GRCh38.98.chr.gtf"))

# Merge K562 replicate TSRs
merged_k562_tsrs <- GenomicRanges::reduce(c(exp@experiment$TSRs$K562_100ng_1,
                                            exp@experiment$TSRs$K562_100ng_2,
                                            exp@experiment$TSRs$K562_100ng_3))

rtracklayer::export(merged_k562_tsrs, con = file.path(distal_tsr_dir, "merged_TSRs.bed"), format = "bed")

# Determine overlap of each replicate with the merged TSR set
merged_k562_tsrs$rep1_overlap <- merged_k562_tsrs %over% exp@experiment$TSRs$K562_100ng_1
merged_k562_tsrs$rep2_overlap <- merged_k562_tsrs %over% exp@experiment$TSRs$K562_100ng_2
merged_k562_tsrs$rep3_overlap <- merged_k562_tsrs %over% exp@experiment$TSRs$K562_100ng_3

# Get TSRs shared across all three replicates
consensus_tsrs_all_reps <- merged_k562_tsrs[values(merged_k562_tsrs)[, "rep1_overlap"] == TRUE &
                                            values(merged_k562_tsrs)[, "rep2_overlap"] == TRUE &
                                            values(merged_k562_tsrs)[, "rep3_overlap"] == TRUE]

rtracklayer::export(consensus_tsrs_all_reps, con = file.path(distal_tsr_dir, "consensus_tsrs_all_reps.bed"), format = "bed")

# Get lists of proximal (< 1 kb from an annotated TSS) and distal (>= 1 kb from an annotated TSS)
tsr_anno_all_reps <- annotatePeak(consensus_tsrs_all_reps, TxDb = txdb, annoDb = "org.Hs.eg.db")

proximal_tsrs_all_reps <- tsr_anno_all_reps@anno[(values(tsr_anno_all_reps@anno)[, "distanceToTSS"] < 1000 & 
                                                      (values(tsr_anno_all_reps@anno)[, "distanceToTSS"] > -1000))]

rtracklayer::export(proximal_tsrs_all_reps, con = file.path(distal_tsr_dir, "proximal_tsrs_all_reps.bed"), format = "bed")

distal_tsrs_all_reps <- tsr_anno_all_reps@anno[(values(tsr_anno_all_reps@anno)[, "distanceToTSS"] <= -1000 | 
                                                    (values(tsr_anno_all_reps@anno)[, "distanceToTSS"] >= 1000))]
                                               
rtracklayer::export(distal_tsrs_all_reps, con = file.path(distal_tsr_dir, "distal_tsrs_all_reps.bed"), format = "bed")

# Overlap distal TSRs with K562 enhancer atlas regions

# Make an enhancer GRanges object and change seqlevelsStyle
enhancer_gr <- read.delim(file.path("human_data", "K562_enhancers_hg38.bed"), header=F) %>%
    makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3")

seqlevelsStyle(enhancer_gr) <- "ensembl"

# Count overlaps of proximal and distal TSRs with enhancer annotations 
proximal_enhancers <- proximal_tsrs_all_reps %over% enhancer_gr %>%
    sum

proximal_non_enhancers <- (length(proximal_tsrs_all_reps) - proximal_enhancers)

distal_enhancers <- countOverlaps(distal_tsrs_all_reps, enhancer_gr) %>%
    sum

distal_non_enhancers <- (length(distal_tsrs_all_reps) - distal_enhancers)

# Generate table for chi-square 
observed_table <- matrix(c(distal_enhancers, distal_non_enhancers, proximal_enhancers, proximal_non_enhancers), nrow = 2, ncol = 2, byrow = T)
rownames(observed_table) <- c("distal", "proximal")
colnames(observed_table) <- c("enhancer", "non-enhancer")

# Perform chi-square test
X <- chisq.test(observed_table)
pchisq(679.3, df = 1, lower.tail = F)