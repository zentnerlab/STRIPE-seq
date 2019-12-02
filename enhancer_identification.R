library(ChIPseeker)
library(org.Hs.eg.db)

setwd(file.path(baseDir, "human_work"))

# Generate TxDb from Ensembl GTF
txdb <- makeTxDbFromGFF(file.path(baseDir, "Homo_sapiens.GRCh38.98.gtf"))

# Merge K562 replicate TSRs
merged_k562_tsrs <- GenomicRanges::reduce(c(exp@experiment$TSRs$K562_100ng_1,exp@experiment$TSRs$K562_100ng_2,exp@experiment$TSRs$K562_100ng_3))

# Determine overlap of each replicate with the merged TSR set
merged_k562_tsrs$rep1_overlap <- merged_k562_tsrs %over% exp@experiment$TSRs$K562_100ng_1
merged_k562_tsrs$rep2_overlap <- merged_k562_tsrs %over% exp@experiment$TSRs$K562_100ng_2
merged_k562_tsrs$rep3_overlap <- merged_k562_tsrs %over% exp@experiment$TSRs$K562_100ng_3

# Get TSRs shared across at least two replicates
consensus_TSRs <- merged_k562_tsrs[values(merged_k562_tsrs)[, "rep1_overlap"] == TRUE & 
                                   values(merged_k562_tsrs)[, "rep2_overlap"] == TRUE |
                                   values(merged_k562_tsrs)[, "rep1_overlap"] == TRUE & 
                                   values(merged_k562_tsrs)[, "rep3_overlap"] == TRUE |    
                                   values(merged_k562_tsrs)[, "rep2_overlap"] == TRUE &
                                   values(merged_k562_tsrs)[, "rep3_overlap"] == TRUE]

# Get a list of consensus TSRs greater than 3 kb from an annotated TSS
tsr_anno <- annotatePeak(consensus_TSRs, TxDb = txdb, annoDb = "org.Hs.eg.db")

distal_tsrs <- tsr_anno@anno[values(tsr_anno@anno)[, "distanceToTSS"] <= -3000 | values(tsr_anno@anno)[, "distanceToTSS"] >= 3000]

export.bed(distal_tsrs, con = "distal_TSRs.bed")
