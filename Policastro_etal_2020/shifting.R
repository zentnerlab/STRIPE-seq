library(tsrexplorer)
library(GenomicRanges)
library(tidyverse)
library(CAGEr)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Note that the newest version of TSRexploreR is required for the syntax in this first part to work
# devtools::install_github("rpolicastro/tsrexplorer", ref = "clean", force = TRUE)

# Read YPD and diamide 100 ng STRIPE-seq TSSs into variables

# STRIPE-seq YPD TSSs
YPD_TSSs <- map(list.files("yeast_data/YPD_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                  makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                           start.field = "TSS", end.field = "TSS")) %>%
  set_names(c("S288C_50ng_1","S288C_50ng_2","S288C_50ng_3",
              "S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
              "S288C_250ng_1","S288C_250ng_2","S288C_250ng_3"))

# STRIPE-seq diamide TSSs
diamide_TSSs <- map(list.files("yeast_data/diamide_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                      makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                               start.field = "TSS", end.field = "TSS")) %>%
  set_names(c("S288C_diamide_100ng_1","S288C_diamide_100ng_2","S288C_diamide_100ng_3"))

# Generate list of all TSS objects
full_TSS_set <- c(YPD_TSSs, diamide_TSSs)

# Add a "score" column to TSS GRanges (duplication and column renaming of TSRchitect "nTAGs" column)
full_TSS_set <- lapply(full_TSS_set, function(x) {x$score <- x$nTAGs; return(x)})

# Read TSSs and TSRs into TSRexploreR and prepare them for analysis
yeast_exp <- tsr_explorer(full_TSS_set)

yeast_exp <- format_counts(yeast_exp, data_type = "tss")

# Write STRIPE-seq TSSs to CTSS format
iwalk(yeast_exp@experiment$TSSs, function(x,y) {x %>% as.data.frame %>% select(-width, -end, -nTAGs, -isreal) %>% 
    mutate(seqnames = str_c("chr", seqnames)) %>%
    write.table(str_c(y, ".txt"), col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)})

# Group YPD 100 ng and YPD + diamide STRIPE-seq samples
diamide <- c("S288C_100ng_1.txt","S288C_100ng_2.txt","S288C_100ng_3.txt",
             "S288C_diamide_100ng_1.txt","S288C_diamide_100ng_2.txt","S288C_diamide_100ng_3.txt")

# Copy 100 ng control and diamide CTSS files to a new directory
if (!dir.exists("cager_tss")){
  message("Creating directory 'cager_tss'...")
  dir.create("cager_tss")
} else {
  message("Directory 'cager_tss' already exists...")
}

file.copy(diamide, "cager_tss/", overwrite = TRUE)

# Read CTSS files into a CAGEr CAGEset object
TSSs_for_CAGEr <- list.files("cager_tss", full.names = TRUE)

samples <- make.names(TSSs_for_CAGEr)

myCAGEset <- new("CAGEset", genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3",
                 inputFiles = TSSs_for_CAGEr, inputFilesType = "ctss",
                 sampleLabels = c("ctrl_1", "ctrl_2", "ctrl_3",
                                  "diamide_1", "diamide_2", "diamide_3"))
                 
getCTSS(myCAGEset)

# Merge control and diamide replicate samples
mergeSamples(myCAGEset, mergeIndex = c(1,1,1,2,2,2),
             mergedSampleLabels = c("control", "diamide"))

# Determine alpha for power-law normalization
plotReverseCumulatives(myCAGEset, fitInRange = c(5, 1000), onePlot = TRUE)
                                      
# Normalize samples
normalizeTagCount(myCAGEset, method = "powerLaw",
                  fitInRange = c(5, 1000), alpha = 1.93, T = 1e6)                 

# Generate CTSS clusters
clusterCTSS(object = myCAGEset, threshold = 1, thresholdIsTpm = TRUE,
            nrPassThreshold = 1, method = "distclu", maxDist = 40,
            removeSingletons = TRUE, keepSingletonsAbove = 3)

# Get signal within clusters
cumulativeCTSSdistribution(myCAGEset, clusters = "tagClusters")
quantilePositions(myCAGEset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)

# Combine tag clusters
aggregateTagClusters(myCAGEset, tpmThreshold = 3,
                     qLow = 0.1, qUp = 0.9, maxDist = 100)

consensusCl <- consensusClusters(myCAGEset)

cumulativeCTSSdistribution(myCAGEset, clusters = "consensusClusters")

# Calculate shifting score
scoreShift(myCAGEset, groupX = "control", groupY = "diamide",
           testKS = TRUE, useTpmKS = FALSE)

# Filter shifted tag clusters and write to a table
shifting.promoters <- getShiftingPromoters(myCAGEset,
                                           tpmThreshold = 5,
                                           fdrThreshold = 0.01,
                                           scoreThreshold = 0.6)

write.table(shifting.promoters, "shifting_tsrs.txt", sep = "\t", row.names = F, col.names = T, quote = F)
