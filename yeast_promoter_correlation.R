library(tsrexplorer)
library(tidyverse)
library(GenomicFeatures)
library(AnnotationDbi)
library(rtracklayer)
library(ComplexHeatmap)
library(viridis)
library(ComplexHeatmap)

# This script generates a list of yeast mRNA TSSs by filtering a TxDb object, then gets their 
# promoter coordinates and sums TSSs in each promoter for each sample. It then normalizes 
# promoter counts and performs Spearman correlation analysis to assess concordance of promoter 
# signals between methods.

prom_corr_dir <- file.path("yeast_work", "promoter_correlation")

if (!dir.exists(prom_corr_dir)) {
    message("Creating directory 'yeast_work/promoter_correlation'...")
    dir.create(prom_corr_dir)
} else {
    message("Directory 'yeast_work/promoter_correlation' already exists.")
}

# Generate promoter annotation file
annotation <- system.file("extdata", "yeast_annotation.gtf", package="tsrexplorer")
txdb <- makeTxDbFromGFF(annotation)

mrna_ids <- AnnotationDbi::select(txdb, keys = keys(txdb, "GENEID"), keytype = "GENEID", 
                              columns = columns(txdb)) %>% 
    as_tibble %>% 
    filter(str_detect(TXNAME, "mRNA"), TXCHROM != "Mito") %>%
    pull(TXNAME) 

### Retrieve promoter windows
promoters <- promoters(txdb, upstream = 250, downstream = 100)

mrna_promoters <- subset(promoters, tx_name %in% mrna_ids)

mrna_promoters$tx_name <- str_replace(mrna_promoters$tx_name, "-", "_")

export(mrna_promoters, file.path(prom_corr_dir, "mRNA_promoters.bed"), format = "bed")

# STRIPE-seq YPD TSSs
YPD_TSSs <- map(list.files("yeast_data/YPD_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                        makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                                 start.field = "TSS", end.field = "TSS")) %>%
        set_names(c("S288C_50ng_1","S288C_50ng_2","S288C_50ng_3",
                    "S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
                    "S288C_250ng_1","S288C_250ng_2","S288C_250ng_3"))
    
# SLIC-CAGE and nanoCAGE TSSs
CAGE_TSSs <- map(list.files("yeast_data/CAGE_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                         makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                                  start.field = "TSS", end.field = "TSS")) %>%
        set_names(c("SLIC_CAGE_100ng_1","SLIC_CAGE_100ng_2",
                    "nanoCAGE_500ng_1","nanoCAGE_500ng_2",
                    "nanoCAGE_25ng_1","nanoCAGE_25ng_2"))

ypd_cage_tss <- c(YPD_TSSs,CAGE_TSSs)

promoter_ntags <- map(ypd_cage_tss, function(x) {
    pairs <- findOverlapPairs(x, mrna_promoters, type = "any", ignore.strand = FALSE)
    tss_tibble <- as_tibble(pairs@first)
    promoter_tibble <- as_tibble(pairs@second)
    join <- bind_cols(tss_tibble, promoter_tibble)
    result <- group_by(join, seqnames, start1, end1, strand1, tx_name) %>%
        summarize(nTAGs = sum(nTAGs)) %>%
        makeGRangesFromDataFrame(start.field = "start1", end.field = "end1",
                               strand.field = "strand1", keep.extra.columns = TRUE)
    return(result)
})

# Store sample names for correlation
all <- c("S288C_50ng_1","S288C_50ng_2","S288C_50ng_3",
         "S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
         "S288C_250ng_1","S288C_250ng_2","S288C_250ng_3",
         "SLIC_CAGE_100ng_1","SLIC_CAGE_100ng_2",
         "nanoCAGE_500ng_1","nanoCAGE_500ng_2",
         "nanoCAGE_25ng_1","nanoCAGE_25ng_2")

# Read data into TSRexploreR
exp <- tsr_explorer(ypd_cage_tss, promoter_ntags)

# Normalize TSS and TSR counts
exp <- count_normalization(exp, data_type = "tss", threshold = 3, n_samples = 1, samples = all)

exp <- count_normalization(exp, data_type = "tsr", threshold = 3, n_samples = 1, samples = all)

# Generate a hierarchically clustered TSR heatmap with correlation values displayed
corr_matrix <- find_correlation(exp, data_type = "tsr", correlation_metric = "spearman")

cairo_pdf(file = file.path(prom_corr_dir, "promoter_correlation_hierarchical.pdf"), width = 22, height = 22)
Heatmap(corr_matrix, col = viridis(256), heatmap_legend_param = list(title = "Spearman"), 
        layer_fun = function(j, i, x, y, width, height, fill)
        {
            grid.text(sprintf("%.3f", pindex(corr_matrix, i, j)), x, y, gp = gpar(fontsize = 28, col = "white"))
        }
)
dev.off()