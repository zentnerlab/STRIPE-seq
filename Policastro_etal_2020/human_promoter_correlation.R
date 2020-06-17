library(tsrexplorer)
library(tidyverse)
library(GenomicFeatures)
library(AnnotationDbi)
library(rtracklayer)
library(ComplexHeatmap)
library(viridis)
library(ComplexHeatmap)

# This script generates a list of human mRNA TSSs by filtering a TxDb object, then gets their 
# promoter coordinates and sums TSSs in each promoter for each sample. It then normalizes 
# promoter counts and performs Spearman correlation analysis to assess concordance of promoter 
# signals between methods.

prom_corr_dir <- file.path("human_work", "promoter_correlation")

if (!dir.exists(prom_corr_dir)) {
    message("Creating directory 'human_work/promoter_correlation'...")
    dir.create(prom_corr_dir)
} else {
    message("Directory 'human_work/promoter_correlation' already exists.")
}

# Generate promoter annotation file
txdb <- makeTxDbFromGFF(file.path("human_data", "Homo_sapiens.GRCh38.98.chr.gtf"))
tx_list <- read.delim("human_data/protein_coding_ens_tx.txt", stringsAsFactors = FALSE) %>%
    pull

mrna_ids <- unique(AnnotationDbi::select(txdb, keys = keys(txdb), keytype = "GENEID", 
                                         columns = columns(txdb))) %>% 
    as_tibble %>% 
    filter(TXNAME %in% tx_list, TXCHROM != "MT", TXCHROM != "Y") %>%
    pull(TXNAME) 

### Retrieve promoter windows
promoters <- promoters(txdb, upstream = 500, downstream = 500)

mrna_promoters <- subset(promoters, tx_name %in% mrna_ids)

export(mrna_promoters, file.path(prom_corr_dir, "mRNA_promoters.bed"), format = "bed")

# STRIPE-seq TSSs
STRIPE_TSSs <- map(list.files("human_data/STRIPE_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                    makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                             start.field = "TSS", end.field = "TSS")) %>%
    set_names(c("K562_100ng_1","K562_100ng_2","K562_100ng_3"))

# SLIC-CAGE and nanoCAGE TSSs
CAGE_TSSs <- map(list.files("human_data/CAGE_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                     makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                              start.field = "TSS", end.field = "TSS")) %>%
    set_names(c("nanoCAGE_XL_7.5ug_1",
                "CAGE_10ug_1","CAGE_10ug_2",
                "RAMPAGE_5ug_1","RAMPAGE_5ug_2"))

stripe_cage_tss <- c(STRIPE_TSSs,CAGE_TSSs)

promoter_ntags <- map(stripe_cage_tss, function(x) {
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
all <- c("K562_100ng_1","K562_100ng_2","K562_100ng_3",
         "CAGE_10ug_1","CAGE_10ug_2",
         "RAMPAGE_5ug_1","RAMPAGE_5ug_2",
         "nanoCAGE_XL_7.5ug_1")

# Read data into TSRexploreR
exp <- tsr_explorer(stripe_cage_tss, promoter_ntags)

# Normalize TSS and TSR counts
exp <- count_normalization(exp, data_type = "tss", threshold = 3, n_samples = 1, samples = all)

exp <- count_normalization(exp, data_type = "tsr", threshold = 3, n_samples = 1, samples = all)

# Generate a hierarchically clustered TSR heatmap with correlation values displayed
corr_matrix <- find_correlation(exp, data_type = "tsr", correlation_metric = "spearman")

cairo_pdf(file = file.path(prom_corr_dir, "promoter_correlation_hierarchical.pdf"), width = 16, height = 16)
Heatmap(corr_matrix, col = viridis(256), heatmap_legend_param = list(title = "Spearman"), 
        layer_fun = function(j, i, x, y, width, height, fill)
        {
            grid.text(sprintf("%.3f", pindex(corr_matrix, i, j)), x, y, gp = gpar(fontsize = 28, col = "white"))
        }
)
dev.off()