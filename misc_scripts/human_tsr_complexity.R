library(ComplexHeatmap)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(multcomp)

# Set base directory for analyses
baseDir <- "/Users/gzentner/Desktop/tsrexplorer/STRIPE-seq/"
setwd(baseDir)

if (!dir.exists(file.path(baseDir, "human_work/TSR_complexity"))){
  print("Creating directory 'human_work/TSR_complexity' and changing working directory...")
  dir.create(file.path(baseDir, "human_work/TSR_complexity"))
  setwd(file.path(baseDir, "human_work/TSR_complexity/"))
} else {
  print("Directory 'human_work/TSR_complexity' already exists, changing working directory...")
  setwd(file.path(baseDir, "human_work/TSR_complexity/"))
}

# Combine CAGE TSRs
# min.gapwidth = 40 merges TSRs less than 40 bp apart
# cage_tsrs_reduced_40 <- GenomicRanges::reduce(c(exp@experiment$TSRs$CAGE_10ug_1,exp@experiment$TSRs$CAGE_10ug_2),
#                                            drop.empty.ranges = FALSE, min.gapwidth = 40, ignore.strand = FALSE, with.revmap = FALSE)
# 
# # Filter out TSRs less than 5 bp in length
# cage_tsrs_reduced_40_filtered_5 <- cage_tsrs_reduced_40[GenomicRanges::width(cage_tsrs_reduced_40) >= 5]
# 
# export.bed(cage_tsrs_reduced_40_filtered_5, con = "cage_tsrs_reduced_40_filtered_5.bed")

# Get promoter windows
promoters <- promoters(genes(TxDb.Hsapiens.UCSC.hg38.knownGene), upstream = 1000, downstream = 1000)
promoters <- trim(promoters)
chr_names <- seqlevels(promoters)
chr_names <- stringr::str_replace(chr_names, "chr", "")
seqlevels(promoters) <- chr_names

samples_for_complexity <- list(exp@experiment$TSSs$K562_100ng_1,
                               exp@experiment$TSSs$K562_100ng_2,
                               exp@experiment$TSSs$K562_100ng_3,
                               exp@experiment$TSSs$nanoCAGE_XL_7.5ug_1,
                               exp@experiment$TSSs$CAGE_10ug_1,
                               exp@experiment$TSSs$CAGE_10ug_2,
                               exp@experiment$TSSs$RAMPAGE_5ug_1,
                               exp@experiment$TSSs$RAMPAGE_5ug_2) %>%
    set_names(c("K562_100ng_1","K562_100ng_2","K562_100ng_3",
                "nanoCAGE_XL_7.5ug_1",
                "CAGE_10ug_1","CAGE_10ug_2",
                "RAMPAGE_5ug_1","RAMPAGE_5ug_2"))

tsr_complexity <- map(samples_for_complexity, ~ countOverlaps(promoters, .x) %>%
                          as.data.frame %>%
                          dplyr::rename(., nTSSs = .))

nTSSs <- bind_rows(tsr_complexity, .id = "sample") %>%
    gather(key = "sample", value = "nTSSs")
    
nTSSs <- nTSSs %>%
    mutate(log2 = log2(nTSSs + 1)) %>%
    mutate(sample = factor(sample, levels=c("K562_100ng_1","K562_100ng_2","K562_100ng_3",
                                            "nanoCAGE_XL_7.5ug_1",
                                            "CAGE_10ug_1","CAGE_10ug_2",
                                            "RAMPAGE_5ug_1","RAMPAGE_5ug_2")))

p <- ggplot(nTSSs, aes(x = sample, y = log2)) + 
    geom_jitter(color = "lightgrey", size = 0.5) +
    geom_boxplot(fill = NA, aes(color = sample), outlier.shape = NA) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1)) +
    scale_fill_viridis_d() +
    scale_color_viridis_d()

ggsave("tsr_complexity.png", plot = p, device = "png", type = "cairo", height = 12, width = 12)

# Convert complexities to a log2-transformed data frame and plot a correlation heatmap
tsr_complexity_df <- as.data.frame(tsr_complexity)
tsr_complexity_df_log2 <- log2(tsr_complexity_df + 1) %>%
    set_names(c("K562_100ng_1","K562_100ng_2","K562_100ng_3",
                "nanoCAGE_XL_7.5ug_1",
                "CAGE_10ug_1","CAGE_10ug_2",
                "RAMPAGE_5ug_1","RAMPAGE_5ug_2"))

corr_matrix <- cor(tsr_complexity_df_log2, method = "spearman")

pdf(file = "tsr_complexity_hierarchical.pdf", width = 16, height = 16)
Heatmap(corr_matrix, col = viridis(256), heatmap_legend_param = list(title = "Spearman"), 
        layer_fun = function(j, i, x, y, width, height, fill)
            {
                grid.text(sprintf("%.3f", pindex(corr_matrix, i, j)), x, y, gp = gpar(fontsize = 16, col = "white"))
            }
        )
dev.off()

# ANOVA
# res.aov <- aov(log2 ~ sample, data = nTSSs)
# 
# summary(glht(res.aov, linfct = mcp(sample = "Tukey")))