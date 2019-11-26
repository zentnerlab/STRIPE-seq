library(ComplexHeatmap)
library(multcomp)

setwd("/Users/gzentner/Desktop/tsrexplorer/yeast/STRIPE-seq/Gabe_yeast_work/TSR_complexity/")

# Combine SLIC-CAGE TSRs
# min.gapwidth = 40 merges TSRs less than 40 bp apart
slic_tsrs_reduced_40 <- GenomicRanges::reduce(c(exp@experiment$TSRs$SLIC_CAGE_100ng_1,exp@experiment$TSRs$SLIC_CAGE_100ng_2),
                                           drop.empty.ranges = FALSE, min.gapwidth = 40, ignore.strand = FALSE, with.revmap = FALSE)

# Filter out TSRs less than 5 bp in length
slic_tsrs_reduced_40_filtered_5 <- slic_tsrs_reduced_40[GenomicRanges::width(slic_tsrs_reduced_40) >= 5]

export.bed(slic_tsrs_reduced_40_filtered_5, con = "slic_tsrs_reduced_40_filtered_5.bed")

samples_for_complexity <- list(exp@experiment$TSSs$S288C_50ng_1,exp@experiment$TSSs$S288C_50ng_2,exp@experiment$TSSs$S288C_50ng_3,
                               exp@experiment$TSSs$S288C_100ng_1,exp@experiment$TSSs$S288C_100ng_2,exp@experiment$TSSs$S288C_100ng_3,
                               exp@experiment$TSSs$S288C_250ng_1,exp@experiment$TSSs$S288C_250ng_2,exp@experiment$TSSs$S288C_250ng_3,
                               exp@experiment$TSSs$SLIC_CAGE_100ng_1,exp@experiment$TSSs$SLIC_CAGE_100ng_2,
                               exp@experiment$TSSs$nanoCAGE_500ng_1,exp@experiment$TSSs$nanoCAGE_500ng_2,
                               exp@experiment$TSSs$nanoCAGE_25ng_1,exp@experiment$TSSs$nanoCAGE_25ng_2) %>%
    set_names(c("S288C_50ng_1","S288C_50ng_2","S288C_50ng_3",
                "S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
                "S288C_250ng_1","S288C_250ng_2","S288C_250ng_3",
                "SLIC_CAGE_100ng_1","SLIC_CAGE_100ng_2",
                "nanoCAGE_500ng_1","nanoCAGE_500ng_2",
                "nanoCAGE_25ng_1","nanoCAGE_25ng_2"))

tsr_complexity <- map(samples_for_complexity, ~ countOverlaps(slic_tsrs_reduced_40_filtered_5, .x) %>%
                          as.data.frame %>%
                          dplyr::rename(., nTSSs = .))

nTSSs <- bind_rows(tsr_complexity, .id = "sample") %>%
    gather(key = "sample", value = "nTSSs")
    
nTSSs <- nTSSs %>%
    mutate(log2 = log2(nTSSs + 1)) %>%
    mutate(sample = factor(sample, levels=c("S288C_50ng_1","S288C_50ng_2","S288C_50ng_3",
                                            "S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
                                            "S288C_250ng_1","S288C_250ng_2","S288C_250ng_3",
                                            "SLIC_CAGE_100ng_1","SLIC_CAGE_100ng_2",
                                            "nanoCAGE_500ng_1","nanoCAGE_500ng_2",
                                            "nanoCAGE_25ng_1","nanoCAGE_25ng_2")))

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
    set_names(c("S288C_50ng_1","S288C_50ng_2","S288C_50ng_3",
                "S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
                "S288C_250ng_1","S288C_250ng_2","S288C_250ng_3",
                "SLIC_CAGE_100ng_1","SLIC_CAGE_100ng_2",
                "nanoCAGE_500ng_1","nanoCAGE_500ng_2",
                "nanoCAGE_25ng_1","nanoCAGE_25ng_2"))

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
res.aov <- aov(log2 ~ sample, data = nTSSs)

summary(glht(res.aov, linfct = mcp(sample = "Tukey")))