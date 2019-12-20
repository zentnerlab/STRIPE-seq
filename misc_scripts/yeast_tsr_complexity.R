library(GenomicFeatures)
library(AnnotationDbi)
library(ComplexHeatmap)
library(multcomp)

complexity_dir <- file.path("yeast_work", "TSR_complexity")

if (!dir.exists(complexity_dir)) {
    message("Creating directory 'yeast_work/TSR_complexity'...")
    dir.create(complexity_dir)
} else {
    message("Directory 'yeast_work/TSR_complexity' already exists.")
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

tsr_complexity <- map(samples_for_complexity, ~ countOverlaps(mrna_promoters, .x) %>%
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

# Fetch median TSR complexities for each sample
nTSSs %>% group_by(sample) %>% summarize(medians = median(log2))

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

#ANOVA
res.aov <- aov(log2 ~ sample, data = nTSSs)

tukey <- summary(glht(res.aov, linfct = mcp(sample = "Tukey")))

tukey_tibble <- tibble::enframe(tukey$test$coefficients)
tukey_tibble <- tukey_tibble %>%
    add_column(tukey$test$pvalues)
