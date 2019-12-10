library(tsrexplorer)
library(tidyverse)
library(GenomicRanges)
library(viridis)
library(ComplexHeatmap)
library(rtracklayer)
library(org.Sc.sgd.db)
library(enrichplot)

# Pull latest master version of tsrexplorer
# devtools::install_github("rpolicastro/tsrexplorer", ref = "master", force = TRUE)

# Pull latest dev version of tsrexplorer
# devtools::install_github("rpolicastro/tsrexplorer", ref = "dev", force = TRUE)

if (!dir.exists("yeast_work/")){
  message("Creating directory 'yeast_work'...")
  dir.create("yeast_work/")
} else {
  message("Directory 'yeast_work' already exists...")
}

####################
### Read in TSSs ###
####################

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

# STRIPE-seq diamide TSSs
diamide_TSSs <- map(list.files("yeast_data/diamide_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                        makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                                 start.field = "TSS", end.field = "TSS")) %>%
    set_names(c("S288C_diamide_100ng_1","S288C_diamide_100ng_2","S288C_diamide_100ng_3"))

# Generate list of all TSS objects
full_TSSs_set <- c(YPD_TSSs,CAGE_TSSs,diamide_TSSs)

####################
### Read in TSRs ###
####################

# STRIPE-seq YPD TSRs
YPD_TSRs <- map(list.files("yeast_data/YPD_TSRs/", full.names = TRUE), ~ read.delim(.x) %>%
                      makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                               start.field = "start", end.field = "end")) %>%
    set_names(c("S288C_50ng_1","S288C_50ng_2","S288C_50ng_3",
                "S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
                "S288C_250ng_1","S288C_250ng_2","S288C_250ng_3"))

# SLIC-CAGE and nanoCAGE TSRs
CAGE_TSRs <- map(list.files("yeast_data/CAGE_TSRs/", full.names = TRUE), ~ read.delim(.x) %>%
                     makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                              start.field = "start", end.field = "end")) %>%
    set_names(c("SLIC_CAGE_100ng_1","SLIC_CAGE_100ng_2",
                "nanoCAGE_500ng_1","nanoCAGE_500ng_2",
                "nanoCAGE_25ng_1","nanoCAGE_25ng_2"))

# STRIPE-seq diamide TSRs
diamide_TSRs <- map(list.files("yeast_data/diamide_TSRs/", full.names = TRUE), ~ read.delim(.x) %>%
                        makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                                 start.field = "start", end.field = "end")) %>%
    set_names(c("S288C_diamide_100ng_1","S288C_diamide_100ng_2","S288C_diamide_100ng_3"))

# Generate list of all TSR objects
full_TSR_set <- c(YPD_TSRs,CAGE_TSRs,diamide_TSRs)

# Read in annotation files (included with package)
annotation <- system.file("extdata", "yeast_annotation.gtf", package="tsrexplorer")
assembly <- system.file("extdata", "yeast_assembly.fasta", package="tsrexplorer")

# Assign sample sets of interest to variables 

# YPD STRIPE-seq samples
stripe <- c("S288C_50ng_1","S288C_50ng_2","S288C_50ng_3",
            "S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
            "S288C_250ng_1","S288C_250ng_2","S288C_250ng_3")

# YPD 100 ng STRIPE-seq, SLIC-CAGE, and nanoCAGE samples
cage <- c("S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
          "SLIC_CAGE_100ng_1","SLIC_CAGE_100ng_2",
          "nanoCAGE_500ng_1","nanoCAGE_500ng_2",
          "nanoCAGE_25ng_1","nanoCAGE_25ng_2")

# All YPD STRIPE-seq and CAGE samples
all <- c("S288C_50ng_1","S288C_50ng_2","S288C_50ng_3",
         "S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
         "S288C_250ng_1","S288C_250ng_2","S288C_250ng_3",
         "SLIC_CAGE_100ng_1","SLIC_CAGE_100ng_2",
         "nanoCAGE_500ng_1","nanoCAGE_500ng_2",
         "nanoCAGE_25ng_1","nanoCAGE_25ng_2")

# YPD 100 ng and YPD + diamide STRIPE-seq samples
diamide <- c("S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
             "S288C_diamide_100ng_1","S288C_diamide_100ng_2","S288C_diamide_100ng_3")

# Read TSSs and TSRs into TSRexploreR
exp <- tsr_explorer(full_TSSs_set,full_TSR_set)

##############################################
### YPD input variation replicate analysis ###
##############################################

if (!dir.exists("yeast_work/YPD_STRIPE/")){
  message("Creating directory 'yeast_work/YPD_STRIPE' and changing working directory...")
  dir.create("yeast_work/YPD_STRIPE/")
  setwd("yeast_work/YPD_STRIPE/")
} else {
  message("Directory 'yeast_work/YPD_STRIPE' already exisits, changing working directory...")
  setwd("yeast_work/YPD_STRIPE/")
}

# Normalize TSS counts
exp <- count_normalization(exp, data_type = "tss", threshold = 3, n_samples = 1, samples = stripe)

# Generate a combinbed TSS correlation plot
p <- plot_correlation(exp, data_type = "tss", font_size = 2, pt_size = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 4))

ggsave("tss_correlation.png", plot = p, device = "png", type = "cairo", height = 5, width = 5)

# Generate a hierarchically clustered TSS heatmap with correlation values displayed
corr_matrix <- find_correlation(exp, data_type = "tss", correlation_metric = "pearson")

cairo_pdf(file = "tss_correlation_hierarchical.pdf", width = 13.5, height = 13.5)
Heatmap(corr_matrix, col = viridis(256), heatmap_legend_param = list(title = "PCC"), 
        layer_fun = function(j, i, x, y, width, height, fill)
        {
            grid.text(sprintf("%.3f", pindex(corr_matrix, i, j)), x, y, gp = gpar(fontsize = 28, col = "white"))
        }
)
dev.off()

# Annotate TSSs relative to genomic features
exp <- annotate_features(exp, annotation_file = annotation, data_type = "tss", feature_type = "transcript", 
                         upstream = 250, downstream = 100)

# Explore TSS read thresholds for promoter fraction and plot
thresh <- explore_thresholds(exp, annotation_file = annotation, feature_type = "transcript", max_threshold = 25, 
                                 upstream = 250, downstream = 100, samples = stripe)

p <- plot_threshold_exploration(thresh, ncol = 3, point_size = 1.5, sample_order = stripe) +
    geom_vline(xintercept = 3, lty = 2) +
    theme(legend.key.size = unit(0.8, "cm"), text = element_text(size = 12))

ggsave("tss_thresholds.pdf", plot = p, device = cairo_pdf, height = 5, width = 10)

# Determine TSS distribution relative to genomic features
tss_distribution <- genomic_distribution(exp, data_type = "tss", threshold = 3, samples = stripe)

p <- plot_genomic_distribution(tss_distribution, sample_order = stripe) +
    ggplot2::theme(text = element_text(size = 14), legend.key.size = unit(0.8, "cm"))

ggsave("tss_genomic_distribution.pdf", plot = p, device = cairo_pdf, height = 6, width = 6)

genomic_dist <- genomic_distribution(exp, data_type = "tss", threshold = 3, quantiles = 5, 
                                     samples = "S288C_100ng_1")

p <- plot_genomic_distribution(genomic_dist, sample_order = stripe) +
    ggplot2::theme(text = element_text(size = 14), legend.key.size = unit(0.8, "cm"))

ggsave("tss_genomic_distribution_quantiles.pdf", plot = p, device = cairo_pdf, height = 4.5, width = 6)

# Plot number of promoter-proximal features with a TSS
features <- detect_features(exp, data_type = "tss", feature_type = "transcript", threshold = 3, 
                            samples = stripe)

p <- plot_detected_features(features, ncol = 3, width = 0.75) +
    ggplot2::theme(text = element_text(size = 14), legend.key.size = unit(0.8, "cm"))

ggsave("tss_feature_plot.pdf", plot = p, device = cairo_pdf, height = 4, width = 8)

# Generate TSS density plots
p <- plot_average(exp, data_type = "tss", threshold = 3, samples = "S288C_100ng_1", upstream = 1000, downstream = 1000) +
    ggplot2::theme(text = element_text(size = 13))

ggsave("tss_average_plot.pdf", plot = p, cairo_pdf, height = 2.5, width = 3.5)

# Generate TSS sequence logos
seqs <- tss_sequences(exp, genome_assembly = assembly, threshold = 3, samples = stripe)

p <- plot_sequence_logo(seqs, ncol = 3, font_size = 10)

ggsave("tss_seq_logo.pdf", plot = p, device = cairo_pdf, height = 5, width = 12)

seqs <- tss_sequences(exp, genome_assembly = assembly, threshold = 3, quantiles = 5, samples = "S288C_100ng_1")

p <- plot_sequence_logo(seqs, ncol = 1, font_size = 10)

ggsave("tss_seq_logo_quantiles.pdf", plot = p, device = cairo_pdf, height = 7, width = 5)

# Generate TSS color plot
seqs <- tss_sequences(exp, genome_assembly = assembly, threshold = 3, samples = "S288C_100ng_1")

p <- plot_sequence_colormap(seqs, ncol = 3) +
    ggplot2::theme(text = element_text(size = 6), legend.key.size = unit(0.4, "cm"))

ggsave("tss_seq_colormap.png", plot = p, device = "png", type = "cairo", height = 2.5, width = 2)

# Assess TSS dinucleotide frequencies
frequencies <- dinucleotide_frequencies(exp, genome_assembly = assembly, threshold = 3, samples = stripe)

p <- plot_dinucleotide_frequencies(frequencies, ncol = 3, sample_order = stripe) +
    ggplot2::theme(text = element_text(size = 6), legend.key.size = unit(0.4, "cm"))

ggsave("tss_dinucleotide_frequencies.pdf", plot = p, device = cairo_pdf, height = 7, width = 8)

# Plot distance of dominant TSS to annotated start codon
dominant <- dominant_tss(exp, threshold = 3, feature_type = "geneId", samples = "S288C_100ng_1")

p <- plot_dominant_tss(dominant)

ggsave("dominant_tss.pdf", plot = p, device = cairo_pdf, height = 4, width = 4)

# Plot hypothetical maximum 5'UTR length
max <- max_utr(exp, threshold = 3, feature_type = "geneId", samples = "S288C_100ng_1")

p <- plot_max_utr(max)

ggsave("max_utr.pdf", plot = p, device = cairo_pdf, height = 4, width = 4)

# Export normalized TSS bedGraphs

if (!dir.exists(file.path("yeast_work", "bedgraphs"))){
  message("Creating directory 'yeast_work/bedgraphs'...")
  dir.create(file.path("yeast_work", "bedgraphs"))
} else {
  message("Directory 'yeast_work/bedgraphs' already exists...")
}

iwalk(exp@counts$TSSs$cpm, function(counts, sample) {
	pos <- counts[strand(counts) == "+"]
	min <- counts[strand(counts) == "-"]

	export(pos, file.path("yeast_work", "bedgraphs", paste(sample, "pos.bedgraph", sep = "_")), format = "bedgraph")
	export(min, file.path("yeast_work", "bedgraphs", paste(sample, "min.bedgraph", sep = "_")), format = "bedgraph")
})

# Normalize TSR counts
exp <- count_normalization(exp, data_type = "tsr", threshold = 3, n_samples = 1, samples = stripe)

# Generate a combinbed TSR correlation plot
p <- plot_correlation(exp, data_type = "tsr", font_size = 2, pt_size = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 4))

ggsave("tsr_correlation.png", plot = p, device = "png", type = "cairo", height = 5, width = 5)

# Generate a hierarchically clustered TSR heatmap with correlation values displayed
corr_matrix <- find_correlation(exp, data_type = "tsr", correlation_metric = "pearson")

cairo_pdf(file = "tsr_correlation_hierarchical.pdf", width = 13.5, height = 13.5)
Heatmap(corr_matrix, col = viridis(256), heatmap_legend_param = list(title = "PCC"), 
        layer_fun = function(j, i, x, y, width, height, fill)
        {
            grid.text(sprintf("%.3f", pindex(corr_matrix, i, j)), x, y, gp = gpar(fontsize = 28, col = "white"))
        }
)
dev.off()

# Annotate TSRs
exp <- annotate_features(exp, annotation_file = annotation, data_type = "tsr", feature_type = "transcript")

# Determine TSR distribution relative to genomic features
tsr_distribution <- genomic_distribution(exp, data_type = "tsr", threshold = 3, samples = "S288C_100ng_1")

p <- plot_genomic_distribution(tsr_distribution) +
    ggplot2::theme(text = element_text(size = 12), legend.key.size = unit(0.6, "cm"))

ggsave("tsr_genomic_distribution.pdf", plot = p, device = cairo_pdf, height = 2, width = 6)

tsr_distribution <- genomic_distribution(exp, data_type = "tsr", threshold = 3, quantiles = 5, 
                                         samples = stripe)

p <- plot_genomic_distribution(tsr_distribution, sample_order = stripe) +
    ggplot2::theme(text = element_text(size = 14), legend.key.size = unit(0.8, "cm"))

ggsave("tsr_genomic_distribution_quantiles.pdf", plot = p, device = cairo_pdf, height = 12, width = 6)

# Plot number of promoter-proximal features with a TSR
features <- detect_features(exp, data_type = "tsr", feature_type = "transcript", samples = stripe)

p <- plot_detected_features(features, ncol = 3, width = 0.75) +
    ggplot2::theme(text = element_text(size = 14), legend.key.size = unit(0.8, "cm"))

ggsave("tsr_feature_plot.pdf", plot = p, device = cairo_pdf, height = 4.5, width = 8)

# # Plot selected TSR metrics
# p <- plot_tsr_metric(exp, tsr_metrics = "nTSSs", log2_transform = TRUE, ncol = 1, plot_type = "boxjitter", 
#                      size = 0.5, samples = "S288C_100ng_1") +
#     ggplot2::theme(text = element_text(size = 6))
# 
# ggsave("tsr_metrics.pdf", plot = p, device = cairo_pdf, width = 7, height = 7)

# Generate TSR density plot
p <- plot_average(exp, data_type = "tsr", samples = "S288C_100ng_1", upstream = 1000, downstream = 1000) +
    ggplot2::theme(text = element_text(size = 13))

ggsave("tsr_density_plot.pdf", plot = p, device = cairo_pdf, height = 2.5, width = 3.5)

####################################
### STRIPE-seq vs. CAGE analysis ###
####################################

setwd("..")

if (!dir.exists(file.path("yeast_work", "CAGE"))){
  message("Creating directory 'yeast_work/CAGE' and changing working directory...")
  dir.create(file.path("yeast_work", "CAGE"))
  setwd(file.path("yeast_work", "CAGE"))
} else {
  message("Directory 'yeast_work/CAGE' already exists, changing working directory...")
  setwd(file.path("yeast_work", "CAGE"))
}

# Normalize TSS counts
exp <- count_normalization(exp, data_type = "tss", threshold = 3, n_samples = 1, samples = cage)

# Annotate TSSs relative to genomic features
exp <- annotate_features(exp, annotation_file = annotation, data_type = "tss", feature_type = "transcript", 
                         upstream = 250, downstream = 100)

# Explore TSS read thresholds for promoter fraction and plot
thresh <- explore_thresholds(exp, annotation_file = annotation, feature_type = "transcript", max_threshold = 25, 
                             upstream = 250, downstream = 100, samples = cage)

p <- plot_threshold_exploration(thresh, ncol = 3, point_size = 2, sample_order = cage) +
    ggplot2::geom_vline(xintercept = 3, lty = 2) + 
    ggplot2::theme(legend.key.size = unit(0.8, "cm"), text = element_text(size = 12))

ggsave("tss_thresholds.pdf", plot = p, device = cairo_pdf, height = 5, width = 10)

# Determine TSS distribution relative to genomic features
tss_distribution <- genomic_distribution(exp, data_type = "tss", threshold = 3, samples = cage)

p <- plot_genomic_distribution(tss_distribution, sample_order = cage) +
    ggplot2::theme(text = element_text(size = 14), legend.key.size = unit(0.8, "cm"))

ggsave("tss_genomic_distribution.pdf", plot = p, device = cairo_pdf, height = 6, width = 6)

# Assess TSS dinucleotide frequencies
frequencies <- dinucleotide_frequencies(exp, genome_assembly = assembly, threshold = 3, 
                                        samples = cage)

p <- plot_dinucleotide_frequencies(frequencies, ncol = 3, sample_order = cage) +
  ggplot2::theme(text = element_text(size = 6), legend.key.size = unit(0.4, "cm"))

ggsave("tss_dinucleotide_frequencies.pdf", plot = p, device = cairo_pdf, height = 7, width = 8)

# Export normalized TSS bedGraphs

# SLIC-CAGE
export.bedGraph(exp@counts$TSSs$cpm$SLIC_CAGE_100ng_1[strand(exp@counts$TSSs$cpm$SLIC_CAGE_100ng_1) == "+"], 
                file.path(baseDir, "yeast_work/bedgraphs/SLIC_CAGE_100ng_1_+.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$SLIC_CAGE_100ng_1[strand(exp@counts$TSSs$cpm$SLIC_CAGE_100ng_1) == "-"], 
                file.path(baseDir, "yeast_work/bedgraphs/SLIC_CAGE_100ng_1_-.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$SLIC_CAGE_100ng_2[strand(exp@counts$TSSs$cpm$SLIC_CAGE_100ng_2) == "+"], 
                file.path(baseDir, "yeast_work/bedgraphs/SLIC_CAGE_100ng_2_+.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$SLIC_CAGE_100ng_2[strand(exp@counts$TSSs$cpm$SLIC_CAGE_100ng_2) == "-"], 
                file.path(baseDir, "yeast_work/bedgraphs/SLIC_CAGE_100ng_2_-.bedgraph"))

# nanoCAGE 500 ng
export.bedGraph(exp@counts$TSSs$cpm$nanoCAGE_500ng_1[strand(exp@counts$TSSs$cpm$nanoCAGE_500ng_1) == "+"], 
                file.path(baseDir, "yeast_work/bedgraphs/nanoCAGE_500ng_1_+.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$nanoCAGE_500ng_1[strand(exp@counts$TSSs$cpm$nanoCAGE_500ng_1) == "-"], 
                file.path(baseDir, "yeast_work/bedgraphs/nanoCAGE_500ng_1_-.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$nanoCAGE_500ng_2[strand(exp@counts$TSSs$cpm$nanoCAGE_500ng_2) == "+"], 
                file.path(baseDir, "yeast_work/bedgraphs/nanoCAGE_500ng_2_+.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$nanoCAGE_500ng_2[strand(exp@counts$TSSs$cpm$nanoCAGE_500ng_2) == "-"], 
                file.path(baseDir, "yeast_work/bedgraphs/nanoCAGE_500ng_2_-.bedgraph"))

# nanoCAGE 25 ng
export.bedGraph(exp@counts$TSSs$cpm$nanoCAGE_25ng_1[strand(exp@counts$TSSs$cpm$nanoCAGE_25ng_1) == "+"], 
                file.path(baseDir, "yeast_work/bedgraphs/nanoCAGE_25ng_1_+.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$nanoCAGE_25ng_1[strand(exp@counts$TSSs$cpm$nanoCAGE_25ng_1) == "-"], 
                file.path(baseDir, "yeast_work/bedgraphs/nanoCAGE_25ng_1_-.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$nanoCAGE_25ng_2[strand(exp@counts$TSSs$cpm$nanoCAGE_25ng_2) == "+"], 
                file.path(baseDir, "yeast_work/bedgraphs/nanoCAGE_25ng_2_+.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$nanoCAGE_25ng_2[strand(exp@counts$TSSs$cpm$nanoCAGE_25ng_2) == "-"], 
                file.path(baseDir, "yeast_work/bedgraphs/nanoCAGE_25ng_2_-.bedgraph"))

# Normalize TSR counts
exp <- count_normalization(exp, data_type = "tsr", threshold = 3, n_samples = 1, samples = all)

# Generate a combinbed TSR correlation plot
p <- plot_correlation(exp, data_type = "tsr", correlation_metric = "spearman") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_correlation.png", plot = p, device = "png", type = "cairo", height = 30, width = 30)

# Generate a hierarchically clustered TSR heatmap with correlation values displayed
corr_matrix <- find_correlation(exp, data_type = "tsr", correlation_metric = "spearman")

cairo_pdf(file = "tsr_correlation_hierarchical.pdf", width = 22.5, height = 22.5)
Heatmap(corr_matrix, col = viridis(256), heatmap_legend_param = list(title = "Spearman"), 
        layer_fun = function(j, i, x, y, width, height, fill)
        {
            grid.text(sprintf("%.3f", pindex(corr_matrix, i, j)), x, y, gp = gpar(fontsize = 28, col = "white"))
        }
)
dev.off()

# Annotate TSRs
exp <- annotate_features(exp, annotation_file = annotation, data_type = "tsr", feature_type = "transcript")

# Determine TSR distribution relative to genomic features
tsr_distribution <- genomic_distribution(exp, data_type = "tsr", threshold = 3, 
                                         samples = cage)

p <- plot_genomic_distribution(tsr_distribution, sample_order = cage) +
    ggplot2::theme(text = element_text(size = 6), legend.key.size = unit(0.4, "cm"))

ggsave("tsr_genomic_distribution.pdf", plot = p, device = cairo_pdf, height = 3, width = 4)

# Plot number of promoter-proximal features with a TSR
features <- detect_features(exp, data_type = "tsr", feature_type = "transcript", samples = cage)

p <- plot_detected_features(features, ncol = 3) +
    ggplot2::theme(text = element_text(size = 5), legend.key.size = unit(0.4, "cm"))

ggsave("tsr_feature_plot.pdf", plot = p, device = cairo_pdf, height = 5, width = 6)

########################
### Diamide analysis ###
########################

if (!dir.exists(file.path(baseDir, "yeast_work/diamide/"))){
  print("Creating directory 'yeast_work/diamide' and changing working directory...")
  dir.create(file.path(baseDir, "yeast_work/diamide/"))
  setwd(file.path(baseDir, "yeast_work/diamide/"))
} else {
  print("Directory 'yeast_work/diamide' already exists, changing working directory...")
  setwd(file.path(baseDir, "yeast_work/diamide/"))
}
# Normalize TSS counts
exp <- count_normalization(exp, data_type = "tss", threshold = 3, n_samples = 1, samples = diamide)

# Generate a hierarchically clustered TSS heatmap with correlation values displayed
corr_matrix <- find_correlation(exp, data_type = "tss", correlation_metric = "pearson")

cairo_pdf(file = "tss_correlation_hierarchical.pdf", width = 9, height = 9)
Heatmap(corr_matrix, col = viridis(256), heatmap_legend_param = list(title = "PCC"), 
        layer_fun = function(j, i, x, y, width, height, fill)
        {
            grid.text(sprintf("%.3f", pindex(corr_matrix, i, j)), x, y, gp = gpar(fontsize = 24, col = "white"))
        }
)
dev.off()

# Annotate TSSs relative to genomic features
exp <- annotate_features(exp, annotation_file = annotation, data_type = "tss", feature_type = "transcript", 
                         upstream = 250, downstream = 100)

# Determine TSS distribution relative to genomic features
tss_distribution <- genomic_distribution(exp, data_type = "tss", threshold = 3, 
                                         samples = diamide)

p <- plot_genomic_distribution(tss_distribution, sample_order = diamide) +
    ggplot2::theme(text = element_text(size = 6), legend.key.size = unit(0.4, "cm"))

ggsave("tss_genomic_distribution.pdf", plot = p, device = cairo_pdf, height = 3, width = 4)

# Plot number of promoter-proximal features with a TSS
features <- detect_features(exp, data_type = "tss", feature_type = "transcript", threshold = 3, 
                            samples = c("S288C_100ng_1","S288C_diamide_100ng_1"))

p <- plot_detected_features(features, ncol = 3) +
    ggplot2::theme(text = element_text(size = 5), legend.key.size = unit(0.4, "cm"))

ggsave("tss_feature_plot.pdf", plot = p, device = cairo_pdf, height = 2, width = 3)

# Generate TSS sequence logos
seqs <- tss_sequences(exp, genome_assembly = assembly, threshold = 3, 
                      samples = c("S288C_100ng_1","S288C_diamide_100ng_1"))

p <- plot_sequence_logo(seqs, ncol = 1) +
    ggplot2::theme(text = element_text(size = 5))

ggsave("tss_seq_logo.pdf", plot = p, device = cairo_pdf, height = 2, width = 3)

# Plot distance of dominant TSS to annotated start codon
dominant <- dominant_tss(exp, threshold = 3, feature_type = "geneId", samples = c("S288C_100ng_1","S288C_diamide_100ng_1"))

p <- plot_dominant_tss(dominant)

ggsave("dominant_tss.pdf", plot = p, device = cairo_pdf, height = 4, width = 4)

# Plot hypothetical maximum 5'UTR length
max <- max_utr(exp, threshold = 3, feature_type = "geneId", samples = c("S288C_100ng_1","S288C_diamide_100ng_1"))

p <- plot_max_utr(max)

ggsave("max_utr.pdf", plot = p, device = cairo_pdf, height = 4, width = 4)

# Export normalized TSS bedGraphs
export.bedGraph(exp@counts$TSSs$cpm$S288C_diamide_100ng_1[strand(exp@counts$TSSs$cpm$S288C_diamide_100ng_1) == "+"], 
                file.path(baseDir, "yeast_work/bedgraphs/S288C_diamide_100ng_1_+.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$S288C_diamide_100ng_1[strand(exp@counts$TSSs$cpm$S288C_diamide_100ng_1) == "-"], 
                file.path(baseDir, "yeast_work/bedgraphs/S288C_diamide_100ng_1_-.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$S288C_diamide_100ng_2[strand(exp@counts$TSSs$cpm$S288C_diamide_100ng_2) == "+"], 
                file.path(baseDir, "yeast_work/bedgraphs/S288C_diamide_100ng_2_+.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$S288C_diamide_100ng_2[strand(exp@counts$TSSs$cpm$S288C_diamide_100ng_2) == "-"], 
                file.path(baseDir, "yeast_work/bedgraphs/S288C_diamide_100ng_2_-.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$S288C_diamide_100ng_3[strand(exp@counts$TSSs$cpm$S288C_diamide_100ng_3) == "+"], 
                file.path(baseDir, "yeast_work/bedgraphs/S288C_diamide_100ng_3_+.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$S288C_diamide_100ng_3[strand(exp@counts$TSSs$cpm$S288C_diamide_100ng_3) == "-"], 
                file.path(baseDir, "yeast_work/bedgraphs/S288C_diamide_100ng_3_-.bedgraph"))

# Normalize TSR counts
exp <- count_normalization(exp, data_type = "tsr", threshold = 3, n_samples = 1, samples = diamide)

# Generate a combinbed TSR correlation plot
p <- plot_correlation(exp, data_type = "tsr") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_correlation.png", plot = p, device = "png", type = "cairo", height = 8, width = 8)

# Generate a hierarchically clustered TSS heatmap with correlation values displayed
corr_matrix <- find_correlation(exp, data_type = "tsr", correlation_metric = "pearson")

cairo_pdf(file = "tsr_correlation_hierarchical.pdf", width = 9, height = 9)
Heatmap(corr_matrix, col = viridis(256), heatmap_legend_param = list(title = "PCC"), 
        layer_fun = function(j, i, x, y, width, height, fill)
        {
            grid.text(sprintf("%.3f", pindex(corr_matrix, i, j)), x, y, gp = gpar(fontsize = 24, col = "white"))
        }
)
dev.off()

# Annotate TSRs
exp <- annotate_features(exp, annotation_file = annotation, data_type = "tsr", feature_type = "transcript")

# Determine TSR distribution relative to genomic features
tsr_distribution <- genomic_distribution(exp, data_type = "tsr", threshold = 3, 
                                         samples = diamide)

p <- plot_genomic_distribution(tsr_distribution, sample_order = diamide) +
    ggplot2::theme(text = element_text(size = 6), legend.key.size = unit(0.4, "cm"))

ggsave("tsr_genomic_distribution.pdf", plot = p, device = cairo_pdf, height = 3, width = 4)

# Differential TSR (dTSR) analysis
edger_model <- fit_edger_model(
    exp,
    data_type = "tsr", 
    samples = c(
        "S288C_100ng_1",
        "S288C_100ng_2",
        "S288C_100ng_3",
        "S288C_diamide_100ng_1",
        "S288C_diamide_100ng_2",
        "S288C_diamide_100ng_3"
    ),
    groups = c(1, 1, 1, 2, 2, 2)
)

diff_tsrs <- differential_expression(edger_model, data_type = "tsr", compare_groups = c(1, 2))

# Write dTSRs to a table
diff_tsrs %>% 
    filter(log2FC <= -1 & FDR < 0.05 | log2FC >= 1 & FDR < 0.05) %>%
    write.table(., "diff_tsrs.tsv", sep="\t", col.names=T, row.names=F, quote=F)

# Annotate dTSRs
annotated_diff_tsrs <- annotate_differential_tsrs(diff_tsrs, annotation_file = annotation, 
                                                  feature_type = "transcript", 
                                                  upstream = 250, downstream = 100)

# Write annotated significant dTSRs to a table
annotated_diff_tsrs %>% 
    filter(annotation == "Promoter" & log2FC <= -1 & FDR < 0.05 | 
           annotation == "Promoter" & log2FC >= 1 & FDR < 0.05) %>%
    write.table(., "promoter_annotated_diff_tsrs.tsv", sep="\t", col.names=T, row.names=F, quote=F)

# Make a volcano plot of dTSRs
p <- plot_volcano(diff_tsrs, size = 0.1) + 
    ggplot2::theme(legend.key.size = unit(0.4, "cm"))

ggsave("diff_tsrs_volcano_plot.pdf", plot = p, device = cairo_pdf, height = 2.5, width = 4)

# Perform GO analysis (work in progress)
enrichment_data <- export_for_enrichment(annotated_diff_tsrs)

library(org.Sc.sgd.db)

go_enrichment <- clusterProfiler::compareCluster(
    geneId ~ change,
    data = enrichment_data,
    fun = "enrichGO",
    OrgDb = "org.Sc.sgd.db",
    pAdjustMethod = "fdr"
)

##########################################
### RNA-seq analysis (YPD and diamide) ###
##########################################

if (!dir.exists(file.path(baseDir, "yeast_work/RNA_seq/"))){
    print("Creating directory 'yeast_work/RNA_seq' and changing working directory...")
    dir.create(file.path(baseDir, "yeast_work/RNA_seq/"))
    setwd(file.path(baseDir, "yeast_work/RNA_seq/"))
} else {
    print("Directory 'yeast_work/RNA_seq' already exisits, changing working directory...")
    setwd(file.path(baseDir, "yeast_work/RNA_seq/"))
}

# Get feature counts for STRIPE-seq and RNA-seq
stripe_counts <- read_tsv(file.path(baseDir, "yeast_data/RNA_seq/cleaned_S288C_feature_counts_STRIPEseq.tsv"))
rnaseq_counts <- read_tsv(file.path(baseDir, "yeast_data/RNA_seq/cleaned_yeast_feature_counts.tsv"))

stripe_counts <- column_to_rownames(stripe_counts, "Geneid") %>% as.matrix
rnaseq_counts <- column_to_rownames(rnaseq_counts, "Geneid") %>% as.matrix

exp <- add_feature_counts(exp, five_prime_feature_counts = stripe_counts, rnaseq_feature_counts = rnaseq_counts)
exp <- count_normalization(exp, data_type = "features")

# Generate sample sets of interest for feature count correlation analysis

# YPD STRIPE-seq and RNA-seq
stripe_rnaseq_ypd <- c("RNASEQ001_S288C_untreated_r1", "RNASEQ002_S288C_untreated_r1", "RNASEQ003_S288C_untreated_r1",
                       "GSF2268_s_SP01_S288C_WT_50ng", "GSF2268_s_SP02_S288C_WT_50ng", "GSF2268_s_SP03_S288C_WT_50ng",
                       "GSF2268_s_SP04_S288C_WT_100ng",	"GSF2268_s_SP05_S288C_WT_100ng", "GSF2268_s_SP06_S288C_WT_100ng",
                       "GSF2268_s_SP07_S288C_WT_250ng", "GSF2268_s_SP08_S288C_WT_250ng", "GSF2268_s_SP09_S288C_WT_250ng")

# Diamide STRIPE-seq and RNA-seq
stripe_rnaseq_diamide <- c("RNASEQ001_S288C_untreated_r1", "RNASEQ002_S288C_untreated_r1", "RNASEQ003_S288C_untreated_r1",
                           "GSF2268_s_SP04_S288C_WT_100ng",	"GSF2268_s_SP05_S288C_WT_100ng", "GSF2268_s_SP06_S288C_WT_100ng",
                           "GSF2268_s_SP28_S288C_Diamide_100ng", "GSF2268_s_SP29_S288C_Diamide_100ng", "GSF2268_s_SP30_S288C_Diamide_100ng",
                           "RNASEQ007_S288C_diamide_r2", "RNASEQ008_S288C_diamide_r2", "RNASEQ009_S288C_diamide_r2")

# Find correlation of YPD STRIPE-seq and RNA-seq feature counts
corr_matrix <- find_correlation(exp, data_type = "features", correlation_metric = "spearman", samples = stripe_rnaseq_ypd)

cairo_pdf(file = "tss_rnaseq_correlation_ypd_hierarchical.pdf", width = 18, height = 18)
Heatmap(corr_matrix, col = viridis(256), heatmap_legend_param = list(title = "spearman"), 
        layer_fun = function(j, i, x, y, width, height, fill)
        {
            grid.text(sprintf("%.3f", pindex(corr_matrix, i, j)), x, y, gp = gpar(fontsize = 26, col = "white"))
        }
)
dev.off()

# Find correlation of diamide STRIPE-seq and RNA-seq feature counts
corr_matrix <- find_correlation(exp, data_type = "features", correlation_metric = "spearman", samples = stripe_rnaseq_diamide)

cairo_pdf(file = "tss_rnaseq_correlation_diamide_hierarchical.pdf", width = 18, height = 18)
Heatmap(corr_matrix, col = viridis(256), heatmap_legend_param = list(title = "spearman"), 
        layer_fun = function(j, i, x, y, width, height, fill)
        {
            grid.text(sprintf("%.3f", pindex(corr_matrix, i, j)), x, y, gp = gpar(fontsize = 26, col = "white"))
        }
)
dev.off()

# Perform YPD vs. diamide differential expression for STRIPE-seq and RNA-seq

# Identify differential STRIPE-seq features
edger_model_stripe <- fit_edger_model(
    exp,
    data_type = "features",
    samples = c(
        "GSF2268_s_SP04_S288C_WT_100ng",
        "GSF2268_s_SP05_S288C_WT_100ng",
        "GSF2268_s_SP06_S288C_WT_100ng",
        "GSF2268_s_SP28_S288C_Diamide_100ng",
        "GSF2268_s_SP29_S288C_Diamide_100ng",
        "GSF2268_s_SP30_S288C_Diamide_100ng"
    ),
    groups = c(1, 1, 1, 2, 2, 2)
)

diff_features_stripe <- differential_expression(edger_model_stripe, data_type = "features", compare_groups = c(1, 2)) %>%
    dplyr::rename("geneId" = "position")

# Write differential STRIPE-seq features to a table
diff_features_stripe %>% 
    filter(log2FC <= -1 & FDR < 0.05 | log2FC >= 1 & FDR < 0.05) %>%
    write.table(., "diff_features_stripe.tsv", sep="\t", col.names=T, row.names=F, quote=F)

# Count differential STRIPE-seq features
diff_features_stripe %>%
    count(log2FC >= 1 & FDR < 0.05)
diff_features_stripe %>%
    count(log2FC <= -1 & FDR < 0.05)

# Export STRIPE-seq results for GO enrichment, adding an identifying column for compareCluster
enrichment_data_stripe <- export_for_enrichment(diff_features_stripe, log2fc_cutoff = 1, fdr_cutoff = 0.05) %>%
    add_column(technology = "STRIPE")

# Identify differential RNA-seq features
edger_model_rnaseq <- fit_edger_model(
    exp,
    data_type = "features",
    samples = c(
        "RNASEQ001_S288C_untreated_r1", 
        "RNASEQ002_S288C_untreated_r1", 
        "RNASEQ003_S288C_untreated_r1",
        "RNASEQ007_S288C_diamide_r2", 
        "RNASEQ008_S288C_diamide_r2", 
        "RNASEQ009_S288C_diamide_r2"
    ),
    groups = c(1, 1, 1, 2, 2, 2)
)

diff_features_rnaseq <- differential_expression(edger_model_rnaseq, data_type = "features", compare_groups = c(1, 2)) %>%
    dplyr::rename("geneId" = "position")

# Write differential RNA-seq features to a table
diff_features_rnaseq %>% 
    filter(log2FC <= -1 & FDR < 0.05 | log2FC >= 1 & FDR < 0.05) %>%
    write.table(., "diff_features_rnaseq.tsv", sep="\t", col.names=T, row.names=F, quote=F)

# Count differential RNA-seq features
diff_features_rnaseq %>%
    count(log2FC >= 1 & FDR < 0.05)
diff_features_rnaseq %>%
    count(log2FC <= -1 & FDR < 0.05)

# Export RNA-seq results for GO enrichment, adding an identifying column for compareCluster
enrichment_data_rnaseq <- export_for_enrichment(diff_features_rnaseq, log2fc_cutoff = 1, fdr_cutoff = 0.05) %>%
    add_column(technology = "RNA-seq")

# Combine enrichment datasets and perform GO analysis
enrichment_data_combined <- bind_rows(enrichment_data_stripe, enrichment_data_rnaseq)

go_enrichment_combined <- clusterProfiler::compareCluster(
    geneId ~ technology + change,
    data = enrichment_data_combined,
    fun = "enrichGO",
    OrgDb = "org.Sc.sgd.db",
    pAdjustMethod = "fdr",
    keyType = "ENSEMBL",
    ont = "BP"
)

cairo_pdf(file = "go_analysis_combined.pdf", height = 10, width = 8)
dotplot(go_enrichment_combined, showCategory = 10) + scale_color_viridis()
dev.off()
