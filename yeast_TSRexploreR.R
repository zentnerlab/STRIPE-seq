library(tsrexplorer)
library(tidyverse)
library(GenomicRanges)
library(viridis)
library(ComplexHeatmap)
library(rtracklayer)
library(eulerr)
library(org.Sc.sgd.db)
library(enrichplot)

# This script runs the majority of TSRexploreR functions and was used to generate the majority of images for the STRIPE-seq manuscript. 
# It is divided into sections for analysis of YPD input replicates, comparison of STRIPE-seq and CAGE, analysis of differential TSR 
# usage upon diamide treatment, and comparison of STRIPE-seq and RNA-seq for transcript abundance and differential expression analysis.

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

cage_slim <- c("S288C_100ng_1",
               "SLIC_CAGE_100ng_1",
               "nanoCAGE_500ng_1",
               "nanoCAGE_25ng_1")

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

ypd_stripe_dir <- file.path("yeast_work", "YPD_STRIPE")

if (!dir.exists(ypd_stripe_dir)) {
  message("Creating directory 'yeast_work/YPD_STRIPE'...")
  dir.create(ypd_stripe_dir)
} else {
  message("Directory 'yeast_work/YPD_STRIPE' already exists.")
}

# Normalize TSS counts
exp <- count_normalization(exp, data_type = "tss", threshold = 3, n_samples = 1, samples = stripe)

# Generate a hierarchically clustered TSS heatmap with correlation values displayed
corr_matrix <- find_correlation(exp, data_type = "tss", correlation_metric = "pearson", samples = stripe)

cairo_pdf(file = file.path(ypd_stripe_dir, "tss_correlation_hierarchical.pdf"), width = 13.5, height = 13.5)
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
    theme(legend.key.size = unit(0.8, "cm"), text = element_text(size = 12), axis.text = element_text(color="black"))

ggsave(file.path(ypd_stripe_dir, "tss_thresholds.pdf"), plot = p, device = cairo_pdf, height = 5, width = 10)

# Export threshold plot for one sample for main figure
thresh <- explore_thresholds(exp, annotation_file = annotation, feature_type = "transcript", max_threshold = 25, 
                             upstream = 250, downstream = 100, samples = "S288C_100ng_1")

p <- plot_threshold_exploration(thresh, ncol = 1, point_size = 4.5) +
    geom_vline(xintercept = 3, lty = 2) +
    theme(legend.key.size = unit(0.8, "cm"), text = element_text(size = 12), axis.text = element_text(color="black"))

ggsave(file.path(ypd_stripe_dir, "tss_thresholds_S288C_100ng_1.pdf"), plot = p, device = cairo_pdf, height = 5, width = 10)

# Determine TSS distribution relative to genomic features
tss_distribution <- genomic_distribution(exp, data_type = "tss", threshold = 3, samples = stripe)

p <- plot_genomic_distribution(tss_distribution, sample_order = stripe) +
    ggplot2::theme(text = element_text(size = 14), legend.key.size = unit(0.8, "cm"), axis.text = element_text(color="black"))

ggsave(file.path(ypd_stripe_dir, "tss_genomic_distribution.pdf"), plot = p, device = cairo_pdf, height = 6, width = 6)

genomic_dist <- genomic_distribution(exp, data_type = "tss", threshold = 3, quantiles = 5, 
                                     samples = "S288C_100ng_1")

p <- plot_genomic_distribution(genomic_dist, sample_order = stripe) +
    ggplot2::theme(text = element_text(size = 14), legend.key.size = unit(0.8, "cm"), axis.text = element_text(color="black"))

ggsave(file.path(ypd_stripe_dir, "tss_genomic_distribution_quantiles.pdf"), plot = p, device = cairo_pdf, height = 4.5, width = 6)

# Plot number of promoter-proximal features with a TSS
features <- detect_features(exp, data_type = "tss", feature_type = "transcript", threshold = 3, 
                            samples = stripe)

p <- plot_detected_features(features, ncol = 3, width = 0.75) +
    ggplot2::theme(text = element_text(size = 14), legend.key.size = unit(0.8, "cm"), axis.text = element_text(color="black"))

ggsave(file.path(ypd_stripe_dir, "tss_feature_plot.pdf"), plot = p, device = cairo_pdf, height = 4, width = 8)

# Generate TSS density plots
p <- plot_average(exp, data_type = "tss", threshold = 3, samples = "S288C_100ng_1", upstream = 1000, downstream = 1000) +
    ggplot2::theme(text = element_text(size = 13), axis.text = element_text(color="black"))

ggsave(file.path(ypd_stripe_dir, "tss_average_plot.pdf"), plot = p, cairo_pdf, height = 2.5, width = 3.5)

# Generate TSS sequence logos
seqs <- tss_sequences(exp, genome_assembly = assembly, threshold = 3, samples = stripe)

p <- plot_sequence_logo(seqs, ncol = 3, font_size = 10)

ggsave(file.path(ypd_stripe_dir, "tss_seq_logo.pdf"), plot = p, device = cairo_pdf, height = 5, width = 12)

seqs <- tss_sequences(exp, genome_assembly = assembly, threshold = 3, quantiles = 5, samples = "S288C_100ng_1")

p <- plot_sequence_logo(seqs, ncol = 1, font_size = 10)

ggsave(file.path(ypd_stripe_dir, "tss_seq_logo_quantiles.pdf"), plot = p, device = cairo_pdf, height = 7, width = 5)

# Generate TSS color plot
seqs <- tss_sequences(exp, genome_assembly = assembly, threshold = 3, samples = "S288C_100ng_1")

p <- plot_sequence_colormap(seqs, ncol = 3) +
    ggplot2::theme(text = element_text(size = 6), legend.key.size = unit(0.4, "cm"), axis.text = element_text(color="black"))

ggsave(file.path(ypd_stripe_dir, "tss_seq_colormap.png"), plot = p, device = "png", type = "cairo", height = 2.5, width = 2)

# Assess TSS dinucleotide frequencies
frequencies <- dinucleotide_frequencies(exp, genome_assembly = assembly, threshold = 3, samples = stripe)

p <- plot_dinucleotide_frequencies(frequencies, ncol = 3, sample_order = stripe) +
    ggplot2::theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"), axis.text = element_text(color="black"))

ggsave(file.path(ypd_stripe_dir, "tss_dinucleotide_frequencies.pdf"), plot = p, device = cairo_pdf, height = 7, width = 7)

# Plot single-sample TSS dinucleotide frequencies
frequencies <- dinucleotide_frequencies(exp, genome_assembly = assembly, threshold = 3, samples = "S288C_100ng_1")

p <- plot_dinucleotide_frequencies(frequencies, ncol = 1) +
    ggplot2::theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"), axis.text = element_text(color="black"))

ggsave(file.path(ypd_stripe_dir, "tss_dinucleotide_frequencies_100ng.pdf"), plot = p, device = cairo_pdf, height = 7, width = 7)

# Plot distance of dominant TSS to annotated start codon
dominant <- dominant_tss(exp, threshold = 3, feature_type = "geneId", samples = "S288C_100ng_1")

p <- plot_dominant_tss(dominant) +
	theme(text = element_text(size = 12), axis.text = element_text(color="black"))

ggsave(file.path(ypd_stripe_dir, "dominant_tss.pdf"), plot = p, device = cairo_pdf, height = 3, width = 3)

# Plot hypothetical maximum 5'UTR length
max <- max_utr(exp, threshold = 3, feature_type = "geneId", samples = "S288C_100ng_1")

p <- plot_max_utr(max) +
	theme(text = element_text(size = 12), axis.text = element_text(color="black"))

ggsave(file.path(ypd_stripe_dir, "max_utr.pdf"), plot = p, device = cairo_pdf, height = 3, width = 3)

# Export normalized TSS bedGraphs

bedgraph_dir <- file.path("yeast_work", "bedgraphs")

if (!dir.exists(bedgraph_dir)){
  message("Creating directory 'yeast_work/bedgraphs'...")
  dir.create(file.path(bedgraph_dir))
} else {
  message("Directory 'yeast_work/bedgraphs' already exists...")
}

iwalk(exp@counts$TSSs$cpm, function(counts, sample) {
	plus <- counts[strand(counts) == "+"]
	minus <- counts[strand(counts) == "-"]

	export(plus, file.path(bedgraph_dir, paste(sample, "plus.bedgraph", sep = "_")), format = "bedgraph")
	export(minus, file.path(bedgraph_dir, paste(sample, "minus.bedgraph", sep = "_")), format = "bedgraph")
})

# Normalize TSR counts
exp <- count_normalization(exp, data_type = "tsr", threshold = 3, n_samples = 1, samples = stripe)

# Annotate TSRs
exp <- annotate_features(exp, annotation_file = annotation, data_type = "tsr", feature_type = "transcript")

# Determine TSR distribution relative to genomic features
tsr_distribution <- genomic_distribution(exp, data_type = "tsr", threshold = 3, samples = "S288C_100ng_1")

p <- plot_genomic_distribution(tsr_distribution) +
    ggplot2::theme(text = element_text(size = 12), legend.key.size = unit(0.6, "cm"), axis.text = element_text(color="black"))

ggsave(file.path(ypd_stripe_dir, "tsr_genomic_distribution.pdf"), plot = p, device = cairo_pdf, height = 2, width = 6)

tsr_distribution <- genomic_distribution(exp, data_type = "tsr", threshold = 3, quantiles = 5, 
                                         samples = stripe)

p <- plot_genomic_distribution(tsr_distribution, sample_order = stripe) +
    ggplot2::theme(text = element_text(size = 14), legend.key.size = unit(0.8, "cm"), axis.text = element_text(color="black"))

ggsave(file.path(ypd_stripe_dir, "tsr_genomic_distribution_quantiles.pdf"), plot = p, device = cairo_pdf, height = 12, width = 6)

# Plot number of promoter-proximal features with a TSR
features <- detect_features(exp, data_type = "tsr", feature_type = "transcript", samples = stripe)

p <- plot_detected_features(features, ncol = 3, width = 0.75) +
    ggplot2::theme(text = element_text(size = 14), legend.key.size = unit(0.8, "cm"), axis.text = element_text(color="black"))

ggsave(file.path(ypd_stripe_dir, "tsr_feature_plot.pdf"), plot = p, device = cairo_pdf, height = 4.5, width = 8)

# Generate TSR density plot
p <- plot_average(exp, data_type = "tsr", samples = "S288C_100ng_1", upstream = 1000, downstream = 1000) +
    ggplot2::theme(text = element_text(size = 13), axis.text = element_text(color="black"))

ggsave(file.path(ypd_stripe_dir, "tsr_density_plot.pdf"), plot = p, device = cairo_pdf, height = 2.5, width = 3.5)

# Export TSR beds
iwalk(exp@counts$TSRs$cpm, function(counts, sample) {
    export(counts, file.path(bedgraph_dir, paste(sample, "TSRs.bed", sep = "_")), format = "bed")
})

####################################
### STRIPE-seq vs. CAGE analysis ###
####################################

cage_dir <- file.path("yeast_work", "CAGE")

if (!dir.exists(cage_dir)){
  message("Creating directory 'yeast_work/CAGE'...")
  dir.create(cage_dir)
} else {
  message("Directory 'yeast_work/CAGE' already exists.")
}

# Normalize TSS counts
exp <- count_normalization(exp, data_type = "tss", threshold = 3, n_samples = 1, samples = all)

# Annotate TSSs relative to genomic features
exp <- annotate_features(exp, annotation_file = annotation, data_type = "tss", feature_type = "transcript", 
                         upstream = 250, downstream = 100)

# Explore TSS read thresholds for promoter fraction and plot
thresh <- explore_thresholds(exp, annotation_file = annotation, feature_type = "transcript", max_threshold = 25, 
                             upstream = 250, downstream = 100, samples = all)

p <- plot_threshold_exploration(thresh, ncol = 3, point_size = 2, sample_order = all) +
    ggplot2::geom_vline(xintercept = 3, lty = 2) + 
    ggplot2::theme(legend.key.size = unit(0.8, "cm"), text = element_text(size = 12), axis.text = element_text(color="black"))

ggsave(file.path(cage_dir, "tss_thresholds.pdf"), plot = p, device = cairo_pdf, height = 8, width = 11)

# Determine TSS distribution relative to genomic features
tss_distribution <- genomic_distribution(exp, data_type = "tss", threshold = 3, samples = all)

p <- plot_genomic_distribution(tss_distribution, sample_order = cage) +
    ggplot2::theme(text = element_text(size = 14), legend.key.size = unit(0.8, "cm"), axis.text = element_text(color="black"))

ggsave(file.path(cage_dir, "tss_genomic_distribution.pdf"), plot = p, device = cairo_pdf, height = 6, width = 6)

# Print all TSS dinucleotide frequencies
frequencies <- dinucleotide_frequencies(exp, genome_assembly = assembly, threshold = 3, samples = all)

p <- plot_dinucleotide_frequencies(frequencies, ncol = 3, sample_order = all) +
    ggplot2::theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"), axis.text = element_text(color="black"))

ggsave(file.path(cage_dir, "tss_dinucleotide_frequencies_all.pdf"), plot = p, device = cairo_pdf, height = 12, width = 9)

# Assess TSS dinucleotide frequencies (1 sample each method)
frequencies <- dinucleotide_frequencies(exp, genome_assembly = assembly, threshold = 3, samples = cage_slim)

p <- plot_dinucleotide_frequencies(frequencies, ncol = 2, sample_order = cage_slim) +
  ggplot2::theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"), axis.text = element_text(color="black"))

ggsave(file.path(cage_dir, "tss_dinucleotide_frequencies_slim.pdf"), plot = p, device = cairo_pdf, height = 5, width = 5)

# Export normalized TSS bedGraphs
iwalk(exp@counts$TSSs$cpm, function(counts, sample) {
	plus <- counts[strand(counts) == "+"]
	minus <- counts[strand(counts) == "-"]

	export(plus, file.path(bedgraph_dir, paste(sample, "plus.bedgraph", sep = "_")), format = "bedgraph")
	export(minus, file.path(bedgraph_dir, paste(sample, "minus.bedgraph", sep = "_")), format = "bedgraph")
})

# Normalize TSR counts
exp <- count_normalization(exp, data_type = "tsr", threshold = 3, n_samples = 1, samples = all)

# Annotate TSRs
exp <- annotate_features(exp, annotation_file = annotation, data_type = "tsr", feature_type = "transcript")

# Determine TSR distribution relative to genomic features
tsr_distribution <- genomic_distribution(exp, data_type = "tsr", threshold = 3, samples = cage)

p <- plot_genomic_distribution(tsr_distribution, sample_order = cage) +
    ggplot2::theme(text = element_text(size = 12), legend.key.size = unit(0.6, "cm"), axis.text = element_text(color="black"))

ggsave(file.path(cage_dir, "tsr_genomic_distribution.pdf"), plot = p, device = cairo_pdf, height = 6, width = 6)

# Plot number of promoter-proximal features with a TSR (1 sample each method)
features <- detect_features(exp, data_type = "tsr", feature_type = "transcript", samples = all)

p <- plot_detected_features(features, ncol = 3, width = 0.75) +
    ggplot2::theme(text = element_text(size = 10), legend.key.size = unit(0.8, "cm"), axis.text = element_text(color="black"))

ggsave(file.path(cage_dir, "tsr_feature_plot.pdf"), plot = p, device = cairo_pdf, height = 10, width = 10)

# Determine average promoter feature counts and plot
p <- features %>% 
    add_column(technology = c(rep("nanoCAGE25", 2), rep("nanoCAGE500", 2), 
                              rep("STRIPE100", 3), rep("STRIPE250", 3), rep("STRIPE50", 3), 
                              rep("SLIC100", 2))) %>%
    mutate(technology = factor(technology, levels=c("STRIPE50", "STRIPE100", "STRIPE250",
                                                    "SLIC100", "nanoCAGE500", "nanoCAGE25"))) %>%
    ggplot(aes(x = technology, y = promoter_proximal_features)) +
    theme_bw() +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), size = 0.4, geom = "errorbar", width = 0.25,
                 alpha = 0.5) +
    geom_jitter(aes(color = technology), position = position_jitter(height = 0, width = 0.2), size = 2) +
    scale_color_viridis_d(begin = 0, end = 0.7) +
    ylim(2000, 6000)

ggsave(file.path(cage_dir, "promoter_feature_jitter_plot.pdf"), plot = p, device = cairo_pdf, height = 4, width = 6)

# Export TSR beds
iwalk(exp@counts$TSRs$cpm, function(counts, sample) {
    export(counts, file.path(bedgraph_dir, paste(sample, "TSRs.bed", sep = "_")), format = "bed")
})

########################
### Diamide analysis ###
########################

diamide_dir <- file.path("yeast_work", "diamide")

if (!dir.exists(diamide_dir)){
  message("Creating directory 'yeast_work/diamide'...")
  dir.create(diamide_dir)
} else {
  message("Directory 'yeast_work/diamide' already exists.")
}

# Normalize TSS counts
exp <- count_normalization(exp, data_type = "tss", threshold = 3, n_samples = 1, samples = diamide)

# Generate a hierarchically clustered TSS heatmap with correlation values displayed
corr_matrix <- find_correlation(exp, data_type = "tss", correlation_metric = "pearson", samples = diamide)

cairo_pdf(file = file.path(diamide_dir, "tss_correlation_hierarchical.pdf"), width = 9, height = 9)
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
tss_distribution <- genomic_distribution(exp, data_type = "tss", threshold = 3, samples = diamide)

p <- plot_genomic_distribution(tss_distribution, sample_order = diamide) +
    ggplot2::theme(text = element_text(size = 14), legend.key.size = unit(0.8, "cm"), axis.text = element_text(color="black"))

ggsave(file.path(diamide_dir, "tss_genomic_distribution.pdf"), plot = p, device = cairo_pdf, height = 4.5, width = 6)

# Generate TSS sequence logos
seqs <- tss_sequences(exp, genome_assembly = assembly, threshold = 3, 
                      samples = c("S288C_100ng_1","S288C_diamide_100ng_1"))

p <- plot_sequence_logo(seqs, ncol = 1, font_size = 10)

ggsave(file.path(diamide_dir, "tss_seq_logo.pdf"), plot = p, device = cairo_pdf, height = 2.5, width = 4)

# Plot distance of dominant TSS to annotated start codon
dominant <- dominant_tss(exp, threshold = 3, feature_type = "geneId", samples = c("S288C_100ng_1","S288C_diamide_100ng_1"))

p <- plot_dominant_tss(dominant) +
	theme(text = element_text(size = 12), axis.text = element_text(color="black"))

ggsave(file.path(diamide_dir, "dominant_tss.pdf"), plot = p, device = cairo_pdf, height = 3, width = 3)

# Plot hypothetical maximum 5'UTR length
max <- max_utr(exp, threshold = 3, feature_type = "geneId", samples = c("S288C_100ng_1","S288C_diamide_100ng_1"))

p <- plot_max_utr(max) +
	theme(text = element_text(size = 12), axis.text = element_text(color="black"))

ggsave(file.path(diamide_dir, "max_utr.pdf"), plot = p, device = cairo_pdf, height = 3, width = 3)

# Export normalized TSS bedGraphs
iwalk(exp@counts$TSSs$cpm, function(counts, sample) {
	plus <- counts[strand(counts) == "+"]
	minus <- counts[strand(counts) == "-"]

	export(plus, file.path(bedgraph_dir, paste(sample, "plus.bedgraph", sep = "_")), format = "bedgraph")
	export(minus, file.path(bedgraph_dir, paste(sample, "minus.bedgraph", sep = "_")), format = "bedgraph")
})

# Normalize TSR counts
exp <- count_normalization(exp, data_type = "tsr", threshold = 3, n_samples = 1, samples = diamide)

# Generate a hierarchically clustered TSR heatmap with correlation values displayed
corr_matrix <- find_correlation(exp, data_type = "tsr", correlation_metric = "pearson", samples = diamide)

cairo_pdf(file = file.path(diamide_dir, "tsr_correlation_hierarchical.pdf"), width = 9, height = 9)
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
    ggplot2::theme(text = element_text(size = 12), legend.key.size = unit(0.6, "cm"), axis.text = element_text(color="black"))

ggsave(file.path(diamide_dir, "tsr_genomic_distribution.pdf"), plot = p, device = cairo_pdf, height = 4, width = 6)

# Export TSR beds
iwalk(exp@counts$TSRs$cpm, function(counts, sample) {
    export(counts, file.path(bedgraph_dir, paste(sample, "TSRs.bed", sep = "_")), format = "bed")
})

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
    filter(abs(log2FC) >= 1, FDR <= 0.05) %>%
    write.table(file.path(diamide_dir, "diff_tsrs.tsv"),
                sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Annotate dTSRs
annotated_diff_tsrs <- annotate_differential(diff_tsrs, annotation_file = annotation, 
                                                  feature_type = "transcript", 
                                                  upstream = 250, downstream = 100)

# Write annotated significant dTSRs to a table
annotated_diff_tsrs %>% 
    filter(annotation == "Promoter", abs(log2FC) >= 1, FDR <= 0.05) %>%
    write.table(file.path(diamide_dir, "promoter_annotated_diff_tsrs.tsv"),
                sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

promoter_diff_tsrs <- annotated_diff_tsrs %>% 
    filter(annotation == "Promoter", abs(log2FC) >= 1, FDR <= 0.05)

# Count dTSRs
diff_tsrs %>%
    mutate(change = case_when(
        log2FC >= 1 & FDR <= 0.05 ~ "increased",
        log2FC <= -1 & FDR <= 0.05 ~ "decreased",
        TRUE ~ "unchanged")) %>%
    dplyr::count(change)

# Count promoter dTSRs
annotated_diff_tsrs %>% 
    filter(annotation == "Promoter") %>%
    mutate(change = case_when(
        log2FC >= 1 & FDR <= 0.05 ~ "increased",
        log2FC <= -1 & FDR <= 0.05 ~ "decreased",
        TRUE ~ "unchanged")) %>%
    dplyr::count(change)

# Make a volcano plot of dTSRs
p <- plot_volcano(diff_tsrs, size = 0.25) + 
    ggplot2::theme(text = element_text(size = 12), axis.text = element_text(color="black"))

ggsave(file.path(diamide_dir, "diff_tsrs_volcano_plot.pdf"), plot = p, device = cairo_pdf, height = 2.5, width = 5)

# Perform GO analysis
enrichment_data <- export_for_enrichment(promoter_diff_tsrs)

go_enrichment <- clusterProfiler::compareCluster(
    geneId ~ change,
    data = enrichment_data,
    fun = "enrichGO",
    OrgDb = "org.Sc.sgd.db",
    pAdjustMethod = "fdr",
    keyType = "ENSEMBL",
    ont = "BP"
)

cairo_pdf(file = file.path(diamide_dir, "go_analysis.pdf"), height = 8, width = 10)
dotplot(go_enrichment, showCategory = 10) + scale_color_viridis()
dev.off()

##########################################
### RNA-seq analysis (YPD and diamide) ###
##########################################

rnaseq_dir <- file.path("yeast_work", "RNA_seq")

if (!dir.exists(rnaseq_dir)){
    message("Creating directory 'yeast_work/RNA_seq'...")
    dir.create(rnaseq_dir)
} else {
    message("Directory 'yeast_work/RNA_seq' already exisits.")
}

# Get feature counts for STRIPE-seq and RNA-seq
stripe_counts <- read.delim(file.path("yeast_data", "RNA_seq", "cleaned_S288C_feature_counts_STRIPEseq.tsv"),
                            sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rnaseq_counts <- read.delim(file.path("yeast_data", "RNA_seq", "cleaned_yeast_feature_counts.tsv"),
                             sep = "\t", header = TRUE, stringsAsFactors = FALSE)

stripe_counts <- column_to_rownames(stripe_counts, "Geneid") %>% 
    as.matrix
rnaseq_counts <- column_to_rownames(rnaseq_counts, "Geneid") %>% 
    as.matrix

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
corr_matrix <- find_correlation(exp, data_type = "features", correlation_metric = "spearman", 
                                samples = stripe_rnaseq_ypd)

cairo_pdf(file = file.path(rnaseq_dir, "tss_rnaseq_correlation_ypd_hierarchical.pdf"), width = 18, height = 18)
Heatmap(corr_matrix, col = viridis(256), heatmap_legend_param = list(title = "spearman"), 
        layer_fun = function(j, i, x, y, width, height, fill)
        {
            grid.text(sprintf("%.3f", pindex(corr_matrix, i, j)), x, y, gp = gpar(fontsize = 26, col = "white"))
        }
)
dev.off()

# Find correlation of diamide STRIPE-seq and RNA-seq feature counts
corr_matrix <- find_correlation(exp, data_type = "features", correlation_metric = "spearman", 
                                samples = stripe_rnaseq_diamide)

cairo_pdf(file = file.path(rnaseq_dir, "tss_rnaseq_correlation_diamide_hierarchical.pdf"), width = 18, height = 18)
Heatmap(corr_matrix, col = viridis(256), heatmap_legend_param = list(title = "spearman"), 
        layer_fun = function(j, i, x, y, width, height, fill)
        {
            grid.text(sprintf("%.3f", pindex(corr_matrix, i, j)), x, y, gp = gpar(fontsize = 26, col = "white"))
        }
)
dev.off()

# Perform YPD vs. diamide RNA-seq-like differential expression analysis for STRIPE-seq and RNA-seq

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

diff_features_stripe <- differential_expression(edger_model_stripe, data_type = "features", 
                                                compare_groups = c(1, 2)) %>%
    dplyr::rename("geneId" = "position")

# Filter differential STRIPE-seq features
diff_features_stripe_filtered <- diff_features_stripe %>% 
    filter(abs(log2FC) >= 1, FDR <= 0.05)

# Write differential STRIPE-seq features to a table
write.table(diff_features_stripe_filtered, file.path(rnaseq_dir, "diff_features_stripe_filtered.tsv"),
            sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Count differential STRIPE-seq features
diff_features_stripe %>%
    dplyr::count(log2FC >= 1 & FDR < 0.05)
diff_features_stripe %>%
    dplyr::count(log2FC <= -1 & FDR < 0.05)

# Add a column for identifying DE STRIPE-seq features in merged DE table
diff_features_stripe_filtered_id <- diff_features_stripe_filtered %>%
    mutate(change = case_when(
        log2FC >= 1 & FDR <= 0.05 ~ "increased",   
        log2FC <= -1 & FDR <= 0.05 ~ "decreased",
        TRUE ~ "unchanged"))

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

diff_features_rnaseq <- differential_expression(edger_model_rnaseq, data_type = "features", 
                                                compare_groups = c(1, 2)) %>%
    dplyr::rename("geneId" = "position")

# Filter differential RNA-seq features
diff_features_rnaseq_filtered <- diff_features_rnaseq %>% 
    filter(abs(log2FC) >= 1, FDR <= 0.05)

# Write differential RNA-seq features to a table
write.table(diff_features_rnaseq_filtered, file.path(rnaseq_dir, "diff_features_rnaseq_filtered.tsv"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Count differential RNA-seq features
diff_features_rnaseq %>%
    dplyr::count(log2FC >= 1 & FDR < 0.05)
diff_features_rnaseq %>%
    dplyr::count(log2FC <= -1 & FDR < 0.05)

# Add a column for identifying DE RNA-seq features in merged DE table
diff_features_rnaseq_filtered_id <- diff_features_rnaseq_filtered %>%
    mutate(change = case_when(
        log2FC >= 1 & FDR <= 0.05 ~ "increased",   
        log2FC <= -1 & FDR <= 0.05 ~ "decreased",
        TRUE ~ "unchanged"))

# Export RNA-seq results for GO enrichment, adding an identifying column for compareCluster
enrichment_data_rnaseq <- export_for_enrichment(diff_features_rnaseq, log2fc_cutoff = 1, fdr_cutoff = 0.05) %>%
    add_column(technology = "RNAseq")

# Combine edgeR results tables and count numbers of shared and distinct DE genes
merged_de <- full_join(diff_features_stripe_filtered_id, diff_features_rnaseq_filtered_id, 
                       by = "geneId", suffix = c("_STRIPE","_RNAseq")) %>%
    as_tibble(.name_repair = "unique") %>%
    replace_na(list(
    FDR_RNAseq = 1, FDR_STRIPE = 1,
    log2FC_RNAseq = 0, log2FC_STRIPE = 0,
    change_RNAseq = "unchanged", change_STRIPE = "unchanged"
    ))

merged_de_counts <- merged_de %>%
    dplyr::count(change_STRIPE, change_RNAseq)
    
# Count various combinations of expression change in STRIPE-seq and RNA-seq and plot Euler diagrams

# Upregulated
up_shared <- merged_de_counts %>%
    dplyr::filter(change_STRIPE == "increased" & change_RNAseq == "increased") %>%
    dplyr::select(n) %>%
    pull

up_stripe <- merged_de_counts %>%
    dplyr::filter(change_STRIPE == "increased" & change_RNAseq != "increased") %>%
    dplyr::select(n) %>%
    sum

up_rnaseq <- merged_de_counts %>%
    dplyr::filter(change_STRIPE != "increased" & change_RNAseq == "increased") %>%
    dplyr::select(n) %>%
    sum

upregulated_fit <- euler(c("A" = up_stripe, "B" = up_rnaseq, "A&B" = up_shared))

euler_fill = viridis(4)

cairo_pdf(file = file.path(rnaseq_dir, "upregulated_euler.pdf"), width = 7, height = 3)
plot(upregulated_fit, fill = euler_fill, alpha=c(0.5,0.5), lab=c("STRIPE-seq","RNA-seq"))
dev.off()

# Count total number of upregulated genes
merged_de_counts %>% 
    filter(change_STRIPE == "increased" | change_RNAseq == "increased") %>% 
    dplyr::select(n) %>% 
    sum

# Downregulated
down_shared <- merged_de_counts %>%
    dplyr::filter(change_STRIPE == "decreased" & change_RNAseq == "decreased") %>%
    dplyr::select(n) %>%
    pull

down_stripe <- merged_de_counts %>%
    dplyr::filter(change_STRIPE == "decreased" & change_RNAseq != "decreased") %>%
    dplyr::select(n) %>%
    sum

down_rnaseq <- merged_de_counts %>%
    dplyr::filter(change_STRIPE != "decreased" & change_RNAseq == "decreased") %>%
    dplyr::select(n) %>%
    sum

downregulated_fit <- euler(c("A" = down_stripe, "B" = down_rnaseq, "A&B" = down_shared))

cairo_pdf(file = file.path(rnaseq_dir, "downregulated_euler.pdf"), width = 7, height = 3)
plot(downregulated_fit, fill = euler_fill, alpha=c(0.5,0.5), lab=c("STRIPE-seq","RNA-seq"))
dev.off()

# Count total number of downregulated genes
merged_de_counts %>% 
    filter(change_STRIPE == "decreased" | change_RNAseq == "decreased") %>% 
    dplyr::select(n) %>% 
    sum

# Determine cumulative averages for STRIPE-seq vs. RNA-seq feature counts
merged_de <- merged_de %>%
    dplyr::arrange(desc(abs(log2FC_RNAseq))) %>%
    mutate(
        match = ifelse(change_RNAseq == change_STRIPE, 1, 0),
        cumfrac = cummean(match)
    )

plot_data <- merged_de %>%
    dplyr::select(geneId, log2FC_RNAseq, cumfrac) %>%
    rowid_to_column %>%
    filter(log2FC_RNAseq > 0) %>%
    gather(key = "metric", value = "value", -rowid, -geneId)

p <- ggplot(plot_data, aes(x = rowid, y = abs(value))) +
    theme_bw() +
    geom_line() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 12)
    ) +
    facet_grid(metric ~ ., scales = "free")

ggsave(file.path(rnaseq_dir, "cumulative_plot.pdf"), plot = p, device = cairo_pdf, height = 3, width = 5)

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

cairo_pdf(file = file.path(rnaseq_dir, "go_analysis_combined.pdf"), height = 10, width = 8)
dotplot(go_enrichment_combined, showCategory = 10) + scale_color_viridis()
dev.off()