library(tsrexplorer)
library(tidyverse)
library(GenomicRanges)
library(viridis)

# Pull latest version of tsrexplorer
# devtools::install_github("rpolicastro/tsrexplorer", ref = "normalization", force = TRUE)

# Read in TSSs
S288C_TSSs <- map(list.files("../yeast_data/YPD_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                      as.data.frame %>%
                      makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                             start.field = "TSS", end.field = "TSS"))

names(S288C_TSSs) <- c("S288C_50ng_1","S288C_50ng_2","S288C_50ng_3",
                       "S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
                       "S288C_250ng_1","S288C_250ng_2","S288C_250ng_3")
                                                 
# Read in TSRs
S288C_TSRs <- map(list.files("../yeast_data/YPD_TSRs/", full.names = TRUE), ~ read.delim(.x) %>%
                      as.data.frame %>%
                      makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                               start.field = "start", end.field = "end"))

names(S288C_TSRs) <- c("S288C_50ng_1","S288C_50ng_2","S288C_50ng_3",
                       "S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
                       "S288C_250ng_1","S288C_250ng_2","S288C_250ng_3")

# Read in annotation files (included with package)
annotation <- system.file("extdata", "yeast_annotation.gtf", package="tsrexplorer")
assembly <- system.file("extdata", "yeast_assembly.fasta", package="tsrexplorer")

# Read in all TSSs and TSRs
exp <- tsr_explorer(S288C_TSSs, S288C_TSRs)

####################
### TSS analysis ###
####################

# Normalize counts
exp <- count_normalization(exp, data_type = "tss", n_samples = 1, threshold = 3)

# Generate a combination correlation plot
p <- plot_correlation(exp, data_type = "tss") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tss_correlation.png", plot = p, device = "png", type = "cairo", height = 3, width = 3)

# Generate a hierarchically clustered heatmap
p <- plot_correlation(exp, data_type = "tss", correlation_plot = "hierarchical", col = viridis(256))

pdf("tss_correlation_hierarchical.pdf", height = 3.5, width = 4)
p
dev.off()

# Annotate TSSs relative to genomic features
exp <- annotate_features(exp, annotation_file = annotation, data_type = "tss", feature_type = "transcript", upstream = 250, downstream = 100)

# Determine TSS distribution relative to genomic features
tss_distribution <- genomic_distribution(exp, data_type = "tss", threshold = 3, samples = "S288C_100ng_1")

p <- plot_genomic_distribution(tss_distribution) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tss_genomic_distribution.pdf", plot = p, device = cairo_pdf, height = 1.5, width = 4)

genomic_dist <- genomic_distribution(exp, data_type = "tss", threshold = 3, quantiles = 5, samples = "S288C_100ng_1")

p <- plot_genomic_distribution(genomic_dist) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tss_genomic_distribution_quantiles.pdf", plot = p, device = cairo_pdf, height = 4, width = 4)

# Plot number of promoter-proximal features with a TSS
features <- detect_features(exp, data_type = "tss", feature_type = "transcript", threshold = 3, samples = "S288C_100ng_1")

p <- plot_detected_features(features, ncol = 3) +
    ggplot2::theme(text = element_text(size = 5))

ggsave("tss_feature_plot.pdf", plot = p, device = cairo_pdf, height = 2, width = 3)

# Generate TSS density plots
p <- plot_average(exp, data_type = "tss", threshold = 3, ncol = 3, samples = "S288C_100ng_1") +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tss_average_plot.pdf", plot = p, cairo_pdf, height = 4, width = 4)

# Generate TSS sequence logos
seqs <- tss_sequences(exp, genome_assembly = assembly, threshold = 3, samples = "S288C_100ng_1")

p <- plot_sequence_logo(seqs, ncol = 3) +
    ggplot2::theme(text = element_text(size = 5))

ggsave("tss_seq_logo.pdf", plot = p, device = cairo_pdf, height = 1, width = 2)

# Generate TSS color plot
seqs <- tss_sequences(exp, genome_assembly = assembly, threshold = 3, samples = "S288C_100ng_1")

p <- plot_sequence_colormap(seqs, ncol = 3) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tss_seq_colormap.png", plot = p, device = "png", type = "cairo", height = 2.5, width = 2)

# Assess TSS dinucleotide frequencies
frequencies <- dinucleotide_frequencies(exp, genome_assembly = assembly, threshold = 3, samples = "S288C_100ng_1")

p <- plot_dinucleotide_frequencies(frequencies, ncol = 3) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tss_dinucleotide_frequencies.pdf", plot = p, device = cairo_pdf, height = 1.7, width = 2.5)

# Plot distance of dominant TSS to annotated start codon
dominant <- dominant_tss(exp, threshold = 3, feature_type = "geneId", samples = "S288C_100ng_1")

p <- plot_dominant_tss(dominant, upstream = 500, downstream = 500)

ggsave("dominant_tss.pdf", plot = p, device = cairo_pdf, height = 4, width = 4)

# Plot hypothetical maximum 5'UTR length
max <- max_utr(exp, threshold = 3, feature_type = "geneId", samples = "S288C_100ng_1")

p <- plot_max_utr(max)

ggsave("max_utr.pdf", plot = p, device = cairo_pdf, height = 4, width = 4)
    
####################
### TSR analysis ###
####################

# Normalize TSR counts
exp <- count_normalization(exp, data_type = "tsr", threshold = 3)

## Generate a combination TSR correlation plot
p <- plot_correlation(exp, data_type = "tsr") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_correlation.png", plot = p, device = "png", type = "cairo", height = 9, width = 9)

# Generate a hierarchically clustered heatmap
p <- plot_correlation(exp, data_type = "tsr", correlation_plot = "hierarchical")

pdf("tsr_correlation_hierarchical.pdf", height = 3.5, width = 4)
p
dev.off()

# Plot selected TSR metrics
p <- plot_tsr_metric(exp, tsr_metrics = c("nTAGs", "nTSSs"), log2_transform = TRUE, ncol = 2) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_metrics.pdf", plot = p, device = cairo_pdf, width = 4, height = 2)

# Annotate TSRs
exp <- annotate_features(exp, annotation_file = annotation, data_type = "tsr", feature_type = "transcript")

# Determine TSR distribution relative to genomic features
tsr_distribution <- genomic_distribution(exp, data_type = "tsr", samples = "S288C_100ng_1")

p <- plot_genomic_distribution(tsr_distribution) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_genomic_distribution.png", plot = p, device = "png", type = "cairo", height = 1.5, width = 4)

# Plot number of promoter-proximal features with a TSR
features <- detect_features(exp, data_type = "tsr", feature_type = "transcript")

p <- plot_detected_features(features) +
    ggplot2::theme(text = element_text(size = 5))

ggsave("tsr_feature_plot.png", plot = p, device = "png", type = "cairo", height = 2, width = 4)

# Generate TSR average plot
p <- plot_average(exp, data_type = "tsr", samples = "S288C_100ng_1") +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_average_plot.png", plot = p, device = "png", type = "cairo", height = 1.5, width = 4)