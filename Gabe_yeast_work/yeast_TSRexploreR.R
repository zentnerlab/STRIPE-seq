library(tsrexplorer)
library(tidyverse)
library(GenomicRanges)
library(viridis)
library(rtracklayer)

# Pull latest version of tsrexplorer
# devtools::install_github("rpolicastro/tsrexplorer", ref = "dev", force = TRUE)

setwd("/Users/gzentner/Desktop/tsrexplorer/yeast/STRIPE-seq/Gabe_yeast_work/")

####################
### Read in TSSs ###
####################

# STRIPE-seq YPD TSSs
YPD_TSSs <- map(list.files("../yeast_data/YPD_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                      as.data.frame %>%
                      makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                             start.field = "TSS", end.field = "TSS")) %>%
    set_names(c("S288C_50ng_1","S288C_50ng_2","S288C_50ng_3",
                "S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
                "S288C_250ng_1","S288C_250ng_2","S288C_250ng_3"))

# SLIC-CAGE and nanoCAGE TSSs
CAGE_TSSs <- map(list.files("../yeast_data/CAGE_TSSs", full.names = TRUE), ~ read.delim(.x) %>%
                     as.data.frame %>%
                     makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                              start.field = "TSS", end.field = "TSS")) %>%
    set_names(c("SLIC_CAGE_100ng_1","SLIC_CAGE_100ng_2",
                "nanoCAGE_500ng_1","nanoCAGE_500ng_2",
                "nanoCAGE_25ng_1","nanoCAGE_25ng_2"))

# STRIPE-seq diamide TSSs
diamide_TSSs <- map(list.files("../yeast_data/diamide_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                        as.data.frame %>%
                        makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                                 start.field = "TSS", end.field = "TSS")) %>%
    set_names(c("S288C_diamide_100ng_1","S288C_diamide_100ng_2","S288C_diamide_100ng_3"))

# Generate list of all TSS objects
full_TSSs_set <- c(YPD_TSSs,CAGE_TSSs,diamide_TSSs)

####################
### Read in TSRs ###
####################

# STRIPE-seq YPD TSRs
YPD_TSRs <- map(list.files("../yeast_data/YPD_TSRs/", full.names = TRUE), ~ read.delim(.x) %>%
                      as.data.frame %>%
                      makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                               start.field = "start", end.field = "end")) %>%
    set_names(c("S288C_50ng_1","S288C_50ng_2","S288C_50ng_3",
                     "S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
                     "S288C_250ng_1","S288C_250ng_2","S288C_250ng_3"))

# SLIC-CAGE and nanoCAGE TSRs
CAGE_TSRs <- map(list.files("../yeast_data/CAGE_TSRs/", full.names = TRUE), ~ read.delim(.x) %>%
                     as.data.frame %>%
                     makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                              start.field = "start", end.field = "end")) %>%
    set_names(c("SLIC_CAGE_100ng_1","SLIC_CAGE_100ng_2",
                "nanoCAGE_500ng_1","nanoCAGE_500ng_2",
                "nanoCAGE_25ng_1","nanoCAGE_25ng_2"))

# STRIPE-seq diamide TSRs
diamide_TSRs <- map(list.files("../yeast_data/diamide_TSRs/", full.names = TRUE), ~ read.delim(.x) %>%
                        as.data.frame %>%
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

# YPD 100 ng and YPD + diamide STRIPE-seq samples
diamide <- c("S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
             "S288C_diamide_100ng_1","S288C_diamide_100ng_2","S288C_diamide_100ng_3")

# Read TSSs and TSRs into TSRexploreR
exp <- tsr_explorer(full_TSSs_set,full_TSR_set)

##############################################
### YPD input variation replicate analysis ###
##############################################

setwd("/Users/gzentner/Desktop/tsrexplorer/yeast/STRIPE-seq/Gabe_yeast_work/YPD_STRIPE/")

# Normalize TSS counts
exp <- count_normalization(exp, data_type = "tss", n_samples = 1, threshold = 3, samples = stripe)

# Generate a hierarchically clustered heatmap
p <- plot_correlation(exp, data_type = "tss", correlation_plot = "hierarchical", col = viridis(256), correlation_metric = "pearson")

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

p <- plot_sequence_logo(seqs, ncol = 1) +
    ggplot2::theme(text = element_text(size = 5))

ggsave("tss_seq_logo.pdf", plot = p, device = cairo_pdf, height = 3, width = 2)

seqs <- tss_sequences(exp, genome_assembly = assembly, threshold = 3, quantiles = 5, samples = "S288C_100ng_1")

p <- plot_sequence_logo(seqs, ncol = 1) +
    ggplot2::theme(text = element_text(size = 5))

ggsave("tss_seq_logo_quantiles.pdf", plot = p, device = cairo_pdf, height = 5, width = 5)

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

# Export normalized TSS bedGraphs
export.bedGraph(exp@counts$TSSs$cpm$S288C_100ng_1[strand(exp@counts$TSSs$cpm$S288C_100ng_1) == "+"], "S288C_100ng_1_+.bedgraph")
export.bedGraph(exp@counts$TSSs$cpm$S288C_100ng_1[strand(exp@counts$TSSs$cpm$S288C_100ng_1) == "-"], "S288C_100ng_1_-.bedgraph")

# Normalize TSR counts
exp <- count_normalization(exp, data_type = "tsr", threshold = 3, n_samples = 1, samples = stripe)

# Generate a hierarchically clustered heatmap
p <- plot_correlation(exp, data_type = "tsr", correlation_plot = "hierarchical", col = viridis(256), correlation_metric = "pearson")

pdf("tsr_correlation_hierarchical.pdf", height = 7, width = 7)
p
dev.off()

# Annotate TSRs
exp <- annotate_features(exp, annotation_file = annotation, data_type = "tsr", feature_type = "transcript")

# Determine TSR distribution relative to genomic features
tsr_distribution <- genomic_distribution(exp, data_type = "tsr", threshold = 3, samples = "S288C_100ng_1")

p <- plot_genomic_distribution(tsr_distribution) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_genomic_distribution.png", plot = p, device = "png", type = "cairo", height = 1.5, width = 4)

tsr_distribution <- genomic_distribution(exp, data_type = "tsr", threshold = 3, quantiles = 5, samples = "S288C_100ng_1")

p <- plot_genomic_distribution(tsr_distribution) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_genomic_distribution_quantiles.png", plot = p, device = "png", type = "cairo", height = 1.5, width = 4)

# Plot number of promoter-proximal features with a TSR
features <- detect_features(exp, data_type = "tsr", feature_type = "transcript", samples = "S288C_100ng_1")

p <- plot_detected_features(features) +
    ggplot2::theme(text = element_text(size = 5))

ggsave("tsr_feature_plot.png", plot = p, device = "png", type = "cairo", height = 2, width = 4)

# Plot selected TSR metrics
p <- plot_tsr_metric(exp, tsr_metrics = "nTSSs", log2_transform = TRUE, ncol = 1, plot_type = "boxjitter", samples = "S288C_100ng_1", 
                     size = 0.5) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_metrics.pdf", plot = p, device = cairo_pdf, width = 7, height = 7)

# Generate TSR average plot
p <- plot_average(exp, data_type = "tsr", samples = "S288C_100ng_1") +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_average_plot.png", plot = p, device = "png", type = "cairo", height = 3, width = 3)

####################################
### STRIPE-seq vs. CAGE analysis ###
####################################

setwd("/Users/gzentner/Desktop/tsrexplorer/yeast/STRIPE-seq/Gabe_yeast_work/CAGE/")

# Normalize TSS counts
exp <- count_normalization(exp, data_type = "tss", n_samples = 1, threshold = 3, samples = cage)

# Export normalized TSS bedGraphs
export.bedGraph(exp@counts$TSSs$cpm$SLIC_CAGE_100ng_1[strand(exp@counts$TSSs$cpm$SLIC_CAGE_100ng_1) == "+"], "SLIC_CAGE_100ng_1_+.bedgraph")
export.bedGraph(exp@counts$TSSs$cpm$SLIC_CAGE_100ng_1[strand(exp@counts$TSSs$cpm$SLIC_CAGE_100ng_1) == "-"], "SLIC_CAGE_100ng_1_-.bedgraph")

export.bedGraph(exp@counts$TSSs$cpm$nanoCAGE_500ng_1[strand(exp@counts$TSSs$cpm$nanoCAGE_500ng_1) == "+"], "nanoCAGE_500ng_1_+.bedgraph")
export.bedGraph(exp@counts$TSSs$cpm$nanoCAGE_500ng_1[strand(exp@counts$TSSs$cpm$nanoCAGE_500ng_1) == "-"], "nanoCAGE_500ng_1_-.bedgraph")

export.bedGraph(exp@counts$TSSs$cpm$nanoCAGE_25ng_1[strand(exp@counts$TSSs$cpm$nanoCAGE_25ng_1) == "+"], "nanoCAGE_25ng_1_+.bedgraph")
export.bedGraph(exp@counts$TSSs$cpm$nanoCAGE_25ng_1[strand(exp@counts$TSSs$cpm$nanoCAGE_25ng_1) == "-"], "nanoCAGE_25ng_1_-.bedgraph")

# Normalize TSR counts
exp <- count_normalization(exp, data_type = "tsr", threshold = 3, n_samples = 1, samples = cage)

# Generate a hierarchically clustered heatmap
p <- plot_correlation(exp, data_type = "tsr", correlation_plot = "hierarchical", col = viridis(256), correlation_metric = "spearman")

pdf("tsr_correlation_hierarchical.pdf", height = 7, width = 7)
p
dev.off()

# Annotate TSRs
exp <- annotate_features(exp, annotation_file = annotation, data_type = "tsr", feature_type = "transcript")

# Determine TSR distribution relative to genomic features
tsr_distribution <- genomic_distribution(exp, data_type = "tsr", threshold = 3, samples = c("S288C_100ng_1","SLIC_CAGE_100ng_1",
                                                                                            "nanoCAGE_500ng_1","nanoCAGE_25ng_1"))

p <- plot_genomic_distribution(tsr_distribution) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_genomic_distribution.png", plot = p, device = "png", type = "cairo", height = 1.5, width = 4)

########################
### Diamide analysis ###
########################

setwd("/Users/gzentner/Desktop/tsrexplorer/yeast/STRIPE-seq/Gabe_yeast_work/diamide/")

# Normalize TSS counts
exp <- count_normalization(exp, data_type = "tss", n_samples = 1, threshold = 3, samples = diamide)

# Generate a hierarchically clustered heatmap
p <- plot_correlation(exp, data_type = "tss", correlation_plot = "hierarchical", col = viridis(256), correlation_metric = "pearson")

pdf("tss_correlation_hierarchical.pdf", height = 3.5, width = 4)
p
dev.off()

# Annotate TSSs relative to genomic features
exp <- annotate_features(exp, annotation_file = annotation, data_type = "tss", feature_type = "transcript", upstream = 250, downstream = 100)

# Determine TSS distribution relative to genomic features
tss_distribution <- genomic_distribution(exp, data_type = "tss", threshold = 3, samples = c("S288C_100ng_1","S288C_diamide_100ng_1"))

p <- plot_genomic_distribution(tss_distribution) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tss_genomic_distribution.pdf", plot = p, device = cairo_pdf, height = 1.5, width = 4)

# Plot number of promoter-proximal features with a TSS
features <- detect_features(exp, data_type = "tss", feature_type = "transcript", threshold = 3, samples = c("S288C_100ng_1","S288C_diamide_100ng_1"))

p <- plot_detected_features(features, ncol = 3) +
    ggplot2::theme(text = element_text(size = 5))

ggsave("tss_feature_plot.pdf", plot = p, device = cairo_pdf, height = 2, width = 3)

# Generate TSS sequence logos
seqs <- tss_sequences(exp, genome_assembly = assembly, threshold = 3, samples = c("S288C_100ng_1","S288C_diamide_100ng_1"))

p <- plot_sequence_logo(seqs, ncol = 1) +
    ggplot2::theme(text = element_text(size = 5))

ggsave("tss_seq_logo.pdf", plot = p, device = cairo_pdf, height = 3, width = 2)

# Plot distance of dominant TSS to annotated start codon
dominant <- dominant_tss(exp, threshold = 3, feature_type = "geneId", samples = c("S288C_100ng_1","S288C_diamide_100ng_1"))

p <- plot_dominant_tss(dominant, upstream = 500, downstream = 500)

ggsave("dominant_tss.pdf", plot = p, device = cairo_pdf, height = 4, width = 4)

# Plot hypothetical maximum 5'UTR length
max <- max_utr(exp, threshold = 3, feature_type = "geneId", samples = c("S288C_100ng_1","S288C_diamide_100ng_1"))

p <- plot_max_utr(max)

ggsave("max_utr.pdf", plot = p, device = cairo_pdf, height = 4, width = 4)

# Export normalized TSS bedGraphs
export.bedGraph(exp@counts$TSSs$cpm$S288C_diamide_100ng_1[strand(exp@counts$TSSs$cpm$S288C_diamide_100ng_1) == "+"], "S288C_diamide_100ng_1_+.bedgraph")
export.bedGraph(exp@counts$TSSs$cpm$S288C_diamide_100ng_1[strand(exp@counts$TSSs$cpm$S288C_diamide_100ng_1) == "-"], "S288C_diamide_100ng_1_-.bedgraph")

# Normalize TSR counts
exp <- count_normalization(exp, data_type = "tsr", threshold = 3, n_samples = 1, samples = diamide)

# Generate a hierarchically clustered heatmap
p <- plot_correlation(exp, data_type = "tsr", correlation_plot = "hierarchical", col = viridis(256), correlation_metric = "pearson")

pdf("tsr_correlation_hierarchical.pdf", height = 7, width = 7)
p
dev.off()

# Annotate TSRs
exp <- annotate_features(exp, annotation_file = annotation, data_type = "tsr", feature_type = "transcript")

# Determine TSR distribution relative to genomic features
tsr_distribution <- genomic_distribution(exp, data_type = "tsr", threshold = 3, samples = c("S288C_100ng_1","S288C_diamide_100ng_1"))

p <- plot_genomic_distribution(tsr_distribution) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_genomic_distribution.png", plot = p, device = "png", type = "cairo", height = 1.5, width = 4)

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
annotated_diff_tsrs <- annotate_differential_tsrs(diff_tsrs, annotation_file = annotation, feature_type = "transcript", 
                                                  upstream = 250, downstream = 100)

# Write annotated significant dTSRs to a table
annotated_diff_tsrs %>% 
    filter(annotation == "Promoter" & log2FC <= -1 & FDR < 0.05 | annotation == "Promoter" & log2FC >= 1 & FDR < 0.05) %>%
    write.table(., "promoter_annotated_diff_tsrs.tsv", sep="\t", col.names=T, row.names=F, quote=F)

# Make a volcano plot of dTSRs
p <- plot_volcano(diff_tsrs, size = 0.1)

ggsave("diff_tsrs_volcano_plot.png", plot = p, device = "png", type = "cairo", height = 2.5, width = 4)

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