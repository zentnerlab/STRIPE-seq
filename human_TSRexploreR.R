library(tsrexplorer)
library(tidyverse)
library(GenomicRanges)
library(Biostrings)
library(viridis)
library(rtracklayer)

# Pull latest version of tsrexplorer
# devtools::install_github("rpolicastro/tsrexplorer", ref = "dev", force = TRUE)

# Set base directory for analysis
baseDir <- "/Users/gzentner/Desktop/tsrexplorer/STRIPE-seq/"
setwd(baseDir)

if (!dir.exists(file.path(baseDir, "human_work/"))){
  print("Creating directory 'human_work'...")
  dir.create(file.path(baseDir, "human_work/"))
} else {
  print("Directory 'human_work' already exists...")
}

####################
### Read in TSSs ###
####################

# STRIPE-seq TSSs
STRIPE_TSSs <- map(list.files("human_data/STRIPE_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                      as.data.frame %>%
                      makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                             start.field = "TSS", end.field = "TSS")) %>%
    set_names(c("K562_100ng_1","K562_100ng_2","K562_100ng_3"))

# nanoCAGE-XL, CAGE, and RAMPAGE TSSs
CAGE_TSSs <- map(list.files("human_data/CAGE_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                     as.data.frame %>%
                     makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                              start.field = "TSS", end.field = "TSS")) %>%
    set_names(c("nanoCAGE_XL_7.5ug_1",
                "CAGE_10ug_1","CAGE_10ug_2",
                "RAMPAGE_5ug_1","RAMPAGE_5ug_2"))

# Generate list of all TSS sets
full_TSS_set <- c(STRIPE_TSSs,CAGE_TSSs)

####################
### Read in TSRs ###
####################

# STRIPE-seq TSRs
STRIPE_TSRs <- map(list.files("human_data/STRIPE_TSRs/", full.names = TRUE), ~ read.delim(.x) %>%
                      as.data.frame %>%
                      makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                               start.field = "start", end.field = "end")) %>%
  set_names(c("K562_100ng_1","K562_100ng_2","K562_100ng_3"))

# nanoCAGE-XL, CAGE, and RAMPAGE TSRs
CAGE_TSRs <- map(list.files("human_data/CAGE_TSRs/", full.names = TRUE), ~ read.delim(.x) %>%
                     as.data.frame %>%
                     makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                              start.field = "start", end.field = "end")) %>%
  set_names(c("nanoCAGE_XL_7.5ug_1",
              "CAGE_10ug_1","CAGE_10ug_2",
              "RAMPAGE_5ug_1","RAMPAGE_5ug_2"))

# Generate list of all TSR sets
full_TSR_set <- c(STRIPE_TSRs,CAGE_TSRs)

# Assign sample sets of interest to variables 

# STRIPE-seq samples
stripe <- c("K562_100ng_1","K562_100ng_2","K562_100ng_3")

# nanoCAGE-XL, CAGE, and RAMPAGE samples
cage <- c("nanoCAGE_XL_7.5ug_1",
          "CAGE_10ug_1","CAGE_10ug_2",
          "RAMPAGE_5ug_1","RAMPAGE_5ug_2")

# All samples
all <- c("K562_100ng_1","K562_100ng_2","K562_100ng_3",
         "nanoCAGE_XL_7.5ug_1",
         "CAGE_10ug_1","CAGE_10ug_2",
         "RAMPAGE_5ug_1","RAMPAGE_5ug_2")

# Read TSSs and TSRs into TSRexploreR
exp <- tsr_explorer(full_TSS_set,full_TSR_set)

###############################
### K562 replicate analysis ###
###############################

if (!dir.exists(file.path(baseDir, "human_work/STRIPE/"))){
  print("Creating directory 'human_work/STRIPE' and changing working directory...")
  dir.create(file.path(baseDir, "human_work/STRIPE/"))
  setwd(file.path(baseDir, "human_work/STRIPE/"))
} else {
  print("Directory 'human_work/STRIPE' already exists, changing working directory...")
  setwd(file.path(baseDir, "human_work/STRIPE/"))
}

# Normalize TSS counts
exp <- count_normalization(exp, data_type = "tss", threshold = 3, n_samples = 1, samples = stripe)

# Generate a combinbed TSS correlation plot
p <- plot_correlation(exp, data_type = "tss") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tss_correlation.png", plot = p, device = "png", type = "cairo", height = 12, width = 12)

# Generate a hierarchically clustered TSS heatmap
p <- plot_correlation(exp, data_type = "tss", correlation_plot = "hierarchical", col = viridis(256), 
                      correlation_metric = "pearson")

pdf("tss_correlation_hierarchical.pdf", height = 4, width = 4)
p
dev.off()

# Annotate TSSs relative to genomic features
exp <- annotate_features(exp, annotation_file = file.path(baseDir, "human_data/Homo_sapiens.GRCh38.98.chr.gtf"),
                         data_type = "tss", feature_type = "transcript", upstream = 250, downstream = 100)

# Determine TSS distribution relative to genomic features
tss_distribution <- genomic_distribution(exp, data_type = "tss", threshold = 3, samples = "K562_100ng_1")

p <- plot_genomic_distribution(tss_distribution) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tss_genomic_distribution.pdf", plot = p, device = cairo_pdf, height = 1.5, width = 4)

genomic_dist <- genomic_distribution(exp, data_type = "tss", threshold = 3, quantiles = 5, 
                                     samples = "K562_100ng_1")

p <- plot_genomic_distribution(genomic_dist) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tss_genomic_distribution_quantiles.pdf", plot = p, device = cairo_pdf, height = 4, width = 4)

# Plot number of promoter-proximal features with a TSS
features <- detect_features(exp, data_type = "tss", feature_type = "transcript", threshold = 3, 
                            samples = "K562_100ng_1")

p <- plot_detected_features(features, ncol = 3) +
    ggplot2::theme(text = element_text(size = 5))

ggsave("tss_feature_plot.pdf", plot = p, device = cairo_pdf, height = 2, width = 3)

# Generate TSS density plots
p <- plot_average(exp, data_type = "tss", threshold = 3, ncol = 3, samples = "K562_100ng_1") +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tss_average_plot.pdf", plot = p, cairo_pdf, height = 4, width = 4)

# Generate TSS sequence logos
seqs <- tss_sequences(exp, genome_assembly = file.path(baseDir, "human_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
                      threshold = 3, samples = "K562_100ng_1")

p <- plot_sequence_logo(seqs, ncol = 1) +
    ggplot2::theme(text = element_text(size = 5))

ggsave("tss_seq_logo.pdf", plot = p, device = cairo_pdf, height = 1, width = 2)

seqs <- tss_sequences(exp, genome_assembly = file.path(baseDir, "human_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
                      threshold = 3, quantiles = 5, samples = "K562_100ng_1")

p <- plot_sequence_logo(seqs, ncol = 1) +
    ggplot2::theme(text = element_text(size = 5))

ggsave("tss_seq_logo_quantiles.pdf", plot = p, device = cairo_pdf, height = 5, width = 5)

# Generate TSS color plot
seqs <- tss_sequences(exp, genome_assembly = file.path(baseDir, "human_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
                      threshold = 3, samples = "K562_100ng_1")

p <- plot_sequence_colormap(seqs, ncol = 3) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tss_seq_colormap.png", plot = p, device = "png", type = "cairo", height = 2.5, width = 2)

# Assess TSS dinucleotide frequencies
frequencies <- dinucleotide_frequencies(exp, genome_assembly = file.path(baseDir, "human_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
                                        threshold = 3, samples = "K562_100ng_1")

p <- plot_dinucleotide_frequencies(frequencies, ncol = 3) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tss_dinucleotide_frequencies.pdf", plot = p, device = cairo_pdf, height = 1.7, width = 2.5)

# Plot distance of dominant TSS to annotated start codon
dominant <- dominant_tss(exp, threshold = 3, feature_type = "geneId", samples = "K562_100ng_1")

p <- plot_dominant_tss(dominant, upstream = 500, downstream = 500)

ggsave("dominant_tss.pdf", plot = p, device = cairo_pdf, height = 4, width = 4)

# Plot hypothetical maximum 5'UTR length
max <- max_utr(exp, threshold = 3, feature_type = "geneId", samples = "K562_100ng_1")

p <- plot_max_utr(max)

ggsave("max_utr.pdf", plot = p, device = cairo_pdf, height = 4, width = 4)

# Export normalized TSS bedGraphs

if (!dir.exists(file.path(baseDir, "human_work/bedgraphs/"))){
  print("Creating directory 'human_work/bedgraphs'...")
  dir.create(file.path(baseDir, "human_work/bedgraphs/"))
} else {
  print("Directory 'human_work/bedgraphs' already exists...")
}

export.bedGraph(exp@counts$TSSs$cpm$K562_100ng_1[strand(exp@counts$TSSs$cpm$K562_100ng_1) == "+"], 
                file.path(baseDir, "human_work/bedgraphs/K562_100ng_1_+.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$K562_100ng_1[strand(exp@counts$TSSs$cpm$K562_100ng_1) == "-"], 
                file.path(baseDir, "human_work/bedgraphs/K562_100ng_1_-.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$K562_100ng_2[strand(exp@counts$TSSs$cpm$K562_100ng_2) == "+"], 
                file.path(baseDir, "human_work/bedgraphs/K562_100ng_2_+.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$K562_100ng_2[strand(exp@counts$TSSs$cpm$K562_100ng_2) == "-"], 
                file.path(baseDir, "human_work/bedgraphs/K562_100ng_2_-.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$K562_100ng_3[strand(exp@counts$TSSs$cpm$K562_100ng_3) == "+"], 
                file.path(baseDir, "human_work/bedgraphs/K562_100ng_3_+.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$K562_100ng_3[strand(exp@counts$TSSs$cpm$K562_100ng_3) == "-"], 
                file.path(baseDir, "human_work/bedgraphs/K562_100ng_3_-.bedgraph"))

# Normalize TSR counts
exp <- count_normalization(exp, data_type = "tsr", threshold = 3, n_samples = 1, samples = stripe)

# Generate a combinbed TSR correlation plot
p <- plot_correlation(exp, data_type = "tsr") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_correlation.png", plot = p, device = cairo_pdf, height = 12, width = 12)

# Generate a hierarchically clustered TSR heatmap
p <- plot_correlation(exp, data_type = "tsr", correlation_plot = "hierarchical", col = viridis(256), 
                      correlation_metric = "pearson")

pdf("tsr_correlation_hierarchical.pdf", height = 7, width = 7)
p
dev.off()

# Annotate TSRs
exp <- annotate_features(exp, annotation_file = file.path(baseDir, "human_data/Homo_sapiens.GRCh38.98.chr.gtf"),
                         data_type = "tsr", feature_type = "transcript")

# Determine TSR distribution relative to genomic features
tsr_distribution <- genomic_distribution(exp, data_type = "tsr", threshold = 3, samples = "K562_100ng_1")

p <- plot_genomic_distribution(tsr_distribution) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_genomic_distribution.pdf", plot = p, device = cairo_pdf, height = 1.5, width = 4)

tsr_distribution <- genomic_distribution(exp, data_type = "tsr", threshold = 3, quantiles = 5, 
                                         samples = "K562_100ng_1")

p <- plot_genomic_distribution(tsr_distribution) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_genomic_distribution_quantiles.pdf", plot = p, device = cairo_pdf, height = 1.5, width = 4)

# Plot number of promoter-proximal features with a TSR
features <- detect_features(exp, data_type = "tsr", feature_type = "transcript", samples = "K562_100ng_1")

p <- plot_detected_features(features) +
    ggplot2::theme(text = element_text(size = 5))

ggsave("tsr_feature_plot.pdf", plot = p, device = cairo_pdf, height = 2, width = 4)

# Plot selected TSR metrics
p <- plot_tsr_metric(exp, tsr_metrics = "nTSSs", log2_transform = TRUE, ncol = 1, plot_type = "boxjitter",
                     size = 0.5, samples = "K562_100ng_1") +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_metrics.pdf", plot = p, device = cairo_pdf, width = 7, height = 7)

# Generate TSR average plot
p <- plot_average(exp, data_type = "tsr", samples = "K562_100ng_1") +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_average_plot.pdf", plot = p, device = cairo_pdf, height = 3, width = 3)

####################################
### STRIPE-seq vs. CAGE analysis ###
####################################

if (!dir.exists(file.path(baseDir, "human_work/CAGE/"))){
  print("Creating directory 'human_work/CAGE' and changing working directory...")
  dir.create(file.path(baseDir, "human_work/CAGE/"))
  setwd(file.path(baseDir, "human_work/CAGE/"))
} else {
  print("Directory 'human_work/CAGE' already exists, changing working directory...")
  setwd(file.path(baseDir, "human_work/CAGE/"))
}

# Normalize TSS counts
exp <- count_normalization(exp, data_type = "tss", threshold = 3, n_samples = 1, samples = all)

# Assess TSS dinucleotide frequencies
frequencies <- dinucleotide_frequencies(exp, genome_assembly = file.path(baseDir, "human_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
                                        threshold = 3, samples = c("K562_100ng_1",
                                                                   "nanoCAGE_XL_7.5ug_1",
                                                                   "CAGE_10ug_1",
                                                                   "RAMPAGE_5ug_1"))

p <- plot_dinucleotide_frequencies(frequencies, ncol = 2) +
  ggplot2::theme(text = element_text(size = 6))

ggsave("tss_dinucleotide_frequencies.pdf", plot = p, device = cairo_pdf, height = 3, width = 6)

# Export normalized TSS bedGraphs

# nanoCAGE-XL
export.bedGraph(exp@counts$TSSs$cpm$nanoCAGE_XL_7.5ug_1[strand(exp@counts$TSSs$cpm$nanoCAGE_XL_7.5ug_1) == "+"], 
                file.path(baseDir, "human_work/bedgraphs/nanoCAGE_XL_7.5ug_1_+.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$nanoCAGE_XL_7.5ug_1[strand(exp@counts$TSSs$cpm$nanoCAGE_XL_7.5ug_1) == "-"], 
                file.path(baseDir, "human_work/bedgraphs/nanoCAGE_XL_7.5ug_1_-.bedgraph"))

# CAGE
export.bedGraph(exp@counts$TSSs$cpm$CAGE_10ug_1[strand(exp@counts$TSSs$cpm$CAGE_10ug_1) == "+"], 
                file.path(baseDir, "human_work/bedgraphs/CAGE_10ug_1_+.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$CAGE_10ug_1[strand(exp@counts$TSSs$cpm$CAGE_10ug_1) == "-"], 
                file.path(baseDir, "human_work/bedgraphs/CAGE_10ug_1_-.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$CAGE_10ug_2[strand(exp@counts$TSSs$cpm$CAGE_10ug_2) == "+"], 
                file.path(baseDir, "human_work/bedgraphs/CAGE_10ug_2_+.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$CAGE_10ug_2[strand(exp@counts$TSSs$cpm$CAGE_10ug_2) == "-"], 
                file.path(baseDir, "human_work/bedgraphs/CAGE_10ug_2_-.bedgraph"))

# RAMPAGE
export.bedGraph(exp@counts$TSSs$cpm$RAMPAGE_5ug_1[strand(exp@counts$TSSs$cpm$RAMPAGE_5ug_1) == "+"], 
                file.path(baseDir, "human_work/bedgraphs/RAMPAGE_5ug_1_+.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$RAMPAGE_5ug_1[strand(exp@counts$TSSs$cpm$RAMPAGE_5ug_1) == "-"], 
                file.path(baseDir, "human_work/bedgraphs/RAMPAGE_5ug_1_-.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$RAMPAGE_5ug_2[strand(exp@counts$TSSs$cpm$RAMPAGE_5ug_2) == "+"], 
                file.path(baseDir, "human_work/bedgraphs/RAMPAGE_5ug_2_+.bedgraph"))
export.bedGraph(exp@counts$TSSs$cpm$RAMPAGE_5ug_2[strand(exp@counts$TSSs$cpm$RAMPAGE_5ug_2) == "-"], 
                file.path(baseDir, "human_work/bedgraphs/RAMPAGE_5ug_2_-.bedgraph"))

# Normalize TSR counts
exp <- count_normalization(exp, data_type = "tsr", threshold = 3, n_samples = 1, samples = all)

# Generate a combinbed TSR correlation plot
p <- plot_correlation(exp, data_type = "tsr", correlation_metric = "spearman") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_correlation.png", plot = p, device = "png", type = "cairo", height = 20, width = 20)

# Generate a hierarchically clustered TSR heatmap
p <- plot_correlation(exp, data_type = "tsr", correlation_plot = "hierarchical", col = viridis(256), 
                      correlation_metric = "spearman")

pdf("tsr_correlation_hierarchical.pdf", height = 7, width = 7)
p
dev.off()

# Annotate TSRs
exp <- annotate_features(exp, annotation_file = file.path(baseDir, "human_data/Homo_sapiens.GRCh38.98.chr.gtf"),
                         data_type = "tsr", feature_type = "transcript")

# Determine TSR distribution relative to genomic features
tsr_distribution <- genomic_distribution(exp, data_type = "tsr", threshold = 3, 
                                         samples = c("K562_100ng_1","K562_100ng_2","K562_100ng_3",
                                                     "nanoCAGE_XL_7.5ug_1",
                                                     "CAGE_10ug_1","CAGE_10ug_2",
                                                     "RAMPAGE_5ug_1","RAMPAGE_5ug_2"))

p <- plot_genomic_distribution(tsr_distribution) +
    ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_genomic_distribution.pdf", plot = p, device = cairo_pdf, height = 1.5, width = 4)

# Plot number of promoter-proximal features with a TSR
features <- detect_features(exp, data_type = "tsr", feature_type = "transcript", 
                            samples = c("K562_100ng_1",
                                        "nanoCAGE_XL_7.5ug_1",
                                        "CAGE_10ug_1",
                                        "RAMPAGE_5ug_1"))

p <- plot_detected_features(features) +
  ggplot2::theme(text = element_text(size = 5))

ggsave("tsr_feature_plot.pdf", plot = p, device = cairo_pdf, height = 3, width = 4)