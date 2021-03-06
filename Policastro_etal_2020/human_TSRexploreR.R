library(tsrexplorer)
library(tidyverse)
library(GenomicRanges)
library(ComplexHeatmap)
library(viridis)
library(rtracklayer)

# Pull latest version of tsrexplorer
# devtools::install_github("rpolicastro/tsrexplorer", ref = "dev", force = TRUE)

# Set base directory for analysis
if (!dir.exists("human_work")){
  message("Creating directory 'human_work'...")
  dir.create("human_work")
} else {
  message("Directory 'human_work' already exists.")
}

####################
### Read in TSSs ###
####################

# STRIPE-seq TSSs
STRIPE_TSSs <- map(list.files("human_data/STRIPE_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                      makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                             start.field = "TSS", end.field = "TSS")) %>%
    set_names(c("K562_100ng_1","K562_100ng_2","K562_100ng_3"))

# nanoCAGE-XL, CAGE, and RAMPAGE TSSs
CAGE_TSSs <- map(list.files("human_data/CAGE_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
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
                      makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                               start.field = "start", end.field = "end")) %>%
  set_names(c("K562_100ng_1","K562_100ng_2","K562_100ng_3"))

# nanoCAGE-XL, CAGE, and RAMPAGE TSRs
CAGE_TSRs <- map(list.files("human_data/CAGE_TSRs/", full.names = TRUE), ~ read.delim(.x) %>%
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

# All samples
all <- c("K562_100ng_1","K562_100ng_2","K562_100ng_3",
         "CAGE_10ug_1","CAGE_10ug_2",
         "RAMPAGE_5ug_1","RAMPAGE_5ug_2",
         "nanoCAGE_XL_7.5ug_1")

# Read TSSs and TSRs into TSRexploreR
exp <- tsr_explorer(full_TSS_set,full_TSR_set)

###############################
### K562 replicate analysis ###
###############################

stripe_dir <- file.path("human_work", "STRIPE")

if (!dir.exists(stripe_dir)){
  print("Creating directory 'human_work/STRIPE'")
  dir.create(stripe_dir)
} else {
  print("Directory 'human_work/STRIPE' already exists.")
}

# Normalize TSS counts
exp <- count_normalization(exp, data_type = "tss", threshold = 3, n_samples = 1, samples = stripe)

# Generate a hierarchically clustered TSS heatmap with correlation values displayed
corr_matrix <- find_correlation(exp, data_type = "tss", correlation_metric = "pearson", samples = stripe)

cairo_pdf(file = file.path(stripe_dir, "tss_correlation_hierarchical.pdf"), width = 5, height = 5)
Heatmap(corr_matrix, col = viridis(256), heatmap_legend_param = list(title = "PCC"), 
        layer_fun = function(j, i, x, y, width, height, fill)
        {
            grid.text(sprintf("%.3f", pindex(corr_matrix, i, j)), x, y, gp = gpar(fontsize = 22, col = "white"))
        }
)
dev.off()

# Annotate TSSs relative to genomic features
exp <- annotate_features(exp, annotation_file = "human_data/Homo_sapiens.GRCh38.98.chr.gtf",
                         data_type = "tss", feature_type = "transcript", upstream = 500, downstream = 500)

# Explore TSS read thresholds for promoter fraction and plot
thresh <- explore_thresholds(exp, annotation_file = "human_data/Homo_sapiens.GRCh38.98.chr.gtf", 
                             feature_type = "transcript", max_threshold = 25, 
                             upstream = 500, downstream = 500, samples = stripe)

p <- plot_threshold_exploration(thresh, ncol = 3, point_size = 1.5, sample_order = stripe) +
    ggplot2::geom_vline(xintercept = 3, lty = 2) +
    ggplot2::theme(legend.key.size = unit(0.6, "cm"), text = element_text(size = 12))

ggsave(file.path(stripe_dir, "tss_thresholds.pdf"), plot = p, device = cairo_pdf, height = 2, width = 10)

# Determine TSS distribution relative to genomic features
tss_distribution <- genomic_distribution(exp, data_type = "tss", threshold = 3, samples = stripe)

p <- plot_genomic_distribution(tss_distribution, sample_order = stripe) +
    ggplot2::theme(text = element_text(size = 14), legend.key.size = unit(0.8, "cm"))

ggsave(file.path(stripe_dir, "tss_genomic_distribution.pdf"), plot = p, device = cairo_pdf, height = 2.5, width = 6)

genomic_dist <- genomic_distribution(exp, data_type = "tss", threshold = 3, quantiles = 5, 
                                     samples = "K562_100ng_1")

p <- plot_genomic_distribution(genomic_dist) +
    ggplot2::theme(text = element_text(size = 14), legend.key.size = unit(0.8, "cm"))

ggsave(file.path(stripe_dir, "tss_genomic_distribution_quantiles.pdf"), plot = p, device = cairo_pdf, height = 4, width = 6)

# Generate TSS density plot
p <- plot_average(exp, data_type = "tss", threshold = 3, samples = "K562_100ng_1", upstream = 1000, downstream = 1000) +
    ggplot2::theme(text = element_text(size = 13))

ggsave(file.path(stripe_dir, "tss_average_plot.pdf"), plot = p, cairo_pdf, height = 2.5, width = 3.5)

# Generate TSS sequence logos
seqs <- tss_sequences(exp, genome_assembly = "human_data/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
                      threshold = 3, samples = stripe)

p <- plot_sequence_logo(seqs, ncol = 1, font_size = 10)

ggsave(file.path(stripe_dir, "tss_seq_logo.pdf"), plot = p, device = cairo_pdf, height = 5, width = 4)

seqs <- tss_sequences(exp, genome_assembly = "human_data/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
                      threshold = 3, quantiles = 5, samples = "K562_100ng_1")

p <- plot_sequence_logo(seqs, ncol = 1, font_size = 10)

ggsave(file.path(stripe_dir, "tss_seq_logo_quantiles.pdf"), plot = p, device = cairo_pdf, height = 6, width = 4)

# Generate TSS color plot
seqs <- tss_sequences(exp, genome_assembly = "human_data/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
                      threshold = 3, samples = "K562_100ng_1")

p <- plot_sequence_colormap(seqs, ncol = 3) +
    ggplot2::theme(text = element_text(size = 6), legend.key.size = unit(0.4, "cm"))

ggsave(file.path(stripe_dir, "tss_seq_colormap.png"), plot = p, device = "png", type = "cairo", height = 2.5, width = 2)

# Assess TSS dinucleotide frequencies
frequencies <- dinucleotide_frequencies(exp, genome_assembly = "human_data/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
                                        threshold = 3, samples = stripe)

p <- plot_dinucleotide_frequencies(frequencies, ncol = 1, sample_order = stripe) +
    ggplot2::theme(text = element_text(size = 12), legend.key.size = unit(0.6, "cm"))

ggsave(file.path(stripe_dir, "tss_dinucleotide_frequencies.pdf"), plot = p, device = cairo_pdf, height = 7, width = 3.5)

frequencies <- dinucleotide_frequencies(exp, genome_assembly = "human_data/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
                                        threshold = 3, samples = "K562_100ng_1")

p <- plot_dinucleotide_frequencies(frequencies, ncol = 1) +
    ggplot2::theme(text = element_text(size = 12), legend.key.size = unit(0.6, "cm"))

ggsave(file.path(stripe_dir, "tss_dinucleotide_frequencies_rep1.pdf"), plot = p, device = cairo_pdf, height = 3, width = 3.8)

# Export normalized TSS bedGraphs

bedgraph_dir <- file.path("human_work", "bedgraphs")

if (!dir.exists(bedgraph_dir)){
  message("Creating directory 'human_work/bedgraphs'...")
  dir.create(bedgraph_dir)
} else {
  message("Directory 'human_work/bedgraphs' already exists.")
}

iwalk(exp@counts$TSSs$cpm, function(counts, sample) {
	plus <- counts[strand(counts) == "+"]
	minus <- counts[strand(counts) == "-"]

	export(plus, file.path(bedgraph_dir, paste(sample, "plus.bedgraph", sep = "_")), format = "bedgraph")
	export(minus, file.path(bedgraph_dir, paste(sample, "minus.bedgraph", sep = "_")), format = "bedgraph")
})

# Normalize TSR counts
exp <- count_normalization(exp, data_type = "tsr", threshold = 3, n_samples = 1, samples = stripe)

# Export TSR beds
iwalk(exp@counts$TSRs$cpm, function(counts, sample) {
    export(counts, file.path(bedgraph_dir, paste(sample, "TSRs.bed", sep = "_")), format = "bed")
})

####################################
### STRIPE-seq vs. CAGE analysis ###
####################################

cage_dir <- file.path("human_work", "CAGE")

if (!dir.exists(cage_dir)){
  message("Creating directory 'human_work/CAGE'...")
  dir.create(cage_dir)
} else {
  message("Directory 'human_work/CAGE' already exists.")
}

# Normalize TSS counts
exp <- count_normalization(exp, data_type = "tss", threshold = 3, n_samples = 1, samples = all)

# Annotate TSSs
exp <- annotate_features(exp, annotation_file = "human_data/Homo_sapiens.GRCh38.98.chr.gtf",
                         data_type = "tss", feature_type = "transcript", upstream = 500, downstream = 500)

# Explore TSS read thresholds for promoter fraction and plot
thresh <- explore_thresholds(exp, annotation_file = "human_data/Homo_sapiens.GRCh38.98.chr.gtf", 
                             feature_type = "transcript", max_threshold = 25, 
                             upstream = 500, downstream = 500, samples = all)

p <- plot_threshold_exploration(thresh, ncol = 3, point_size = 1.5, sample_order = all) +
    ggplot2::geom_vline(xintercept = 3, lty = 2) +
    ggplot2::theme(legend.key.size = unit(0.8, "cm"), text = element_text(size = 12))

ggsave(file.path(cage_dir, "tss_thresholds.pdf"), plot = p, device = cairo_pdf, height = 5, width = 10)

# Assess TSS dinucleotide frequencies
frequencies <- dinucleotide_frequencies(exp, genome_assembly = "human_data/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
                                        threshold = 3, samples = all)

p <- plot_dinucleotide_frequencies(frequencies, ncol = 3, sample_order = all) +
    ggplot2::theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"))

ggsave(file.path(cage_dir, "tss_dinucleotide_frequencies.pdf"), plot = p, device = cairo_pdf, height = 8, width = 7)

# Export normalized TSS bedGraphs

bedgraph_dir <- file.path("human_work", "bedgraphs")

iwalk(exp@counts$TSSs$cpm, function(counts, sample) {
	plus <- counts[strand(counts) == "+"]
	minus <- counts[strand(counts) == "-"]

	export(plus, file.path(bedgraph_dir, paste(sample, "plus.bedgraph", sep = "_")), format = "bedgraph")
	export(minus, file.path(bedgraph_dir, paste(sample, "minus.bedgraph", sep = "_")), format = "bedgraph")
})

# Normalize TSR counts
exp <- count_normalization(exp, data_type = "tsr", threshold = 3, n_samples = 1, samples = all)

# Annotate TSRs
exp <- annotate_features(exp, annotation_file = "human_data/Homo_sapiens.GRCh38.98.chr.gtf",
                         data_type = "tsr", feature_type = "transcript", upstream = 500, downstream = 500)

# Determine average promoter feature counts and plot
features <- detect_features(exp, data_type = "tsr", feature_type = "transcript", 
                            samples = all, downstream = 500, upstream = 500)

p <- features %>% 
    add_column(technology = c(rep("CAGE", 2), rep("STRIPE", 3), "nanoCAGEXL", rep("RAMPAGE", 2))) %>%
    mutate(technology = factor(technology, levels=c("STRIPE", "CAGE", "RAMPAGE", "nanoCAGEXL"))) %>%
    ggplot(aes(x = technology, y = promoter_proximal_features)) +
    theme_bw() +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), size = 0.4, geom = "errorbar", width = 0.25,
                 alpha = 0.5) +
    geom_jitter(aes(color = technology), position = position_jitter(height = 0, width = 0.2), size = 2) +
    scale_color_viridis_d(begin = 0, end = 0.7) +
    ylim(5000, 15000)

ggsave(file.path(cage_dir, "promoter_feature_jitter_plot.pdf"), plot = p, device = cairo_pdf, height = 4, width = 6)

# Export TSR beds
iwalk(exp@counts$TSRs$cpm, function(counts, sample) {
    export(counts, file.path(bedgraph_dir, paste(sample, "TSRs.bed", sep = "_")), format = "bed")
})

########################
### RNA-seq analysis ###
########################

rnaseq_dir <- file.path("human_work", "RNA_seq")

if (!dir.exists(rnaseq_dir)){
    message("Creating directory 'human_work/RNA_seq'...")
    dir.create(rnaseq_dir)
} else {
    message("Directory 'human_work/RNA_seq' already exists.")
}

# Get feature counts for STRIPE-seq and RNA-seq
stripe_counts <- read.delim(
	file.path("human_data", "RNA_seq", "cleaned_K562_feature_counts_STRIPEseq.tsv"),
	sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

rnaseq_counts <- read.delim(
	file.path("human_data", "RNA_seq", "cleaned_human_feature_counts.tsv"),
	sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

stripe_counts <- column_to_rownames(stripe_counts, "Geneid") %>% as.matrix
rnaseq_counts <- column_to_rownames(rnaseq_counts, "Geneid") %>% as.matrix

exp <- add_feature_counts(exp, five_prime_feature_counts = stripe_counts, rnaseq_feature_counts = rnaseq_counts)
exp <- count_normalization(exp, data_type = "features")

# Group samples for analysis
stripe_rnaseq_k562 <- c("GSF2268_s_SP52_K562_WT_100ng",	"GSF2268_s_SP53_K562_WT_100ng", "GSF2268_s_SP54_K562_WT_100ng",
                        "RNASEQ004_K562_untreated_r1", "RNASEQ005_K562_untreated_r1", "RNASEQ006_K562_untreated_r1")

# Find correlation of YPD STRIPE-seq and RNA-seq feature counts
corr_matrix <- find_correlation(exp, data_type = "features", correlation_metric = "spearman", samples = stripe_rnaseq_k562)

cairo_pdf(file = file.path(rnaseq_dir, "tss_rnaseq_correlation_hierarchical.pdf"), width = 10, height = 10)
Heatmap(corr_matrix, col = viridis(256), heatmap_legend_param = list(title = "spearman"), 
        layer_fun = function(j, i, x, y, width, height, fill)
        {
            grid.text(sprintf("%.3f", pindex(corr_matrix, i, j)), x, y, gp = gpar(fontsize = 22, col = "white"))
        }
)
dev.off()