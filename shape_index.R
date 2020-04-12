library(tsrexplorer)
library(tidyverse)
library(GenomicRanges)
library(viridis)

# Note that the newest version of TSRexploreR is required for SI calculation and plotting
# devtools::install_github("rpolicastro/tsrexplorer", ref = "clean", force = TRUE)

#################
### S288C YPD ###
#################

# Read STRIPE-seq and CAGE TSSs and TSRs into variables
yeast_STRIPE_TSSs <- map(list.files("yeast_data/YPD_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                  makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                           start.field = "TSS", end.field = "TSS")) %>%
  set_names(c("S288C_50ng_1","S288C_50ng_2","S288C_50ng_3",
              "S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
              "S288C_250ng_1","S288C_250ng_2","S288C_250ng_3"))

yeast_CAGE_TSSs <- map(list.files("yeast_data/CAGE_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                   makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                            start.field = "TSS", end.field = "TSS")) %>%
  set_names(c("SLIC_CAGE_100ng_1","SLIC_CAGE_100ng_2",
              "nanoCAGE_500ng_1","nanoCAGE_500ng_2",
              "nanoCAGE_25ng_1","nanoCAGE_25ng_2"))

yeast_STRIPE_TSRs <- map(list.files("yeast_data/YPD_TSRs/", full.names = TRUE), ~ read.delim(.x) %>%
                  makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                           start.field = "start", end.field = "end")) %>%
  set_names(c("S288C_50ng_1","S288C_50ng_2","S288C_50ng_3",
              "S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
              "S288C_250ng_1","S288C_250ng_2","S288C_250ng_3"))

yeast_CAGE_TSRs <- map(list.files("yeast_data/CAGE_TSRs/", full.names = TRUE), ~ read.delim(.x) %>%
                   makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                            start.field = "start", end.field = "end")) %>%
  set_names(c("SLIC_CAGE_100ng_1","SLIC_CAGE_100ng_2",
              "nanoCAGE_500ng_1","nanoCAGE_500ng_2",
              "nanoCAGE_25ng_1","nanoCAGE_25ng_2"))

# Combine lists of STRIPE-seq and CAGE TSSs and TSRs
yeast_TSSs <- c(yeast_STRIPE_TSSs, yeast_CAGE_TSSs)

yeast_TSRs <- c(yeast_STRIPE_TSRs, yeast_CAGE_TSRs)

# Add a "score" column to TSS and TSR GRanges (duplication and column renaming of TSRchitect "nTAGs" column)
yeast_TSSs <- lapply(yeast_TSSs, function(x) {x$score <- x$nTAGs; return(x)})

yeast_TSRs <- lapply(yeast_TSRs, function(x) {x$score <- x$nTAGs; return(x)})

# Generate list of all samples for later ordering
yeast_STRIPE_CAGE <- c("S288C_50ng_1","S288C_50ng_2","S288C_50ng_3",
                       "S288C_100ng_1","S288C_100ng_2","S288C_100ng_3",
                       "S288C_250ng_1","S288C_250ng_2","S288C_250ng_3",
                       "SLIC_CAGE_100ng_1","SLIC_CAGE_100ng_2",
                       "nanoCAGE_500ng_1","nanoCAGE_500ng_2",
                       "nanoCAGE_25ng_1","nanoCAGE_25ng_2")

# Read data into TSRexploreR
yeast_exp <- tsr_explorer(yeast_TSSs, yeast_TSRs)

# Format TSS and TSR data for analysis
yeast_exp <- format_counts(yeast_exp, data_type = "tss")

yeast_exp <- format_counts(yeast_exp, data_type = "tsr")

# Normalize TSS counts
yeast_exp <- cpm_normalize(yeast_exp, data_type = "tss")

# Associate TSSs with TSRs
yeast_exp <- associate_with_tsr(
  yeast_exp, use_sample_sheet = FALSE,
  sample_list <- list(
    "S288C_50ng_1" = "S288C_50ng_1",
    "S288C_50ng_2" = "S288C_50ng_2",
    "S288C_50ng_3" = "S288C_50ng_3",
    "S288C_100ng_1" = "S288C_100ng_1",
    "S288C_100ng_2" = "S288C_100ng_2",
    "S288C_100ng_3" = "S288C_100ng_3",
    "S288C_250ng_1" = "S288C_250ng_1",
    "S288C_250ng_2" = "S288C_250ng_2",
    "S288C_250ng_3" = "S288C_250ng_3",
    "SLIC_CAGE_100ng_1" = "SLIC_CAGE_100ng_1",
    "SLIC_CAGE_100ng_2" = "SLIC_CAGE_100ng_2",
    "nanoCAGE_500ng_1" = "nanoCAGE_500ng_1",
    "nanoCAGE_500ng_2" = "nanoCAGE_500ng_2",
    "nanoCAGE_25ng_1" = "nanoCAGE_25ng_1",
    "nanoCAGE_25ng_2" = "nanoCAGE_25ng_2"
  )
)

# Calculate TSR metrics
yeast_exp <- tsr_metrics(yeast_exp)

# Plot shape index
ysi <- plot_tsr_metric(yeast_exp, tsr_metrics = "shape_index", log2_transform = FALSE, ncol = 1, 
                       threshold = 10, samples = yeast_STRIPE_CAGE) +
  theme_bw() +
  ggplot2::theme(text = element_text(size = 10), legend.key.size = unit(0.5, "cm")) +
  geom_hline(yintercept = -1, size = 0.4, lty = 2)

ggsave("yeast_shape_index.pdf", plot = ysi, device = cairo_pdf, height = 3.5, width = 8.5)

# Count TSRs in each shape class, only considering those with at least 10 counts
yeast_df <- bind_rows(yeast_exp@counts$TSRs$raw, .id = "sample")

yeast_df <- filter(yeast_df, score >= 10)

count(yeast_df, sample, shape_class) %>%
  group_by(sample) %>%
  mutate(total = sum(n)) %>%
  mutate(percent = round(n / total * 100, 2)) %>%
  print(n = Inf)

############
### K562 ###
############

# Read STRIPE-seq and CAGE TSSs and TSRs into variables
human_STRIPE_TSSs <- map(list.files("human_data/STRIPE_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                     makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                              start.field = "TSS", end.field = "TSS")) %>%
  set_names(c("K562_100ng_1","K562_100ng_2","K562_100ng_3"))

human_CAGE_TSSs <- map(list.files("human_data/CAGE_TSSs/", full.names = TRUE), ~ read.delim(.x) %>%
                   makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                            start.field = "TSS", end.field = "TSS")) %>%
  set_names(c("nanoCAGE_XL_7.5ug_1",
              "CAGE_10ug_1","CAGE_10ug_2",
              "RAMPAGE_5ug_1","RAMPAGE_5ug_2"))

human_STRIPE_TSRs <- map(list.files("human_data/STRIPE_TSRs/", full.names = TRUE), ~ read.delim(.x) %>%
                     makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                              start.field = "start", end.field = "end")) %>%
  set_names(c("K562_100ng_1","K562_100ng_2","K562_100ng_3"))

human_CAGE_TSRs <- map(list.files("human_data/CAGE_TSRs/", full.names = TRUE), ~ read.delim(.x) %>%
                   makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "seq", 
                                            start.field = "start", end.field = "end")) %>%
  set_names(c("nanoCAGE_XL_7.5ug_1",
              "CAGE_10ug_1","CAGE_10ug_2",
              "RAMPAGE_5ug_1","RAMPAGE_5ug_2"))

# Combine lists of STRIPE-seq and CAGE TSSs and TSRs
human_TSSs <- c(human_STRIPE_TSSs, human_CAGE_TSSs)

human_TSRs <- c(human_STRIPE_TSRs, human_CAGE_TSRs)

# Add a "score" column to TSS and TSR GRanges (duplication and column renaming of TSRchitect "nTAGs" column)
human_TSSs <- lapply(human_TSSs, function(x) {x$score <- x$nTAGs; return(x)})

human_TSRs <- lapply(human_TSRs, function(x) {x$score <- x$nTAGs; return(x)})

# Generate list of all samples for later ordering
human_STRIPE_CAGE <- c("K562_100ng_1","K562_100ng_2","K562_100ng_3",
                       "CAGE_10ug_1","CAGE_10ug_2",
                       "RAMPAGE_5ug_1","RAMPAGE_5ug_2",
                       "nanoCAGE_XL_7.5ug_1")

# Read data into TSRexploreR
human_exp <- tsr_explorer(human_TSSs, human_TSRs)

# Format TSS and TSR data for analysis
human_exp <- format_counts(human_exp, data_type = "tss")

human_exp <- format_counts(human_exp, data_type = "tsr")

# Normalize TSS counts
human_exp <- cpm_normalize(human_exp, data_type = "tss")

# Associate TSSs with TSRs
human_exp <- associate_with_tsr(
  human_exp, use_sample_sheet = FALSE,
  sample_list <- list(
    "K562_100ng_1" = "K562_100ng_1",
    "K562_100ng_2" = "K562_100ng_2",
    "K562_100ng_3" = "K562_100ng_3",
    "nanoCAGE_XL_7.5ug_1" = "nanoCAGE_XL_7.5ug_1",
    "CAGE_10ug_1" = "CAGE_10ug_1",
    "CAGE_10ug_2" = "CAGE_10ug_2",
    "RAMPAGE_5ug_1" = "RAMPAGE_5ug_1",
    "RAMPAGE_5ug_2" = "RAMPAGE_5ug_2"
  )
)

# Calculate TSR metrics
human_exp <- tsr_metrics(human_exp)

# Plot shape index
hsi <- plot_tsr_metric(human_exp, tsr_metrics = "shape_index", log2_transform = FALSE, ncol = 1, 
                       threshold = 10, samples = human_STRIPE_CAGE) +
  theme_bw() +
  ggplot2::theme(text = element_text(size = 10), legend.key.size = unit(0.5, "cm")) + 
  geom_hline(yintercept = -1, size = 0.4, lty = 2)

ggsave("human_shape_index.pdf", plot = hsi, device = cairo_pdf, height = 3.5, width = 6.5)

# Count TSRs in each shape class, only considering those with at least 10 counts=
human_df <- bind_rows(human_exp@counts$TSRs$raw, .id = "sample")

human_df <- filter(human_df, score >= 10)

count(human_df, sample, shape_class) %>%
  group_by(sample) %>%
  mutate(total = sum(n)) %>%
  mutate(percent = round(n / total * 100, 2)) %>%
  print(n = Inf)
