library(GenomicFeatures)
library(Gviz)
library(viridis)
library(scales)

# This script imports signal tracks and displays them at the specified regions of the human 
# genome, listed at the end of the script.

gviz_dir <- file.path("human_work", "Gviz")
bedgraph_dir <- file.path("human_work", "bedgraphs")
rnaseq_dir <- file.path("human_data", "RNA_seq")
cutntag_dir <- file.path("human_data", "CUTnTag")

if (!dir.exists(gviz_dir)) {
    message("Creating directory 'human_work/Gviz'...")
    dir.create(gviz_dir)
} else {
    message("Directory 'human_work/Gviz' already exists.")
}

# Create genomic axis track
axis.track <- GenomeAxisTrack(col = "black", scale = 0.1, col.range = "black")

options(ucscChromosomeNames = FALSE)

# Create gene annotation track
txdb <- makeTxDbFromGFF(file.path("human_data", "Homo_sapiens.GRCh38.98.chr.gtf"))

genome.track <- GeneRegionTrack(txdb, genome = "GRCh38", shape = "arrow", names = "Genes", col = "black",
                                showId = TRUE, fill = "black", trancriptAnnotation = "gene_symbol")

# Create data tracks

# Get colors from the viridis palette for the number of tracks to be plotted
show_col(viridis_pal()(20))

# Set y-axis limits
stripe_pos_lim <- c(0,20)
stripe_neg_lim <- c(0,20)

cage_pos_lim <- c(0,5)
cage_neg_lim <- c(0,5)

rampage_pos_lim <- c(0,20)
rampage_neg_lim <- c(0,20)

nanocage_pos_lim <- c(0,5)
nanocage_neg_lim <- c(0,5)

rnaseq_pos_lim <- c(0,1)
rnaseq_neg_lim <- c(0,1)

h3k4cnt_lim <- c(0,50)
h3k27cnt_lim <- c(0,0.5)

# STRIPE-seq 
K562_100ng_1_pos <- DataTrack(range = file.path(bedgraph_dir, "K562_100ng_1_plus.bedgraph"), genome = "GRCh38", 
                              name = "100 ng 1 plus", col.histogram = "#440154FF", fill.histogram = "#440154FF", ylim = stripe_pos_lim)
K562_100ng_1_neg <- DataTrack(range = file.path(bedgraph_dir, "K562_100ng_1_minus.bedgraph"), genome = "GRCh38", 
                              name = "100 ng 1 minus", col.histogram = "#440154FF", fill.histogram = "#440154FF", ylim = stripe_neg_lim)
K562_100ng_2_pos <- DataTrack(range = file.path(bedgraph_dir, "K562_100ng_2_plus.bedgraph"), genome = "GRCh38", 
                              name = "100 ng 2 plus", col.histogram = "#440154FF", fill.histogram = "#440154FF", ylim = stripe_pos_lim)
K562_100ng_2_neg <- DataTrack(range = file.path(bedgraph_dir, "K562_100ng_2_minus.bedgraph"), genome = "GRCh38", 
                              name = "100 ng 2 minus", col.histogram = "#440154FF", fill.histogram = "#440154FF", ylim = stripe_neg_lim)
K562_100ng_3_pos <- DataTrack(range = file.path(bedgraph_dir, "K562_100ng_3_plus.bedgraph"), genome = "GRCh38", 
                              name = "100 ng 3 plus", col.histogram = "#440154FF", fill.histogram = "#440154FF", ylim = stripe_pos_lim)
K562_100ng_3_neg <- DataTrack(range = file.path(bedgraph_dir, "K562_100ng_3_minus.bedgraph"), genome = "GRCh38", 
                              name = "100 ng 3 minus", col.histogram = "#440154FF", fill.histogram = "#440154FF", ylim = stripe_neg_lim)

# CAGE
CAGE_10ug_1_pos <- DataTrack(range = file.path(bedgraph_dir, "CAGE_10ug_1_plus.bedgraph"), genome = "GRCh38", 
                             name = "CAGE 10 1 plus", col.histogram = "#56C667FF", fill.histogram = "#56C667FF", ylim = cage_pos_lim)
CAGE_10ug_1_neg <- DataTrack(range = file.path(bedgraph_dir, "CAGE_10ug_1_minus.bedgraph"), genome = "GRCh38", 
                             name = "CAGE 10 1 minus", col.histogram = "#56C667FF", fill.histogram = "#56C667FF", ylim = cage_neg_lim)

# RAMPAGE
RAMPAGE_5ug_1_pos <- DataTrack(range = file.path(bedgraph_dir, "RAMPAGE_5ug_1_plus.bedgraph"), genome = "GRCh38", 
                             name = "RAMPAGE 5 ug 1 plus", col.histogram = "#94D840FF", fill.histogram = "#94D840FF", ylim = rampage_pos_lim)
RAMPAGE_5ug_1_neg <- DataTrack(range = file.path(bedgraph_dir, "RAMPAGE_5ug_1_minus.bedgraph"), genome = "GRCh38", 
                             name = "RAMPAGE 5 ug 1 minus", col.histogram = "#94D840FF", fill.histogram = "#94D840FF", ylim = rampage_neg_lim)

# nanoCAGE-XL
nanoCAGE_XL_7.5ug_1_pos <- DataTrack(range = file.path(bedgraph_dir, "nanoCAGE_XL_7.5ug_1_plus.bedgraph"), genome = "GRCh38", 
                               name = "nanoCAGE-XL 7.5 ug plus", col.histogram = "#DCE318FF", fill.histogram = "#DCE318FF", ylim = nanocage_pos_lim)
nanoCAGE_XL_7.5ug_1_neg <- DataTrack(range = file.path(bedgraph_dir, "nanoCAGE_XL_7.5ug_1_minus.bedgraph"), genome = "GRCh38", 
                               name = "nanoCAGE-XL 7.5 ug minus", col.histogram = "#DCE318FF", fill.histogram = "#DCE318FF", ylim = nanocage_neg_lim)

# RNA-seq
K562_rnaseq_1_pos <- DataTrack(range = file.path(rnaseq_dir, "RNASEQ004_K562_untreated_r1.CPM.bs1.smooth25.plus.bw"), genome = "GRCh38", 
                              name = "K562 RNA-seq 1 plus", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = rnaseq_pos_lim)
K562_rnaseq_1_neg <- DataTrack(range = file.path(rnaseq_dir, "RNASEQ004_K562_untreated_r1.CPM.bs1.smooth25.minus.bw"), genome = "GRCh38", 
                               name = "K562 RNA-seq 1 minus", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = rnaseq_neg_lim)

# K562 histone marks 
H3K27ac <- DataTrack(range = file.path(cutntag_dir, "H3K27ac_merged.spikenorm.bw"), genome = "GRCh38",
                     name = "H3K27ac CUT&Tag", col.histogram = "#3CBC75FF", fill.histogram = "#3CBC75FF", ylim = h3k27cnt_lim)

H3K4me1 <- DataTrack(range = file.path(cutntag_dir, "H3K4me1_merged.spikenorm.bw"), genome = "GRCh38",
                     name = "H3K4me1 CUT&Tag", col.histogram = "#39558CFF", fill.histogram = "#39558CFF", ylim = h3k4cnt_lim) 

H3K4me2 <- DataTrack(range = file.path(cutntag_dir, "H3K4me2_merged.spikenorm.bw"), genome = "GRCh38",
                     name = "H3K4me2 CUT&Tag", col.histogram = "#287D8EFF", fill.histogram = "#287D8EFF", ylim = h3k4cnt_lim) 

H3K4me3 <- DataTrack(range = file.path(cutntag_dir, "H3K4me3_merged.spikenorm.bw"), genome = "GRCh38",
                     name = "H3K4me3 CUT&Tag", col.histogram = "#1F968BFF", fill.histogram = "#1F968BFF", ylim = h3k4cnt_lim)

IgG <- DataTrack(range = file.path(cutntag_dir, "IgG_merged.spikenorm.bw"), genome = "GRCh38",
                     name = "IgG CUT&Tag", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = h3k27cnt_lim) 

# K562 enhancer annotations
enhancers <- AnnotationTrack(range = file.path("human_data", "K562_enhancers_hg38.bed"), genome = "GRCh38",
                       name = "Enhancers", col = "black", col.line = "white", fill = "black")

# STRIPE-seq - TBPL1
cairo_pdf(file = file.path(gviz_dir, "STRIPE_TBPL1.pdf"), height = 8, width = 12)
plotTracks(list(axis.track, 
                K562_100ng_1_pos, 
                K562_100ng_1_neg, 
                K562_100ng_2_pos, 
                K562_100ng_2_neg,
                K562_100ng_3_pos, 
                K562_100ng_3_neg,
                K562_rnaseq_1_pos,
                K562_rnaseq_1_neg,
                genome.track), 
           chromosome = "6", from = 133953000, to = 133953500, 
           background.title = "white", 
           col.title = "black", 
           col.axis = "black", 
           type = "histogram", 
           baseline = 0, 
           col.baseline = "black")
dev.off()

# STRIPE-seq - MALAT1
cairo_pdf(file = file.path(gviz_dir, "STRIPE_MALAT1.pdf"), height = 8, width = 12)
plotTracks(list(axis.track, 
                K562_100ng_1_pos, 
                K562_100ng_1_neg, 
                K562_100ng_2_pos, 
                K562_100ng_2_neg,
                K562_100ng_3_pos, 
                K562_100ng_3_neg,
                K562_rnaseq_1_pos,
                K562_rnaseq_1_neg,
                genome.track), 
           chromosome = "11", from = 65498800, to = 65499300, 
           background.title = "white", 
           col.title = "black", 
           col.axis = "black", 
           type = "histogram", 
           baseline = 0, 
           col.baseline = "black")
dev.off()

# STRIPE-seq + CAGE - RARS2/ORC3
cairo_pdf(file = file.path(gviz_dir, "STRIPE_RARS2_ORC3.pdf"), height = 12, width = 12)
plotTracks(list(axis.track, 
                K562_100ng_1_pos, 
                K562_100ng_1_neg, 
                CAGE_10ug_1_pos,
                CAGE_10ug_1_neg,
                RAMPAGE_5ug_1_pos,
                RAMPAGE_5ug_1_neg,
                nanoCAGE_XL_7.5ug_1_pos,
                nanoCAGE_XL_7.5ug_1_neg,
                K562_rnaseq_1_pos,
                K562_rnaseq_1_neg,
                genome.track), 
           chromosome = "6", from = 87589800, to = 87590300, 
           background.title = "white", 
           col.title = "black", 
           col.axis = "black", 
           type = "histogram", 
           baseline = 0, 
           col.baseline = "black")
dev.off()

# STRIPE-seq + CAGE + CUT&Tag - LDLRAD3 enhancer
cairo_pdf(file = file.path(gviz_dir, "STRIPE_CnT_LDLRAD3_enhancer.pdf"), height = 12, width = 12)
plotTracks(list(axis.track, 
                enhancers,
                K562_100ng_1_pos, 
                K562_100ng_1_neg, 
                CAGE_10ug_1_pos,
                CAGE_10ug_1_neg,
                RAMPAGE_5ug_1_pos,
                RAMPAGE_5ug_1_neg,
                nanoCAGE_XL_7.5ug_1_pos,
                nanoCAGE_XL_7.5ug_1_neg,
                H3K4me1,
                H3K4me2,
                H3K4me3,
                H3K27ac,
                IgG,
                genome.track), 
           chromosome = "11", from = 36075800, to = 36141000, 
           background.title = "white", 
           col.title = "black", 
           col.axis = "black", 
           type = "histogram", 
           baseline = 0, 
           col.baseline = "black")
dev.off()

# STRIPE-seq + CAGE + CUT&Tag - chr2 intergenic enhancer
cairo_pdf(file = file.path(gviz_dir, "STRIPE_CnT_chr2_intergenic_enhancer.pdf"), height = 12, width = 12)
plotTracks(list(axis.track, 
                enhancers,
                K562_100ng_1_pos, 
                K562_100ng_1_neg, 
                CAGE_10ug_1_pos,
                CAGE_10ug_1_neg,
                RAMPAGE_5ug_1_pos,
                RAMPAGE_5ug_1_neg,
                nanoCAGE_XL_7.5ug_1_pos,
                nanoCAGE_XL_7.5ug_1_neg,
                H3K4me1,
                H3K4me2,
                H3K4me3,
                H3K27ac,
                IgG,
                genome.track), 
           chromosome = "2", from = 65805000, to = 65835600, 
           background.title = "white", 
           col.title = "black", 
           col.axis = "black", 
           type = "histogram", 
           baseline = 0, 
           col.baseline = "black")
dev.off()