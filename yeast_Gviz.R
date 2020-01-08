library(tsrexplorer)
library(GenomicFeatures)
library(Gviz)
library(viridis)
library(scales)

# This script imports signal tracks and displays them at the specified regions of the yeast genome, listed at the end of the script.

gviz_dir <- file.path("yeast_work", "Gviz")
bedgraph_dir <- file.path("yeast_work", "bedgraphs")
rnaseq_dir <- file.path("yeast_work", "bedgraphs")

if (!dir.exists(gviz_dir)) {
    message("Creating directory 'yeast_work/Gviz'...")
    dir.create(gviz_dir)
} else {
    message("Directory 'yeast_work/Gviz' already exists.")
}

# Create genomic axis track
axis.track <- GenomeAxisTrack(col = "black", scale = 0.1, col.range = "black")

options(ucscChromosomeNames = FALSE)

# Create gene annotation track
txdb <- makeTxDbFromGFF(annotation)

genome.track <- GeneRegionTrack(txdb, genome = "sacCer3", shape = "arrow", names = "Genes", col = "black",
                                showId = TRUE, fill = "black", trancriptAnnotation = "gene_symbol", collapseTranscripts = "meta")

# Create data tracks

# Get colors from the viridis palette for the number of tracks to be plotted
show_col(viridis_pal()(20))

# Set y-axis limits
stripe_pos_lim <- c(0,250)
stripe_neg_lim <- c(0,250)

slic_pos_lim <- c(0,300)
slic_neg_lim <- c(0,300)

nano500_pos_lim <- c(0,3000)
nano500_neg_lim <- c(0,3000)

nano25_pos_lim <- c(0,3500)
nano25_neg_lim <- c(0,3500)

rnaseq_pos_lim <- c(0,150)
rnaseq_neg_lim <- c(0,250)

riboseq_pos_lim <- c(0,6)
riboseq_neg_lim <- c(0,5)

# STRIPE-seq 50 ng
S288C_50ng_1_pos <- DataTrack(range = file.path(bedgraph_dir, "S288C_50ng_1_plus.bedgraph"), genome = "sacCer3", 
                            name = "50 ng 1 plus", col.histogram = "#440154FF", fill.histogram = "#440154FF", ylim = stripe_pos_lim)
S288C_50ng_1_neg <- DataTrack(range = file.path(bedgraph_dir, "S288C_50ng_1_minus.bedgraph"), genome = "sacCer3", 
                            name = "50 ng 1 minus", col.histogram = "#440154FF", fill.histogram = "#440154FF", ylim = stripe_neg_lim)

# STRIPE-seq 100 ng
S288C_100ng_1_pos <- DataTrack(range = file.path(bedgraph_dir, "S288C_100ng_1_plus.bedgraph"), genome = "sacCer3", 
                              name = "100 ng 1 plus", col.histogram = "#31688EFF", fill.histogram = "#31688EFF", ylim = stripe_pos_lim)
S288C_100ng_1_neg <- DataTrack(range = file.path(bedgraph_dir, "S288C_100ng_1_minus.bedgraph"), genome = "sacCer3", 
                              name = "100 ng 1 minus", col.histogram = "#31688EFF", fill.histogram = "#31688EFF", ylim = stripe_neg_lim)

# STRIPE-seq 250 ng
S288C_250ng_1_pos <- DataTrack(range = file.path(bedgraph_dir, "S288C_250ng_1_plus.bedgraph"), genome = "sacCer3", 
                              name = "250 ng 1 plus", col.histogram = "#35B779FF", fill.histogram = "#35B779FF", ylim = stripe_pos_lim)
S288C_250ng_1_neg <- DataTrack(range = file.path(bedgraph_dir, "S288C_250ng_1_minus.bedgraph"), genome = "sacCer3", 
                              name = "250 ng 1 minus", col.histogram = "#35B779FF", fill.histogram = "#35B779FF", ylim = stripe_neg_lim)

# STRIPE-seq diamide
S288C_diamide_100ng_1_pos <- DataTrack(range = file.path(bedgraph_dir, "S288C_diamide_100ng_1_plus.bedgraph"), genome = "sacCer3", 
                                       name = "Diamide 1 plus", col.histogram = "#414487FF", fill.histogram = "#414487FF", ylim = stripe_pos_lim)
S288C_diamide_100ng_1_neg <- DataTrack(range = file.path(bedgraph_dir, "S288C_diamide_100ng_1_minus.bedgraph"), genome = "sacCer3", 
                                       name = "Diamide 1 minus", col.histogram = "#414487FF", fill.histogram = "#414487FF", ylim = stripe_neg_lim)

# SLIC-CAGE
SLIC_CAGE_100ng_1_pos <- DataTrack(range = file.path(bedgraph_dir, "SLIC_CAGE_100ng_1_plus.bedgraph"), genome = "sacCer3", 
                               name = "SLIC-CAGE 100 ng 1 plus", col.histogram = "#56C667FF", fill.histogram = "#56C667FF", ylim = slic_pos_lim)
SLIC_CAGE_100ng_1_neg <- DataTrack(range = file.path(bedgraph_dir, "SLIC_CAGE_100ng_1_minus.bedgraph"), genome="sacCer3", 
                               name = "SLIC-CAGE 100 ng 1 minus", col.histogram = "#56C667FF", fill.histogram = "#56C667FF", ylim = slic_neg_lim)

# 500 ng nanoCAGE
nanoCAGE_500ng_1_pos <- DataTrack(range = file.path(bedgraph_dir, "nanoCAGE_500ng_1_plus.bedgraph"), genome = "sacCer3", 
                               name = "nanoCAGE 500 ng 1 plus", col.histogram = "#94D840FF", fill.histogram = "#94D840FF", ylim = nano500_pos_lim)
nanoCAGE_500ng_1_neg <- DataTrack(range = file.path(bedgraph_dir, "nanoCAGE_500ng_1_minus.bedgraph"), genome = "sacCer3", 
                                  name = "nanoCAGE 500 ng 1 minus", col.histogram = "#94D840FF", fill.histogram = "#94D840FF", ylim = nano500_neg_lim)

# 25 ng nanoCAGE
nanoCAGE_25ng_1_pos <- DataTrack(range = file.path(bedgraph_dir, "nanoCAGE_25ng_1_plus.bedgraph"), genome = "sacCer3", 
                                  name = "nanoCAGE 25 ng 1 plus", col.histogram = "#DCE318FF", fill.histogram = "#DCE318FF", ylim = nano25_pos_lim)
nanoCAGE_25ng_1_neg <- DataTrack(range = file.path(bedgraph_dir, "nanoCAGE_25ng_1_minus.bedgraph"), genome = "sacCer3", 
                                  name = "nanoCAGE 25 ng 1 minus", col.histogram = "#DCE318FF", fill.histogram = "#DCE318FF", ylim = nano25_neg_lim)

# YPD RNA-seq
ypd_rnaseq_1_pos <- DataTrack(range = file.path(rnaseq_dir, "RNASEQ001_S288C_untreated_r1.CPM.bs1.smooth25.plus.bw"), genome = "sacCer3", 
                              name = "YPD RNA-seq 1 plus", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = rnaseq_pos_lim)
ypd_rnaseq_1_neg <- DataTrack(range = file.path(rnaseq_dir, "RNASEQ001_S288C_untreated_r1.CPM.bs1.smooth25.minus.bw"), genome = "sacCer3", 
                              name = "YPD RNA-seq 1 minus", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = rnaseq_neg_lim)

# Diamide RNA-seq
diamide_rnaseq_1_pos <- DataTrack(range = file.path(rnaseq_dir, "RNASEQ009_S288C_diamide_r2.CPM.bs1.smooth25.plus.bw"), genome = "sacCer3", 
                              name = "YPD RNA-seq 1 plus", col.histogram = "#22A884FF", fill.histogram = "#22A884FF", ylim = rnaseq_pos_lim)
diamide_rnaseq_1_neg <- DataTrack(range = file.path(rnaseq_dir, "RNASEQ009_S288C_diamide_r2.CPM.bs1.smooth25.minus.bw"), genome = "sacCer3", 
                              name = "YPD RNA-seq 1 minus", col.histogram = "#22A884FF", fill.histogram = "#22A884FF", ylim = rnaseq_neg_lim)

# Ribosome profiling
S288C_riboseq_1_pos <- DataTrack(range = file.path(bedgraph_dir, "GSM1949550_rep1_positive.wig"), genome = "sacCer3",
                                 name = "S288C Ribo-seq 1 plus", col.histogram = "#3CBC75FF", fill.histogram = "#3CBC75FF", 
                                 ylim = riboseq_pos_lim)
S288C_riboseq_1_neg <- DataTrack(range = file.path(bedgraph_dir, "GSM1949550_rep1_negative.wig"), genome = "sacCer3", 
                                 name = "S288C Ribo-seq 1 minus", col.histogram = "#3CBC75FF", fill.histogram = "#3CBC75FF", 
                                 ylim = riboseq_neg_lim)
S288C_riboseq_2_pos <- DataTrack(range = file.path(bedgraph_dir, "GSM1949551_rep2_positive.wig"), genome = "sacCer3",
                                 name = "S288C Ribo-seq 2 plus", col.histogram = "#3CBC75FF", fill.histogram = "#3CBC75FF", 
                                 ylim = riboseq_pos_lim)
S288C_riboseq_2_neg <- DataTrack(range = file.path(bedgraph_dir, "GSM1949551_rep2_negative.wig"), genome = "sacCer3", 
                                 name = "S288C Ribo-seq 2 minus", col.histogram = "#3CBC75FF", fill.histogram = "#3CBC75FF", 
                                 ylim = riboseq_neg_lim)

# STRIPE-seq - GCN4
cairo_pdf(file = file.path(gviz_dir, "STRIPE_GCN4.pdf"), width = 12, height = 10)
plotTracks(list(axis.track, 
                S288C_50ng_1_pos,
                S288C_50ng_1_neg,
                S288C_100ng_1_pos, 
                S288C_100ng_1_neg, 
                S288C_250ng_1_pos,
                S288C_250ng_1_neg,
                ypd_rnaseq_1_pos,
                ypd_rnaseq_1_neg,
                genome.track), 
           chromosome = "V",from = 138000, to = 141000, 
           background.title = "white", 
           col.title = "black", 
           col.axis= "black", 
           type="histogram", 
           baseline = 0, 
           col.baseline = "black"
)
dev.off()

# STRIPE-seq - AIM39
cairo_pdf(file = file.path(gviz_dir, "STRIPE_AIM39.pdf"), width = 12, height = 15)
plotTracks(list(axis.track, 
                S288C_50ng_1_pos,
                S288C_50ng_1_neg,
                S288C_100ng_1_pos, 
                S288C_100ng_1_neg, 
                S288C_250ng_1_pos,
                S288C_250ng_1_neg,
                ypd_rnaseq_1_pos,
                ypd_rnaseq_1_neg,
                genome.track), 
           chromosome = "XV", from = 230000, to = 232000, 
           background.title = "white", 
           col.title = "black", 
           col.axis = "black", 
           type = "histogram", 
           baseline = 0, 
           col.baseline = "black"
)
dev.off()

# STRIPE-seq + Ribo-seq - AIM39
highlight.track <- HighlightTrack(trackList = list(S288C_100ng_1_pos, 
                                                   S288C_100ng_1_neg,
                                                   ypd_rnaseq_1_pos,
                                                   ypd_rnaseq_1_neg,
                                                   S288C_riboseq_1_pos,
                                                   S288C_riboseq_1_neg,
                                                   S288C_riboseq_2_pos,
                                                   S288C_riboseq_2_neg),
                                  chromosome = "XV", start = 230100, end = 230200, 
                                  inBackground = FALSE, 
                                  col = "black", fill = "lightblue", alpha = 0.5)

cairo_pdf(file = file.path(gviz_dir, "STRIPE_AIM39_riboseq.pdf"), width = 12, height = 12)
plotTracks(list(axis.track,
                highlight.track,
                genome.track), 
           chromosome = "XV",from = 230000, to = 232000, 
           background.title = "white", 
           col.title = "black", 
           col.axis = "black", 
           type ="histogram", 
           baseline = 0, 
           col.baseline = "black"
)
dev.off()

# STRIPE-seq + CAGE - YGR250C
cairo_pdf(file = file.path(gviz_dir, "STRIPE_CAGE_RIE1.pdf"), width = 12, height = 16)
plotTracks(list(axis.track, 
                S288C_50ng_1_pos,
                S288C_50ng_1_neg,
                S288C_100ng_1_pos, 
                S288C_100ng_1_neg, 
                S288C_250ng_1_pos,
                S288C_250ng_1_neg,
                SLIC_CAGE_100ng_1_pos,
                SLIC_CAGE_100ng_1_neg,
                nanoCAGE_500ng_1_pos,
                nanoCAGE_500ng_1_neg,
                nanoCAGE_25ng_1_pos,
                nanoCAGE_25ng_1_neg,
                ypd_rnaseq_1_pos,
                ypd_rnaseq_1_neg,
                genome.track), 
           chromosome = "VII",from = 993400, to = 994000, 
           background.title = "white", 
           col.title = "black", 
           col.axis= "black", 
           type="histogram", 
           baseline = 0, 
           col.baseline = "black"
)
dev.off()

# STRIPE-seq + CAGE - RPL17B
cairo_pdf(file = file.path(gviz_dir, "STRIPE_CAGE_RPL17B.pdf"), width = 12, height = 16)
plotTracks(list(axis.track, 
                S288C_50ng_1_pos,
                S288C_50ng_1_neg,
                S288C_100ng_1_pos, 
                S288C_100ng_1_neg, 
                S288C_250ng_1_pos,
                S288C_250ng_1_neg,
                SLIC_CAGE_100ng_1_pos,
                SLIC_CAGE_100ng_1_neg,
                nanoCAGE_500ng_1_pos,
                nanoCAGE_500ng_1_neg,
                nanoCAGE_25ng_1_pos,
                nanoCAGE_25ng_1_neg,
                ypd_rnaseq_1_pos,
                ypd_rnaseq_1_neg,
                genome.track), 
           chromosome = "X",from = 90500, to = 91000, 
           background.title = "white", 
           col.title = "black", 
           col.axis= "black", 
           type="histogram", 
           baseline = 0, 
           col.baseline = "black"
)
dev.off()

# YPD and diamide STRIPE-seq (HSP150/CIS3 region)
cairo_pdf(file = file.path(gviz_dir, "diamide_HSP150_CIS3.pdf"), width = 12, height = 12)
plotTracks(list(axis.track, 
                S288C_100ng_1_pos, 
                S288C_100ng_1_neg, 
                S288C_diamide_100ng_1_pos,
                S288C_diamide_100ng_1_neg,
                ypd_rnaseq_1_pos,
                ypd_rnaseq_1_neg,
                diamide_rnaseq_1_pos,
                diamide_rnaseq_1_neg,
                genome.track), 
           chromosome = "X", from = 120000, to = 123500, 
           background.title = "white", 
           col.title = "black", 
           col.axis = "black", 
           type = "histogram", 
           baseline = 0, 
           col.baseline = "black")
dev.off()

# YPD and diamide STRIPE-seq (RPS16B/RPL13A region)
cairo_pdf(file = file.path(gviz_dir, "diamide_RPS16B_RPL13A.pdf"), width = 12, height = 12)
plotTracks(list(axis.track, 
                S288C_100ng_1_pos, 
                S288C_100ng_1_neg, 
                S288C_diamide_100ng_1_pos,
                S288C_diamide_100ng_1_neg,
                ypd_rnaseq_1_pos,
                ypd_rnaseq_1_neg,
                diamide_rnaseq_1_pos,
                diamide_rnaseq_1_neg,
                genome.track), 
           chromosome = "IV",from = 306600, to = 309750, 
           background.title = "white", 
           col.title = "black", 
           col.axis= "black", 
           type="histogram", 
           baseline = 0, 
           col.baseline = "black")
dev.off()