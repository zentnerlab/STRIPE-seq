library(GenomicFeatures)
library(Gviz)

baseDir <- "/Users/gzentner/Desktop/tsrexplorer/STRIPE-seq/"
setwd(baseDir)

if (!dir.exists(file.path(baseDir, "yeast_work/Gviz/"))){
    print("Creating directory 'yeast_work/Gviz' and changing working directory...")
    dir.create(file.path(baseDir, "yeast_work/Gviz/"))
    setwd(file.path(baseDir, "yeast_work/Gviz/"))
} else {
    print("Directory 'yeast_work/Gviz' already exisits, changing working directory...")
    setwd(file.path(baseDir, "yeast_work/Gviz/"))
}

# Create genomic axis track
axis.track <- GenomeAxisTrack(col="black", scale=0.1, col.range="black")

options(ucscChromosomeNames=FALSE)

# Create gene annotation track
txdb <- makeTxDbFromGFF(annotation)

genome.track <- GeneRegionTrack(txdb,genome="sacCer3", shape="arrow", names="Genes", col="black",
                                showId=TRUE, fill="black", trancriptAnnotation="gene_symbol", collapseTranscripts="meta")

# Create data tracks

# STRIPE-seq 50 ng
S288C_50ng_1_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_50ng_1_+.bedgraph"), genome="sacCer3", 
                            name="50 ng 1 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_50ng_1_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_50ng_1_-.bedgraph"), genome="sacCer3", 
                            name="50 ng 1 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_50ng_2_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_50ng_2_+.bedgraph"), genome="sacCer3", 
                              name="50 ng 2 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_50ng_2_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_50ng_2_-.bedgraph"), genome="sacCer3", 
                              name="50 ng 2 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_50ng_3_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_50ng_3_+.bedgraph"), genome="sacCer3", 
                              name="50 ng 3 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_50ng_3_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_50ng_3_-.bedgraph"), genome="sacCer3", 
                              name="50 ng 3 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))

# STRIPE-seq 100 ng
S288C_100ng_1_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_100ng_1_+.bedgraph"), genome="sacCer3", 
                              name="100 ng 1 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_100ng_1_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_100ng_1_-.bedgraph"), genome="sacCer3", 
                              name="100 ng 1 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_100ng_2_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_100ng_2_+.bedgraph"), genome="sacCer3", 
                               name="100 ng 2 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_100ng_2_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_100ng_2_-.bedgraph"), genome="sacCer3", 
                               name="100 ng 2 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_100ng_3_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_100ng_3_+.bedgraph"), genome="sacCer3", 
                               name="100 ng 3 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_100ng_3_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_100ng_3_-.bedgraph"), genome="sacCer3", 
                               name="100 ng 3 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))

# STRIPE-seq 250 ng
S288C_250ng_1_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_250ng_1_+.bedgraph"), genome="sacCer3", 
                              name="250 ng 1 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_250ng_1_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_250ng_1_-.bedgraph"), genome="sacCer3", 
                              name="250 ng 1 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_250ng_2_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_250ng_2_+.bedgraph"), genome="sacCer3", 
                               name="250 ng 2 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_250ng_2_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_250ng_2_-.bedgraph"), genome="sacCer3", 
                               name="250 ng 2 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_250ng_3_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_250ng_3_+.bedgraph"), genome="sacCer3", 
                               name="250 ng 3 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_250ng_3_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_250ng_3_-.bedgraph"), genome="sacCer3", 
                               name="250 ng 3 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))

# STRIPE-seq diamide
S288C_diamide_100ng_1_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_diamide_100ng_1_+.bedgraph"), genome="sacCer3", 
                                       name="Diamide 1 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_diamide_100ng_1_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_diamide_100ng_1_-.bedgraph"), genome="sacCer3", 
                                       name="Diamide 1 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_diamide_100ng_2_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_diamide_100ng_2_+.bedgraph"), genome="sacCer3", 
                                       name="Diamide 2 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_diamide_100ng_2_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_diamide_100ng_2_-.bedgraph"), genome="sacCer3", 
                                       name="Diamide 2 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_diamide_100ng_3_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_diamide_100ng_3_+.bedgraph"), genome="sacCer3", 
                                       name="Diamide 3 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
S288C_diamide_100ng_3_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/S288C_diamide_100ng_3_-.bedgraph"), genome="sacCer3", 
                                       name="Diamide 3 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))

# SLIC-CAGE
SLIC_CAGE_100ng_1_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/SLIC_CAGE_100ng_1_+.bedgraph"), genome="sacCer3", 
                               name="SLIC_CAGE_100ng_1_+", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
SLIC_CAGE_100ng_1_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/SLIC_CAGE_100ng_1_-.bedgraph"), genome="sacCer3", 
                               name="SLIC_CAGE_100ng_1_-", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
SLIC_CAGE_100ng_2_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/SLIC_CAGE_100ng_2_+.bedgraph"), genome="sacCer3", 
                               name="SLIC_CAGE_100ng_2_+", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
SLIC_CAGE_100ng_2_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/SLIC_CAGE_100ng_2_-.bedgraph"), genome="sacCer3", 
                               name="SLIC_CAGE_100ng_2_-", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))

# 500 ng nanoCAGE
nanoCAGE_500ng_1_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/nanoCAGE_500ng_1_+.bedgraph"), genome="sacCer3", 
                               name="nanoCAGE_500ng_1_+", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
nanoCAGE_500ng_1_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/nanoCAGE_500ng_1_-.bedgraph"), genome="sacCer3", 
                                  name="nanoCAGE_500ng_1_-", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
nanoCAGE_500ng_2_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/nanoCAGE_500ng_2_+.bedgraph"), genome="sacCer3", 
                                  name="nanoCAGE_500ng_2_+", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
nanoCAGE_500ng_2_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/nanoCAGE_500ng_2_-.bedgraph"), genome="sacCer3", 
                                  name="nanoCAGE_500ng_2_-", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))

# 25 ng nanoCAGE
nanoCAGE_25ng_1_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/nanoCAGE_25ng_1_+.bedgraph"), genome="sacCer3", 
                                  name="nanoCAGE_25ng_1_+", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
nanoCAGE_25ng_1_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/nanoCAGE_25ng_1_-.bedgraph"), genome="sacCer3", 
                                  name="nanoCAGE_25ng_1_-", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
nanoCAGE_25ng_2_pos <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/nanoCAGE_25ng_2_+.bedgraph"), genome="sacCer3", 
                                  name="nanoCAGE_25ng_2_+", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))
nanoCAGE_25ng_2_neg <- DataTrack(range = file.path(baseDir, "yeast_work/bedgraphs/nanoCAGE_25ng_2_-.bedgraph"), genome="sacCer3", 
                                  name="nanoCAGE_25ng_2_-", col.histogram = "black", fill.histogram = "black", ylim=c(0, 350))

# STRIPE-seq
cairo_pdf(file = "STRIPE_HHF1_HHT1.pdf", width = 12, height = 10)
plotTracks(list(axis.track, 
                S288C_50ng_1_pos,
                S288C_50ng_1_neg,
                S288C_100ng_1_pos, 
                S288C_100ng_1_neg, 
                S288C_250ng_1_pos,
                S288C_250ng_1_neg,
                genome.track), 
           chromosome = "II",from = 254000, to = 259000, 
           background.title = "white", 
           col.title = "black", 
           col.axis= "black", 
           type="histogram", 
           baseline = 0, 
           col.baseline = "black")
dev.off()

# STRIPE-seq and CAGE (HHF1/HHT1 region)
cairo_pdf(file = "STRIPE_CAGE_HHF1_HHT1.pdf", width = 12, height = 15)
plotTracks(list(axis.track, 
                S288C_100ng_1_pos, 
                S288C_100ng_1_neg, 
                SLIC_CAGE_100ng_1_pos,
                SLIC_CAGE_100ng_1_neg,
                nanoCAGE_500ng_1_pos,
                nanoCAGE_500ng_1_neg,
                nanoCAGE_25ng_1_pos,
                nanoCAGE_25ng_2_neg,
                genome.track), 
           chromosome = "II",from = 254000, to = 259000, 
           background.title = "white", 
           col.title = "black", 
           col.axis= "black", 
           type="histogram", 
           baseline = 0, 
           col.baseline = "black")
dev.off()

# YPD and diamide STRIPE-seq (JEN1/URA1 region)
cairo_pdf(file = "diamide_JEN1_URA1.pdf", width = 12, height = 12)
plotTracks(list(axis.track, 
                S288C_100ng_1_pos, 
                S288C_100ng_1_neg, 
                S288C_diamide_100ng_1_pos,
                S288C_diamide_100ng_1_neg,
                genome.track), 
           chromosome = "XI",from = 18000, to = 30000, 
           background.title = "white", 
           col.title = "black", 
           col.axis= "black", 
           type="histogram", 
           baseline = 0, 
           col.baseline = "black")
dev.off()