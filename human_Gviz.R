library(GenomicFeatures)
library(Gviz)

baseDir <- "/Users/gzentner/Desktop/tsrexplorer/STRIPE-seq/"
setwd(baseDir)

if (!dir.exists(file.path(baseDir, "human_work/Gviz/"))){
    print("Creating directory 'human_work/Gviz' and changing working directory...")
    dir.create(file.path(baseDir, "human_work/Gviz/"))
    setwd(file.path(baseDir, "human_work/Gviz/"))
} else {
    print("Directory 'human_work/Gviz' already exisits, changing working directory...")
    setwd(file.path(baseDir, "human_work/Gviz/"))
}

# Create genomic axis track
axis.track <- GenomeAxisTrack(col="black", scale=0.1, col.range="black")

options(ucscChromosomeNames=FALSE)

# Create gene annotation track
txdb <- makeTxDbFromGFF(file.path(baseDir, "Homo_sapiens.GRCh38.98.gtf"))

genome.track <- GeneRegionTrack(txdb, genome="GRCh38", shape="arrow", names="Genes", col="black",
                                showId=TRUE, fill="black", trancriptAnnotation="gene_symbol", collapseTranscripts="meta")

# Create data tracks

# STRIPE-seq 
K562_100ng_1_pos <- DataTrack(range = file.path(baseDir, "human_work/bedgraphs/K562_100ng_1_+.bedgraph"), genome="GRCh38", 
                              name = "K562 100 ng 1 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 60))
K562_100ng_1_neg <- DataTrack(range = file.path(baseDir, "human_work/bedgraphs/K562_100ng_1_-.bedgraph"), genome="GRCh38", 
                              name = "K562 100 ng 1 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 60))
K562_100ng_2_pos <- DataTrack(range = file.path(baseDir, "human_work/bedgraphs/K562_100ng_2_+.bedgraph"), genome="GRCh38", 
                              name = "K562 100 ng 2 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 60))
K562_100ng_2_neg <- DataTrack(range = file.path(baseDir, "human_work/bedgraphs/K562_100ng_2_-.bedgraph"), genome="GRCh38", 
                              name = "K562 100 ng 2 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 60))
K562_100ng_3_pos <- DataTrack(range = file.path(baseDir, "human_work/bedgraphs/K562_100ng_3_+.bedgraph"), genome="GRCh38", 
                              name = "K562 100 ng 3 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 60))
K562_100ng_3_neg <- DataTrack(range = file.path(baseDir, "human_work/bedgraphs/K562_100ng_3_-.bedgraph"), genome="GRCh38", 
                              name = "K562 100 ng 3 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 60))

# CAGE
CAGE_10ug_1_pos <- DataTrack(range = file.path(baseDir, "human_work/bedgraphs/CAGE_10ug_1_+.bedgraph"), genome="GRCh38", 
                             name = "CAGE 10 ug 1 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 60))
CAGE_10ug_1_neg <- DataTrack(range = file.path(baseDir, "human_work/bedgraphs/CAGE_10ug_1_-.bedgraph"), genome="GRCh38", 
                             name = "CAGE 10 ug 1 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 60))
CAGE_10ug_2_pos <- DataTrack(range = file.path(baseDir, "human_work/bedgraphs/CAGE_10ug_2_+.bedgraph"), genome="GRCh38", 
                             name = "CAGE 10 ug 2 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 60))
CAGE_10ug_2_neg <- DataTrack(range = file.path(baseDir, "human_work/bedgraphs/CAGE_10ug_2_-.bedgraph"), genome="GRCh38", 
                             name = "CAGE 10 ug 2 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 60))

# RAMPAGE
RAMPAGE_5ug_1_pos <- DataTrack(range = file.path(baseDir, "human_work/bedgraphs/RAMPAGE_5ug_1_+.bedgraph"), genome="GRCh38", 
                             name = "RAMPAGE 10 ug 1 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 60))
RAMPAGE_5ug_1_neg <- DataTrack(range = file.path(baseDir, "human_work/bedgraphs/RAMPAGE_5ug_1_-.bedgraph"), genome="GRCh38", 
                             name = "RAMPAGE 10 ug 1 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 60))
RAMPAGE_5ug_2_pos <- DataTrack(range = file.path(baseDir, "human_work/bedgraphs/RAMPAGE_5ug_2_+.bedgraph"), genome="GRCh38", 
                               name = "RAMPAGE 10 ug 2 +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 60))
RAMPAGE_5ug_2_neg <- DataTrack(range = file.path(baseDir, "human_work/bedgraphs/RAMPAGE_5ug_2_-.bedgraph"), genome="GRCh38", 
                               name = "RAMPAGE 10 ug 2 -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 60))

# nanoCAGE-XL
nanoCAGE_XL_7.5ug_1_pos <- DataTrack(range = file.path(baseDir, "human_work/bedgraphs/nanoCAGE_XL_7.5ug_1_+.bedgraph"), genome="GRCh38", 
                               name = "nanoCAGE-XL 7.5 ug +", col.histogram = "black", fill.histogram = "black", ylim=c(0, 60))
nanoCAGE_XL_7.5ug_1_neg <- DataTrack(range = file.path(baseDir, "human_work/bedgraphs/nanoCAGE_XL_7.5ug_1_-.bedgraph"), genome="GRCh38", 
                               name = "nanoCAGE-XL 7.5 ug -", col.histogram = "black", fill.histogram = "black", ylim=c(0, 60))

# STRIPE-seq and CAGE (GAPDH region)
cairo_pdf(file = "tracks.pdf", height = 12, width = 12)
plotTracks(list(axis.track, 
                K562_100ng_1_pos, 
                K562_100ng_1_neg, 
                CAGE_10ug_1_pos, 
                CAGE_10ug_1_neg, 
                RAMPAGE_5ug_1_pos, 
                RAMPAGE_5ug_1_neg,
                nanoCAGE_XL_7.5ug_1_pos,
                nanoCAGE_XL_7.5ug_1_neg,
                genome.track), 
           chromosome = "12", from = 6484201, to = 6595900, 
           background.title = "white", 
           col.title = "black", 
           col.axis = "black", 
           type = "histogram", 
           baseline = 0, 
           col.baseline = "black")
dev.off()