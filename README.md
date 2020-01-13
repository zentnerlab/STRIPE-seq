# STRIPE-seq analysis
Policastro et al, 2020

## About
The scripts and data contained in this repository can be used to reproduce the majority of the images found in the STRIPE-seq publication. TSS bedgraphs, TSR beds, and RNA-seq bigwigs are available from GEO (GSE142524).

## Cloning Repository

Clone the ChIPseq automation repository. Navigate to a directory you would like to clone the repository to and enter  
`git clone https://github.com/zentnerlab/STRIPE-seq.git`

## Script overviews

1. *_TSRexplorer.R scripts load TSSs and TSRs detected in STRIPE-seq and CAGE datasets and perform the analyses described in the paper (threshold analysis, correlation, genomic distribution, etc). The scripts also analyze transcript abundances using STRIPE-seq and RNA-seq data. 

2. *_promoter_correlation.R scripts quantify signal in promoter windows and perform correlation analysis.

3. *_Gviz.R scripts load signal tracks for TSS mapping methods as well as RNA-seq, Ribo-seq (for yeast), and CUT&Tag (for human) and displays them at the specified genomic regions.

4. distal_TSR_analysis.R retrieves TSS-proximal and TSS-distal TSRs a specified distance from an annotated human TSS. It also calculates the overlap of proximal and distal TSRs with enhancer annotations from EnhancerAtlas 2.0 and performs statistical comparison of overlaps with a chi-square test.

5. distal_TSR_heatmaps.sh generates heatmaps of TSS and histone modification CUT&Tag signal around distal TSRs.

## A note regarding annotation files
For human analyses, the user must place a genome sequence file in FASTA format and a genome annotation file in GTF/GFF format in the human_data subfolder, as these are too large to be hosted on GitHub (we used Ensembl GRCh38.p13 for the analyses in the paper). The yeast genome sequence and annotation files (Ensembl V64.1.1) are built into TSRexploreR.