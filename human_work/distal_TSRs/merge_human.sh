#!/bin/bash

# Script to generate spike-in-normalized human CUT&Tag bigwig files to assess histone modification at distal TSRs in K562 STRIPE-seq data. Two replicate datasets for each antibody are merged and indexed with SAMTools. The spike-in factors determined by merge_ecoli.sh are then used in the generation of normalized bigwig files with deepTools bamCoverage.

# Merge replicate GRCh38 alignment BAM files

cd /N/dc2/scratch/gzentner/K562_histone/human_results/aligned

samtools merge \
-f \
-@ 8 \
H3K27ac_merged.bam \
SRR8383507_Millipore_cat_MABE647_1.bam \
SRR8383508_Millipore_cat_MABE647_2.bam 

samtools merge \
-f \
-@ 8 \
H3K4me1_merged.bam \
SRR8383512_H3K4me1_Abcam_cat_ab8895_1.bam \
SRR8383513_H3K4me1_Abcam_cat_ab8895_2.bam

samtools merge \
-f \
-@ 8 \
H3K4me2_merged.bam \
SRR8383514_H3K4me2_Upstate_cat_07_030_1.bam \
SRR8383515_H3K4me2_Upstate_cat_07_030_2.bam

samtools merge \
-f \
-@ 8 \
H3K4me3_merged.bam \
SRR8383516_H3K4me3_Active_Motif_cat_39159_1.bam \
SRR8383517_H3K4me3_Active_Motif_cat_39159_2.bam

samtools merge \
-f \
-@ 8 \
IgG_merged.bam \
SRR8435051.bam \
SRR8435052.bam

for i in *merged.bam 
do 
	samtools index $i 
done 

# Generate spike-in-normalized bigWig files

bamCoverage \
-b H3K27ac_merged.bam \
-of bigwig \
-o H3K27ac_merged.spikenorm.bw \
-bs 25 \
--scaleFactor .00981 \
-p 8 \
-e

bamCoverage \
-b H3K4me1_merged.bam \
-of bigwig \
-o H3K4me1_merged.spikenorm.bw \
-bs 25 \
--scaleFactor .51813 \
-p 8 \
-e

bamCoverage \
-b H3K4me2_merged.bam \
-of bigwig \
-o H3K4me2_merged.spikenorm.bw \
-bs 25 \
--scaleFactor .12680 \
-p 8 \
-e

bamCoverage \
-b H3K4me3_merged.bam \
-of bigwig \
-o H3K4me3_merged.spikenorm.bw \
-bs 25 \
--scaleFactor .05466 \
-p 8 \
-e

bamCoverage \
-b IgG_merged.bam \
-of bigwig \
-o IgG_merged.spikenorm.bw \
-bs 25 \
--scaleFactor .00354 \
-p 8 \
-e