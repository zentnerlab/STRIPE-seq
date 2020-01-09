#!/bin/bash

# Script to determine the spike-in normalization for CUT&Tag datasets used to assess 
# histone modification enrichment at distal TSRs in K562 STRIPE-seq data. Two replicate
# datasets for each antibody are merged and the number of reads mapping to the E. coli
# genome is counted for each. The spike-in factor is then calculated as 1000 / the
# number of reads mapped to the E. coli genome and output to the command line.

# Merge replicate EB1 E. coli alignment BAM files 

cd /N/dc2/scratch/gzentner/K562_histone/ecoli_results/aligned

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
SRR8383516_H3K4me3_Active_Motif_cat_39159_1.bam \ SRR8383517_H3K4me3_Active_Motif_cat_39159_2.bam

samtools merge \
-f \
-@ 8 \
IgG_merged.bam \
SRR8435051.bam \
SRR8435052.bam

# Count reads mapped to E. coli genome in each sample and calculate spike-in normalization factors

for i in *merged.bam
do
	NREADS=$( samtools view -c -F 260 $i )
	SCALE_NUMERATOR=1000
	echo -e ""
	echo -e $NREADS "reads in" $i "mapped to the E. coli genome"
	SPIKE=$( echo "scale=5 ; $SCALE_NUMERATOR / $NREADS" | bc )
	echo -e "Spike-in factor ( 1000 / E. coli mapped reads ) for" $i "=" $SPIKE
done
