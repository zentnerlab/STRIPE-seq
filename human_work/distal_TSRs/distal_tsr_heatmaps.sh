#!/bin/bash

#################################
### TSS signal at distal TSRs ###
#################################

# Compute TSS mapping signal around distal TSRs
computeMatrix \
reference-point \
-R distal_tsrs_all_reps.bed \
--referencePoint center \
-S tss_bigwigs/K562_100ng_1_plus.bw \
tss_bigwigs/K562_100ng_1_minus.bw \
tss_bigwigs/CAGE_10ug_1_plus.bw \
tss_bigwigs/CAGE_10ug_1_minus.bw \
tss_bigwigs/RAMPAGE_5ug_1_plus.bw \
tss_bigwigs/RAMPAGE_5ug_1_minus.bw \
tss_bigwigs/nanoCAGE_XL_7.5ug_1_plus.bw \
tss_bigwigs/nanoCAGE_XL_7.5ug_1_minus.bw \
-b 100 \
-a 100 \
-bs 1 \
-p 4 \
-o distal_tsr_tss_signal_matrix \
--verbose

# Heatmap TSS mapping signal around distal TSRs
plotHeatmap \
-m distal_tsr_tss_signal_matrix \
-o distal_tsr_tss_signal_matrix.pdf \
--colorMap Greys \
--missingDataColor '#ffffff' \
--sortUsingSamples 1 \
--sortRegions ascend \
--heatmapHeight 8 \
--zMin 0 \
--zMax 0.1 \
--whatToShow 'heatmap and colorbar'

#####################################
### CUT&Tag signal at distal TSRs ###
#####################################

# Compute histone mark signal around distal TSRs
computeMatrix \
reference-point \
-R distal_tsrs_all_reps.bed \
--referencePoint center \
-S cutntag_bigwigs/H3K4me1_merged.spikenorm.bw \
cutntag_bigwigs/H3K4me2_merged.spikenorm.bw \
cutntag_bigwigs/H3K4me3_merged.spikenorm.bw \
cutntag_bigwigs/H3K27ac_merged.spikenorm.bw \
cutntag_bigwigs/IgG_merged.spikenorm.bw \
-b 1000 \
-a 1000 \
-bs 25 \
-p 4 \
-o distal_tsr_histone_matrix_cutntag \
--verbose

# Heatmap histone mark signal around distal TSRs
plotHeatmap \
-m distal_tsr_histone_matrix_cutntag \
-o distal_tsr_histone_matrix_cutntag.pdf \
--colorMap Greys \
--sortUsingSamples 1 \
--heatmapHeight 8 \
--colorMap viridis \
--zMax 20 20 20 0.5 0.5 \
--whatToShow 'heatmap and colorbar' \
--outFileSortedRegions distal_tsr_histone_clusters.txt