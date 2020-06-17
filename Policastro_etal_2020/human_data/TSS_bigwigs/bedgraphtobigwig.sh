#!/bin/bash

BGDIR='/Users/gzentner/Desktop/tsrexplorer/STRIPE-seq/human_work/bedgraphs'

CHROMSIZES='/Users/gzentner/Desktop/tsrexplorer/STRIPE-seq/human_data/TSS_bigwigs/hg38.nochr.sizes'

OUTDIR='/Users/gzentner/Desktop/tsrexplorer/STRIPE-seq/human_data/TSS_bigwigs'

cd $BGDIR

mkdir -p tss_bigwigs

for i in *bedgraph
do
	sort -k1,1 -k2,2n $i -o $i
	bedGraphToBigWig $i $CHROMSIZES ${i/bedgraph/bw}
done

mv *bw $OUTDIR