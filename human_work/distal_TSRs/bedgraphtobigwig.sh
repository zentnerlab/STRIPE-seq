#!/bin/bash

cd /Users/gzentner/Desktop/tsrexplorer/STRIPE-seq/human_work/bedgraphs

mkdir -p tss_bigwigs

for i in *bedgraph
do
	sort -k1,1 -k2,2n $i -o $i
	bedGraphToBigWig $i hg38.nochr.sizes ${i/bedgraph/bw}
done

mv *bw tss_bigwigs