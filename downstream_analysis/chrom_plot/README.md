## Chromosome-wide plotting of DSBs

1. Creating sliding windows using bedtools to reduce the size of the data
`bedtools makewindow -g chr_length.tsv -w 1000 -s 500 > chr_sw.bed`
This case the window size is 1kb and the overlap is 500bp

2. "Expanding" bedgraph file 
The output of the current mapping pipeline is a bedgraph file.
Bedtools coverage does not take the count into account, only the number of entries (lines).
The following command will repeat the lines as many times as the count is in the 4th column.
`cat $file | awk '{for(i=1;i<=$4;i++) print $0}' >$ID"_exp.bed" `

3. Calculating the coverage for each window, with `-counts` flag it will only report the number of overlap
`bedtools coverage -counts -a chr_sw.bed -b $ID"_exp.bed" > $ID"_coverage.bed"`
This can be for for each chromosome separately.

4. Plotting with chrom_plotter.R
This one is a custom-made R script based on [TitanCNA](http://bioconductor.org/packages/release/bioc/vignettes/TitanCNA/inst/doc/TitanCNA.pdf)
