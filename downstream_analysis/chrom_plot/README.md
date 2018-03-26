## Chromosome-wide plotting of DSBs

1. Creating sliding windows using bedtools to reduce the size of the data

`bedtools makewindows -g chr_length.tsv -w 1000 -s 500 > chr_sw.bed`

This case the window size is 1kb and the overlap is 500bp. So the result is:

```
chr10	0	1000
chr10	500	1500
chr10	1000	2000
chr10	1500	2500
chr10	2000	3000
chr10	2500	3500
chr10	3000	4000
chr10	3500	4500
chr10	4000	5000
chr10	4500	5500
...
```

Note: Depending on the `chr_lenght.tsv` input it can cover one or more chromosome. I create this per chromosome, so the input has only one line: `chr10	135534747`.

2. "Expanding" bedgraph file 
The output of the current mapping pipeline is a bedgraph file.
Bedtools coverage does not take the forth column into account, only the number of lines.
The following command will repeat the lines as many times as the count is in the 4th column.

`cat $file | awk '{for(i=1;i<=$4;i++) print $0}' >$ID"_exp.bed" `

so the result will look like this:

```
chr1	115748	115749	1
chr1	162488	162489	1
chr1	172554	172555	2
chr1	172554	172555	2
chr1	534242	534243	1
chr1	542170	542171	1
chr1	542175	542176	1
chr1	564879	564880	1
chr1	565853	565854	1
chr1	566032	566033	1
...
```

3. Calculating the coverage for each window
I use `-counts` flag to have only the number of overlap reported.

`bedtools coverage -counts -a chr_sw.bed -b $ID"_exp.bed" > $ID"_coverage.bed"`

This can be for for each chromosome separately. And the result is:

```
chr10	0	1000	0
chr10	500	1500	0
chr10	1000	2000	0
chr10	1500	2500	0
chr10	2000	3000	0
chr10	2500	3500	0
chr10	3000	4000	0
chr10	3500	4500	0
chr10	4000	5000	0
chr10	4500	5500	0
...
```

Note: You get the following warning if your input bed file doesn't follow the same naming convention, ie. if chromosome 1 is labelled as "1" instead of "chr1".

***** WARNING: File chr10_w1000_s500.bed has inconsistent naming convention for record

To fix this, you can run this command on your input file (it will replace the expanded file):
```
for file in *exp.bed;\
do name=$(echo $file | cut -d"_" -f1-2);\
echo $name;\
cat $file | awk '{print "chr"$0}' >$name"_mcn.bed";\
rm $name"_exp.bed";\
mv $name"_mcn.bed" $name"_exp.bed";\
done
```

4. Removing regions with zero counts
this command filters out all lines with 0 counts and adds the ID at the end of the line
(that will be useful when plotting different datasets):

```
mkdir -p ../nz
for file in *.bed;\
do name=$(echo $file | cut -d"_" -f1-3);\
echo $name;\
cat $file | awk -v name="$name" '{if($4>0){print $0 "\t" name}}' > ../nz/$name"_nz_cov.bed" ;\
done
```


And the result is:

```
chr10	84500	85500	1 ID1
chr10	85000	86000	1 ID1
chr10	93500	94500	2 ID1
chr10	94000	95000	2 ID1
chr10	94500	95500	2 ID1
chr10	95000	96000	2 ID1
chr10	102000	103000	2 ID1
chr10	102500	103500	2 ID1
chr10	104000	105000	3 ID1
chr10	104500	105500	3 ID1
...
```

Note: the code expects the file naming convention such as chrID_libraryID_experimentID_coverage.bed
The chrID was placed there in the previous step. If you don't have both library ID and experiment ID, run this code instead:
```
mkdir -p ../nz
for file in *.bed;\
do name=$(echo $file | cut -d"_" -f1-2);\
echo $name;\
cat $file | awk -v name="$name" '{if($4>0){print $0 "\t" name}}' > ../nz/$name"_nz_cov.bed" ;\
done
```

5. Plotting with chrom_plotter function (part of BLISS package) that is an R script based on [TitanCNA](http://bioconductor.org/packages/release/bioc/vignettes/TitanCNA/inst/doc/TitanCNA.pdf)

The input is the path to the files (one per chromosome) with the reported count per region. The resolution depends on the previous steps.
