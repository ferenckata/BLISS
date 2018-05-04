### 1. Download TCGA dataset

My latest information (04.05.2018) is to download from here: https://portal.gdc.cancer.gov/

### 2. Copy number variation

The open datasets are based on SNPs and already segmented. It means that the breakpoints are rather "breakregions" with the length of approx 1-10kb between consecutive regions with diffrent copy numbers.

For comparison, the resolution of the DSB distribution should be also in this range. that can be done like this:
```
for file in chrlength/*.bed;\
do chrn=$(echo $file | cut -d"/" -f2 | cut -d"_" -f1);\
echo $chrn;\
bedtools makewindows -g $file -w 10000 >w1kb_allhg19/$chrn"_w10kb.bed";\
done 


```

