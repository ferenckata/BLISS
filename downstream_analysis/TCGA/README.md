### 1. Download TCGA dataset

My latest information (04.05.2018) is to download from here: https://portal.gdc.cancer.gov/

### 2. Copy number variation

The open datasets are based on SNPs and already segmented. It means that the breakpoints are rather "breakregions" with the length of approx 1-10kb between consecutive regions with diffrent copy numbers. To get the breakregions, `breakpoint_finder.py` code is used that extracts the regions between consecutive regions with diffrent copy numbers from all datasts.

Than, `bedtools makewindws` and `bedtools coverage` are used to find the number of breakregions across datasets per 10kb, that is to find the recurrent breaklocations. 


### 3. Prepare BLISS output

For comparison, the resolution of the DSB distribution should be also in this range. that can be done like this:

```
for file in hg38_chrlength/*.bed;\
do chrn=$(echo $file | cut -d"/" -f2 | cut -d"_" -f1);\
echo $chrn;\
bedtools makewindows -g $file -w 10000 >w10kb_allhg38/$chrn"_w10kb.bed";\
done
```

The corrdinates should be lifted over to hg38 assembly. Here the local liftover tool is used with the corresponding chain file downloaded from [UCSC](http://hgdownload.cse.ucsc.edu/downloads.html). Note: before using liftover tool read [this](https://genome-store.ucsc.edu/)!

First, the chromosome names should get"chr" prefix to match with the UCSC format. Than they can be overlifted.
```
for file in *UMI.bed;\
do name=$(echo $file | cut -d"_" -f1);\
echo $name;\
cat $file | awk '{print "chr" $0}' > $name"_ucscformat.bed";\
~/Applications/liftOver $name"_ucscformat.bed" ~/Documents/bliss/hg19ToHg38.over.chain $name"_hg38.bed" $name"_unmapped.bed";\
done
```

The file should be "expanded" and moved to a different folder to not to interfere with hg19 data
```
for file in *hg38.bed;\
do name=$(echo $file | cut -d"_" -f1);\
echo $name;\
cat $file | awk '{for(i=1;i<=$4;i++) print $0}' >$name"_hg38_exp.bed";\
done

mkdir hg38
mv *hg38_* hg38/
mv *hg38.bed hg38/
```

DSB coverage per chromosome is calculated
```
for file in ../*hg38_exp.bed;\
do name=$(echo $file | cut -d"/" -f4 | cut -d"_" -f1);\
echo $name;\
for chrc in ../w10kb_allhg38/*.bed;\
do chrnum=$(echo $chrc| cut -d"_" -f1);\
echo $chrnum;\
$bedpath"bedtools" coverage -counts -a $chrc -b $file > ../nz/$name"_"$chrnum"_10kb_cov.bed";\
done;\
done
```

Zero DSB counts are removed
```
for file in *.bed;\
do name=$(echo $file | cut -d"_" -f1-2);\
echo $name;\
ID=$(echo $name | cut -d"_" -f1);\
cat $file | awk -v ID="$ID" '{if($4>0){print $0 "\t" ID}}' > ../nz/$name"_nz_cov.bed";\
done

```

