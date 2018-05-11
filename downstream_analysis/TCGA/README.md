### 1. Download TCGA dataset

My latest information (04.05.2018) is to download from here: https://portal.gdc.cancer.gov/

### 2. Copy number variation

The open datasets are based on SNPs and already segmented. It means that the breakpoints are rather "breakregions" with the length of approx 1-10kb between consecutive regions with diffrent copy numbers. To get the breakregions, `breakpoint_finder.py` code is used that extracts the regions between consecutive regions with diffrent copy numbers from all datasts.

Than, `bedtools makewindws` and `bedtools coverage` are used to find the number of breakregions across datasets per 10kb, that is to find the recurrent breaklocations. 


### 3. Prepare BLISS output

### 4. Cookbook below

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

Arrange conveniently each dataset to a corresponding folder (ID) and merge files of different chromosomes from the same dataset to one file, while "expanding" it again for bedtools. Also removing the ID because this time it is inconveniet to have.

` for dr in ID*;do cat $dr"/"*"nz_cov.bed" |awk '{for(i=1;i<=$4;i++) {print $1 "\t" $2 "\t" $3 "\t" $4}}' > "all_"$dr".bed";done `

Use `breakpoint_finder.ipynb` to extract the breakregions from the downloaded TCGA samples. A file with the hg38 chromosome lengths is needed.

Calculate coverage for all TCGA file on all chromsomes. First, 10kb long windows are made of hg38.

` $bedpath"bedtools" makewindows -g hg38.txt -w 10000 > allchr_hg38_w10kb.bed `

If it starts with chr10 for whatever reason, do a sorting step.

` cat allchr_hg38_w10kb.bed | sort -k1,1 -k2,2n > hg38_srt_10kb.bed `

Than calculate the coverage with bedtools ...

` cat breakregions/* | $bedpath"bedtools" coverage -counts -a hg38_srt_10kb.bed -b stdin > alltcga_hg38_cov_10kb.bed `

... and finally remove zero counts.

` cat alltcga_hg38_cov_10kb.bed | awk '{if($4>0) print}' >nz_tcga_hg38_10kb_cov.bed `

Open `TCGA_breakpoint_enrichment.R` script and runt he first few command to have an overview on the data.

Recurrent breakpoint is identified here as a 10kb region where at least 0.5% of the samples have breakpoint (in this case it is 7). 

` cat nz_tcga_hg38_10kb_cov.bed | awk '{if($4>7) print}' > tcga_min7_hg38_10kb_cov.bed `

Create 10kb bins around the recurrent CN breakpoints.

```
cat tcga_min7_hg38_10kb_cov.bed | awk '{print $1 "\t" $2-10000 "\t" $2 "\t" $4}' > tcga_min7_hg38_10kb_m0-10_cov.bed
cat tcga_min7_hg38_10kb_cov.bed | awk '{print $1 "\t" $2-20000 "\t" $2-10000 "\t" $4}' > tcga_min7_hg38_10kb_m10-20_cov.bed
cat tcga_min7_hg38_10kb_cov.bed | awk '{print $1 "\t" $2-30000 "\t" $2-20000 "\t" $4}' > tcga_min7_hg38_10kb_m20-30_cov.bed
cat tcga_min7_hg38_10kb_cov.bed | awk '{print $1 "\t" $2-40000 "\t" $2-30000 "\t" $4}' > tcga_min7_hg38_10kb_m30-40_cov.bed
cat tcga_min7_hg38_10kb_cov.bed | awk '{print $1 "\t" $2-50000 "\t" $2-40000 "\t" $4}' > tcga_min7_hg38_10kb_m40-50_cov.bed
cat tcga_min7_hg38_10kb_cov.bed | awk '{print $1 "\t" $3+40000 "\t" $3+50000 "\t" $4}' > tcga_min7_hg38_10kb_p40-50_cov.bed
cat tcga_min7_hg38_10kb_cov.bed | awk '{print $1 "\t" $3+30000 "\t" $3+40000 "\t" $4}' > tcga_min7_hg38_10kb_p30-40_cov.bed
cat tcga_min7_hg38_10kb_cov.bed | awk '{print $1 "\t" $3+20000 "\t" $3+30000 "\t" $4}' > tcga_min7_hg38_10kb_p20-30_cov.bed
cat tcga_min7_hg38_10kb_cov.bed | awk '{print $1 "\t" $3+10000 "\t" $3+20000 "\t" $4}' > tcga_min7_hg38_10kb_p10-20_cov.bed
cat tcga_min7_hg38_10kb_cov.bed | awk '{print $1 "\t" $3 "\t" $3+10000 "\t" $4}' > tcga_min7_hg38_10kb_p0-10_cov.bed
```

Intersect these with DSBs.

```
for file in ../../hg38_deepseqcov/all*.bed;\
do name=$(echo $file | cut -d"/" -f4 | cut -d"." -f1); \
echo $name;\
for tcga in ../*.bed;\
do tcgn=$(echo $tcga | cut -d"_" -f5 | cut -d"." -f1);\
echo $tcgn;\
$bedpath"bedtools" coverage -counts -a $tcga -b $file > $name"_"$tcgn"_cov.bed";\
done;\
done
```

The `TCGA_breakpoint_enrichment.R` script does the plotting as well. It expects the data to be in different folders for different datasets.

To generate a random bed file use this command (here the same amount of entry is generated as many breakpoint were identified in TCGA):

` $bedpath"bedtools" random -l 10000 -n 30478 -g hg38_chrlength.tsv >random_points.bed `

And to the same as before:

```
cat random_points.bed | awk '{print $1 "\t" $2-10000 "\t" $2 "\t" $4}' > random_m0-10_cov.bed
cat random_points.bed | awk '{print $1 "\t" $2-20000 "\t" $2-10000 "\t" $4}' > random_m10-20_cov.bed
cat random_points.bed | awk '{print $1 "\t" $2-30000 "\t" $2-20000 "\t" $4}' > random_m20-30_cov.bed
cat random_points.bed | awk '{print $1 "\t" $2-40000 "\t" $2-30000 "\t" $4}' > random_m30-40_cov.bed
cat random_points.bed | awk '{print $1 "\t" $2-50000 "\t" $2-40000 "\t" $4}' > random_m40-50_cov.bed
cat random_points.bed | awk '{print $1 "\t" $3+40000 "\t" $3+50000 "\t" $4}' > random_p40-50_cov.bed    
cat random_points.bed | awk '{print $1 "\t" $3+30000 "\t" $3+40000 "\t" $4}' > random_p30-40_cov.bed
cat random_points.bed | awk '{print $1 "\t" $3+20000 "\t" $3+30000 "\t" $4}' > random_p20-30_cov.bed
cat random_points.bed | awk '{print $1 "\t" $3+10000 "\t" $3+20000 "\t" $4}' > random_p10-20_cov.bed
cat random_points.bed | awk '{print $1 "\t" $3 "\t" $3+10000 "\t" $4}' > random_p0-10_cov.bed
```

Note: sometimes it generates invalid records since the breakpoint or the random bed can be close to the end of the genome.
You can remove them with this line (save to another folder):
` for file in random_*;do echo $file; cat $file | awk '{if($2>=0 && $3>=0) print}' > forplot/$file;done `


Intersect these with DSBs.

```
 for file in hg38_deepseqcov/all*.bed;\
 do name=$(echo $file | cut -d"/" -f5 | cut -d"." -f1);\
 echo $name;\
 for rnd in ../*.bed;\
 do rndn=$(echo $rnd | cut -d"_" -f2);\
 echo $rndn;\
 $bedpath"bedtools" coverage -counts -a $rnd -b $file > $name"_"$rndn"_cov.bed";\
 done;\
 done
```
