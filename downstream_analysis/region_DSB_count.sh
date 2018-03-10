# Gencode database can be found here for hg19:
# https://www.gencodegenes.org/releases/19.html
# (note: make very sure you're using the same assembly that has been used for mapping!)

# Unzip the downloaded file
gunzip gencode.v19.annotation.gtf.gz

# Exon, gene and TSS coordinates are extracted and merged respectively using bedtools 2.27 merge.
cat gencode.v19.annotation.gtf | awk '{if($3=="exon"){print $1 "\t" $4 "\t" $5 "\t" $28 "\t0\t" $7}}' > exon_gencode19.bed
cat gencode.v19.annotation.gtf | awk '{if($3=="gene"){print $1 "\t" $4 "\t" $5 "\t" $10 "\t0\t" $7}}' > gene_gencode19.bed
cat gencode.v19.annotation.gtf | awk '{if($20=="\"protein_coding\";" && $3=="transcript")\
{if($7=="-"){print $1 "\t" $5-2500 "\t" $5+2500 "\t" $18 "\t0\t" $7}else{print $1 "\t" $4-2500 "\t" $4+2500 "\t" $18 "\t0\t" $7}}}'\
> tss_gencode19.bed

# Sorting is always useful before using bedtools (even if you suspect that your file is already sorted).
sort -k1,1 -k2,2n exon_gencode19.bed > exon_srt_gencode19.bed
sort -k1,1 -k2,2n gene_gencode19.bed > gene_srt_gencode19.bed
sort -k1,1 -k2,2n tss_gencode19.bed > tss_srt_gencode19.bed

# Merge with respect of strandedness and keeping the gene info.
#### TODO: change to gene symbol!! ####
bedtools merge -i exon_srt_gencode19.bed -s -c 6 -o distinct > exon_srt_m_gencode19.bed
bedtools merge -i gene_srt_gencode19.bed -s -c 6 -o distinct > gene_srt_m_gencode19.bed

bedtools merge -i tss_srt_gencode19.bed -s -c 4 -o distinct > tss_srt_g_gencode19.bed

# DSBs are counted in these exonic regions. DSB number in intron is calculated as [(DSB# in gene) - (DSB# in exon)].
# DSB number in intergenic region is calculated as [(total # DSB) - (DSB# in gene)].

# Bedtools intersect could be used, -wa and -wb options to have both strand and DSB count reported.
# But bedtools coverage is better for further analysis
# coverage cannot take the DSB counts into account, so the files should be "expanded"
# meaning that the number of lines should correspond to the number of DSBs in the region
for file in *UMI.bed;\
do name=$(echo $file | cut -d"_" -f1);\
echo $name;\
cat $file | awk '{for(i=1;i<=$4;i++) print $0}' >$name"_exp.bed";\
done
# now coverage can be used
for file in *exp.bed;do name=$(echo $file | cut -d"_" -f1);echo $name;\
bedtools coverage -counts -a gene_srt_m_gencode19.bed -b $file >$name"_gene.bed";\
echo "exon";\
bedtools coverage -counts -a exon_srt_m_gencode19.bed -b $file >$name"_exon.bed";\
echo "tss" ;\
bedtools coverage -counts -a tss_srt_g_gencode19.bed -b $file >$name"_tss.bed" ;\
done

# From this point R is used, because it’s good for managing data frames and plotting.
DSBcountR.r $path1 $path2
