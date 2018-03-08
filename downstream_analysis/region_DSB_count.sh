# Gencode database can be found here for hg19:
# https://www.gencodegenes.org/releases/19.html
# (note: make very sure you're using the same assembly that has been used for mapping!)

# Unzip the downloaded file
gunzip gencode.v19.annotation.gtf.gz

# Exon and gene coordinates are extracted and merged respectively using bedtools 2.27 merge.
cat gencode.v19.annotation.gtf | awk '{if($3=="exon"){print $1 "\t" $4 "\t" $5 "\t" $28 "\t0\t" $7}}' > exon_gencode19.bed
cat gencode.v19.annotation.gtf | awk '{if($3=="gene"){print $1 "\t" $4 "\t" $5 "\t" $10 "\t0\t" $7}}' > gene_gencode19.bed

# Sorting is always useful before using bedtools (even if you suspect that your file is already sorted).
sort -k1,1 -k2,2n exon_gencode19.bed > exon_srt_gencode19.bed
sort -k1,1 -k2,2n gene_gencode19.bed > gene_srt_gencode19.bed

# Merge with respect of strandedness and keeping the strand info.
bedtools merge -i exon_srt_gencode19.bed -s -c 6 -o distinct > exon_srt_m_gencode19.bed
bedtools merge -i gene_srt_gencode19.bed -s -c 6 -o distinct > gene_srt_m_gencode19.bed

# DSBs are counted in these exonic regions. DSB number in intron is calculated as [(DSB# in gene) - (DSB# in exon)].
# DSB number in intergenic region is calculated as [(total # DSB) - (DSB# in gene)].

# Bedtools intersect is used, -wa and -wb options to have both strand and DSB count reported.
bedtools intersect -wa -wb -a exon_srt_m_gencode19.bed -b $file > $id"_exon.bed"
bedtools intersect -wa -wb -a gene_srt_m_gencode19.bed -b $file > $id"_gene.bed"

# From this point R is used, because it’s good for managing data frames and plotting.
DSBcountR.r $path