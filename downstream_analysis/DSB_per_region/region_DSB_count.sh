# ------------------- first part: region information from gencode database ----------------------

# Gencode database can be found here for hg19:
# https://www.gencodegenes.org/releases/19.html
# (note: make very sure you're using the same assembly that has been used for mapping!)

# create and step in the folder where you intend to keep the gene, exon and TSS regions
cd ~/Documents # or wherever you wish to keep your files
mkdir gene_db
genecodepath=~/Documents/gene_db # this can be different for you!!
cd $genecodepath
# Move here the downloaded file
mv ~/Downloads/gencode.v19.annotation.gtf.gz .
# Unzip the downloaded file
gunzip gencode.v19.annotation.gtf.gz

# ----- REGIONAL COUNTS --------

# Exon, gene and TSS coordinates are extracted and merged by strand respectively using bedtools 2.27 merge. The gene symbol is kept.
# Sorting is always useful before using bedtools (even if you suspect, but not 100% sure that your file is already sorted).
cat gencode.v19.annotation.gtf | awk '{if($20=="\"protein_coding\";" && $3=="exon"){print $1 "\t" $4 "\t" $5 "\t" $18 "\t0\t" $7}}' |\
sort -k1,1 -k2,2n > exon_srt_gencode19.bed
cat gencode.v19.annotation.gtf | awk '{if($20=="\"protein_coding\";" && $3=="gene"){print $1 "\t" $4 "\t" $5 "\t" $18 "\t0\t" $7}}' |\
sort -k1,1 -k2,2n > gene_srt_gencode19.bed
# for transcript start site, the strandedness is important:
# if the transcript is on the - strand, the higher coordinate is the TSS on the reference genome
cat gencode.v19.annotation.gtf | awk '{if($20=="\"protein_coding\";" && $3=="transcript")\
{if($7=="-"){print $1 "\t" $5-2500 "\t" $5+2500 "\t" $18 "\t0\t" $7}else{print $1 "\t" $4-2500 "\t" $4+2500 "\t" $18 "\t0\t" $7}}}'\
| sort -k1,1 -k2,2n > tss_srt_gencode19.bed

# Merge without respect of strandedness and keeping the gene symbol for the total DSB counts per region
$bedpath"bedtools" merge -i exon_srt_gencode19.bed -c 4 -o distinct > exon_srt_m_gencode19.bed
$bedpath"bedtools" merge -i gene_srt_gencode19.bed -c 4 -o distinct > gene_srt_m_gencode19.bed
$bedpath"bedtools" merge -i tss_srt_gencode19.bed -c 4 -o distinct > tss_srt_m_gencode19.bed

# Additionally the total length of genes, exons and TSS is caculated and stored for further analysis
echo "total" >> sumbp.tsv
# this number is the total number of non-N basepairs in the hg19
# downloaded from here:http://genomewiki.ucsc.edu/index.php/Hg19_Genome_size_statistics
echo "2897310462" >> sumbp.tsv
echo "gene" >> sumbp.tsv 
cat gene_srt_m_gencode19.bed | awk '{print $3-$2}' | paste -sd+ - | bc >> sumbp.tsv
echo "exon" >> sumbp.tsv 
cat exon_srt_m_gencode19.bed | awk '{print $3-$2}' | paste -sd+ - | bc >> sumbp.tsv
echo "TSS" >> sumbp.tsv 
cat tss_srt_m_gencode19.bed | awk '{print $3-$2}' | paste -sd+ - | bc >> sumbp.tsv
cat sumbp.tsv | paste - - > total_length.tsv
rm sumbp.tsv


# ----- ENRICHMENT ANALYSIS --------

# Create dataset of the +/-3kb region from the genebody for enrichment analysis
# remove the genes that are closer than 3kb to the chromosome end, because it would mess up the calculations downstream
cat gencode.v19.annotation.gtf | awk '{if($20=="\"protein_coding\";" && $3=="gene"){print $1 "\t" $4-3000 "\t" $4 "\t" $18 "\t0\t" $7}}' |\
awk '{if($2>=0){print $0}}' | sort -k1,1 -k2,2n > upstr_gene_srt_gencode19.bed
cat gencode.v19.annotation.gtf | awk '{if($20=="\"protein_coding\";" && $3=="gene"){print $1 "\t" $5 "\t" $5+3000 "\t" $18 "\t0\t" $7}}' |\
sort -k1,1 -k2,2n > dwnstr_gene_srt_gencode19.bed

# Make different files for - or + strand
# plus upstream = upstream
# plus downstream = downstream
# minus upstream = downstream & reversed
# minus downstream = upstream & reversed
cat upstr_gene_srt_gencode19.bed | awk '{if($6=="-"){print $0}}' > minus_upstr_gencode19.bed
cat upstr_gene_srt_gencode19.bed | awk '{if($6=="+"){print $0}}' > plus_upstr_gencode19.bed
cat dwnstr_gene_srt_gencode19.bed | awk '{if($6=="-"){print $0}}' > minus_dwnstr_gencode19.bed
cat dwnstr_gene_srt_gencode19.bed | awk '{if($6=="+"){print $0}}' > plus_dwnstr_gencode19.bed
cat gene_srt_gencode19.bed | awk '{if($6=="-"){print $0}}' > minus_gene_gencode19.bed
cat gene_srt_gencode19.bed | awk '{if($6=="+"){print $0}}' > plus_gene_gencode19.bed

# Create fixed size (100nt) sliding windows in these up- and downstream regions respecting the strandedness
# bedpath= ... in case you have no sudo right and bedtools is in a user folder
$bedpath"bedtools" makewindows -b plus_upstr_gencode19.bed -w 100 > plus_upstr_w100_g19.bed
$bedpath"bedtools" makewindows -b minus_upstr_gencode19.bed -w 100 > minus_upstr_w100_g19.bed
$bedpath"bedtools" makewindows -b plus_dwnstr_gencode19.bed -w 100 > plus_dwnstr_w100_g19.bed
$bedpath"bedtools" makewindows -b minus_dwnstr_gencode19.bed -w 100 > minus_dwnstr_w100_g19.bed

# Create fixed number (100) sliding window in gene bodies
$bedpath"bedtools" makewindows -b plus_gene_gencode19.bed -n 100 > plus_gene_n100_g19.bed
$bedpath"bedtools" makewindows -b minus_gene_gencode19.bed -n 100 > minus_gene_n100_g19.bed

# WARNING: Interval chr1:147706573-147706607 is smaller than the number of windows requested. Skipping.
# These are very short genes, so I just ignored them

# ------------------- second part: region information intersection with DSB data ----------------------


# step in the folder where you have your UMI-filtered bed files
cd $myfolder

# Bedtools intersect could be used for the counts of topX% of genes and TSS
# with -wa and -wb options to have both gene name and DSB count reported
# But bedtools coverage is better for further analysis
# coverage cannot take the DSB counts into account, so the files should be "expanded"
# meaning that the number of lines should correspond to the number of DSBs in the region
for file in *UMI.bed;\
do name=$(echo $file | cut -d"_" -f1);\
echo $name;\
cat $file | awk '{for(i=1;i<=$4;i++) print "chr" $0}' >$name"_exp.bed";\
done

# NOTE: this script above add a "chr" prefix to the chromsome name, if it already has one, the 4th line should look like this instead:
# cat $file | awk '{for(i=1;i<=$4;i++) print $0}' >$name"_exp.bed";\

# now coverage can be used
# gencodepath= ... the path to the previously prepared genecode stuff
# bedpath= ... in case you have no sudo right and bedtools is in a user folder

mkdir -p ../exon_gene_counts

# ----- REGIONAL COUNTS --------

for file in *exp.bed;\
do name=$(echo ../exon_gene_counts/$file | cut -d"_" -f1-4);\
echo $name;\
echo "gene";\
$bedpath"bedtools" coverage -counts -a $gencodepath"gene_srt_m_gencode19.bed" -b $file >$name"_gene.bed";\
echo "exon";\
$bedpath"bedtools" coverage -counts -a $gencodepath"exon_srt_m_gencode19.bed" -b $file >$name"_exon.bed";\
echo "tss";\
$bedpath"bedtools" coverage -counts -a $gencodepath"tss_srt_m_gencode19.bed" -b $file >$name"_tss.bed";\
done

# If you got this error message:
# Error: Unable to open file gene_srt_m_gencode19.bed. Exiting.
# That probably means that you have a spelling mistake in $gencodepath. You can test is with echo $gencodepath.

# The following piece cuts the filename in the second "_"
# name=$(echo ../exon_gene_counts/$file | cut -d"_" -f1-4)
# if you have only one ID in the file name, you should use this:
# name=$(echo ../exon_gene_counts/$file | cut -d"_" -f1-3)
# same goes for the piece below

# Do the same, but on non-merged gene and tss regions for later analysis on the top1% fragile genes/TSS

for file in *exp.bed;\
do name=$(echo ../exon_gene_counts/$file | cut -d"_" -f1-4);\
echo $name;\
echo "gene";\
$bedpath"bedtools" coverage -counts -a $gencodepath"gene_srt_gencode19.bed" -b $file >$name"_c_gene.bed";\
echo "tss";\
$bedpath"bedtools" coverage -counts -a $gencodepath"tss_srt_gencode19.bed" -b $file >$name"_c_tss.bed";\
done


# ----- ENRICHMENT ANALYSIS --------

for file in *exp.bed;\
do name=$(echo ../exon_gene_counts/$file | cut -d"_" -f1-4);\
echo $name;\
echo "upstream";\
$bedpath"bedtools" coverage -counts -a $gencodepath"plus_upstr_w100_g19.bed" -b $file >$name"_p_upstr_cov.bed";\
$bedpath"bedtools" coverage -counts -a $gencodepath"minus_upstr_w100_g19.bed" -b $file >$name"_m_upstr_cov.bed";\
echo "downstream";\
$bedpath"bedtools" coverage -counts -a $gencodepath"plus_dwnstr_w100_g19.bed" -b $file >$name"_p_dwnstr_cov.bed";\
$bedpath"bedtools" coverage -counts -a $gencodepath"minus_dwnstr_w100_g19.bed" -b $file >$name"_m_dwnstr_cov.bed";\
echo "genebody";\
$bedpath"bedtools" coverage -counts -a $gencodepath"plus_gene_n100_g19.bed" -b $file >$name"_p_gene_cov.bed";\
$bedpath"bedtools" coverage -counts -a $gencodepath"minus_gene_n100_g19.bed" -b $file >$name"_m_gene_cov.bed";\
done

# ------------------- third part: plotting ----------------------

# From this point R is used, because it’s good for managing data frames and plotting.
# DSBs are counted in these exonic regions. DSB number in intron is calculated as [(DSB# in gene) - (DSB# in exon)].
# DSB number in intergenic region is calculated as [(total # DSB) - (DSB# in gene)].
# the corresponding codes are described in the README

