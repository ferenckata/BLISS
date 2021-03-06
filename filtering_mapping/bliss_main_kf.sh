# Processing of the BLISS sequencing results
# this script can be called as
# bliss_main_kf.sh infastq inputdir whitelist refgenome pydir quality

# INPUTS (it will be in a config file soon...)
inputdir = $1		# full path to the fastq file
infastq = $2		# the name of the gz compressed fastq file without the ".gz" ending (e.i. rm31.fastq)
WL = $3		        # the name of the file containing a list of the experiment barcodes (one barcode per line, no header)
refgenome = $4		# full path to reference genome fastq file or fasta file if it has been indexed already
pydir = $5		# full path to the mismatchmaker.py file
q = $6			# mapping quality threshold for filtering

# DEPENDENCIES
# umi_tools
# samtools
# bwa mem
# python3
# mismatchmaker.py

# CODE
# unzipping the fastq file of interest
gunzip $inputdir"/"$infastq".gz"
# retrieving the run ID
name = $(echo $infastq | cut -d'.' -f1)

# creating all possible mismatches of the experiment barcode
# each mismatched barcode can be saved to different files, and the barcode filtering could be run on each file
cat $WL | while read LINE; do python3 $pydir"/"mismatchmaker.py $LINE >> mmWL ;done

# experiment barcode filtering with umi_tools
umi_tools extract -p NNNNNNNNCCCCCCCC -I $inputdir"/"$infastq --filter-cell-barcode --whitelist=mmWL > $name"_ft.fastq"
# remove the comment lines from the beginning of the file to map (otherwise samtools crashes)
tail -n+35 $name"_ft.fastq" > $name"_nhft.fastq"

# mapping the barcode filtered fastq files
# bwa index $refseqfq      # only for the first time
# bwa mem $refseqfa $outfq > $outsam
bwa mem $refgenome $name"_nhft.fastq" > $name".sam"

# converting, quality filtering, sorting and indexing sam to bam for later use
samtools view -Sb -q $q $name".sam" > $name".bam"
samtools sort $name".bam" -o $name"_srt.bam"
samtools index $name"_srt.bam"

# deduplicating the UMIs by cells with umi_tools
umi_tools dedup -I $name"_srt.bam" --edit-distance-threshold 2 --per-cell -S $name"_dedup.bam" -L $name"_dedup.log"

# retrieving the strand information from the groupped BAM file
samtools view -Xf 0x10 $name"_srt_gp.bam" | awk '{print $1}' > $name"_rev.tsv"
samtools view -XF 0x10 $name"_srt_gp.bam" | awk '{print $1}' > $name"_fwd.tsv"

# converting final bam to sam
# samtools view -h -o $name"_ft_mp_srt_gp.sam" $name"_ft_mp_srt_gp.bam"
