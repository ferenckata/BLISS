# Processing of the BLISS sequencing results
# this script can be called as
# bliss_main_kf.sh infastq inputdir whitelist refgenome pydir quality

# INPUTS
inputdir = $1		# full path to the fastq file
infastq = $2		# the name of the gz compressed fastq file without the ".gz" ending (e.i. rm31.fastq)
WL = $3		    # the name of the file containing a list of the experiment barcodes (one barcode per line, no header)
refgenome = $4		# full path to reference genome fastq file or fasta file if it has been indexed already
pydir = $5		# full path to the mismatchmaker.py file
q = $6			# mapping quality threshold for filtering

# DEPENDENCIES
# umi_tools
# samtools
# bwa mem
# python3
# matchmaker.py

# CODE
# unzipping the fastq file of interest
gunzip $inputdir"/"$infastq".gz"
# retrieving the run ID
name = $(echo $infastq | cut -d'.' -f1)

# creating all possible mismatches of the experiment barcode
cat $WL | while read LINE; do python3 $pydir"/"mismatchmaker.py $LINE >> mmWL ;done

# experiment barcode filtering with umi_tools
umi_tools extract -p NNNNNNNNCCCCCCCC -I $inputdir"/"$infastq --filter-cell-barcode --whitelist=mmWL > $name"_ft.fastq"

# mapping the barcode filtered fastq files
# bwa index $refseqfq      # only for the first time
# bwa mem $refseqfa $outfq > $outsam
bwa mem $refgenome $name"_ft.fastq" > $name"_ft_mp.sam"

# quality filtering, converting, sorting and indexing sam to bam for later use
samtools view -Sb -q $q $name"_ft_mp.sam" > $name"_ft_mp.bam"
samtools sort $name"_ft_mp.bam" -o $name"_ft_mp_srt.bam"
samtools index $name"_ft_mp_srt.bam"

# grouping the UMIs with umi_tools
umi_tools group -I $name"_ft_mp_srt.bam" --group-out=groups.tsv --per-cell --output-bam -S $name"_ft_mp_srt_gp.bam"

# converting final bam to sam
samtools view -h -o $name"_ft_mp_srt_gp.sam" $name"_ft_mp_srt_gp.bam"
