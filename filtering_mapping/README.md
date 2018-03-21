## BLISS mapping and UMI and barcode filtering

_Dependences_
- [UMI tools](https://github.com/CGATOxford/UMI-tools)
- [mismathcmaker.py](https://github.com/ferenckata/BLISS/blob/master/filtering_mapping/mismatchmaker.py)
- [BWA MEM](http://bio-bwa.sourceforge.net/bwa.shtml)
- [samtools](http://www.htslib.org/doc/samtools.html)

_How to run the scripts_
It is better to run `bliss_main_kf.sh` line by line, at least for the first few times.
The script calls several tools. Read the manual for each tool carefully!
The file structure is not fixed, you should optimize it for yourself.

_What does it do?_
`bliss_main_kf.sh` is the main script.
After unzipping the fastq file, it prepares the whitelist of barcodes for `UMItools extract`, which will do two things:
 - filter for the experiment barcode
 - relocate the experiment barcode and UMI from the sequence to the header
 
 Next step is the mapping that is done by the BWA MEM algorithm.
 The output is converted to the binary version with samtools for further analysis.
 `UMItools dedup` will group the UMIs and filter out pcr duplicates. The current version groups the UMIs with respet to the experiment barcode. If the whitelist has all the barcodes, it is needed. If you prefer to keep the experiments separated and run everything X-times, the `--per-cell` option can be removed.
