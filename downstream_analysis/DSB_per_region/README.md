## Analysis of DSBs in different genomic regions

Gencode database download and file processing to gain three type of data
- gene regions
- exon regions
- TSS regions

After retrieving the start and end coordinates of different regions, bedtools can be used to calculate the total number of DSBs in each of these regions. Importing the number of breaks per gene/exon/TSS region to R let us perform various analysis and plotting of the data.

The `DSBcountR.r` script
- calculates the total number of DSBs in the whole genome, in genes, in intergenic regions, in exons, in introns and in the pre-defined TSS region
- plots the number of DSBs in these regions across datasets
- finds the top1% of the genes and TSS regions with the most number of breaks
- plots the overlap between top1% fragile genes/TSSs across datasets
- plots DSB counts versus gene length
