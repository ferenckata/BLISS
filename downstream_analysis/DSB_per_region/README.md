## Analysis of DSBs in different genomic regions

1. Gencode database download and file processing to gain three type of data
- gene regions
- exon regions
- TSS regions

After retrieving the start and end coordinates of different regions, bedtools can be used to calculate the total number of DSBs in each of these regions. Importing the number of breaks per gene/exon/TSS region to R let us perform various analysis and plotting of the data.

The `DSBcountR.r` script calculates the total number of DSBs in the whole genome, in genes, in intergenic regions, in exons, in introns and in the pre-defined TSS region. It also finds the top10% of the genes and TSS regions with the most number of breaks.
