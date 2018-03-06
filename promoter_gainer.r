library(biomaRt)
grch37 = useEnsembl(biomart="regulation",GRCh=37)
listDatasets(grch37)
regulatory=useEnsembl(biomart = "regulation",
                     dataset="hsapiens_regulatory_feature",
                     GRCh=37)
listFilters(regulatory)
listAttributes(regulatory)
promoters = getBM(attributes = c('chromosome_name','chromosome_start','chromosome_end','feature_type_name'),
                  filters='regulatory_feature_type_name',values='promoter',mart=regulatory)
write.table(promoters,file='promoter_regions.tsv',col.names = F,row.names = F,quote = F,sep="\t")
