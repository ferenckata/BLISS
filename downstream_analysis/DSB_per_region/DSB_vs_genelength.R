library(biomaRt)
library(plyr)

# plot dsb number vs gene length in four color for four different files

userpath = ""
gepath = ""
normfile = ""
outpath = ""

files = list.files(path=userpath,pattern = "exp.bed")
idList = c("ID1","ID2","ID3","ID4","ID5")

i=0
par(mfrow=c(1,1))
for(infile in files){
  i = i+1
  # read non merged files and store the counts per gene
  runid = unlist(strsplit(infile,"_",fixed = T))[1]
  nm = paste(gepath,runid,"_c_gene.bed",sep="")
  
  nmc = read.delim(nm,header = F,sep="\t",as.is = T)
  
  print(summary(nmc$V7))
  thrh = quantile(nmc$V7,probs = 0.75)
  
  pdf(paste(outpath,runid,"_boxplot.pdf",sep=""))
  # plotting the number of DSBs vs the gene length
  boxplot(nmc$V7/abs((nmc$V3-nmc$V2)),main=runid,
          data=nmc,names=runid)
  dev.off()
  
  # find outliers, plot them and calculate linear regression for them
  outlier = nmc[which(nmc$V7>thrh),]
  outlier$V4 = gsub(";","",outlier$V4)
  sprintf("There are %s outliers in %s dataset.",as.character(dim(outlier)[1]),runid)
  
  print(summary(outlier$V7))
  minth = quantile(outlier$V7,probs = 0)
  
  pdf(paste(outpath,runid,"_scatterplot.pdf",sep = ""))
  
  scatter.smooth(nmc$V3-nmc$V2,nmc$V7,ylab = "DBS counts per gene",
                 xlab = "gene length",main = runid,pch=20,cex=0.5,
                 lpars=list(col="red",lwd=2))
  
  # plot(nmc$V3-nmc$V2,nmc$V7,ylab = "DBS counts per gene",xlab = "gene length",main = runid,pch=20,cex=0.5)
  # abline(lm(nmc[which(nmc$V7<minth),3]-nmc[which(nmc$V7<minth),2]~nmc[which(nmc$V7<minth),7]),col="red")
  # NOTE : the abline seems a good idea, but the values around 0 are usually so many that it looks weird
  # looks better to use scatter.smooth that also plots the dots
  # BUT it may not be the statistically sound approach
  
  dev.off()
  
  # calculate correlation
  sprintf("The correlation between DSB number and gene length in the genes is %f",
          cor(nmc$V3-nmc$V2,nmc$V7))
  # report linear regression model
  print(summary(lm(nmc$V3-nmc$V2 ~ nmc$V7)))
  
  # save outlier genes
  write.table(outlier,file = paste(outpath,runid,"_fragile_gene.tsv",sep=""),
              quote = F,sep = "\t",
              row.names = F,col.names = c("chromosome","start","end","gene symbols","NA","strand","DSB count"))
  
  # get GO terms for the outlier genes
  mymart = useMart("ensembl","hsapiens_gene_ensembl")
  listAttributes(mymart)
  goforols = getBM(attributes = c('hgnc_symbol','name_1006'),
                   filters = 'hgnc_symbol',
                   values = outlier$V4,
                   mart = mymart)
  write.table(goforols,file = paste(outpath,runid,"_fragile_gene_GO.tsv",sep=""),quote = F,sep = "\t",
              row.names = F,col.names = c("gene symbols","GO term"))
  # sort by GO frequency
  gofreq = count(goforols,vars = 'name_1006')
  gofreq = gofreq[order(-gofreq$freq),]
  # save GO counts
  write.table(gofreq,file=paste(outpath,runid,"_fragile_gene_GOfreq.tsv",sep=""),quote = F,sep='\t',row.names = F,col.names = T)
}

# # # If you prefer colored dots on top of each other for comparison  # # #

# one can include more colors here and more IDs in the idList if there are more than 5 datasets
dsbCol = c("darkgreen","palegreen","gold","steelblue3","yellow")
names(dsbCol) = idList

for(infile in files){
  i = i+1
  
  # read non merged files and store the counts per gene
  runid = unlist(strsplit(infile,"_",fixed = T))[1]
  nm = paste(gepath,runid,"_c_gene.bed",sep="")
  
  nmc = read.delim(nm,header = F,sep="\t",as.is = T)
  colnm = dsbCol[runid]
  
  # plotting the number of DSBs vs the gene length
  
  if(i==1){
    pdf(paste(outpath,"all_dsb_vs_genelength.pdf",sep=""))
    plot(nmc[which(nmc$V7>0),3]-nmc[which(nmc$V7>0),2],nmc[which(nmc$V7>0),7],
         ylab = "DBS counts per gene",xlab = "gene length",main = "DSB vs gene length",
         ylim=c(0,24000),col= colnm,pch=20,cex=0.5)
  }else{
    points(nmc[which(nmc$V7>0),3]-nmc[which(nmc$V7>0),2],nmc[which(nmc$V7>0),7],
           col= colnm,pch=20,cex=0.5)
  }
}  
dev.off()

  
  
