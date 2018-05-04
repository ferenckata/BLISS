library(biomaRt)
library(plyr)

# plot dsb number vs gene length in four color for four different files

userpath = "~/Documents/bliss/MCF7/MCF7_BICRO63/mapped_filtered/"
gepath = "~/Documents/bliss/MCF7/MCF7_BICRO63/rtest/exon_gene_counts/"
normfile = "~/Documents/gene_db/hg19/total_length.tsv"
outpath = "~/Documents/bliss/MCF7/MCF7_BICRO63/rtest/"

files = list.files(path=userpath,pattern = "exp.bed")
idList = c("RM118","RM119","RM120","RM121")

i = 0

par(mfrow=c(1,1))
i=0
for(infile in files){
  i = i+1
  # read non merged files and store the counts per gene
  runid = unlist(strsplit(infile,"_",fixed = T))[1]
  nm = paste(gepath,runid,"_c_gene.bed",sep="")
  
  nmc = read.delim(nm,header = F,sep="\t",as.is = T)
  colnm = dsbCol[runid]
  
  print(summary(nmc$V7))
  thrh = quantile(nmc$V7,probs = 0.75)
  
  pdf(paste(runid,"_boxplot.pdf",sep=""))
  par(mfrow=c(1,1))
  # plotting the number of DSBs vs the gene length
  boxplot(nmc$V7/abs((nmc$V3-nmc$V2)),main=runid,
          data=nmc,names=runid,col = colnm)
  dev.off()
  
  # find outliers, plot them and calculate linear regression for them
  outlier = nmc[which(nmc$V7>thrh),]
  sprintf("There are %s outliers in %s dataset.",as.character(dim(outlier)[1]),runid)
  
  print(summary(outlier$V7))
  minth = quantile(outlier$V7,probs = 0)
  
  pdf(paste(runid,"_scatterplot.pdf",sep = ""))
  par(mfrow=c(2,1))
  
  plot(nmc[which(nmc$V7<minth),3]-nmc[which(nmc$V7<minth),2],nmc[which(nmc$V7<minth),7],
       ylab = "DBS counts per gene",xlab = "gene length",main = runid,pch=20,cex=0.5)
  
  # abline(lm(nmc[which(nmc$V7<minth),3]-nmc[which(nmc$V7<minth),2]~nmc[which(nmc$V7<minth),7]),col="red")
  # NOTE : the abline seems a good idea, but the values around 0 are usually so many that it looks weird
  # better to use scatter.smooth that also plots the dots
  # 
  
  plot(outlier$V3-outlier$V2,outlier$V7,
       ylab = "DBS counts per gene",xlab = "gene length",main = runid,pch=20,cex=0.5)
  abline(lm(outlier$V3-outlier$V2 ~ outlier$V7),col="red")
  
  dev.off()
  
  # calculate correlation
  sprintf("The correlation between DSB number and gene length in the outlier genes is %f",
         cor(outlier$V3-outlier$V2,outlier$V7))
  # report linear regression model
  print(summarylm(outlier$V3-outlier$V2 ~ outlier$V7))
  
  # get GO terms for the outlier genes
  mymart = useMart("ensembl","hsapiens_gene_ensembl")
  listAttributes(mymart)
  goforols = getBM(attributes = c('hgnc_symbol','name_1006'),
                  filters = 'hgnc_symbol',
                   values = outlier$V4,
                   mart = mymart)
  write.table(goforols,file = "fragile_gene_GO.tsv",quote = F,sep = "\t",row.names = F,col.names = c("gene symbols","GO term"))
  # sort by GO frequency
  gofreq = count(goforols,vars = 'name_1006')
  gofreq = gofreq[order(-gofrec$freq),]
  # save GO counts
  write.table(gofreq,file="fragile_gene_GOfreq.tsv",quote = F,sep='\t',row.names = F,col.names = T)
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

  
  
