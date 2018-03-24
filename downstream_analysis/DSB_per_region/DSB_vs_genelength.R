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

# # # # # # # # INDIVIDUAL DATASETS # # # # # # #

# one can include more colors here and more IDs in the idList if there are more than 4 datasets
dsbCol = c("darkgreen","palegreen","gold","steelblue3")
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
    png(paste(outpath,"all_dsb_vs_genelength.png",sep=""))
    plot(nmc[which(nmc$V7>0),3]-nmc[which(nmc$V7>0),2],nmc[which(nmc$V7>0),7],
         ylab = "DBS counts per gene",xlab = "gene length",main = "DSB vs gene length",
         ylim=c(0,24000),col= colnm,pch=16)
  }else{
    points(nmc[which(nmc$V7>0),3]-nmc[which(nmc$V7>0),2],nmc[which(nmc$V7>0),7],
           col= colnm,pch=16)
  }
}  
dev.off()


# ---------- boxplot ---------------
par(mfrow=c(2,2))
i=0
for(infile in files){
  i = i+1
  # read non merged files and store the counts per gene
  runid = unlist(strsplit(infile,"_",fixed = T))[1]
  nm = paste(gepath,runid,"_c_gene.bed",sep="")
  
  nmc = read.delim(nm,header = F,sep="\t",as.is = T)
  colnm = dsbCol[runid]
  
  # plotting the number of DSBs vs the gene length
  # png(paste(userpath,"all_dsb_vs_genelength.png",sep=""))
  boxplot(nmc[which(nmc$V7>0),7]/abs((nmc[which(nmc$V7>0),3]-nmc[which(nmc$V7>0),2])),data=nmc,names=runid,col = colnm,ylim=c(0,15))
}

# # # # # # POOLED DATASETS # # # # # # #

userpath = "~/Documents/bliss/MCF7/MCF7_BICRO63/rtest/exon_gene_counts/"
i=0
# read in all files
for(infile in list.files(userpath,pattern = "_c_gene.bed")){
  i = i+1
  if(i==1){
    nm1 = read.delim(paste(userpath,infile,sep=""),header = F,sep="\t",as.is = T)
    print("hey")
  }else if(i==2){
    nm2 = read.delim(paste(userpath,infile,sep=""),header = F,sep="\t",as.is = T)
  }else if(i==3){
    nm3 = read.delim(paste(userpath,infile,sep=""),header = F,sep="\t",as.is = T)
  }else if(i==4){
    nm4 = read.delim(paste(userpath,infile,sep=""),header = F,sep="\t",as.is = T)
  }
}
# create an empty data frame to store the location, gene symbol and DSB counts
nmall = cbind(nm1$V1,nm1$V2,nm1$V3,nm1$V4,0)

# merge the counts from each
for(n in 1:dim(nm1)[1]){
  nmall[n,5] = sum(nm1$V7[n],nm2$V7[n],nm3$V7[n],nm4$V7[n])
}

nmalldf = data.frame(nmall)

# make scatterplot
png("all_DSB_vs_genelength.png",width=600,height=700)
plot(as.integer(nmall[which(as.integer(nmall[,5])>0),3])-as.integer(nmall[which(as.integer(nmall[,5])>0),2]),
     as.integer(nmall[which(as.integer(nmall[,5])>0),5]),ylab = "DBS counts per gene",
     xlab = "gene length",main = "DSB vs gene length",pch=16)
dev.off()
# plot without outliers (genes with >20 000 DSB)
png("all_DSB_vs_genelength_no_outl.png",width=600,height=700)
plot(as.integer(nmall[which(as.integer(nmall[,5])>0 & as.integer(nmall[,5])<20000),3])-as.integer(nmall[which(as.integer(nmall[,5])>0 & as.integer(nmall[,5])<20000),2]),
     as.integer(nmall[which(as.integer(nmall[,5])>0 & as.integer(nmall[,5])<20000),5]),ylab = "DBS counts per gene",
     xlab = "gene length",main = "DSB vs gene length",pch=16)
dev.off()

# set a limit and pick all genes above that
outliers = nmall[which(as.integer(nmall[,5])>20000),4]

for(i in 1:length(outliers)){
  outliers[i] = unlist(strsplit(outliers[i],";",fixed=T)[1])[1]
}

# get GO terms for these
mymart = useMart("ensembl","hsapiens_gene_ensembl")
listAttributes(mymart)
goforols = getBM(attributes = c('hgnc_symbol','name_1006'),
                 filters = 'hgnc_symbol',
                 values = outliers,
                 mart = mymart)
write.table(goforols,file = "fragile_gene_GO.tsv",quote = F,sep = "\t",row.names = F,col.names = c("gene symbols","GO term"))
# sort by GO frequency
gofreq = count(goforols,vars = 'name_1006')
gofreq = gofreq[order(-gofrec$freq),]
# save GO counts
write.table(gofreq,file="fragile_gene_GOfreq.tsv",quote = F,sep='\t',row.names = F,col.names = T)

  
  
