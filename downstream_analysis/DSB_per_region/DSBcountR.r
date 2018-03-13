# count and plot DSBs in exons, introns and intergenic regions

library(ggplot2)
library(reshape2)
library(wesanderson)

userpath = "~/Documents/bliss/MCF7/MCF7_BICRO63/rtest/originals/"
gepath = "~/Documents/bliss/MCF7/MCF7_BICRO63/rtest/exon_gene_counts/"
normfile = "~/Documents/gene_db/hg19/total_length.tsv"

dsbcount = function(infile,i,outfile,dsbcounts){
  
  # read files and store the counts
  runid = unlist(strsplit(infile,"_",fixed = T))[1]
  genefile = paste(gepath,runid,"_gene.bed",sep="")
  exonfile = paste(gepath,runid,"_exon.bed",sep="")
  tssfile = paste(gepath,runid,"_tss.bed",sep="")
  
  totalcounts = read.delim(paste(userpath,infile,sep=""),header = F,sep="\t",as.is = T)
  genecounts = read.delim(genefile,header = F,sep="\t",as.is = T)
  exoncounts = read.delim(exonfile,header = F,sep="\t",as.is = T)
  tsscounts = read.delim(tssfile,header = F,sep="\t",as.is = T)
  
  # calculating the number of breaks in each region
  totaldsb = sum(totalcounts$V4)
  genedsb = sum(genecounts$V5)
  exondsb = sum(exoncounts$V5)
  tssdsb = sum(tsscounts$V5)
  introndsb = genedsb - exondsb
  intergendsb = totaldsb - genedsb
  
  dsbcounts[i,"ID"] = runid
  dsbcounts[i,"total"] = as.numeric(totaldsb)
  dsbcounts[i,"gene"] = as.numeric(genedsb)
  dsbcounts[i,"exon"] = as.numeric(exondsb)
  dsbcounts[i,"intron"] = as.numeric(introndsb)
  dsbcounts[i,"tss"] = as.numeric(tssdsb)
  dsbcounts[i,"intergenic"] = as.numeric(intergendsb)
  
  # calculating the number of breaks per bp in each region
  totalnorm = totaldsb/normcounts[which(normcounts$V1=="total"),2]
  genenorm = genedsb/normcounts[which(normcounts$V1=="gene"),2]
  exonnorm = exondsb/normcounts[which(normcounts$V1=="exon"),2]
  tssnorm = tssdsb/normcounts[which(normcounts$V1=="TSS"),2]
  intronnorm = introndsb/(normcounts[which(normcounts$V1=="gene"),2]-normcounts[which(normcounts$V1=="exon"),2])
  intergenenorm = intergendsb/(normcounts[which(normcounts$V1=="total"),2]-normcounts[which(normcounts$V1=="gene"),2])
  
  dsbcounts[i,"totalnorm"] = as.numeric(totalnorm)
  dsbcounts[i,"genenorm"] = as.numeric(genenorm)
  dsbcounts[i,"exonnorm"] = as.numeric(exonnorm)
  dsbcounts[i,"intronnorm"] = as.numeric(intronnorm)
  dsbcounts[i,"tssnorm"] = as.numeric(tssnorm)
  dsbcounts[i,"intergenicnorm"] = as.numeric(intergenenorm)
  
  # calculating the TSS and GeneBody regions with the highest number of DSBs
  toptent = floor(dim(tsscounts)[1]/100)
  srttssc = tsscounts[order(-tsscounts$V5),]
  toptentss = srttssc[1:toptent,]
  
  write.table(toptentss,file=paste(userpath,runid,"_top1_TSS.bed",sep=""),quote = F,
              row.names = F,col.names = F,sep = '\t')
  
  topteng = floor(dim(genecounts)[1]/100)
  srtgenec = genecounts[order(-genecounts$V5),]
  toptengene = srtgenec[1:topteng,]
  
  write.table(toptengene,file=paste(userpath,runid,"_top1_gene.bed",sep=""),quote = F,
              row.names = F,col.names = F,sep = '\t')
  
  return(dsbcounts)
  
}

#### RUN ####
normcounts = read.delim(normfile,header = F,sep = '\t')
outfile = paste(userpath,"DSB_distribution.tsv",sep="")
files = list.files(path=userpath,pattern = "exp.bed")
dsbcounts = matrix(data=NA,nrow=length(files),ncol=13)
colnames(dsbcounts) = c("ID","total","totalnorm","intergenic","intergenicnorm",
                        "gene","genenorm","exon","exonnorm","intron","intronnorm",
                        "tss","tssnorm")
i = 0
for(f in files){
  i = i+1
  dsbcounts = dsbcount(f,i,outfile,dsbcounts)
}

#----- save data ------
dsbframe = data.frame(dsbcounts)
write.csv(dsbframe,file=paste(userpath,'DSB_counts.csv',sep=""))

#### PLOT ####
# ---- original counts -------
orig = cbind(as.character(dsbframe$ID),as.character(dsbframe$total),
             as.character(dsbframe$intergenic),as.character(dsbframe$gene),
             as.character(dsbframe$exon),as.character(dsbframe$intron),
             as.character(dsbframe$tss))
colnames(orig) = c("ID","total","intergenic","gene","exon","intron","tss")
origdf = data.frame(orig)
meltedsb = melt(origdf,id.vars = 'ID')
meltedsb$value <-as.numeric(meltedsb$value)
pdf(paste(userpath,'DSB_count_distribution.pdf',sep=""))
ggplot(meltedsb,aes(variable,value)) + geom_bar(aes(fill=ID), stat = "identity",position = "dodge") + scale_fill_brewer(palette = "Paired")
dev.off()

# ---- normalized counts -----
normed = cbind(as.character(dsbframe$ID),as.character(dsbframe$totalnorm),
             as.character(dsbframe$intergenicnorm),as.character(dsbframe$genenorm),
             as.character(dsbframe$exonnorm),as.character(dsbframe$intronnorm),
             as.character(dsbframe$tssnorm))
colnames(normed) = c("ID","total/bp","intergenic/bp","gene/bp","exon/bp",
                   "intron/bp","tss/bp")
normeddf = data.frame(normed)
nmeltedsb = melt(normeddf,id.vars = 'ID')
nmeltedsb$value <-as.numeric(nmeltedsb$value)
pdf(paste(userpath,'DSB_norm_distribution.pdf',sep=""))
ggplot(nmeltedsb,aes(variable,value))+ geom_bar(aes(fill=ID),stat = "identity",position = "dodge") + scale_fill_brewer(palette = "Paired")
dev.off()


