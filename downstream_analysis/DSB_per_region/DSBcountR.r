# count and plot DSBs in exons, introns and intergenic regions

library(ggplot2)
library(reshape2)
library(wesanderson)

userpath = "~/Documents/bliss/MCF7/MCF7_BICRO63/rtest/originals"
gepath = "~/Documents/bliss/MCF7/MCF7_BICRO63/rtest/exon_gene_counts"

dsbcount = function(infile,i,outfile,dsbcounts){
  
  # read files and store the counts
  runid = unlist(strsplit(infile,"_",fixed = T))[1]
  genefile = paste(gepath,'/',runid,"_gene.bed",sep="")
  exonfile = paste(gepath,'/',runid,"_exon.bed",sep="")
  # for total number of DSB per TSS this file
  tssfile = paste(gepath,'/',runid,"_tss.bed",sep="")
  # for topten TSS of genes with DSB
  gtssfile = paste(gepath,'/',runid,"_g_tss.bed",sep="")
  
  totalcounts = read.delim(paste(userpath,infile,sep="/"),header = F,sep="\t",as.is = T)
  genecounts = read.delim(genefile,header = F,sep="\t",as.is = T)
  exoncounts = read.delim(exonfile,header = F,sep="\t",as.is = T)
  tsscounts = read.delim(tssfile,header = F,sep="\t",as.is = T)
  gtsscounts = read.delim(gtssfile,header = F,sep="\t",as.is = T)
  
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
  
  return(dsbcounts)
  
  wline = sprintf("The total number of DSB in %s:\t%d\nin intergenic regions:\t%d\nin genes:\t%d\nin exons:\t%d\nin introns:\t%d\nin +/-2.5kb of coding TSS:\t%d\n",runid,totaldsb,intergendsb,genedsb,exondsb,introndsb,tssdsb)
  cat(text=wline,file=outfile,sep='\n',append=T)
  
  # calculating the TSS and GeneBody regions with the highest number of DSBs
  toptent = floor(dim(gtsscounts)[1]/10)
  srttssc = gtsscounts[order(-gtsscounts$V5),]
  toptentss = srttssc[1:toptent,]
  tttss[i,"ID"] = runid
  tttss[i,"max"] = toptentss$V5[1]
  tttss[i,"median"] = median(toptentss$V5)
  tttss[i,"min"] = min(toptentss$V5)
  
  fragtss = "?"
  
  topteng = floor(dim(genecounts)[1]/10)
  srtgenec = genecounts[order(-genecounts$V5),]
  toptengene = srtgenec[1:topteng,]
  ttgene[i,"ID"] = runid
  ttgene[i,"max"] = toptengene$V5[1]
  ttgene[i,"median"] = median(toptengene$V5)
  ttgene[i,"min"] = min(toptengene$V5)
  
}

outfile = paste(userpath,'/',"DSB_distribution.tsv",sep="")
files = list.files(path=userpath,pattern = "exp.bed")
dsbcounts = matrix(data=NA,nrow=length(files),ncol=7)
colnames(dsbcounts) = c("ID","total","intergenic","gene","exon","intron","tss")
tttss = matrix(data=NA,nrow=length(files),ncol=4)
colnames(tttss) = c("ID","max","median","min")
ttgene = matrix(data=NA,nrow=length(files),ncol=4)
colnames(ttgene) = c("ID","max","median","min")
i = 0
for(f in files){
  i = i+1
  dsbcounts = dsbcount(f,i,outfile,dsbcounts)
}

dsbframe = data.frame(dsbcounts)
meltedsb = melt(dsbframe,id.vars = 'ID')
meltedsb$value <-as.numeric(meltedsb$value)
ggplot(meltedsb,aes(variable,value)) + geom_bar(aes(fill=ID), stat = "identity",position = "dodge") + scale_fill_brewer(palette = "Paired")



