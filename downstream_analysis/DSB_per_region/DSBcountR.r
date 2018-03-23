# count and plot DSBs in exons, introns and intergenic regions

# ----------- DEPENDENCIES -------------
# install.packages('VennDiagram')
library(ggplot2)
library(reshape2)
library(wesanderson)
library(plyr)
library(RColorBrewer)
library(VennDiagram)

# ---------- INPUTS -------------
userpath = ""
gepath = ""
normfile = ""

files = list.files(path=userpath,pattern = "exp.bed")

# ----------- FUNCTIONS --------------

# There are four functions:
# 1) for counting the DSBs per region 
# 2) for finding the genes with the most DSBs and plotting the DSB count vs gene length
# 3) for plotting Venn diagram for 4 datasets
# 4) for plotting fractions

# 1) 
dsbcount = function(infile,i,outfile,dsbcounts){
  
  # ------------ overview on counts ----------------  
  
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
  
  # calculating the fraction of the total number of DSBs in each region
  frtotal = 1
  frgene = genedsb/totaldsb
  frintergene = intergendsb/totaldsb
  frexon = exondsb/totaldsb
  frintron = introndsb/totaldsb
  frtss = tssdsb/totaldsb
  
  dsbcounts[i,"totalfr"] = frtotal
  dsbcounts[i,"genefr"] = as.numeric(frgene)
  dsbcounts[i,"exonfr"] = as.numeric(frexon)
  dsbcounts[i,"intronfr"] = as.numeric(frintron)
  dsbcounts[i,"tssfr"] = as.numeric(frtss)
  dsbcounts[i,"intergenicfr"] = as.numeric(frintergene)
  
  return(dsbcounts)
  
}

# 2)
fragility = function(infile,i,fragdb,endng){
  
  # -------------- topX% fragile gene/TSS --------------------------
  
  # read non merged files and store the counts per gene
  runid = unlist(strsplit(infile,"_",fixed = T))[1]
  nm = paste(gepath,runid,endng,sep="")
  
  nmc = read.delim(nm,header = F,sep="\t",as.is = T)
  
  # plotting the number of DSBs vs the gene length
  if(endng=="_c_gene.bed"){
    png(paste(userpath,runid,"_dsb_vs_genelength.png",sep=""))
    plot(nmc[which(nmc$V7>0),3]-nmc[which(nmc$V7>0),2],nmc[which(nmc$V7>0),7],
         ylab = "DBS counts per gene",xlab = "gene length",main = runid)
    dev.off()
  }
  
  # calculating the TSS regions with the highest number of DSBs
  tt = floor(dim(nmc)[1]/100)
  srt = nmc[order(-nmc$V7),]
  topten = srt[1:tt,]
  # removing duplicated genes
  fr = unique(topten$V4, incomparables = FALSE, MARGIN = 1,fromLast = FALSE)
  fr = c(runid,fr)
  
  # concat to the database
  fragdb=rbind.fill(as.data.frame(t(fr)),fragdb)
  
  # write to file the top1% of TSS with the counts
  write.table(topten,file=paste(userpath,runid,endng,sep=""),quote = F,
              row.names = F,col.names = F,sep = '\t')
  
  return(fragdb)
  
}

# 3)
venner = function(fragdb,outname){
  # intersect the data
  fragm = as.matrix(fragdb)
  
  # only works with four datasets!!
  areas = matrix(data = NA,nrow = 15,ncol=2)
  
  for(k in 1:length(files)){
    areas[k,1] = fragm[k,1]
    areas[k,2] = sum(!is.na(fragm[k,]))-1
    if(k ==1){
      areas[k+4,1] = paste(fragm[k,1],fragm[k+1,1],sep=":")
      areas[k+4,2] = length(intersect(fragm[k,],fragm[k+1,]))-1
      areas[k+5,1] = paste(fragm[k,1],fragm[k+2,1],sep=":")
      areas[k+5,2] = length(intersect(fragm[k,],fragm[k+2,]))-1
      areas[k+6,1] = paste(fragm[k,1],fragm[k+3,1],sep=":")
      areas[k+6,2] = length(intersect(fragm[k,],fragm[k+3,]))-1
      areas[k+10,1] = paste(fragm[k,1],fragm[k+1,1],fragm[k+2,1],sep=":")
      areas[k+10,2] = length(intersect(intersect(fragm[k,],fragm[k+1,]),fragm[k+2,]))-1
      areas[k+11,1] = paste(fragm[k,1],fragm[k+1,1],fragm[k+3,1],sep=":")
      areas[k+11,2] = length(intersect(intersect(fragm[k,],fragm[k+1,]),fragm[k+3,]))-1
      areas[k+12,1] = paste(fragm[k,1],fragm[k+2,1],fragm[k+3,1],sep=":")
      areas[k+12,2] = length(intersect(intersect(fragm[k,],fragm[k+2,]),fragm[k+3,]))-1
      areas[k+14,1] = paste(fragm[k,1],fragm[k+1,1],fragm[k+2,1],fragm[k+3,1],sep=":")
      areas[k+14,2] = length(intersect(intersect(intersect(fragm[k,],fragm[k+1,]),fragm[k+2,]),fragm[k+3,]))-1
    }else if(k == 2){
      areas[k+6,1] = paste(fragm[k,1],fragm[k+1,1],sep=":")
      areas[k+6,2] = length(intersect(fragm[k,],fragm[k+1,]))-1
      areas[k+7,1] = paste(fragm[k,1],fragm[k+2,1],sep=":")
      areas[k+7,2] = length(intersect(fragm[k,],fragm[k+2,]))-1
      areas[k+12,1] = paste(fragm[k,1],fragm[k+1,1],fragm[k+2,1],sep=":")
      areas[k+12,2] = length(intersect(intersect(fragm[k,],fragm[k+1,]),fragm[k+2,]))-1
    }else if(k ==3){
      areas[k+7,1] = paste(fragm[k,1],fragm[k+1,1],sep=":")
      areas[k+7,2] = length(intersect(fragm[k,],fragm[k+1,]))-1
    }
  }
  
  # write out the list of genes broken in all datasets
  write.table(intersect(intersect(intersect(fragm[1,],fragm[2,]),fragm[3,]),fragm[4,]),
              file=paste(outname,'.txt',sep=""),quote = F,col.names = F)
  
  # save Venn diagram as jpeg file
  jpeg(paste(outname,".jpeg",sep = ""))
  draw.quad.venn(as.integer(areas[1,2]),as.integer(areas[2,2]),as.integer(areas[3,2]),
                 as.integer(areas[4,2]),as.integer(areas[5,2]),as.integer(areas[6,2]),
                 as.integer(areas[7,2]),as.integer(areas[8,2]),as.integer(areas[9,2]),
                 as.integer(areas[10,2]),as.integer(areas[11,2]),as.integer(areas[12,2]),
                 as.integer(areas[13,2]),as.integer(areas[14,2]),as.integer(areas[15,2]),
                 category = c(areas[1,1],areas[2,1],areas[3,1],areas[4,1]),
                 col = brewer.pal(4,"Paired"),fill = brewer.pal(4,"Paired"),
                 alpha = 0.75, cex=1.25,cat.cex = 1.25)
  dev.off()
}

# 4)
fracplot = function(pdfname,orig){
  origdf = data.frame(orig)
  meltedsb = melt(origdf,id.vars = 'ID')
  meltedsb$value <-as.numeric(meltedsb$value)
  pdf(paste(userpath,pdfname,sep=""))
  ggplot(meltedsb,aes(variable,value)) + geom_bar(aes(fill=ID), stat = "identity",position = "dodge") + scale_fill_brewer(palette = "Paired")
  dev.off()
}

# ----------------------------------- RUN ----------------------------------- 


# ----- OVERALL COUNTS -------

normcounts = read.delim(normfile,header = F,sep = '\t',as.is = T)
outfile = paste(userpath,"DSB_distribution.tsv",sep="")

dsbcounts = matrix(data=NA,nrow=length(files),ncol=19)
colnames(dsbcounts) = c("ID","total","totalnorm","totalfr",
                        "intergenic","intergenicnorm","intergenicfr",
                        "gene","genenorm","genefr","exon","exonnorm","exonfr",
                        "intron","intronnorm","intronfr",
                        "tss","tssnorm","tssfr")

i = 0
for(f in files){
  i = i+1
  dsbcounts = dsbcount(f,i,outfile,dsbcounts)
}

#-------------- save data ---------------
  
dsbframe = data.frame(dsbcounts)
write.csv(dsbframe,file=paste(userpath,'DSB_counts.csv',sep=""))

# -------------- plot ---------------------

# original counts 
orig = cbind(as.character(dsbframe$ID),as.character(dsbframe$total),
             as.character(dsbframe$intergenic),as.character(dsbframe$gene),
             as.character(dsbframe$exon),as.character(dsbframe$intron),
             as.character(dsbframe$tss))
colnames(orig) = c("ID","total","intergenic","gene","exon","intron","tss")
  
pdfname = 'DSB_count_distribution.pdf'
fracplot(pdfname,orig)

# normalized counts 
normed = cbind(as.character(dsbframe$ID),as.character(dsbframe$totalnorm),
             as.character(dsbframe$intergenicnorm),as.character(dsbframe$genenorm),
             as.character(dsbframe$exonnorm),as.character(dsbframe$intronnorm),
             as.character(dsbframe$tssnorm))
colnames(normed) = c("ID","total/bp","intergenic/bp","gene/bp","exon/bp",
                   "intron/bp","tss/bp")
  
pdfname = 'DSB_norm_distribution.pdf'
fracplot(pdfname,normed)

# fractions
frac = cbind(as.character(dsbframe$ID),as.character(dsbframe$totalfr),
               as.character(dsbframe$intergenicfr),as.character(dsbframe$genefr),
               as.character(dsbframe$exonfr),as.character(dsbframe$intronfr),
               as.character(dsbframe$tssfr))
colnames(frac) = c("ID","total/fraction","intergenic/fraction","gene/fraction",
                     "exon/fraction","intron/fraction","tss/fraction")
  
pdfname = 'DSB_fraction_distribution.pdf'
fracplot(pdfname,frac)


# ------------ FRAGILITY  ------------------

fragm = matrix(data=NA,nrow=1,ncol=1)
tfragdb = data.frame(fragm)
gfragdb = data.frame(fragm)

i = 0
for(f in files){
  i = i+1
  gfragdb = fragility(f,i,gfragdb,endng = "_c_gene.bed")
  tfragdb = fragility(f,i,tfragdb,endng = "_c_tss.bed")
}

# create Venn diagram

venner(gfragdb,paste(userpath,"gene_overlap",sep=""))
venner(tfragdb,paste(userpath,"tss_overlap",sep=""))




