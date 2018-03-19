# dsb number vs gene length

userpath = ""
gepath = ""

files = list.files(path=userpath,pattern = "exp.bed")
idList = c("RM118","RM119","RM120","RM121")

i = 0
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
    png(paste(userpath,"all_dsb_vs_genelength.png",sep=""))
    plot(nmc[which(nmc$V7>0),3]-nmc[which(nmc$V7>0),2],nmc[which(nmc$V7>0),7],
         ylab = "DBS counts per gene",xlab = "gene length",main = "DSB vs gene length",
         ylim=c(0,24000),col= colnm,pch=16)
  }else{
    points(nmc[which(nmc$V7>0),3]-nmc[which(nmc$V7>0),2],nmc[which(nmc$V7>0),7],
           col= colnm,pch=16)
  }
}  
dev.off()
