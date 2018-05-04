# DSB vs chromosome length

library(ggplot2)
library(gtools)

inpath = ''
outpath = ''
chromsize = ''

# store chromosome sizes
print("loading chromosome file...")
chrom_size <- read.delim(chromsize,header = F)
files=mixedsort(list.files(inpath,pattern = "UMI.bed"))

outfile = paste(outpath,"DSB_vs_chrlength.pdf",sep="")
DSBnum = matrix(data = NA,nrow = 24,ncol = 4)
i = 0
for(infile in files){
  i = i+1
  chrname = unlist(strsplit(unlist(strsplit(infile,'_',fixed=T))[1],"r",fixed = T)[1])[2]
  print(chrname)
  # store coverage file (per chromosome)
  print("loading coverage file...")
  dataByChr <- read.delim(paste(inpath,infile,sep=""),header = F,sep = "\t",as.is = T)
  chrlen = chrom_size[which(as.character(chrom_size$V1)==paste("chr",chrname,sep = "")),2]
  DSBnum[i,1] = chrlen
  DSBnum[i,2] = sum(dataByChr$V4)
  DSBnum[i,3] = chrname
}
DSBdf = data.frame(DSBnum)
# for numerical axis
DSBdf$X1<-as.integer(DSBnum[,1])
DSBdf$X2<-as.integer(DSBnum[,2])

print("plotting...")
pdf(outfile,width=600,height = 600)
ggplot(data=DSBdf,aes(DSBdf$X1,DSBdf$X2)) + geom_point() + geom_text(aes(label=DSBdf$X3),hjust=0.7,vjust=-0.9) + labs(x="Chromosome length",y="DSB number") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
dev.off()
