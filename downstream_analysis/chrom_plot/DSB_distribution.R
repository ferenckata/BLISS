library("SNPchip")
library(ggplot2)
library(gtools)

inpath = ''
outpath = ''
chromsize = '~/Documents/bliss/data/chr_lengths.tsv'

# store chromosome sizes
print("loading chromosome file...")
chrom_size <- read.delim(chromsize,header = F)
files=list.files(inpath,pattern = ".bed")

# save all chromosome plots in the same file
outfile = paste(outpath,"all_DSB_distr.png",sep="")
png(outfile,width=4000,height = 5000)
par(mfrow = c(8,3),mar=c(6,4,4,2),cex.main=5)

for(infile in mixedsort(files)){
  chrname = unlist(strsplit(infile,'_',fixed=T))[1]
  print(chrname)
  if(chrname=="chrY" || chrname=="chrX"){
    chrnum = unlist(strsplit(chrname,"r",fixed=T))[2]
  }else{
    chrnum = as.integer(unlist(strsplit(chrname,"r",fixed=T))[2])
  }
  # store coverage file (per chromosome)
  print("loading coverage file...")
  dataByChr <- read.delim(paste(inpath,infile,sep=""),header = F,sep = "\t",as.is = T)
  chrlen = chrom_size[which(as.character(chrom_size$V1)==chrname),2]
  # plot and save to png
  print("Plotting...")
  # plotting DSBs
  plot(dataByChr$V2,dataByChr$V4,pch=16,xaxt = "n",
       las = 1,ylab = "DSB counts",xlim = c(1,chrlen),
       cex = 1,ylim=c(-40,500),
       xlab="",main=chrname)
  # plotting karyogram
  pI <- plotIdiogram(chromosome = chrnum,build = "hg19",unit = "bp",label.y = -0.35,new=F,ylim=c(-15,-30),label.cytoband = F)
  sprintf("Plotting of %s is done", chrname)
}

dev.off()
