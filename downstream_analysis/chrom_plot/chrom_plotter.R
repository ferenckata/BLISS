# plot DSB data in sliding windows
# http://bioconductor.org/packages/release/bioc/vignettes/TitanCNA/inst/doc/TitanCNA.pdf
# based on library(TitanCNA)
# library(argparser)
library("SNPchip")

inpath = '' # with all 
outpath = '' # should be different from inpath
chromsize = ''
ID = '' #

# store chromosome sizes
print("loading chromosome file...")
chrom_size <- read.delim(chromsize,header = F)
files=list.files(inpath,pattern = ID) 

for(infile in list.files(inpath)){
  chrname = unlist(strsplit(infile,'_',fixed=T))[2]
  print(chrname)
  # store chromosome number for plotIdiogram
  if(chrname=="chrY" || chrname=="chrX"){
    chrnum = unlist(strsplit(chrname,"r",fixed=T))[2]
  }else{
    chrnum = as.integer(unlist(strsplit(chrname,"r",fixed=T))[2])
  }
  # store coverage file (per chromosome)
  print("loading coverage file...")
  dataByChr <- read.delim(paste(inpath,infile,sep=""),header = F,sep = "\t",as.is = T)
  spacing=4
  chrlen = chrom_size[which(as.character(chrom_size$V1)==chrname),2]
  # plot and save to png
  outfile = paste(oupath,chrname,"_",ID,"_DSB_distr.png",sep="")
  print("Plotting...")
  png(outfile)
  par(mar = c(spacing, 8, 2, 2))
  # plotting DSBs, keeping the same y-axis for all
  plot(dataByChr$V2,dataByChr$V4,pch=16,xaxt = "n",
      las = 1,ylab = "DSB counts",xlim = c(1,chrlen),
      cex = 0.25,ylim=c(-50,1500),
      xlab="",main=chrname)
  # plotting karyogram
  pI <- plotIdiogram(chromosome = "1",build = "hg19",unit = "bp",label.y = -0.35,new=F,ylim=c(-15,-40),label.cytoband = F)
  dev.off()
  sprintf("Plotting of %s is done", chrname)
}
