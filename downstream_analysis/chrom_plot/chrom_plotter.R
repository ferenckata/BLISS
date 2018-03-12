# plot DSB data in sliding windows
# http://bioconductor.org/packages/release/bioc/vignettes/TitanCNA/inst/doc/TitanCNA.pdf
# based on library(TitanCNA)
# library(argparser)

inpath = ''
chromsize = ''
# store chromosome sizes
print("loading chromosome file...")
chrom_size <- read.delim(chromsize,header = F)
for(infile in list.files(inpath)){
  chrname = unlist(strsplit(infile,'_',fixed=T))[1]
  print(chrname)
  print("loading coverage file...")
  dataByChr <- read.delim(paste(inpath,infile,sep=""),header = F,sep = "\t",as.is = T)
  spacing=4
  chrlen = chrom_size[which(as.character(chrom_size$V1)==chrname),2]
  outfile = paste(inpath,chrname,"_DSB_distr.png",sep="")
  print("Plotting...")
  png(outfile,width=10,height=6)
  par(mar = c(spacing, 8, 2, 2))
  plot(dataByChr$V2,dataByChr$V4,pch=16,xaxt = "n",
      las = 1,ylab = "DSB counts",xlim = c(1,chrlen),
      cex = 0.25,ylim=c(0,max(dataByChr$V4)),
      xlab="",main=chrname)
  dev.off()
  sprintf("Plotting of % is done", chrname)
}
