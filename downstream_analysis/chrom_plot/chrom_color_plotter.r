# to compare results from replicates

inpath = ''
outpath = ''
chromsize = ''

idList = c("ID1","ID2","ID3","ID4") # or less or more, up to 6 at the moment

print("loading chromosome file...")
chrom_size <- read.delim(chromsize,header = F)

infile='chr1_allreplicates.bed'

spacing=4
i=0

chrname = unlist(strsplit(infile,'_',fixed=T))[1]
print(chrname)

# assigning colors to different datasets
dsbCol = c("darkgreen","red3","gold","steelblue3","khaki4","purple")
# naming the colors by the IDs
names(dsbCol) = idList

print("loading coverage file...")
dataByChr <- read.delim(paste(inpath,infile,sep=""),header = F,sep = "\t",as.is = T)

# plot all data from the same chromosome on the same plot, with different color respectively to the dataset
pdf('test.pdf')
plot(dataByChr$V2,dataByChr$V4,pch=16,xaxt = "n",
     las = 1,ylab = "DSB counts",xlim = c(1,chrlen),
     col= dsbCol[dataByChr$V5],cex = 0.25,ylim=c(0,max(dataByChr$V4)),
     xlab="",main=chrname)
dev.off()
