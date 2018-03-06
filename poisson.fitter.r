args <- commandArgs(trailingOnly=TRUE)
inputdir=args[1]
bedfile=args[2]
experiment=args[3]
pvalue=args[4]
outdir=args[5]

setwd(inputdir)
data = read.delim(bedfile,header = FALSE)
lamb = mean(data$V4)
setwd(outdir)
data$V5= rpois(length(data$V4),lamb)
hr=hist(data$V5)
jpeg(paste(experiment,'_break_distr.jpg',sep=''))
par(mfrow=c(1,2))
mxbreak = max(hr$breaks)
ho = hist(data$V4[data$V4<mxbreak],breaks = seq(0,mxbreak),xlab="break number per location",main=paste("Poisson distribution\n(lambda=",toString(lamb),')',sep=''))
hr = hist(data$V5,breaks = seq(0,mxbreak),xlab="break number per location",main=experiment)
dev.off()
jpeg(paste(experiment,'_poisson_fitting.jpg',sep=''))
plot(log(hr$counts),log(ho$counts),ylab='actual break counts (log)',xlab='fitted break counts (log)')
dev.off()
# the threshold is 11, because there is nothing in the created dataset that is bigger than 11
# so the chance to have in a value greater than 11 in this number of events is basically 0
# the poisson distribution has been used, because it is a counting stuff and it represents it fairly nicely
p05 = qpois(p=as.numeric(pvalue),lambda = lamb)
enrichment = matrix(data=NA,nrow=144,ncol=2)
k=0
for (i in p05:max(data$V4)){
  if (length(which(data$V4==i))>0){
    k = k+1
    enrichment[k,1] = i
    enrichment[k,2] = length(which(data$V4==i))/(length(which(data$V5==i))+length(which(data$V4==i)))*100
  }
}
jpeg(paste(experiment,'_',pvalue,'_enrichment.jpg',sep=''))
plot(enrichment[,2]~enrichment[,1],ylab='enrichment (%)',xlab='number of breaks')
dev.off()

#everything is a peak that is above this threshold
#saving it to a bed file (with the number of breaks, so actually it's a bedgraph)
peaks=matrix(data=NA,nrow=dim(data)[1],ncol=4)
peaks=data[which(data$V4>p05),1:4]
write.table(peaks,paste(experiment,'_',pvalue,'_peaks.bed',sep=''),append=FALSE,quote=FALSE,sep="\t",row.names = FALSE,col.names = FALSE,na="")
