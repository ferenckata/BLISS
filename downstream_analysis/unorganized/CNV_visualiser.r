# do something with the tcga ovarian cancer cnv dataset
# downloaded from TCGA database on 27th of April, 2018

setwd("C:/Users/GG/Desktop/Kata/KI/vancouver/tcga/cnv/")
dirs = list.dirs(path=".",recursive = F)

plot(0,0,ylim = c(-10,10),xlim = c(0,25000000))
for(d in 1:length(dirs)){
  print(dirs[d])
  cnvfile = list.files(path = dirs[d],pattern = ".seg.txt")
  cnv_data = read.delim(paste(dirs[d],"/",cnvfile,sep=""),
                      sep = "\t",header = T,as.is = T)
  max(cnv_data$Segment_Mean)
  min(cnv_data$Segment_Mean)

  dataByChr = cnv_data[which(cnv_data$Chromosome=="1"),]
  for(i in 1:dim(dataByChr)[1]){
      points(dataByChr$Start[i+1],dataByChr$Segment_Mean[i+1]-dataByChr$Segment_Mean[i]) 
    if(dataByChr$Segment_Mean[i]>=-1 && dataByChr$Segment_Mean[i] <1){
      dataByChr$Color[i] = "green"
    }else if(dataByChr$Segment_Mean[i]>=1 && dataByChr$Segment_Mean[i] <2){
      dataByChr$Color[i] = "yellow"
    }else if(dataByChr$Segment_Mean[i]>=2 && dataByChr$Segment_Mean[i] <3){
      dataByChr$Color[i] = "orange"
    }else if(dataByChr$Segment_Mean[i]>=3){
      dataByChr$Color[i] = "red"
    }else if(dataByChr$Segment_Mean[i]>=-2 && dataByChr$Segment_Mean[i]<(-1)){
      dataByChr$Color[i] = "blue"
    }else if(dataByChr$Segment_Mean[i]<(-2)){
      dataByChr$Color[i] = "purple"
    }
  }
  #if(a==1){
   # plot(dataByChr$Start,dataByChr$Segment_Mean,type='p',xaxt = "n",
    #     las = 1,ylab = "Segment Mean",col=dataByChr$Color,
     #    cex = 1,xlab="",main="chr1")
  #}else{
   # points(dataByChr$Start,dataByChr$Segment_Mean,type='p',xaxt = "n",
    #       las = 1,col=dataByChr$Color)
 # }
  a = 2
}

