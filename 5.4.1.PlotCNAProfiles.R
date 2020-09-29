# Script to plot CNA profiles of all WGS and lpWGS samples together
# To do the following:
# 1) Load chromosome and sample information
# 2) Load CNA and ploidy data for all samples
# 3) Plot WGS CNA profiles
# 4) Plot lpWGS CNA profiles

# Functions ####
# Classify WGS CNAs in terms of gain/loss
classifySeq <- function(segs) {
  segs <- segs[which(segs$chromosome %in% paste0('chr',c(1:22))),]
  segs <- segs[which(!is.na(segs$A) & !is.na(segs$B)),]
  segs$CNA <- 'Dip';segs[segs$CNt==3,'CNA'] <- 'Gain';segs[segs$CNt>3,'CNA'] <- 'HighGain'
  segs[segs$CNt==1,'CNA'] <- 'Mono';segs[segs$CNt<1,'CNA'] <- 'Loss';segs[segs$A==2 & segs$B==0,'CNA'] <- 'cnLOH'
  segs <- segs[,c("chromosome","start.pos","end.pos","CNt","CNA")]
  return(segs)
}
# Classify lpWGS CNAs in terms of gain/loss
classifyqdna <- function(qdna,bins) {
  segs <- bins[,c("chromosome","start.pos","end.pos")]
  segs$CNt <- qdna[,1] 
  segs$CNA <- 'Dip';segs[segs$CNt==3,'CNA'] <- 'Gain';segs[segs$CNt>3,'CNA'] <- 'HighGain'
  segs[segs$CNt==1,'CNA'] <- 'Mono';segs[segs$CNt<1,'CNA'] <- 'Loss'
  return(segs)
}

# Main code ####

# 0. Prepare environment
setwd("~/Documents/EPICC/Data/CopyNumber");library(readr);library(ggplot2);library(gridExtra)
library(data.table);library(stringr);'%ni%' <- Negate('%in%')
regcols <- c(A='#E31A1C',B='#377DB8',C='#4DAE49',D='#904A9A')

# 1. Get chromosome plotting and sample information ####
chroms <- paste0('chr',c(1:22))
chrInfo <- read.csv("centromerePositions_hg38.csv");chrInfo <- chrInfo[c(1:22),]
matTab <- data.frame(chrInfo[c("chromosome","plotStart")])

glandPass <- readRDS('~/Documents/EPICC/Data/QTLs/nodesamples/wgs_used_in_phylo.rds')
wgssam <- gsub('EPICC_(C\\d+_\\S\\d+_\\S\\d+)_.+','\\1',glandPass)
wgssam <- wgssam[grep('C\\d+_(A|B|C|D)\\d+_G\\d+',wgssam)]
keep <- c()
for(i in c(1:length(wgssam))) {
  if(file.exists(paste0('~/Documents/EPICC/Data/CopyNumber/Segments/EPICC_',wgssam[i],'_D1_GRCh38_segments.txt'))) {
    keep <- c(keep,i)
  }
}
wgssam <- wgssam[keep]

lpsam <- read.table('~/Documents/EPICC/Data/CopyNumber/LowPass/lp_samplelist.txt',header=T)
lpsam <- lpsam[grep('C\\d+_(A|B|C|D)\\d+_G\\d+',lpsam$Sample),];row.names(lpsam) <- c(1:nrow(lpsam))

samples <- sort(unique(c(wgssam,lpsam$Sample)))
patients <- as.factor(gsub('(C\\d+).+','\\1',samples))

# 2. Load CNA and ploidy data for both WGS (Sequenza) and lpWGS (QDNAseq) ####

# Get bin info
segments <- read.table('~/Documents/EPICC/Data/CopyNumber/LowPass/Segments/C516_EPICC_C516_A1_G3_L1_500kb_GRCh38_multiregion_segmentation_calls.txt')
rowbins <- unlist(str_split(row.names(segments),':'));rowpos <- unlist(str_split(rowbins[seq.int(2,nrow(segments)*2,by=2)],'-'))
bins <- data.frame(chromosome=paste0('chr',rowbins[seq.int(1,nrow(segments)*2,by=2)]),start.pos=as.integer(rowpos[seq.int(1,nrow(segments)*2,by=2)]),end.pos=as.integer(rowpos[seq.int(2,nrow(segments)*2,by=2)]))

# Add classification of CNAs to segment dataframes
ploidys <- c()
for(i in c(1:length(samples))) {
  sam <- samples[i]
  if(sam %in% wgssam) {
    segments <- read.table(paste0('~/Documents/EPICC/Data/CopyNumber/Segments/EPICC_',sam,'_D1_GRCh38_segments.txt'),header=T)
    segmuts <- classifySeq(segments)
    assign(paste0(samples[i],'Segs'),segmuts)
    candp <- read.table(paste0('~/Documents/EPICC/Data/CopyNumber/CellandPloidy/EPICC_',samples[i],'_D1_GRCh38_confints_CP.txt'),header=T)
    ploidys <- c(ploidys,candp$ploidy.estimate[2])
  } else {
    segments <- read.table(paste0('~/Documents/EPICC/Data/CopyNumber/LowPass/Segments/',as.character(patients[i]),'_EPICC_',sam,'_L1_500kb_GRCh38_multiregion_segmentation_calls.txt'))
    segmuts <- classifyqdna(segments,bins)
    assign(paste0(samples[i],'Segs'),segmuts)
    ploidys <- c(ploidys,lpsam[which(lpsam$Sample==samples[i]),'Ploidy'])
  }
}
patploid <-colMeans(boxplot(ploidys~patients,plot=F)$stats)
names(patploid) <- unique(patients)#;patploid <- sort(patploid,decreasing = T)

# 3. Plot WGS CNA profiles ####
ploidys <- c();wgspat <- gsub('(C\\d+)_\\S+','\\1',wgssam)
for(i in c(1:length(wgssam))) {
  sam <- wgssam[i]
  candp <- read.table(paste0('~/Documents/EPICC/Data/CopyNumber/CellandPloidy/EPICC_',wgssam[i],'_D1_GRCh38_confints_CP.txt'),header=T)
  ploidys <- c(ploidys,candp$ploidy.estimate[2])
}
patploid <-colMeans(boxplot(ploidys~wgspat,plot=F)$stats)
names(patploid) <- unique(wgspat)

sortsam <- c();for(i in c(1:length(patploid))) { sortsam <- c(sortsam,wgssam[which(grepl(names(patploid)[i],wgssam))]) }
patients <- as.factor(gsub('(C\\d+).+','\\1',sortsam));patCol <- patients;levels(patCol) <- 1:length(levels(patCol))
patCol <- as.numeric(patCol)
regs <- gsub('C\\d+_(\\S).+','\\1',sortsam);depth <- rep('WGS',length(sortsam));depth[which(sortsam %ni% wgssam)] <- 'LP'

pdf('~/Documents/EPICC/Data/CopyNumber/LowPass/onlydeep.pdf',width=15,height=20)
topy <- 277.4
plot(x=c(0,3.2e9),y=c(0,topy+0.5), col='white',xaxt='n',yaxt='n',xlab='Chromosomes',ylab=paste0('Samples (n=',length(sortsam),')'),bty='n',yaxs='i',xaxs='i',main='',cex.lab=1,font.lab=2)
ypos <- topy+1;patpos <- c();gappos <- c()
for(i in c(1:length(sortsam))) {
  sample <- sortsam[i];segmuts <- get(paste0(sample,'Segs'))
  if(i!=1) {
    if(patCol[i]!=patCol[i-1]) {
      ypos <- ypos-3;gappos <- c(gappos,ypos)
      patpos <- c(patpos,ypos-.5-length(which(gsub('(C\\d+).+','\\1',sortsam)==patients[i]))/2)
    }
  } else {
    patpos <- c(patpos,ypos-1-length(which(gsub('(C\\d+).+','\\1',sortsam)==patients[i]))/2)
  }
  
  ypos <- ypos-1
  # Convert start and end positions to genomic coordinates
  chrIDs <- segmuts$chromosome
  gPositions <- matTab[['plotStart']][match(chrIDs,matTab$chromosome)]
  startPos <- segmuts$start.pos+gPositions
  endPos <- segmuts$end.pos+gPositions
  
  # Select colours to plot gains, dip and loss
  colour = rep(NA,length=length(sortsam))
  colour[which(segmuts$CNA=='HighGain')] = 'red4'
  colour[which(segmuts$CNA=='Gain')] = 'red1'
  colour[which(segmuts$CNA=='Dip')] = 'gray90'
  colour[which(segmuts$CNA=='cnLOH')] = 'mediumpurple'
  colour[which(segmuts$CNA=='Mono')] = 'steelblue1'
  colour[which(segmuts$CNA=='Loss')] = 'darkblue'
  
  # Plot segments
  rect(xleft=startPos, xright=endPos, ybottom=(ypos-.4), ytop=(ypos+.4),col=colour,border=NA)
  
  rect(xleft=2.9e9,xright=2.93e9,ybottom=(ypos-.5), ytop=(ypos+.5),col=regcols[regs[i]],border=NA)
}
axis(side=1,at=c(chrInfo$plotCentroStart,2.915e9),labels=c(gsub('chr(\\S+)','\\1',chrInfo$chromosome),'Region'),cex.axis=0.8,las=1,font=2,line=0)
axis(side=2,at=patpos,labels=unique(patients),cex.axis=.8,las=1,font=2,line=0)
abline(v=chrInfo$plotCentroStart,lty=2,lwd=0.3);abline(v=chrInfo$plotCentroEnd,lty=2,lwd=0.3)
abline(v=chrInfo$plotStart,lty=1,lwd=0.6);abline(v=chrInfo$plotEnd,lty=1,lwd=0.6)
rect(xleft=0,xright=max(chrInfo$plotEnd),ytop=topy+0.52,ybottom=topy+0.48,lwd=0.6)
rect(xleft=0,xright=max(chrInfo$plotEnd),ytop=0.02,ybottom=-0.02,lwd=0.6)
rect(xleft=1e7,xright=max(chrInfo$plotEnd)-1e7,ytop=gappos+2.5,ybottom=gappos-.5,col='white',border=NA)
rect(xleft=0,xright=max(chrInfo$plotEnd),ytop=gappos+2.6,ybottom=gappos-.6,lwd=0.3)
par(xpd=T,font=2)
legend(x=3e9,y=topy,legend=c('4+','3','2','2:0','1','0'),fill=c('red4','red1','gray90','mediumpurple','steelblue2','darkblue'),border=NA,cex=1,box.lwd=0.5,title='CNA')
legend(x=3e9,y=topy-32,legend=c('A','B','C','D'),fill=regcols,border=NA,cex=1.14,box.lwd=0.5,title='Region')
par(xpd=F,font=1)
dev.off()

# 4. Plot lpWGS CNA profiles ####
ploidys <- lpsam$Ploidy
patploid <-colMeans(boxplot(ploidys~lppat,plot=F)$stats)
names(patploid) <- unique(lppat)

sortsam <- c();for(i in c(1:length(patploid))) { sortsam <- c(sortsam,lpsam$Sample[which(grepl(names(patploid)[i],lpsam$Sample))]) }
patients <- as.factor(gsub('(C\\d+).+','\\1',sortsam));patCol <- patients;levels(patCol) <- 1:length(levels(patCol))
patCol <- as.numeric(patCol)
regs <- gsub('C\\d+_(\\S).+','\\1',sortsam);depth <- rep('LP',length(sortsam))

pdf('~/Documents/EPICC/Data/CopyNumber/LowPass/onlylowpass.pdf',width=15,height=20)
topy <- 382.4
plot(x=c(0,3.2e9),y=c(0,topy+0.5), col='white',xaxt='n',yaxt='n',xlab='Chromosomes',ylab=paste0('Samples (n=',length(sortsam),')'),bty='n',yaxs='i',xaxs='i',main='',cex.lab=1,font.lab=2)
ypos <- topy+1;patpos <- c();gappos <- c()
for(i in c(1:length(sortsam))) {
  sample <- sortsam[i];curpat <- gsub('(C\\d+)_\\S+','\\1',sample)
  segments <- read.table(paste0('~/Documents/EPICC/Data/CopyNumber/LowPass/Segments/',curpat,'_EPICC_',sample,'_L1_500kb_GRCh38_multiregion_segmentation_calls.txt'))
  segmuts <- classifyqdna(segments,bins)
  if(i!=1) {
    if(patCol[i]!=patCol[i-1]) {
      ypos <- ypos-3;gappos <- c(gappos,ypos)
      patpos <- c(patpos,ypos-.5-length(which(gsub('(C\\d+).+','\\1',sortsam)==patients[i]))/2)
    }
  } else {
    patpos <- c(patpos,ypos-1-length(which(gsub('(C\\d+).+','\\1',sortsam)==patients[i]))/2)
  }
  
  ypos <- ypos-1
  # Convert start and end positions to genomic coordinates
  chrIDs <- segmuts$chromosome
  gPositions <- matTab[['plotStart']][match(chrIDs,matTab$chromosome)]
  startPos <- segmuts$start.pos+gPositions
  endPos <- segmuts$end.pos+gPositions
  
  # Select colours to plot gains, dip and loss
  colour = rep(NA,length=length(sortsam))
  colour[which(segmuts$CNA=='HighGain')] = 'red4'
  colour[which(segmuts$CNA=='Gain')] = 'red1'
  colour[which(segmuts$CNA=='Dip')] = 'gray90'
  colour[which(segmuts$CNA=='cnLOH')] = 'mediumpurple'
  colour[which(segmuts$CNA=='Mono')] = 'steelblue1'
  colour[which(segmuts$CNA=='Loss')] = 'darkblue'
  
  # Plot segments
  rect(xleft=startPos, xright=endPos, ybottom=(ypos-.4), ytop=(ypos+.4),col=colour,border=NA)
  
  rect(xleft=2.9e9,xright=2.93e9,ybottom=(ypos-.5), ytop=(ypos+.5),col=regcols[regs[i]],border=NA)
}
axis(side=1,at=c(chrInfo$plotCentroStart,2.915e9),labels=c(gsub('chr(\\S+)','\\1',chrInfo$chromosome),'Region'),cex.axis=0.8,las=1,font=2,line=0)
axis(side=2,at=patpos,labels=unique(patients),cex.axis=.8,las=1,font=2,line=0)
abline(v=chrInfo$plotCentroStart,lty=2,lwd=0.3);abline(v=chrInfo$plotCentroEnd,lty=2,lwd=0.3)
abline(v=chrInfo$plotStart,lty=1,lwd=0.6);abline(v=chrInfo$plotEnd,lty=1,lwd=0.6)
rect(xleft=0,xright=max(chrInfo$plotEnd),ytop=topy+0.52,ybottom=topy+0.48,lwd=0.6)
rect(xleft=0,xright=max(chrInfo$plotEnd),ytop=0.02,ybottom=-0.02,lwd=0.6)
rect(xleft=1e7,xright=max(chrInfo$plotEnd)-1e7,ytop=gappos+2.5,ybottom=gappos-.5,col='white',border=NA)
rect(xleft=0,xright=max(chrInfo$plotEnd),ytop=gappos+2.6,ybottom=gappos-.6,lwd=0.3)
par(xpd=T,font=2)
legend(x=3e9,y=topy,legend=c('4+','3','2','2:0','1','0'),fill=c('red4','red1','gray90','mediumpurple','steelblue2','darkblue'),border=NA,cex=1,box.lwd=0.5,title='CNA')
legend(x=3e9,y=topy-37,legend=c('A','B','C','D'),fill=regcols,border=NA,cex=1.14,box.lwd=0.5,title='Region')
par(xpd=F,font=1)
dev.off()





