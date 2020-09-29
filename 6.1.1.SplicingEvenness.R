# Script to calculate evenness of splicing and compare to gene expression
# To do the following:
# 1) Load sample and gene expression data
# 2) Load, reformat and filter MISO data
# 3) Calculate evenness and plot all evenness vs expression data
# 4) Output results and summary file
# 5) Plot examples of evenness vs expression
# 6) Plot evenness in normal samples vs tumour samples

# Main code ####

# 0. Prepare environment
setwd('~/Documents/EPICC/Data/Splicing/');library(dplyr);library(data.table)
library(vegan);library(e1071);library(wesanderson)

# 1. Get sample and gene expression data ####
# Get sample list
rnasam <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]
sampleList <- readRDS('miso_list_files.rds')
sampleList <- sampleList[which(sampleList$File!='missing'),];row.names(sampleList) <- c(1:nrow(sampleList))

# Load gene expression
geneexp <- as.data.frame(fread('~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_tpm.txt'))

# 2. Read in MISO data, reformat and save if sample has more than 2000 events ####
keep <- c();misoList <- list()
for(i in c(1:nrow(sampleList))) {
  cursam <- sampleList[i,]
  tmpsum <- as.data.frame(fread(cursam$File))
  tmpsum$event_name <- gsub('(ENSG\\d+)\\.\\d+','\\1',tmpsum$event_name)
  tmpsum <- merge(tmpsum,geneexp[,c('GeneID',cursam$Sample)],by.x='event_name',by.y='GeneID')
  colnames(tmpsum) <- c(colnames(tmpsum[,c(1:ncol(tmpsum)-1)]),'Expression')
  tmpsum <- tmpsum[which(tmpsum$Expression>=10),]
  if(nrow(tmpsum)>=2000) {
    keep <- c(keep,i)
    row.names(tmpsum) <- c(1:nrow(tmpsum))
    tmpsum$LogExp <- log(tmpsum$Expression)
    misoList[[cursam$Sample]] <- tmpsum
  }
  print(paste0('Loaded and filtered ',cursam$Sample,': ',nrow(tmpsum),' events'))
}
sampleList <- sampleList[keep,];row.names(sampleList) <- c(1:nrow(sampleList))

# Save miso list
saveRDS(misoList,file='~/Documents/ThesisOther/ScriptsForThesis/RData/miso_list.rds')

# 3. Calculate metrics and plot evenness vs expression and evenness histogram ####
misoList <- readRDS(file='~/Documents/ThesisOther/ScriptsForThesis/RData/miso_list.rds')
pdf('~/Documents/EPICC/Data/Splicing/AllSamplePlots.pdf',height=14,width=10)
misodf <- data.frame(Sample=sampleList$Sample,Patient=gsub('(C\\d+)_\\S+','\\1',sampleList$Sample),Type='Tumour',NumEvents=0,MeanEven=0,MedianEven=0,MeanExp=0,MedianExp=0,Cor=0,pCor=0,Skew=0)
misodf[grep('^C\\d+_E',misodf$Sample),'Type'] <- 'Normal'
par(mar=c(5,5,3,1),font=2);annoList <- list()
layout(matrix(c(1:8),ncol=2,byrow=T))
for(i in c(1:nrow(sampleList))) {
  cursam <- sampleList[i,]
  shans <- c();numisos <- c();maxshans <- c();psis <- c();evens <- c()
  sam <- misoList[[cursam$Sample]]
  for(j in c(1:nrow(sam))) {
    mis <- sam[j,]
    means <- as.numeric(strsplit(mis$miso_posterior_mean,",")[[1]])
    if(length(means)==1) {
      means <- c(means,1-means)
    }
    psis <- c(psis,paste0(means,collapse=','))
    numiso <- length(means)
    shan <- diversity(means, index = "shannon", MARGIN = 1, base = exp(1));maxshan <- log(numiso,base=exp(1))
    even <- shan/maxshan
    numisos <- c(numisos,numiso);maxshans <- c(maxshans,maxshan);evens <- c(evens,even)
    shans <- c(shans,shan)
  }
  sam$PSI <- psis
  sam$Shannon <- shans;sam$nIso <- numisos
  sam$nIso <- numisos;sam$MaxShannon <- maxshans;sam$Evenness <- evens
  annoList[[cursam$Sample]] <- sam
  print(paste0('Analysed ',cursam$Sample,': ',nrow(sam),' events'))
  
  misodf[i,'NumEvents'] <- nrow(sam)
  misodf[i,'MeanEven'] <- mean(sam$Evenness);misodf[i,'MeanExp'] <- mean(sam$LogExp)
  misodf[i,'MedianEven'] <- median(sam$Evenness);misodf[i,'MedianExp'] <- median(sam$LogExp)
  rescor <- suppressWarnings(cor.test(sam$Evenness,sam$LogExp,method='spearman'))
  misodf[i,'Cor'] <- rescor$estimate;misodf[i,'pCor'] <- rescor$p.value
  misodf[i,'Skew'] <- skewness(sam$Evenness)
  
  hist(sam$Evenness,xlim=c(0,1),col='skyblue3',xlab='',ylab='',axes=F,breaks=30,border=NA,main='')
  abline(v=mean(sam$Evenness),lty=2,col=scales::alpha('firebrick3',0.65),lwd=2)
  axis(side=1,font=2,cex.axis=1.6);mtext(side=1,line=3,text='Evenness (H/Hmax)',cex=1.3,font=2)
  axis(side=2,font=2,cex.axis=1.6,las=2);mtext(side=2,line=3,text='Frequency',cex=1.3,font=2)
  mtext(side=3,line=0,text=paste0(gsub('_',' ',cursam$Sample),': Skewness=',signif(skewness(sam$Evenness),2)),cex=1.3,font=2)
  
  plot(sam$Evenness,sam$LogExp,pch=16,xlab='',ylab='',col=scales::alpha('skyblue3',0.4),axes=F,xlim=c(0,1),ylim=c(2,10))
  abline(lm(sam$LogExp~sam$Evenness),lty=2,col=scales::alpha('firebrick3',0.65),lwd=2)
  axis(side=1,font=2,cex.axis=1.6);mtext(side=1,line=3,text='Evenness (H/Hmax)',cex=1.3,font=2)
  axis(side=2,font=2,cex.axis=1.6,las=2);mtext(side=2,line=3,text='Expression (logTPM)',cex=1.3,font=2)
  mtext(side=3,line=0,text=paste0(gsub('_',' ',cursam$Sample),': rho=',signif(rescor$estimate,2),', p=',formatC(rescor$p.value, format = "e", digits = 2)),cex=1.3,font=2)
}
dev.off()
misodf$fdrCor <- p.adjust(misodf$pCor,method='fdr')

# 4. Output summary file and annotated data ####
saveRDS(misodf,file='~/Documents/ThesisOther/ScriptsForThesis/RData/evenness_results.rds')
saveRDS(annoList,file='~/Documents/ThesisOther/ScriptsForThesis/RData/miso_anno_list.rds')

# 5. Plot examples of evenness vs expression data ####
misodf <- readRDS(file='~/Documents/ThesisOther/ScriptsForThesis/RData/evenness_results.rds')
annoList <- readRDS(file='~/Documents/ThesisOther/ScriptsForThesis/RData/miso_anno_list.rds')

pdf('~/Documents/Thesis/figures/Splicing_Median_Dists.pdf',width = 10,height=10)
figsam <- c('C559_A1_G5','C559_E1_B1')
par(mar=c(5,5.5,4,2),font=2)
layout(matrix(c(1:4),ncol=2,byrow=T))
sam <- annoList[[figsam[1]]];curdf <- misodf[which(misodf$Sample==figsam[1]),]
plot(sam$Evenness,sam$LogExp,pch=16,xlab='',ylab='',col=scales::alpha('skyblue3',0.3),axes=F,xlim=c(0,1),ylim=c(2,10))
abline(lm(sam$LogExp~sam$Evenness),lty=2,col=scales::alpha('firebrick3',0.65),lwd=2)
axis(side=1,font=2,cex.axis=1.6);mtext(side=1,line=3,text='Evenness (H/Hmax)',cex=1.3,font=2)
axis(side=2,font=2,cex.axis=1.6,las=2);mtext(side=2,line=3,text='Expression (logTPM)',cex=1.3,font=2)
mtext(side=3,line=-0.5,text=paste0(gsub('_',' ',figsam[1]),': rho=',signif(curdf$Cor,2),', adj-p=',formatC(curdf$fdrCor, format = "e", digits = 2)),cex=1.3,font=2)
text(x=-0.1,y=11.2,label='(a)',font=2,xpd=T,cex=2.5)

hist(sam$Evenness,xlim=c(0,1),col='skyblue3',xlab='',ylab='',axes=F,breaks=30,border=NA,main='')
abline(v=median(sam$Evenness),lty=2,col=scales::alpha('firebrick3',0.65),lwd=2)
axis(side=1,font=2,cex.axis=1.6);mtext(side=1,line=3,text='Evenness (H/Hmax)',cex=1.3,font=2)
axis(side=2,font=2,cex.axis=1.6,las=2);mtext(side=2,line=3.5,text='Frequency',cex=1.3,font=2)
mtext(side=3,line=0,text=paste0(gsub('_',' ',figsam[1]),': Median=',signif(median(sam$Evenness),2)),cex=1.3,font=2)
text(x=-0.1,y=890,label='(b)',font=2,xpd=T,cex=2.5)

plot(sam$Evenness,sam$LogExp,pch=16,xlab='',ylab='',col=scales::alpha('skyblue3',0.3),axes=F,xlim=c(0,1),ylim=c(2,10))
abline(lm(sam$LogExp~sam$Evenness),lty=2,col=scales::alpha('firebrick3',0.65),lwd=2)
axis(side=1,font=2,cex.axis=1.6);mtext(side=1,line=3,text='Evenness (H/Hmax)',cex=1.3,font=2)
axis(side=2,font=2,cex.axis=1.6,las=2);mtext(side=2,line=3.5,text='Expression (logTPM)',cex=1.3,font=2)
mtext(side=3,line=-0.5,text=paste0(gsub('_',' ',figsam[2]),': rho=',signif(curdf$Cor,2),', adj-p=',formatC(curdf$fdrCor, format = "e", digits = 2)),cex=1.3,font=2)
text(x=-0.1,y=11.2,label='(c)',font=2,xpd=T,cex=2.5)

sam <- annoList[[figsam[2]]];curdf <- misodf[which(misodf$Sample==figsam[2]),]
hist(sam$Evenness,xlim=c(0,1),col='skyblue3',xlab='',ylab='',axes=F,breaks=30,border=NA,main='')
abline(v=median(sam$Evenness),lty=2,col=scales::alpha('firebrick3',0.65),lwd=2)
axis(side=1,font=2,cex.axis=1.6);mtext(side=1,line=3,text='Evenness (H/Hmax)',cex=1.3,font=2)
axis(side=2,font=2,cex.axis=1.6,las=2);mtext(side=2,line=3.5,text='Frequency',cex=1.3,font=2)
mtext(side=3,line=0,text=paste0(gsub('_',' ',figsam[2]),': Median=',signif(median(sam$Evenness),2)),cex=1.3,font=2)
text(x=-0.1,y=405,label='(d)',font=2,xpd=T,cex=2.5)
dev.off()

# 6. Plot median evenness normal vs tumour ####
pdf('~/Documents/Thesis/figures/Splicing_MedianEven_Boxplot.pdf')
par(mar=c(4.6,5.1,3.6,0.6),font.lab=2,cex.axis=1.5,font.axis=2,font=2)
boxplot(MedianEven~Type,data=misodf,las=1,frame=F,ylim=c(0,1),xlab='',ylab='',names=paste0(c('Normal','Tumour'),' (n=',as.numeric(table(misodf$Type)),')'),col=wes_palette("Royal1")[c(1,2)],boxcol=wes_palette("Royal1")[c(1,2)],outpch=20,border=wes_palette("BottleRocket1")[c(4,3)],staplecol=wes_palette("BottleRocket1")[c(4,3)],staplelwd=3,outcol=wes_palette("BottleRocket1")[c(4,3)],medcol=wes_palette("BottleRocket1")[c(4,3)],medlwd=3,cex.lab=1.8,cex.axis=1.8)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="#EBEBEB",border=NA)
grid(lty=1,col='white');par(new=T)
boxplot(MedianEven~Type,data=misodf,las=1,frame=F,ylim=c(0,1),xlab='',ylab='',names=paste0(c('Normal','Tumour'),' (n=',as.numeric(table(misodf$Type)),')'),col=wes_palette("Royal1")[c(1,2)],boxcol=wes_palette("Royal1")[c(1,2)],outpch=20,border=wes_palette("BottleRocket1")[c(4,3)],staplecol=wes_palette("BottleRocket1")[c(4,3)],staplelwd=3,outcol=wes_palette("BottleRocket1")[c(4,3)],medcol=wes_palette("BottleRocket1")[c(4,3)],medlwd=3,cex.lab=1.8,cex.axis=1.8)
segments(1,1.1,x1=2,lwd=3,xpd=T);segments(c(1,2),1.1,y1=1.05,lwd=3,xpd=T);text(x=1.5,y=1.15,labels = "*",cex=1.5,xpd=T)
mtext(side=2,'Median splicing evenness',line=3,cex=2)
mtext(side=1,text='Sample type',line=3,cex=2,font=2)
dev.off()



