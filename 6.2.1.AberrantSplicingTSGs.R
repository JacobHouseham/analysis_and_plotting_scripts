# Script to assess aberrant splicing and whether it is more common in tumour suppressor genes (TSGs)
# To do the following:
# 1) Load expression data and pre-filter MISO compare data
# 2) Load, reformat and save miso_compare data
# 3) Assess if more aberrant splicing occurs in TSGs
# 4) Plots results of TSG analysis

# Functions ####
# Read in a comparison and output the significantly differentially spliced isoforms
readinandfilter <- function(sam1,sam2) {
  comp <- read.delim(paste0("~/Documents/Splicing/allCompare/",sam1,"_vs_",sam2,".miso_bf"),sep="\t")
  colnames(comp) <- c('GeneID','Sample1Mean','Sample1CILow','Sample1CIHigh','Sample2Mean','Sample2CILow','Sample2CIHigh','Diff','BayesFactor','Isoforms','Sample1Counts','Sample1Assigned','Sample2Counts','Sample2Assigned','Chr','Strand','Start','End')
  
  maxdiff <- c();maxbayes <- c()
  for(i in c(1:nrow(comp))) {
    curgen <- comp[i,]
    cdiff <- abs(as.numeric(strsplit(curgen$Diff,',')[[1]]))
    maxdiff <- c(maxdiff,max(abs(cdiff)))
    
    cbayes <- abs(as.numeric(strsplit(curgen$BayesFactor,',')[[1]]))
    maxbayes <- c(maxbayes,max(cbayes))
  }
  comp$MaxDiff <- maxdiff;comp$MaxBayes <- maxbayes
  comp$GeneID <- gsub('(ENSG\\d+).\\d+','\\1',comp$GeneID)
  return(comp)
}

# Main code ####

# 0. Prepare environment
library(gplots);library(data.table);setwd("~/Documents/Splicing")
library(moments);"%ni%" <- Negate("%in%");library(caroline);library(wesanderson)
geneinfo <- read.table("~/Documents/EPICC/Data/Expression/compiledGeneInfo.txt",header=T)
genehg38 <- merge(geneinfo,EPICC,by='GeneID')
cosmic <- read.csv('~/Documents/EPICC/Data/Mutations/NewCGC_Tier1.csv')
tsgs <- cosmic[grep('TSG',cosmic$Role.in.Cancer),'Gene.Symbol']

# 1. Load gene expression and find the MISO comparisons that have at least 1000 events ####
EPICC <- as.data.frame(fread(paste0("~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_tpm.txt"),header=T))
samlist <- read.table(file='~/Documents/Splicing/NewSplicing/sample_list.txt',header=T)
rnasam <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]
samlist <- samlist[which(samlist$Samples %in% rnasam & samlist$MatchedNorm %in% rnasam),]
nrows <- c() # Only analyse samples with a decent number of events
for(i in c(1:nrow(samlist))) {
  tmp <- read.delim(paste0("~/Documents/Splicing/allCompare/",samlist[i,'MatchedNorm'],"_vs_",samlist[i,'Samples'],".miso_bf"),sep="\t")
  nrows <- c(nrows,nrow(tmp))
}
samlist <- samlist[which(nrows>=1000),];row.names(samlist) <- c(1:nrow(samlist))

allsam <- sort(unique(c(samlist$Samples,samlist$MatchedNorm)))
patients <- sort(unique(gsub('(C\\d+)_.+','\\1',samlist$MatchedNorm)))
EPICC <- EPICC[,c('GeneID',allsam)]

# 2. Load MISO comparison data and record the events in TSGs ####
allevents <- c();usesam <- c();spliceList <- tsgList <- list()
for(i in c(1:nrow(samlist))) {
  curcom <- samlist[i,]
  expsams <- genehg38[,c('GeneID','Name',curcom$MatchedNorm,curcom$Samples)]
  cursig <- readinandfilter(curcom$MatchedNorm,curcom$Samples)
  cursig <- merge(cursig,expsams,by='GeneID')
  cursig <- cursig[which(cursig$MaxDiff>0.2 & cursig$MaxBayes>20 & cursig[,curcom$MatchedNorm]>=10 & cursig[,curcom$Samples]>=10),];rownames(cursig) <- c(1:nrow(cursig))
  cursigTSG <- cursig[which(cursig$Name %in% tsgs),]
  if(nrow(cursig)>=1 & nrow(cursigTSG)>=1) {
    print(paste0('Sample ',curcom$Samples,': ',nrow(cursig),' sig events (',nrow(cursigTSG),' in TSGs)'))
    rownames(cursigTSG) <- c(1:nrow(cursigTSG))
    spliceList[[curcom$Samples]] <- cursig
    tsgList[[curcom$Samples]] <- cursigTSG
    allevents <- sort(unique(c(allevents,cursig$GeneID)));usesam <- c(usesam,curcom$Samples)
  }
}

save(samlist,allevents,usesam,spliceList,tsgList,file='~/Documents/ThesisOther/ScriptsForThesis/RData/aberrant_splicing_data.Rdata')

# 3. Assess if TSGs are more commonly differentially spliced than all other genes? ####
load(file='~/Documents/ThesisOther/ScriptsForThesis/RData/aberrant_splicing_data.Rdata')
minexp <- 10;pvals <- c();chis <- c()
asRes <- matrix(nrow=2,ncol=nrow(samlist));row.names(asRes) <- c('non-TSGs','TSGs');colnames(asRes) <- samlist$Samples
for(i in c(1:nrow(samlist))) {
  sam <- samlist$Samples[i];normal <- samlist$MatchedNorm[i]
  cursig <- spliceList[[sam]]
  expgenes <- genehg38[which(genehg38[,sam]>=minexp & genehg38[,normal]>=minexp),];exptsgs <- tsgs[which(tsgs %in% expgenes$Name)]
  expnontsgs <- genehg38[which(genehg38[,sam]>=minexp & genehg38$Name %ni% tsgs & genehg38[,normal]>=minexp),]
  astsgs <- nrow(cursig[which(cursig$Name %in% tsgs),]);asnontsgs <- nrow(cursig[which(cursig$Name %ni% tsgs),])
  nonastsgs <- length(which(exptsgs %ni% cursig$Name));nonasnotsgs <- length(which(expnontsgs$Name %ni% cursig$Name))
  
  tab <- rbind(c(nonasnotsgs,nonastsgs),c(asnontsgs,astsgs))
  row.names(tab) <- c('Not AS','AS');colnames(tab) <- c('Not TSG','TSG')
  
  res <- chisq.test(tab)
  pvals <- c(pvals,res$p.value)
  chis <- c(chis,res$statistic)
  
  newtab <- c(tab[2,1]/sum(tab[,1])*100,tab[2,2]/sum(tab[,2])*100)
  asRes[,i] <- newtab
}

adjps <- p.adjust(pvals,method='fdr')
sigval <- sapply(adjps,getstarsNS)

# 4. Plot result of aberrant splicing in TSGs ####
pdf('~/Documents/Thesis/figures/ASinTSGs.pdf',height=8,width=8)
par(mar=c(4.6,5.1,1.1,0.6),font.lab=2,font.axis=2,font=2)
boxplot(asRes[1,],asRes[2,],las=1,frame=F,ylim=c(0,100),xlab='',ylab='',names=c('non-TSGs','TSGs'),col=c(wes_palette("Zissou1")[1],wes_palette("Royal1")[2]),boxcol=c(wes_palette("Zissou1")[1],wes_palette("Royal1")[2]),outpch=20,border=c(wes_palette("BottleRocket2")[3],wes_palette("BottleRocket1")[3]),staplecol=c(wes_palette("BottleRocket2")[3],wes_palette("BottleRocket1")[3]),staplelwd=3,outcol=c(wes_palette("BottleRocket2")[3],wes_palette("BottleRocket1")[3]),medcol=c(wes_palette("BottleRocket2")[3],wes_palette("BottleRocket1")[3]),medlwd=3,cex.lab=1.5,cex.axis=1.5)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="#EBEBEB",border=NA)
grid(lty=1,col='white');par(new=T)
segments(x0=1,y0=asRes[1,],x1=2,y1=asRes[2,],col=scales::alpha('gray40',0.75),lty=2,lwd=1.5)
boxplot(asRes[1,],asRes[2,],las=1,frame=F,ylim=c(0,100),xlab='',ylab='',names=c('non-TSGs','TSGs'),col=c(wes_palette("Zissou1")[1],wes_palette("Royal1")[2]),boxcol=c(wes_palette("Zissou1")[1],wes_palette("Royal1")[2]),outpch=20,border=c(wes_palette("BottleRocket2")[3],wes_palette("BottleRocket1")[3]),staplecol=c(wes_palette("BottleRocket2")[3],wes_palette("BottleRocket1")[3]),staplelwd=3,outcol=c(wes_palette("BottleRocket2")[3],wes_palette("BottleRocket1")[3]),medcol=c(wes_palette("BottleRocket2")[3],wes_palette("BottleRocket1")[3]),medlwd=3,cex.lab=1.5,cex.axis=1.5,add=T)
mtext(side=2,'% of available genes that are abberantly spliced',line=3,cex=1.5)
mtext(side=1,text='Gene type',line=3,cex=1.5,font=2)
dev.off()
