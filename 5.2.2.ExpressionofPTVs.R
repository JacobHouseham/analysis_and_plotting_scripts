# Script for looking at the expression of genes which contain protein truncating variants (PTVs)
# To do the following:
# 1) Load PTV and expression data
# 2) Combine and analyse PTV and expression data
# 3) Plots the results

# Main code ####

# 0. Prepare environment
setwd("~/Documents/EPICC/Data/Mutations");library(dplyr);library(data.table)
library(stringr);'%ni%' <- Negate('%in%');library(sos);library(DESeq2)
cols <- c('Chr','Start','End','Ref','Alt','Location','Gene','Detail','ExonicFunc','AAChange','dbsnpID','cosmicCoding','cosmicNonCoding','Other','Depth','Mut','VAF','Locus','Context')
rnasam <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]

# 1. Load PTV mutations and gene expression ####
stopmat <- readRDS('~/Documents/EPICC/Data/QTLs/allvariants/PTV_mutburdengene_allDNAsamples.rds')
wgssam <- readRDS('~/Documents/EPICC/Data/QTLs/nodesamples/wgs_used_in_phylo.rds')
wgssam <- gsub('EPICC_(C\\d+_\\S\\d+_\\S\\d+)\\S+','\\1',wgssam)
stopmat <- stopmat[,wgssam]
samples <- rnasam[which(rnasam %in% colnames(stopmat))]
stopmat <- stopmat[,samples]
stopmat <- stopmat[which(rowSums(stopmat)>0),]

# Load gene expression
vsd <- readRDS('~/Documents/EPICC/Data/Expression/ProcessedCounts/allgenes.vsd.ensembl.rds')
geneexp <- as.data.frame(assay(vsd));geneexp$GeneID <- row.names(geneexp)
geneexp <- geneexp[,c('GeneID',samples)];row.names(geneexp) <- c(1:nrow(geneexp))

# 2. Analyse the expression of PTVs ####
meanstop <- c();meanwt <- c()
for(i in c(1:nrow(stopmat))) {
  curstop <- as.numeric(stopmat[i,])
  curgene <- as.numeric(geneexp[which(geneexp$GeneID==row.names(stopmat)[i]),samples])
  
  meanstop <- c(meanstop,mean(curgene[which(curstop!=0)]))
  meanwt <- c(meanwt,mean(curgene[which(curstop==0)]))
}

ptvexp <- data.frame(MeanExpression=c(meanwt,meanstop),Type=c(rep('non-PTV',length(meanwt)),rep('PTV',length(meanstop))))
ptvexp$Type <- factor(ptvexp$Type,levels=c('non-PTV','PTV'))

res <- wilcox.test(MeanExpression~Type,alternative='greater',data=ptvexp,paired=T)
formatC(res$p.value, format = "e", digits = 2)

# 3. Plot the expression of PTVs ####

pdf('~/Documents/Thesis/figures/STOPmutVSexp_boxplot.pdf')
par(font=2,font.axis=2,cex.axis=1.5,mar=c(5,5,3.5,0.6))
boxplot(MeanExpression~Type,ylim=c(6,16),las=1,frame=F,data=ptvexp,outline=T,xlab='',ylab='',boxcol=c('gray55','firebrick3'),boxlwd=3,col=c('gray70','#e88989'),medcol=c('gray40','firebrick4'),medlwd=3.5,staplecol=c('gray40','firebrick4'),staplelwd=3,border=c('gray40','firebrick4'),outcol=c('gray40','firebrick4'),outpch=20)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="#EBEBEB",border=NA)
grid(lty=1,col='white');par(new=T)
boxplot(MeanExpression~Type,ylim=c(6,16),las=1,frame=F,data=ptvexp,outline=T,xlab='',ylab='',boxcol=c('gray55','firebrick3'),boxlwd=3,col=c('gray70','#e88989'),medcol=c('gray40','firebrick4'),medlwd=3.5,staplecol=c('gray40','firebrick4'),staplelwd=3,border=c('gray40','firebrick4'),outcol=c('gray40','firebrick4'),outpch=20)
mtext(side=2,'Mean expression of genes (DESeq2 VST)',line=3,cex=1.5)
segments(x0=1,x1=2,y0=17,lwd=4,xpd=T);segments(x0=c(1,2),y0=c(17,17),y=c(16.5,16.5),lwd=4,xpd=T)
text(1.5,17.3,labels = "***",xpd=T,cex=1.5,font=2)
mtext(side=1,text='Gene mutation status',line=3,cex=1.5,font=2)
dev.off()



