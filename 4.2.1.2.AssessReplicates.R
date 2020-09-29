# Script to assess whether to combine sequencing replicates (R1/R2)
# To do the following:
# 1) Identify replicate samples
# 2) Load and process their raw gene counts
# 3) Assess the QC of these duplicates
# 4) Look at how well these samples cluster on plots
# 5) Normalise gene expression (including all other samples)
# 6) Assess the correlation of replicates compared to all other samples
# 7) Plot examples to illustrate

# Functions ####
# Plot the correlation of two samples normalised gene expression
plotcors <- function(gen1,gen2,colplot,sam1,sam2,let) {
  gen1 <- gen1+1;gen2 <- gen2+1
  spear <- signif(cor(log(gen1),log(gen2),method=cormethod),2)
  plot(log(gen1),log(gen2),pch=16,col=scales::alpha(colplot,0.3),cex=0.8,xlab=paste0(sam1,' (log(TPM+1))'),ylab=paste0(sam2,' (log(TPM+1))'),font=2,cex.lab=1.2)
  abline(a=0,b=1,col='gray',lty=2);legend('topleft',legend=paste0('rho: ',spear),bty='n',cex=1)
  mtext(side=3,at=1,text=let,font=2,line=1.5,cex=1.5)
  return(spear)
}

# Main code ####

# 0. Prepare environment
setwd("~/Documents/EPICC/Data/Expression/")
'%ni%' <- Negate('%in%');nucs <- c('A','C','G','T')
library(dplyr);library(data.table);library(stringr)

# 1. Load sample information and identify duplicates ####
# Load Rdata detailing the location and replication type of each sample
filelist <- readRDS('~/Documents/ThesisOther/ScriptsForThesis/RData/tmp.currentreadcountfileloc.rds')

# Generate a list of duplicated samples
duplist <- data.frame(Sample=filelist$Sample[which(duplicated(filelist$Sample))],Old=rep('',length(which(duplicated(filelist$Sample)))),OldFile=rep('',length(which(duplicated(filelist$Sample)))),New=rep('',length(which(duplicated(filelist$Sample)))),NewFile=rep('',length(which(duplicated(filelist$Sample)))))
for(i in c(1:nrow(duplist))) {
  curdup <- filelist[which(filelist$Sample==duplist[i,'Sample']),]
  if('S1' %in% curdup$ID) {
    duplist[i,c(2:5)] <- c('S1',curdup[which(curdup$ID=='S1'),'FileUse'],curdup[which(curdup$ID!='S1'),'ID'],curdup[which(curdup$ID!='S1'),'FileUse'])
  } else {
    duplist[i,c(2:5)] <- c('R1',curdup[which(curdup$ID=='R1'),'FileUse'],curdup[which(curdup$ID!='R1'),'ID'],curdup[which(curdup$ID!='R1'),'FileUse'])
  }
}
row.names(duplist) <- duplist$Sample

# 2. Load and process the raw gene counts of these duplicates ####
# Compile gene counts into one data frame
firstdup <- as.data.frame(fread(duplist[1,'OldFile']))
geneexp <- firstdup[c(5:(nrow(firstdup))),c('V1','V4')];colnames(geneexp) <- c('GeneID',paste0(row.names(duplist)[1],'_',duplist[1,'Old']))
seconddup <- as.data.frame(fread(duplist[1,'NewFile']))
geneexp[,paste0(row.names(duplist)[1],'_',duplist[1,'New'])] <- seconddup[c(5:(nrow(seconddup))),'V4']
for(i in c(2:nrow(duplist))) {
  print(row.names(duplist)[i])
  curdup <- duplist[i,]
  olddup <- as.data.frame(fread(curdup$OldFile))
  geneexp[,paste0(row.names(duplist)[i],'_',curdup$Old)] <- olddup[c(5:(nrow(olddup))),'V4']
  newdup <- as.data.frame(fread(curdup$NewFile))
  geneexp[,paste0(row.names(duplist)[i],'_',curdup$New)] <- newdup[c(5:(nrow(newdup))),'V4']
}

# Remove PAR_Y genes (never expressed) and make gene IDs stable
geneexp <- geneexp[-grep("PAR_Y",geneexp$GeneID),]
geneexp$GeneID <- gsub('(ENSG\\d+)\\.\\d+','\\1',geneexp$GeneID)
geneexp <- geneexp[order(geneexp$GeneID),]

# Save gene counts for duplicated samples
saveRDS(geneexp,file='~/Documents/ThesisOther/ScriptsForThesis/RData/reranGeneCounts.rds')

# 3. Generate QC information on duplicated samples ####
# Load gene information
ribo_genes <- read.table('~/Documents/EPICC/Data/Expression/GSEA/ribosomal_genes_hsa03010.txt')[,1]
geneexp <- readRDS('~/Documents/ThesisOther/ScriptsForThesis/RData/reranGeneCounts.rds')
samnames <- row.names(duplist);fullnames <- colnames(geneexp[,c(2:ncol(geneexp))])
geneinfo <- read.table("~/Documents/EPICC/Data/Expression/compiledGeneInfo.txt",header=T)
geneexpinfo <- merge(geneinfo,geneexp,by='GeneID')

# Assess the number of reads assigned to different categories of genes
mitos <- c();prots <- c();immos <- c();snos <- c();others <- c()
for(sam in fullnames) {
  sortCount <- geneexpinfo[order(geneexpinfo[,sam],decreasing=T),]
  mitos <- c(mitos,sum(sortCount[which(sortCount$Type=='Mt_rRNA'),sam]))
  prots <-c(prots,sum(sortCount[which(sortCount$Type=='protein_coding'),sam]))
  others <- c(others,sum(sortCount[which(sortCount$Type %ni% c('Mt_rRNA','protein_coding')),sam]))
}
QCtable <- geneexp[,fullnames]
QCtable <- as.data.frame(t(QCtable[c(1:2),]));colnames(QCtable) <- c('Assigned','ProteinCoding')
QCtable$Assigned <- colSums(as.matrix(geneexp[,fullnames]))/1e6
QCtable$UsableReads <- colSums(as.matrix(geneexpinfo[which(geneexpinfo$Chr %in% paste0('chr',c(1:22,'X','Y')) & geneexpinfo$Name %ni% ribo_genes & geneexpinfo$Type=='protein_coding'),fullnames]))/1e6
QCtable$ProteinCoding <- prots/1e6;QCtable$PercentProteinCoding <- QCtable$ProteinCoding/QCtable$Assigned*100
QCtable$MTrRNA <- mitos/1e6;QCtable$PercentMTrRNA <- QCtable$MTrRNA/QCtable$Assigned*100
QCtable$Other <- others/1e6;QCtable$PercentOther <- QCtable$Other/QCtable$Assigned*100

# Assess whether to use originals, repeats or (possibly) combine
olds <- c();news <- c()
for(r in c(1:nrow(duplist))) { olds <- c(olds,QCtable[paste0(row.names(duplist)[r],'_',duplist[r,'Old']),'UsableReads']);news <- c(news,QCtable[paste0(row.names(duplist)[r],'_',duplist[r,'New']),'UsableReads']) }
OldNew <- data.frame(Sample=samnames,OldID=duplist$Old,NewID=duplist$New,Old=olds,New=news,stringsAsFactors = F);OldNew$Total <- OldNew$Old+OldNew$New
OldNew$QC <- rep('Fail',nrow(OldNew))
for(i in c(1:nrow(OldNew))) {
  if(OldNew[i,'Old']>=5) {
    if(OldNew[i,'New']<5) { OldNew[i,'QC'] <- paste0('UseOld')
    } else { OldNew[i,'QC'] <- 'PossCombine' }
  } else if(OldNew[i,'Old']/OldNew[i,'New']>=0.5) {
    if(OldNew[i,'New']>=5) { OldNew[i,'QC'] <- 'PossCombine'
    } else if(OldNew[i,'Total']>=5) { OldNew[i,'QC'] <- 'Combine' }
  } else {
    if(OldNew[i,'New']>=5) { OldNew[i,'QC'] <- paste0('UseNew')
    } else if(OldNew[i,'Total']>=5) { OldNew[i,'QC'] <- 'Combine' }
  }
}

# 4. For samples which can benefit from being combined, assess their similarity by PCA and heatmap ####
# Filter for samples to combine
geneexp <- geneexp[which(geneexp$GeneID %in% geneinfo[which(geneinfo$Chr %in% paste0('chr',c(1:22,'X','Y')) & geneinfo$Name %ni% ribo_genes & geneinfo$Type=='protein_coding'),'GeneID']),];row.names(geneexp) <- c(1:nrow(geneexp))
combsam <- paste0(OldNew[which(OldNew$QC %in% c('Combine','PossCombine')),'Sample'],'_',OldNew[which(OldNew$QC %in% c('Combine','PossCombine')),'OldID'])
combsam <- c(combsam,paste0(OldNew[which(OldNew$QC %in% c('Combine','PossCombine')),'Sample'],'_',OldNew[which(OldNew$QC %in% c('Combine','PossCombine')),'NewID']))

# Reformat and run DESeq2
RERUN <- geneexp[,combsam];row.names(RERUN) <- geneexp$GeneID
RERUNdata <- as.data.frame(t(RERUN[c(1:2),]));colnames(RERUNdata) <- c('Patient','Type')
RERUNdata$Patient <- gsub('(C\\d+)_\\S\\d+_\\S\\d+_\\S+','\\1',colnames(RERUN))
RERUNdata$Type <- c(rep('Original',nrow(RERUNdata)/2),rep('Rerun',nrow(RERUNdata)/2))
dds <- DESeqDataSetFromMatrix(countData = RERUN,colData = RERUNdata,design = ~ Type)
dds <-DESeq(dds);vsd <- varianceStabilizingTransformation(dds, blind=T)
pca <- prcomp(t(assay(vsd)));sumpca <- summary(pca)$importance
mycols <- c(brewer.pal(12,'Paired')[c(1:10,12)],brewer.pal(6,'Dark2')[1:6])

# Plot PCA
pdf('~/Documents/ThesisOther/ScriptsForThesis/Plots/rerun_pca.pdf',width=10,height = 6)
par(font=2,font.axis=2,font.lab=2,mar=c(4.5,4.5,2,1))
plot(pca$x[,1],pca$x[,2],lwd=1,xlim=c(-150,750),cex.axis=1.2,cex.lab=1.3,ylim=c(-450,250),pch=c(rep(22,17),rep(24,17)),bg = rep(scales::alpha(mycols,0.3),2),cex=1.5,col=rep(mycols,2),xlab=paste0('PC1: ',signif(sumpca[2,1]*100,digits=4),'% variance'),ylab=paste0('PC2: ',signif(sumpca[2,2]*100,digits=4),'% variance'))
legend('topleft',legend=c('Original','Rerun'),pt.lwd = 2,pt.cex=2,pch=c(22,24),pt.bg=scales::alpha('black',0.3),border='black',bty='o',cex=1.2)
legend('bottomright',legend=gsub('(C\\S+)_\\S\\d+$','\\1',row.names(RERUNdata[1:17,])),fill=scales::alpha(mycols,0.75),border=NA,bty='o',cex=0.85)
dev.off()

# Plot sample-sample heatmap
pdf('~/Documents/ThesisOther/ScriptsForThesis/Plots/rerun_heatmap.pdf',width=10,height = 6)
mat <- as.matrix(dist(t(assay(vsd))))
patcols <- brewer.pal(12,'Paired')[c(1:length(unique(RERUNdata$Patient)))];names(patcols) <- unique(RERUNdata$Patient)
my_colour <- list(Patient=patcols,Type=c(Original=brewer.pal(5,'Pastel2')[1],Rerun=brewer.pal(5,'Pastel2')[4]))
pheatmap(mat,show_rownames=T,treeheight_col = 50,cluster_rows = T,treeheight_row = 50,annotation_col = RERUNdata[,c(2,1)],annotation_colors = my_colour,fontsize_row = 8,fontsize_col = 8)
dev.off()

# 5. Convert gene counts to TPM  ####
# Load pre-compiled gene length data and reformat
load(file="~/Documents/EPICC/Data/Expression/outputforNormalisation.RData");alllens <- as.data.frame(output)
alllens$GeneID <- row.names(alllens);row.names(alllens) <- c(1:nrow(alllens));alllens <- alllens[,c(3,1)]
alllens <- alllens[-grep("PAR_Y",alllens$GeneID),]
alllens$GeneID <- gsub('(ENSG\\d+)\\.\\d+','\\1',alllens$GeneID)
alllens <- alllens[order(alllens$GeneID),]

# Convert possible combining replicates to TPM
geneexplen <- merge(geneexp,alllens,by='GeneID')
normgeneexp <- geneexplen[,grep('C',colnames(geneexplen))];row.names(normgeneexp) <- geneexplen$GeneID
for(i in c(1:ncol(normgeneexp))) { normgeneexp[,i] <- normgeneexp[,i]/(geneexplen$Length/1000) }
for(i in c(1:ncol(normgeneexp))) {
  sfactor <- sum(normgeneexp[,i])/1000000
  normgeneexp[,i] <- normgeneexp[,i]/sfactor
}
tpmexp <- normgeneexp

# Have to run first 55 lines of 4.2.1.3.CompileReadCountsQC.R to get this file
# Then rerun those lines once this script has finished
# Load and normalise to TPM all gene count samples
allgenecounts <- as.data.frame(fread('~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_counts.allgenes.txt'))
allgenecounts <- allgenecounts[which(allgenecounts$GeneID %in% geneinfo[which(geneinfo$Chr %in% paste0('chr',c(1:22,'X','Y')) & geneinfo$Name %ni% ribo_genes & geneinfo$Type=='protein_coding'),'GeneID']),];row.names(allgenecounts) <- c(1:nrow(allgenecounts))
allgenelen <- merge(allgenecounts,alllens,by='GeneID')
normallgene <- allgenelen[,grep('C',colnames(allgenelen))];row.names(normallgene) <- allgenelen$GeneID
for(i in c(1:ncol(normallgene))) { normallgene[,i] <- normallgene[,i]/(allgenelen$Length/1000) }
for(i in c(1:ncol(normallgene))) {
  sfactor <- sum(normallgene[,i])/1000000
  normallgene[,i] <- normallgene[,i]/sfactor
}
tpmall <- normallgene

# Remove the 'duplicate' samples to be assessed from the all sample gene count data
notrep <- colnames(tpmall)[which(colnames(tpmall) %ni% OldNew$Sample)]
tpmall <- tpmall[,notrep]

# 6. Assess if correlation between replicates is higher than original to all other samples ####
# Setup matrix to store results in
decision <- c();cormethod='spearman';matres <- as.data.frame(matrix(0L,nrow=length(which(OldNew$QC %in% c('Combine','PossCombine'))),ncol=4))
row.names(matres) <- OldNew[which(OldNew$QC %in% c('Combine','PossCombine')),'Sample'];colnames(matres) <- c('RepCor','MeanOtherCor','Pvalue','Padj')

# For each replicate, get their correlation and compare to mean correlation of all other samples
for(i in c(1:nrow(OldNew))) {
  if(OldNew[i,'QC'] %in% c('UseNew','UseOld','Fail')) {
    decision <- c(decision,OldNew[i,'QC'])
  } else {
    currep <- OldNew[i,]
    print(currep$Sample)
    gen1 <- tpmexp[,paste0(currep$Sample,'_',currep$OldID)]+1
    gen2 <- tpmexp[,paste0(currep$Sample,'_',currep$NewID)]+1
    spear <- signif(cor(log(gen1),log(gen2),method=cormethod),2)
    
    distcors <- c()
    for(j in c(2:ncol(tpmall))) {
      gen2 <- tpmall[,j]
      distcors <- c(distcors,cor(log(gen1),log(gen2),method=cormethod))
    }
    res <- wilcox.test(distcors,mu=spear,alternative = 'less')
    # If correlation significantly higher, replicate counts can be combined
    if(p.adjust(res$p.value,method = 'fdr',n=17)<0.01) {
      decision <- c(decision,'Combine')
    } else {
      if(OldNew[i,'QC']=='PossCombine') {
        decision <- c(decision,'UseNew')
      } else if(OldNew[i,'QC']=='Combine') {
        decision <- c(decision,'Fail')
      }
    }
    matres[currep$Sample,] <- c(spear,mean(distcors),res$p.value,p.adjust(res$p.value,method = 'fdr',n=17))
  }
}
OldNew$Decision <- decision

# Save RData for use in future scripts
saveRDS(matres,file='~/Documents/ThesisOther/ScriptsForThesis/RData/combining_correlation_results.rds')
saveRDS(OldNew,file='~/Documents/ThesisOther/ScriptsForThesis/RData/reranQCDecision.rds')

# 7. Plot examples of replicates that pass and fail the combination test ####
pdf('~/Documents/Thesis/figures/ReplicateCorrelations.pdf')
layout(matrix(c(1:4),nrow=2,ncol=2,byrow = T))
par(font=2,font.lab=2,font.axis=2,mar=c(4.5,4.5,3.5,1))

# Correlation plot of good (PASS/COMBINE) sample
sam <- 'C552_B1_G9';currep <- OldNew[which(OldNew$Sample==sam),]
gen1 <- tpmexp[,paste0(currep$Sample,'_',currep$OldID)]
gen2 <- tpmexp[,paste0(currep$Sample,'_',currep$NewID)]
spear <- plotcors(gen1,gen2,'forestgreen',paste0(currep$Sample,'_',currep$OldID),paste0(currep$Sample,'_',currep$NewID),'(a)')
# Histogram comparing to all other rho values for good sample
distcors <- c()
for(j in c(2:ncol(tpmall))) {
  if(colnames(tpmall)[j]!=currep$Sample) {
    gen2 <- tpmall[,j]
    gen1 <- gen1+1;gen2 <- gen2+1
    distcors <- c(distcors,cor(log(gen1),log(gen2),method=cormethod))
  }
}
hist(distcors,col='dimgray',border=NA,breaks=50,xlim=c(0,1),main='',xlab='Rho values');abline(v=spear,col='forestgreen',lwd=3)
mtext(side=3,text='one-sided wilcoxon test: FDR<0.001',line=0.5,cex=0.8)
mtext(side=3,at=0,text='(b)',font=2,line=1.5,cex=1.5)

# Correlation plot of bad (PASS/COMBINE) sample
sam <- 'C544_A1_G9';currep <- OldNew[which(OldNew$Sample==sam),]
gen1 <- tpmexp[,paste0(currep$Sample,'_',currep$OldID)]
gen2 <- tpmexp[,paste0(currep$Sample,'_',currep$NewID)]
spear <- plotcors(gen1,gen2,'firebrick3',paste0(currep$Sample,'_',currep$OldID),paste0(currep$Sample,'_',currep$NewID),'(c)')

# Histogram comparing to all other rho values for bad sample
distcors <- c()
for(j in c(2:ncol(tpmall))) {
  gen2 <- tpmall[,j]
  gen1 <- gen1+1;gen2 <- gen2+1
  distcors <- c(distcors,cor(log(gen1),log(gen2),method=cormethod))
}
hist(distcors,col='dimgray',border=NA,breaks=50,xlim=c(0,1),main='',xlab='Rho values');abline(v=spear,col='firebrick3',lwd=3)
mtext(side=3,text='one-sided wilcoxon test: FDR=1',line=0.5,cex=0.8)
mtext(side=3,at=0,text='(d)',font=2,line=1.5,cex=1.5)
dev.off()
