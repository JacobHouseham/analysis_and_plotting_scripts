# Script to assess the relationship of splicing evenness to mutations in splicing factors
# To do the following:
# 1) Load splicing evenness and splicing factor data
# 2) Load mutation data and filter for splicing factor gene mutations
# 3) Assess the association of splicing factor mutations and evenness
# 4) Determine the recurrence within and across tumours of splicing factor mutations

# Main code ####

# 0. Prepare environment
setwd('~/Documents/EPICC/Data/Splicing/')
EnsemblConseq <- read.csv('~/Documents/EPICC/Data/Mutations/VCFtoTSV/ensembl_consequence_ordered.csv')
orderConseq <- EnsemblConseq$Term
nsvs <- orderConseq[c(1:12)]

# 1. Load splicing evenness data and splicing factor information ####
# Load splicing factors and find ensg ids
splicing <- read.table('splicing_factors_seiler_cell_2018.txt')[,1]
geneinfo <- as.data.frame(fread('~/Documents/EPICC/Data/Expression/compiledGeneInfo.txt'))
splice_ensg <- sort(unique(geneinfo[which(geneinfo$Name %in% splicing),'GeneID']))

# Load results of evenness analysis
misodf <- readRDS('~/Documents/ThesisOther/ScriptsForThesis/RData/evenness_results.rds')

# 2. Load non-synonymous mutation data and filter for those in splicing factor genes ####
splicemat <- readRDS('~/Documents/EPICC/Data/QTLs/allvariants/NS_mutburdengene_allDNAsamples.rds')
samples <- misodf$Sample[which(misodf$Sample %in% colnames(splicemat))]
matchsplice <- splice_ensg[which(splice_ensg %in% row.names(splicemat))]
splicemat <- splicemat[matchsplice,samples]
splicemat <- splicemat[which(rowSums(splicemat)>0),]
splicemat <- merge(splicemat,geneinfo[,c('GeneID','Name')],by.x=0,by.y='GeneID')
row.names(splicemat) <- splicemat$Name;splicemat <- as.matrix(splicemat[,samples])

misodf <- misodf[which(misodf$Sample %in% samples),]

# 3. Compare evenness when samples have mutations in splicing factors ####
# Overall for any splicing factor
resSF <- data.frame(Sample=samples,Mutated=ifelse(colSums(splicemat)==0,'WT','Mut'))
resSF <- merge(resSF,misodf[,c('Sample','MedianEven','MedianEven')],by='Sample');resSF$Mutated <- factor(resSF$Mutated,levels=c('WT','Mut'))
wilcox.test(MedianEven~Mutated,data=resSF,alternative='less');boxplot(MedianEven~Mutated,data=resSF)

# For one splicing factor at a time
pvals <- c()
for(i in c(1:nrow(splicemat))) {
  resSF <- data.frame(Sample=samples,Mutated=ifelse(splicemat[i,]==0,'WT','Mut'))
  resSF <- merge(resSF,misodf[,c('Sample','MedianEven','Skew')],by='Sample');resSF$Mutated <- factor(resSF$Mutated,levels=c('WT','Mut'))
  res <- wilcox.test(MedianEven~Mutated,data=resSF,alternative='less')
  boxplot(MedianEven~Mutated,main=row.names(splicemat)[i],xlab=res$p.value,data=resSF)
  pvals <- c(pvals,res$p.value)
}
names(pvals) <- row.names(splicemat)
names(which(p.adjust(pvals,method='fdr')<0.01))

library(pheatmap)
pheatmap(splicemat,fontsize = 5,cluster_cols = F)

# 4. Look at recurrence of splicing factor mutations and their effect on eveness ####
srecpat <- matrix(0L,nrow=nrow(splicemat),ncol=length(unique(misodf$Patient)))
row.names(srecpat) <- row.names(splicemat);colnames(srecpat) <- unique(misodf$Patient)
for(i in c(1:ncol(srecpat))) {
  pat <- colnames(srecpat)[i]
  patsplice <- splicemat[,grep(pat,colnames(splicemat))]
  if(is.matrix(patsplice)) {
    srecpat[rowSums(patsplice)>0,i] <- 1
  } else {
    srecpat[patsplice>0,i] <- 1
  }
}
srecpat <- srecpat[,which(colSums(srecpat)!=0)]
pheatmap(srecpat,fontsize = 5)

