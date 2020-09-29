# Script to make phylogenetic like clustering dendrograms
# To do the following:
# 1) Load normalised gene expression and epithelial expression information
# 2) Plot dendrograms based on all expressed genes
# 3) Plot dendrograms based on expressed epithelial genes
# 4) Plot dendrograms based on expressed stromal genes
# 5) Save sample-sample distances for use in mixing analysis

# Functions ####
# Make a dendrogram based on gene expression - return the sample-sample distances
make_dendro <- function(patsam,pat,filgeneexp,my_colour) {
  EPICC <- filgeneexp[,patsam];row.names(EPICC) <- filgeneexp$GeneID
  colnames(EPICC) <- gsub('C\\d+_(\\S\\d+_\\S\\d+)','\\1',colnames(EPICC))
  EPICCdata <- as.data.frame(t(EPICC[c(1:2),]));colnames(EPICCdata) <- c('Region','Type')
  EPICCdata$Region <- as.factor(gsub("(\\S)\\d+_.+$","\\1",row.names(EPICCdata)))
  EPICCdata$Type <- as.factor(sapply(row.names(EPICCdata),function(x) if(gsub("\\S\\d+_(\\S).+$","\\1",x)=='B') {x<-'Bulk'} else {x<-'Gland'}))
  
  cexplot <- ifelse(length(patsam)>=30,1,1.35)
  dd <- dist(t(EPICC),method='euclidean')
  hc <- hclust(dd, method = "ward.D2")
  plot(as.phylo(hc),cex = cexplot, label.offset = 2,edge.width=2,tip.color=my_colour$Region[EPICCdata$Region],font=2)
  mtext(text=pat,side=3,at=0,cex=1.5)
  add.scale.bar(1,0.75,lwd=2)
  return(dd)
}

# Main code ####

# 0. Prepare environment
setwd('~/Documents/ThesisOther/ScriptsForThesis/');library(dplyr);library(dendextend)
library(pheatmap);library(pvclust);library(fpc);library(phangorn);library(DESeq2)
my_colour <- list(Region=c(A='#E31A1C',B='#377DB8',C='#4DAE49',D='#904A9A'),Type=c(Bulk='gray30',Gland="#e89829"))

# 1. Load expression data and other information ####
# Load in gene expression 
vsd <- readRDS('~/Documents/EPICC/Data/Expression/ProcessedCounts/filgenes.vsd.ensembl.rds')
geneexp <- as.data.frame(assay(vsd));geneexp$GeneID <- row.names(assay(vsd))
geneinfo <- read.table('~/Documents/EPICC/Data/Expression/compiledGeneInfo.txt',header=T)

# Load in sample data. filter for patients that have at least 5 samples
samples <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]
samples <- samples[which(!grepl('E1',samples))]
geneexp <- geneexp[,c('GeneID',samples)]
allpat <- gsub('(C\\d+)\\S+','\\1',samples);patients <- unique(allpat)
numpat <- table(allpat)
multipat <- names(numpat[which(numpat>=5)])
newsam <- samples[which(allpat %in% multipat)]
newallpat <- gsub('(C\\d+)\\S+','\\1',newsam)
geneexp <- geneexp[,c('GeneID',newsam)]

# Load in supplementary Table 11a from Isella (2015)
isellasupp <- as.data.frame(readxl::read_xlsx('~/Documents/EPICC/Data/Expression/Classification/CRIS/IsellaSupp11.xlsx')[,c(1:3)])
colnames(isellasupp) <- c('Mouse','GeneSymbol','PercMurine')

# 2. Plot all (expressed) gene expression clustering ####
filgeneexp <- geneexp;distlist <- list()

pdf('~/Documents/ThesisOther/ScriptsForThesis/Plots/ExpPhylogenies.pdf',width=18,height=10)
layout(matrix(c(1:18),nrow=3,ncol=6,byrow = T))
par(mar=c(2.1,2.1,2.1,0.1),xpd=T)
for(p in c(1:length(multipat))) {
  pat <- multipat[p]
  patsam <- newsam[which(newallpat==pat)]
  patsam <- patsam[which(!grepl('_E',patsam))]
  distlist[[pat]] <- make_dendro(patsam,pat,filgeneexp,my_colour)
}
dev.off()

# 3. Plot epithelial gene expression clustering ####
# Get expressed genes that are of non-stromal origin
epithelial <- isellasupp[which(isellasupp$PercMurine<=0.25),'GeneSymbol']
ens_epithelial <- geneinfo[which(geneinfo$Name %in% epithelial),'GeneID']
filgeneexp <- geneexp[which(geneexp$GeneID %in% unique(ens_epithelial)),]
distlist_epi <- list()

pdf('~/Documents/ThesisOther/ScriptsForThesis/Plots/ExpPhylogenies_Epithelial.pdf',width=18,height=10)
layout(matrix(c(1:18),nrow=3,ncol=6,byrow = T))
par(mar=c(2.1,2.1,2.1,0.1),xpd=T)
for(p in c(1:length(multipat))) {
  pat <- multipat[p]
  patsam <- newsam[which(newallpat==pat)]
  patsam <- patsam[which(!grepl('_E',patsam))]
  distlist_epi[[pat]] <- make_dendro(patsam,pat,filgeneexp,my_colour)
}
dev.off()

# 4. Plot stromal gene expression clustering ####
# Get expressed genes that are of non-stromal origin
stromal <- isellasupp[which(isellasupp$PercMurine>=0.5),'GeneSymbol']
ens_stromal <- geneinfo[which(geneinfo$Name %in% stromal),'GeneID']
filgeneexp <- geneexp[which(geneexp$GeneID %in% unique(ens_stromal)),]
distlist_str <- list()

pdf('~/Documents/ThesisOther/ScriptsForThesis/Plots/ExpPhylogenies_Stromal.pdf',width=18,height=10)
layout(matrix(c(1:18),nrow=3,ncol=6,byrow = T))
par(mar=c(2.1,2.1,2.1,0.1),xpd=T)
for(p in c(1:length(multipat))) {
  pat <- multipat[p]
  patsam <- newsam[which(newallpat==pat)]
  patsam <- patsam[which(!grepl('_E',patsam))]
  distlist_str[[pat]] <- make_dendro(patsam,pat,filgeneexp,my_colour)
}
dev.off()

# 5. Save list of sample-sample distances for each analysis ####
save(distlist,distlist_epi,distlist_str,file='~/Documents/ThesisOther/ScriptsForThesis/RData/expression_clustering_distances.Rdata')

