# Script to look for recurrent aberrant splicing patterns
# To do the following:
# 1) Load expression data and aberrant splicing data
# 2) Get matrix of aberrant splicing sample by gene
# 3) Get matrix of aberrant splicing tumour by gene
# 4) Perform KEGG enrichment of recurrent genes

# Main code ####

# 0. Prepare environment
library(data.table);setwd("~/Documents/Splicing");load("~/Documents/EPICC/Data/Splicing/refinedGFF.RData")
"%ni%" <- Negate("%in%")
geneinfo <- read.table("~/Documents/EPICC/Data/Expression/compiledGeneInfo.txt",header=T)
genehg38 <- merge(geneinfo,EPICC,by='GeneID')
tsgs <- read.table('~/Documents/EPICC/Data/Splicing/TSGs300.txt')[,1]
ensg_tsgs <- geneinfo[which(geneinfo$Name %in% tsgs),'GeneID']

# 1. Load expression and aberrant splicing data ####
EPICC <- as.data.frame(fread(paste0("~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_tpm.txt"),header=T))
rnasam <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]

load(file='~/Documents/ThesisOther/ScriptsForThesis/RData/aberrant_splicing_data.Rdata')

allsam <- sort(unique(c(samlist$Samples,samlist$MatchedNorm)))
patients <- sort(unique(gsub('(C\\d+)_.+','\\1',samlist$MatchedNorm)))
EPICC <- EPICC[,c('GeneID',allsam)]

# 2. Find recurrent aberrant splicing events ####
# A '2' signifies high expression and aberrant splicing in that sample
recsplice <- matrix(0,nrow=length(allevents),ncol=length(usesam));colnames(recsplice) <- usesam
row.names(recsplice) <- allevents
expsplice <- recsplice
for(i in c(1:nrow(samlist))) {
  sam <- samlist$Samples[i];normal <- samlist$MatchedNorm[i]
  cursig <- spliceList[[sam]]
  expgenes <- EPICC$GeneID[EPICC[,sam]>=10 & EPICC[,normal]>=10]
  
  recsplice[row.names(recsplice) %ni% cursig$GeneID & row.names(recsplice) %in% expgenes,sam] <- 1
  recsplice[row.names(recsplice) %in% cursig$GeneID & row.names(recsplice) %in% expgenes,sam] <- 2
}
head(sort(rowSums(recsplice==2),decreasing = T))

# 3. Find aberrant splicing events that occur in at least half of the samples of each tumour ####
# Here, a '2' signifies that at least half the samples of that tumour had expression and aberrant splicing of that gene
allpats <- gsub('(C\\d+)_\\S+','\\1',samlist$Samples);patients <- unique(allpats)
recpat <- recsplice[,c(1:length(patients))];colnames(recpat) <- patients
for(pat in patients) {
  patsam <- samlist$Samples[grep(pat,samlist$Samples)]
  patsplice <- recsplice[,patsam]
  recpat[,pat] <- ifelse(rowSums(patsplice==2)>=(length(patsam)/2),2,ifelse(rowSums(patsplice==1)>=(length(patsam)/2),1,0))
}
length(which(rowSums(recpat)==12))
namesplice <- recsplice;row.names(namesplice) <- geneinfo[match(row.names(recsplice),geneinfo$GeneID),'Name']
# So with 6 tumours, recurrently aberrantly spliced genes will had 2 in every row (6*2=12)
symrec <- row.names(namesplice)[which(rowSums(recpat)==12)]

# 4. Perform KEGG enrichment analysis on recurrently aberrantly spliced genes ####
library(clusterProfiler)
keggpaths <- read.gmt('~/Downloads/c2.cp.kegg.v7.1.symbols.gmt')
keggres <- enricher(symrec, TERM2GENE=keggpaths,pAdjustMethod = 'fdr',pvalueCutoff = 0.01,qvalueCutoff = 0.01)
View(keggres@result[which(keggres@result$qvalue<0.01),])
keggres@result$Description <- gsub('KEGG_(\\S+)','\\1',keggres@result$Description)

pdf('~/Documents/Thesis/figures/Splicing_Recurrent_KEGG.pdf',width=10)
dotplot(keggres,font.size=10)
dev.off()

