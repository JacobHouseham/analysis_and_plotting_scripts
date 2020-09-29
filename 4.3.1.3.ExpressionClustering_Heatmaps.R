# Script to produce a heatmap clustering by gene expression for all samples
# To do the following:
# 1) Load data
# 2) Select tumours with 5 or more samples
# 2) Plot heatmap

# Main code ####

# 0. Prepare environment
setwd('~/Documents/ThesisOther/ScriptsForThesis/');library(dplyr);library(dendextend)
library(pheatmap);library(pvclust);library(fpc);library(RColorBrewer);library(made4)

# 1. Load expression data and filter for epithelial expression ####
# Load in supplementary Table 11a from Isella (2015)
isellasupp <- as.data.frame(readxl::read_xlsx('~/Documents/EPICC/Data/Expression/Classification/CRIS/IsellaSupp11.xlsx')[,c(1:3)])
colnames(isellasupp) <- c('Mouse','GeneSymbol','PercMurine')
nonstromal <- isellasupp[which(isellasupp$PercMurine<=0.25),'GeneSymbol']

# Load in gene expression
vsd <- readRDS('~/Documents/EPICC/Data/Expression/ProcessedCounts/filgenes.vsd.ensembl.rds')
geneexp <- as.data.frame(assay(vsd));geneexp$GeneID <- row.names(assay(vsd))
geneinfo <- read.table('~/Documents/EPICC/Data/Expression/compiledGeneInfo.txt',header=T)

# Get filtered (so expressed genes) that are of non-stromal origin
ens_nonstromal <- geneinfo[which(geneinfo$Name %in% nonstromal),'GeneID']
geneexp <- geneexp[which(geneexp$GeneID %in% unique(ens_nonstromal)),]

# 2. Load sample data and filter for patients that have at least 5 samples ####
samples <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]
samples <- samples[which(!grepl('E1',samples))]
geneexp <- geneexp[,c('GeneID',samples)]
allpat <- gsub('(C\\d+)\\S+','\\1',samples);patients <- unique(allpat)
numpat <- table(allpat)
multipat <- names(numpat[which(numpat>=5)])
newsam <- samples[which(allpat %in% multipat)]
newallpat <- gsub('(C\\d+)\\S+','\\1',newsam)
geneexp <- geneexp[,c('GeneID',newsam)]

# 3. Plot all samples together, is ITH greater than inter-patient heterogeneity? ####
patsam <- newsam
EPICC <- geneexp[,patsam];row.names(EPICC) <- geneexp$GeneID
EPICCdata <- as.data.frame(t(EPICC[c(1:3),]));colnames(EPICCdata) <- c('Patient','Region','Type')
EPICCdata$Patient <- newallpat[which(newsam %in% patsam)]
EPICCdata$Region <- as.factor(gsub("C\\d\\d\\d_(\\S)\\d+_.+$","\\1",row.names(EPICCdata)))
EPICCdata$Type <- ifelse(grepl('C\\d+_\\S\\d+_G',patsam),'Gland','Bulk')

patcols <- getcol(c(1:18), palette="colours2")
names(patcols) <- unique(EPICCdata$Patient)
my_colour <- list(Patient=patcols,Region=c(A='#E31A1C',B='#377DB8',C='#4DAE49',D='#904A9A'),Type=c(Bulk='gray30',Gland="#e89829"))
pdf('~/Documents/Thesis/figures/EpithelialExp_AllSamples_Clustering.pdf',width = 14,height = 14)
pheatmap(EPICC,show_rownames=F,treeheight_col = 100,cluster_rows = F,treeheight_row = 0,fontsize_col = 3,annotation_col = EPICCdata[,c(3,2,1)],annotation_colors = my_colour,clustering_method = "ward.D2")
dev.off()

