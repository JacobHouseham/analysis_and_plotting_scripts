# Script to assess how closely matched `top-up` RNA-seq samples are to original (data from C518)
# To do the following:
# 1) Load pre-compiled gene counts for original and topup samples
# 2) Filter for protein coding and non-ribosomal genes
# 3) Run DESeq2 to get normalised counts
# 4) Plot PCA and sample-sample heatmap

# Main code ####

# 0. Prepare environment ####
setwd("~/Documents/EPICC/QC/TOPUP_QC/")
'%ni%' <- Negate('%in%');library(DESeq2);library(RColorBrewer)
library(dplyr);library(pheatmap)

# 1. Load raw counts and other information ####
counts <- read.table('TOPUPcounts.txt',header=T)

# Remove suffixes so that GeneIDs represent stable ensembl gene IDs
counts$GeneID <- gsub('(ENSG\\d+)\\.\\d+','\\1',counts$GeneID)

# Load information on gene location and biotype
geneinfo <- read.table('~/Documents/EPICC/Data/Expression/compiledGeneInfo.txt',header=T)

# Load KEGG ribosomal genes
ribo_genes <- read.table('~/Documents/EPICC/Data/Expression/GSEA/ribosomal_genes_hsa03010.txt')[,1]

# 2. Filter gene counts ####
# Filter for canonical chromosomes and protein coding genes
counts <- counts[which(counts$GeneID %in% geneinfo[which(geneinfo$Chr %in% paste0('chr',c(1:22,'X','Y'))),'GeneID'] & geneinfo$Type=='protein_coding'),];row.names(counts) <- c(1:nrow(counts))

# Filter out ribosomal protein coding genes
counts <- counts[which(counts$GeneID %in% geneinfo[which(geneinfo$Name %ni% ribo_genes),'GeneID']),];row.names(counts) <- c(1:nrow(counts))

# Filter out genes with 0 read counts in any sample (for speed)
nonzero <- which(rowSums(counts[,c(2:ncol(counts))])!=0)
counts <- counts[nonzero,]
row.names(counts) <- counts$GeneID;counts <- counts[,c(2:ncol(counts))]

TOPUP <- counts[,c(grep('RNA',colnames(counts)),grep('TOP',colnames(counts)))]
samples <- unique(gsub('(C\\d+_\\S\\d+_\\S\\d+)_\\S+','\\1',colnames(counts)))

# 3. Prepare and then run DESeq2 ####
# Make metadata dataframe for use in DEseq2
TOPUPdata <- as.data.frame(t(TOPUP[c(1:2),]));colnames(TOPUPdata) <- c('Region','Type')
TOPUPdata$Region <- gsub('C\\d+_(\\S)\\d+_\\S\\d+_\\S+','\\1',colnames(TOPUP))
TOPUPdata$Type <- c(rep('Original',12),rep('TopUp',12))

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = TOPUP,colData = TOPUPdata,design = ~ Type)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind=F)

# 4. Plot the normalised counts ####
# Make a sample-sample heatmap
pdf('~/Documents/Thesis/figures/topup_heatmap.pdf',width=10,height = 6)
mat <- as.matrix(dist(t(assay(vsd))))
my_colour <- list(Region=c(A='#E31A1C',B='#377DB8',C='#4DAE49',D='#904A9A'),Type=c(Original=brewer.pal(5,'Pastel2')[1],TopUp=brewer.pal(5,'Pastel2')[4]))
pheatmap(mat,show_rownames=T,treeheight_col = 50,cluster_rows = T,treeheight_row = 50,annotation_col = TOPUPdata[,c(2,1)],annotation_colors = my_colour,fontsize_row = 8,fontsize_col = 8)
dev.off()

# Make a PCA plot
pca <- prcomp(t(assay(vsd)));sumpca <- summary(pca)$importance
mycols <- c(brewer.pal(12,'Paired')[c(1:10,12)],brewer.pal(3,'Dark2')[1])

pdf('~/Documents/Thesis/figures/topup_pca.pdf',width=10,height = 6)
par(font=2,font.axis=2,font.lab=2,mar=c(4.5,4.5,2,1))
plot(pca$x[,1],pca$x[,2],lwd=2,xlim=c(-225,250),cex.axis=1.2,cex.lab=1.3,ylim=c(-80,200),pch=c(rep(22,12),rep(24,12)),bg = rep(scales::alpha(mycols,0.3),2),cex=3,col=rep(mycols,2),xlab=paste0('PC1: ',signif(sumpca[2,1]*100,digits=4),'% variance'),ylab=paste0('PC2: ',signif(sumpca[2,2]*100,digits=4),'% variance'))
legend('topleft',legend=c('Original','Topup'),pt.lwd = 2,pt.cex=2,pch=c(22,24),pt.bg=scales::alpha('black',0.3),border='black',bty='o',cex=1.3)
legend('topright',legend=samples,fill=scales::alpha(mycols,0.75),border=NA,bty='o',cex=0.9)
dev.off()







