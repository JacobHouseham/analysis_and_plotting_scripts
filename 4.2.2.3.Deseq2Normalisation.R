# Script to convert raw gene counts into DESeq2 variance stabilising
# normalised counts
# To do the following:
# 1) Load ensembl raw counts for expressed genes and PASS samples
# 2) Run DESeq2 and output normalised counts for ensembl gene IDs
# 3) Load symbol raw counts for expressed genes and PASS samples
# 4) Run DESeq2 and output normalised counts for symbol gene IDs

# Script run on HPC as it takes a long time

# Main code ####

# 0. Prepare environment
setwd('~/Documents/EPICC/Data/Expression/');'%ni%' <- Negate('%in%')
library(DESeq2);library(data.table)

# 1) Load ensembl raw gene counts and filter for expressed genes ####
EPICC <- as.data.frame(fread('ProcessedCounts/All_EPICC_counts.txt'))
passsam <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]
filteredGenes <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/filteredgenes.20perAL1TPM.ensembl.txt')[,1]
EPICC <- EPICC[which(EPICC$GeneID %in% filteredGenes),c('GeneID',passsam)]

# 2) Run DEseq2 on ensembl raw counts ####
# Convert ready for DESeq2
row.names(EPICC) <- EPICC[,1];EPICC <- EPICC[,c(2:ncol(EPICC))]
EPICCdata <- as.data.frame(t(EPICC[c(1:2),]));colnames(EPICCdata) <- c('Patient','Type')
EPICCdata$Patient <- gsub('(C\\d+)\\S+','\\1',row.names(EPICCdata))
regs <- gsub('C\\d+_(\\S)\\d+_\\S+','\\1',row.names(EPICCdata))
EPICCdata$Type <- ifelse(regs=='E','Normal','Tumour')

# Run Deseq2 to get dds and normalised counts
dds <- DESeqDataSetFromMatrix(countData = EPICC,colData = EPICCdata,design = ~ Patient + Type)
dds <-DESeq(dds);vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

# Output vsd assay as normalised read counts
saveRDS(vsd,file='~/Documents/ThesisOther/ScriptsForThesis/Rdata/filgenes.vsd.ensembl.rds')

# 1) Load symbol raw gene counts and filter for expressed genes ####
EPICC <- as.data.frame(fread('ProcessedCounts/All_EPICC_symbol_counts.txt'))
passsam <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]
filteredGenes <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/filteredgenes.20perAL1TPM.symbol.txt')[,1]
EPICC <- EPICC[which(EPICC$GeneID %in% filteredGenes),c('GeneID',passsam)]

# 2) Run DEseq2 on symbol raw counts ####
# Convert ready for DESeq2
row.names(EPICC) <- EPICC[,1];EPICC <- EPICC[,c(2:ncol(EPICC))]
EPICCdata <- as.data.frame(t(EPICC[c(1:2),]));colnames(EPICCdata) <- c('Patient','Type')
EPICCdata$Patient <- gsub('(C\\d+)\\S+','\\1',row.names(EPICCdata))
regs <- gsub('C\\d+_(\\S)\\d+_\\S+','\\1',row.names(EPICCdata))
EPICCdata$Type <- ifelse(regs=='E','Normal','Tumour')

# Run Deseq2 to get dds and normalise
dds_sym <- DESeqDataSetFromMatrix(countData = EPICC,colData = EPICCdata,design = ~ Patient + Region)
dds_sym <-DESeq(dds_sym);vsd_sym <- varianceStabilizingTransformation(dds_sym, blind=TRUE)

# Output vsd assay as normalised read counts
saveRDS(vsd_sym,file='~/Documents/ThesisOther/ScriptsForThesis/Rdata/filgenes.vsd.symbol.rds')


