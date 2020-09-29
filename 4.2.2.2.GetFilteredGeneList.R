# Script to generate expression filtered gene list
# To do the following:
# 1) Load ensembl TPM data
# 2) Filter for expressed genes and output result
# 3) Load symbol TPM data
# 4) Filter for expressed genes and output result

# Functions ####

# Main code ####

# 0. Prepare environment
setwd('~/Documents/EPICC/Data/Expression/');'%ni%' <- Negate('%in%')
library(DESeq2);library(data.table)

# 1. Load ensembl TPM data and filter for PASS samples ####
tpm <- as.data.frame(fread('ProcessedCounts/All_EPICC_tpm.txt'))
passsam <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]
tpm <- tpm[,c('GeneID',passsam)]

# 2. Find which genes are >=1TPM in >=20% of PASS samples ####
perc20 <- ceiling(length(passsam)*.2)
atl1 <- apply(tpm[,passsam],1,function(x) sum(x>=1))
filens <- tpm$GeneID[which(atl1>=perc20)]
write.table(filens,file='~/Documents/ThesisOther/ScriptsForThesis/Tables/filteredgenes.20perAL1TPM.ensembl.txt',row.names = F,quote = F,col.names = F)

# 3. Load symbol TPM data and filter for PASS samples ####
tpm <- as.data.frame(fread('ProcessedCounts/All_EPICC_symbol_tpm.txt'))
passsam <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]
tpm <- tpm[,c('GeneID',passsam)]

# 4. Find which genes are >=1TPM in >=20% of PASS samples ####
perc20 <- ceiling(length(passsam)*.2)
atl1 <- apply(tpm[,passsam],1,function(x) sum(x>=1))
filsym <- tpm$GeneID[which(atl1>=perc20)]
write.table(filsym,file='~/Documents/ThesisOther/ScriptsForThesis/Tables/filteredgenes.20perAL1TPM.symbol.txt',row.names = F,quote = F,col.names = F)

