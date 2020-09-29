# Script for loading in all count data and assessing QC 
# To do the following:
# 1) Load all gene count data for all samples
# 2) Convert to TPM and log(TPM+1)
# 3) Convert Ensembl gene IDs to symbol
# 4) Convert gene symbol counts to TPM

# Functions ####
# Convert a dataframe of TPM into log(tpm+1)
TPMTologTpm <- function(tpm) {
  for(i in c(1:ncol(tpm))) { tpm[,i] <- log(tpm[,i]+1) }
  return(tpm)
}

# Main code ####

# 0. Prepare environment
setwd("~/Documents/EPICC/Data/Expression");library(data.table)
library(dplyr);'%ni%' <- Negate('%in%')

# 1. Load data and reformat ready for normalisation ####
# Load gene count matrix from 4.2.1.3
EPICC <- as.data.frame(fread('ProcessedCounts/All_EPICC_counts.txt'))

# Load pre-compiled gene length data
load(file="outputforNormalisation.RData");alllens <- as.data.frame(output)
alllens$GeneID <- row.names(alllens);row.names(alllens) <- c(1:nrow(alllens));alllens <- alllens[,c(3,1)]
alllens <- alllens[-grep("PAR_Y",alllens$GeneID),]
alllens$GeneID <- gsub('(ENSG\\d+)\\.\\d+','\\1',alllens$GeneID)
alllens <- alllens[order(alllens$GeneID),]

# Merge data together
EPICC <- merge(EPICC,alllens,by='GeneID')

# 2. Convert raw gene counts to TPM and log(TPM+1) ####
# a. Normalise for gene length: gene counts / gene length (in kb)
normEPICC <- EPICC[,grep('C',colnames(EPICC))];row.names(normEPICC) <- EPICC$GeneID
for(i in c(1:ncol(normEPICC))) { normEPICC[,i] <- normEPICC[,i]/(EPICC$Length/1000) }

# b. Normalise for sequencing depth: sum the normalised gene counts and divide by a million
# then divide each normalised gene count by that samples scaling factor
for(i in c(1:ncol(normEPICC))) {
  sfactor <- sum(normEPICC[,i])/1000000
  normEPICC[,i] <- normEPICC[,i]/sfactor
}
epiccTPM <- normEPICC

# TPM to log(TPM+1)
epicclogTPM <- TPMTologTpm(epiccTPM)

# Output TPM and log(TPM+1) files
# TPM
epiccTPM$GeneID <- row.names(epiccTPM);epiccTPM <- epiccTPM[,c(ncol(epiccTPM),1:(ncol(epiccTPM)-1))]
write.table(epiccTPM,"~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_tpm.txt",sep='\t',quote=F,row.names = F)
# logTPM
epicclogTPM$GeneID <- row.names(epicclogTPM);epicclogTPM <- epicclogTPM[,c(ncol(epicclogTPM),1:(ncol(epicclogTPM)-1))]
write.table(epicclogTPM,"~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_logtpm.txt",sep='\t',quote=F,row.names = F)

# 3. Convert Ensembl raw gene counts to gene symbols ####
# Load in data mapping ensembl gene IDs to gene symbols
geneinfo <- read.table("~/Documents/EPICC/Data/Expression/compiledGeneInfo.txt",header=T)
counts <- as.data.frame(fread("~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_counts.txt"))
merged <- merge(geneinfo,counts,by='GeneID');merged <- merge(merged,alllens,by='GeneID')
merged <- merged[,c('GeneID','Name','Length',colnames(merged)[grep('C\\d\\d\\d',colnames(merged))])]

# Assess duplicate gene names
alldups <- unique(merged[which(duplicated(merged$Name)),'Name']);epiccTMP <- merged
merged <- merged[which(merged$Name %ni% alldups),]
for(i in c(1:length(alldups))) {
  dupped <- epiccTMP[which(epiccTMP$Name==alldups[i]),]
  counts <- colSums(dupped[,c(3:ncol(merged))])
  len <- dupped[which(dupped$Length==max(dupped$Length)),'Length'][1]
  merged <- rbind(merged,c('Dupped',alldups[i],len,counts))
}
merged <- merged[,c(2:ncol(merged))];merged <- merged[order(merged$Name),]
for(col in c(2:ncol(merged))) { merged[,col] <- as.integer(merged[,col]) }

# Output 
symcounts <- merged;symcounts$GeneID <- symcounts$Name
symcounts <- symcounts[,c(ncol(symcounts),3:(ncol(symcounts)-1))]
write.table(symcounts,"~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_symbol_counts.txt",sep='\t',quote=F,row.names = F)

# 4. Convert gene symbol counts to TPM and log(TPM+1) ####
# a. Normalise for gene length: do gene counts / gene length (in kb)
symEPICC <- merged[,colnames(merged)[grep('C\\d\\d\\d',colnames(merged))]];row.names(symEPICC) <- merged$Name
for(i in c(1:ncol(symEPICC))) { symEPICC[,i] <- symEPICC[,i]/(merged$Length/1000) }

# b. Normalise for sequencing depth: sum the normalised gene counts and divide by a million
# then divide each normalised gene count by that samples scaling factor
for(i in c(1:ncol(symEPICC))) {
  sfactor <- sum(symEPICC[,i])/1000000
  symEPICC[,i] <- symEPICC[,i]/sfactor
}
symTPM <- symEPICC

# TPM to log(TPM+1)
symlogTPM <- TPMTologTpm(symTPM)

# Output TPM and logTPM files
# TPM
symTPM$GeneID <- row.names(symTPM);symTPM <- symTPM[,c(ncol(symTPM),1:(ncol(symTPM)-1))]
write.table(symTPM,"~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_symbol_tpm.txt",sep='\t',quote=F,row.names = F)
# logTPM
symlogTPM$GeneID <- row.names(symlogTPM);symlogTPM <- symlogTPM[,c(ncol(symlogTPM),1:(ncol(symlogTPM)-1))]
write.table(symlogTPM,"~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_symbol_logtpm.txt",sep='\t',quote=F,row.names = F)



