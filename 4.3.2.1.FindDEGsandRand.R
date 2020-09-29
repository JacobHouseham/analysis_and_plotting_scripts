# Script to get the number of DEGs between regions
# To do the following:
# 1) Load expression data and filter for expressed genes
# 2) Assess clonal mixing with different filtered gene lists
# 3) For regions with multiple samples find the number of DEGs
# and also randomise sample names and repeat analysis
# 4) Save results
# 5) Do Reactome enrichment on highly heterogeneous and homogenous genes

# Main code ####

# 0. Prepare environment
setwd("~/Documents/EPICC/Data/Expression")
library(DESeq2);library(data.table);library(biomaRt);library(dplyr)

# 1. Load expression data and filter for expressed genes ####
# Load raw gene counts
counts <- as.data.frame(fread("ProcessedCounts/All_EPICC_counts.txt",header=T))
# Get PASS samples
samples <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]
tsam <- samples[which(!grepl('E1',samples))];nsam <- samples[which(grepl('E1',samples))]
counts <- counts[,c('GeneID',tsam)];row.names(counts) <- counts$GeneID
# Filter for expressed genes and only consider tumour samples
filgenes <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/filteredgenes.20perAL1TPM.ensembl.txt')[,1]
counts <- counts[which(counts$GeneID %in% filgenes),tsam]

# 2. Set up dataframes prior to analysis ####
cNames <- c();patients <- unique(gsub('(C\\d+)_\\S+','\\1',tsam))
for(patient in patients) {
  patsam <- colnames(counts)[grep(patient,names(counts))]
  patregs <- gsub("C\\d+_(\\S)\\d+_.+$","\\1",patsam)
  regions <- table(patregs)
  multireg <- regions[which(regions>=2)]
  if(length(multireg)>=2) {
    for(r1 in c(1:length(multireg))) {
      for(r2 in c(1:length(multireg))) {
        if(r1<r2) {
          reg1 <- names(multireg[r1]);reg2 <- names(multireg[r2])
          cNames <- c(cNames,paste0(patient,'_',reg1,multireg[r1],'vs',reg2,multireg[r2]))
        }
      }
    }
  }
}
degDiff <- data.frame(matrix(0L,nrow=nrow(counts),ncol=length(cNames)));row.names(degDiff) <- row.names(counts);colnames(degDiff) <- cNames
randDiff <- degDiff

# 3. Use DESeq2 to get DEGs between regions + reorder sample names to get random ####
setwd('~/Documents/EPICC/Data/Expression/GSEA')
patients <- unique(gsub('(C\\d+)_\\S+','\\1',cNames))
for(patient in patients) {
  write(paste0('Analysing expression heterogeneity of patient ',patient), stdout())
  
  # Filter for patient and convert ready for 'count matrix' input for Deseq2 analysis
  EPICC <- counts[,grep(patient,names(counts))]
  EPICCdata <- as.data.frame(t(EPICC[c(1:2),]));colnames(EPICCdata) <- c('Region','Type')
  EPICCdata$Region <- as.factor(gsub("C\\d+_(\\S)\\d+_.+$","\\1",row.names(EPICCdata)))
  EPICCdata$Type <- as.factor(sapply(row.names(EPICCdata),function(x) if(gsub("C\\d+_\\S\\d+_\\S.+$","\\1",x)=='B') {x<-'Bulk'} else {x<-'Gland'}))
  
  # DEGs ####
  regions <- table(EPICCdata$Region)
  multireg <- regions[which(regions>=2)]
  if(length(multireg)>=2) {
    # Run Deseq2 to get dds and normalise
    write('Performing DESeq2 analysis', stdout())
    dds <- DESeqDataSetFromMatrix(countData = EPICC,colData = EPICCdata,design = ~ Region)
    suppressMessages(dds <-DESeq(dds))
  
    # Get differentially expressed genes for any regions with 2 or more glands
    for(r1 in c(1:length(multireg))) {
      for(r2 in c(1:length(multireg))) {
        if(r1<r2) {
          reg1 <- names(multireg[r1]);reg2 <- names(multireg[r2]) # Define regions
          compname <- paste0(patient,'_',reg1,multireg[r1],'vs',reg2,multireg[r2])
          
          write(paste0('Finding GSEA of DEGs between regions ',reg1,' (n=',multireg[r1],') and ',reg2,' (n=',multireg[r2],')'), stdout())
          # Get DEG results
          res <- results(dds,contrast=c('Region',reg1,reg2))
          # Filter for significant DEGs
          sigres <- res[which(res$padj<0.05),]
          
          # Record the DEGs per comparison
          write(paste0(compname,' - ',nrow(sigres),' DEGs'), stdout())
          degDiff[row.names(degDiff) %in% row.names(sigres),compname] <- 1
        }
      }
    }
      
    Random <- counts[,grep(patient,names(counts))]
    Randomdata <- as.data.frame(t(Random[c(1:2),]));colnames(Randomdata) <- c('Region','Type')
    Randomdata$Type <- as.factor(sapply(row.names(Randomdata),function(x) if(gsub("C\\d+_\\S\\d+_\\S.+$","\\1",x)=='B') {x<-'Bulk'} else {x<-'Gland'}))
    Randomdata$Region <- as.factor(sample(gsub("C\\d+_(\\S)\\d+_.+$","\\1",row.names(EPICCdata))))
      
    while(length(which(EPICCdata$Region==Randomdata$Region))==nrow(EPICCdata) | length(which(EPICCdata$Region==rev(Randomdata$Region)))==nrow(EPICCdata)) {
      Randomdata$Region <- as.factor(sample(gsub("C\\d+_(\\S)\\d+_.+$","\\1",row.names(EPICCdata))))
    }
      
    write('Performing DESeq2 analysis on random samples', stdout())
    dds <- DESeqDataSetFromMatrix(countData = Random,colData = Randomdata,design = ~ Region)
    suppressMessages(dds <-DESeq(dds))
      
    for(r1 in c(1:length(multireg))) {
      for(r2 in c(1:length(multireg))) {
        if(r1<r2) {
          reg1 <- names(multireg[r1]);reg2 <- names(multireg[r2]) # Define regions
          compname <- paste0(patient,'_',reg1,multireg[r1],'vs',reg2,multireg[r2])
          res <- results(dds,contrast=c('Region',reg1,reg2)) # Get DEG results
          sigres <- res[which(res$padj<0.05),]
          randDiff[row.names(randDiff) %in% row.names(sigres),compname] <- 1
          write(paste0(compname,' randomized - ',nrow(sigres),' DEGs'), stdout())
        }
      }
    }
  } else {
    write('No regions with enough samples', stdout())
  }
}

# 4. Save RData for use in subsequent scripts ####
saveRDS(degDiff,file='~/Documents/EPICC/Data/Expression/GSEA/degDiff.rds')
saveRDS(randDiff,file='~/Documents/EPICC/Data/Expression/GSEA/randDiff.rds')

