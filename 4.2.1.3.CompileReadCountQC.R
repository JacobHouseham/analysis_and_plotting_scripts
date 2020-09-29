# Script for loading in all count data and assessing QC 
# To do the following:
# 1) Load raw gene count and replicate QC data to assess which files should be used
# 2) Load and process all raw gene counts
# 3) Assess the QC of all samples
# 4) Output a filtered gene count matrix for all samples
# 5) Output a QC table for all samples, and a list of only PASS samples

# Main code ####

# 0. Prepare environment
setwd("~/Documents/EPICC/Data/Expression/")
'%ni%' <- Negate('%in%');nucs <- c('A','C','G','T')
library(dplyr);library(data.table);library(stringr)

# 1. Load in data from 4.2.1.2 and use it to assess which raw gene count files to use ####
# Get file list and duplicate sample list
filelist <- readRDS('~/Documents/ThesisOther/ScriptsForThesis/RData/tmp.currentreadcountfileloc.rds')
oldnew <- readRDS('~/Documents/ThesisOther/ScriptsForThesis/RData/reranQCDecision.rds')

# Based on replicate QC decision generate a final sample list
sampleList <- data.frame(Sample=unique(filelist$Sample),File=rep('',length(unique(filelist$Sample))))
for(i in c(1:nrow(sampleList))) {
  cursam <- sampleList[i,]
  if(cursam$Sample %in% oldnew$Sample) {
    curoldnew <- oldnew[which(oldnew$Sample==cursam$Sample),]
    if(curoldnew$Decision=='Fail') {
      maxsam <- names(which.max(curoldnew[,c(4,5)]))
      chooseSam <- paste0(cursam$Sample,'_',curoldnew[,paste0(maxsam,'ID')])
      sampleList[i,'File'] <- filelist[which(filelist$SampleID==chooseSam),'FileUse']
    } else if(curoldnew$Decision=='UseNew') {
      chooseSam <- paste0(cursam$Sample,'_',curoldnew[,'NewID'])
      sampleList[i,'File'] <- filelist[which(filelist$SampleID==chooseSam),'FileUse']
    } else if(curoldnew$Decision=='UseOld') {
      chooseSam <- paste0(cursam$Sample,'_',curoldnew[,'OldID'])
      sampleList[i,'File'] <- filelist[which(filelist$SampleID==chooseSam),'FileUse']
    } else if(curoldnew$Decision=='Combine') {
      sam1 <- paste0(cursam$Sample,'_',curoldnew[,'OldID'])
      sam2 <- paste0(cursam$Sample,'_',curoldnew[,'NewID'])
      sampleList[i,'File'] <- paste0(filelist[which(filelist$SampleID==sam1),'FileUse'],'&&',filelist[which(filelist$SampleID==sam2),'FileUse'])
    }
  } else {
    sampleList[i,'File'] <- filelist[which(filelist$Sample==cursam$Sample),'FileUse']
  }
}

# 2. Compile all sample raw gene counts into one data frame ####
# Load the first sample
firstsam <- as.data.frame(fread(gsub('\\S+(EPICC_\\S+)','rawCounts\\/\\1',sampleList$File[1])))
EPICC <- firstsam[c(5:(nrow(firstsam))),c('V1','V4')];colnames(EPICC) <- c('GeneID',sampleList$Sample[1])

# Append all other samples onto this data frame
for(i in c(2:nrow(sampleList))) {
  print(paste0('Now loading sample ',sampleList[i,'Sample'],' ',i,'/',nrow(sampleList),' (',signif(i/nrow(sampleList)*100,3),'%)'))
  # If there are multiple files for one sample (i.e. replicates to be combined)
  # Add the counts from both files together
  if(length(grep('&&',sampleList[i,'File']))>0) {
    print(paste0('Combining counts for sample ',sampleList[i,'Sample']))
    file1 <- str_split(sampleList[i,'File'],'&&')[[1]][1];file2 <- str_split(sampleList[i,'File'],'&&')[[1]][2]
    tmpsam <- as.data.frame(fread(gsub('\\S+(EPICC_\\S+)','rawCounts\\/\\1',file1)));tmpsam2 <- as.data.frame(fread(gsub('\\S+(EPICC_\\S+)','rawCounts\\/\\1',file2)))
    tmpsam$V4 <- tmpsam$V4+tmpsam2$V4
  # Load normally if just one file for current sample
  } else {
    tmpsam <- as.data.frame(fread(gsub('\\S+(EPICC_\\S+)','rawCounts\\/\\1',sampleList[i,'File'])))
  }
  EPICC[,sampleList[i,'Sample']] <- tmpsam[c(5:(nrow(tmpsam))),'V4']
}

# Filter and reorder gene list
EPICC <- EPICC[-grep("PAR_Y",EPICC$GeneID),]
EPICC$GeneID <- gsub('(ENSG\\d+)\\.\\d+','\\1',EPICC$GeneID)
EPICC <- EPICC[order(EPICC$GeneID),]

# Save raw gene count data frame containing all samples
write.table(EPICC,'~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_counts.allgenes.txt',sep='\t',row.names = F,quote = F)

# 3. Load gene info and make QC table of total counts and counts by biotype ####
EPICC <- as.data.frame(fread('~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_counts.allgenes.txt'))
geneinfo <- read.table("~/Documents/EPICC/Data/Expression/compiledGeneInfo.txt",header=T)
ribo_genes <- read.table('~/Documents/EPICC/Data/Expression/GSEA/ribosomal_genes_hsa03010.txt')[,1]
EPICCcounts <- merge(geneinfo,EPICC,by='GeneID')
mitos <- c();prots <- c();immos <- c();snos <- c();others <- c()
for(sam in sampleList$Sample) {
  sortCount <- EPICCcounts[order(EPICCcounts[,sam],decreasing=T),]
  mitos <- c(mitos,sum(sortCount[which(sortCount$Type=='Mt_rRNA' | sortCount$Name %in% ribo_genes),sam]))
  prots <-c(prots,sum(sortCount[which(sortCount$Type=='protein_coding'),sam]))
  others <- c(others,sum(sortCount[which(sortCount$Type %ni% c('Mt_rRNA','protein_coding') & sortCount$Name %ni% ribo_genes),sam]))
}
QCtable <- EPICC[,sampleList$Sample]
QCtable <- as.data.frame(t(QCtable[c(1:2),]));colnames(QCtable) <- c('Assigned','ProteinCoding')
QCtable$Assigned <- colSums(as.matrix(EPICC[,sampleList$Sample]))/1e6
QCtable$UsableReads <- colSums(as.matrix(EPICCcounts[which(EPICCcounts$Chr %in% paste0('chr',c(1:22,'X','Y')) & EPICCcounts$Name %ni% ribo_genes & EPICCcounts$Type=='protein_coding'),sampleList$Sample]))/1e6
QCtable$ProteinCoding <- prots/1e6;QCtable$PercentProteinCoding <- QCtable$ProteinCoding/QCtable$Assigned*100
QCtable$MTrRNA <- mitos/1e6;QCtable$PercentMTrRNA <- QCtable$MTrRNA/QCtable$Assigned*100
QCtable$Other <- others/1e6;QCtable$PercentOther <- QCtable$Other/QCtable$Assigned*100

# 4. Write processed raw gene counts file ####
# Write file with only protein-coding non-ribosomal genes from canonical chromosomes
EPICC <- EPICC[which(EPICC$GeneID %in% geneinfo[which(geneinfo$Type=='protein_coding' & geneinfo$Chr %in% paste0('chr',c(1:22,'X','Y'))& geneinfo$Name %ni% ribo_genes),'GeneID']),];row.names(EPICC) <- c(1:nrow(EPICC))
write.table(EPICC,'~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_counts.txt',sep='\t',row.names = F,quote = F)

# 5. Reformat and output QC table and list of only PASS samples ####
QCtable$Sample <- row.names(QCtable)
qcFinal <- QCtable[,c('Sample','Assigned','ProteinCoding','PercentProteinCoding','MTrRNA','PercentMTrRNA','Other','PercentOther','UsableReads')]

# Samples PASS QC if Number of million usable reads >=5 
qcFinal$QC <- ifelse(qcFinal$UsableReads>=5,'PASS','FAIL')

# Print out QCtable as .csv file
write.csv(qcFinal,file='~/Documents/ThesisOther/ScriptsForThesis/Tables/Final_EPICC_QC_ReadCounts.csv',row.names = F)
# Write file of PASS samples
write.table(qcFinal[which(qcFinal$QC=='PASS'),'Sample'],file='~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt',row.names = F,col.names = F,quote=F)


