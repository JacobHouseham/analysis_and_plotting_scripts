# Script to perform GSEA on DEGs between regions
# To do the following:
# 1) Load expression and sample data
# 2) Pre set up dataframes for analysis input
# 3) Find DEGs and perform GSEA
# 4) Output the most commonly heterogeneous pathways

# Functions ####
getGAGE <- function(v,nl) {
  selG <- v$greater[, "q.val"] < 0.1 & !is.na(v$greater[, "q.val"]);v1 <- rownames(v$greater)[selG]
  selL <- v$less[, "q.val"] < 0.1 & !is.na(v$less[, "q.val"]);v2 <- rownames(v$less)[selL]
  path <- substr(c(v1, v2), 1, nl);return(path)
}

# Main code ####

# 0. Prepare environment
setwd("~/Documents/EPICC/Data/Expression")
suppressMessages(library(DESeq2));suppressMessages(library(pathview));suppressMessages(library(pcaExplorer));suppressMessages(library(gage));suppressMessages(library(gageData));data(go.sets.hs);data(go.subs.hs);data(kegg.gs)
suppressMessages(library(data.table));suppressMessages(library(ReactomePA));suppressMessages(library(biomaRt))
suppressMessages(library(clusterProfiler));suppressMessages(library(dplyr))
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(clusterProfiler)

# 1. Load sample info and raw counts ####
counts <- as.data.frame(fread("ProcessedCounts/All_EPICC_symbol_counts.txt",header=T))
samples <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]
tsam <- samples[which(!grepl('E1',samples))];nsam <- samples[which(grepl('E1',samples))]
counts <- counts[,c('GeneID',tsam)]
filgenes <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/filteredgenes.20perAL1TPM.symbol.txt')[,1]
counts <- counts[which(counts$GeneID %in% filgenes),];row.names(counts) <- c(1:nrow(counts))

# Convert to entrezid, since all GSEA seems to be in them
eg <- bitr(counts$GeneID,fromType = 'SYMBOL',toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg <- eg[-which(eg$ENTREZID=='100124696'),];row.names(eg) <- c(1:nrow(eg))

counts <- merge(eg,counts,by.x='SYMBOL',by.y='GeneID');row.names(counts) <- counts$ENTREZID

# Save conversion for later
conmat <- counts[c('SYMBOL','ENTREZID')]
saveRDS(conmat,file='~/Documents/EPICC/Data/Expression/GSEA/symbol_entrez_conversion.rds')

row.names(counts) <- counts$ENTREZID
counts <- counts[,tsam]
patients <- unique(gsub('(C\\d+)_\\S+','\\1',tsam))

# 2. Set up dataframes before running analysis ####
cNames <- c()
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
          cNames <- c(cNames,paste0(patient,'_',reg1,'vs',reg2))
        }
      }
    }
  }
}
kgDF <- data.frame(matrix(0L,nrow=length(names(kegg.gs)),ncol=length(cNames)));row.names(kgDF) <- sort(substr(names(kegg.gs),1,8));colnames(kgDF) <- cNames
mfDF <- data.frame(matrix(0L,nrow=length(names(go.sets.hs[go.subs.hs$MF])),ncol=length(cNames)));row.names(mfDF) <- sort(substr(names(go.sets.hs[go.subs.hs$MF]),1,10));colnames(mfDF) <- cNames
bpDF <- data.frame(matrix(0L,nrow=length(names(go.sets.hs[go.subs.hs$BP])),ncol=length(cNames)));row.names(bpDF) <- sort(substr(names(go.sets.hs[go.subs.hs$BP]),1,10));colnames(bpDF) <- cNames

patients <- unique(gsub('(C\\d+)_\\S+','\\1',cNames))

# 3. Find DEGs and do 3 separate GSEAs on them ####
setwd('~/Documents/EPICC/Data/Expression/GSEA')
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
    #vsd <- varianceStabilizingTransformation(dds,blind=T)
    #geneexp <- assay(vsd)
    
    # Get differentially expressed genes for any regions with 3 or more glands
    for(r1 in c(1:length(multireg))) {
      for(r2 in c(1:length(multireg))) {
        if(r1<r2) {
          reg1 <- names(multireg[r1]);reg2 <- names(multireg[r2]) # Define regions
          compname <- paste0(patient,'_',reg1,'vs',reg2)
          
          write(paste0('Finding GSEA of DEGs between regions ',reg1,' (n=',multireg[r1],') and ',reg2,' (n=',multireg[r2],')'), stdout())
          
          res <- results(dds,contrast=c('Region',reg1,reg2)) # Get DEG results
          res.fc <- res$log2FoldChange;names(res.fc) <- row.names(res);res.fc  <- sort(res.fc,decreasing = T) # Get fold changes
          sigres <- res[which(res$padj<0.05),]
          
          # Record the DEGs per comparison
          write(paste0(compname,' - ',nrow(sigres),' DEGs'), stdout())
          
          gse <- gseGO(geneList=res.fc,ont ="ALL",keyType = "ENTREZID",nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05,verbose=F,OrgDb = organism,pAdjustMethod = "fdr")
          # Get Enriched Up/Down molecular functions
          MF <- gse@result[which(gse@result$ONTOLOGY=='MF'),'ID']
          mfDF[row.names(mfDF) %in% MF,compname] <- 1 
          # Get Enriched Up/Down biological processes
          BP <- gse@result[which(gse@result$ONTOLOGY=='BP'),'ID']
          bpDF[row.names(bpDF) %in% BP,compname] <- 1
          
          # Get Enriched Up/Down kegg pathways
          kegg_organism = "hsa"
          kk2 <- gseKEGG(geneList=res.fc,organism=kegg_organism,nPerm=10000,minGSSize=3,maxGSSize=800,pvalueCutoff=0.05,verbose = F,pAdjustMethod="fdr",keyType="ncbi-geneid")
          kegg <- kk2@result$ID
          kgDF[row.names(kgDF) %in% kegg,compname] <- 1
          
          write(paste0(compname,' - ',length(kegg),' KEGG, ',length(MF),' MF, ',length(BP),' BP and ',length(pathsig),' Reactome significant pathways found'), stdout())
        }
      }
    }
  } else {
    write('No regions with enough samples', stdout())
  }
}

gseaDFs <- list(kgDF,mfDF,bpDF,reactDF);names(gseaDFs) <- c('KEGG','MF','BP','Reactome')
saveRDS(gseaDFs,file='~/Documents/EPICC/Data/Expression/GSEA/gsea_res_list.newtry.rds')

# 4. Determine most common GSEAs and output as csvs ####
gseaDFs <- readRDS(file='~/Documents/EPICC/Data/Expression/GSEA/gsea_res_list.newtry.rds')

commonKG <- data.frame(ID=substr(names(kegg.gs),1,8),Description=substr(names(kegg.gs),9,nchar(names(kegg.gs))),Common=rowSums(gseaDFs$KEGG))
commonKG <- arrange(commonKG,desc(Common))
write.csv(commonKG,file='~/Documents/Thesis/mainText/commonKG.csv',row.names=F,quote=T)
# Molecular Function
commonMF <- data.frame(ID=substr(names(go.sets.hs[go.subs.hs$MF]),1,10),Description=substr(names(go.sets.hs[go.subs.hs$MF]),11,nchar(names(go.sets.hs[go.subs.hs$MF]))),Common=rowSums(gseaDFs$MF))
commonMF <- arrange(commonMF,desc(Common))
commonMF[which(commonMF$ID=='GO:0000976'),'Description'] <- 'transcription regulatory region DNA binding'
commonMF[which(commonMF$ID=='GO:0000977'),'Description'] <- 'RNA Pol II regulatory region DNA binding'
commonMF[which(commonMF$ID=='GO:0000978'),'Description'] <- 'RNA Pol II core promoter proximal region DNA binding'
write.csv(commonMF,file='~/Documents/Thesis/mainText/commonMF.csv',row.names=F,quote=T)
# Biological Processes
commonBP <- data.frame(ID=substr(names(go.sets.hs[go.subs.hs$BP]),1,10),Description=substr(names(go.sets.hs[go.subs.hs$BP]),11,nchar(names(go.sets.hs[go.subs.hs$BP]))),Common=rowSums(gseaDFs$BP))
commonBP <- arrange(commonBP,desc(Common))
write.csv(commonBP,file='~/Documents/Thesis/mainText/commonBP.csv',row.names=F,quote=T)





