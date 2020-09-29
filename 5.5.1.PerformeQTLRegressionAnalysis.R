# Script to perform eQTL analysis by integrating genetic, epigenetic and expression data
# To do the following:
# 1) Load mutation and expression data
# 2) Load atac and copy number data
# 3) Get symbol and entrez ids for all genes
# 4) Filter, reformat and save matrices
# 5) Perform eQTL (linear regression) analysis on each gene

# Functions ####
# Extract the p-value for the overall F-statistic of a linear regression model
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# Main code ####

# 0. Prepare environment
setwd('~/Documents/EPICC/Data/QTLs/allvariants/')
library(readxl);library(stringr);library(data.table);library(MASS);library(rms)
library(dplyr);'%ni%' <- Negate('%in%');library(DESeq2);library(limma)
library(wesanderson)
convcn <- c('0'='Zero','1'='One','2'='Two','3'='Three','4'='Four')

# 1. Load mutation and expression data ####
# Load mut matrix of mutated genes
genemut <- readRDS('~/Documents/EPICC/Data/QTLs/allvariants/NS_mutburdengene_allDNAsamples.rds')
mutsam <- colnames(genemut)

# Load mut matrix of mutated enhancers
enhmut <- readRDS('~/Documents/EPICC/Data/QTLs/allvariants/mutburdenenhancer_allDNAsamples.rds')

# Load expression and get genes+samples that are shared
rnasam <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]
tsam <- rnasam[-grep('C\\d+_E',rnasam)];tsam <- tsam[-grep('C306',tsam)]
allpats <- gsub('(C\\d+)_\\S+','\\1',tsam);patients <- unique(allpats)
normsam <- rnasam[grep('^C\\d+_E',rnasam)];normsam <- normsam[which(gsub('(C\\d+)_\\S+','\\1',normsam) %in% patients)]

# vsd is used by TRACERx 
vsd <- readRDS('~/Documents/EPICC/Data/Expression/ProcessedCounts/filgenes.vsd.ensembl.rds')
geneexp <- as.data.frame(assay(vsd));normexp <- rowMeans(geneexp[,normsam]);normsd <- rowSds(as.matrix(geneexp[,normsam]))
allexp <- rowMeans(geneexp);allsd <- rowSds(as.matrix(geneexp))
for(i in c(1:ncol(geneexp))) {
  geneexp[,i] <- (geneexp[,i]-allexp)/allsd
}

# Add columns of 0s for normal samples
for(norm in normsam) { tmpnam <- colnames(genemut);genemut <- cbind(genemut,rep(0,nrow(genemut)));colnames(genemut) <- c(tmpnam,norm)}
for(norm in normsam) { tmpnam <- colnames(enhmut);enhmut <- cbind(enhmut,rep(0,nrow(enhmut)));colnames(enhmut) <- c(tmpnam,norm)}

matchsam <- sort(colnames(genemut)[which(colnames(genemut) %in% colnames(geneexp))])

# 2. Load ATAC-seq and CNA data ####
geneatac <- readRDS('~/Documents/EPICC/Data/QTLs/nodesamples/atacbygene_accessible.rds')
matchsam <- sort(matchsam[which(matchsam %in% colnames(geneatac))])
matchgene <- sort(row.names(geneatac)[which(row.names(geneatac) %in% row.names(geneexp))])

# Load CNAs and get genes+samples that are shared
genecna <- readRDS('~/Documents/EPICC/Data/QTLs/revisedmodel/cnabygene_numeric.rds')

# Add columns of 2s for normal samples
for(norm in normsam) { tmpnam <- colnames(genecna);genecna <- cbind(genecna,rep(2,nrow(genecna)));colnames(genecna) <- c(tmpnam,norm)}
cnasam <- colnames(genecna)

# Remove NA columns from genecna
matchsam <- sort(matchsam[which(matchsam %in% cnasam)])
genecna <- genecna[,matchsam]
nagenes <- apply(genecna, 1, function(x) sum(is.na(x)))
genecna <- genecna[which(nagenes==0),]
matchgene <- sort(matchgene[which(matchgene %in% row.names(genecna))])

# 3. Get symbol and entrez ids for every gene ####
# Filter out ENSGs that don't appear in biomart
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
geneMap = getBM(attributes=c("ensembl_gene_id","entrezgene_id","hgnc_symbol"), mart = ensembl,filters="ensembl_gene_id",values = matchgene)
dupens <- geneMap[which(duplicated(geneMap$ensembl_gene_id)),'ensembl_gene_id']
uniqmap <- geneMap[which(geneMap$ensembl_gene_id %ni% dupens),];row.names(uniqmap) <- c(1:nrow(uniqmap))
for(gene in dupens) {
  curdup <- geneMap[which(geneMap$ensembl_gene_id==gene),]
  uniqmap[nrow(uniqmap)+1,] <- c(curdup$ensembl_gene_id[1],paste(sort(unique(curdup$entrezgene_id)),collapse = ','),paste(sort(unique(curdup$hgnc_symbol)),collapse = ','))
}
geneDF <- uniqmap[order(uniqmap$ensembl_gene_id),];row.names(geneDF) <- c(1:nrow(geneDF))
matchgene <- matchgene[which(matchgene %in% geneDF$ensembl_gene_id)]

# 4. Filter and reformat all data type matrices and save as RData ####
# Filter all data types for genes and samples
segPat <- tmpPat <- genecna[matchgene,matchsam]
atacPat <- geneatac[matchgene,matchsam]
expPat <- as.matrix(geneexp[matchgene,matchsam])

addgens <- matrix(0L,nrow=length(matchgene[which(matchgene %ni% row.names(enhmut))]),ncol=ncol(enhmut))
row.names(addgens) <- matchgene[which(matchgene %ni% row.names(enhmut))];colnames(addgens) <- colnames(enhmut)
enhmut <- rbind(enhmut,addgens)
enhPat <- enhmut[matchgene,matchsam]

addgens <- matrix(0L,nrow=length(matchgene[which(matchgene %ni% row.names(genemut))]),ncol=ncol(genemut))
row.names(addgens) <- matchgene[which(matchgene %ni% row.names(genemut))];colnames(addgens) <- colnames(genemut)
genemut <- rbind(genemut,addgens)
mutPat <- genemut[matchgene,matchsam]

# Convert to characters for linear regression
for(i in c(1:ncol(mutPat))) { mutPat[,i] <- ifelse(mutPat[,i]==0,'WT','Mut')}
for(i in c(1:ncol(enhPat))) { enhPat[,i] <- ifelse(enhPat[,i]==0,'WT','Mut')}
for(i in c(1:ncol(atacPat))) { atacPat[,i] <- ifelse(atacPat[,i]==0,'Closed','Open')}
typetissue <- ifelse(gsub('C\\d+_(\\S)\\S+','\\1',colnames(atacPat))=='E','Normal','Tumour')

listmat <- list(expPat,mutPat,enhPat,segPat,atacPat,typetissue)
names(listmat) <- c('Expression','Mut','Enh','CNA','ATAC','Tissue')
saveRDS(listmat,'~/Documents/EPICC/Data/QTLs/allvariants/matricesOfModel.latest.cnanum.nolog.allvars.zscore.rds')

# 5. Perform eQTL analysis and save results as RData ####
tumsam <- matchsam[-grep('C\\d+_E',matchsam)]
patients <- unique(gsub('(C\\d+)_\\S+','\\1',tumsam))
allpats <- gsub('(C\\d+)_\\S+','\\1',colnames(expPat))

pvals <- c();modeltype <- c();varmod <- c();trackmod <- c();rsquared <- c();varps <- c()
for(gene in matchgene) {
  geneData <- data.frame(Expression=as.numeric(expPat[gene,]),Mut=mutPat[gene,],Enh=enhPat[gene,],CNA=segPat[gene,],ATAC=atacPat[gene,],Tissue=typetissue)
  keep <- c();for(i in c(2:(ncol(geneData)-1))) { if(length(table(geneData[,i]))>1) { keep <- c(keep,i) } }
  geneData <- geneData[,c(1,keep,ncol(geneData))]
  if('Mut' %in% colnames(geneData)) { geneData$Mut <- factor(geneData$Mut,levels=c('WT','Mut')) }
  if('Enh' %in% colnames(geneData)) { geneData$Enh <- factor(geneData$Enh,levels=c('WT','Mut')) }
  trackmod <- c(trackmod,paste(sort(colnames(geneData)[c(2:ncol(geneData))]),collapse = ' + '))
  res.lm <- lm(Expression ~., data = geneData)
  
  step <- stepAIC(res.lm, direction = 'backward', trace = F,k=3.8415)
  rsquared <- c(rsquared,summary(step)$adj.r.squared)
  
  if(as.character(step$terms)[3]==1) {
    pvals <- c(pvals,1)
    modeltype <- c(modeltype,'None')
    varmod <- c(varmod,'None')
    varps <- c(varps,'None')
  } else {
    pvals <- c(pvals,lmp(step))
    preds <- sort(str_split(as.character(step$terms)[3]," \\+ ")[[1]])
    modeltype <- c(modeltype,paste(preds,collapse = ' + '))
    
    # Record the regression coefficients for each linear regression separately
    covars <- rep(0,length(preds));names(covars) <- preds;pvars <- covars
    coefstep <- coef(step);indps <- summary(step)$coefficients[,'Pr(>|t|)']
    for(pred in preds) {
      covars[pred] <- as.numeric(coefstep[grep(paste0('^',pred),names(coefstep))])
      pvars[pred] <- indps[grep(paste0('^',pred),names(indps))]
    }
    varmod <- c(varmod,paste(paste(preds,covars,sep=':'),collapse = ' + '))
    varps <- c(varps,paste(paste(preds,pvars,sep=':'),collapse = ' + '))
    
  }
  
  if((which(matchgene==gene) %% 100)==0) {
    print(paste0('Analysing gene ',gene,' - ',which(matchgene==gene),'/',length(matchgene),' (',signif(which(matchgene==gene)/length(matchgene)*100,digits = 2),'%)'))
  }
}
modelres <- geneDF;modelres$MaxModel <- trackmod;modelres$FinalModel <- modeltype;modelres$VarMod <- varmod
modelres$VarP <- varps
modelres$Rsq <- rsquared;modelres$pval <- pvals;modelres$padj <- p.adjust(pvals,method='fdr')

saveRDS(modelres,file='~/Documents/EPICC/Data/QTLs/allvariants/qtl_model_res_latest.cnanum.nolog.allvars.zscore.rds')


