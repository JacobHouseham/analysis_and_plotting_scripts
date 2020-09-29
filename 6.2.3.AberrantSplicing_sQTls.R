# Script to assess if mutations in splice sites affect the splicing of those genes
# To do the following:
# 1) Load expression data and aberrant splicing data
# 2) Load splice site mutation data and combine with aberrant splicing data
# 3) C554 analysis of splicing mutations and splicing
# 4) Plot C554 analysis
# 5) C516 analysis of splicing mutations and splicing
# 6) Plot C516 analysis
# 7) All sample analysis of splicing mutations and splicing
# 8) Plot analysis of all samples together

# Main code ####

# 0. Prepare environment
library(data.table);setwd("~/Documents/Splicing");load("~/Documents/EPICC/Data/Splicing/refinedGFF.RData")
"%ni%" <- Negate("%in%");library(pheatmap);library(RColorBrewer)
geneinfo <- read.table("~/Documents/EPICC/Data/Expression/compiledGeneInfo.txt",header=T)
genehg38 <- merge(geneinfo,EPICC,by='GeneID')
tsgs <- read.table('~/Documents/EPICC/Data/Splicing/TSGs300.txt')[,1]

# 1. Load expression and splicing data ####
EPICC <- as.data.frame(fread(paste0("~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_tpm.txt"),header=T))
rnasam <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]

load(file='~/Documents/ThesisOther/ScriptsForThesis/RData/aberrant_splicing_data.Rdata')

allsam <- sort(unique(c(samlist$Samples,samlist$MatchedNorm)))
patients <- sort(unique(gsub('(C\\d+)_.+','\\1',samlist$MatchedNorm)))
EPICC <- EPICC[,c('GeneID',allsam)]

# 2. Load pre-compiled splice-site mutation data and combine with aberrant splicing data ####
splicemut <- readRDS('~/Documents/EPICC/Data/QTLs/allvariants/splice_mutburdengene_allDNAsamples.rds')
usesam <- c();allevents <- c()
for(i in c(1:nrow(samlist))) {
  sam <- samlist$Samples[i];normal <- samlist$MatchedNorm[i]
  cursig <- spliceList[[sam]]
  
  if(sam %in% colnames(splicemut)) {
    print(paste0(sam,' : ',nrow(cursig),' sig splicing events, ',
                length(which(splicemut[,sam]!=0)),' splice Muts, ',
                length(which(cursig$GeneID %in% row.names(splicemut[which(splicemut[,sam]!=0),]))),' sig splicing events have Muts'))
    if(length(which(cursig$GeneID %in% row.names(splicemut[which(splicemut[,sam]!=0),])))>0) {
      print(cursig$Name[which(cursig$GeneID %in% row.names(splicemut[which(splicemut[,sam]!=0),]))])
    }
    usesam <- c(usesam,sam);allevents <- sort(unique(c(allevents,cursig$Name[which(cursig$GeneID %in% row.names(splicemut[which(splicemut[,sam]!=0),]))])))
  }
}

# 3. For tumour C554, assess the association of splice site mutations and aberrant splicing ####
pat <- 'C554'
sqtl <- matrix(-1,nrow=112,ncol=length(grep(pat,usesam)));row.names(sqtl) <- allevents
colnames(sqtl) <- usesam[grep(pat,usesam)]
ensgs <- row.names(sqtl);names(ensgs)  <- geneinfo[match(ensgs,geneinfo$Name),'GeneID']

for(i in c(1:ncol(sqtl))) {
  sam <- colnames(sqtl)[i];normal <- samlist[which(samlist$Samples==sam),'MatchedNorm']
  cursig <- spliceList[[sam]]
  expgenes <- row.names(sqtl)[which(row.names(sqtl) %in% unique(genehg38[which(genehg38[,sam]>=10 & genehg38[,normal]>=10),'Name']))]
  mutgenes <- as.character(ensgs[names(which(splicemut[,sam]!=0 & row.names(splicemut) %in% names(ensgs)))])
  asgenes <- row.names(sqtl)[which(row.names(sqtl) %in% cursig$Name)]
  # Give 0 to genes that are available
  sqtl[expgenes,sam] <- 0
  # Give 1 to genes that are aberrantly spliced
  sqtl[asgenes,sam] <- 1
  # Give 2 to genes that are only mutated
  sqtl[mutgenes[which(mutgenes %ni% expgenes)],sam] <- 2
  # Give 3 to genes that are mutated and available but not AS
  sqtl[mutgenes[which(mutgenes %in% expgenes & mutgenes %ni% asgenes)],sam] <- 3
  # Give 4 to genes that are mutated and AS
  sqtl[mutgenes[which(mutgenes %in% asgenes)],sam] <- 4
}
sqtl <- sqtl[which(rowSums(sqtl>=2)>=1),]

which(rowSums(sqtl==4)>=1 & (rowSums(sqtl==-1)+rowSums(sqtl==0)+rowSums(sqtl==2))>=1 & rowSums(sqtl==1)==0 & rowSums(sqtl==3)==0)

# 4. Plot heatmap of C554 splice mutation vs. splicing analysis ####
pdf('~/Documents/Thesis/figures/sQTLs_C554.pdf')
sqtldata <- data.frame(matrix(nrow=ncol(sqtl),ncol=1));colnames(sqtldata) <- c('Region')
row.names(sqtldata) <- colnames(sqtl)
sqtldata$Region <- as.factor(gsub('C\\d+_(\\S)\\d+\\S+','\\1',row.names(sqtldata)))
my_colour <- list(Region=c(A='#E31A1C',B='#377DB8',C='#4DAE49',D='#904A9A'))
breaksList = seq(-1, 4, by = 1)
pheatmap(sqtl,show_rownames=T,treeheight_col = 25,legend_breaks = c(-0.5833,0.25,1.0833,1.9167,2.7500,3.5833,4),legend_labels = c('NA','Avail (Not Mut)','AS (Not Mut)','Mut (Not Avail)','Mut (Avail)','Mut + AS','Event'),color=c('gray95',wes_palette("Darjeeling1")[c(5,2)],wes_palette("BottleRocket2")[1],wes_palette("Darjeeling1")[c(3,1)]),cluster_rows = T,treeheight_row = 25,annotation_col = sqtldata,annotation_colors = my_colour)
dev.off()

# 5. For tumour C516, assess the association of splice site mutations and aberrant splicing ####
pat <- 'C516'
sqtl <- matrix(-1,nrow=112,ncol=length(grep(pat,usesam)));row.names(sqtl) <- allevents
colnames(sqtl) <- usesam[grep(pat,usesam)]
ensgs <- row.names(sqtl);names(ensgs)  <- geneinfo[match(ensgs,geneinfo$Name),'GeneID']

for(i in c(1:ncol(sqtl))) {
  sam <- colnames(sqtl)[i];normal <- samlist[which(samlist$Samples==sam),'MatchedNorm']
  cursig <- spliceList[[sam]]
  expgenes <- row.names(sqtl)[which(row.names(sqtl) %in% unique(genehg38[which(genehg38[,sam]>=10 & genehg38[,normal]>=10),'Name']))]
  mutgenes <- as.character(ensgs[names(which(splicemut[,sam]!=0 & row.names(splicemut) %in% names(ensgs)))])
  asgenes <- row.names(sqtl)[which(row.names(sqtl) %in% cursig$Name)]
  # Give 0 to genes that are available
  sqtl[expgenes,sam] <- 0
  # Give 1 to genes that are aberrantly spliced
  sqtl[asgenes,sam] <- 1
  # Give 2 to genes that are only mutated
  sqtl[mutgenes[which(mutgenes %ni% expgenes)],sam] <- 2
  # Give 3 to genes that are mutated and available but not AS
  sqtl[mutgenes[which(mutgenes %in% expgenes & mutgenes %ni% asgenes)],sam] <- 3
  # Give 4 to genes that are mutated and AS
  sqtl[mutgenes[which(mutgenes %in% asgenes)],sam] <- 4
}
sqtl <- sqtl[which(rowSums(sqtl>=2)>=1),]

which(rowSums(sqtl==4)>=1 & (rowSums(sqtl==-1)+rowSums(sqtl==0)+rowSums(sqtl==2))>=1 & rowSums(sqtl==1)==0 & rowSums(sqtl==3)==0)

# 6. Plot heatmap of C516 splice mutation vs. splicing analysis ####
pdf('~/Documents/Thesis/figures/sQTLs_C516.pdf')
sqtldata <- data.frame(matrix(nrow=ncol(sqtl),ncol=1));colnames(sqtldata) <- c('Region')
row.names(sqtldata) <- colnames(sqtl)
sqtldata$Region <- as.factor(gsub('C\\d+_(\\S)\\d+\\S+','\\1',row.names(sqtldata)))
my_colour <- list(Region=c(A='#E31A1C',B='#377DB8',C='#4DAE49',D='#904A9A'))
breaksList = seq(-1, 4, by = 1)
pheatmap(sqtl,show_rownames=T,treeheight_col = 25,legend_breaks = c(-0.5833,0.25,1.0833,1.9167,2.7500,3.5833,4),legend_labels = c('NA','Avail (Not Mut)','AS (Not Mut)','Mut (Not Avail)','Mut (Avail)','Mut + AS','Event'),color=c('gray95',wes_palette("Darjeeling1")[c(5,2)],wes_palette("BottleRocket2")[1],wes_palette("Darjeeling1")[c(3,1)]),cluster_rows = T,treeheight_row = 25,annotation_col = sqtldata,annotation_colors = my_colour)
dev.off()

# 7. For all samples across tumours, assess the association of splice site mutations and aberrant splicing ####
sqtl <- matrix(-1,nrow=112,ncol=length(usesam));row.names(sqtl) <- allevents
colnames(sqtl) <- usesam
ensgs <- row.names(sqtl);names(ensgs)  <- geneinfo[match(ensgs,geneinfo$Name),'GeneID']

for(i in c(1:ncol(sqtl))) {
  sam <- colnames(sqtl)[i];normal <- samlist[which(samlist$Samples==sam),'MatchedNorm']
  cursig <- spliceList[[sam]]
  expgenes <- row.names(sqtl)[which(row.names(sqtl) %in% unique(genehg38[which(genehg38[,sam]>=10 & genehg38[,normal]>=10),'Name']))]
  mutgenes <- as.character(ensgs[names(which(splicemut[,sam]!=0 & row.names(splicemut) %in% names(ensgs)))])
  asgenes <- row.names(sqtl)[which(row.names(sqtl) %in% cursig$Name)]
  # Give 0 to genes that are available
  sqtl[expgenes,sam] <- 0
  # Give 1 to genes that are aberrantly spliced
  sqtl[asgenes,sam] <- 1
  # Give 2 to genes that are only mutated
  sqtl[mutgenes[which(mutgenes %ni% expgenes)],sam] <- 2
  # Give 3 to genes that are mutated and available but not AS
  sqtl[mutgenes[which(mutgenes %in% expgenes & mutgenes %ni% asgenes)],sam] <- 3
  # Give 4 to genes that are mutated and AS
  sqtl[mutgenes[which(mutgenes %in% asgenes)],sam] <- 4
}
sqtl <- sqtl[which(rowSums(sqtl>=2)>=1),]

which(rowSums(sqtl==4)>=1 & (rowSums(sqtl==-1)+rowSums(sqtl==0)+rowSums(sqtl==2))>=1 & rowSums(sqtl==1)==0 & rowSums(sqtl==3)==0)

# 8. Plot heatmap of all splice mutation vs. splicing results ####
pdf('~/Documents/Thesis/figures/sQTLs_allsamples.pdf',width=8,height=10)
sqtldata <- data.frame(matrix(nrow=ncol(sqtl),ncol=1));colnames(sqtldata) <- c('Patient')
row.names(sqtldata) <- colnames(sqtl)
sqtldata$Patient <- as.factor(gsub('(C\\d+)_\\S\\d+\\S+','\\1',row.names(sqtldata)))
patcols <- c(wes_palette("Moonrise3"),wes_palette("Moonrise2")[2]);names(patcols) <- unique(sqtldata$Patient)
my_colour <- list(Region=c(A='#E31A1C',B='#377DB8',C='#4DAE49',D='#904A9A'),Patient=patcols)
breaksList = seq(-1, 4, by = 1)
pheatmap(sqtl,show_rownames=T,fontsize_row = 6,fontsize_col = 6,cluster_cols = F,cluster_rows = T,treeheight_col = 25,legend_breaks = c(-0.5833,0.25,1.0833,1.9167,2.7500,3.5833,4),legend_labels = c('NA','Avail (Not Mut)','AS (Not Mut)','Mut (Not Avail)','Mut (Avail)','Mut + AS','Event'),color=c('gray95',wes_palette("Darjeeling1")[c(5,2)],wes_palette("BottleRocket2")[1],wes_palette("Darjeeling1")[c(3,1)]),treeheight_row = 25,annotation_col = sqtldata,annotation_colors = my_colour)
dev.off()






