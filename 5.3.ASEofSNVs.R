# Script to assess the ASE of SNVs by comparing the allele counts in WGS and RNA-seq
# To do the following:
# 1) Load expression and sample data
# 2) Determine the ASE of SNVs within tumours
# 3) Plot the ASE of SNVs
# 4) Look at case studies and output summary file

# Main code ####

# 0. Prepare environment
setwd("~/Documents/EPICC/Data/Mutations/VCFtoTSV")
'%ni%' <- Negate('%in%');nucs <- c('A','C','G','T')
library(dplyr);library(data.table);library(stringr);library(wesanderson)
colrm <- c('Feature_type','cDNA_position','CDS_position','Protein_position','miRNA','DISTANCE','FLAGS','VARIANT_CLASS','SYMBOL_SOURCE','HGNC_ID','TSL','APPRIS','CCDS','SWISSPROT','TREMBL','UNIPARC','GENE_PHENO','DOMAINS','HGVS_OFFSET',
           'AFR_AF','AMR_AF','EAS_AF','EUR_AF','SAS_AF','AA_AF','EA_AF','gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_ASJ_AF','gnomAD_EAS_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_OTH_AF','gnomAD_SAS_AF','MAX_AF','MAX_AF_POPS','AF','SOMATIC','PHENO','PUBMED')

# 1. Load samples names and gene expression ####
wgssam <- readRDS('~/Documents/EPICC/Data/QTLs/nodesamples/wgs_used_in_phylo.rds')
wgssam <- gsub('EPICC_(C\\d+_\\S\\d+_\\S\\d+)\\S+','\\1',wgssam)
rnasam <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]
samples <- wgssam[which(wgssam %in% rnasam)]

geneexp <- as.data.frame(fread('~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_symbol_tpm.txt',stringsAsFactors = F))

# 2. Determine ASE of WGS-based SNVs comparing allele counts in WGS and RNA-seq ####
snvASEres <- matrix(0L,nrow=3,ncol=length(samples))#;colnames(snvASEres) <- samples
row.names(snvASEres) <- c('LowCov','Non-ASE','ASE')
combList <- depthList <- sigList <- list()
par(mar=c(4.6,4.6,3.6,5.1),font=2,font.axis=2,font.lab=2,cex.axis=1.5,cex.lab=1.5)
for(s in c(1:length(samples))) { 
  sam <- samples[s]
  combmut <- as.data.frame(fread(paste0('~/Documents/EPICC/Data/Mutations/VCFtoTSV/rna_matched/',sam,'.txt')))
  combmut <- combmut[,which(colnames(combmut) %ni% colrm)]
  combmut$Locus <- paste0(combmut$Chr,':',combmut$Pos,'_',combmut$Ref,'/',combmut$Alt)
  depthmut <- combmut[which(combmut$NR_RNA>=10),]
  
  if(nrow(depthmut)!=0) {
    fishp <- c()
    for(i in c(1:nrow(depthmut))) {
      mut <- depthmut[i,]
      mat <- rbind(c(mut$NR-mut$NV,mut$NV),c(mut$NR_RNA-mut$NV_RNA,mut$NV_RNA))
      fishp <- c(fishp,fisher.test(mat)$p.value)
    }
    depthmut$pval_fish <- fishp
    depthmut$fdr_fish <- p.adjust(depthmut$pval_fish,method='fdr')
    
    sigmut <- depthmut[which(depthmut$fdr_fish<0.01),]
    
  } else {
    sigmut <- depthmut
  }
  snvASEres[,s] <- c(length(which(combmut$Locus %ni% depthmut$Locus)),length(which(depthmut$Locus %ni% sigmut$Locus)),nrow(sigmut))/nrow(combmut)*100
  combList[[sam]] <- combmut;depthList[[sam]] <- depthmut;sigList[[sam]] <- sigmut
}

# 3. Plot summary of ASE of SNVs ####

pdf('~/Documents/Thesis/figures/ASEinSNVs.pdf',height=9)
par(mar=c(3.1,4.1,0.1,6.1),font=2,font.axis=2,font.lab=2,las=1)
yy <- barplot(snvASEres[c(3:1),rev(which(snvASEres[1,]!=100))],horiz = T,axes=F,cex.names=0.5,col=wes_palette("Zissou1")[c(5,3,1)],border=NA)
axis(side=1,at=seq.int(0,100,by=10),labels=seq.int(0,100,by=10),cex=0.75,line=-.75,lwd=2)
mtext(side=1,line=1.5,text='% of matched DNA-RNA SNVs',cex=1)
mtext(side=2,at=yy,text=samples[rev(which(snvASEres[1,]!=100))],cex=0.65)
legend(x=98,y=55,legend=c('Low Cov','Non-ASE','ASE'),fill=wes_palette("Zissou1")[c(1,3,5)],border=NA,cex=1.1,xpd=T,bty='n')
dev.off()

# 4. Examine case studies and output a csv summary ####

# Case study C552
sigList[["C552_D1_G7"]]
View(combList[["C552_C1_G3"]])
View(combList[["C552_D1_G5"]])

# Case study C559
sigList[["C559_A1_G10"]]
grep('ZC3H13',depthList[["C559_B1_G10"]]$SYMBOL)
head(depthList[["C559_B1_G10"]])

# Case study
View(sigList[["C516_B3_B1"]])
View(combList[["C516_B1_G7"]])
geneexp[which(geneexp$GeneID=='MMP9'),samples[grep('C516',samples)]]

# Make a table of selected mutations and patients
tabres <- as.data.frame(matrix(0L,nrow=14,ncol=10));row.names(tabres) <- samples[c(grep('C516',samples)[c(1,2)],grep('C552',samples),grep('C559',samples)[c(1:9)])]
colnames(tabres) <- c('Sample','Mutation','Gene','Conseq','NRD','NVD','NRR','NVR','FDR','Decision')
tabres[1,] <- c('C516_B1_G7',depthList[["C516_B1_G7"]][which(depthList[["C516_B1_G7"]]$SYMBOL=='MMP9'),c('Locus','SYMBOL','Consequence','NR','NV','NR_RNA','NV_RNA','fdr_fish')],'Non-ASE')
tabres[2,] <- c('C516_B3_B1',sigList[["C516_B3_B1"]][which(sigList[["C516_B3_B1"]]$SYMBOL=='MMP9'),c('Locus','SYMBOL','Consequence','NR','NV','NR_RNA','NV_RNA','fdr_fish')],'ASE')
tabres[3,] <- c('C552_C1_G3',sigList[["C552_C1_G3"]][which(sigList[["C552_C1_G3"]]$SYMBOL=='MUC5B')[2],c('Locus','SYMBOL','Consequence','NR','NV','NR_RNA','NV_RNA','fdr_fish')],'ASE')
tabres[4,] <- c('C552_D1_G5',combList[["C552_D1_G5"]][which(combList[["C552_D1_G5"]]$SYMBOL=='MUC5B')[2],c('Locus','SYMBOL','Consequence','NR','NV','NR_RNA','NV_RNA')],NA,'Low Cov')
tabres[5,] <- c('C552_D1_G7',sigList[["C552_D1_G7"]][which(sigList[["C552_D1_G7"]]$SYMBOL=='MUC5B'),c('Locus','SYMBOL','Consequence','NR','NV','NR_RNA','NV_RNA','fdr_fish')],'ASE')
tabres[6,] <- c('C559_A1_G10',sigList[["C559_A1_G10"]][which(sigList[["C559_A1_G10"]]$SYMBOL=='ZC3H13'),c('Locus','SYMBOL','Consequence','NR','NV','NR_RNA','NV_RNA','fdr_fish')],'ASE')
tabres[7,] <- c('C559_A1_G1',sigList[["C559_A1_G1"]][which(sigList[["C559_A1_G1"]]$SYMBOL=='ZC3H13'),c('Locus','SYMBOL','Consequence','NR','NV','NR_RNA','NV_RNA','fdr_fish')],'ASE')
tabres[8,] <- c('C559_B1_G10',depthList[["C559_B1_G10"]][which(depthList[["C559_B1_G10"]]$SYMBOL=='ZC3H13'),c('Locus','SYMBOL','Consequence','NR','NV','NR_RNA','NV_RNA','fdr_fish')],'Non-ASE')
tabres[9,] <- c('C559_B1_G7',sigList[["C559_B1_G7"]][which(sigList[["C559_B1_G7"]]$SYMBOL=='ZC3H13'),c('Locus','SYMBOL','Consequence','NR','NV','NR_RNA','NV_RNA','fdr_fish')],'ASE')
tabres[10,] <- c('C559_B1_G8',sigList[["C559_B1_G8"]][which(sigList[["C559_B1_G8"]]$SYMBOL=='ZC3H13'),c('Locus','SYMBOL','Consequence','NR','NV','NR_RNA','NV_RNA','fdr_fish')],'ASE')
tabres[11,] <- c('C559_C1_G4',sigList[["C559_C1_G4"]][which(sigList[["C559_C1_G4"]]$SYMBOL=='ZC3H13'),c('Locus','SYMBOL','Consequence','NR','NV','NR_RNA','NV_RNA','fdr_fish')],'ASE')
tabres[12,] <- c('C559_C1_G7',sigList[["C559_C1_G7"]][which(sigList[["C559_C1_G7"]]$SYMBOL=='ZC3H13'),c('Locus','SYMBOL','Consequence','NR','NV','NR_RNA','NV_RNA','fdr_fish')],'ASE')
tabres[13,] <- c('C559_C1_G8',depthList[["C559_C1_G8"]][which(depthList[["C559_C1_G8"]]$SYMBOL=='ZC3H13'),c('Locus','SYMBOL','Consequence','NR','NV','NR_RNA','NV_RNA','fdr_fish')],'Non-ASE')
tabres[14,] <- c('C559_D1_G10',combList[["C559_D1_G10"]][which(combList[["C559_D1_G10"]]$SYMBOL=='ZC3H13'),c('Locus','SYMBOL','Consequence','NR','NV','NR_RNA','NV_RNA')],NA,'Low Cov')
tabres$Sample <- gsub('(C\\d+)_(\\S\\d+)_(\\S\\d+)','\\1 \\2 \\3',tabres$Sample)
tabres$Mutation <- gsub('(\\S+:\\S+)_(\\S+)','\\1:\\2',tabres$Mutation)
tabres$Conseq <- gsub('(\\S+)_(\\S+)','\\1 \\2',tabres$Conseq)

tabres[which(!is.na(tabres$FDR)),'FDR'] <- signif(tabres[which(!is.na(tabres$FDR)),'FDR'],digits = 2)
tabres[which(tabres$Conseq=='missense variant'),'Conseq'] <- 'missense'
tabres[which(tabres$Conseq=='stop gained'),'Conseq'] <- 'stop gain'

write.csv(tabres,file='~/Documents/Thesis/mainText/ase_of_snv_results.csv',row.names=F,quote=F)




