# Script for looking at the expression (in RNA-seq bams) of WGS-based SNVs
# To do the following:
# 1) Load data
# 2) Combine WGS SNVs with RNA-seq read count data
# 3) Analyse results of above
# 4) Plot results and output summary csv
# 5) Compare results to the depth of RNA-seq samples

# Main code ####

# 0. Prepare environment
setwd("~/Documents/EPICC/Data/Mutations/VCFtoTSV")
'%ni%' <- Negate('%in%');nucs <- c('A','C','G','T')
library(dplyr);library(data.table);library(stringr)

# 1. Get sample and expression data ####
wgssam <- readRDS('~/Documents/EPICC/Data/QTLs/nodesamples/wgs_used_in_phylo.rds')
wgssam <- gsub('EPICC_(C\\d+_\\S\\d+_\\S\\d+)\\S+','\\1',wgssam)
rnasam <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]
samples <- wgssam[which(wgssam %in% rnasam)]

geneexp <- as.data.frame(fread('~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_tpm.txt',stringsAsFactors = F))

# 2. Load WGS SNVs and combine with RNA-seq bam read count data ####
sumMuts <- matrix(0L,nrow=length(samples),ncol=5);row.names(sumMuts) <- samples
rerun <- c();useable <- c()
colnames(sumMuts) <- c('percNonExome','percExome','percNonExpExome','percExpExome','percHighExome')
for(s in c(1:length(samples))) { 
  sam <- samples[s]
  
  nonexpgene <-  geneexp[which(geneexp[,sam]<1),'GeneID']
  expgene <- geneexp[which(geneexp[,sam]>=1 & geneexp[,sam]<10),'GeneID']
  expgenehigh <- geneexp[which(geneexp[,sam]>=10),'GeneID']
  
  wgsmut <- as.data.frame(fread(paste0('~/Documents/EPICC/Data/Mutations/VCFtoTSV/only_snvs/EPICC_',sam,'_D1.txt'),sep='\t',header=T))
  if(file.exists(paste0('~/Documents/COLCC/EPICC/epicc_rna_data/other/tools_tried/bamreadcount/recal/results/',sam,'.txt'))) {
    print(paste0('Looking up RNA bam-readcount results for sample ',sam,' (',s,'/',length(samples),')'))
    rnacounts <- read.table(paste0('~/Documents/COLCC/EPICC/epicc_rna_data/other/tools_tried/bamreadcount/recal/results/',sam,'.txt'),fill=T)
    if(nrow(rnacounts)>1) {
      colnames(rnacounts) <- c('Chr','Pos','Ref','Depth','None','Ainfo','Cinfo','Ginfo','Tinfo','N')
      rnacounts <- rnacounts[which(rnacounts$Depth>0),];row.names(rnacounts) <- c(1:nrow(rnacounts))
      rnacounts$A <- as.numeric(sapply(str_split(rnacounts$Ainfo,':'), "[[", 2));rnacounts$C <- as.numeric(sapply(str_split(rnacounts$Cinfo,':'), "[[", 2))
      rnacounts$G <- as.numeric(sapply(str_split(rnacounts$Ginfo,':'), "[[", 2));rnacounts$T <- as.numeric(sapply(str_split(rnacounts$Tinfo,':'), "[[", 2))
        
      nrefs <- c();nvs <- c()
      for(i in c(1:nrow(rnacounts))) {
        curalt <- wgsmut[which(wgsmut$Chr==rnacounts[i,'Chr'] & wgsmut$Pos==rnacounts[i,'Pos'] & wgsmut$Ref==rnacounts[i,'Ref']),'Alt']
        nrefs <- c(nrefs,rnacounts[i,rnacounts[i,'Ref']])
        nvs <- c(nvs,rnacounts[i,curalt])
      }
      rnacounts$NR_RNA <- nvs+nrefs;rnacounts$NV_RNA <- nvs
      rnacounts <- rnacounts[which(rnacounts$NR_RNA>0),];row.names(rnacounts) <- c(1:nrow(rnacounts))
      rnacounts$VAF_RNA <- rnacounts$NV/rnacounts$NR
      rnacounts <- rnacounts[,c('Chr','Pos','Ref','NR_RNA','NV_RNA','VAF_RNA')]
        
      combmut <- merge(wgsmut,rnacounts,by.x=c('Chr','Pos','Ref'))
      assign(paste0(sam,'_muts'),combmut)
      
      write.table(combmut,file=paste0('~/Documents/EPICC/Data/Mutations/VCFtoTSV/rna_matched_recal/',sam,'.txt'),row.names = F,quote=F,sep='\t')
      filmut <- combmut[which(combmut$VAF_RNA>0),]
      
      sumMuts[s,] <- c(length(which(filmut$EXON==''))/length(which(wgsmut$EXON==''))*100,
                       length(which(filmut$EXON!=''))/length(which(wgsmut$EXON!=''))*100,
                       length(which(filmut$EXON!='' & filmut$Gene %in% nonexpgene))/length(which(wgsmut$EXON!='' & wgsmut$Gene %in% nonexpgene))*100,
                       length(which(filmut$EXON!='' & filmut$Gene %in% expgene))/length(which(wgsmut$EXON!='' & wgsmut$Gene %in% expgene))*100,
                       length(which(filmut$EXON!='' & filmut$Gene %in% expgenehigh))/length(which(wgsmut$EXON!='' & wgsmut$Gene %in% expgenehigh))*100)
      
    } else {
      print(paste0('No mutations for sample ',sam,' - rerun'))
      rerun <- c(rerun,s)
    }
  } else {
    print(paste0("RNA bam-readcount file doesn't exist for ",sam," (",s,"/",length(samples),")"))
    rerun <- c(rerun,s)
  }
}

saveRDS(sumMuts,'~/Documents/ThesisOther/ScriptsForThesis/RData/expressed_mutations.rds')

# 3. Analyse SNV expression results ####
sumMuts <- readRDS('~/Documents/ThesisOther/ScriptsForThesis/RData/expressed_mutations.rds')

sumtmp <- as.data.frame(sumMuts)

toplot <- sumtmp[,c(3:5)];colnames(toplot) <- c('Non/Low Expressed','Mid Expressed','High Expressed')

library(FSA)
aovdf <- data.frame(Cat=rep(colnames(sumtmp[,c(3:5)]),each=nrow(sumtmp)),Perc=0)
aovdf[aovdf$Cat=='percNonExpExome','Perc'] <- sumtmp$percNonExpExome
aovdf[aovdf$Cat=='percExpExome','Perc'] <- sumtmp$percExpExome
aovdf[aovdf$Cat=='percHighExome','Perc'] <- sumtmp$percHighExome
kw <- kruskal.test(Perc ~ Cat,data = aovdf)
dt = dunn.test::dunn.test(aovdf$Perc,aovdf$Cat,method="bh",altp=T)

# 4. Plot summary figure and output data of SNV expression ####
pdf('~/Documents/Thesis/figures/SNVsDNAtoRNA.pdf',width=11,height=8)
library(vioplot);par(mar=c(4,4,4,0))
vioplot(toplot,font.axis=2,cex.axis=1.2,col=c('firebrick4','firebrick3','firebrick2','firebrick1'),ylim=c(0,70),border = c('firebrick3','firebrick2','firebrick1','pink'),axes=F,las=1)
mtext(side=2,text='% of WGS SNVs found in matched RNA',font=2,line=2.4,cex=1.2)
segments(1,45,x1=2,lwd=2);segments(c(1,2),45,y1=43,lwd=2);text(x=1.5,y=47,labels = "**",cex=1.5,xpd=T)
segments(1,69,x1=3,lwd=2);segments(c(1,3),69,y1=67,lwd=2);text(x=2,y=71,labels = "***",cex=1.5,xpd=T)
segments(2,65,x1=3,lwd=2);segments(c(2,3),65,y1=63,lwd=2);text(x=2.5,y=67,labels = "***",cex=1.5,xpd=T)
mtext(side=3,at=1.85,text=bquote('Kruskal-Wallis Test '*chi^2*'('*.(kw$parameter)*')='*.(signif(kw$statistic,digits=5))*', p='*.(formatC(kw$p.value, format = "e", digits = 2))*',n='*.(nrow(sumtmp))*'         Post-hoc test: Dunn (FDR)'),cex=1.1,font=2,line=1)
dev.off()

tocsv <- sumtmp[,c(2:5)];tocsv$Sample <- row.names(tocsv)
row.names(tocsv) <- c(1:nrow(tocsv));tocsv <- tocsv[,c(5,1:4)]
colnames(tocsv) <- c('Samples','AllExonic','Non/Low Expressed','Mid Expressed','High Expressed')
for(i in c(2:ncol(tocsv))) {
  tocsv[,i] <- signif(tocsv[,i],4)
}
tocsv$Samples <- gsub('(C\\d+)_(\\S\\d+)_(\\S\\d+)','\\1 \\2 \\3',tocsv$Samples)

write.csv(tocsv,file='~/Documents/Thesis/mainText/supplementary_expressedmuts.csv',quote=F,row.names = F)

# 5. Compare to SNV expression results to RNA depth ####
sumMuts <- as.data.frame(readRDS('~/Documents/ThesisOther/ScriptsForThesis/RData/expressed_mutations.recal.rds'))
qc <- read.csv('~/Documents/ThesisOther/ScriptsForThesis/Tables/Final_EPICC_QC_ReadCounts.csv')
qc <- qc[which(qc$QC=='PASS'),];row.names(qc) <- c(1:nrow(qc))

covmuts <- merge(qc[,c('Sample','UsableReads')],sumMuts,by.x='Sample',by.y=0)

pdf('~/Documents/Thesis/figures/correlation_depth_vs_matchedSNVs.pdf')
par(font=2,font.axis=2,font.lab=2,mar=c(5,5,1,1))
plot(covmuts$UsableReads,covmuts$percExome,pch=16,las=1,cex.lab=1.2,cex.axis=1.2,cex=1.2,col=scales::alpha('firebrick3',0.5),log='x',xlab='RNA-seq Depth (million reads)',ylab='% WGS SNVs found in matched RNA')
abline(lm(covmuts$percExome~log(covmuts$UsableReads,base=10)),lwd=2,lty=2)
legend('topright',legend=c(paste0('r = ',signif(cor(log(covmuts$UsableReads,base=10),covmuts$percExome),digits=3)),
                           paste0('p = ',signif(cor.test(log(covmuts$UsableReads,base=10),covmuts$percExome)$p.value,digits=3))),bty='n')
dev.off()
