# Script to assess the heterogeneity of SNVs found in both DNA and RNA
# To do the following:
# 1) Load expression, SNV and sample data
# 2) Analyse the heterogeneity of SNVs within tumours DNA vs. RNA
# 3) Plot the heterogeneity of SNVs

# Main code ####

# 0. Prepare environment
setwd("~/Documents/EPICC/Data/Mutations/VCFtoTSV")
'%ni%' <- Negate('%in%');library(dplyr);library(data.table)

# 1. Load SNV and gene expression data ####
annosnv <- readRDS('~/Documents/EPICC/Data/Mutations/VCFtoTSV/allsnvs_annotated.rds')
wgssam <- readRDS('~/Documents/EPICC/Data/QTLs/nodesamples/wgs_used_in_phylo.rds')
wgssam <- gsub('EPICC_(C\\d+_\\S\\d+_\\S\\d+)\\S+','\\1',wgssam)
rnasam <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]
samples <- wgssam[which(wgssam %in% rnasam)]

allpats <- gsub('(C\\d+)\\S+','\\1',samples);patients <- unique(allpats)
patients <- names(which(table(allpats)>=5))
samples <- samples[which(allpats %in% patients)]

vsd <- readRDS('~/Documents/EPICC/Data/Expression/ProcessedCounts/filgenes.vsd.ensembl.rds')
geneexp <- as.data.frame(assay(vsd))

# 2. Analyse the heterogeneity of SNVs within tumours in DNA and RNA data ####
snvnum <- c()
for(p in c(1:length(patients))) {
  pat <- patients[p];patsam <- samples[grep(pat,samples)]
  snvmat <- cbind(as.matrix(annosnv),matrix(0L,nrow=nrow(annosnv),ncol=length(patsam)))
  colnames(snvmat) <- c(colnames(annosnv),patsam);row.names(snvmat) <- annosnv$Locus
  
  # Get clonality distribution for DNA
  for(i in c(1:length(patsam))) {
    cursam <- patsam[i]
    print(paste0('Adding mutations from ',cursam,' - ',i,'/',length(patsam),' (',signif(i/length(patsam)*100,digits = 2),'%)'))
    
    if(file.exists(paste0('~/Documents/EPICC/Data/Mutations/VCFtoTSV/only_snvs/EPICC_',cursam,'_D1.txt'))) {
      tmpsnv <- as.data.frame(fread(paste0('~/Documents/EPICC/Data/Mutations/VCFtoTSV/only_snvs/EPICC_',cursam,'_D1.txt'),sep='\t',header=T))
      tmpsnv$Locus <- paste0(tmpsnv$Chr,':',tmpsnv$Pos,'_',tmpsnv$Ref,'/',tmpsnv$Alt)
      snvmat[row.names(snvmat) %in% tmpsnv$Locus,cursam] <- 1
    } else {
      print(paste0('SNV file for ',cursam,' is missing'))
    }
  }
  
  dnadf <- as.data.frame(snvmat)
  for(s in which(colnames(dnadf) %in% patsam)) { dnadf[,s] <- as.numeric(dnadf[,s]) }
  
  # Get clonality distribution for RNA
  snvmat <- cbind(as.matrix(annosnv),matrix(0L,nrow=nrow(annosnv),ncol=length(patsam)))
  colnames(snvmat) <- c(colnames(annosnv),patsam);row.names(snvmat) <- annosnv$Locus
  
  for(i in c(1:length(patsam))) {
    cursam <- patsam[i]
    print(paste0('Adding mutations from ',cursam,' - ',i,'/',length(patsam),' (',signif(i/length(patsam)*100,digits = 2),'%)'))
    tmprna <- fread(paste0('~/Documents/EPICC/Data/Mutations/VCFtoTSV/rna_matched/',cursam,'.txt'))
    tmprna <- tmprna[which(tmprna$VAF_RNA>0),]
    tmprna$Locus <- paste0(tmprna$Chr,':',tmprna$Pos,'_',tmprna$Ref,'/',tmprna$Alt)
    snvmat[row.names(snvmat) %in% tmprna$Locus,cursam] <- 1
  }
  
  rnadf <- as.data.frame(snvmat)
  for(s in which(colnames(rnadf) %in% patsam)) { rnadf[,s] <- as.numeric(rnadf[,s]) }
  
  rnadf <- rnadf[which(rowSums(rnadf[,patsam])>0),]
  #rnadf <- rnadf[which(rnadf$Gene %in% row.names(geneexp)),] # Barely changes result
  abovevaf <- as.numeric(apply(rnadf[,patsam],1, function(x) sum(x==1)))
  rnadf$Clonality <- 'notintumour';rnadf[abovevaf==1,'Clonality'] <- 'private'
  rnadf[abovevaf>1 & abovevaf<length(patsam),'Clonality'] <- 'subclonal'
  rnadf[abovevaf==length(patsam),'Clonality'] <- 'clonal'
  
  dnadf <- dnadf[which(rowSums(dnadf[,patsam])>0),]
  dnadf <- dnadf[which(dnadf$Locus %in% rnadf$Locus),]
  #dnadf <- dnadf[which(dnadf$Gene %in% row.names(geneexp)),] # Barely changes result
  abovevaf <- as.numeric(apply(dnadf[,patsam],1, function(x) sum(x==1)))
  dnadf$Clonality <- 'notintumour';dnadf[abovevaf==1,'Clonality'] <- 'private'
  dnadf[abovevaf>1 & abovevaf<length(patsam),'Clonality'] <- 'subclonal'
  dnadf[abovevaf==length(patsam),'Clonality'] <- 'clonal'
  
  snvnum <- c(snvnum,nrow(rnadf))
  
  clonal <- matrix(0L,nrow=3,ncol=2);row.names(clonal) <- c('private','subclonal','clonal');colnames(clonal) <- c('DNA','RNA')
  clonal['private','DNA'] <- length(which(dnadf$Clonality=='private'));clonal['subclonal','DNA'] <- length(which(dnadf$Clonality=='subclonal'));clonal['clonal','DNA'] <- length(which(dnadf$Clonality=='clonal'))
  clonal['private','RNA'] <- length(which(rnadf$Clonality=='private'));clonal['subclonal','RNA'] <- length(which(rnadf$Clonality=='subclonal'));clonal['clonal','RNA'] <- length(which(rnadf$Clonality=='clonal'))
  assign(paste0(pat,'_clonal'),clonal)
}

# 3. Plot the heterogeneity of SNVs in DNA and RNA ####
pdf('~/Documents/Thesis/figures/RNAvsDNAclonality.pdf',width=8,height=8)
layout(matrix(c(1:6),nrow=2,ncol=3,byrow = T));par(font=2)
for(p in c(1:length(patients))) {
  pat <- patients[p];patsam <- samples[grep(pat,samples)]
  clonal <- get(paste0(pat,'_clonal'));colnames(clonal) <- c('','')
  clonal[,1] <- clonal[,1]/snvnum[p]*100;clonal[,2] <- clonal[,2]/snvnum[p]*100
  
  barplot(clonal,col=RColorBrewer::brewer.pal(3,'Oranges'),font=2,las=1,border=NA,cex.axis = 1.5)
  mtext(side=1,at=c(.7,1.85),text=c('DNA','RNA'),font=2,line=1,cex=1.5)
  mtext(side=1,at=c(1.275),text=paste0('n=',snvnum[p]),font=2,line=2.5)
  mtext(side=2,text='% of SNVs',font=2,line=2.5,cex=1.2)
  mtext(side=3,text=paste0(pat,' (',length(patsam),' samples)'),line=1,cex=1.5)
}
plot.new()
legend(x=-.1,y=0.7,legend=c('All','Some','One'),bty='n',fill=rev(RColorBrewer::brewer.pal(3,'Oranges')),xpd=F,border=NA,cex=2.2)
dev.off()
