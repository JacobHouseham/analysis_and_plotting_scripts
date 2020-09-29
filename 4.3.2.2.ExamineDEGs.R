# Script to analyse the results of DEG between regions analysis
# To do the following:
# 1) Load DEG data
# 2) Assess the bias of sample number of DEG numbers and compare to random analysis
# 3) Define heterogeneity groups based on DEG numbers
# 4) Plot the number of DEGs per patient
# 5) Do Reactome enrichment on highly heterogeneous and homogenous genes

# Main code ####

# 0. Prepare environment
setwd('~/Documents/EPICC/Data/Expression/GSEA')
hg38genes <- read.table("~/Documents/EPICC/Data/Expression/compiledGeneInfo.txt",header=T)
library(dplyr);library(RColorBrewer)

# 1. Load DEG data ####
degDiff <- readRDS('~/Documents/EPICC/Data/Expression/GSEA/degDiff.rds')

numcomp <- table(gsub('(C\\d+)_\\S+','\\1',colnames(degDiff)))
multcomp <- names(which(numcomp>1))
samples <- colnames(degDiff)[which(gsub('(C\\d+)_\\S+','\\1',colnames(degDiff)) %in% multcomp)]

degDiff <- degDiff[,samples]
commonDEGs <- arrange(data.frame(GeneID=row.names(degDiff),Common=rowSums(degDiff)),desc(Common))
write.table(commonDEGs[which(commonDEGs$Common>=10),'Name'],file='commonDEGs.txt',sep='\n',row.names=F,col.names=F,quote=F)

randDiff <- readRDS('~/Documents/EPICC/Data/Expression/GSEA/randDiff.rds')
randDiff <- randDiff[,samples]
randDEGs <- arrange(data.frame(GeneID=row.names(randDiff),Common=rowSums(randDiff)),desc(Common))

# 2. Assess whether the number of samples in comparisons correlates with number of DEGs ####
plot(colSums(randDiff),colSums(degDiff),xlim=c(0,5e3),ylim=c(0,5e3),pch=16,col='darkblue')
text(x=colSums(randDiff),y=colSums(degDiff),colnames(degDiff),cex=0.4,pos=4)
boxplot(data.frame(Random=colSums(randDiff),Comparisons=colSums(degDiff)))

rep1 <- as.numeric(gsub("C\\d+_\\S(\\d+)vs\\S\\d+","\\1",colnames(degDiff)))
rep2 <- as.numeric(gsub("C\\d\\d\\d_\\S\\d+vs\\S(\\d+)","\\1",colnames(degDiff)))
reg1 <- gsub("C\\d\\d\\d_(\\S)\\d+vs\\S\\d+","\\1",colnames(degDiff));reg2 <- gsub("C\\d\\d\\d_\\S\\d+vs(\\S)\\d+","\\1",colnames(degDiff))
reg1[which(reg1=='A')] <- rep(1,length(which(reg1=='A')));reg1[which(reg1=='B')] <- rep(2,length(which(reg1=='B')));reg1[which(reg1=='C')] <- rep(3,length(which(reg1=='C')));reg1[which(reg1=='D')] <- rep(4,length(which(reg1=='D')))
reg2[which(reg2=='A')] <- rep(1,length(which(reg2=='A')));reg2[which(reg2=='B')] <- rep(2,length(which(reg2=='B')));reg2[which(reg2=='C')] <- rep(3,length(which(reg2=='C')));reg2[which(reg2=='D')] <- rep(4,length(which(reg2=='D')))
reg1 <- as.numeric(reg1);reg2 <- as.numeric(reg2);regDiff <- reg2-reg1

tot <- rep1+rep2;lowest <- pmin(rep1,rep2);highest <- pmax(rep1,rep2);diffNum <- abs(rep2-rep1)
ndeg <- colSums(degDiff)
allDEGs <- data.frame(DEGs=ndeg,NumSam1=rep1,NumSam2=rep2,Tot=tot,Lowest=lowest,Highest=highest,Diff=diffNum)
par(mfrow=c(2,3))
for(i in c(2:ncol(allDEGs))) {
  res <- cor.test(ndeg,allDEGs[,i],method='spearman',exact=F)
  plot(ndeg,allDEGs[,i],pch=16,cex=0.5,xlab='Num DEGs',ylab=colnames(allDEGs)[i],main='',sub=paste0('Rho: ',signif(res$estimate,digits=3),' p-val: ',signif(res$p.value,digits=3)),col='darkblue')
  abline(lm(allDEGs[,i]~ndeg),col='darkred',lty=2)
}

forano <- data.frame(DEGs=ndeg,regDiff=as.factor(regDiff))
res.aov <- aov(ndeg ~ regDiff, data = forano)
summary(res.aov);TukeyHSD(res.aov)

# 3. Look at DEG numbers per patient and perform k-means clustering ####
allpats <- gsub('^(C\\d+)_.+','\\1',colnames(degDiff))
patients <- unique(allpats)
medians <- c()
for(pat in patients) {
  curpat <- grep(pat,colnames(degDiff))
  if(length(curpat)>1) {
    medians <- c(medians,median(colSums(degDiff[,curpat])))
  } else {
    medians <- c(medians,sum(degDiff[,curpat]))
  }
}
allDEGs$Patient <- allpats
clust <- kmeans(medians,2)


if(median(medians[which(clust$cluster==1)])<median(medians[which(clust$cluster==2)])) {
  groups <- rep(1,length(clust$cluster))
  groups[clust$cluster==2] <- 2
} else {
  groups <- rep(2,length(clust$cluster))
  groups[clust$cluster==1] <- 1
}
allDEGs$Heterogeneity <- rep(1,nrow(allDEGs));allDEGs[allDEGs$Patient %in% patients[groups==2],'Heterogeneity'] <- 2

# 4. Plot DEGs per patient with heterogeneity groups highlighted ####
colbox <- c('skyblue2','firebrick2')[groups]
colbord <- c('skyblue4','firebrick4')[groups]

pdf('~/Documents/Thesis/figures/DEGsPerPatient.pdf',height=7,width=12)
par(mar=c(5.1,7.1,2.1,1.1),font.lab=2,font.axis=2,font=2,cex.axis=1.25)
boxplot(DEGs~Patient,data=allDEGs,ylab='',ylim=c(0,3500),las=1,col=colbox,bty='n',frame.plot=F,lwd=1.5,pch=20,border=colbord,xlab='')
mtext(side=c(1,2),text=c('Patients','Number of DEGs between two regions'),line=c(3,4),cex=c(1.7,1.7))
legend(x=.5,y=3500,legend=c('Low ITH Group','High ITH Group'),fill=c('skyblue2','firebrick2'),cex=1.7,bty='n')
dev.off()

# 5. Perform Reactome pathway enrichment on het and hom genes ####
numcomp <- table(gsub('(C\\d+)_\\S+','\\1',colnames(degDiff)))
multcomp <- names(which(numcomp>1))
samples <- colnames(degDiff)[which(gsub('(C\\d+)_\\S+','\\1',colnames(degDiff)) %in% multcomp)]

degDiff <- degDiff[,samples]
commonDEGs <- arrange(data.frame(GeneID=row.names(degDiff),Common=rowSums(degDiff)),desc(Common))

hetdegs <- commonDEGs$GeneID[which(commonDEGs$Common>=15)]
nondegs <- commonDEGs$GeneID[which(commonDEGs$Common==0)]

library(ReactomePA);library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
symbolandentrez <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "entrezgene_id", values = commonDEGs$GeneID, bmHeader = T, mart = mart)
colnames(symbolandentrez) <- c('Name','Entrez')

enrichHet <- enrichPathway(gene=hetdegs,pvalueCutoff=0.01, readable=T)
enrichNon <- enrichPathway(gene=nondegs,pvalueCutoff=0.01, readable=T)

barplot(enrichHet)
barplot(enrichNon)



