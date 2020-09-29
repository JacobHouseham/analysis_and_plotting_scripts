# Script to perform quadrant gene expression heterogeneity analysis
# To do the following:
# 1) Load expression and sample data
# 2) Determine gene-ITH scores
# 3) Determine gene inter-tumour heterogeneity scores
# 4) Save those scores
# 5) Determine quadrants and plot
# 6) Plot comparison of quadrants to original paper
# 7) Do Reactome pathway analysis
# 8) Plot Reactome results
# 9) Compare CMS and CRIS genes to quadrants

# Functions ####
CV <- function(x){
  (sd(x)/mean(x))*100
}

# Main code ####

# 0. Prepare environment
library(data.table);library(DESeq2)
quadcols <- c(Q1='#E43B38',Q2='#8E55A2',Q3='#F8D30E',Q4='#5FC5D5')

# 1. Load normalised gene expression and sample info ####
vsd <- readRDS('~/Documents/EPICC/Data/Expression/ProcessedCounts/filgenes.vsd.ensembl.rds')
geneexp <- as.data.frame(assay(vsd));samples <- colnames(geneexp)
geneexp$GeneID <- row.names(assay(vsd));geneexp <- geneexp[,c('GeneID',samples)]

tsam <- samples[which(!grepl('E1',samples))];nsam <- samples[which(grepl('E1',samples))]
patnum <- table(gsub('(C\\d+)\\S+','\\1',tsam))
tsam <- tsam[which(gsub('(C\\d+)\\S+','\\1',tsam) %in% names(patnum[which(patnum>2)]))]

# 2. Get RNA-ITH scores ####
patients <- unique(gsub('(C\\d+)\\S+','\\1',tsam))
ithSD <- matrix(0L,ncol=length(patients),nrow=nrow(geneexp));colnames(ithSD) <- patients;row.names(ithSD) <- geneexp$GeneID
ithMAD <- ithCV <- ithSD
for(i in c(1:length(patients))) {
  pat <- patients[i];patsam <- tsam[grepl(pat,tsam)]
  metrics <- transform(geneexp[,patsam], SD=apply(geneexp[,patsam],1, sd, na.rm = TRUE),MAD=apply(geneexp[,patsam],1, mad, na.rm = TRUE),CV=apply(geneexp[,patsam],1, CV))
  
  ithSD[,pat] <- metrics$SD
  ithMAD[,pat] <- metrics$MAD
  ithCV[,pat] <- metrics$CV
}
numpat <- as.numeric(table(gsub('(C\\d+)\\S+','\\1',tsam)))
geneith <- apply(ithSD,1, median, na.rm = TRUE);names(geneith) <- row.names(ithSD)
geneMAD <- apply(ithMAD,1, median, na.rm = TRUE);names(geneMAD) <- row.names(ithSD)
geneCV <- apply(ithCV,1, median, na.rm = TRUE);names(geneCV) <- row.names(ithSD)
patith <- apply(ithSD,2, median, na.rm = TRUE);names(patith) <- colnames(ithSD)

barplot(patith,col='darkred',border=NA,las=2,font=2)

# 3. Get inter-tumour heterogeneity scores ####
interdf <- matrix(0L,ncol=10,nrow=nrow(geneexp));colnames(interdf) <- paste0('iter',c(1:ncol(interdf)));row.names(interdf) <- geneexp$GeneID
for(i in c(1:ncol(interdf))) {
  cursams <- c()
  for(pat in patients) {  patsam <- tsam[grepl(pat,tsam)];cursams <- c(cursams,sample(patsam,1))  }
  sds <- transform(geneexp[,cursams], SD=apply(geneexp[,cursams],1, sd))$SD
  interdf[,i] <- sds
}
geneinter <- apply(interdf,1, median, na.rm = TRUE);names(geneinter) <- row.names(interdf)

# 4. Save (and reload) ITH and inter heterogeneity scores ####
# Save geneith and geneinter so don't have to rerun everything every time
saveRDS(geneith,file='~/Documents/ThesisOther/ScriptsForThesis/RData/quadrant_intraSD.rds')
saveRDS(geneinter,file='~/Documents/ThesisOther/ScriptsForThesis/RData/quadrant_interSD.rds')

# Load those saved files
geneith <- readRDS(file='~/Documents/ThesisOther/ScriptsForThesis/RData/quadrant_intraSD.rds')
geneinter <- readRDS(file='~/Documents/ThesisOther/ScriptsForThesis/RData/quadrant_interSD.rds')

# 5. Extract quadrants and plot ####
# New get correlation of inter and intra
cor.test(geneinter,geneith,method='pearson')

# Identify quadrants
q1 <- names(geneith)[which(geneith>=mean(geneith) & geneinter<mean(geneinter))]
q2 <- names(geneith)[which(geneith<mean(geneith) & geneinter<mean(geneinter))]
q3 <- names(geneith)[which(geneith>=mean(geneith) & geneinter>=mean(geneinter))]
q4 <- names(geneith)[which(geneith<mean(geneith) & geneinter>=mean(geneinter))]
plotext <- c();for(i in c(1:4)) { plotext <- c(plotext,paste0('Q',i,' (',length(get(paste0('q',i))),' genes)'))}

pdf('~/Documents/Thesis/figures/DhruvaQuadrants.pdf')
layout(matrix(c(3,3,3,0,1,1,1,2,1,1,1,2,1,1,1,2),ncol=4))
par(mar=c(0,0,2.1,2.1),font=2)
x <- geneinter;y <- geneith
colpoints <- rep(quadcols['Q1'],length(x));colpoints[x<mean(x) & y<mean(y)] <- quadcols['Q2']
colpoints[x>=mean(x) & y>=mean(y)] <- quadcols['Q3'];colpoints[x>=mean(x) & y<mean(y)] <- quadcols['Q4']
plot(x,y,xaxt="none",yaxt="none",xlab='',ylab='',pch=16,col=colpoints,cex=0.6)
abline(v=mean(x),col='#67A5CB',lty=2);abline(h=mean(y),col='#F48A42',lty=2)
text(x=c(1.15,1.15,max(x)-.5,max(x)-.5),y=c(max(y)-.03,0.03,max(y)-.03,0.03),cex=1.5,labels = plotext,col=quadcols)
d.x <- density(x);d.y <- density(y)
par(mar=c(3.1,0,0,2.1))
plot(d.x$x, 1-d.x$y, xlim=range(x),axes=F,type='l',col='#67A5CB',xlab='',font.lab=2,lwd=3)
mtext(side=1,text='Between-tumour diversity (median sd)',cex=1.2,line = 0.5)
polygon(d.x$x, 1-d.x$y,col=scales::alpha('#67A5CB',0.4),border=NA)
par(mar=c(0,3.1,2.1,0))
plot(d.y$y, d.y$x, ylim=range(y),axes=F ,xlim=rev(range(d.y$y)),ylab='',type='l',col='#F48A42',font.lab=2,lwd=3)
mtext(side=2,text='Within-tumour diversity (median sd)',cex=1.2,line = 0.5)
polygon(d.y$y, d.y$x ,col=scales::alpha('#F48A42',0.4),border=NA)
dev.off()

# 6. Plot comparison of EPICC and TRACERx quadrants ####
epicc <- c(length(q1),length(q2),length(q3),length(q4))
epiccpercs <- epicc/sum(epicc)*100
# Input quadrant numbers from Dhruva paper
dhruva <- c(798,4766,9642,1080);dhruvapercs <- dhruva/sum(dhruva)*100

pdf('~/Documents/Thesis/figures/quadrants_compare.pdf')
par(font=2,mar=c(3,4.5,1,1))
x <-barplot(cbind(epiccpercs,dhruvapercs),col=quadcols,border=NA,cex.names = 1.5,cex.axis = 1.5,names.arg = c('EPICC','TRACERx'),font=2,las=1)
mtext(side=2,text='Genes per heterogeneity quadrant (%)',line=2.75,cex=1.5)
text(x=0.7,y=c(epiccpercs[1]/2,
               (epiccpercs[1]+epiccpercs[2])-(epiccpercs[2]/2),
               (epiccpercs[1]+epiccpercs[2]+epiccpercs[3])-(epiccpercs[3]/2),
               (epiccpercs[1]+epiccpercs[2]+epiccpercs[3]+epiccpercs[4])-(epiccpercs[4]/2)),labels=paste0('Q',c(1:4),': ',round(epiccpercs),'%'),col='white',cex=1.35)
text(x=1.9,y=c(dhruvapercs[1]/2,
               (dhruvapercs[1]+dhruvapercs[2])-(dhruvapercs[2]/2),
               (dhruvapercs[1]+dhruvapercs[2]+dhruvapercs[3])-(dhruvapercs[3]/2),
               (dhruvapercs[1]+dhruvapercs[2]+dhruvapercs[3]+dhruvapercs[4])-(dhruvapercs[4]/2)),labels=paste0('Q',c(1:4),': ',round(dhruvapercs),'%'),col='white',cex=1.35)
dev.off()

# 7. Perform Reactome pathway analysis ####
library(ReactomePA);library(biomaRt)

# Convert ensembl IDs to entrez IDs for each quadrant gene set
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",mirror = "useast")
convertgen <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), filters = "ensembl_gene_id", values =geneexp$GeneID, bmHeader = T, mart = mart)
colnames(convertgen) <- c('Ensembl','Entrez')
q1entrez <- convertgen[which(convertgen$Ensembl %in% q1),'Entrez'];q1entrez <- unique(q1entrez[which(!is.na(q1entrez))])
q2entrez <- convertgen[which(convertgen$Ensembl %in% q2),'Entrez'];q2entrez <- unique(q2entrez[which(!is.na(q2entrez))])
q3entrez <- convertgen[which(convertgen$Ensembl %in% q3),'Entrez'];q3entrez <- unique(q3entrez[which(!is.na(q3entrez))])
q4entrez <- convertgen[which(convertgen$Ensembl %in% q4),'Entrez'];q4entrez <- unique(q4entrez[which(!is.na(q4entrez))])

# Save these results
quad_ensembl <- list(q1,q2,q3,q4);names(quad_ensembl) <- paste0('q',c(1:4))
saveRDS(quad_ensembl,file='~/Documents/ThesisOther/ScriptsForThesis/RData/quadrant_genes_ensembl.rds')
quad_entrez <- list(q1entrez,q2entrez,q3entrez,q4entrez);names(quad_entrez) <- paste0('q',c(1:4))
saveRDS(quad_entrez,file='~/Documents/ThesisOther/ScriptsForThesis/RData/quadrant_genes_entrez.rds')

# Load these results back in and perform Reactome pathway enrichment analysis
quad_entrez <- readRDS(file='~/Documents/ThesisOther/ScriptsForThesis/RData/quadrant_genes_entrez.rds')
enrichq1 <- enrichPathway(gene=quad_entrez$q1,pvalueCutoff=0.05, readable=T)
enrichq2 <- enrichPathway(gene=quad_entrez$q2,pvalueCutoff=0.05, readable=T)
enrichq3 <- enrichPathway(gene=quad_entrez$q3,pvalueCutoff=0.05, readable=T)
enrichq4 <- enrichPathway(gene=quad_entrez$q4,pvalueCutoff=0.05, readable=T)

# 8. Plot enrichment of Reactome pathways ####
pdf('~/Documents/Thesis/figures/quadrants_reactome.pdf',width=12,height=8)
maxxs <- c(0,150,200,60);bys <- c(10,25,25,10)
lets <- c('(a) Q1 pathways','(c) Q2 pathways','(b) Q3 pathways','(d) Q4 pathways')
xlet <- c(-0.6,-85,-110,-35)
layout(matrix(c(1:4),nrow=2,ncol=2))
par(mar=c(4.1,16.6,2.1,2.1))
for(i in c(1:4)) {
  if(i==1) {
    plot.new()
    text(x=-0.3,y=0.5,labels='No pathways significantly enriched',cex=1.5,xpd=T,font=2)
    text(x=xlet[i],y=1.08,labels=lets[i],cex=2,xpd=T,font=2)
  } else {
    curreact <- head(get(paste0('enrichq',i))@result,n=10)
    if(i==2) {
      curreact$Description[which(curreact$ID=='R-HSA-2029482')] <- 'Regulation of actin dynamics'
      curreact$Description[which(curreact$ID=='R-HSA-72203')] <- 'Processing of Pre-mRNA'
    } else if(i==3) {
      curreact$Description[which(curreact$ID=='R-HSA-3700989')] <- 'TP53 Transcriptional Regulation'
      curreact$Description[which(curreact$ID=='R-HSA-5693538')] <- 'Homology Directed Repair (HDR)'
      curreact$Description[which(curreact$ID=='R-HSA-5693567')] <- 'HDR through HRR or SSA'
    } else if(i==4) {
      curreact$Description[which(curreact$ID=='R-HSA-1428517')] <- 'TCA cycle and RET'
      curreact$Description[which(curreact$ID=='R-HSA-163200')] <- 'RET and ATP synthesis'
      curreact$Description[which(curreact$ID=='R-HSA-69656')] <- 'Cdk2 events at S phase'
      curreact$Description[which(curreact$ID=='R-HSA-8939236')] <- 'RUNX1 regulates HSCs'
      curreact$Description[which(curreact$ID=='R-HSA-8939902')] <- 'Regulation of RUNX2'
      curreact$Description[which(curreact$ID=='R-HSA-75815')] <- 'Degradation of Cyclin D'
      curreact$Description[which(curreact$ID=='R-HSA-187577')] <- 'Degradation of p27/p21'
    }
    yy <- barplot(rev(curreact$Count),col=quadcols[i],border=NA,horiz=T,axes=F,xlim=c(0,maxxs[i]))
    axis(side=1,at=seq(0,maxxs[i],by=bys[i]),labels=seq(0,maxxs[i],by=bys[i]),font=2,xpd=T,cex=3.5,lwd=2)
    mtext(side=1,text='Number of genes',line=2.5,font=2,xpd=T)
    mtext(side=2,at=yy,text=rev(curreact$Description),las=2,line=0.5,font=2,xpd=T)
    text(x=xlet[i],y=13,labels=lets[i],cex=2,xpd=T,font=2)
  }
}
dev.off()





