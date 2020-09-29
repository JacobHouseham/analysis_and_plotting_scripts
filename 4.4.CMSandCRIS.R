# Script to analyse the heterogeneity of CMS and CRIS classifications
# To do the following:
# 1) Load data and define CMS and CRIS classes
# 2) Determine low and high confidence results
# 3) Reformat results to by-patient
# 4) Plot CMS and CRIS results
# 5) Compare heterogeneity of both classifications
# 6) Analyse overlap of CMS/CRIS and quadrant genes
# 7) Plot enrichment of CMS/CRIS in quadrant subset genes

# Main code ####

# 0. Prepare environment
setwd("~/Documents/EPICC/Data/Expression/")
library(data.table);library(CMScaller)
colcris <- c(A='#F7A143',B='#EE1F25',C='#3753A4',D='#6BBC42',E='#71C69F',Unk='#D4D4D4')
colcms <- c(CMS1='#E4E515',CMS2='#8789C2',CMS3='#C191C3',CMS4='#33B577',Unk='#D4D4D4')
samples <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1];samples <- samples[which(!grepl('^.+E1.+',samples))]
allpats <- gsub('(C\\d+)_\\S+','\\1',samples);patients <- unique(allpats)

# 1. Load expression data and use CMScaller to call CMS and CRIS ####
counts <- as.data.frame(fread("~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_symbol_counts.txt"))
row.names(counts) <- counts$GeneID
entrezids <- mapIds(org.Hs.eg.db, keys=row.names(counts), column="ENTREZID", keytype="SYMBOL", multiVals="first")
inds <- which(!is.na(entrezids));found_ids <- entrezids[inds]
entrezcounts <- counts[names(found_ids), ];rownames(entrezcounts) <- found_ids
forcris <- as.matrix(entrezcounts[,samples])
cms <- CMScaller(forcris,templates = templates.CMS,rowNames = 'entrez',RNAseq=TRUE,nPerm=1000,FDR=0.05,doPlot = T,verbose = T)
cris <- CMScaller(forcris,templates = templates.CRIS,rowNames = 'entrez',RNAseq=TRUE,nPerm=1000,FDR=0.05,doPlot = T,verbose = T)

# 2. Determine low and high confidence calls ####
callingRes <- data.frame(Samples=samples,CMS='',Con_CMS='Low',CRIS='',Con_CRIS='Low')
for(i in c(1:nrow(callingRes))) {
  curcms <- cms[which(row.names(cms)==callingRes[i,'Samples']),]
  callingRes[i,'CMS'] <- gsub('d.(CMS\\d)','\\1',names(which.min(curcms[,grep('CMS',colnames(curcms))])))
  if(curcms$FDR<0.05) { callingRes[i,'Con_CMS'] <- 'High' }
  curcris <- cris[which(row.names(cris)==callingRes[i,'Samples']),]
  callingRes[i,'CRIS'] <- gsub('d.(CRIS\\S)','\\1',names(which.min(curcris[,grep('CRIS',colnames(curcris))])))
  if(curcris$FDR<0.05) { callingRes[i,'Con_CRIS'] <- 'High' }
}

# 3. Organise results by patient ####
patRes <- data.frame(CMS1=rep(0,length(patients)),CMS2=0L,CMS3=0L,CMS4=0L,CRISA=0L,CRISB=0L,CRISC=0L,CRISD=0L,CRISE=0L,NumSam=as.numeric(table(allpats)),row.names = patients)
patHigh <- patRes
for(i in c(1:length(patients))) {
  pat <- patients[i]
  patsam <- samples[grep(pat,samples)]
  patHigh[i,'NumSam'] <- length(which(callingRes$Samples %in% patsam & callingRes$Con_CMS=='High'))
  for(cat in paste0('CMS',c(1:4))) {
    patRes[i,cat] <- length(which(callingRes[which(callingRes$Samples %in% patsam),'CMS']==cat))
    patHigh[i,cat] <- length(which(callingRes[which(callingRes$Samples %in% patsam & callingRes$Con_CMS=='High'),'CMS']==cat))
  }
  for(cat2 in paste0('CRIS',c('A','B','C','D','E'))) {
    patRes[i,cat2] <- length(which(callingRes[which(callingRes$Samples %in% patsam),'CRIS']==cat2))
    patHigh[i,cat2] <- length(which(callingRes[which(callingRes$Samples %in% patsam & callingRes$Con_CRIS=='High'),'CRIS']==cat2))
  }
}

# 4. Plot all CMS and CRIS (low and high confidence) results ####
pdf('~/Documents/Thesis/figures/CMSandCRIS.pdf',width=14,height=8)
layout(matrix(c(1:4),nrow=2,ncol=2))
par(mar=c(4.1,4.1,2.6,1.1),font=2,font.axis=2)
cmsplot <- t(patRes[,grep('CMS',colnames(patRes))])
barplot(cmsplot,col=colcms,border=NA,ylim=c(0,40),las=2)
mtext(side=2,text='Number of samples',line=2.5,font=2,xpd=T)
text(x=4,y=43,labels='(a) CMS - All samples',font=2,xpd=T,cex=1.4)
legend(x=-1,y=42,bty='n',pt.cex=2,cex=1.2,fill=colcms[1:4],legend=c('CMS1','CMS2','CMS3','CMS4'),border=NA,xpd=T)

cmshighplot <- t(patHigh[,grep('CMS',colnames(patHigh))])
barplot(cmshighplot,col=colcms,border=NA,ylim=c(0,40),las=2)
mtext(side=2,text='Number of samples',line=2.5,font=2,xpd=T)
text(x=8.8,y=43,labels='(c) CMS - High confidence samples',font=2,xpd=T,cex=1.4)
legend(x=-1,y=42,bty='n',pt.cex=2,cex=1.2,fill=colcms[1:4],legend=c('CMS1','CMS2','CMS3','CMS4'),border=NA,xpd=T)

crisplot <- t(patRes[,grep('CRIS',colnames(patRes))])
barplot(crisplot,col=colcris,border=NA,ylim=c(0,40),las=2)
mtext(side=2,text='Number of samples',line=2.5,font=2,xpd=T)
text(x=4,y=43,labels='(b) CRIS - All samples',font=2,xpd=T,cex=1.4)
legend(x=-1,y=42,bty='n',pt.cex=2,cex=1.1,fill=colcris[1:5],legend=c('CRIS-A','CRIS-B','CRIS-C','CRIS-D','CRIS-E'),border=NA,xpd=T)

crishighplot <- t(patHigh[,grep('CRIS',colnames(patHigh))])
barplot(crishighplot,col=colcris,border=NA,ylim=c(0,40),las=2)
mtext(side=2,text='Number of samples',line=2.5,font=2,xpd=T)
text(x=8.8,y=43,labels='(d) CRIS - High confidence samples',font=2,xpd=T,cex=1.4)
legend(x=-1,y=42,bty='n',pt.cex=2,cex=1.1,fill=colcris[1:5],legend=c('CRIS-A','CRIS-B','CRIS-C','CRIS-D','CRIS-E'),border=NA,xpd=T)
dev.off()


# 5. Compare CRIS and CMS heterogenity ####

multpat <- row.names(patRes)[which(patRes$NumSam>=5)]
multpat_high <- names(patHigh)[which(patHigh$NumSam>=5)]

# Get diversity scores
maxCRIS <- log(5,base=exp(1));maxCMS <- log(4,base=exp(1))

evenCMS <- diversity(t(cmsplot[,multpat]))/maxCMS
evenCRIS <- diversity(t(crisplot[,multpat]))/maxCRIS

wilcox.test(evenCMS,evenCRIS,paired = T)

pdf('~/Documents/Thesis/figures/ShanEvennessCMSCRIS.pdf')
par(mar=c(4.6,5.1,3.6,0.6),font.lab=2,cex.axis=1.5,font.axis=2,font=2)
boxplot(evenCMS,evenCRIS,ylim=c(0,1),frame=F,las=1,names=c('CMS','CRIS'),ylab='',xlab='',outline=T,boxcol=c('skyblue3','firebrick3'),boxlwd=3,col=c('skyblue2','firebrick2'),medcol=c('skyblue4','firebrick4'),medlwd=3.5,staplecol=c('skyblue4','firebrick4'),staplelwd=3,border=c('skyblue4','firebrick4'),outcol=c('skyblue4','firebrick4'),outpch=20)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="#EBEBEB",border=NA)
grid(lty=1,col='white');par(new=T)
boxplot(evenCMS,evenCRIS,ylim=c(0,1),frame=F,las=1,names=c('CMS','CRIS'),ylab='',xlab='',outline=T,boxcol=c('skyblue3','firebrick3'),boxlwd=3,col=c('skyblue2','firebrick2'),medcol=c('skyblue4','firebrick4'),medlwd=3.5,staplecol=c('skyblue4','firebrick4'),staplelwd=3,border=c('skyblue4','firebrick4'),outcol=c('skyblue4','firebrick4'),outpch=20)
segments(1,1.1,x1=2,lwd=3,xpd=T);segments(c(1,2),1.1,y1=1.05,lwd=3,xpd=T);text(x=1.5,y=1.13,labels = "n.s.",cex=1.5,xpd=T)
mtext(side=2,"Shannon's Evenness (H/Hmax)",line=3,cex=1.5)
mtext(side=1,text='Classification Signatures',line=3,cex=1.5,font=2)
dev.off()

rescor <- cor.test(evenCMS,evenCRIS)

pdf('~/Documents/Thesis/figures/CMSCRIScorrelation.pdf')
par(mar=c(4.1,4.6,1.1,1.1),font.lab=2,cex.axis=1.5,font.axis=2,font=2)
plot(evenCMS,evenCRIS,xlim=c(0,1),ylim=c(0,1),cex=1.5,pch=16,col='firebrick3',axes=F,xlab='',ylab='')
axis(side=1,at=seq(0,1,by=0.2),labels=seq(0,1,by=0.2),font=2,lwd=2)
mtext(side=1,"Evenness of CMS classification",line=2.5,cex=1.3)
axis(side=2,at=seq(0,1,by=0.2),labels=seq(0,1,by=0.2),font=2,lwd=2,las=2)
mtext(side=2,"Evenness of CRIS classification",line=2.9,cex=1.3)
abline(lm(evenCRIS~evenCMS),lty=2,lwd=2,col='dimgray')
text(x=0.6,y=0.95,col='dimgray',labels=paste0('Pearson correlation = ',signif(rescor$estimate,3),' (p=',signif(rescor$p.value,2),')'),cex=1.3)
text(x=evenCMS,y=evenCRIS,labels=multpat,offset=0.5,pos=c(1,1,1,1,3,2,1,3,2,1,4,1,4,2,4,1,1,4),cex=1.3)
dev.off()

# 6. Look at enrichment of Quadrant genes in CMS and CRIS genes ####
quadcols <- c(Q1='#E43B38',Q2='#8E55A2',Q3='#F8D30E',Q4='#5FC5D5')
geneinfo <- fread('~/Documents/EPICC/Data/Expression/compiledGeneInfo.txt',data.table = F)
# Load data from 4.3.4
quad_ensembl <- readRDS('~/Documents/ThesisOther/ScriptsForThesis/RData/quadrant_genes_ensembl.rds')
quad_symbol <- list()
for(i in c(1:4)) {
  curq <- quad_ensembl[[paste0('q',i)]]
  quad_symbol[[paste0('q',i)]] <- geneinfo[which(geneinfo$GeneID %in% curq),'Name']
}
allq <- unlist(quad_symbol);rm(quad_ensembl)

percSig <- data.frame(CMS=rep(0,4),CRIS=rep(0,4),All=rep(0,4),stringsAsFactors=F)
row.names(percSig) <- c(paste0('Q',c(1:4)))
for(i in c(1:4)) {
  curq <- quad_symbol[[paste0('q',i)]]
  percSig[i,1] <- length(which(unique(templates.CMS$symbol) %in% curq))/length(unique(templates.CMS$symbol[which(templates.CMS$symbol %in% allq)]))*100
  percSig[i,2] <- length(which(unique(templates.CRIS$symbol) %in% curq))/length(unique(templates.CRIS$symbol[which(templates.CRIS$symbol %in% allq)]))*100
  percSig[i,3] <- length(curq)/length(allq)*100
}
perSig <- as.matrix(percSig[c(1:4),])
barplot(perSig,beside=T,col=quadcols[c(1:4)],border='white')

enrichCMS <- as.data.frame(matrix(0L,nrow=4,ncol=8));row.names(enrichCMS) <- paste0('Q',c(1:4))
colnames(enrichCMS) <- c('NotCMSNotQ','NotCMSQ','CMSNotQ','CMSQ','OR','lowOR','highOR','pval')
enrichCRIS <- enrichCMS;colnames(enrichCRIS) <- c('NotCRISNotQ','NotCRISQ','CRISNotQ','CRISQ','OR','lowOR','highOR','pval')
for(i in c(1:4)) {
  test <- i;nontest <- c(1:4)[which(c(1:4)!=i)]
  curquad <- quad_symbol[[paste0('q',i)]];curnonquad <- unique(c(quad_symbol[[paste0('q',nontest[1])]],quad_symbol[[paste0('q',nontest[2])]],quad_symbol[[paste0('q',nontest[3])]]))
  testmat <- rbind(c(length(which(curnonquad %ni% unique(templates.CMS$symbol))),length(which(curquad %ni% unique(templates.CMS$symbol)))),c(length(which(curnonquad %in% unique(templates.CMS$symbol))),length(which(curquad %in% unique(templates.CMS$symbol)))))
  enrichCMS[i,'NotCMSNotQ'] <-length(which(curnonquad %ni% unique(templates.CMS$symbol)));enrichCMS[i,'NotCMSQ'] <-length(which(curquad %ni% unique(templates.CMS$symbol)))
  enrichCMS[i,'CMSNotQ'] <-length(which(curnonquad %in% unique(templates.CMS$symbol)));enrichCMS[i,'CMSQ'] <-length(which(curquad %in% unique(templates.CMS$symbol)))
  res <- fisher.test(testmat)
  enrichCMS[i,'OR'] <- res$estimate;enrichCMS[i,'lowOR'] <- res$conf.int[1];enrichCMS[i,'highOR'] <- res$conf.int[2]
  enrichCMS[i,'pval'] <- res$p.value
  
  testmat <- rbind(c(length(which(curnonquad %ni% unique(templates.CRIS$symbol))),length(which(curquad %ni% unique(templates.CRIS$symbol)))),c(length(which(curnonquad %in% unique(templates.CRIS$symbol))),length(which(curquad %in% unique(templates.CRIS$symbol)))))
  enrichCRIS[i,'NotCRISNotQ'] <-length(which(curnonquad %ni% unique(templates.CRIS$symbol)));enrichCRIS[i,'NotCRISQ'] <-length(which(curquad %ni% unique(templates.CRIS$symbol)))
  enrichCRIS[i,'CRISNotQ'] <-length(which(curnonquad %in% unique(templates.CRIS$symbol)));enrichCRIS[i,'CRISQ'] <-length(which(curquad %in% unique(templates.CRIS$symbol)))
  res <- fisher.test(testmat,conf.level = 0.95)
  enrichCRIS[i,'OR'] <- res$estimate;enrichCRIS[i,'lowOR'] <- res$conf.int[1];enrichCRIS[i,'highOR'] <- res$conf.int[2]
  enrichCRIS[i,'pval'] <- res$p.value
}
sigs <- c();for(p in p.adjust(c(enrichCMS$pval,enrichCRIS$pval),method='fdr')) { sigs <- c(sigs,getstars(p)) }

# 7. Plot enrichment of CMS and CRIS genes in quadrants ####

pdf('~/Documents/Thesis/figures/QuadEnrichmentCMSCRIS.pdf')
par(mar=c(4.6,5.1,3.6,0.6),font.lab=2,cex.axis=1.5,font.axis=2,font=2)
xx <- barplot(rbind(enrichCMS$OR,enrichCRIS$OR),beside=T,axes=F,las=2,col=scales::alpha(c('skyblue3','firebrick3'),0.9),ylim=c(0,3),border=NA);abline(h=1,lty=3,lwd=2)
axis(side=2,at=c(0,1,2,3),labels=c(0,1,2,3),lwd=2,las=2)
mtext(side=2,line=3,text='Odds ratio: CMS/CRIS genes in quadrants',cex=1.5)
axis(side=1,at=(xx[2,]-((xx[2,]-xx[1,])/2)),labels = paste0('Q',c(1:4)),lwd=2,line=1)
segments(x0=xx[1,],y0=enrichCMS$lowOR,y1=enrichCMS$highOR,lwd=4)
segments(x0=xx[2,],y0=enrichCRIS$lowOR,y1=enrichCRIS$highOR,lwd=4)
text(x=c(xx[1,],xx[2,]),y=0.1,labels=sigs,cex=1.5,font=2)
legend(x=0.5,y=3.1,legend=c('CMS','CRIS'),fill=scales::alpha(c('skyblue3','firebrick3'),0.9),cex=2,bty='n',border = NA)
dev.off()
