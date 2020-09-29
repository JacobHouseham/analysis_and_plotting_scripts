# Script to get expression distance and compare to genetic distance
# To do the following:
# 1) Load expression data and calculate expression distance between samples
# 2) Plot expression distance by pair-type per patient
# 3) Model expression distance
# 4) Combine genetic and expression distance
# 5) Plot comparison of genetic vs. expression distance
# 6) Model expression distance against genetic distance

# Main code ####

# 0. Prepare environment
setwd('~/Documents/EPICC/Data/Expression/Expression_Distance/')
library(data.table);'%ni%' <- Negate('%in%');library(RColorBrewer)
options(scipen = -1);library(DESeq2)

# 1. Calculate expression distance ####
vsd <- readRDS('~/Documents/EPICC/Data/Expression/ProcessedCounts/filgenes.vsd.ensembl.rds')
geneexp <- as.data.frame(assay(vsd));rnasam <- colnames(geneexp)
rnaallpat <- gsub('(C\\d+)_\\S+','\\1',rnasam);rnapat <- unique(rnaallpat)

# Get expression distance for all pairwise comparisons
disexp <- data.frame(matrix(ncol=6,nrow=0));colnames(disexp) <- c('s1','s2','patient','pair','pair_type','rna_d')
for(pat in rnapat) {
  curpat <- rnasam[grep(pat,rnasam)]
  for(s1 in c(1:length(curpat))) {
    sam1 <- curpat[s1]
    for(s2 in c(1:length(curpat))) {
      if(s1<s2) {
        sam2 <- curpat[s2]
        
        # square root(sum of the squared differences)
        expdiff <- sqrt(sum((geneexp[,sam1]-geneexp[,sam2])^2))
        
        reg1 <- gsub('C\\d+_(\\S)\\S+','\\1',sam1);reg2 <- paste0(gsub('C\\d+_(\\S)\\S+','\\1',sam2))
        type <- ifelse(reg1==reg2,'within-regions','between-regions')
        
        disexp[nrow(disexp)+1,] <- c(paste0('EPICC_',sam1),paste0('EPICC_',sam2),pat,paste0(reg1,'-',reg2),type,expdiff)
      }
    }
  }
}
disexp$rna_d <- as.numeric(disexp$rna_d)
disexp$pair_type <- factor(disexp$pair_type,levels=c('within-regions','between-regions'))
disexp <- disexp[-grep('E',disexp$pair),]
disexp <- disexp[-grep('C\\d+_\\S\\d+_B',disexp$s1),]
disexp <- disexp[-grep('C\\d+_\\S\\d+_B',disexp$s2),]

boxplot(rna_d~pair_type,names=F,xaxt='n',axes=F,las=2,cex.axis=1.2,frame=F,outline=FALSE,data=disexp,xlab='',ylab='',boxcol=scales::alpha(c('#E31A1C','#377DB8'),.9),boxlwd=2,col=scales::alpha(c('#E31A1C','#377DB8'),.4),medcol=scales::alpha(c('#E31A1C','#377DB8'),.9),medlwd=2.5,staplecol=scales::alpha(c('#E31A1C','#377DB8'),.9),staplelwd=1.5)
axis(2,font=2,las=2,cex.axis=1.4)
stripchart(rna_d~pair_type,vertical=T,data=disexp,method="jitter",add=TRUE,pch=20,col=scales::alpha(c('#E31A1C','#377DB8'),0.5),cex=2)

saveRDS(disexp,file='~/Documents/ThesisOther/ScriptsForThesis/RData/expression_distance.rds')

# 2. Plot boxplots of expression distance per patient ####
disexp <- readRDS('~/Documents/ThesisOther/ScriptsForThesis/RData/expression_distance.rds')
pdf('~/Documents/Thesis/figures/expression_distance_by_patient.pdf',width=20,height=10)
par(font=2,mar=c(2,5.5,1,0));pvals <- c();meddiff <- c()
layout(matrix(c(1:21),nrow=3,ncol=7,byrow = T))
for(pat in unique(disexp$patient)) {
  patdist <- disexp[which(disexp$patient==pat),]
  patdist$pair_type <- factor(patdist$pair_type,levels=c('within-regions','between-regions'))
  numpairs <- table(patdist$pair_type)
  if(numpairs[1]>1 & numpairs[2]>1) {
    meddiff <- c(meddiff,median(patdist[which(patdist$pair_type=='between-regions'),'rna_d'])-median(patdist[which(patdist$pair_type=='within-regions'),'rna_d']))
    test <- p.adjust(wilcox.test(rna_d~pair_type,alternative='less',data=patdist)$p.value,method='fdr',n=length(unique(disexp$patient)))
    pvals <- c(pvals,test)
    
    boxplot(rna_d~pair_type,names=F,yaxt='n',xaxt='n',ylim=c(0,800),las=2,cex.axis=1.2,frame=F,outline=FALSE,data=patdist,xlab='',ylab='',boxcol=scales::alpha(c('#E31A1C','#377DB8'),.9),boxlwd=2,col=scales::alpha(c('#E31A1C','#377DB8'),.4),medcol=scales::alpha(c('#E31A1C','#377DB8'),.9),medlwd=2.5,staplecol=scales::alpha(c('#E31A1C','#377DB8'),.9),staplelwd=1.5)
    axis(2,font=2,las=2,labels=c(0,200,400,600,800),at=c(0,200,400,600,800),cex.axis=1.4)
    mtext(side=3,line=-2,cex=1.4,text=pat)
    par(font=1);mtext(side=1,line=-1,cex=0.75,text=paste0('FDR adjusted p-value = ',signif(test,digits=2)));par(font=2)
    stripchart(rna_d~pair_type,vertical=T,data=patdist,method="jitter",add=TRUE,pch=20,col=scales::alpha(c('#E31A1C','#377DB8'),0.5),cex=2)
    
    if(pat=='C306') {
      legend(x=0.4,y=500,legend=c('within-regions','between-regions'),col=c('#E31A1C','#377DB8'),pch=16,bty='n',cex=2)
    }
    if(pat=='C537') {
      mtext(side=2,text='RNA-seq Distance',line=3.5,cex=1.4)
    }
  }
}
dev.off()

# 3. Apply mixed random effects model to expression distance ####
library(lme4);library(lmerTest)
pairres <- lmer(rna_d~pair_type+(1|patient),data=disexp)
sumres <- summary(pairres)
nullres <- lmer(rna_d~(1|patient),data=disexp)
pairnull_rna <- anova(pairres,nullres)

matres <- matrix(0L,nrow=1,ncol=6);colnames(matres) <- c('Fixed effect','Estimate','Standard Error','df','t value','$Pr (>|t|)$')
matres[1,] <- c('Pair',signif(sumres$coefficients[2,1],5),signif(sumres$coefficients[2,2],5),signif(sumres$coefficients[2,3],5),signif(sumres$coefficients[2,4],5),'$\\textless2e^{-16}$')
write.csv(matres,file='~/Documents/Thesis/mainText/expressiondistance_lmer.csv',row.names = F,quote=F)

compres <- matrix('',nrow=2,ncol=9);colnames(compres) <- c('Model','N','AIC','BIC','logLik','deviance','Chisq','Df','$Pr (>Chisq)$')
compres[1,] <- c('$R \\sim (1|Patient)$',pairnull_rna$npar[1],signif(pairnull_rna$AIC[1],digits=5),signif(pairnull_rna$BIC[1],digits=5),signif(pairnull_rna$logLik[1],digits=5),signif(pairnull_rna$deviance[1],digits=5),'','','')
compres[2,] <- c('$R \\sim Pair + (1|Patient)$',pairnull_rna$npar[2],signif(pairnull_rna$AIC[2],digits=5),signif(pairnull_rna$BIC[2],digits=5),signif(pairnull_rna$logLik[2],digits=5),signif(pairnull_rna$deviance[2],digits=5),signif(pairnull_rna$Chisq[2],digits=5),pairnull_rna$Df[2],'$\\textless2e^{-16}$')
write.csv(compres,file='~/Documents/Thesis/mainText/expressiondistance_lmer_chisq.csv',row.names = F,quote=F)


# 4. Compare genetic and expression distance ####
genomicdist <- readRDS('~/Documents/ThesisOther/ScriptsForThesis/RData/genetic_distance.rds')
colnames(genomicdist) <- c('s1','s2','patient','pair','pair_type','wgs_d')

dnarnadist <- merge(genomicdist,disexp,by=c('s1','s2','patient','pair','pair_type'))
paircols <- ifelse(dnarnadist$pair_type=='within-regions','#E31A1C','#377DB8')

resall <- lm(dnarnadist$rna_d~dnarnadist$wgs_d)
resbet <- lm(dnarnadist$rna_d[dnarnadist$pair_type=='between-regions']~dnarnadist$wgs_d[dnarnadist$pair_type=='between-regions'])
reswith <- lm(dnarnadist$rna_d[which(dnarnadist$pair_type=='within-regions')]~dnarnadist$wgs_d[dnarnadist$pair_type=='within-regions'])

corall <- cor.test(dnarnadist$wgs_d,dnarnadist$rna_d,method='spearman')
corbet <- cor.test(dnarnadist$wgs_d[dnarnadist$pair_type=='between-regions'],dnarnadist$rna_d[dnarnadist$pair_type=='between-regions'],method='spearman')
corwith <- cor.test(dnarnadist$wgs_d[dnarnadist$pair_type=='within-regions'],dnarnadist$rna_d[dnarnadist$pair_type=='within-regions'],method='spearman')

# 5. Plot genetic vs expression distance ####

pdf('~/Documents/Thesis/figures/genetic_vs_expression_distance_allsamples.pdf',width=17,height=15)
par(font=2,mar=c(6.5,6.1,6.6,2.1),cex.lab=2)
layout(matrix(c(rep(1,6),2:10),nrow=5,byrow = T))
plot(dnarnadist$wgs_d,dnarnadist$rna_d,pch=16,col=scales::alpha(rev(paircols),.6),cex.axis=2,cex=2,cex.lab=2.5,font=2,ylim=c(0,600),xlim=c(0,1),ylab='Expression Distance (sum of squares)',xlab='Genetic Distance (proportion of SNVs not shared)',font.lab=2)
title(main = paste0('All DNA-RNA samples (n=',length(unique(c(dnarnadist$s1,dnarnadist$s2))),')'),cex.main=3)
abline(resbet$coefficients[1],resbet$coefficients[2],col='skyblue3',lwd=2)
abline(reswith$coefficients[1],reswith$coefficients[2],col='firebrick3',lwd=2)
legend('bottom',legend=c('within-regions','between-regions'),col=c('#E31A1C','#377DB8'),pch=16,bty='n',cex=2)
text(x=-0.04,y=700,labels='(a)',xpd=T,font=2,cex=4)
par(mar=c(6.5,6.1,5.1,1.1))
for(pat in unique(dnarnadist$patient)) {
  patdist <- dnarnadist[which(dnarnadist$patient==pat),]
  numpairs <- table(patdist$pair_type)
  if(length(numpairs)>1) {
    if(numpairs[1]>1 & numpairs[2]>1) {
      paircols <- ifelse(patdist$pair_type=='within-regions','#E31A1C','#377DB8')
      resbettmp <- lm(patdist$rna_d[patdist$pair_type=='between-regions']~patdist$wgs_d[patdist$pair_type=='between-regions'])
      reswithtmp <- lm(patdist$rna_d[which(patdist$pair_type=='within-regions')]~patdist$wgs_d[patdist$pair_type=='within-regions'])
      
      corbettmp <- cor.test(patdist$wgs_d[patdist$pair_type=='between-regions'],patdist$rna_d[patdist$pair_type=='between-regions'],method='spearman')
      corwithtmp <- cor.test(patdist$wgs_d[patdist$pair_type=='within-regions'],patdist$rna_d[patdist$pair_type=='within-regions'],method='spearman')
      
      plot(patdist$wgs_d,patdist$rna_d,pch=16,xlim=c(0,1),cex.axis=2,cex=2,cex.lab=2,ylim=c(0,600),col=scales::alpha(paircols,.6),font=2,ylab='Expression Distance',xlab='Genetic Distance',font.lab=2)
      mtext(side=3,cex=2,line=1,text=pat)
      abline(resbettmp$coefficients[1],resbettmp$coefficients[2],col='#377DB8')
      abline(reswithtmp$coefficients[1],reswithtmp$coefficients[2],col='#E31A1C')
      text(c(0.8,0.8),c(100,20),labels=c(paste0('rho: ',signif(corwithtmp$estimate,2),' p: ',signif(corwithtmp$p.value,2)),
                                                    paste0('rho: ',signif(corbettmp$estimate,2),' p: ',signif(corbettmp$p.value,2))),cex=2,col=c('#E31A1C','#377DB8'))
      if(pat=='C524') {
        text(x=-0.04,y=750,labels='(b)',xpd=T,font=2,cex=4)
      }
      
    }
  }
  if(pat==unique(dnarnadist$patient)[length(unique(dnarnadist$patient))]) {
    plot.new()
    legend('center',legend=c('within-regions','between-regions'),col=c('#E31A1C','#377DB8'),pch=16,bty='n',cex=4)
  }
}
dev.off()

saveRDS(dnarnadist,file='~/Documents/ThesisOther/ScriptsForThesis/RData/genetic_vs_expression_distance.rds')

# 6. Apply mixed effects model to genetic vs expression distance ####
library(lme4);library(lmerTest)
wgsres <- lmer(rna_d~wgs_d+pair_type+(1|patient),data=dnarnadist)
sumres <- summary(wgsres)
wgsonly <- lmer(rna_d~wgs_d+(1|patient),data=dnarnadist)
pairres <- lmer(rna_d~pair_type+(1|patient),data=dnarnadist)
nullres <- lmer(rna_d~(1|patient),data=dnarnadist)

wgspair <- anova(wgsres,pairres)
wgsonlywgs <- anova(wgsonly,wgsres)
wgsonlynull <- anova(wgsonly,nullres)
wgsnull <- anova(wgsres,nullres)
pairnull_both <- anova(pairres,nullres)

matres <- matrix(0L,nrow=2,ncol=6);colnames(matres) <- c('Fixed effect','Estimate','Standard Error','df','t value','$Pr (>|t|)$')
matres[1,] <- c('D',signif(sumres$coefficients[2,1],5),signif(sumres$coefficients[2,2],5),signif(sumres$coefficients[2,3],5),signif(sumres$coefficients[2,4],5),signif(sumres$coefficients[2,5],5))
matres[2,] <- c('Pair',signif(sumres$coefficients[3,1],5),signif(sumres$coefficients[3,2],5),signif(sumres$coefficients[3,3],5),signif(sumres$coefficients[3,4],5),signif(sumres$coefficients[3,5],5))
write.csv(matres,file='~/Documents/Thesis/mainText/expvsgendis_lmer.csv',row.names = F,quote=F)

compres <- matrix('',nrow=2,ncol=9);colnames(compres) <- c('Model','N','AIC','BIC','logLik','deviance','Chisq','Df','$Pr (>Chisq)$')
compres[1,] <- c('$R \\sim (1|Patient)$',wgsonlynull$npar[1],signif(wgsonlynull$AIC[1],digits=5),signif(wgsonlynull$BIC[1],digits=5),signif(wgsonlynull$logLik[1],digits=5),signif(wgsonlynull$deviance[1],digits=5),'','','')
compres[2,] <- c('$R \\sim D + (1|Patient)$',wgsonlynull$npar[2],signif(wgsonlynull$AIC[2],digits=5),signif(wgsonlynull$BIC[2],digits=5),signif(wgsonlynull$logLik[2],digits=5),signif(wgsonlynull$deviance[2],digits=5),signif(wgsonlynull$Chisq[2],digits=5),wgsonlynull$Df[2],signif(wgsonlynull$`Pr(>Chisq)`[2],digits=5))
write.csv(compres,file='~/Documents/Thesis/mainText/modeldistanceresults.csv',row.names = F,quote=F)


