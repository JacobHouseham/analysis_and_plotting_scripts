# Script to assess genetic distance within and between regions of tumours
# To do the following:
# 1) Load data
# 2) Determine genetic distance for all within-tumour sample comparisons
# 3) Plot results of genetic distance by patient
# 4) Do mixed effects model of genetic distance by pair-type

# Main code ####

# 0. Prepare environment
regcol <- c(A='#E31A1C',B='#377DB8',C='#4DAE49',D='#904A9A',Root='#808080');'%ni%' <- Negate('%in%')
setwd("~/Documents/EPICC/Data/Mutations");library(dplyr);library(data.table);library(phangorn)
allpats <- gsub('(C\\d+)_\\S+','\\1',multisam);patients <- unique(allpats)
options(scipen = -1)

# 1. Load SNV and sample data ####
snvmat <- readRDS('~/Documents/EPICC/Data/QTLs/nodesamples/snvmatrix_allDNAsamples.rds')
wgssam <- readRDS('~/Documents/EPICC/Data/QTLs/nodesamples/wgs_used_in_phylo.rds')
wgssam <- gsub('EPICC_(C\\d+_\\S\\d+_\\S\\d+)\\S+','\\1',wgssam)
samples <- colnames(snvmat);samples <- samples[which(samples %in% wgssam)]
samples <- samples[grep('C\\d+_(A|B|C|D)\\d+_G',samples)]
wgspat <- unique(gsub('(C\\d+)_\\S+','\\1',samples))

# 2. Get genetic distance between samples within tumours ####
disgen <- data.frame(matrix(ncol=6,nrow=0));colnames(disgen) <- c('s1','s2','patient','pair','pair_type','wgs_d')
for(pat in wgspat) {
  curpat <- samples[grep(pat,samples)]
  for(s1 in c(1:length(curpat))) {
    sam1 <- curpat[s1]
    print(sam1)
    for(s2 in c(1:length(curpat))) {
      if(s1<s2) {
        sam2 <- curpat[s2]
        
        snvs1 <- row.names(snvmat[which(snvmat[,sam1]==1),])
        snvs2 <- row.names(snvmat[which(snvmat[,sam2]==1),])
        
        gendiff <- as.numeric((length(which(snvs2 %ni% snvs1))+length(which(snvs1 %ni% snvs2)))/sum(length(unique(c(snvs1,snvs2)))))
        
        reg1 <- gsub('C\\d+_(\\S)\\S+','\\1',sam1);reg2 <- paste0(gsub('C\\d+_(\\S)\\S+','\\1',sam2))
        type <- ifelse(reg1==reg2,'within-regions','between-regions')
        
        disgen[nrow(disgen)+1,] <- c(paste0('EPICC_',sam1),paste0('EPICC_',sam2),pat,paste0(reg1,'-',reg2),type,gendiff)
      }
    }
  }
}
disgen$wgs_d <- as.numeric(disgen$wgs_d)
saveRDS(disgen,file='~/Documents/ThesisOther/ScriptsForThesis/RData/genetic_distance.rds')

# 3. Plot boxplots of genetic distance per patient ####
disgen <- readRDS('~/Documents/ThesisOther/ScriptsForThesis/RData/genetic_distance.rds')

pdf('~/Documents/Thesis/figures/genetic_distance_by_patient.pdf',width=20,height=10)
par(font=2,mar=c(4,6,1.1,0));pvals <- c()
layout(matrix(c(1:28),nrow=4,ncol=7,byrow = T))
for(pat in unique(disgen$patient)) {
  patdist <- disgen[which(disgen$patient==pat),]
  
  patdist$pair_type <- factor(patdist$pair_type,levels=c('within-regions','between-regions'))
  numpairs <- table(patdist$pair_type)
  if(numpairs[1]>=1 & numpairs[2]>=1) {
    test <- p.adjust(wilcox.test(wgs_d~pair_type,alternative='less',data=patdist)$p.value,method='fdr',n=length(unique(disgen$patient)))
    pvals <- c(pvals,test)
    
    boxplot(wgs_d~pair_type,names=F,yaxt='n',xaxt='n',ylim=c(0,1),las=2,cex.axis=1.2,frame=F,outline=FALSE,data=patdist,xlab='',ylab='',boxcol=scales::alpha(c('#E31A1C','#377DB8'),.9),boxlwd=2,col=scales::alpha(c('#E31A1C','#377DB8'),.4),medcol=scales::alpha(c('#E31A1C','#377DB8'),.9),medlwd=2.5,staplecol=scales::alpha(c('#E31A1C','#377DB8'),.9),staplelwd=1.5)
    axis(2,font=2,las=2,labels=seq(0,1,by=0.2),at=seq(0,1,by=0.2),cex.axis=1.4)
    mtext(side=3,line=-2,cex=1.4,text=pat)
    par(font=1);mtext(side=1,line=-.5,cex=0.75,text=paste0('FDR adjusted p-value = ',signif(test,digits=2)));par(font=2)
    stripchart(wgs_d~pair_type,vertical=T,data=patdist,method="jitter",add=TRUE,pch=20,col=scales::alpha(c('#E31A1C','#377DB8'),0.5),cex=2)
    
    if(pat=='C543') {
      mtext(side=2,text='WGS Distance',line=4,cex=1.4,xpd=T,at=1.2)
    }
  }
}
plot.new()
legend(x=0,y=1,legend=c('within-regions','between-regions'),col=c('#E31A1C','#377DB8'),pch=16,bty='n',cex=1.65)
dev.off()

# 4. Analyse mixed effect model and output results ####
library(lme4);library(lmerTest)

pairres <- lmer(wgs_d~pair_type+(1|patient),data=disgen)
sumres <- summary(pairres)
nullres <- lmer(wgs_d~(1|patient),data=disgen)
pairnull_wgs <- anova(pairres,nullres)

matres <- matrix(0L,nrow=1,ncol=6);colnames(matres) <- c('Fixed effect','Estimate','Standard Error','df','t value','$Pr (>|t|)$')
matres[1,] <- c('Pair',signif(sumres$coefficients[2,1],5),signif(sumres$coefficients[2,2],5),signif(sumres$coefficients[2,3],5),signif(sumres$coefficients[2,4],5),'$\\textless2e^{-16}$')
write.csv(matres,file='~/Documents/Thesis/mainText/geneticdistance_lmer.csv',row.names = F,quote=F)

compres <- matrix('',nrow=2,ncol=9);colnames(compres) <- c('Model','N','AIC','BIC','logLik','deviance','Chisq','Df','$Pr (>Chisq)$')
compres[1,] <- c('$D \\sim (1|Patient)$',pairnull_wgs$npar[1],signif(pairnull_wgs$AIC[1],digits=5),signif(pairnull_wgs$BIC[1],digits=5),signif(pairnull_wgs$logLik[1],digits=5),signif(pairnull_wgs$deviance[1],digits=5),'','','')
compres[2,] <- c('$D \\sim Pair + (1|Patient)$',pairnull_wgs$npar[2],signif(pairnull_wgs$AIC[2],digits=5),signif(pairnull_wgs$BIC[2],digits=5),signif(pairnull_wgs$logLik[2],digits=5),signif(pairnull_wgs$deviance[2],digits=5),signif(pairnull_wgs$Chisq[2],digits=5),pairnull_wgs$Df[2],'$\\textless2e^{-16}$')
write.csv(compres,file='~/Documents/Thesis/mainText/geneticdistance_lmer_chisq.csv',row.names = F,quote=F)

