# Script to analyse and plot the results of eQTL analysis
# To do the following:
# 1) Load results of eQTL analysis
# 2) Plots pie charts of inital and final models
# 3) Plot bar plot showing proportion of data types in significant genes
# 4) Extract the coefficient and p-value for each gene and each data type
# 5) Plot the coefficient size of significant genes by data type
# 6) Plots volcano plots of Enh and Mut significant genes
# 7) Plot CXCL11 example data

# Functions ####
# Plot 'heatmap' of expression and all other factors involved in analysis for a particular gene
plotvar <- function(geneData,gene,namcol,patcol) {
  typecol <- rep('black',5)
  rescols <- wes_palette("Zissou1",100, type = "continuous");cnacols <- rep("gray85",nrow(geneData))
  if(is.numeric(geneData$CNA)) {;
    cnacols[geneData$CNA==3] <- wes_palette("BottleRocket1")[3];cnacols[geneData$CNA>3] <- wes_palette("BottleRocket1")[6]
    cnacols[geneData$CNA==1] <- wes_palette("Darjeeling1")[5];cnacols[geneData$CNA==0] <- wes_palette("Darjeeling2")[2]
  } else {
    cnacols[geneData$CNA=='Three'] <-wes_palette("BottleRocket1")[3];cnacols[geneData$CNA %in% c('Four','Amp')] <-wes_palette("BottleRocket1")[6]
    cnacols[geneData$CNA=='One'] <- wes_palette("Darjeeling1")[5];cnacols[geneData$CNA=='Zero'] <- wes_palette("Darjeeling2")[2]
  }
  samnegs <- rescale(c(1:nrow(geneData)))-rescale(c(1:nrow(geneData)))[2]/2;sampos <- rescale(c(1:nrow(geneData)))+rescale(c(1:nrow(geneData)))[2]/2
  geneData <- geneData[c(nrow(geneData):1),];row.names(geneData) <- c(1:nrow(geneData))
  
  ints <- seq.int(-100,100,by=30)/100
  image(t(scale(geneData$Expression)),axes=F,col=rescols)#;title(main=gene,cex=0.75,col.main='black')
  #if(namcol=='red2') { mtext(side=3,line=1.5,at=0.25,text='***',xpd=T,cex=2)}
  rect(rep(ints[3],nrow(geneData)),samnegs,rep(ints[2],nrow(geneData)),sampos,border=NA,col=ifelse(geneData$ATAC=='Closed','gray85',datacol[1]))
  rect(rep(ints[4],nrow(geneData)),samnegs,rep(ints[3],nrow(geneData)),sampos,border=NA,col=cnacols)
  rect(rep(ints[5],nrow(geneData)),samnegs,rep(ints[4],nrow(geneData)),sampos,border=NA,col=ifelse(geneData$Enh=='WT','gray85',datacol[3]))
  rect(rep(ints[6],nrow(geneData)),samnegs,rep(ints[5],nrow(geneData)),sampos,border=NA,col=ifelse(geneData$Mut=='WT','gray85',datacol[4]))
  rect(rep(ints[7],nrow(geneData)),samnegs,rep(ints[6],nrow(geneData)),sampos,border=NA,col=ifelse(geneData$Tissue=='Normal','gray85',datacol[5]))
  rect(rep(1,nrow(geneData)),samnegs,rep(ints[7],nrow(geneData)),sampos,border=NA,col=rev(patcol))
  
  mtext(side=2,at=rescale(c(1:nrow(geneData))),text=rev(matchsam),las=2,font=2,cex=0.28,line=0.2)
  mtext(side=3,at=c(ints[c(1:6)]+.15,0.9),text=c('Exp','ATAC','CNA','Enh','Mut','Tis','Pat'),line=0.1,font=2,cex=0.7,col=c('black',typecol,'black'))
  abline(h=c(samnegs,1.00365),col='white',lwd=0.2);abline(v=ints)
}

# Plot boxplots of expression vs various factors for a particular gene
plotboxes <- function(geneData) {
  cnacols <- rep("gray85",nrow(geneData));cnacols[geneData$CNA==3] <-"firebrick2";cnacols[geneData$CNA>3] <-"firebrick4"
  cnacols[geneData$CNA==1] <- "blue2";cnacols[geneData$CNA==0] <- "blue4"
  maxexp <- ceiling(max(geneData$Expression));minexp <- floor(min(geneData$Expression))
  
  if(length(table(geneData$ATAC))>1) {
    boxplot(Expression~ATAC,ylim=c(minexp,maxexp),las=1,cex.axis=1.125,frame=F,outline=F,data=geneData,xlab='',ylab='',boxcol=c('gray85',datacol[1]),boxlwd=1,col='white',medcol=c('gray70',datacol[1]),medlwd=1.5,staplecol=c('gray70',datacol[1]),staplelwd=1,outcol=c('gray70',datacol[1]),outpch=20)
    stripchart(Expression~ATAC,vertical=T,data=geneData,method="jitter",add=TRUE,pch=20,col=scales::alpha(c('gray70',datacol[1]),0.5),cex=1.5)
    mtext(side=3,line=0.2,cex=0.75,font=2,text=paste0('ATAC, wilcoxon p: ',signif(wilcox.test(Expression~ATAC,data=geneData)$p.value,2)))
    text(x=0.5,y=maxexp+(maxexp*.25),labels='(b)',cex=2,xpd=T,font=2)
  } else {
    plot.new()
    legend('center',legend = 'No variation in ATAC',bty='n',cex=1.125)
    text(x=0,y=1.1,labels='(b)',cex=2,xpd=T,font=2)
  }
  if(length(table(geneData$CNA))>1) {
    cnacols <- rep('gray70',length(table(geneData$CNA)));cnacols[as.numeric(names(table(geneData$CNA)))==3] <- wes_palette("BottleRocket1")[3];cnacols[as.numeric(names(table(geneData$CNA)))>3] <- wes_palette("BottleRocket1")[6]
    cnacols[as.numeric(names(table(geneData$CNA)))==1] <- wes_palette("Darjeeling1")[5];cnacols[as.numeric(names(table(geneData$CNA)))==0] <- wes_palette("Darjeeling2")[2]
    boxplot(Expression~CNA,ylim=c(minexp,maxexp),las=1,cex.axis=1.125,frame=F,outline=FALSE,data=geneData,xlab='',ylab='',boxcol=cnacols,boxlwd=1,col='white',medcol=cnacols,medlwd=1.5,staplecol=cnacols,staplelwd=1)
    stripchart(Expression~CNA,vertical=T,data=geneData,method="jitter",add=TRUE,pch=20,col=scales::alpha(cnacols,0.5),cex=1.5)
    mtext(side=3,line=0.2,cex=0.75,font=2,text=paste0('CNA, kruskal-wallis p: ',signif(kruskal.test(Expression~CNA,data=geneData)$p.value,3)))
    text(x=0.5,y=maxexp+(maxexp*.25),labels='(c)',cex=2,xpd=T,font=2)
  } else {
    plot.new()
    legend('center',legend = 'No variation in CNA',bty='n',cex=1.125)
    text(x=0,y=1.1,labels='(c)',cex=2,xpd=T,font=2)
  }
  if(length(table(geneData$Enh))>1) {
    geneData$Enh <- factor(geneData$Enh,levels=c('WT','Mut'))
    boxplot(Expression~Enh,ylim=c(minexp,maxexp),las=1,cex.axis=1.125,frame=F,outline=F,data=geneData,xlab='',ylab='',boxcol=c('gray85',datacol[3]),boxlwd=1,col='white',medcol=c('gray70',datacol[3]),medlwd=1.5,staplecol=c('gray70',datacol[3]),staplelwd=1,outcol=c('gray70',datacol[3]),outpch=20)
    stripchart(Expression~Enh,vertical=T,data=geneData,method="jitter",add=TRUE,pch=20,col=scales::alpha(c('gray70',datacol[3]),0.5),cex=1.5)
    mtext(side=3,line=0.2,cex=0.75,font=2,text=paste0('Enhancer Mutation, wilcoxon p: ',signif(wilcox.test(Expression~Enh,data=geneData)$p.value,2)))
    text(x=0.5,y=maxexp+(maxexp*.25),labels='(d)',cex=2,xpd=T,font=2)
  } else {
    plot.new()
    legend('center',legend = 'No variation in Enhancer Mutation',bty='n',cex=1.125)
    text(x=0,y=1.1,labels='(d)',cex=2,xpd=T,font=2)
  }
  if(length(table(geneData$Mut))>1) {
    geneData$Mut <- factor(geneData$Mut,levels=c('WT','Mut'))
    boxplot(Expression~Mut,ylim=c(minexp,maxexp),las=1,cex.axis=1.125,frame=F,outline=F,data=geneData,xlab='',ylab='',boxcol=c('gray85',datacol[4]),boxlwd=1,col='white',medcol=c('gray70',datacol[4]),medlwd=1.5,staplecol=c('gray70',datacol[4]),staplelwd=1,outcol=c('gray70',datacol[4]),outpch=20)
    stripchart(Expression~Mut,vertical=T,data=geneData,method="jitter",add=TRUE,pch=20,col=scales::alpha(c('gray70',datacol[4]),0.5),cex=1.5)
    mtext(side=3,line=0.2,cex=0.75,font=2,text=paste0('Mutation, wilcoxon p: ',signif(wilcox.test(Expression~Mut,data=geneData)$p.value,2)))
    text(x=0.5,y=maxexp+(maxexp*.25),labels='(e)',cex=2,xpd=T,font=2)
  } else {
    plot.new()
    legend('center',legend = 'No variation in Genic Mutation',bty='n',cex=1.125)
    text(x=0,y=1.1,labels='(e)',cex=2,xpd=T,font=2)
  }
  
  boxplot(Expression~Tissue,ylim=c(minexp,maxexp),las=1,cex.axis=1.125,frame=F,outline=F,data=geneData,xlab='',ylab='',boxcol=c('gray85',datacol[5]),boxlwd=1,col='white',medcol=c('gray70',datacol[5]),medlwd=1.5,staplecol=c('gray70',datacol[5]),staplelwd=1,outcol=c('gray70',datacol[5]),outpch=20)
  stripchart(Expression~Tissue,vertical=T,data=geneData,method="jitter",add=TRUE,pch=20,col=scales::alpha(c('gray70',datacol[5]),0.5),cex=1.5)
  mtext(side=3,line=0.2,cex=0.75,font=2,text=paste0('Tissue, wilcoxon p: ',signif(wilcox.test(Expression~Tissue,data=geneData)$p.value,2)))
  text(x=0.5,y=maxexp+(maxexp*.25),labels='(f)',cex=2,xpd=T,font=2)
  
}

# Main code ####

# 0. Prepare environment
setwd('~/Documents/EPICC/Data/QTLs/allvariants/')
library(data.table);library(dplyr);'%ni%' <- Negate('%in%');library(DESeq2)
library(RColorBrewer);library(scales);library(wesanderson)
datacol <- wes_palette("Cavalcanti1")


# 1. Load the results of the eQTL analysis from 5.5.1 ####
modelres <- readRDS('qtl_model_res_latest.cnanum.nolog.allvars.zscore.rds')
sigres <- modelres[which(modelres$padj<0.01),];row.names(sigres) <- c(1:nrow(sigres))
listmat <- readRDS('matricesOfModel.latest.cnanum.nolog.allvars.zscore.rds')
matchsam <- colnames(listmat$Expression)

# 2. Plot pies of the starting models and final significant models ####
pdf('~/Documents/Thesis/figures/eQTL_pies.pdf')
layout(matrix(c(1,2),nrow=2))
par(font=2,mar=c(0,2,2,21.2),xpd=T)
pie(table(modelres$MaxModel),clockwise = T,labels = NA,radius=1.2,col=RColorBrewer::brewer.pal(8,'Set1'),border = NA)
legend(x=1.2,y=0.6,legend=paste0(unique(modelres$MaxModel),', n=',as.numeric(table(modelres$MaxModel)[unique(modelres$MaxModel)])),xpd = T,fill=RColorBrewer::brewer.pal(8,'Set1'),border = NA,bty='n')
text(x=1.5,y=1.35,labels=paste0('Starting models (n = ',nrow(modelres),' genes)'),cex=1.5)
text(x=-1,y=1.3,'(a)',cex=2,font=2)

colmod <- rev(c(brewer.pal(12,'Paired'),brewer.pal(12,'Set3'),brewer.pal(6,'Dark2')))
pie(table(sigres$FinalModel),clockwise = T,labels = NA,radius=1.2,col=colmod,border = NA)
legend(x=1.2,y=1.1,cex=.7,legend=paste0(sort(unique((sigres$FinalModel)))[c(1:15)],', n=',as.numeric(table(sigres$FinalModel)[sort(unique(sigres$FinalModel))][c(1:15)])),xpd = T,fill=colmod[c(1:15)],border = NA,bty='n')
legend(x=3.05,y=1.1,cex=.7,legend=paste0(sort(unique((sigres$FinalModel)))[c(16:30)],', n=',as.numeric(table(sigres$FinalModel)[sort(unique(sigres$FinalModel))][c(16:30)])),xpd = T,fill=colmod[c(16:30)],border = NA,bty='n')
text(x=1.5,y=1.35,labels=paste0('Significant models (n = ',nrow(sigres),' genes)'),cex=1.5)
text(x=-1,y=1.3,'(b)',cex=2,font=2)
dev.off()

# 3. Plot barplot of significant gene proportions by data type ####
pdf('~/Documents/Thesis/figures/eQTL_results_by_data_type.pdf',width=8)
par(font=2,mar=c(3.6,5.6,3.1,1.1))
startnum <- rep(0,5);names(startnum) <- c('ATAC','CNA','Enh','Mut','Tissue');for(dat in names(startnum)) { startnum[dat] <- length(grep(dat,modelres$MaxModel))}
signum <- rep(0,5);names(signum) <- c('ATAC','CNA','Enh','Mut','Tissue');for(dat in names(signum)) { signum[dat] <- length(grep(dat,sigres$FinalModel))}
sigprop <- signum/startnum*100
xx<- barplot(sigprop,ylim=c(0,30),col=datacol,border=NA,las=1,font=2,cex.axis=1.8,cex.names = 2.1)
text(x=xx,y=sigprop+1,label=paste0(signum,'/',startnum),cex=1.5,font=2,col='black',xpd=T)
mtext(side=2,text='Percentage significant genes (%)',line=3.3,cex=1.8,xpd=T)
dev.off()

# 4. Get the size of coefficients and p-values per data type ####
coefmat <- matrix(0L,nrow=nrow(sigres),ncol=5);colnames(coefmat) <- c('ATAC','CNA','Enh','Mut','Tissue')
for(var in colnames(coefmat)) {
  for(i in c(1:nrow(sigres))) {
    coefmat[i,var] <- ifelse(length(grep(var,sigres$VarMod[i]))==1,as.numeric(gsub(paste0('.*',var,':(-\\d+\\.\\d+|\\d+\\.\\d+|-\\d+|\\d+).*'),'\\1',sigres$VarMod[i])),NA)
  }
}
pvalmat <- matrix(0L,nrow=nrow(sigres),ncol=5);colnames(pvalmat) <- c('ATAC','CNA','Enh','Mut','Tissue')
for(var in colnames(pvalmat)) {
  for(i in c(1:nrow(sigres))) {
    pvalmat[i,var] <- ifelse(length(grep(var,sigres$VarP[i]))==1,as.numeric(gsub(paste0('.*',var,':(-\\d+\\.\\d+|\\d+\\.\\d+|-\\d+|\\d+).*'),'\\1',sigres$VarP[i])),NA)
  }
}

# 5. Plot the size of coefficients of significant genes by data type ####
pdf('~/Documents/Thesis/figures/eQTL_effect_size.pdf')
library(vioplot);par(font=2,font.axis=2,cex.axis=1.8,mar=c(3.1,4.5,1.1,1.1))
vioplot(coefmat,col=datacol,ylim=c(-4,6),border=datacol,axes=F,las=1)
abline(h=0,lty=2)
vioplot(coefmat,col=datacol,ylim=c(-4,6),border=datacol,axes=F,las=1,add=T)
stripchart(list(coefmat[,1],coefmat[,2],coefmat[,3],coefmat[,4],coefmat[,5]),vertical=T,method="jitter",add=TRUE,pch=20,col=scales::alpha('dimgray',0.15),cex=0.5)
mtext(side=2,text='Regression coefficient',line=2.8,cex=1.8)
dev.off()

# 6. Plot volcano plots of Mutations and Enhancer significant genes ####
pdf('~/Documents/Thesis/figures/eQTL_volcano_plots.pdf',height=11,width=7)
layout(matrix(c(1,2),nrow=2))
par(mar=c(4.6,4.6,3,1),font.lab=2)
var <- 'Mut';cutp <- 2;cutreg <- 1.96
tmpco <- coefmat[which(!is.na(coefmat[,var])),var]
tmpp <- log(p.adjust(pvalmat[which(!is.na(pvalmat[,var])),var],method='fdr'),base=10)*-1

plotcol <- rep('dimgray',length(tmpp))
plotcol[tmpp>=cutp & abs(tmpco)>=cutreg] <- datacol[which(colnames(pvalmat)==var)]
genes <- rep('',length(tmpp))
tmpgen <- sigres[grep(var,sigres$FinalModel),];row.names(tmpgen) <- c(1:nrow(tmpgen))
genes <- tmpgen$hgnc_symbol[which(tmpp>=cutp & abs(tmpco)>=cutreg)]
labtab <- data.frame(Gene=genes,xLab=tmpco[which(abs(tmpco)>=cutreg & tmpp>=cutp)],yLab=tmpp[which(abs(tmpco)>=cutreg & tmpp>=cutp)])

plot(tmpco,tmpp,pch=16,xlim=c(-6,6),las=1,ylim=c(0,8),main=var,cex.main=1.75,col=scales::alpha(plotcol,0.85),bty='n',xlab='',ylab='',cex.axis=1.7,font.axis=2)
mtext(side=2,text='-log10(pval)',line=2.8,cex=1.7,font=2)
mtext(side=1,text='Regression coefficient',line=2.7,cex=1.7,font=2)
abline(h=cutp,lty=2,col=datacol[which(colnames(pvalmat)==var)])
abline(v=c(cutreg*-1,cutreg),lty=2,col=datacol[which(colnames(pvalmat)==var)])
abline(v=0,lty=3,col='dimgray')

# Manually position gene labels and lines
labtab$newX <- labtab$xLab+c(-1,-2.5,1.5,-0.5,-2.5,-2.5,-1,-1,0.2,-2.5,-2,1.5,1,-2.5,-0.2,2.5,-2,2.5,1.5,-3,2,3,1,0.7,1.5)
labtab$newY <- labtab$yLab+c(2.5,0.5,-1.5,-2,-1.5,0.8,-1.5,1.5,1,-0.2,0.5,1,-2,1,1,-1,1.2,-0.2,1,-1,0.2,0.8,0.5,2,0.2)
text(x=labtab$newX,y=labtab$newY,labels=labtab$Gene,cex=1.1,xpd=T,font=2)
segments(x0=labtab$xLab,x1=labtab$newX+c(0,0.7,0,0,0.6,0.7,0,0,0,0.8,0.6,-0.8,0,0.7,0,-0.7,0.6,-0.7,-0.7,0.7,-0.7,-0.7,-0.8,0,-0.6),
         y0=labtab$yLab,y1=labtab$newY+c(-0.15,0,0.15,0.15,0,0,0.15,-0.15,-0.15,0,0,0,0.15,0,-0.15,0,0,0,0,0,0,0,0,-0.15,0),
         col='black',lty=2)
text(x=-7,y=9,labels='(a)',cex=2,font=2,xpd=T)

par(mar=c(4.6,4.6,3,1),font.lab=2)
var <- 'Enh';cutp <- 2;cutreg <- 1.96
tmpco <- coefmat[which(!is.na(coefmat[,var])),var]
tmpp <- log(p.adjust(pvalmat[which(!is.na(pvalmat[,var])),var],method='fdr'),base=10)*-1

plotcol <- rep('dimgray',length(tmpp))
plotcol[tmpp>=cutp & abs(tmpco)>=cutreg] <- datacol[which(colnames(pvalmat)==var)]
genes <- rep('',length(tmpp))
tmpgen <- sigres[grep(var,sigres$FinalModel),];row.names(tmpgen) <- c(1:nrow(tmpgen))
genes <- tmpgen$hgnc_symbol[which(tmpp>=cutp & abs(tmpco)>=cutreg)]
labtab <- data.frame(Gene=genes,xLab=tmpco[which(abs(tmpco)>=cutreg & tmpp>=cutp)],yLab=tmpp[which(abs(tmpco)>=cutreg & tmpp>=cutp)])

plot(tmpco,tmpp,pch=16,xlim=c(-6,6),las=1,ylim=c(0,8),main=var,cex.main=1.75,col=scales::alpha(plotcol,0.85),bty='n',xlab='',ylab='',cex.axis=1.7,font.axis=2)
mtext(side=2,text='-log10(pval)',line=2.8,cex=1.7,font=2)
mtext(side=1,text='Regression coefficient',line=2.7,cex=1.7,font=2)
abline(h=cutp,lty=2,col=datacol[which(colnames(pvalmat)==var)])
abline(v=c(cutreg*-1,cutreg),lty=2,col=datacol[which(colnames(pvalmat)==var)])
abline(v=0,lty=3,col='dimgray')

# Manually position gene labels and lines
labtab$newX <- labtab$xLab+c(1.5,2.5,3,1,1.1,0.5,-1.5,3,0.5,1,2.8,2,1.5)
labtab$newY <- labtab$yLab+c(-2,-0.5,-2,0.5,0.1,1,1,-0.2,1,-0.4,-1.4,0.2,0.5)
text(x=labtab$newX,y=labtab$newY,labels=labtab$Gene,cex=1.1,xpd=T,font=2)
segments(x0=labtab$xLab,x1=labtab$newX+c(0,-0.7,0,0,-0.7,0,0,-0.7,0,0,-0.7,-0.7,-0.7),
         y0=labtab$yLab,y1=labtab$newY+c(0.15,0,0.15,0.15,0,-0.15,-0.15,0,-0.15,0.15,0,0,0),
         col='black',lty=2)
text(x=-7,y=9,labels='(b)',cex=2,font=2,xpd=T)
dev.off()

# 7. Plot heatmap and boxplots of CXCL11 ####
pdf('~/Documents/Thesis/figures/eQTL_CXCL11_example.pdf')
layout(matrix(c(1,1,1,1,1,2,3,4,5,6),ncol=2))
otherres <- sigres[which(sigres$hgnc_symbol=='CXCL11'),]
for(i in c(1:nrow(otherres))) {
  curmod <- otherres[i,]
  genename <- curmod$hgnc_symbol;gen <- curmod$ensembl_gene_id
  
  geneData <- data.frame(Expression=listmat[["Expression"]][gen,],ATAC=listmat[["ATAC"]][gen,],CNA=listmat[["CNA"]][gen,],Enh=listmat[["Enh"]][gen,],Mut=listmat[["Mut"]][gen,],Tissue=listmat[["Tissue"]],Patient=as.factor(gsub('(C\\d+)_\\S\\d+_\\S\\d+','\\1',matchsam)))
  
  twocol <- c('gray55','gray80');matchpat <- gsub('(C\\d+)_\\S+','\\1',matchsam)
  patcol <- c('gray80');for(j in c(2:length(matchsam))) { patcol <- c(patcol,ifelse(matchpat[j]==matchpat[j-1],patcol[j-1],twocol[which(twocol!=patcol[j-1])])) }
  
  namcol <- ifelse(curmod$padj<0.01,'red2','black')
  
  if(curmod$padj<0.01) {
    par(mar=c(.5,3.5,2.5,.5),font=2)
    plotvar(geneData,genename,namcol,patcol)
    mtext(side=3,text=paste0(genename,': Final model - Exp ~ ',curmod$FinalModel,', padj=',formatC(curmod$padj, format = "e", digits = 2)),font=2,line=1,cex=.7)
    text(x=-1.2,y=1.03,labels='(a)',cex=2,xpd=T,font=2)
    par(mar=c(2.5,2.5,2.5,1),font=2,font.axis=2,font.lab=2)
    plotboxes(geneData)
  }
}
dev.off()


