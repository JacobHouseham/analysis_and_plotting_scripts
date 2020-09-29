# Script to assess the correlation of CNAs with gene expression
# To do the following:
# 1) Load sample data
# 2) Load gene expression data and reformat
# 3) Load matrix of CNAs by gene per sample
# 4) Assess correlation of CNAs and expression
# 5) Do inference of CNAs by gene expression 
# 6) Plot CNA profiles and gene expression rolling averages together

# Functions ####
# Plot the copy number with gene expression rolling average
plotcn <- function(sam,segs,chrInfo,regionhg38,plotx,colreg) {
  
  segs <- segs[which(segs$chromosome %in% chrInfo$chromosome),]
  segs <- segs[which(!is.na(segs$Status)),];row.names(segs) <- c(1:nrow(segs))
  chrIDs <- segs$chromosome;gPositions <- chrInfo[['plotStart']][match(chrIDs,chrInfo$chromosome)]
  startPos <- segs$start.pos+gPositions;endPos <- segs$end.pos+gPositions
  
  plot(x=c(0,max(endPos)),y=c(0,5),font.axis=2,col='white',xaxt='n',main='',xlab='',ylab='',cex.main=1.1,xaxs='i',las=1)
  grid(nx=NA,ny=NULL,lty=1)
  
  scaled <- rescale(regionhg38[,sam],c(0,5))
  toplot <- ((regionhg38$NewEnd-regionhg38$NewStart)/2)+regionhg38$NewStart
  rect(xleft=toplot-5e6,xright=toplot+5e6,ybottom=(scaled-0.08),ytop=(scaled+0.08),col=alpha(colreg,0.8),border=NA)
  
  rect(xleft=startPos,xright=endPos,ybottom=(segs$Status-0.1),ytop=(segs$Status+0.1),col='blue3',border='blue3')
  
  if(max(segs$Status)>5) {
    above5 <- which(segs$Status>5)
    rect(xleft=startPos[above5],xright=endPos[above5],ybottom=4.9,ytop=5.1,col='blue4',border='blue4')
    text(x=startPos[above5][1],y=4.75,label=paste0('Max CN: ',max(segs$Status)),cex=0.75,font=3)
  }
  text(x=1e8,y=4.5,labels=sam,cex=1)
  abline(v=chrInfo$plotCentroStart,lty=2,lwd=0.8,col='lightgray')
  abline(v=chrInfo$plotCentroEnd,lty=2,lwd=0.8,col='lightgray')
  abline(v=chrInfo$plotStart,lty=1,lwd=0.8)
  abline(v=chrInfo$plotEnd,lty=1,lwd=0.8)
  if(plotx) { axis(side=1,at=chrInfo$plotCentroStart,labels=gsub('chr(\\S+)','\\1',chrInfo$chromosome),cex=0.75,font=2.5) }
  
}

# Main code ####

# 0. Prepare environment
convcn <- c('0'='Zero','1'='One','2'='Two','3'='Three','4'='Four','5'='Five','6'='Six','7'='Seven',
            '8'='Eight','9'='Nine','10'='Ten','11'='Eleven','12'='Twelve','13'='Thirteen','14'='Fourteen',
            '15'='Fourteen','16'='Sixteen','17'='Seventeen','18'='Eighteen','19'='Nineteen','20'='Twenty',
            '22'='TwentyTwo','24'='TwentyFour','25'='TwentyFive','27'='TwentySeven','29'='TwentyNine',
            '32'='ThirtyTwo','41'='FortyOne')
numcn <- as.numeric(names(convcn));names(numcn) <- unname(convcn)
regcols <- c(A='#E31A1C',B='#377DB8',C='#4DAE49',D='#904A9A')

setwd("~/Documents/EPICC/Data/Mutations");library(dplyr);library(data.table)
library(stringr);'%ni%' <- Negate('%in%');library(zoo);library(scales)
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
geneMap = getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), mart = ensembl)

# 1. Load sample data  ####

rnasam <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]
tsam <- rnasam[-grep('_E\\d+',rnasam)]
wgssam <- readRDS('~/Documents/EPICC/Data/QTLs/nodesamples/wgs_used_in_phylo.rds')
wgssam <- gsub('EPICC_(C\\S+)_D1','\\1',wgssam)

# Add in low pass
lpsam <- read.table('~/Documents/EPICC/Data/CopyNumber/LowPass/lp_samplelist.txt',header=T)

cnainput <- readRDS('~/Documents/EPICC/Data/QTLs/revisedmodel/cnabygene_numeric.rds')
samples <- tsam[which(tsam %in% colnames(cnainput))]
allpats <- gsub('(C\\d+).+','\\1',samples);patients <- unique(allpats)

# 2. Load gene expression data, filter and transform for somatic epithelial expression ####
geneexp <- as.data.frame(fread('~/Documents/EPICC/Data/Expression/ProcessedCounts/All_EPICC_tpm.txt'))
filgenes <- read.table('~/Documents/ThesisOther/ScriptsForThesis/Tables/filteredgenes.20perAL1TPM.ensembl.txt')[,1]
row.names(geneexp) <- geneexp$GeneID
geneexp <- geneexp[filgenes,rnasam]

isellasupp <- as.data.frame(readxl::read_xlsx('~/Documents/EPICC/Data/Expression/Classification/CRIS/IsellaSupp11.xlsx')[,c(1:3)])
colnames(isellasupp) <- c('Mouse','GeneSymbol','PercMurine')
isellasupp <- isellasupp[order(isellasupp$GeneSymbol),]
nonstromal <- isellasupp[which(isellasupp$PercMurine<=0.25),'GeneSymbol']

nonstromal_ens <- geneMap[which(geneMap$hgnc_symbol %in% nonstromal),'ensembl_gene_id']
geneexp <- geneexp[which(row.names(geneexp) %in% nonstromal_ens),]
for(s in c(1:length(rnasam))) { geneexp[,rnasam[s]] <- log2(geneexp[,rnasam[s]]+1)}

tpats <- unique(gsub('(C\\d+)_\\S+','\\1',tsam))
normsam <- rnasam[grep('_E\\d+',rnasam)];normsam <- normsam[which(gsub('(C\\d+)_\\S+','\\1',normsam) %in% tpats)]

normexp <- rowMeans(geneexp[,normsam])
for(s in c(1:length(rnasam))) { geneexp[,rnasam[s]] <- geneexp[,rnasam[s]]-normexp }
geneexp$GeneID <- row.names(geneexp);geneexp <- geneexp[,c(ncol(geneexp),1:(ncol(geneexp)-1))]
row.names(geneexp) <- c(1:nrow(geneexp))

# 3. Load pre-compiled matrix of CNAs by gene and PGA data ####
cnamat <- cnainput[geneexp$GeneID,samples]
nagenes <- apply(cnamat, 1, function(x) sum(is.na(x)))
matchgene <- sort(row.names(cnamat)[which(nagenes==0)])

geneexp <- geneexp[which(geneexp$GeneID %in% matchgene),]
cnamat <- as.data.frame(cnamat[matchgene,]);cnamat$GeneID <- row.names(cnamat)
cnamat <- cnamat[,c(ncol(cnamat),1:(ncol(cnamat)-1))];row.names(cnamat) <- c(1:nrow(cnamat))

pgas <- readRDS('~/Documents/ThesisOther/ScriptsForThesis/Tables/copynumberPGA.rds')
pgas <- pgas[samples]
matchsam <- samples
matchsam <- matchsam[which(sapply(cnamat[,matchsam],function(x) length(table(x)))>1)]
cnamat <- cnamat[,c('GeneID',matchsam)]
geneexp <- geneexp[,c('GeneID',matchsam)]

# 4. Get correlation of CNA and expression by matched DNA-RNA sample ####
expcors <- data.frame(Sample=matchsam,Type=ifelse(matchsam %in% wgssam,'Deep','LP'),PGA=pgas[matchsam],Correlation=0L,Gradient=0L,Intercept=0L,Rsquared=0L,Pval=0L)
for(s in c(1:length(matchsam))) {
  sam <- matchsam[s]
  ExpCN <- data.frame(Status=cnamat[,sam],Exp=geneexp[,sam])
  
  expcors[s,"Correlation"] <-  cor(ExpCN$Status,ExpCN$Exp,method='pearson')
  res <-lm(Exp~Status,data=ExpCN);ressum <- summary(res)
  expcors[s,"Gradient"] <- as.numeric(res$coefficients[2])
  expcors[s,"Intercept"] <- as.numeric(res$coefficients[1])
  expcors[s,"Rsquared"] <- ressum$r.squared
  if(nrow(ressum$coefficients)>1) {
    expcors[s,"Pval"] <- ressum$coefficients[2,4]
  } else {
    expcors[s,"Pval"] <- 1
  }
}
expcors$Padj <- p.adjust(expcors$Pval)

barplot(expcors[,'Gradient'],main='Red means significant correlation (FDR<0.01)',las=2,cex.names = 0.4,col=ifelse(expcors[,'Padj']<0.01,'red2','dimgray'),ylab='Correlation of CN and Expression (ρ)',border=NA)

plot(1, type="n", xlab="", ylab="", xlim=c(min(expcors$Intercept), max(expcors$Intercept)), ylim=c(min(expcors$Gradient), max(expcors$Gradient)))
for(s in c(1:nrow(expcors))) {
  sam <- expcors[s,'Sample']
  abline(a=expcors[s,'Intercept'],b=expcors[s,'Gradient'],lty=2,col=ifelse(expcors[s,'Padj']<0.01,'red2','dimgray'))
}

boxplot(expcors[which(expcors$Type=='Deep'),'Gradient'],expcors[which(expcors$Type=='LP'),'Gradient'],names=c('Deep (n=72)','LP (n=81)'),col=c('gray40','gray85'));abline(h=0,lty=2)

msi_positiv = c("C536","C548","C516","C518","C552")

# 5. Get rolling average of gene expression across genome ####
chrInfo <- read.csv('~/Documents/EPICC/Data/CopyNumber/centromerePositions_hg38.csv')
chrInfo <- chrInfo[which(chrInfo$chromosome %in% paste0('chr',c(1:22))),]
geneinfo <- as.data.frame(fread('~/Documents/EPICC/Data/Expression/compiledGeneInfo.txt'))
expinfo <- merge(geneinfo[,c('GeneID','Chr','Start','End','NewStart','NewEnd')],geneexp,by='GeneID',sort=F)

regionhg38 <- data.frame(NewStart=as.numeric(),NewEnd=as.numeric(),stringsAsFactors = FALSE)
rollmax <- 200
for(sam in matchsam) {
  if(sam == matchsam[1]) {
    means <- c()
    expinfo <- expinfo[order(expinfo$NewStart),]
    means <- c(means,rollmean(expinfo[,sam],rollmax))
    starts <- c();ends <- c()
    for(i in c(1:(nrow(expinfo)-(rollmax-1)))) {
      starts <- c(starts,as.numeric(expinfo[i,'NewStart']))
      ends <- c(ends,as.numeric(expinfo[i+(rollmax-1),'NewEnd']))
    }
    regionhg38[c((nrow(regionhg38)+1):(nrow(regionhg38)+length(starts))),] <- c(starts,ends)
    regionhg38[,sam] <- means
    regionhg38$NewStart <- as.numeric(regionhg38$NewStart);regionhg38$NewEnd <- as.numeric(regionhg38$NewEnd)
  } else {
    means <- c()
    expinfo <- expinfo[order(expinfo$NewStart),]
    means <- c(means,rollmean(expinfo[,sam],rollmax))
    regionhg38[,sam] <- means
  }
}

# 6. Plot rolling average of gene expression with CNA data  ####

segments <- read.table('~/Documents/EPICC/Data/CopyNumber/LowPass/Segments/C516_EPICC_C516_A1_G3_L1_500kb_GRCh38_multiregion_segmentation_calls.txt')
rowbins <- unlist(str_split(row.names(segments),':'));rowpos <- unlist(str_split(rowbins[seq.int(2,nrow(segments)*2,by=2)],'-'))
bins <- data.frame(chromosome=paste0('chr',rowbins[seq.int(1,nrow(segments)*2,by=2)]),start.pos=as.integer(rowpos[seq.int(1,nrow(segments)*2,by=2)]),end.pos=as.integer(rowpos[seq.int(2,nrow(segments)*2,by=2)]))

allpats <- gsub('(C\\d+)\\S+','\\1',matchsam);patients <- unique(allpats)

for(pat in patients) {
  patsam <- matchsam[grep(pat,matchsam)];numsam <- length(patsam)
  pdf(paste0('~/Documents/EPICC/Data/Expression/InferredCN/',pat,'.pdf'),width = 14,height=(1.5*numsam))
  layout(matrix(c(1:numsam),nrow=numsam,ncol=1))
  for(i in c(1:numsam)) {
    if(i==1) {
      par(mar=c(0,2,3,1),xpd=F,font=2);plotx <- FALSE
    } else if(i==numsam) {
      par(mar=c(3,2,0,1),xpd=F,font=2);plotx <- TRUE
    } else {
      par(mar=c(0,2,0,1),xpd=F,font=2);plotx <- FALSE
    }
    sam <- patsam[i];colreg <- as.character(regcols[gsub('C\\d+_(\\S)\\S+','\\1',sam)])
    if(sam %in% wgssam) { 
      segments <- as.data.frame(fread(paste0('~/Documents/EPICC/Data/CopyNumber/Segments/EPICC_',sam,'_D1_GRCh38_segments.txt')))
      segments$Status <- segments$CNt;segments$Length <- segments$end.pos-segments$start.pos
      segs <- segments[,c('chromosome','start.pos','end.pos','Status')]
    } else {
      segments <- read.table(paste0('~/Documents/EPICC/Data/CopyNumber/LowPass/Segments/',pat,'_EPICC_',sam,'_L1_500kb_GRCh38_multiregion_segmentation_calls.txt'))
      segs <- bins;segs$Status <- segments[,1]
    }
    
    plotcn(sam,segs,chrInfo,regionhg38,plotx,colreg)
  }
  dev.off()
}
