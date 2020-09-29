# Script to assess the clonal mixing of phylogenetic trees
# To do the following:
# 1) Load phylogenetic data
# 2) Assess clonal mixing of phylogenetic trees

# Functions ####
# Takes a distance matrix and returns the sum of
# dists between samples of the same group
getsumdists <- function(dm,sams) {
  regs <- gsub('^(\\S)\\d+_\\S+','\\1',sams)
  dists <- c();for(i in c(1:length(sams))) { for(j in c(1:length(sams))) { if(i<j) {if(regs[i]==regs[j]) { dists <- c(dists,dm[i,j]) } } }}
  return(sum(dists))
}

# Takes a distance matrix and returns a vector of all 
# dists between samples of the same group
getalldists <- function(dm,sams,alldist,nodedist,pat) {
  regs <- gsub('^(\\S)\\d+_\\S+','\\1',sams)
  dists <- c();nams <- c();tots <- c();for(i in c(1:length(sams))) { for(j in c(1:length(sams))) { if(i<j) {if(regs[i]==regs[j]) { dists <- c(dists,dm[i,j]);nams <- c(nams,paste0(pat,':',sams[i],'vs',sams[j]));tots <- c(tots,nodedist[getMRCA(MPtree,tip=c("Root",sams[i],sams[j])),getMRCA(MPtree,tip=c(sams[i],sams[j]))]) } } }}
  tmpnam <- names(alldist)
  alldist <- c(alldist,dists)
  names(alldist) <- c(tmpnam,nams)
  return(alldist)
}

# Main code ####

# 0. Prepare environment
setwd("~/Documents/EPICC/Data/Mutations");library(dplyr);library(data.table);library(phangorn)

# 1. Load phylogenetic data ####
DistList <- readRDS('~/Documents/ThesisOther/ScriptsForThesis/RData/DNAPhyDists.rds')
load(file='~/Documents/ThesisOther/ScriptsForThesis/RData/DNAphylogenies.Rdata')

# Get the patients with no mixing according to just looking at the phylogenetic trees
nomixing <- c('C516','C518','C525','C528','C532','C536','C537','C538','C539','C543','C544','C548','C549','C550','C552','C554','C560','C561')

patients <- names(DistList)

# 2. Analyse and plot clonal mixing in phylogenetic trees of each patient ####

pdf('~/Documents/ThesisOther/ScriptsForThesis/Plots/mixing_permutations.pdf')
layout(matrix(c(1:6),nrow=3,ncol=2,byrow=T))
mixps <- c();usable <- c();nsam <- c()
for(p in c(1:length(patients))) {
  print(p)
  pat <- patients[p]
  MPtree <- MPtreesList[[pat]]
  
  distmat <- cophenetic.phylo(MPtree);nodedists <- dist.nodes(MPtree)
  samples <- colnames(distmat)[grep('^(A|B|C|D)',colnames(distmat))]
  regpat <- gsub('^(\\S)\\d+_\\S+','\\1',samples)
  sumdist <- getsumdists(distmat[samples,samples],samples)
  
  if(length(table(regpat))>1 & max(table(regpat))>1) {
    usable <- c(usable,pat)
    nsam <- c(nsam,length(samples))
    
    permdists <- c()
    for(i in c(1:10000)) {
      tmpmat <- distmat[samples,samples]
      colnames(tmpmat) <- row.names(tmpmat) <- sample(samples)
      permdists <- c(permdists,getsumdists(tmpmat,colnames(tmpmat)))
    } 
    pval <- (length(which(permdists<=sumdist))+1)/(10000+1)
    mixps <- c(mixps,pval)
    if(sumdist<min(permdists)) { 
      hist(permdists,breaks=30,col='skyblue3',xlim=c(sumdist-1000,max(permdists)),border=NA,main=pat,xlab=paste0('p=',pval));abline(v=sumdist,lty=2,lwd=2,col='firebrick3')
    } else {
      hist(permdists,breaks=30,col='skyblue3',border=NA,main=pat,xlab=paste0('p=',pval));abline(v=sumdist,lty=2,lwd=2,col='firebrick3')
    }
  }
}
dev.off()
names(mixps) <- names(nsam) <- usable
signomix <- names(mixps)[which(p.adjust(mixps,method='fdr')<0.05)]

