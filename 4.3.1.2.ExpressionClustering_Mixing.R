# Script to assess the degree of clonal mixing in clustering by gene expression
# To do the following:
# 1) Load data and assess which tumour are multi-regional
# 2) Assess clonal mixing with different filtered gene lists

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
  dists <- c();nams <- c();for(i in c(1:length(sams))) { for(j in c(1:length(sams))) { if(i<j) {if(regs[i]==regs[j]) { dists <- c(dists,dm[i,j]);nams <- c(nams,paste0(pat,':',sams[i],'vs',sams[j])) } } }}
  tmpnam <- names(alldist)
  alldist <- c(alldist,dists)
  names(alldist) <- c(tmpnam,nams)
  return(alldist)
}

# Main code ####

# 0. Prepare environment
setwd("~/Documents/EPICC/Data/Expression");library(dplyr);library(data.table);library(phangorn)

# 1. Assess which tumours are multi-sample ####
# Load data from 4.3.1.1
load(file='~/Documents/ThesisOther/ScriptsForThesis/RData/expression_clustering_distances.Rdata')

# Filter for tumours with more than one region, and at least on of those regions has multiple samples
usable <- c();nsam <- c()
patients <- names(distlist)
for(p in c(1:length(patients))) {
  pat <- patients[p]
  distmat <- as.matrix(distlist[[pat]])
  samples <- colnames(distmat)[grep('^(A|B|C|D)',colnames(distmat))]
  regpat <- gsub('^(\\S)\\d+_\\S+','\\1',samples)
  if(length(table(regpat))>1 & max(table(regpat))>1) {
    usable <- c(usable,pat)
    nsam <- c(nsam,length(samples))
  }
}

# 2. For the different filter types (All, Epithelial or Stromal) assess clonal mixing ####
resdist <- matrix(0L,nrow=length(usable),ncol=3);colnames(resdist) <- c('All','Epithelial','Stromal')
row.names(resdist) <- usable;sigdist <- resdist
types <- c('','_epi','_str')
for(t in c(1:3)) {
  typ <- types[t]
  pdf(paste0('~/Documents/ThesisOther/ScriptsForThesis/Plots/mixing_permutations_expression',typ,'.pdf'))
  layout(matrix(c(1:6),nrow=3,ncol=2,byrow=T))
  for(p in c(1:length(usable))) {
    pat <- usable[p]
    distmat <- as.matrix(get(paste0('distlist',typ))[[pat]])
    samples <- colnames(distmat)[grep('^(A|B|C|D)',colnames(distmat))]
    sumdist <- getsumdists(distmat,samples)
    resdist[pat,t] <- sumdist#/mean(colSums(distmat))
    
    permdists <- c()
    for(i in c(1:10000)) {
      tmpmat <- distmat
      colnames(tmpmat) <- row.names(tmpmat) <- sample(samples)
      permdists <- c(permdists,getsumdists(tmpmat,colnames(tmpmat)))#/mean(colSums(tmpmat)))
    } 
    sigdist[pat,t] <- (length(which(permdists<=resdist[pat,t]))+1)/(10000+1)
    if(resdist[pat,t]<min(permdists)) {
      hist(permdists,breaks=30,xlim=c(resdist[pat,t],max(permdists)),col='skyblue3',border=NA,main=pat,xlab=paste0('p=',sigdist[pat,t]));abline(v=resdist[pat,t],lty=2,lwd=2,col='firebrick3')
    } else {
      hist(permdists,breaks=30,col='skyblue3',border=NA,main=pat,xlab=paste0('p=',sigdist[pat,t]));abline(v=resdist[pat,t],lty=2,lwd=2,col='firebrick3')
    }
  }
  dev.off()
}

