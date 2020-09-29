# Script to plot phylogenetic trees for each patient
# To do the following:
# 1) Manually colour the nodes of each tree
# 2) Plot the trees with colour coding of leaves and nodes

# Main code ####

# 0. Prepare environment
regcol <- c(A='#E31A1C',B='#377DB8',C='#4DAE49',D='#904A9A',Root='#808080')
setwd("~/Documents/EPICC/Data/Mutations");library(dplyr);library(data.table);library(phangorn)
load(file='~/Documents/ThesisOther/ScriptsForThesis/RData/DNAphylogenies.Rdata')
allpats <- gsub('(C\\d+)_\\S+','\\1',multisam);patients <- unique(allpats)

# 1. Manual colour coding of branches to fit with clustering of regions ####
C516para <- c(2,2,2,1,3,3,4,3,5,5,5,5)
C518para <- c(4,4,2,2,1,1,1,2,5,4,5,5)
C522para <- c(1,1,1,1,3,4,5,3,4,4,4,5,5,1,5,5)
C524para <- c(3,3,3,3,3,4,5,4,5,3,2,2,2,2,5,3,5,5,5,5)
C525para <- c(4,4,3,3,2,2,2,2,2,3,5,4,5,5)
C528para <- c(3,3,3,2,2,5,5,1,5,1,5,5)
C530para <- c(4,4,4,4,4,4,3,2,1,1,1,1,1,1,1,5,5,4,5,2,5,5)
C531para <- c(2,2,2,2,1,1,3,3,3,3,3,3,3,1,4,4,4,5,5,2,5,5)
C532para <- c(3,3,3,3,1,1,1,1,1,3,5,2,5,4,5,5)
C536para <- c(2,2,4,4,4,4,4,4,4,2,5,5)
C537para <- c(1,1,3,3,3,3,3,4,5,2,5,1,5,5)
C538para <- c(4,4,2,2,3,3,3,3,3,2,5,2,5,4,5,4,5,1,5,5)
C539para <- c(1,1,1,1,1,2,5,2,4,4,3,3,3,4,5,5,5,5)
C542para <- c(3,3,3,3,3,2,2,2,4,1,5,1,5,2,5,5,5,5)
C543para <- c(1,1,2,2,2,2,2,2,2,1,5,5)
C544para <- c(2,2,4,4,1,1,1,1,1,3,5,3,5,4,5,2,5,5)
C548para <- c(2,2,1,1,1,2,3,3,4,4,4,3,5,5,5,5)
C549para <- c(2,2,2,2,2,3,5,1,5,5)
C550para <- c(3,3,3,4,5,4,5,5)
C551para <- c(3,3,3,3,1,1,4,1,4,1,5,5,5,3,5,2,5,2,5,5)
C552para <- c(3,3,4,4,4,3,5,5)
C554para <- c(3,3,3,3,4,4,4,3,2,2,2,2,2,5,5,5)
C555para <- c(4,4,4,4,4,5)
C559para <- c(2,2,2,2,1,1,4,4,4,4,3,4,5,3,5,3,5,5,5,2,5,5)
C560para <- c(1,1,3,3,3,1,5,2,5,2,5,5)
C561para <- c(2,2,2,2,2,3,5,3,1,1,4,4,4,4,4,1,5,5,5,5)
C562para <- c(1,1,1,1,1,5)

# 2. Plot phylogenetic trees per patient, colouring leaves and nodes by region of origin ####

pdf('~/Documents/Thesis/figures/EPICC_DNA_Phylogenies.pdf',width = 18,height=10)
layout(matrix(c(1:30),nrow=4,ncol=7,byrow = T))
par(mar=c(3.1,1.1,2.1,1.1),xpd=T)
for(p in c(1:length(patients))) {
  pat <- patients[p]
  patsam <- multisam[grep(pat,multisam)]
  
  regions <- gsub('C\\d+_(\\S).+','\\1',patsam)
  MPtree <- MPtreesList[[pat]];BStrees <- BStreesList[[pat]];hi <- his[p]
  
  suppressWarnings(plotBS(MPtree,BStrees,type.tree='phylogram',bs.adj=c(1.1,1.2),label.offset=sum(MPtree$edge.length)*0.01,p=0,type="phylogram",edge.color=regcol[get(paste0(pat,'para'))],edge.width=1.5,font=2,cex=1,tip.color=regcol[c(regions,'Root')]))
  add.scale.bar(0,0.5,lwd=1.5,font=2);mtext(paste0('HI: ',hi),side=1,font=1,cex=0.75,ask=T)
  mtext(side=3,text=pat,line=0.5,cex.main=1,font=2)
}
dev.off()

