# Script to get phylogenetic trees inferred by SNVs for each tumour
# To do the following:
# 1) Load sample information
# 2) Perform phylogenetic analysis
# 3) Save results

# Main code ####

# 0. Prepare environment
setwd("~/Documents/EPICC/Data/Mutations")
library(dplyr);library(data.table);library(phangorn)

# 1. Load sample data ####
samples <- readRDS('~/Documents/EPICC/Data/QTLs/nodesamples/wgs_used_in_phylo.rds')
samples <- gsub('EPICC_(C\\d+_\\S\\d+_\\S\\d+)\\S+','\\1',samples)
samples <- samples[grep('C\\d+_(A|B|C|D)\\d+_G',samples)]
allpats <- gsub('(C\\d+).+','\\1',samples);patients <- unique(allpats)

# 2. For each patient, perform SNV phylogenetics ####
MPtreesList <- list();BStreesList <- list();DistList <- list()
his <- c();multisam <- c();snvsam <- c();snvpat <- c()
for(p in c(1:length(patients))) {
  pat <- patients[p]
  patsam <- samples[which(allpats==pat)]
  if(length(patsam)<=2) {
    write(paste0('Not enough samples for phylogenetics for patient ',pat,' (',length(patsam),' sample(s))'), stdout())
  } else {
    write(paste0('Loading SNVs and outputting phylogeny of patient ',pat,' (',length(patsam),' samples)'), stdout())
    multisam <- c(multisam,patsam)
    regions <- gsub('C\\d+_(\\S).+','\\1',patsam)
    
    dnalist <- list()
    for(i in c(1:length(patsam))) {
      sam <- patsam[i]
      # Load pre-comiled SNV data for each sample
      wgs <- as.data.frame(fread(paste0('~/Documents/EPICC/Data/Mutations/VCFtoTSV/only_snvs/EPICC_',sam,'_D1.txt'),sep='\t',header=T))
      wgs$Locus <- paste0(wgs$Chr,':',wgs$Pos,'_',wgs$Ref,'/',wgs$Alt)
      dnalist[[i]] <- wgs
      snvsam <- c(snvsam,nrow(wgs))
      names(snvsam) <- c(names(snvsam)[c(1:length(snvsam)-1)],sam)
    }
    
    dnamerge <- dnalist %>%
     Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Locus"), .)
    dnamerge <- dnamerge[,grep('Ref',names(dnamerge))] 
    for(i in c(1:length(patsam))) { dnamerge[,patsam[i]] <- ifelse(is.na(dnamerge[,i]), 0, 1) }
    
    DNAbin <- dnamerge[,patsam];DNAbin$Root <- rep(0,nrow(DNAbin))
    colnames(DNAbin) <- c(gsub('C\\d+_(\\S+)','\\1',patsam),'Root')
    DNAphy <- as.phyDat(DNAbin,type="USER",levels=c(0,1))
    DNAtreeMP <- pratchet(DNAphy);DNAtreeMP <- root(DNAtreeMP,'Root',resolve.root=T)
    DNAtreeMP <- acctran(DNAtreeMP, DNAphy);set.seed(123)
    DNABStrees <- bootstrap.phyDat(DNAphy, pratchet, bs = 100,jumble = F)
    
    # Calculate Consistency Index and Homoplasy Index
    ci <- CI(DNAtreeMP,DNAphy) # ci=1 means no homoplasy
    hi <- signif(1-ci,digits=2) # so hi=0 means no homoplasy
    his <- c(his,hi)
    
    BStreesList[[pat]] <- DNABStrees
    MPtreesList[[pat]] <- DNAtreeMP
    DistList[[pat]] <- DNAphy
  }
}

# 3. Save phylogeny results including distance matrices ####
save(multisam,MPtreesList,BStreesList,his,file='~/Documents/ThesisOther/ScriptsForThesis/RData/DNAphylogenies.Rdata')
saveRDS(DistList,file='~/Documents/ThesisOther/ScriptsForThesis/RData/DNAPhyDists.rds')


