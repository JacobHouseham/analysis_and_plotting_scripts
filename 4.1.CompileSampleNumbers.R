# Script to compile number of samples from WGS and RNA-seq
# To do the following:
# 1) Fetch RNA PASS samples
# 2) Fetch WGS PASS samples
# 3) Combine into one table and output
# 4) Run DESeq2 and output normalised counts for symbol gene IDs

# Note: can't be run until after running 4.2.. scripts

# Main code ####

# 0. Prepare environment
library(dplyr);library(data.table)
library(stringr);'%ni%' <- Negate('%in%');library(sos)

# 1. Read in RNA QC ####
rnasam <- read.table(file='~/Documents/ThesisOther/ScriptsForThesis/Tables/ListRNAPass.txt')[,1]
tumrnapat <- table(gsub('(C\\d+)_\\S+','\\1',rnasam[grep('C\\d+_(A|B|C|D)',rnasam)]))
normrnapat <- table(gsub('(C\\d+)_\\S+','\\1',rnasam[grep('C\\d+_(E)',rnasam)]))

# 2. Read in WGS QC ####
wgssam <- readRDS('~/Documents/EPICC/Data/QTLs/nodesamples/wgs_used_in_phylo.rds')
wgssam <- gsub('EPICC_(C\\d+_\\S\\d+_\\S\\d+)\\S+','\\1',wgssam)
tumwgspat <- table(gsub('(C\\d+)_\\S+','\\1',wgssam[grep('C\\d+_(A|B|C|D)',wgssam)]))
tumwgspat['C306'] <- 2;tumwgspat <- tumwgspat[order(names(tumwgspat))]

vcflog <- read.table('~/Documents/EPICC/Data/Mutations/VCFtoTSV/vcf2tsv.log')[,'V2']
vcflog <- unique(vcflog[grep('Normal',vcflog)])
normals <- gsub('.+EPICC_C\\d+_(Z|E)\\S+','\\1',vcflog)
names(normals) <- gsub('.+EPICC_(C\\d+)_\\S+','\\1',vcflog)
normtype <- ifelse(normals=='E','Bulk','Blood')
normtype['C306'] <- 'Bulk'
normtype <- normtype[order(names(normtype))]

# 3. Combine QC and output as a csv file ####
combsam <- rnasam[grep('C\\d+_(A|B|C|D)',rnasam)][which(rnasam[grep('C\\d+_(A|B|C|D)',rnasam)] %in% wgssam)]
matchtum <- table(gsub('(C\\d+)_\\S+','\\1',combsam))
matchtum['C306'] <- 1;matchtum <- matchtum[order(names(matchtum))]

patinfo <- as.data.frame(readxl::read_xlsx('~/Documents/EPICC/Sample Info/epicc_pathology_slim.xlsx'))
patinfo <- patinfo[c(1,which(patinfo$Patient %in% names(normals))),c('Patient','Gender','Age')]
row.names(patinfo) <- c(1:nrow(patinfo))
patinfo$tumRNA <- 0;patinfo[patinfo$Patient %in% names(tumrnapat),'tumRNA'] <- tumrnapat
patinfo$tumDNA <- tumwgspat
patinfo$DNARNA <- 0;patinfo[patinfo$Patient %in% names(matchtum),'DNARNA'] <- matchtum
patinfo$normRNA <- 0;patinfo[patinfo$Patient %in% names(normrnapat),'normRNA'] <- normrnapat
patinfo$normDNA <- normtype
meds <- c(median(as.numeric(patinfo$tumRNA[which(patinfo$tumRNA!=0)])),
          median(as.numeric(patinfo$tumDNA)),
          median(as.numeric(patinfo$DNARNA[which(patinfo$DNARNA!=0)])))
patinfo[(nrow(patinfo)+1),] <- c(rep('',2),'Totals',as.numeric(colSums(patinfo[,c(4:7)])),'')
patinfo[(nrow(patinfo)+1),] <- c('','','Medians',meds,'','')

write.csv(patinfo,'~/Documents/Thesis/mainText/samplenumbers.csv',row.names = F,quote=F)




