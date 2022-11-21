#Note： The extraction process of poly(A) refers to scAPAtrap : http://www.bmibig.cn/mnt/scAPAtrap/Tutorial/scAPAtrap_compare.html
## Call scAPAtrap script  
Usage: ./Def.sh ./scAPAtrapDef.R ./findscAPA.R [sample_name] bam-file   
Example: ./Def.sh ./scAPAtrapDef.R ./findscAPA.R st_mob ./st_mob.bam  
As a result, three files will be generated:   
1. DefExp/st_mob.counts.tsv.gz   
2. DefSAF/Def_st_mob.saf    
3. scPA/ st_mob_scAPApeak.Rda.   
Used in the next step to generate movAPA objects.
## Generate PACds object
library(scAPAtrap)
library(movAPA)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
sample = 'st_mob'
load(paste0('./scPA/', sample, '_scAPApeak.Rda'))
head(bam)
##If the chr column in bam does not have the prefix of chr, add it.
bam$chr<-paste0("chr",bam$chr)
barcodefile <- './barcodes.tsv'
barcode <- read.delim2(barcodefile, header = F)
barcode <- gsub('-[0-9]','',barcode$V1)
countsfile <- paste0('./DefExp/', sample, '.counts.tsv.gz')
countsfile
peaksfile <- paste0('./DefSAF/Def_', sample, '.saf')
peaksfile
expma<- generatescExpMa(countsfile, peaksfile, barcode, bam, min. cells = 2, min. count = 0, d=100)
coldat <- data.frame(group = colnames(expma)[7:ncol(expma)], row.names = colnames(expma)[7:ncol(expma)])
scPACds <- movAPA::readPACds(pacFile = expma, colDataFile = coldata)
## Annotate PA (modify according to your own annotation file)
txdbmmu.10 <- TxDb.Mmusculus.UCSC.mm10.knownGene
txdbmmu.10 <- parseGenomeAnnotation(txdbmmu.10)
scPACds<- annotatePAC(scPACds, txdbmmu.10)
## Apply 1000bp to the 3utr extension, and classify the PA sites located in it as 3utrPA sites
scPACds <- ext3UTRPACds(scPACds, ext3UTRlen = 1000)

