# Identification and quantification of APA sites from ST data

The process of extracting APA sites from ST data is similar to that of scRNA-seq data. The raw ST data are double stranded, including read 1 and read 2. First, barcodes and UMIs were extracted from read 1 and then umi_tools was adopted to append the barcode and UMI information to the sequence header of read 2 to generate a new read 2 FASTQ file. Reads from the read 2 file were aligned to the reference genome using STAR and only uniquely mapped reads were retained. After mapping, PCR duplicates were removed by the dedup function in umi_tools and one read per UMI for each genomic coordinate was remained. Then scAPAtrap (other tools like Sierra can also be used) was used to identify poly(A) sites from these mapped reads. Poly(A) sites were annotated with information, such as genomic region and gene, with the movAPA package using the R annotation package ‘txdb.musculus.ucsc.mm10.knowngene’. 3′ UTRs of annotated genes were extended by 1000 bp to recruit intergenic poly(A) sites which may be originated from authentic 3′ UTRs.

## Using scAPAtrap for identifying poly(A) sites from ST data

Usage: `./Def.sh ./scAPAtrapDef.R ./findscAPA.R [sample_name] bam-file`
Example: `./Def.sh ./scAPAtrapDef.R ./findscAPA.R st_mob ./st_mob.bam`

The pipeline of scAPAtrap generates three files:   
1. DefExp/st_mob.counts.tsv.gz   
2. DefSAF/Def_st_mob.saf    
3. scPA/ st_mob_scAPApeak.Rda.   

Please refer to the [scAPAtrap pipeline](http://www.bmibig.cn/mnt/scAPAtrap/Tutorial/scAPAtrap_compare.html) for identifying poly(A) sites from BAM files. 

## Using movAPA for annotating poly(A) sites

```
library(scAPAtrap)
library(movAPA)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

## load results from scAPAtrap
sample = 'st_mob'
load(paste0('./scPA/', sample, '_scAPApeak.Rda'))

## If the chr column in bam does not have the prefix of chr, add it.
head(bam)
bam$chr<-paste0("chr",bam$chr)
barcodefile <- './barcodes.tsv'
barcode <- read.delim2(barcodefile, header = F)
barcode <- gsub('-[0-9]','',barcode$V1)
countsfile <- paste0('./DefExp/', sample, '.counts.tsv.gz')
countsfile
peaksfile <- paste0('./DefSAF/Def_', sample, '.saf')

expma<- generatescExpMa(countsfile, peaksfile, barcode, bam, min. cells = 2, min. count = 0, d=100)
coldat <- data.frame(group = colnames(expma)[7:ncol(expma)], row.names = colnames(expma)[7:ncol(expma)])

## make the PACdataset for movAPA
scPACds <- movAPA::readPACds(pacFile = expma, colDataFile = coldata)

## Annotate poly(A) sites using the TXDB annotation (please modify it according to your own species)
txdbmmu.10 <- TxDb.Mmusculus.UCSC.mm10.knownGene
txdbmmu.10 <- parseGenomeAnnotation(txdbmmu.10)
scPACds<- annotatePAC(scPACds, txdbmmu.10)

## Apply 1000bp to the 3utr extension, and classify the PA sites located in it as 3utrPA sites
scPACds <- ext3UTRPACds(scPACds, ext3UTRlen = 1000)
```
