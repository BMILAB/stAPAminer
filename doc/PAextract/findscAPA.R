# 加载R包
library(Rsamtools)
library(GenomicRanges)
library(stringi)
library(GenomicAlignments)
library(dplyr)
#导入BAM文件
args <- commandArgs(T)
bamfile <- args[1]
result_dir <- args[2]
sample_name <- args[3]
what <- c("rname","strand","pos","cigar","seq")
#bamfile <- "./temp/dedup.SRR6228889.bam"
param <- ScanBamParam(what=what)
gal1 <- GenomicAlignments::readGAlignments(bamfile, use.names=TRUE, param=param) 
s_1 <- (grepl("[0-9]*M[1-9]{2,}S",gal1@cigar) & as.vector(gal1@strand) == "+")
s_2 <- (grepl("[1-9]{2,}S[0-9]*M",gal1@cigar) & as.vector(gal1@strand) == "-")

bam1 <- gal1[s_1]
bam2 <- gal1[s_2]

bam1 <- bam1[grepl("(A{3,}[^A]{0,2})*A{6,}([^A]{0,2}A{3,})*.{0,2}?$",bam1@elementMetadata@listData$seq)]
bam2 <- bam2[grepl("^.{0,2}?(T{3,}[^T]{0,2})*T{6,}([^T]{0,2}T{3,})*",bam2@elementMetadata@listData$seq)]

final_bam1 <- data.frame(chr=as.vector(seqnames(bam1)),strand=as.vector(strand(bam1)),coord=end(bam1))
final_bam2 <- data.frame(chr=as.vector(seqnames(bam2)),strand=as.vector(strand(bam2)),coord=start(bam2))

bam <- rbind(final_bam1,final_bam2)
bam <- dplyr::group_by(bam,chr,strand,coord) %>% summarise(count = n())
print(dim(bam))
result_path <- paste0(result_dir,"/",sample_name,"_","scAPApeak.Rda")
save(bam,file=result_path)
