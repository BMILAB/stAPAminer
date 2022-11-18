library("derfinder")
library("GenomicRanges")
library("dplyr")
library("scAPA")
library("bumphunter")
args <- commandArgs(T)
# dir_name <- "/media/aa/KESU/scAPA/"
# sampledir <- "AML"
# filetermR <- "dedup.aml035_post_transplant_possorted_genome_bam.reverse.sorted.bam"
# filetermF <- "dedup.aml035_post_transplant_possorted_genome_bam.forward.sorted.bam"
dir_name <- args[1]
sampledir <- args[2]
filetermR <- args[3]
filetermF <- args[4]
result_file <- args[5]

splitPeak <- function(pos_count,chrom,pos){
  chr <- rep(chrom,length(pos_count))
  throshold <- quantile(pos_count,0.25)
  tempdata <- regionFinder(pos_count,chr = chr,pos = pos,cutoff = throshold)
  tempdata <- subset(tempdata,L >= 40)
  
  return(tempdata)
}

getDefAno <-function(dir_name,sampledir,fileterm,Strand){
  files <- rawFiles(dir_name,sampledirs = sampledir,fileterm = fileterm)
  print(files)
  chrs <- paste0('chr',c(as.character(1:19),'X','Y'))
  fullCov <- fullCoverage(files = files,chrs = chrs)
  regionMat <- regionMatrix(fullCov = fullCov,L = 96,cutoff = 0)
  DefPeak <- GenomicRanges::GRanges()
  for(i in names(regionMat)){
    index <- which(width(regionMat[[i]]$regions) > 1000)
    tempDataFram <- data.frame()
    for(j in index){
      #browser()
      rawstart <- start(regionMat[[i]]$regions[j])-1
      temp <- regionMat[[i]]$bpCoverage[[j]]
      temp$pos <- c(1:dim(temp)[1])
      # 进行第一次切分，如果还存在peak大于1000,则继续切分
      first.split.peaks <- splitPeak(temp$value,i,temp$pos)
      pre.split.peaks <- first.split.peaks
      split.peaks <- data.frame()
      iter <- 0
      while(TRUE){
        if (iter >= 3){
          break()
        }
        Ind <- which(pre.split.peaks$L > 1000)
        cur.split.peaks <- data.frame()
        if(length(Ind) > 0){
          # 将小于1000的peak 存储到split.peaks
          split.peaks <- rbind(split.peaks,pre.split.peaks[-Ind,])
          for (k in Ind){
            pos.temp <- temp[pre.split.peaks$start[k]:pre.split.peaks$end[k],]
            cur.split.peaks <- splitPeak(pos.temp$value,i,pos.temp$pos)
          }
        }else{
          split.peaks <- rbind(split.peaks,pre.split.peaks)
          break()
        }
        pre.split.peaks <- cur.split.peaks
        iter <- iter + 1
        
      }
      split.peaks$start <- split.peaks$start + rawstart
      split.peaks$end <- split.peaks$end + rawstart
      tempDataFram <- rbind(tempDataFram,split.peaks)
    }
    DefPeak <- c(DefPeak,regionMat[[i]]$regions[-index])
    if(length(tempDataFram) != 0){
      DefPeak <- c(DefPeak,regioneR::toGRanges(tempDataFram))
    }
  }
  DefPeak <- DefPeak[-which(DefPeak$area == width(DefPeak))]
  DefPeak <- DefPeak[width(DefPeak) >= 40]
  strand(DefPeak) <- Strand
  if(Strand == "+"){
    DefPeak <- as.data.frame(DefPeak,row.names = 1:length(DefPeak))
    DefPeak <- dplyr::mutate(DefPeak,chr=seqnames,coord=end)
  }else{
    DefPeak <- as.data.frame(DefPeak,row.names = 1:length(DefPeak))
    DefPeak <- dplyr::mutate(DefPeak,chr=seqnames,coord=start)
  }
  
  return(DefPeak)
}

DeFF <- getDefAno(dir_name = dir_name,sampledir = sampledir,fileterm = filetermF,"+")
DeFR <- getDefAno(dir_name = dir_name,sampledir = sampledir,fileterm = filetermR,"-")
DeF <- rbind(DeFF,DeFR)
#DeF <- subset(DeF,width >= 98 & width <= 1000)
DeF$PeakID <- paste0("peak","_",1:length(DeF$chr))
print("识别的peak数目")
print(dim(DeF))
DeF.saf <- DeF[,c("PeakID","chr","start","end","strand")]
write.bed(DeF.saf,result_file)
