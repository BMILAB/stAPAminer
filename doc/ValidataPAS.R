## This script shows the code for the validation of poly(A) sites and identification of poly(A) signals, corresponding to Figures 2B and 2C in the manuscript.

library(movAPA)
library(stAPAminer)
library(BSgenome.Mmusculus.UCSC.mm10)

#### Load data ----

## download the data from http://www.bmibig.cn/mnt/stAPAminer/

## reference PolyAsite 2.0
mm10polyAsite2 <- readRDS("mm10polyAsite2.rds")

## poly(A) sites from scRNA-seq of MOB data (SC-MOB)
SCMOB <- readRDS("SCMOB.rds")

## poly(A) sites from ST-MOB Replicate 11
STMOB <- readRDS("STMOB.rds")

#### Plot single nucleotide frequency (Figure 2B) ----
bsgenome <- BSgenome.Mmusculus.UCSC.mm10

faFromPACds(mm10polyAsite2, bsgenome, what='updn', fapre='mm10polyAsite2', up=-100, dn=100)
mm10polyAsite2_plot <- plotATCG(fafile="mm10polyAsite2.fa",refPos=101,title="PolyASite 2.0")

faFromPACds(SCMOB, bsgenome, what='updn', fapre='SCMOB', up=-100, dn=100)
SCMOB_plot <- plotATCG(fafile="SCMOB.fa",refPos=101,title="SCMOB scAPAtrap")

faFromPACds(SCMOB_Sierra, bsgenome, what='updn', fapre='SCMOB_Sierra', up=-100, dn=100)
SCMOB.Sierra_plot <- plotATCG(fafile="SCMOB_Sierra.fa",refPos=101,title="SCMOB Sierra",annotate.y=0.4)

STMOB_3UTR = get3UTRAPAds(APA,sortPA=TRUE)
STMOB_3UTR@counts$tag <- 1
STMOB_3UTR@counts <- as.data.frame(STMOB_3UTR@counts[,ncol(STMOB_3UTR@counts)])
colnames(STMOB_3UTR@counts) <- "tag"
rownames(STMOB_3UTR@counts) <- rownames(STMOB_3UTR@anno) 

faFromPACds(STMOB_3UTR, bsgenome, what='updn', fapre='STMOB', up=-100, dn=100)
STMOB_plot <- plotATCG(fafile="STMOB.fa",refPos=101,title="STMOB")

library(patchwork)
comATCGplot <- mm10polyAsite2_plot+ SCMOB_plot+ SCMOB.Sierra_plot + STMOB_plot
library(eoffice)
topptx(comATCGplot,filename = "comATCGplot.pptx", width = 10,height = 6)

#### validate poly(A) sites (Figure 2C) ----

# adjusting the coordinate (0-based in mm10polyAsite2; 1-based in STMOB)
STMOB_adj <- STMOB
STMOB_adj@anno[STMOB_adj@anno$strand == "+",]$coord <- STMOB_adj@anno[STMOB_adj@anno$strand == "+",]$coord - 1
STMOB_adj@anno[STMOB_adj@anno$strand == "-",]$coord <- STMOB_adj@anno[STMOB_adj@anno$strand == "-",]$coord + 1

Compare1 <- annotateByKnownPAC(STMOB_adj, mm10polyAsite2, labels = "Mm", d = 100)
Compare2 <- annotateByKnownPAC(STMOB, SCMOB, labels = "Mm", d = 100)

length(na.omit(Compare1@anno$Mm_dist))
length(na.omit(Compare2@anno$Mm_dist))
sum(Compare1@anno$Mm_ovp)

sum((Compare1@anno$Mm_ovp + Compare2@anno$Mm_ovp) == 0 )

PolyASite.dist <- as.data.frame(table(na.omit(Compare1@anno$Mm_dist)))
scMOB.dist <- as.data.frame(table(na.omit(Compare2@anno$Mm_dist)))

PolyASite.dist$Freq <- PolyASite.dist$Freq/sum(PolyASite.dist$Freq)
scMOB.dist$Freq <- scMOB.dist$Freq/sum(scMOB.dist$Freq)

PolyASite.dist$method <- rep("PolyAsite 2.0", length(PolyASite.dist$Freq))
scMOB.dist$method <- rep("scMOB", length(scMOB.dist$Freq))

all_dist <- rbind(PolyASite.dist, scMOB.dist)
all_dist$Var1 <- as.numeric(as.character(all_dist$Var1))

all_dist$method <- factor(all_dist$method, levels = c("PolyAsite 2.0", "scMOB"))

ggplot(all_dist, aes(x=Var1, y=Freq,fill=method,colour=method)) +
  geom_line(size=1.5)+labs(title="")+scale_fill_brewer(palette = "Spectral")+
  xlab("Distance(bp)") + ylab("Density")+theme_bw() +
  theme(panel.grid =element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour='black'),
        plot.title = element_text(size = 14, face =  "bold",hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))

topptx(filename = "all_dist.pptx", width = 4,height = 3)


annotateByKnownPAC <- function (pacds, knownPACdss, labels, d=50, verbose=TRUE) {
  
  #library(GenomicRanges, verbose = FALSE)
  
  
  gr1 <- with(pacds@anno, GRanges(seqnames = chr,
                                  ranges =IRanges(start=coord,
                                                  end=coord),
                                  strand = strand) )
  
  if (class(knownPACdss)=='PACdataset') {
    if (length(labels)!=1) stop("knownPACdss is PACdataset, but more than one elements in labels\n")
    knownPACdss=list(knownPACdss)
    names(knownPACdss)=labels
  } else if (!is.list(knownPACdss)) {
    stop("knownPACdss is neither a PACdataset nor a list!\n")
  } else if (is.list(knownPACdss)) {
    if (length(labels)!=length(knownPACdss)) stop("Not equal element# in knownPACdss and labels\n")
    names(knownPACdss)=labels
  }
  
  for (i in 1:length(knownPACdss)) {
    gr2 <- with(knownPACdss[[i]]@anno, GRanges(seqnames = chr,
                                               ranges =IRanges(start=coord,
                                                               end=coord),
                                               strand = strand) )
    lbl=names(knownPACdss)[i]
    
    ov = findOverlaps(gr1, gr2,
                      maxgap=d-1, minoverlap=0L,
                      type=c("any"), select='first',
                      ignore.strand=FALSE)
    isOvp=ov
    im=which(!is.na(ov))
    isOvp[-im]=0
    isOvp[im]=1
    minD=rep(NA, length(ov))
    minD[im]=knownPACdss[[i]]@anno$coord[ov[im]]-pacds@anno$coord[im]
    
    pacds@anno[, paste0(lbl,'_ovp')]=isOvp
    pacds@anno[, paste0(lbl,'_dist')]=minD
    
    if (verbose) {
      cat('query:',length(gr1),'\n')
      cat('known:',length(gr2),'(',lbl,')\n')
      cat(length(im), 'query PACs overlapping known PACs\n')
    }
  }
  return(pacds)
}


plotATCG <- function (fafile ,refPos = NULL, title="",
                      start = NA, end = NA, annotate.x=50,annotate.y=0.5) 
{
  
  seq = readDNAStringSet(fafile, format = "fasta")
  nseq = length(seq)
  if (!is.na(start) | !is.na(end)) {
    seq = subseq(seq, start = start, end = end)
    refPos = NULL
  }
  din <- consensusMatrix(seq, as.prob = TRUE, baseOnly = TRUE)
  din <- as.data.frame(t(din))
  din$other <- NULL
  din$pos <- c(1:nrow(din))
  if (!is.null(refPos)) {
    din$pos = din$pos - refPos
  }
  
  din <- reshape2::melt(din, id.vars = "pos", variable.name = "base", 
                        value.name = "freq")
  annotates <- paste0("-1 : A= ",round(din[din$pos == -1 & din$base == "A",]$freq,2),"\n",
                      "0 : A= ",round(din[din$pos == 0 & din$base == "A",]$freq,2),"\n",
                      "1 : A= ",round(din[din$pos == 1 & din$base == "A",]$freq,2))
  ATCGpicture <- ggplot(din, aes(x = pos, y = freq,colour = base)) + geom_line() + 
    labs(title=title)+
    xlab("Position") + 
    ylab("Base component") + 
    theme_bw()+
    annotate("text", x=annotate.x, y=annotate.y, label=annotates)+
    theme(panel.grid =element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour='black'),
          plot.title = element_text(size = 14, face =  "bold",hjust = 0.5),
          text = element_text(size = 12),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 12))
  
  
  return(ATCGpicture)
}

