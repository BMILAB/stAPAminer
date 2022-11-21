# stAPAminer
Mining Spatial Patterns of Alternative Polyadenylation for Spatially Resolved Transcriptomic Studies

## About
Alternative polyadenylation (APA) contributes to transcriptome complexity and gene expression regulation, which has been implicated in various cellular processes and diseases. Single-cell RNA-seq (scRNA-seq) has led to the profile of APA at the single-cell level, however, the spatial information of cells is not preserved in scRNA-seq. Alternatively, spatial transcriptomics (ST) technologies provide opportunities to decipher the spatial context of the transcriptomic landscape within single cells and/or across tissue sections. Pioneering studies on ST have unveiled potential spatially variable genes and/or splice isoforms, however, the pattern of APA usages in spatial contexts remains unappreciated. 
Here, we developed a toolkit called stAPAminer for mining spatial patterns of APA from spatial barcoded ST data. APA sites were identified and quantified from the ST data. Particularly, an imputation model based on k-nearest neighbors algorithm was designed for recovering APA signals. Then APA genes with spatial patterns of APA usage variation were identified. By analyzing the well-established ST data of mouse olfactory bulb (MOB), we present a detailed view of spatial APA usage across morphological layers of MOB with stAPAminer. We complied a comprehensive list of genes with spatial APA dynamics and obtained several major spatial expression patterns representing spatial APA dynamics in different morphological layers. Extending this analysis to two additional replicates of the MOB ST data, we found that spatial APA patterns of many genes are reproducible among replicates.

* Preparing the input for stAPAminer

The input of stAPAminer is a poly(A) site matrix with each row being a poly(A) site and each column being a spot. APA sites can only be identified from the spatial barcoded ST data for now, using existing scRNA-seq tools like scAPAtrap or Sierra. Please refer to this [document](https://github.com/BMILAB/stAPAminer/blob/main/doc/PAextract/README.md) for detailed pipeline of identifying APA sites from spatial barcoded ST data.

* The stAPAminer package consists of three main modules.
<img src="pic/process.jpg" width="100%" />

A. Quantification and imputation of APA usages

B. Verifying the correctness of the imputation method

C. Identifying genes with differential APA usage between morphological layers and genes with spatial patterns of APA usage variation

## Getting started
### Mandatory
* R (>=3.6.0). [R 3.6.3](https://www.r-project.org/) is recommended.

### Required R Packages
* movAPA, org.Mm.eg.db, Seurat, ggplot2, fdm2id, ClusterR, cluster, clusterSim, SPARK, limma, edgeR, ggpubr

### Installation
* Install the R package using the following commands on the R console:

```
install.packages("devtools")
require(devtools)
install_github("BMILAB/stAPAminer")
library(stAPAminer)
##or you can download ZIP, and then unzip
install.packages("you unzip file path", repos = NULL, type = "source")
```

## Application examples
### Quantification and imputation of APA usages
```
library(stAPAminer)
## the gene-spot count matrix
count<-read.csv(system.file("extdata", "count.csv", package = "stAPAminer"))
## individual spots and coordinates
position<-read.table(system.file("extdata", "position.txt", package = "stAPAminer"))

## the APA-spot count matrix stored as a movAPA PACdataset
load(system.file("extdata", "APA.RDA", package = "stAPAminer"))

## calculating APA ratios
RUD <- computeAPAIndex(APA,rownames(position))

## imputation
RUD <- imputeAPAIndex(RUD_RAW,count,k=10)
```
### Verifying the correctness of the imputation method
```
## clustering
stObj <- createStAPAminerObject(RUD,position)
stObj <- makeStCluster(stObj,resolution = 0.6,dims = 1:10,k=30,nfeatures = 4500)
findLabels(stObj)
stObj <- RenameCluster(stObj,c("GCL","GL","ONL","OPL","MCL"))
drawClusterPlot(stObj,color=distinctColorPalette(5),islegend = T,size = 3)

##Verify by correlation
pearson_RUD <- computePearsonIndex(stObj,"RUD")
index <- computeIndex(stObj)
drawVolcano(stObj,"GCL","GL")
```
### Identification of DEAPA, LSAPA and SVAPA 
```
stObj <- findContrastAPA(stObj)
stObj <- findAllMarkers(stObj,logFC = 0.5)
stObj <- findSVAPA(stObj)
geneSet <- findstAPASet(stObj)

##draw APA Spatial Express picture
drawSpatialExpress(stObj,"Pde1c")
```

## Citation
Ji G, Tang Q, Zhu S, Zhu J, Ye P, Xia S, Wu X: stAPAminer: Mining Spatial Patterns of Alternative Polyadenylation for Spatially Resolved Transcriptomic Studies. bioRxiv 2022:2022.2007.2020.500789.

