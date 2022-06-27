setClass("stAPAminerObj",representation(count = "data.frame",metaData = "data.frame" , normCount = "data.frame" , seurat = "Seurat" , DEAPA = "list",LSAPA = "data.frame",SVAPA = "data.frame"))


#' @title Create a stAPAminer object
#' @description Create a stAPAminer object from raw data.
#' @param count APA index
#' @param metaData Coordinate information for spatial transcriptome data
#'
#' @return
#' @export
#'
#' @examples
createStAPAminerObject <-function(count,metaData){
  count<-count[,rownames(metaData)]
  seurat<- CreateSeuratObject(count)
  seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  normCount <- as.data.frame(apply(as.matrix(seurat@assays$RNA@data),1,relative_func),stringsAsFactors = F)
  normCount[is.na(normCount)]<-0
  obj<-new("stAPAminerObj",count=count,metaData=metaData,normCount=normCount,seurat=seurat,DEAPA=list(),LSAPA=data.frame(),SVAPA=data.frame())
  return(obj)
}

relative_func <- function(expres) {
  maxd = max(expres) - min(expres)
  rexpr = (expres - min(expres))/maxd
  return(rexpr)
}



#' @title Preprocessing and Clustering
#' @description Create a stAPAminer object from raw data.
#' @param stAPAminerObj A stAPAminer object.
#' @param nfeatures Number of features to select as top variable features.
#' @param dims Dimensions of reduction to use as input.
#' @param resolution Value of the resolution parameter.
#' @param k Defines k for the k-nearest neighbor algorithm
#'
#' @return
#' @export
#'
#' @examples
makeStCluster <- function(stAPAminerObj,nfeatures = 2000,dims = 1:5,resolution = 0.5,k=20){
  seurat <- stAPAminerObj@seurat
  seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = nfeatures)
  seurat <- ScaleData(seurat)
  seurat <- RunPCA(seurat,features = VariableFeatures(object = seurat))
  seurat <- FindNeighbors(seurat, dims = dims, k.param = k)
  seurat <- FindClusters(seurat, resolution = resolution)
  stAPAminerObj@seurat <- seurat
  return(stAPAminerObj)
}



#' @title Spatial Scatter Plot
#' @description Graphs the output of a spatial scatter plot where each point is a spot and it's positioned based on the additional coordinate information.
#' By default, spots are colored by the active.ident.
#' @param stAPAminerObj A stAPAminer object.
#' @param size The size of spots.
#'
#' @return
#' @export
#'
#' @examples
findLabels<-function(stAPAminerObj,size=4){
  seurat<-stAPAminerObj@seurat
  cluster <- data.frame(cluster = seurat@meta.data$seurat_clusters,barcode = rownames(seurat@meta.data))
  cluster_nums<-length(unique(seurat@active.ident))
  layer <- stAPAminerObj@metaData
  p<-ggplot(data = cluster,aes(x = layer$x,y = layer$y))+geom_point(aes(color=cluster),size=size,shape=16)+
    labs(x = NULL, y = NULL,color="cluster")+scale_colour_manual(values = rainbow(cluster_nums))+
    theme(panel.grid =element_blank(),plot.title = element_text(hjust = 0.5),legend.position = 'right')
  return(p)
}

#' @title Rename the object's identity classes
#' @description Rename the object's identity classes.
#' @param stAPAminerObj A stAPAminer object.
#' @param layer_names The new object's identity classes
#'
#' @return
#' @export
#'
#' @examples
RenameCluster<-function(stAPAminerObj,layer_names=c("c1","c2","c3","c4","c5")){
  names(layer_names) <- levels(stAPAminerObj@seurat)
  stAPAminerObj@seurat <- RenameIdents(stAPAminerObj@seurat,layer_names)
  return(stAPAminerObj)
}

#' @title Draw the colored spatial scatter plot
#' @description Graphs the output of a spatial scatter plot where each point is a spot and it's positioned based on the additional coordinate information.
#' By default, spots are colored by the active.ident.
#'
#' @param stAPAminerObj A stAPAminer object.
#' @param title the title of the plot.
#' @param labelOrder the order of legend.
#' @param color The color level.
#' @param islegend whether to display legend.
#' @param size The size of spot.
#'
#' @return
#' @export
#'
#' @examples
#'
drawClusterPlot<-function(stAPAminerObj,title = NULL,labelOrder = levels(stAPAminerObj@seurat@active.ident),color = rainbow(cluster_nums) ,size = 3,islegend = TRUE){
  cluster<-as.data.frame(stAPAminerObj@seurat@active.ident)
  cluster_nums<-length(unique(stAPAminerObj@seurat@active.ident))
  colnames(cluster)<-"cluster"
  layer <- stAPAminerObj@metaData
  p <- ggplot(data = cluster,aes(x = layer$x,y = layer$y))+geom_point(aes(color=cluster),size=size,shape=16)+
    labs(x = NULL, y = NULL,title = title,color="")+
    scale_colour_manual(breaks = labelOrder,values = color)+
    theme(panel.grid =element_blank(),plot.title = element_text(hjust = 0.5),legend.position = 'bottom',
          panel.background = element_rect(fill = "#FFFFFF", colour = NULL),axis.text = element_blank(),axis.ticks = element_blank())
  if(islegend == FALSE){
    p<- p + guides(color=F)
  }
  return(p)
}

#' @title Compute Pearson correlation coefficient
#' @description Compute the Pearson correlation coefficient for each spot in the different layers.
#' @param stAPAminerObj A stAPAminer object.
#' @param source Source name for the object.
#'
#' @return
#' @export
#'
#' @examples
computePearsonIndex<-function(stAPAminerObj,source="RUD"){
  ratio <- stAPAminerObj@count
  label <- stAPAminerObj@metaData
  la <- levels(factor(label$label))
  pearson_b<-data.frame()
  for(i in 1:length(la)){
    pearson<-cor(ratio[,which(label$label==la[i])],method="pearson")
    pearson<-pearson[upper.tri(pearson, diag = FALSE)]
    pearson_b<-rbind(pearson_b,data.frame(label = la[i],value = pearson,source=source))
  }
  return(pearson_b)
}

#' @title Cluster evaluation of different indexes
#' @description Cluster evaluation of different indexes,Cluster evaluation of different indexes,
#' including 4 internal indexes(DBI, CH, SC, Denn) and 4 external indexes(ARI, jaccard, purity, NMI).
#' @param stAPAminerObj A stAPAminer object.
#'
#' @return
#' @export
#'
#' @examples
computeIndex<-function(stAPAminerObj){
  layer <- stAPAminerObj@metaData
  data <- t(as.data.frame(stAPAminerObj@seurat@assays$RNA@counts))
  x = as.integer(stAPAminerObj@seurat$seurat_clusters)
  y = as.integer(as.factor(layer$label))
  dbi<-index.DB(data,x)$DB
  chi<-index.G1(data,x)
  sc<-index.S(as.matrix(dist(data)),x)
  denn<-intern.dunn(x,data)
  ari<-external_validation(y,x,method="adjusted_rand_index")
  jaccard<-external_validation(y,x,method="jaccard_index")
  purity<-external_validation(y,x,method="purity")
  nmi<-external_validation(y,x,method="nmi")
  index<-list(DBI = dbi,CH = chi,SC = sc,Denn = denn,ARI = ari,jaccard = jaccard,purity = purity,NMI = nmi)
  return(index)
}


#' @title Visualize individual gene expression levels on a spatial scatter plot
#' @description Colors single spots on a spatial scatter plot according to the gene expression.
#' @param stAPAminerObj A stAPAminer object.
#' @param geneName The gene name
#' @param title the title of the plot.
#' @param size the size of the plot.
#'
#' @return
#' @export
#'
#' @examples
drawSpatialExpress<-function(stAPAminerObj,geneName,title = geneName,size = 2){
  layer<-stAPAminerObj@metaData
  count<-stAPAminerObj@normCount
  if(!geneName%in%colnames(count)){
    p<-ggplot(data = count,aes(x = layer$x,y = layer$y))+geom_point(size=size,shape=16,color="#BBBBBB")+
      labs(x = NULL, y = NULL, title=title ,color="Value")+
      theme(panel.grid =element_blank(),legend.position = 'right',legend.text = NULL,
            panel.background = element_rect(fill = "#FFFFFF", colour = NULL)) + theme(plot.title = element_text(hjust = 0.5))+
      guides(color = guide_colorbar(title = "",title.theme = element_text(angle = 90,size = 10)))+
      scale_y_continuous(breaks = NULL)+scale_x_continuous(breaks = NULL)+
      guides(color=F)
    return(p)
  }
  p<-ggplot(data = count,aes(x = layer$x,y = layer$y,color = get(geneName)))+geom_point(size=size,shape=16)+
    labs(x = NULL, y = NULL, title=title ,color="Value")+
    theme(panel.grid =element_blank(),legend.position = 'right',legend.text = NULL,
          panel.background = element_rect(fill = "#FFFFFF", colour = NULL)) + theme(plot.title = element_text(hjust = 0.5))+
    guides(color = guide_colorbar(title = "",title.theme = element_text(angle = 90,size = 10)))+
    scale_color_viridis_c()+
    scale_y_continuous(breaks = NULL)+scale_x_continuous(breaks = NULL)+
    guides(color=F)
  return(p)
}



#' @title APA index markers for pairs of identity classes
#' @description Finds markers (differentially expressed APA index) for identity classes
#' @param stAPAminerObj A stAPAminer object.
#' @param padj Maximum adjusted p-value.
#' @param logFC Minimum log fold-change of the average expression between the two groups.
#' @param levels Identity classes to be contrasted
#'
#' @return
#' @export
#'
#' @examples
findContrastAPA<-function(stAPAminerObj,padj = 0.05,logFC = 0.5,levels=c("GCL","MCL","OPL","GL","ONL")){
  stAPAminerObj@DEAPA$all <- c()
  stAPAminerObj@seurat@active.ident = factor(stAPAminerObj@seurat@active.ident, levels=levels)
  layersPair<-combn(levels(stAPAminerObj@seurat@active.ident),2)
  for(i in 1:ncol(layersPair)){
    layer1<-layersPair[1,i]
    layer2<-layersPair[2,i]
    print(paste0("Calculating cluster ",layer1,"-",layer2))
    markers <- FindMarkers(stAPAminerObj@seurat, ident.1 = layer1,ident.2 = layer2)
    if("avg_logFC" %in% colnames(markers)){
      markers <- markers[markers$p_val_adj<padj&abs(markers$avg_logFC)>logFC,]
    }
    if("avg_log2FC" %in% colnames(markers)){
      markers <- markers[markers$p_val_adj<padj&abs(markers$avg_log2FC)>logFC,]
    }
    stAPAminerObj@DEAPA[[paste0(layer1,'-',layer2)]]<-markers

    stAPAminerObj@DEAPA$all <- c(stAPAminerObj@DEAPA$all,rownames(markers))
  }

  stAPAminerObj@DEAPA$all <- unique(stAPAminerObj@DEAPA$all)

  return(stAPAminerObj)
}


#' @title APA index markers for all identity classes
#' @description Finds markers (differentially expressed APA index) for identity classes
#' @param stAPAminerObj A stAPAminer object.
#' @param padj Maximum adjusted p-value.
#' @param logFC Minimum log fold-change of the average expression between the two groups.
#'
#' @return
#' @export
#'
#' @examples
findAllMarkers<-function(stAPAminerObj,padj = 0.05,logFC = 0.5){
  LSAPA<-FindAllMarkers(stAPAminerObj@seurat)
  if("avg_logFC" %in% colnames(LSAPA)){
    LSAPA <- LSAPA[LSAPA$p_val_adj<padj&abs(LSAPA$avg_logFC)>logFC,]
  }
  if("avg_log2FC" %in% colnames(LSAPA)){
    LSAPA <- LSAPA[LSAPA$p_val_adj<padj&abs(LSAPA$avg_log2FC)>logFC,]
  }

  stAPAminerObj@LSAPA<-LSAPA
  return(stAPAminerObj)
}


#' @title Find genes with spatial patterns of APA usage variation
#' @description we also adopted SPARK to detect genes with spatial patterns of APA usage variation. SVAPA genes mark distinct spatial APA usage patterns in the global spatial context.
#'
#' @param stAPAminerObj A stAPAminer object.
#' @param padj Maximum adjusted p-value.
#'
#' @return
#' @export
#'
#' @examples
findSVAPA<-function(stAPAminerObj,padj=0.05){
  if(!("SPARK" %in% installed.packages()[, "Package"])){
    devtools::install_github('xzhoulab/SPARK')
  }
  require(SPARK)
  count<-stAPAminerObj@count
  count<-count[which(apply(count,1,var)!=0),]
  count<-10**count
  spark <- CreateSPARKObject(count,
                             location=stAPAminerObj@metaData[,1:2],
                             percentage = 0.1,
                             min_total_counts = 10)
  spark@lib_size <- apply(spark@counts, 2, sum)
  spark <- spark.vc(spark,
                    covariates = NULL,
                    lib_size = spark@lib_size,
                    num_core = 3,
                    verbose = T)
  spark <- spark.test(spark, check_positive = T,verbose = F)
  sv_gene<-spark@res_mtest[which(spark@res_mtest$adjusted_pvalue<padj),c(11,12)]
  stAPAminerObj@SVAPA<-sv_gene
  return(stAPAminerObj)
}


findstAPASet<-function(stAPAminerObj){
  set<-c(stAPAminerObj@DEAPA$all,stAPAminerObj@LSAPA$gene,rownames(stAPAminerObj@SVAPA))
  set<-unique(set)
  return(set)
}


drawVolcano<-function(stAPAminerObj,layer1,layer2,padj = 0.05,logFC = 0.5,tittle=NULL,color=c("#2f5688","#BBBBBB","#CC0000"),legend=FALSE){
  data<- as.data.frame(t(stAPAminerObj@normCount),stringsAsFactors = F)
  data<-data[,colSums(data)>0]
  layer<-stAPAminerObj@metaData
  layer<-layer[colnames(data),]
  group_list = factor(layer$label)
  design <- model.matrix(~0+group_list)
  rownames(design) = colnames(data)
  colnames(design) <- levels(group_list)
  cont.matrix <- makeContrasts(contrasts = c(paste0(layer1,'-',layer2)), levels = design)
  dge <- DGEList(counts = data)
  dge <- calcNormFactors(dge)
  #logCPM<-cpm(dge,log=TRUE,prior.count = 3)
  #fit <- lmFit(logCPM, design)
  v <- voom(dge, design, plot=TRUE) #会自动计算log(cpm)值
  #拟合线性模型
  fit <- lmFit(v, design)
  #针对给定的对比计算估计系数和标准误差
  fit2 <- contrasts.fit(fit, cont.matrix)
  #计算出t统计量，F统计量和差异表达倍数的对数
  fit2 <- eBayes(fit2)
  allDEG <- topTable(fit2, coef = paste0(layer1,'-',layer2), n = Inf)
  allDEG <- na.omit(allDEG)
  deg<-allDEG
  deg$logP<- -log10(deg$adj.P.Val)
  deg$group<-'no-significant'
  if(nrow(deg[which((deg$adj.P.Val<padj)&(deg$logFC > logFC)),])>0){
    deg[which((deg$adj.P.Val<padj)&(deg$logFC > logFC)),]$group="up-regulated"
  }
  if(nrow(deg[which((deg$adj.P.Val<padj)&(deg$logFC < -logFC)),])>0){
    deg[which((deg$adj.P.Val<padj)&(deg$logFC < -logFC)),]$group="down-regulated"
  }
  table(deg$group)
  deg$label = ""
  deg<-deg[order(deg$adj.P.Val),]
  up.genes<-head(rownames(deg[which(deg$group=="up-regulated"),]),10)
  down.genes<-head(rownames(deg[which(deg$group=="down-regulated"),]),10)
  deg.top10.genes<-c(as.character(up.genes),as.character(down.genes))
  deg$label[match(deg.top10.genes,rownames(deg))]<-deg.top10.genes
  p<-ggscatter(deg,x="logFC",y="logP",color = 'group',palette = color,size = 0.5,
               label = deg$label,font.label = 5,repel = T,xlab = "log2FoldChange",ylab = "-log10(Adjust P-value)")+
    #theme_base()+
    ggtitle(tittle)+theme(plot.title = element_text(hjust = 0.5),
                          axis.text = element_text(size = 8,family = "Arial"),axis.title =element_text(size = 8,family = "Arial"))+
    geom_hline(yintercept = -log10(padj),linetype="dashed")+
    geom_vline(xintercept = c(-logFC,logFC),linetype="dashed")
  if(legend==FALSE){
    p<-p+guides(color=F)
  }
  return(p)
}

#' @title Identification of APA genes with spatial patterns
#' @description Identification of APA genes with spatial patterns of APA usage variation
#' @param stAPAminerObj A stAPAminer object.
#' @param features Genes to be tested with spatial patterns of variation in APA usage
#' @param k Number of spatial patterns
#'
#' @return
#' @export
#'
#' @examples
findSpatialPattern<-function(stAPAminerObj,features,k=10){
  features.filter <- features[features %in% rownames(stAPAminerObj@count)]
  count_feature<-stAPAminerObj@count[features.filter,]
  label<-stAPAminerObj@metaData
  uuKmeans<-kmeans(count_feature,k)
  pattern = data.frame()
  for (i in 1:k) {
    pattern<-rbind(pattern,apply(count_feature[uuKmeans[["cluster"]]==i,],2,mean))
  }
  colnames(pattern)<-colnames(count_feature)
  pattern<-as.data.frame(t(pattern))
  colnames(pattern)<-paste0("pattern_",c(1:k))
  uuKmeans$patternCount<-pattern
  return(uuKmeans)
}


#' @title Convert Seurat object to SingleCellExperiment object
#' @description Convert Seurat object to SingleCellExperiment object
#' @param stAPAminerObj A stAPAminer object.
#'
#' @return
#' @export
#'
#' @examples
convertToSCE <- function(stAPAminerObj){
  SCE <- as.SingleCellExperiment(stAPAminerObj@seurat)
  return(SCE)
}


#' @title Convert Sierra result files to movAPA object
#' @description Convert Sierra result files to movAPA object
#' @param path path to Sierra result file. The path needs to contain peak_annotations.txt, counts files
#' @param annodb GTF file annotating poly(A) sites information
#'
#' @return
#' @export
#'
#' @examples
readSierra <- function(path,annodb){
  library(movAPA)
  dirlist <- dir(path)

  peakdir <- dirlist[grep("_peak_annotations.txt",dirlist)]
  countdir <- dirlist[grep("_count$",dirlist)]

  peak_annotations <- read.delim(paste0(path,"/",peakdir))
  peak_annotations$coord <- 0
  peak_annotations[peak_annotations$strand == "+",]$coord <- peak_annotations[peak_annotations$strand == "+",]$end
  peak_annotations[peak_annotations$strand == "-",]$coord <- peak_annotations[peak_annotations$strand == "-",]$start

  peak_annotations$seqnames

  PACds <- peak_annotations[,c("seqnames","start","end","strand","coord","gene_id")]
  colnames(PACds) <- c("chr","start","end","strand","coord","gene_id")
  PACds <- unique(PACds)
  anno <- annotatePAC(PACds, aGFF = annodb) # movAPA annotation
  anno <- arrange(anno, chr,ftr_start)

  count <- ReadPeakCounts(paste0(path,countdir,"/"))
  count <- as.data.frame(count)
  intpas <- rownames(anno)[rownames(anno) %in% rownames(count)]

  count <- count[intpas,]
  anno <- anno[intpas,]
  anno$gene_range <- rownames(anno)
  anno$gene <- gsub("[.].*","",anno$gene)
  anno <- anno[!is.na(anno$gene),]
  p <- readPACds(anno, colDataFile=NULL,noIntergenic=FALSE)
  p@anno$ftr_raw <- p@anno$ftr
  p <- ext3UTRPACds(p, 1000)

  panno <- arrange(p@anno, chr,ftr_start)

  count <- count[panno$gene_range,]
  rownames(count) <- rownames(panno)
  puni <- unique(panno[,c("chr","coord","strand","ftr","gene")])

  colnames(count) <- gsub("-.*","",colnames(count))
  p@counts <- count[rownames(puni),]
  p@anno <- panno[rownames(puni),]
  return(p)

}

ReadPeakCounts <- function(data.dir = NULL,
                           mm.file = NULL,
                           barcodes.file = NULL,
                           sites.file = NULL) {

  if (is.null(data.dir) & is.null(mm.file) & is.null(barcodes.file) & is.null(sites.file)) {
    stop("Please provide either a directory or file names.")
  }

  if (!is.null(data.dir)) {

    ## First check if files are compressed
    file.list <- list.files(data.dir)
    if (sum(endsWith(file.list, ".gz")) == 3) {
      gzipped = TRUE
    } else{
      gzipped = FALSE
    }

    if (gzipped) {
      mm.file <- paste0(data.dir, "/matrix.mtx.gz")
      barcodes.file <- paste0(data.dir, "/barcodes.tsv.gz")
      sites.file <- paste0(data.dir, "/sitenames.tsv.gz")
    } else {
      mm.file <- paste0(data.dir, "/matrix.mtx")
      barcodes.file <- paste0(data.dir, "/barcodes.tsv")
      sites.file <- paste0(data.dir, "/sitenames.tsv")
    }

  } else{

    if (is.null(mm.file) | is.null(barcodes.file) | is.null(sites.file)) {
      stop("No directory provided, but an input file appears to be missing. Please check.")
    }

    ## Individual files proivded - check if compressed
    file.list <- c(mmfile, barcodes.file, sites.file)

    if (sum(endsWith(file.list, ".gz")) == 3) {
      gzipped = TRUE
    } else{
      gzipped = FALSE
    }
  }

  if (gzipped) {
    count.mat <- Matrix::readMM(gzfile(mm.file))

    barcodes.con <- gzfile(barcodes.file)
    barcodes <- readLines(barcodes.con)
    close(barcodes.con)

    peaks.con <- gzfile(sites.file)
    peaks <- readLines(peaks.con)
    close(barcodes.con)
  } else {
    count.mat <- Matrix::readMM(mm.file)

    barcodes <- readLines(barcodes.file)
    peaks <- readLines(sites.file)
  }


  colnames(count.mat) <- barcodes
  rownames(count.mat) <- peaks

  return(count.mat)
}
