# library(movAPA)
# library(org.Mm.eg.db)



#' @title Compute APA index
#' @description Calculates APA index by different metrics for each condition.
#' @param APA A PACdataset.
#' @param barcodes Selected cells or spots.
#' @param method Can be WUL, RUD, SLR. WUL: Weighted UTR length (Jia 2017; Ultisky) (3UTRAPA, choose2PA=NULL/PD/MOST) RUD: Relative Usage of Distal Poly(A) Sites (Ji 2009) (3UTRAPA with/without non-3UTR PACs, choose2PA=NULL/PD/MOST) only calc distal ratio, RUD=distal/gene SLR: short to long ratio (Begik 2017) (3UTRAPA, choose2PA=PD/most) #only for PD, SLR=short/long GPI: geometric proximal index (Shulman 2019) (3UTRAPA, choose2PA=NULL/PD/most) #GPI=proximal/geo_mean_PD
#'
#' @return
#' @export
#'
#' @examples
computeAPAIndex<-function(APA,barcodes = NULL,method="RUD"){
  if(!is.null(barcodes)){
    APA@counts<-APA@counts[,barcodes]
  }
  index<-movAPAindex(APA,method=method,choose2PA=NULL, RUD.includeNon3UTR=FALSE,clearPAT=0)
  index<-index[-nrow(index),]
  # geneSymbol <-  mapIds(org.Mm.eg.db,rownames(index),'SYMBOL','ENTREZID')
  # rownames(index)<-geneSymbol
  return(index)
}


#' @title Imputation APA index
#' @description A KNN-based imputation model for recovering APA signals in the APA index matrix.
#' @param index APA index to be recovered
#' @param gene reference gene expression matrix
#' @param k Defines k for the k-nearest neighbor algorithm
#' @param init Whether to initialize the APA index
#'
#' @return
#' @export
#'
#' @examples
imputeAPAIndex <- function (index, gene, k = 10, init = TRUE)
{
  colNames <- colnames(index)[colnames(index) %in% colnames(gene)]
  index <- index[, colNames]
  print(paste0("NA values before imputed: ", sum(is.na(index))))
  gene <- gene[, colNames]
  if (init == TRUE) {
    gene_some <- gene[rownames(index), ]
    is_zero <- gene_some == 0
    index[is_zero & is.na(index)] = 0
    print(paste0("NA values after initial imputed: ", sum(is.na(index))))
  }
  scaleData <- as.data.frame(t(scale(gene)), stringsAsFactors = F)
  dist <- as.matrix(dist(scaleData, method = "euclidean"))
  weight <- data.frame()
  for (i in 1:nrow(dist)) {
    weight <- rbind(weight, order(dist[i, ])[2:(k + 1)])
  }
  rownames(weight) <- rownames(dist)
  colnames(weight) <- paste0("N", 1:k)
  num = 1
  while (sum(is.na(index)) > 0 & num <= 10) {
    print(paste0("imputing....................", num))

    for (i in 1:ncol(index)) {

      target <- index[, i] # 第i个细胞下每个基因的RUD表达量：
      if(sum(is.na(target)) == 0) next
      tianchong <- index[, as.numeric(weight[i, ])] # 第i个细胞最近的k个细胞的RUD表达量:gene x k
      tianchong.mean <- apply(tianchong,1,mean,na.rm = TRUE)
      target[is.na(target)] <- tianchong.mean[is.na(target)]
      index[, i] <- target
    }
    print(paste0("NA values : ", sum(is.na(index))))
    num = num + 1
  }
  if (num > 10) {
    print("over the max impute times,convert all Na to Zero")
    index[is.na(index)] = 0
  }
  return(index)
}
