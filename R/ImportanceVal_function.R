#' Importance Values (IV) of diversity metrics
#'
#' @param matrix.M Matrix containing diversity metrics (columns) values for each sample (lines)
#' @param scale Logical. If TRUE (default) the diversity metrics are scaled for maximun values of each metric
#' @param method Method used to scale matrix M, default is scaling each diversity metric by their maximmun value
#' @param stopRule Logical. If TRUE (default) a stop rule are used to select only the axis
#'     of the PCA that matched the rule (Kaiser-Gutmann criterion)
#'
#' @return A list of length 3 containing:
#'     1- The IV for each diversity metric (IV.obs_stopRule if stopRule = TRUE)
#'     2- The proportion of variation accounted by each metric in each axis
#'     3- The correlation (Pearson) among all metrics
#' @export
#'
#' @examples
#'     \donotrun{
#'
#'     data(rodents_diversity)
#'     IV_rodents <- ImportanceVal(matrixM = rodents_diversity, method = "max", stopRule = TRUE)
#'     IV_rodents$IV.obs_stopRule # IV results
#'
#'     }
ImportanceVal<- function(matrix.M, scale= TRUE, method= "max", stopRule= TRUE){
  if(is.matrix(matrix.M) == FALSE){
    matrix.M<- as.matrix(matrix.M)
    if(ncol(matrix.M)<3){
      stop("\n matrix M must be at least 3 components of diversity\n")
    }
    if(nrow(matrix.M)<3){
      stop("\n Matrix M must be at least 3 communities\n")
    }
  }
  matrix.M.stand<-vegan::decostand(x = matrix.M, method = method, MARGIN = 2)[1:nrow(matrix.M),]
  if(scale == TRUE){
    metric.sqrt.corr<- (prcomp(x = matrix.M.stand, scale.= FALSE)$rotation ^ 2)
    prop.var<- summary(prcomp(x = matrix.M.stand, scale.= FALSE))$importance[2,]
    IVs.resul<- matrix(nrow= ncol(metric.sqrt.corr), ncol= ncol(matrix.M), dimnames= list(paste("PC",1:ncol(metric.sqrt.corr)), colnames(matrix.M)))
    for(i in 1:nrow(metric.sqrt.corr)){
      IVs.resul[,i]<- metric.sqrt.corr[i,] * as.matrix(prop.var)
    }
    IVs.proportion<-matrix(nrow= ncol(metric.sqrt.corr), ncol= ncol(matrix.M), dimnames= list(paste("PC",1:ncol(metric.sqrt.corr)), colnames(matrix.M)))
    for(i in 1:nrow(IVs.resul)){
      IVs.proportion[i,]<-IVs.resul[i,]/prop.var[i]
    }
    if(stopRule==TRUE){
      sig.eigen<-which(prcomp(matrix.M.stand, scale. = FALSE)$sdev^2>mean(prcomp(matrix.M.stand, scale. = FALSE)$sdev^2))
      IVs.resul.sig<-IVs.resul[sig.eigen,]
      IV.resul.sig<- setNames(list(IVs.resul.sig, prop.var, metric.sqrt.corr), c("IV.obs_stopRule", "Var.by.axis","Metrics_correlation"))
      return(IV.resul.sig)
    } else{
      IV.resul<- setNames(list(IVs.resul, prop.var, metric.sqrt.corr), c("IV.obs", "Var.by.axis","Metrics_correlation"))
      return(IV.resul)
    }
  }
  if(scale == FALSE){
    metric.sqrt.corr<- (prcomp(x = matrix.M, scale.= FALSE)$rotation ^ 2)
    prop.var<- summary(prcomp(x = matrix.M, scale.= FALSE))$importance[2,]
    names.matrix.IV<- list("IV.resu", colnames(matrix.M))
    IVs.resul<- matrix(nrow= ncol(metric.sqrt.corr), ncol= ncol(matrix.M), dimnames= names.matrix.IV)
    for(i in 1:nrow(metric.sqrt.corr)){
      IVs.resul[,i]<- metric.sqrt.corr[i,] * as.matrix(prop.var)
    }
    if(stopRule==TRUE){
      sig.eigen<-which(prcomp(matrix.M, scale. = FALSE)$sdev^2>mean(prcomp(matrix.M, scale. = FALSE)$sdev^2))
      IVs.resul.sig<-IVs.resul[sig.eigen,]
      IV.resul.sig<- setNames(list(IVs.resul.sig, prop.var, metric.sqrt.corr), c("IV.obs_stopRule", "Var.by.axis","Metrics_correlation"))
      return(IV.resul.sig)
    } else{
      IV.resul<- setNames(list(IVs.resul, prop.var, metric.sqrt.corr), c("IV.obs", "Var.by.axis","Metrics_correlation"))
      return(IV.resul)
    }
  }
}
