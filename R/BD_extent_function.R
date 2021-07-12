#' Extension of BD framework for phylogenetic and functional beta diversity
#'
#' @param Y A matrix of species occurrence or abundance per site.
#' @param dist_spp A matrix containing all trait or phylogenetic distances
#'     (upper and lower diagonal elements) among species.
#' @param nperm Scalar indicating the number of permutations to be used in taxa and site shuffle null-models.
#' @param method Character. If "raw" the BD will be computed using directly matrix Y (with reliable transformation),
#'     if "distance" BD will be computed by using distance matrix.
#'
#' @return A list containing the observed values of BD (BDextend.obs), LCBD (LCBDextend.obs), SCBD (if method
#'     is raw) and the respective values of probability calculating accordingly to site shuffle null models (ptaxa.BDextend, ptaxa.LCBDextend)
#'
#' @export
#'
#' @examples
BD.extent <- function(Y,
                      dist_spp,
                      nperm=999,
                      method = "raw"){
  if(dim(Y)[1] <= 1){
    stop("\n matrix Y must be at least three communities\n")
  }
  if(sum(is.na(match(colnames(Y), rownames(dist_spp)))) > 0){
    stop("\n there are species in community matrix no present in phylogeny\n")
  }
  if(is.matrix(Y) == FALSE){
    Y <- as.matrix(Y)
  }

  if(sum(is.na(match(colnames(Y), rownames(dist_spp)))) >= 1){
    stop("\n species names in Y and phylogenetic tree must be the same\n")
  }
  Yorg<- Y[,match(colnames(Y), rownames(dist_spp))]
  fuzzy_mat<- SYNCSA::matrix.p(comm = Yorg, phylodist = dist_spp, notification = FALSE)$matrix.P
  if(method =="raw"){
    BDextend_all<- adespatial::beta.div(Y = fuzzy_mat, method = "chord", sqrt.D = FALSE, samp = TRUE, nperm = 0)
    SCBDextend<- BDextend_all$SCBD #SCBD component
  } else{
    BDextend_all<- adespatial::beta.div(Y = fuzzy_mat, method = "percentdiff", sqrt.D = TRUE, samp = TRUE, nperm = 0) #distance based calculation
  }
  BDextend<- BDextend_all$beta[2] #BD measure
  LCBDextend<- BDextend_all$LCBD #LCBD measure
  BD_allNull<- matrix(NA, nrow = nperm, ncol= 1, dimnames= list(paste("run",1:nperm, sep = ""), "BDextend_null"))
  LCBDextend_allNull<- matrix(NA, nrow = nperm, ncol= nrow(Y), dimnames= list(paste("run", 1:nperm, sep=""), 1:nrow(Y)))
  SAMP <- SYNCSA::permut.vector(ncol(Y), nset = nperm)
  for(i in 1:nperm) {
    distspp_null <- dist_spp[SAMP[i,],SAMP[i,]] #taxa shuffle matrix P
    fuzzyExtend_null<-SYNCSA::matrix.p(comm = Y, phylodist = distspp_null, notification = FALSE)$matrix.P
    if(method =="raw"){
      BDextend_allNull<- adespatial::beta.div(Y = fuzzyExtend_null, method = "chord", sqrt.D = FALSE, samp = TRUE, nperm = 0)
    } else{
      BDextend_allNull<- adespatial::beta.div(Y = fuzzyExtend_null, method = "percentdiff", sqrt.D = TRUE, samp = TRUE, nperm = 0)
    }
    BD_allNull[i,1]<- BDextend_allNull$beta[2]
    LCBDextend_allNull[i,]<- BDextend_allNull$LCBD
  }
  ptaxa.BDextend<- (sum(ifelse(BD_allNull[,"BDextend_null"] >= BDextend, 1, 0)) + 1)/(nperm + 1)
  ptaxa.LCBDextend<- numeric(length = nrow(Y))
  for(j in 1:nrow(Y)){
    ptaxa.LCBDextend[j]<- (sum(ifelse(LCBDextend_allNull[,j] >= LCBDextend[j], 1, 0)) + 1)/(nperm + 1)
  }
  if(method =="raw"){
    return(list(BDextend.obs= BDextend, LCBDextend.obs= LCBDextend, SCBDextend.obs= SCBDextend,
                ptaxa.BDextend= ptaxa.BDextend, ptaxa.LCBDextend= ptaxa.LCBDextend))
  } else{
    return(list(BDextend.obs= BDextend, LCBDextend.obs= LCBDextend,
                ptaxa.BDextend= ptaxa.BDextend, ptaxa.LCBDextend= ptaxa.LCBDextend))
  }
}
