#' Evenness of Eigenvalues for diversity metrics (EE)
#'
#' @param matrix.M A matrix containing diversity metrics values (columns) for each sample (lines)
#' @param scale Logical. If TRUE (default) the diversity metrics are scaled to 0 mean and variance of 1 (standardization)
#' @param method Character, the method used to perform the scaling of diversity metrics.
#'     Default is "standardize".
#' @param evenness Character, the method used to calculate the Evenness of Eigenvalues.
#'     Default is "Camargo". Other options are "Pielou" or "both"
#'
#' @return A scalar (if used "Camargo" or "Pielou") or a numeric vector of length 2
#'     if "both" in evenness argument.
#'
#' @export
#'
#' @examples
EvennessEigen <- function(matrix.M,
                           scale = TRUE,
                           method = "standardize",
                           evenness = "Camargo"){
    library(vegan)
    library(picante)
    #Camargo's index
    camargo.eveness <- function(n_spec, include_zeros = T){
        if(is.vector(n_spec)==FALSE){
            stop("\n n_spec must be a vector of abundance of species \n")
        }
        if (include_zeros){
            n <- n_spec
        }  else{
            n <- n_spec[n_spec > 0]
        }
        S <- length(n)
        camar<-matrix(nrow=length(n), ncol=length(n))
        for (i in 1:S)
        {
            for (j in 1:S)
            {
                p_i <- n[i]/sum(n)
                p_j <- n[j]/sum(n)
                camar[i,j] <- ((abs(p_i - p_j))/S)
            }
        }
        sum.camar<- abs(sum(proxy::as.dist(camar, diag= FALSE, upper= FALSE)))
        return(1-sum.camar)
    }
    #Pielou index
    pielou.evenness<- function(x) {
        H<- vegan::diversity(x)
        nspec<- length(x)
        J<- H/log(nspec)
        return(J)
    }
    if(is.matrix(matrix.M) == FALSE){
        matrix.M<- as.matrix(matrix.M)
        if(ncol(matrix.M)<3){
            stop("\n matrix M must be at least 3 components of diversity\n")
        }
        if(nrow(matrix.M)<3){
            stop("\n Matrix M must be at least 3 communities\n")
        }
    }
    matrix.M.stand <- vegan::decostand(x = matrix.M, method = method, MARGIN = 2)[1:nrow(matrix.M),]
    if(scale==TRUE){
        eingen.obs.stand <- (summary(stats::prcomp(x = matrix.M.stand, scale. = FALSE))$sdev)^2
        eveness.obs.camargo<- camargo.eveness(n_spec= as.vector(eingen.obs.stand))
        eveness.obs.pielou<- pielou.evenness(x = eingen.obs.stand)
        both_evenness<- list(Camargo= eveness.obs.camargo, Pielou= eveness.obs.pielou)
        if(evenness=="Camargo"){
            return(eveness.obs.camargo)
        }
        if(evenness=="Pielou") {
            return(eveness.obs.pielou)
        }
        if(evenness=="both"){
            return(both_evenness)
        }
    }
    if(scale==FALSE){
        eingen.obs<- (summary(stats::prcomp(x = matrix.M, scale. = FALSE))$sdev)^2
        eveness.obs.camargo<- camargo.eveness(n_spec= as.vector(eingen.obs))
        eveness.obs.pielou<- pielou.evenness(x = eingen.obs)
        both_evenness<- list(Camargo= eveness.obs.camargo, Pielou= eveness.obs.pielou)
        if(evenness=="Camargo"){
            return(eveness.obs.camargo)
        }
        if(evenness=="Pielou") {
            return(eveness.obs.pielou)
        }
        if(evenness=="both"){
            return(both_evenness)
        }
    }
}
