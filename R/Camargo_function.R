#' Camargo evenness index
#'
#' @param n_spec numeric vector containing abundances or presence of species in communities
#' @param include_zeros Logical, if TRUE (default) the absences are accounted in the calculation of mean
#'
#' @return Numeric value of Camargo Evenness measure
#' @export
#'
#' @examples
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
    sum.camar<- abs(sum(as.dist(camar, diag= FALSE, upper= FALSE)))
    return(1-sum.camar)
}
