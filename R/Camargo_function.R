#' Camargo´s index of diversity
#'
#' @param n_spec numeric vector with the number of species for each community
#' @param include_zeros Logical indicating if zeroes will be included in the calculation
#'     default is TRUE
#'
#' @return A scalar
#'

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
