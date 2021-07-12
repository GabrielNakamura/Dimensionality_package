#' Pielou´s evenness index
#'
#' @param x numeric vector
#'
#' @return
#'

pielou.evenness<- function(x) {
  H<- vegan::diversity(x)
  nspec<- length(x)
  J<- H/log(nspec)
  return(J)
}
