#' Calculating the coefficient of variation for a numeric vector.

#' @param x Numeric; numeric vector.
#'
#' @return A value for the coefficient of variation for a given numeric vector.
#' @export

cv <- function(x){
  sigma <- sd(x)
  mu <- mean(x)
  val <- sigma/mu*100
  return(val)
}