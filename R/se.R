#' Calculating standard error from a string of numbers
#'
#' @param x Numeric; numeric vector.
#'
#' @return A value for the standard error.
#' @export

se <- function(x){
  n <- length(x) # calculate n
  s <- sd(x) # calculate standard deviation
  se_val <- s/sqrt(n)
  return(se_val)
}