#' Calculating the Z score when the sample mean, population mean, and standard deviation are known.
#'
#' @param xbar Numeric; value representing the sample mean.
#' @param mu Numeric; value representing the population mean.
#' @param sd.x Numeric; value representing the standard deviation of the sample, if known. Can be replaced with population mean in some scenarios.
#' @param n Numeric; value representing the sample size. Defaults to 1 if not defined.
#'
#' @return A value for the Z score.
#' @export

zscore <- function(xbar, mu, sd.x, n = 1){
  z <- (xbar - mu)/(sd.x/sqrt(n))
  return(z)
}