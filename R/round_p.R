#' Calculating the coefficient of variation for a numeric vector.

#' @param p Numeric; numeric vector.
#'
#' @return Correct rounding for p values based on Ornithological Applications.
#' @export

p_rounder <- function(p){
  p <- as.numeric(p)
  if(p > 1){print("¡p must be smaller than 1!")}
  if(p >= 0.01){pval <- round(p, 2)}
  if(p < 0.01 & p >= 0.001){pval <- round(p, 3)}
  if(p < 0.001 & p >= 0.0001){pval <- round(p, 4)}
  if(p < 0.0001){pval <- "p < 0.0001"}
  
  return(pval)
}