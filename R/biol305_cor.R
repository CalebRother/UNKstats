#' Calculating the coefficient of variation for a numeric vector.

#' @param x Numeric; numeric vector.
#' @param y Numeric; numeric vector.
#' @param method Character; method of test to perform. Options are "pearson", "kendall", "spearman". Defaults to "pearson".
#'
#' @return A printout of values needed for correlations, including values for testing them and p values for reporting.
#' @export

biol305_cor <- function(x=NA, y=NA, method = "pearson"){
  if(is.data.frame(x)==T){
    if(ncol(x)==2){
      r <- cor(x[,1], x[,2], method = method)
    }else{
      r <- cor(x, method = method)
    }
    
    r2 <- r
    r[r==1|r==-1] <- 0
    
    n <- 2*nrow(x)
  }else{
    r <- cor(x, y, method = method)
    
    n <- 2*length(x)
  }
  
  t_val <- r*sqrt((n-2)/(1-r^2))
  
  p <- pt(t_val, df = n - 2)
  
  p[p > 0.5] <- 1 - p[p > 0.5]
  p[p > 0.005] <- round(p[p > 0.005],2)
  p[p > 0.0005] <- round(p[p > 0.0005],3)
  p[p > 0.00005] <- round(p[p > 0.00005],4)
  p[p < 0.00005] <- "< 0.0001"
  
  if(is.data.frame(x)==T){
    print("Correlation:")
    print(round(r, 2))
    if(ncol(x) == 2){
      print(paste0("Degrees of freedom: ", n - 2))
      print(paste0("t value: ", round(t_val, 2)))
      print(paste0("P value: ", p))
    }
    if(ncol(x) > 2){
      print(paste0("Degrees of freedom: ", n - 2))
      print("")
      print("t value: ")
      print(round(t_val, 2))
      print("")
      print("P value: ")
      print(p)
    }
  }else{
    print(paste0("Correlation: ", round(r, 2)))
    print(paste0("Degrees of freedom: ", n - 2))
    print(paste0("t value: ", round(t_val, 2)))
    print(paste0("P value: ", p))
  }
}