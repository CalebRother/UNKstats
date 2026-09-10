#' Calculating the confidence interval when mean, sd, n, and CI are known
#'
#' @param percent Numeric; A number between 0 and 1 corresponding to the confidence interval level. If not defined, defaults to 0.95.
#' @param s Numeric; the standard deviation for the given data.
#' @param n Numeric; the number of records in a dataset; if blank, defaults to 1.
#' @param data Numeric vector; a vector of numberic data; if used, arguments for s and n ignored.
#'
#'
#' @return A value for the confidence interval, given the value on either side of the mean that encompasses the interval.
#' @export

confidence_interval <- function(percent = 0.95, s, n=1, data=NULL){
    if(is.null(data)==F){
        if(is.vector(data, mode = "numeric") == F){
            stop("Data must be a numeric vector.")
        }
        xbar <- mean(data)
        print(paste0("Mean = ",xbar))
        n <- length(data)
        s <- sd(data)
    }
    if(percent < 0 | percent > 1){
        stop("Percent must be between 0 and 1.")
    }
    Z <- qnorm(1-((1-percent)/2))
    
    ci <- Z*(s/sqrt(n))
    return(ci)
}
