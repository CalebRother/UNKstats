#' Finding the mode for a numeric vector.
#' Based on Statology function
#' define function to calculate mode
#' works on vectors of data

#' @param x Numeric; numeric vector.
#'
#' @return A value for the mode of a data string. Indicates if no mode is found.
#' @export

find_mode <- function(x) {
  # get unique values from vector
  u <- unique(x)
  # count number of occurrences for each value
  tab <- tabulate(match(x, u))
  
  # if no mode, say so
  if(length(x)==length(u[tab == max(tab)])){
    print("No mode.")
  }else{
    # return the value with the highest count
    mode_x <- u[tab == max(tab)]
    count_mode <- max(tab)
    print(paste0("Mode: ", mode_x))
    print(paste0("Mode count: ", count_mode))
  }
}