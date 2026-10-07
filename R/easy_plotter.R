#' Easy Plotting with Error Bars
#'
#' Creates simple ggplots with data with error bars of your choice shown.
#'
#' @param data A data frame.
#' @param response The response variable.
#' @param group Character; grouping factor.
#' @param error_measure The measurement used for the error bars. Options include "se" for standard error, "sd" for standard deviation, and "ci" for 95% confidence interval.
#' @param theme_base ggplot2 theme.
#'
#' @return A plot with the desired parameters.
#' @export

easy_plotter <- function(data = NULL, 
                         response = NULL, 
                         group = NULL,
                         error_measure = "se",
                         theme_base = theme_classic()){
  if(is.null(data)){stop("Data must be defined.")}
  if(is.null(response)){stop("Response must be defined.")}
  if(is.null(group)){stop("Group must be defined.")}
  
  data |> 
    select(.data[[group]], .data[[response]]) |> 
    na.omit() |> 
    group_by(.data[[group]]) |> 
    summarise(mu = mean(.data[[response]]),
              se = se(.data[[response]]),
              sd = sd(.data[[response]]),
              ci = qnorm(0.975)*(sd(.data[[response]])/sqrt(nrow(data)))) ->
    summary_data
  
  a <- ggplot(summary_data, 
         aes(x = .data[[group]],
             y = mu)) +
    geom_point() +
    theme_base +
    geom_errorbar(aes(ymin = mu - .data[[error_measure]], 
                      ymax = mu + .data[[error_measure]])) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  return(a)
}