#' Two-way factorial ANOVA
#'
#' Fits a 2×(2+) between-subjects factorial ANOVA, optionally with interaction.
#' Performs post-hoc pairwise comparisons for a chosen term via emmeans,
#' with compact letter displays and a simple boxplot.
#'
#' @param data A data frame.
#' @param dv Character; numeric dependent variable.
#' @param factor_a Character; first factor (main effect A).
#' @param factor_b Character; second factor (main effect B).
#' @param include_interaction Logical; if TRUE includes A*B interaction term (default TRUE).
#' @param which_factor Character; which term to perform post-hoc tests on
#'   (e.g. "A", "B", or "A:B"). Defaults to "A" if unspecified.
#' @param adjust Adjustment method for emmeans ("tukey","sidak","holm","bonferroni","BH").
#' @param show_means "point" or "none" for mean overlay.
#' @param theme_base ggplot2 theme.
#'
#' @return A list with test_info, model, anova_table, posthoc, letters, and plot.
#'   Class "teach_anova_result".
#' @export
run_twoway <- function(
    data,
    dv,
    factor_a,
    factor_b,
    include_interaction = TRUE,
    which_factor = NULL,
    adjust = c("tukey","sidak","holm","bonferroni","BH"),
    show_means = c("point","none"),
    theme_base = ggplot2::theme_bw()
) {
  
  .require_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required for this function. Please install it.", call. = FALSE)
  }
  .require_pkg("ggplot2"); .require_pkg("emmeans"); .require_pkg("multcomp"); .require_pkg("broom")
  
  adjust <- match.arg(adjust)
  show_means <- match.arg(show_means)
  
  # --- data prep ----------------------------------------------------------
  df <- data[stats::complete.cases(data[, c(dv, factor_a, factor_b), drop = FALSE]), , drop = FALSE]
  if (nrow(df) < nrow(data)) {
    message("Removed ", nrow(data) - nrow(df), " row(s) with missing values in {",
            dv, ", ", factor_a, ", ", factor_b, "}." )
  }
  
  if (!is.numeric(df[[dv]])) stop("`dv` must be numeric.", call. = FALSE)
  df[[factor_a]] <- as.factor(df[[factor_a]])
  df[[factor_b]] <- as.factor(df[[factor_b]])
  
  # --- formula ------------------------------------------------------------
  rhs <- paste(factor_a, factor_b, sep = if (include_interaction) " * " else " + ")
  fml <- stats::as.formula(paste(dv, "~", rhs))
  
  # --- run ANOVA ----------------------------------------------------------
  fit <- stats::aov(fml, data = df)
  an_tbl <- broom::tidy(stats::anova(fit))
  
  # --- emmeans posthoc ----------------------------------------------------
  if (is.null(which_factor)) which_factor <- factor_a
  emm <- emmeans::emmeans(fit, specs = stats::as.formula(paste0("~", which_factor)))
  pw <- emmeans::contrast(emm, method = "pairwise", adjust = adjust)
  posthoc <- broom::tidy(pw)
  
  cld <- multcomp::cld(emm, Letters = letters, adjust = adjust)
  cld <- as.data.frame(cld)
  cld$.group <- gsub("\\s+", "", cld$.group)
  
  if (grepl(":", which_factor)) {
    facs <- strsplit(which_factor, ":")[[1]]
    cld$.x <- interaction(cld[, facs], drop = TRUE)
    letters_df <- cld[, c(".x", ".group")]
  } else {
    letters_df <- cld[, c(which_factor, ".group")]
  }
  
  # --- plot ---------------------------------------------------------------
  if (grepl(":", which_factor)) {
    df$.x <- interaction(df[, strsplit(which_factor, ":")[[1]], drop = FALSE], drop = TRUE)
    xlab <- which_factor
  } else {
    df$.x <- df[[which_factor]]
    xlab <- which_factor
  }
  
  y_range <- range(df[[dv]], na.rm = TRUE)
  y_pad   <- 0.05 * diff(y_range)
  y_top   <- max(df[[dv]], na.rm = TRUE) + y_pad
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .x, y = .data[[dv]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.5, size = 1.6) +
    {
      if (show_means == "point") {
        ggplot2::stat_summary(fun = "mean", geom = "point", size = 2.6,
                              shape = 21, fill = "white")
      } else NULL
    } +
    ggplot2::geom_text(
      data = letters_df,
      ggplot2::aes(x = if ("x" %in% names(letters_df)) .data$x else .data[[which_factor]],
                   y = y_top, label = .group),
      vjust = 0, size = 5
    ) +
    theme_base +
    ggplot2::labs(
      x = xlab,
      y = dv,
      subtitle = paste0("Two-way ANOVA", if (include_interaction) " (with interaction)" else "")
    )
  
  # --- return -------------------------------------------------------------
  out <- list(
    test_info = list(
      test = "two_way_anova",
      parametric = TRUE,
      include_interaction = include_interaction,
      which_factor = which_factor
    ),
    model = fit,
    anova_table = an_tbl,
    posthoc = posthoc,
    letters = letters_df,
    plot = p
  )
  
  class(out) <- c("teach_anova_result", class(out))
  attr(out, "meta") <- list(
    design = "two_way_between",
    dv = dv,
    factors = list(A = factor_a, B = factor_b),
    adjust = adjust,
    interaction = include_interaction
  )
  
  out
}
