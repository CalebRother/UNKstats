#' Two-way factorial ANOVA (Grouped-Bar Version)
#'
#' Fits a 2×(2+) between-subjects factorial ANOVA, optionally including interaction.
#' Performs post-hoc pairwise comparisons for a chosen term via `emmeans`,
#' produces compact letter displays, and visualizes results with a grouped bar chart.
#'
#' @param data A data frame containing the variables to be analyzed.
#' @param dv Character; name of the numeric dependent variable.
#' @param factor_a Character; first factor (main effect A, used for x-axis grouping).
#' @param factor_b Character; second factor (main effect B, used for color fill).
#' @param include_interaction Logical; if TRUE includes A*B interaction term (default TRUE).
#' @param which_factor Character; which term to perform post-hoc tests on
#'   (e.g., "A", "B", or "A:B"). Defaults to "A" if unspecified.
#' @param adjust Adjustment method for pairwise comparisons
#'   ("tukey","sidak","holm","bonferroni","BH").
#' @param show_means "point" or "none" — controls whether to show mean points on bars.
#' @param theme_base A ggplot2 theme (default = `ggplot2::theme_bw()`).
#'
#' @return A list containing:
#' \itemize{
#'   \item `test_info` — basic metadata on the test.
#'   \item `model` — the ANOVA model object.
#'   \item `anova_table` — ANOVA summary table.
#'   \item `posthoc` — post-hoc pairwise comparison results.
#'   \item `letters` — compact letter display results.
#'   \item `plot` — grouped bar chart (ggplot object).
#' }
#' Class `"teach_anova_result"`.
#'
#' @examples
#' set.seed(123)
#' df <- expand.grid(
#'   A = factor(c("Low","High")),
#'   B = factor(c("X","Y","Z"))
#' )
#' df$Response <- rnorm(nrow(df), mean = rep(c(5,7,9,6,8,10), each = 1), sd = 1)
#' result <- run_twoway(df, dv = "Response", factor_a = "A", factor_b = "B")
#' result$plot
#'
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
  # --- dependencies ----------------------------------------------------------
  .require_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required for this function. Please install it.", call. = FALSE)
  }
  .require_pkg("ggplot2"); .require_pkg("emmeans"); .require_pkg("multcomp")
  .require_pkg("broom"); .require_pkg("dplyr"); .require_pkg("tibble"); .require_pkg("tidyr")

  adjust <- match.arg(adjust)
  show_means <- match.arg(show_means)

  # --- data prep -------------------------------------------------------------
  df <- data[stats::complete.cases(data[, c(dv, factor_a, factor_b), drop = FALSE]), , drop = FALSE]
  if (nrow(df) < nrow(data)) {
    message("Removed ", nrow(data) - nrow(df), " row(s) with missing values in {",
            dv, ", ", factor_a, ", ", factor_b, "}." )
  }
  if (!is.numeric(df[[dv]])) stop("`dv` must be numeric.", call. = FALSE)
  df[[factor_a]] <- as.factor(df[[factor_a]])
  df[[factor_b]] <- as.factor(df[[factor_b]])

  # --- formula & ANOVA -------------------------------------------------------
  rhs <- paste(factor_a, factor_b, sep = if (include_interaction) " * " else " + ")
  fml <- stats::as.formula(paste(dv, "~", rhs))
  fit <- stats::aov(fml, data = df)
  an_tbl <- broom::tidy(stats::anova(fit))

  # --- post-hoc analysis -----------------------------------------------------
  if (is.null(which_factor)) which_factor <- factor_a
  emm <- emmeans::emmeans(fit, specs = stats::as.formula(paste0("~", which_factor)))
  pw <- emmeans::contrast(emm, method = "pairwise", adjust = adjust)
  posthoc <- broom::tidy(pw)

  # --- compact letter display -----------------------------------------------
  cld <- multcomp::cld(emm, Letters = letters, adjust = adjust)
  cld <- as.data.frame(cld)
  cld$.group <- gsub("\\s+", "", cld$.group)
  if (grepl(":", which_factor)) {
    facs <- strsplit(which_factor, ":")[[1]]
    cld$.combo <- interaction(cld[, facs], drop = TRUE)
    letters_df <- cld[, c(".combo", ".group")]
    names(letters_df)[1] <- "combo"
  } else {
    letters_df <- cld[, c(which_factor, ".group")]
    names(letters_df)[1] <- which_factor
  }

  # --- summary data (means + SEs) -------------------------------------------
  summary_df <- df |>
    dplyr::group_by(.data[[factor_a]], .data[[factor_b]]) |>
    dplyr::summarise(mean = mean(.data[[dv]], na.rm = TRUE),
                     se = stats::sd(.data[[dv]], na.rm = TRUE) / sqrt(dplyr::n()),
                     .groups = "drop")

  y_pad <- 0.05 * diff(range(df[[dv]], na.rm = TRUE))

  # --- attach letters for labeling ------------------------------------------
  if (grepl(":", which_factor)) {
    summary_df <- summary_df |>
      dplyr::mutate(combo = interaction(.data[[factor_a]], .data[[factor_b]], drop = TRUE))
    label_df <- dplyr::left_join(summary_df, letters_df, by = "combo") |>
      dplyr::mutate(y = mean + se + y_pad)
  } else {
    join_var <- ifelse(which_factor == factor_a, factor_a, factor_b)
    label_df <- dplyr::left_join(summary_df, letters_df, by = join_var) |>
      dplyr::mutate(y = mean + se + y_pad)
  }

  # --- plot (grouped bar) ----------------------------------------------------
  pd <- ggplot2::position_dodge(width = 0.7)
  p <- ggplot2::ggplot(summary_df, ggplot2::aes(x = .data[[factor_a]], y = mean, fill = .data[[factor_b]])) +
    ggplot2::geom_col(position = pd, width = 0.6) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), position = pd, width = 0.2) +
    ggplot2::geom_text(data = label_df, ggplot2::aes(label = .group, y = y), position = pd, vjust = 0, size = 5) +
    {
      if (show_means == "point")
        ggplot2::geom_point(data = summary_df, ggplot2::aes(y = mean), position = pd, color = "black", size = 2)
      else NULL
    } +
    theme_base +
    ggplot2::labs(
      x = factor_a,
      y = dv,
      fill = factor_b,
      subtitle = paste0("Two-way ANOVA", if (include_interaction) " (with interaction)" else "")
    )

  # --- return ----------------------------------------------------------------
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
