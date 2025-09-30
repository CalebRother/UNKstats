#' @noRd
.as_factor <- function(x) if (is.factor(x)) x else factor(x)

#' @noRd
.stop_miss <- function(df, cols) {
  miss <- cols[!cols %in% names(df)]
  if (length(miss)) stop("Column(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
}

#' One-way ANOVA (student-friendly)
#'
#' Runs one-way ANOVA, Welch, or Kruskal–Wallis automatically and returns
#' test info, tidy tables, and a plot.
#'
#' @export
run_oneway <- function(
    data,
    dv,                     # numeric column name, e.g. "Score"
    group,                  # factor with 3+ levels; if 2 levels, it still works but t-tests may be better
    adjust = c("tukey","sidak","holm","bonferroni","BH"),
    show_means = c("point","point+ci","none"),
    theme_base = ggplot2::theme_bw()
) {
  adjust <- match.arg(adjust)
  show_means <- match.arg(show_means)

  .stop_miss(data, c(dv, group))

  df <- data[stats::complete.cases(data[, c(dv, group), drop = FALSE]), , drop = FALSE]
  if (nrow(df) < nrow(data)) {
    message("Removed ", nrow(data) - nrow(df), " row(s) with missing values in {", dv, ", ", group, "}." )
  }

  # coerce types
  df[[group]] <- .as_factor(df[[group]])
  if (!is.numeric(df[[dv]])) stop("`dv` must be numeric.", call. = FALSE)

  # build formula and call your engine
  fml <- stats::as.formula(paste(dv, "~", group))
  res <- auto_anova(
    formula = fml,
    data = df,
    which_factor = group,
    pairwise_adjust = adjust,
    show_means = show_means,
    theme_base = theme_base
  )

  # Attach a light APA one-liner for printing
  class(res) <- c("teach_anova_result", class(res))
  attr(res, "meta") <- list(
    design = "oneway_between",
    dv = dv,
    factors = list(between = group),
    adjust = adjust
  )
  res
}

#' Factorial ANOVA (student-friendly)
#'
#' Runs factorial ANOVA or ART as needed and returns tidy tables + plot.
#'
#' @export
run_factorial <- function(
    data,
    dv,                       # numeric
    between,                  # character vector of 2+ factor column names
    which_factor = NULL,      # which term to letter/pairwise (e.g., "A", "B", or "A:B")
    adjust = c("sidak","tukey","holm","bonferroni","BH"),
    show_means = c("point","point+ci","none"),
    theme_base = ggplot2::theme_bw()
) {
  adjust <- match.arg(adjust)
  show_means <- match.arg(show_means)

  if (length(between) < 2) stop("`between` must have 2 or more factors for factorial designs.", call. = FALSE)
  .stop_miss(data, c(dv, between))

  df <- data[stats::complete.cases(data[, c(dv, between), drop = FALSE]), , drop = FALSE]
  if (nrow(df) < nrow(data)) {
    message("Removed ", nrow(data) - nrow(df), " row(s) with missing values in {", dv, ", ", paste(between, collapse=", "), "}." )
  }

  # coerce
  for (b in between) df[[b]] <- .as_factor(df[[b]])
  if (!is.numeric(df[[dv]])) stop("`dv` must be numeric.", call. = FALSE)

  # default which_factor = first main effect if not given
  if (is.null(which_factor)) which_factor <- between[1]

  rhs <- paste(between, collapse = " * ")
  fml <- stats::as.formula(paste(dv, "~", rhs))

  res <- auto_anova(
    formula = fml,
    data = df,
    which_factor = which_factor,
    pairwise_adjust = adjust,
    show_means = show_means,
    theme_base = theme_base
  )

  class(res) <- c("teach_anova_result", class(res))
  attr(res, "meta") <- list(
    design = "factorial_between",
    dv = dv,
    factors = list(between = between, which_factor = which_factor),
    adjust = adjust
  )
  res
}

#' @exportS3Method print teach_anova_result
print.teach_anova_result <- function(x, ...) {
  meta <- attr(x, "meta")
  ti   <- x$test_info
  atbl <- x$anova_table

  cat("\nDesign:", meta$design, "\n")
  cat("DV:", meta$dv, "\n")
  if (!is.null(meta$factors$between)) {
    cat("Between factor(s):", paste(meta$factors$between, collapse = ", "), "\n")
  }
  if (!is.null(meta$factors$which_factor)) {
    cat("Letters/posthocs on:", meta$factors$which_factor, "\n")
  }
  cat("Adjustment:", meta$adjust, "\n")
  cat("Assumptions → Normality:", ifelse(ti$normality_ok, "OK", "fail"),
      sprintf("(%s)", ti$normality_note), "; Homogeneity:",
      ifelse(ti$homogeneity_ok, "OK", "fail"), "\n")
  if (isTRUE(ti$welch_preference_applied)) {
    cat("Decision:", ti$welch_preference_reason, "\n")
  }
  cat("Chosen omnibus:", ti$test, "\n")

  # Show a compact ANOVA table if present
  if (!is.null(atbl)) {
    try({
      keep <- intersect(
        c("term","df","df1","df2","statistic","p.value","method",
          "Sum Sq","Mean Sq","F value","Pr(>F)"),
        names(atbl)
      )
      print(utils::head(atbl[keep], 10), row.names = FALSE)
    }, silent = TRUE)
  }

  invisible(x)
}
