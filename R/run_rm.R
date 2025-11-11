#' @noRd
.require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is required for this function. Please install it.", call. = FALSE)
}

# ---------- Repeated-measures ANOVA (within-subjects) ----------
# - 1+ within factors
# - requires id
# - Posthocs via emmeans on the 'within' term or interaction you choose

#' Repeated-measures ANOVA (student-friendly)
#'
#' Fits an RM-ANOVA via `afex::aov_ez()`, reports GG/HF corrections,
#' and provides optional emmeans pairwise comparisons with letters.
#'
#' @param data A data frame.
#' @param dv Character. Numeric DV column name.
#' @param id Character. Subject identifier column.
#' @param within Character vector of within-subject factor(s).
#' @param which_within Character. Term for emmeans/posthocs (e.g. "Time" or "Time:Condition").
#' @param adjust P-adjust for emmeans ("sidak","tukey","holm","bonferroni","BH").
#' @param show_means "point","point+ci","none".
#' @param theme_base ggplot2 theme.
#' @param spaghetti Logical. Draw per-subject lines.
#' @param sphericity_correction "GG","HF","none" for ANOVA table.
#'
#' @return A list with test_info, model, anova_table, emmeans, posthoc, letters, plot.
#' @export
run_rm <- function(
    data,
    dv,                       # numeric dv column
    id,                       # subject id column
    within,                   # character vector of within factors (>=1)
    which_within = NULL,      # e.g., "Time" or "Time:Condition"
    adjust = c("sidak","tukey","holm","bonferroni","BH"),
    show_means = c("point","point+ci","none"),
    theme_base = ggplot2::theme_bw(),
    spaghetti = TRUE,         # draw faint subject lines
    sphericity_correction = c("GG","HF","none")  # report table using this correction
) {
  .require_pkg("afex"); .require_pkg("emmeans"); .require_pkg("ggplot2")
  adjust <- match.arg(adjust)
  show_means <- match.arg(show_means)
  sphericity_correction <- match.arg(sphericity_correction)

  # coerce
  if (!is.numeric(df[[dv]])) stop("`dv` must be numeric.", call. = FALSE)
  df[[id]] <- .as_factor(df[[id]])
  for (w in within) df[[w]] <- .as_factor(df[[w]])

  # afex model: type III SS, sphericity tests included
  afex_mod <- afex::aov_ez(
    id = id, dv = dv, data = df,
    within = within, type = 3, return = "afex_aov"
  )

  # ANOVA table with requested correction
  an_tbl <- afex::nice(afex_mod, es = "pes",
                       correction = sphericity_correction, observed = NULL)
  an_tbl <- as.data.frame(an_tbl)

  # choose term for posthocs/letters
  if (is.null(which_within)) {
    which_within <- if (length(within) == 1) within[1] else paste(within, collapse = ":")
  }

  # emmeans on chosen within term
  emm <- try(emmeans::emmeans(afex_mod, specs = stats::as.formula(paste0("~", which_within))), silent = TRUE)
  pw  <- if (!inherits(emm, "try-error")) emmeans::contrast(emm, "pairwise", adjust = adjust) else NULL
  posthoc_tbl <- if (!is.null(pw)) broom::tidy(pw) else NULL

  # Letters (optional)
  letters_df <- NULL
  if (!is.null(emm)) {
    cld <- try(multcomp::cld(emm, Letters = letters, adjust = adjust), silent = TRUE)
    if (!inherits(cld, "try-error")) {
      cld <- as.data.frame(cld)
      cld$.group <- gsub("\\s+", "", cld$.group)
      if (grepl(":", which_within)) {
        facs <- strsplit(which_within, ":")[[1]]
        cld$.x <- interaction(cld[, facs], drop = TRUE)
        letters_df <- cld[, c(".x", ".group")]
      } else {
        letters_df <- cld[, c(which_within, ".group")]
      }
    }
  }

  # Plot: box+jitter, optional spaghetti, mean±CI
  df$.x <- if (length(within) == 1) df[[within[1]]]
  else interaction(df[, within, drop = FALSE], drop = TRUE)
  .response <- dv

  mean_cl_normal_local <- function(x, conf = 0.95) {
    x <- x[is.finite(x)]
    n  <- length(x); m <- mean(x)
    se <- stats::sd(x) / max(1, sqrt(n))
    mult <- if (n > 1) stats::qt((1 + conf)/2, df = n - 1) else 0
    data.frame(y = m, ymin = m - mult * se, ymax = m + mult * se)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .x, y = .data[[.response]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.5, size = 1.6)

  if (spaghetti) {
    p <- p + ggplot2::geom_line(ggplot2::aes(group = .data[[id]]), alpha = 0.25)
  }

  if (show_means == "point" || show_means == "point+ci") {
    p <- p + ggplot2::stat_summary(fun = "mean", geom = "point", size = 2.6, shape = 21)
    if (show_means == "point+ci") {
      p <- p + ggplot2::stat_summary(fun.data = mean_cl_normal_local, geom = "errorbar", width = 0.2)
    }
  }

  if (!is.null(letters_df)) {
    y_positions <- dplyr::group_by(df, .x) |>
      dplyr::summarize(ymax = max(.data[[.response]], na.rm = TRUE), .groups = "drop")
    y_range <- range(df[[.response]], na.rm = TRUE); y_pad <- 0.05 * diff(y_range)

    letters_pos <- if (".x" %in% names(letters_df)) {
      dplyr::rename(letters_df, x = .x)
    } else {
      dplyr::rename(letters_df, x = dplyr::all_of(names(letters_df)[1]))
    }
    letters_pos <- dplyr::left_join(letters_pos, y_positions, by = c("x" = ".x")) |>
      dplyr::mutate(y = ymax + y_pad)

    p <- p + ggplot2::geom_text(data = letters_pos,
                                ggplot2::aes(x = x, y = y, label = .group),
                                vjust = 0, size = 5)
  }

  p <- p + ggplot2::labs(
    x = if (length(within) == 1) within[1] else paste(within, collapse = ":"),
    y = dv,
    subtitle = paste0("RM-ANOVA (", sphericity_correction, " correction in table); posthocs: ", adjust)
  ) + theme_base

  test_info <- list(
    test = "rm_anova",
    normality_ok = NA, normality_note = "rm_model",
    homogeneity_ok = NA, levene = NULL,
    n_groups_for_letters = if (!is.null(letters_df)) nrow(letters_df) else NA_integer_,
    letters_skipped = is.null(letters_df),
    welch_preference_applied = FALSE,
    welch_preference_reason = NULL
  )

  res <- list(
    test_info = test_info,
    model = afex_mod,
    anova_table = an_tbl,
    emmeans = if (!inherits(emm, "try-error")) emm else NULL,
    posthoc = posthoc_tbl,
    letters = letters_df,
    plot = p
  )
  class(res) <- c("teach_anova_result", class(res))
  attr(res, "meta") <- list(
    design = "rm_within",
    dv = dv,
    factors = list(within = within, which_factor = which_within),
    adjust = adjust
  )
  res
}
