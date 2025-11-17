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
#' @param spaghetti Logical. (Currently unused for plotting; reserved for future.)
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
    spaghetti = TRUE,         # kept for API, not used in current plot
    sphericity_correction = c("GG","HF","none")  # report table using this correction
) {

  adjust <- match.arg(adjust)
  show_means <- match.arg(show_means)
  sphericity_correction <- match.arg(sphericity_correction)

  # work on a local copy
  df <- data

  # coerce/check
  if (!is.numeric(df[[dv]]))
    stop("`dv` must be numeric.", call. = FALSE)

  df[[id]] <- factor(df[[id]])
  for (w in within) df[[w]] <- factor(df[[w]])

  # afex model: type III SS, sphericity tests included
  afex_mod <- afex::aov_ez(
    id    = id,
    dv    = dv,
    data  = df,
    within = within,
    type  = 3,
    return = "afex_aov"
  )

  # ANOVA table with requested correction
  an_tbl <- afex::nice(
    afex_mod,
    es         = "pes",
    correction = sphericity_correction,
    observed   = NULL
  )
  an_tbl <- as.data.frame(an_tbl)

  # choose term for posthocs/letters
  if (is.null(which_within)) {
    which_within <- if (length(within) == 1) {
      within[1]
    } else {
      paste(within, collapse = ":")
    }
  }

  # emmeans on chosen within term
  emm <- try(
    emmeans::emmeans(
      afex_mod,
      specs = stats::as.formula(paste0("~", which_within))
    ),
    silent = TRUE
  )

  pw <- if (!inherits(emm, "try-error")) {
    emmeans::contrast(emm, "pairwise", adjust = adjust)
  } else {
    NULL
  }

  posthoc_tbl <- if (!is.null(pw)) broom::tidy(pw) else NULL

  # Letters (CLD)
  letters_df <- NULL
  if (!inherits(emm, "try-error")) {
    cld <- try(
      multcomp::cld(emm, Letters = letters, adjust = adjust),
      silent = TRUE
    )
    if (!inherits(cld, "try-error")) {
      cld <- as.data.frame(cld)
      cld$.group <- gsub("\\s+", "", cld$.group)

      if (grepl(":", which_within)) {
        facs <- strsplit(which_within, ":")[[1]]
        cld$combo <- interaction(cld[, facs, drop = FALSE], drop = TRUE)
        letters_df <- cld[, c("combo", ".group")]
      } else {
        letters_df <- cld[, c(which_within, ".group")]
        names(letters_df)[1] <- which_within
      }
    }
  }

  # --------- Summary data for plotting (grouped bar style) -------------------

  # group by all within factors and compute mean/SE
  summary_df <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(within))) |>
    dplyr::summarise(
      mean = mean(.data[[dv]], na.rm = TRUE),
      se   = stats::sd(.data[[dv]], na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    )

  # choose x and fill to mirror run_twoway style
  if (length(within) == 1) {
    x_var   <- within[1]
    fill_var <- NULL
  } else {
    x_var   <- within[1]
    fill_var <- within[2]
  }

  # attach letters to summary_df
  y_pad <- 0.05 * diff(range(df[[dv]], na.rm = TRUE))

  if (!is.null(letters_df)) {
    if ("combo" %in% names(letters_df)) {
      # interaction letters
      summary_df <- summary_df |>
        dplyr::mutate(combo = interaction(dplyr::across(dplyr::all_of(strsplit(which_within, ":")[[1]])), drop = TRUE))
      label_df <- dplyr::left_join(summary_df, letters_df, by = "combo") |>
        dplyr::mutate(y = mean + se + y_pad)
    } else {
      join_var <- names(letters_df)[1]  # e.g. which_within
      label_df <- dplyr::left_join(summary_df, letters_df, by = join_var) |>
        dplyr::mutate(y = mean + se + y_pad)
    }
  } else {
    label_df <- NULL
  }

  # --------- Plot: grouped bar + SE + letters -------------------------------

  pd <- ggplot2::position_dodge(width = if (is.null(fill_var)) 0.6 else 0.7)

  if (is.null(fill_var)) {
    p <- ggplot2::ggplot(summary_df, ggplot2::aes(x = .data[[x_var]], y = mean))
  } else {
    p <- ggplot2::ggplot(
      summary_df,
      ggplot2::aes(x = .data[[x_var]], y = mean, fill = .data[[fill_var]])
    )
  }

  p <- p +
    ggplot2::geom_col(position = pd, width = 0.6) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean - se, ymax = mean + se),
      position = pd,
      width = 0.2
    )

  if (!is.null(label_df)) {
    p <- p +
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(label = .group, y = y),
        position = pd,
        vjust = 0,
        size = 5
      )
  }

  if (show_means %in% c("point", "point+ci")) {
    p <- p + ggplot2::geom_point(
      data = summary_df,
      ggplot2::aes(y = mean),
      position = pd,
      color = "black",
      size = 2
    )
  }

  p <- p +
    theme_base +
    ggplot2::labs(
      x = if (length(within) == 1) within[1] else paste(within, collapse = " × "),
      y = dv,
      subtitle = paste0(
        "RM-ANOVA (", sphericity_correction,
        " correction in table); posthocs: ", adjust
      )
    )

  # --------- Return object ---------------------------------------------------

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
    test_info   = test_info,
    model       = afex_mod,
    anova_table = an_tbl,
    emmeans     = if (!inherits(emm, "try-error")) emm else NULL,
    posthoc     = posthoc_tbl,
    letters     = letters_df,
    plot        = p
  )
  class(res) <- c("teach_anova_result", class(res))
  attr(res, "meta") <- list(
    design  = "rm_within",
    dv      = dv,
    factors = list(within = within, which_factor = which_within),
    adjust  = adjust
  )
  res
}
