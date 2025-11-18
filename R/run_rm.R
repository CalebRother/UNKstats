#' Repeated-measures ANOVA (student-friendly)
#'
#' Fits an RM-ANOVA via `afex::aov_ez()`, reports GG/HF corrections,
#' and provides optional emmeans pairwise comparisons with letters.
#'
#' Visualizes results with boxplots and jittered raw points.
#'
#' @param data A data frame.
#' @param dv Character. Numeric DV column name.
#' @param id Character. Subject identifier column.
#' @param within Character vector of within-subject factor(s).
#' @param which_within Character. Term for emmeans/posthocs
#'   (e.g. "Time" or "Time:Condition"). If NULL, uses the single within factor,
#'   or the interaction of all within factors if there are 2+.
#' @param adjust P-adjust for emmeans ("sidak","tukey","holm","bonferroni","BH").
#' @param show_means "point" or "none".
#' @param theme_base ggplot2 theme.
#' @param spaghetti Logical. Reserved for future (currently not used in plot).
#' @param sphericity_correction "GG","HF","none" for ANOVA table.
#'
#' @return A list with test_info, model, anova_table, emmeans, posthoc, letters, plot.
#' @export
run_rm <- function(
    data,
    dv,
    id,
    within,
    which_within = NULL,
    adjust = c("sidak","tukey","holm","bonferroni","BH"),
    show_means = c("point","none"),
    theme_base = ggplot2::theme_bw(),
    spaghetti = TRUE,
    sphericity_correction = c("GG","HF","none")
) {

  adjust     <- match.arg(adjust)
  show_means <- match.arg(show_means)
  sphericity_correction <- match.arg(sphericity_correction)

  # --- Data Prep -------------------------------------------------------------
  df <- data
  if (!is.numeric(df[[dv]])) stop("`dv` must be numeric.", call. = FALSE)

  df[[id]] <- factor(df[[id]])
  for (w in within) df[[w]] <- factor(df[[w]])

  # --- Fit RM-ANOVA via afex -------------------------------------------------
  afex_mod <- afex::aov_ez(
    id     = id,
    dv     = dv,
    data   = df,
    within = within,
    type   = 3,
    return = "afex_aov"
  )

  an_tbl <- afex::nice(
    afex_mod,
    es         = "pes",
    correction = sphericity_correction,
    observed   = NULL
  )
  an_tbl <- as.data.frame(an_tbl)

  # --- Post-hoc / Emmeans ----------------------------------------------------
  if (is.null(which_within)) {
    which_within <- if (length(within) == 1) within[1] else paste(within, collapse = ":")
  }

  emm <- try(
    emmeans::emmeans(
      afex_mod,
      specs = stats::as.formula(paste0("~", which_within))
    ),
    silent = TRUE
  )

  if (!inherits(emm, "try-error")) {
    pw          <- emmeans::contrast(emm, "pairwise", adjust = adjust)
    posthoc_tbl <- broom::tidy(pw)
  } else {
    posthoc_tbl <- NULL
  }

  # --- Letters (CLD) – OLD STYLE with multcomp -------------------------------
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
        # interaction term like "Time:Condition"
        facs <- strsplit(which_within, ":")[[1]]
        cld$combo <- interaction(cld[, facs, drop = FALSE], drop = TRUE)
        letters_df <- cld[, c("combo", ".group")]
      } else {
        # single within factor
        letters_df <- cld[, c(which_within, ".group")]
        names(letters_df)[1] <- which_within
      }
    }
  }

  # --- Summary Data (Means + Max for labeling) -------------------------------
  summary_df <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(within))) |>
    dplyr::summarise(
      mean    = mean(.data[[dv]], na.rm = TRUE),
      max_val = max(.data[[dv]],  na.rm = TRUE),
      .groups = "drop"
    )

  # choose x + fill aesthetics
  if (length(within) == 1) {
    x_var    <- within[1]
    fill_var <- NULL
  } else {
    x_var    <- within[1]
    fill_var <- within[2]
  }

  y_range <- diff(range(df[[dv]], na.rm = TRUE))
  y_pad   <- 0.05 * y_range

  # attach letters → label_df
  if (!is.null(letters_df)) {
    if ("combo" %in% names(letters_df)) {
      # interaction letters
      summary_df <- summary_df |>
        dplyr::mutate(
          combo = interaction(
            dplyr::across(dplyr::all_of(strsplit(which_within, ":")[[1]])),
            drop = TRUE
          )
        )
      label_df <- dplyr::left_join(summary_df, letters_df, by = "combo") |>
        dplyr::mutate(y = max_val + y_pad)
    } else {
      # single factor letters
      join_var <- names(letters_df)[1]
      label_df <- dplyr::left_join(summary_df, letters_df, by = join_var) |>
        dplyr::mutate(y = max_val + y_pad)
    }
  } else {
    label_df <- NULL
  }

  # --- Plotting --------------------------------------------------------------
  pd <- ggplot2::position_dodge(width = if (is.null(fill_var)) 0.6 else 0.7)

  # base ggplot
  if (is.null(fill_var)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = .data[[dv]]))
  } else {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data[[x_var]], y = .data[[dv]], fill = .data[[fill_var]])
    )
  }

  # box + jitter
  p <- p +
    ggplot2::geom_boxplot(
      position      = pd,
      outlier.shape = NA,
      width         = 0.6
    ) +
    ggplot2::geom_jitter(
      position = if (is.null(fill_var)) {
        ggplot2::position_jitter(width = 0.12)
      } else {
        ggplot2::position_jitterdodge(
          jitter.width = 0.12,
          dodge.width  = 0.7
        )
      },
      alpha = 0.5,
      size  = 1.4
    )

  # letters
  if (!is.null(label_df) && !all(is.na(label_df$.group))) {
    if (is.null(fill_var)) {
      p <- p +
        ggplot2::geom_text(
          data = label_df,
          ggplot2::aes(
            x     = .data[[x_var]],
            y     = y,
            label = .group
          ),
          vjust = 0,
          size  = 5
        )
    } else {
      p <- p +
        ggplot2::geom_text(
          data = label_df,
          ggplot2::aes(
            x     = .data[[x_var]],
            y     = y,
            label = .group,
            group = .data[[fill_var]]
          ),
          position = pd,
          vjust    = 0,
          size     = 5
        )
    }
  }

  # optional mean points
  if (show_means == "point") {
    if (is.null(fill_var)) {
      p <- p +
        ggplot2::geom_point(
          data = summary_df,
          ggplot2::aes(
            x = .data[[x_var]],
            y = mean
          ),
          color = "black",
          size  = 2
        )
    } else {
      p <- p +
        ggplot2::geom_point(
          data = summary_df,
          ggplot2::aes(
            x     = .data[[x_var]],
            y     = mean,
            group = .data[[fill_var]]
          ),
          position = pd,
          color    = "black",
          size     = 2
        )
    }
  }

  p <- p +
    theme_base +
    ggplot2::labs(
      x = if (length(within) == 1) within[1] else paste(within, collapse = " × "),
      y = dv,
      subtitle = paste0(
        "RM-ANOVA (", sphericity_correction, " correction); Post-hocs: ", adjust
      )
    )

  # --- Return ----------------------------------------------------------------
  test_info <- list(
    test           = "rm_anova",
    normality_ok   = NA,
    normality_note = "rm_model",
    homogeneity_ok = NA,
    levene         = NULL
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
