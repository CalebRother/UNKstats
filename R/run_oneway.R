#' Calculating the coefficient of variation for a numeric vector.
#'
#' @param x Numeric; numeric vector.
#'
#' @return A value for the coefficient of variation for a given numeric vector.
#' @export
cv <- function(x) {
  stats::sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}


#' One-way comparison (ANOVA or Kruskal–Wallis only)
#'
#' Performs a one-way between-subjects comparison on a numeric outcome.
#' If `parametric = TRUE`, runs a classic one-way ANOVA with Tukey HSD
#' post-hoc tests. If `parametric = FALSE` (default), runs Kruskal–Wallis
#' with Dunn post-hoc tests (BH-adjusted).
#'
#' If `blocking` is supplied (e.g., a block or subject ID), a randomized
#' block ANOVA is fit via `dv ~ group + blocking` when `parametric = TRUE`.
#' Blocking is *not* supported for the Kruskal–Wallis version.
#'
#' @param data A data frame.
#' @param dv Character; numeric outcome column.
#' @param group Character; grouping factor (≥3 levels recommended).
#' @param blocking Character; optional block factor column name for randomized
#'   block ANOVA (parametric only). Use `run_rm()` for repeated-measures designs.
#' @param parametric Logical; if TRUE runs classic one-way ANOVA (+ Tukey HSD),
#'   if FALSE (default) runs Kruskal–Wallis (+ Dunn, BH).
#' @param adjust P-adjust for ANOVA Tukey letters display label only
#'   (TukeyHSD uses Tukey internally). Ignored for Kruskal (BH is used).
#' @param show_means "point","none" — controls mean point display.
#' @param theme_base ggplot2 theme.
#'
#' @return A list with test_info, model (ANOVA only), anova_table,
#'   posthoc, letters, plot. Class "teach_anova_result".
#' @export
run_oneway <- function(
    data,
    dv,
    group,
    blocking   = NULL,
    parametric = FALSE,
    adjust     = c("tukey","sidak","holm","bonferroni","BH"),
    show_means = c("point","none"),
    theme_base = ggplot2::theme_bw()
) {

  adjust     <- match.arg(adjust)
  show_means <- match.arg(show_means)

  # Guard: blocking only implemented for parametric
  if (!is.null(blocking) && !isTRUE(parametric)) {
    stop("`blocking` is currently supported only when `parametric = TRUE` ",
         "(randomized block ANOVA). For nonparametric blocked designs, ",
         "use `run_rm()` or a Friedman-type test.",
         call. = FALSE)
  }

  # --------------------------------------------------------------------
  # Drop incomplete cases, including blocking if present
  # --------------------------------------------------------------------
  if (is.null(blocking)) {
    df <- data[stats::complete.cases(data[, c(dv, group), drop = FALSE]),
               , drop = FALSE]
  } else {
    df <- data[stats::complete.cases(data[, c(dv, group, blocking), drop = FALSE]),
               , drop = FALSE]
    df[[blocking]] <- as.factor(df[[blocking]])
  }

  if (nrow(df) < nrow(data)) {
    msg_cols <- if (is.null(blocking)) {
      paste0(dv, ", ", group)
    } else {
      paste0(dv, ", ", group, ", ", blocking)
    }
    message("Removed ", nrow(data) - nrow(df),
            " row(s) with missing values in {", msg_cols, "}." )
  }

  # --------------------------------------------------------------------
  # Coerce types
  # --------------------------------------------------------------------
  if (!is.numeric(df[[dv]])) {
    stop("`dv` must be numeric.", call. = FALSE)
  }
  df[[group]] <- as.factor(df[[group]])

  # --------------------------------------------------------------------
  # Build formula: with or without blocking
  # --------------------------------------------------------------------
  if (is.null(blocking)) {
    fml <- stats::as.formula(paste(dv, "~", group))
  } else {
    fml <- stats::as.formula(paste(dv, "~", group, "+", blocking))
  }

  # --------------------------------------------------------------------
  # Helper: p-value matrix -> letters
  # --------------------------------------------------------------------
  pmat_to_letters <- function(pairs_df, g1, g2, pcol, alpha = 0.05) {
    groups <- sort(unique(c(as.character(pairs_df[[g1]]),
                            as.character(pairs_df[[g2]]))))
    M <- matrix(1, length(groups), length(groups),
                dimnames = list(groups, groups))
    for (i in seq_len(nrow(pairs_df))) {
      a <- as.character(pairs_df[[g1]][i])
      b <- as.character(pairs_df[[g2]][i])
      p <- as.numeric(pairs_df[[pcol]][i])
      M[a, b] <- p
      M[b, a] <- p
    }
    multcompView::multcompLetters(M < alpha, compare = "<")$Letters
  }

  # --------------------------------------------------------------------
  # Run test
  # --------------------------------------------------------------------
  if (isTRUE(parametric)) {
    # -------- One-way (or randomized block) ANOVA + TukeyHSD --------
    fit    <- stats::aov(fml, data = df)
    an_tbl <- broom::tidy(stats::anova(fit))

    tk <- stats::TukeyHSD(fit)[[group]]
    posthoc <- tibble::tibble(
      contrast = rownames(tk),
      diff     = tk[, "diff"],
      lwr      = tk[, "lwr"],
      upr      = tk[, "upr"],
      p.adj    = tk[, "p adj"]
    )
    spl <- strsplit(posthoc$contrast, "-")
    posthoc$group2 <- vapply(spl, `[`, character(1), 1)
    posthoc$group1 <- vapply(spl, `[`, character(1), 2)

    letters_vec <- pmat_to_letters(posthoc, "group1", "group2", "p.adj", alpha = 0.05)
    letters_df  <- tibble::tibble(!!group := names(letters_vec),
                                  .group   = unname(letters_vec))

    subtitle_txt <- if (is.null(blocking)) {
      "One-way ANOVA + Tukey HSD"
    } else {
      "Randomized block ANOVA + Tukey HSD"
    }
    model <- fit

  } else {
    # -------- Kruskal–Wallis + Dunn (BH) --------
    kw <- stats::kruskal.test(fml, data = df)
    an_tbl <- tibble::tibble(
      term      = group,
      statistic = unname(kw$statistic),
      df        = unname(kw$parameter),
      p.value   = unname(kw$p.value),
      method    = "Kruskal–Wallis"
    )

    dunn <- rstatix::dunn_test(df, fml, p.adjust.method = "BH")
    posthoc <- dplyr::transmute(dunn,
                                group1 = group1,
                                group2 = group2,
                                p.adj  = p.adj)

    letters_vec <- pmat_to_letters(posthoc, "group1", "group2", "p.adj", alpha = 0.05)
    letters_df  <- tibble::tibble(!!group := names(letters_vec),
                                  .group   = unname(letters_vec))

    subtitle_txt <- "Kruskal–Wallis + Dunn (BH)"
    model        <- NULL
  }

  # --------------------------------------------------------------------
  # Plot
  # --------------------------------------------------------------------
  y_range <- range(df[[dv]], na.rm = TRUE)
  y_pad   <- 0.05 * diff(y_range)
  y_top   <- max(df[[dv]], na.rm = TRUE) + y_pad

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[group]], y = .data[[dv]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.5, size = 1.6) +
    {
      if (show_means == "point") {
        ggplot2::stat_summary(
          fun  = "mean",
          geom = "point",
          size = 2.6,
          shape = 21,
          fill  = "white"
        )
      } else {
        NULL
      }
    } +
    ggplot2::geom_text(
      data = letters_df,
      ggplot2::aes(x = .data[[group]], y = y_top, label = .group),
      vjust = 0, size = 5
    ) +
    theme_base +
    ggplot2::labs(x = group, y = dv, subtitle = subtitle_txt)

  # --------------------------------------------------------------------
  # Return
  # --------------------------------------------------------------------
  out <- list(
    test_info = list(
      test       = if (parametric) "oneway_anova_tukey" else "kruskal_dunn",
      parametric = parametric
    ),
    model       = model,
    anova_table = an_tbl,
    posthoc     = posthoc,
    letters     = letters_df,
    plot        = p
  )

  class(out) <- c("teach_anova_result", class(out))
  attr(out, "meta") <- list(
    design    = "oneway_between",
    dv        = dv,
    factors   = list(
      between = group,
      block   = if (!is.null(blocking)) blocking else NULL
    ),
    adjust    = if (parametric) adjust else "BH",
    parametric = parametric
  )

  out
}
