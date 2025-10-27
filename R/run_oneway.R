#' One-way comparison (ANOVA or Kruskal–Wallis only)
#'
#' @param data A data frame.
#' @param dv Character; numeric outcome column.
#' @param group Character; grouping factor (≥3 levels recommended).
#' @param parametric Logical; if TRUE runs classic one-way ANOVA (+ Tukey HSD),
#'   if FALSE (default) runs Kruskal–Wallis (+ Dunn, BH).
#' @param adjust P-adjust for ANOVA Tukey letters display label only
#'   (TukeyHSD uses Tukey internally). Ignored for Kruskal (BH is used).
#' @param show_means "point","point+ci","none".
#' @param theme_base ggplot2 theme.
#'
#' @return A list with test_info, model (ANOVA only), anova_table,
#'   posthoc, letters, plot. Class "teach_anova_result".
#' @export
run_oneway <- function(
    data,
    dv,
    group,
    parametric = FALSE,
    adjust = c("tukey","sidak","holm","bonferroni","BH"),
    show_means = c("point","point+ci","none"),
    theme_base = ggplot2::theme_bw()
) {
  adjust <- match.arg(adjust)
  show_means <- match.arg(show_means)

  .stop_miss(data, c(dv, group))

  df <- data[stats::complete.cases(data[, c(dv, group), drop = FALSE]), , drop = FALSE]
  if (nrow(df) < nrow(data)) {
    message("Removed ", nrow(data) - nrow(df),
            " row(s) with missing values in {", dv, ", ", group, "}." )
  }

  # coerce types
  if (!is.numeric(df[[dv]])) stop("`dv` must be numeric.", call. = FALSE)
  df[[group]] <- .as_factor(df[[group]])

  # helpers ---------------------------------------------------------------
  fml <- stats::as.formula(paste(dv, "~", group))

  mean_cl_normal_local <- function(x, conf = 0.95) {
    x <- x[is.finite(x)]
    n  <- length(x); m <- mean(x)
    se <- stats::sd(x) / max(1, sqrt(n))
    mult <- if (n > 1) stats::qt((1 + conf)/2, df = n - 1) else 0
    data.frame(y = m, ymin = m - mult * se, ymax = m + mult * se)
  }

  pmat_to_letters <- function(pairs_df, g1, g2, pcol, alpha = 0.05) {
    groups <- sort(unique(c(as.character(pairs_df[[g1]]), as.character(pairs_df[[g2]]))))
    M <- matrix(1, length(groups), length(groups), dimnames = list(groups, groups))
    for (i in seq_len(nrow(pairs_df))) {
      a <- as.character(pairs_df[[g1]][i]); b <- as.character(pairs_df[[g2]][i])
      p <- as.numeric(pairs_df[[pcol]][i]); M[a,b] <- p; M[b,a] <- p
    }
    multcompView::multcompLetters(M < alpha, compare = "<")$Letters
  }

  # run -------------------------------------------------------------------
  if (isTRUE(parametric)) {
    # -------- One-way ANOVA + TukeyHSD --------
    fit <- stats::aov(fml, data = df)
    an_tbl <- broom::tidy(stats::anova(fit))

    # TukeyHSD returns matrix with rownames like "B-A"
    tk <- stats::TukeyHSD(fit)[[group]]
    posthoc <- tibble::tibble(
      contrast = rownames(tk),
      diff = tk[, "diff"],
      lwr  = tk[, "lwr"],
      upr  = tk[, "upr"],
      p.adj = tk[, "p adj"]
    )
    # split contrast "B-A" → group1, group2
    spl <- strsplit(posthoc$contrast, "-")
    posthoc$group2 <- vapply(spl, `[`, character(1), 1)
    posthoc$group1 <- vapply(spl, `[`, character(1), 2)

    # letters from Tukey p-values
    letters_vec <- pmat_to_letters(posthoc, "group1", "group2", "p.adj", alpha = 0.05)
    letters_df <- tibble::tibble(!!group := names(letters_vec), .group = unname(letters_vec))

    subtitle_txt <- "One-way ANOVA + Tukey HSD"

    model <- fit

  } else {
    # -------- Kruskal–Wallis + Dunn (BH) --------
    kw <- stats::kruskal.test(fml, data = df)
    an_tbl <- tibble::tibble(
      term = group,
      statistic = unname(kw$statistic),
      df = unname(kw$parameter),
      p.value = unname(kw$p.value),
      method = "Kruskal–Wallis"
    )

    dunn <- rstatix::dunn_test(df, fml, p.adjust.method = "BH")
    posthoc <- dplyr::transmute(dunn, group1 = group1, group2 = group2, p.adj = p.adj)

    letters_vec <- pmat_to_letters(posthoc, "group1", "group2", "p.adj", alpha = 0.05)
    letters_df <- tibble::tibble(!!group := names(letters_vec), .group = unname(letters_vec))

    subtitle_txt <- "Kruskal–Wallis + Dunn (BH)"
    model <- NULL
  }

  # plot ------------------------------------------------------------------
  y_range <- range(df[[dv]], na.rm = TRUE)
  y_pad   <- 0.05 * diff(y_range)
  y_top   <- max(df[[dv]], na.rm = TRUE) + y_pad

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[group]], y = .data[[dv]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.5, size = 1.6) +
    {
      if (show_means == "point") {
        list(ggplot2::stat_summary(fun = "mean", geom = "point", size = 2.6, shape = 21))
      } else if (show_means == "point+ci") {
        list(
          ggplot2::stat_summary(fun = "mean", geom = "point", size = 2.6, shape = 21),
          ggplot2::stat_summary(fun.data = mean_cl_normal_local, geom = "errorbar", width = 0.2)
        )
      } else NULL
    } +
    ggplot2::geom_text(
      data = letters_df,
      ggplot2::aes(x = .data[[group]], y = y_top, label = .group),
      vjust = 0, size = 5
    ) +
    theme_base +
    ggplot2::labs(x = group, y = dv, subtitle = subtitle_txt)

  # return ----------------------------------------------------------------
  out <- list(
    test_info = list(
      test = if (parametric) "oneway_anova_tukey" else "kruskal_dunn",
      parametric = parametric
    ),
    model = model,
    anova_table = an_tbl,
    posthoc = posthoc,
    letters = letters_df,
    plot = p
  )
  class(out) <- c("teach_anova_result", class(out))
  attr(out, "meta") <- list(
    design = "oneway_between",
    dv = dv,
    factors = list(between = group),
    adjust = if (parametric) adjust else "BH",
    parametric = parametric
  )
  out
}
