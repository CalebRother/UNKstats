#' Auto-select ANOVA path and plot
#'
#' @export
auto_anova <- function(
    formula, data,
    which_factor = NULL,                         # for factorials: "A", "B", or "A:B"
    ci_level = 0.95,                             # for mean CI bars (when enabled)
    point_alpha = 0.5,
    jitter_width = 0.12,
    show_means = c("point","point+ci","none"),   # default: mean point, no CI
    theme_base = ggplot2::theme_bw(),
    # normality strategy
    normality_method = c("auto","subsample","ad","jb","none","shapiro"),
    shapiro_limit = 5000,                        # cap for Shapiro subsample
    seed = 42,                                   # reproducible subsampling
    # pairwise p-adjust for emmeans/grid comparisons
    pairwise_adjust = c("sidak","tukey","bonferroni","holm","BH"),
    # performance/control
    max_groups_for_letters = 25,                 # skip letters/posthocs when groups > this
    alpha_letters = 0.05,                        # p-threshold for p-matrix → letters
    prefer_welch_when_large_n = TRUE,            # heuristic when Levene fails
    min_n_for_welch = 15,                        # all groups ≥ this → prefer Welch
    # subtitle verbosity
    subtitle_style = c("full","compact","minimal")
) {
  # -------------------------- Arg normalization ------------------------------
  normality_method <- match.arg(normality_method)
  pairwise_adjust  <- match.arg(pairwise_adjust)
  subtitle_style   <- match.arg(subtitle_style)
  if (is.logical(show_means)) {
    show_means <- if (show_means) "point+ci" else "none"    # legacy TRUE/FALSE
  } else {
    show_means <- match.arg(show_means)
  }

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # -------------------------- Parse formula ----------------------------------
  .response <- all.vars(formula)[1]
  trms <- attr(terms(formula), "term.labels")
  k_terms <- length(trms)

  # infer which factor to letter if not provided
  if (k_terms == 1) {
    which_factor <- trms
  } else if (is.null(which_factor)) {
    which_factor <- trms[1]  # default first main effect
  }

  # -------------------------- Small utilities --------------------------------
  # local t-based CI (replaces Hmisc::mean_cl_normal)
  mean_cl_normal_local <- function(x, conf = 0.95) {
    x <- x[is.finite(x)]
    n  <- length(x); m <- mean(x)
    se <- stats::sd(x) / max(1, sqrt(n))
    mult <- if (n > 1) stats::qt((1 + conf)/2, df = n - 1) else 0
    data.frame(y = m, ymin = m - mult * se, ymax = m + mult * se)
  }

  # Build a single grouping column for plotting when interaction requested
  if (grepl(":", which_factor)) {
    facs <- strsplit(which_factor, ":")[[1]]
    data <- dplyr::mutate(data, .x = interaction(!!!rlang::syms(facs), drop = TRUE))
    x_aes <- ".x"
  } else {
    x_aes <- which_factor
  }

  # Range for letter placement
  y_range <- range(data[[.response]], na.rm = TRUE)
  y_pad   <- 0.05 * diff(y_range)

  # p-value matrix -> CLD letters (for non-emmeans paths)
  pmat_to_letters <- function(pairs_df, g1, g2, pcol, alpha = 0.05) {
    groups <- sort(unique(c(pairs_df[[g1]], pairs_df[[g2]])))
    M <- matrix(1, length(groups), length(groups), dimnames = list(groups, groups))
    for (i in seq_len(nrow(pairs_df))) {
      a <- as.character(pairs_df[[g1]][i]); b <- as.character(pairs_df[[g2]][i])
      p <- as.numeric(pairs_df[[pcol]][i]); M[a,b] <- p; M[b,a] <- p
    }
    multcompView::multcompLetters(M < alpha, compare = "<")$Letters
  }

  # -------------------------- DRY helpers (letters & posthocs) ---------------
  .build_letters_skipped <- function(data, factor_or_x) {
    if (factor_or_x == ".x") {
      tibble::tibble(.x = levels(as.factor(data$.x)), .group = "")
    } else {
      lv <- levels(as.factor(data[[factor_or_x]]))
      tibble::tibble(!!factor_or_x := lv, .group = "")
    }
  }

  .build_letters_from_emm <- function(emm, which_factor, ci_level, adjust) {
    cld <- multcomp::cld(emm, Letters = letters, adjust = adjust,
                         alpha = 1 - ci_level) |> as.data.frame()
    cld$.group <- gsub("\\s+", "", cld$.group)
    if (grepl(":", which_factor)) {
      facs <- strsplit(which_factor, ":")[[1]]
      cld$.x <- interaction(cld[, facs], drop = TRUE)
      cld[, c(".x", ".group")]
    } else {
      cld[, c(which_factor, ".group")]
    }
  }

  .build_letters_from_pmat <- function(posthoc_tbl, g1, g2, pcol, alpha, factor_name) {
    letters_vec <- pmat_to_letters(posthoc_tbl, g1, g2, pcol, alpha = alpha)
    tibble::tibble(!!factor_name := names(letters_vec), .group = unname(letters_vec))
  }

  # -------------------------- Normality check (AUTO) -------------------------
  pick_norm_method <- function(n) {
    if (normality_method != "auto") return(normality_method)
    if (n <= 2000) "shapiro" else if (n <= 50000) "ad" else "jb"
  }
  run_norm_test <- function(vec) {
    n <- length(vec)
    meth <- pick_norm_method(n)
    if (meth == "none") return(list(p = 1, note = "none"))
    if (meth == "shapiro") {
      if (n <= shapiro_limit) {
        return(list(p = stats::shapiro.test(vec)$p.value, note = "shapiro"))
      } else {
        set.seed(seed); idx <- sample.int(n, shapiro_limit)
        return(list(p = stats::shapiro.test(vec[idx])$p.value, note = "shapiro_subsample"))
      }
    }
    if (meth == "subsample") {
      set.seed(seed); idx <- if (n > shapiro_limit) sample.int(n, shapiro_limit) else seq_len(n)
      return(list(p = stats::shapiro.test(vec[idx])$p.value, note = "shapiro_subsample"))
    }
    if (meth == "ad")  return(list(p = nortest::ad.test(vec)$p.value,     note = "anderson_darling"))
    if (meth == "jb")  { jb <- tseries::jarque.bera.test(vec); return(list(p = jb$p.value, note = "jarque_bera")) }
  }

  normality_ok <- TRUE
  normality_note <- "auto"

  if (k_terms == 1) {
    tmp_fit <- stats::lm(formula, data = data)
    res <- stats::residuals(tmp_fit)
    if (length(res) >= 3) {
      r <- run_norm_test(res); normality_ok <- (r$p > 0.05); normality_note <- r$note
    }
  } else {
    # Factorial grouping: use only *real* columns for per-cell checks (no "A:B")
    facs_for_cells <- setdiff(all.vars(stats::terms(formula)), .response)
    facs_for_cells <- intersect(facs_for_cells, names(data))
    grp <- interaction(data[, facs_for_cells, drop = FALSE], drop = TRUE)
    normality_ok <- TRUE; normality_note <- "cells_auto"
    for (lev in levels(grp)) {
      y <- data[[.response]][grp == lev]
      if (length(y) >= 3) {
        r <- run_norm_test(y)
        if (r$p <= 0.05) normality_ok <- FALSE
        normality_note <- paste0("cells_", r$note)
      }
    }
  }

  # -------------------------- Homogeneity (Levene) ---------------------------
  if (k_terms == 1) {
    data[[trms[1]]] <- as.factor(data[[trms[1]]])   # avoid warning
    lev <- car::leveneTest(formula, data = data)
  } else {
    facs <- setdiff(all.vars(stats::terms(formula)), .response)
    facs <- intersect(facs, names(data))
    data$.levgrp <- as.factor(interaction(data[, facs, drop = FALSE], drop = TRUE))
    lev <- car::leveneTest(stats::reformulate(".levgrp", response = .response), data = data)
  }
  homogeneity_ok <- lev$`Pr(>F)`[1] > 0.05

  # -------------------------- Welch preference heuristic ---------------------
  welch_preference_applied <- FALSE
  welch_preference_reason  <- NULL

  normality_ok_eff <- normality_ok
  if (k_terms == 1 && prefer_welch_when_large_n && !homogeneity_ok) {
    grp_sizes <- table(as.factor(data[[trms[1]]]))
    if (length(grp_sizes) > 1 && all(grp_sizes >= min_n_for_welch)) {
      normality_ok_eff <- TRUE
      welch_preference_applied <- TRUE
      nmin <- min(grp_sizes)
      welch_preference_reason <- paste0(
        "Welch preferred over Kruskal because Levene's test failed (heteroscedastic) ",
        "and all group sizes ≥ ", min_n_for_welch, " (min n = ", nmin, ")."
      )
      normality_note <- paste0(normality_note, "_welch_pref")
    }
  }

  # -------------------------- Decide test ------------------------------------
  test_chosen <- NULL
  fit <- NULL
  anova_tbl <- NULL
  emm <- NULL
  posthoc_tbl <- NULL
  letters_df <- NULL
  letter_source <- NULL

  # how many groups are we lettering?
  n_levels_target <- function() {
    if (k_terms == 1) {
      length(levels(as.factor(data[[trms[1]]])))
    } else if (grepl(":", which_factor)) {
      facs <- strsplit(which_factor, ":")[[1]]
      length(levels(interaction(data[, facs], drop = TRUE)))
    } else {
      length(levels(as.factor(data[[which_factor]])))
    }
  }
  n_groups <- n_levels_target()
  skip_letters <- is.finite(max_groups_for_letters) && (n_groups > max_groups_for_letters)

  if (k_terms == 1) {
    # -------- ONE-WAY --------
    grp_var <- trms[1]

    if (normality_ok_eff && homogeneity_ok) {
      test_chosen <- "oneway_anova"
      fit <- stats::aov(formula, data = data)
      anova_tbl <- broom::tidy(stats::anova(fit))

      if (skip_letters) {
        message("Many groups (", n_groups, "). Skipping letters/posthocs for speed.")
        letters_df <- .build_letters_skipped(data, grp_var)
        letter_source <- "skipped_large_k"; posthoc_tbl <- NULL; emm <- NULL
      } else {
        emm <- emmeans::emmeans(fit, specs = stats::as.formula(paste0("~", grp_var)))
        pw  <- emmeans::contrast(emm, method = "pairwise", adjust = pairwise_adjust)
        posthoc_tbl <- broom::tidy(pw)
        letters_df <- .build_letters_from_emm(emm, grp_var, ci_level, pairwise_adjust)
        letter_source <- paste0("emmeans_cld_", pairwise_adjust)
      }

    } else if (normality_ok_eff && !homogeneity_ok) {
      test_chosen <- "oneway_welch_gameshowell"
      welch <- stats::oneway.test(formula, data = data, var.equal = FALSE)
      anova_tbl <- tibble::tibble(
        term = grp_var,
        statistic = unname(welch$statistic),
        df = unname(welch$parameter),
        p.value = unname(welch$p.value),
        method = "Welch ANOVA"
      )

      if (skip_letters) {
        message("Many groups (", n_groups, "). Skipping letters/posthocs for speed.")
        letters_df <- .build_letters_skipped(data, grp_var)
        letter_source <- "skipped_large_k"; posthoc_tbl <- NULL; emm <- NULL
      } else {
        gh <- rstatix::games_howell_test(data, stats::as.formula(paste(.response, "~", grp_var)))
        posthoc_tbl <- dplyr::rename(gh, group1 = group1, group2 = group2, p.adj = p.adj)
        letters_df <- .build_letters_from_pmat(posthoc_tbl, "group1", "group2", "p.adj",
                                               alpha_letters, grp_var)
        letter_source <- "games_howell_letters"
      }

    } else {
      test_chosen <- "kruskal_dunn"
      kw <- rstatix::kruskal_test(data, stats::as.formula(paste(.response, "~", grp_var)))
      anova_tbl <- tibble::tibble(
        term = grp_var,
        statistic = kw$statistic[1],
        p.value  = kw$p[1],
        method   = "Kruskal–Wallis"
      )

      if (skip_letters) {
        message("Many groups (", n_groups, "). Skipping letters/posthocs for speed.")
        letters_df <- .build_letters_skipped(data, grp_var)
        letter_source <- "skipped_large_k"; posthoc_tbl <- NULL; emm <- NULL
      } else {
        dunn <- rstatix::dunn_test(
          data, stats::as.formula(paste(.response, "~", grp_var)),
          p.adjust.method = "BH"
        )
        posthoc_tbl <- dplyr::transmute(dunn, group1 = group1, group2 = group2, p.adj = p.adj)
        letters_df <- .build_letters_from_pmat(posthoc_tbl, "group1", "group2", "p.adj",
                                               alpha_letters, grp_var)
        letter_source <- "dunn_letters"
      }
    }

  } else {
    # -------- TWO-WAY (or higher) --------
    if (normality_ok && homogeneity_ok) {
      test_chosen <- "factorial_anova"
      fit <- stats::aov(formula, data = data)
      anova_tbl <- broom::tidy(stats::anova(fit))

      if (skip_letters) {
        message("Many groups (", n_groups, "). Skipping letters/posthocs for speed.")
        letters_df <- if (grepl(":", which_factor)) .build_letters_skipped(data, ".x")
        else .build_letters_skipped(data, which_factor)
        letter_source <- "skipped_large_k"; posthoc_tbl <- NULL; emm <- NULL
      } else {
        emm <- emmeans::emmeans(fit, specs = stats::as.formula(paste0("~", which_factor)))
        pw  <- emmeans::contrast(emm, method = "pairwise", adjust = pairwise_adjust)
        posthoc_tbl <- broom::tidy(pw)
        letters_df <- .build_letters_from_emm(emm, which_factor, ci_level, pairwise_adjust)
        letter_source <- paste0("emmeans_cld_", pairwise_adjust)
      }

    } else {
      # ----- ART path (use artlm() for emmeans on term of interest) -----
      test_chosen <- "factorial_art"
      art_fit <- ARTool::art(formula, data = data)
      fit <- art_fit

      art_aov <- stats::anova(art_fit)                       # robust tidy for ARTool
      anova_tbl <- tibble::as_tibble(art_aov, rownames = "term")

      if (skip_letters) {
        message("Many groups (", n_groups, "). Skipping letters/posthocs for speed.")
        letters_df <- if (grepl(":", which_factor)) .build_letters_skipped(data, ".x")
        else .build_letters_skipped(data, which_factor)
        letter_source <- "skipped_large_k"; posthoc_tbl <- NULL; emm <- NULL
      } else {
        # emmeans on ART: must use artlm() for the term of interest
        lm_term <- ARTool::artlm(art_fit, which_factor)

        emm <- emmeans::emmeans(lm_term, specs = stats::as.formula(paste0("~", which_factor)))
        pw  <- emmeans::contrast(emm, method = "pairwise", adjust = pairwise_adjust)
        posthoc_tbl <- broom::tidy(pw)

        letters_df <- .build_letters_from_emm(emm, which_factor, ci_level, pairwise_adjust)
        letter_source <- paste0("art_emmeans_cld_", pairwise_adjust)
      }
    }
  }

  # -------------------------- Build plot + letters ----------------------------
  if (!grepl(":", which_factor)) data$.x <- data[[x_aes]]

  y_positions <- dplyr::group_by(data, .x) |>
    dplyr::summarize(ymax = max(.data[[.response]], na.rm = TRUE), .groups = "drop")

  letters_pos <- if (".x" %in% names(letters_df)) {
    dplyr::rename(letters_df, x = .x)
  } else {
    dplyr::rename(letters_df, x = dplyr::all_of(names(letters_df)[1]))
  }
  letters_pos <- dplyr::left_join(letters_pos, y_positions, by = c("x" = ".x")) |>
    dplyr::mutate(y = ymax + y_pad)

  # Subtitle styles (glue if available)
  .subtitle_full <- function() {
    if (requireNamespace("glue", quietly = TRUE)) {
      glue::glue(
        "Auto: {test_chosen}  |  Normality: {ifelse(normality_ok,'OK','fail')} ({normality_note})  |  ",
        "Homogeneity: {ifelse(homogeneity_ok,'OK','fail')}  |  Letters: {letter_source}",
        "{if (welch_preference_applied) paste0('  |  Decision: ', welch_preference_reason) else ''}"
      )
    } else {
      paste0(
        "Auto: ", test_chosen,
        "  |  Normality: ", ifelse(normality_ok,"OK","fail"), " (", normality_note, ")",
        "  |  Homogeneity: ", ifelse(homogeneity_ok,"OK","fail"),
        "  |  Letters: ", letter_source,
        if (welch_preference_applied) paste0("  |  Decision: ", welch_preference_reason) else ""
      )
    }
  }
  .subtitle_compact <- function() {
    if (requireNamespace("glue", quietly = TRUE)) {
      glue::glue(
        "{test_chosen} • Norm {ifelse(normality_ok,'OK','fail')} ({normality_note}) • ",
        "Lev {ifelse(homogeneity_ok,'OK','fail')}",
        "{if (welch_preference_applied) glue::glue(' • Welch: {welch_preference_reason}') else ''}"
      )
    } else {
      paste0(
        test_chosen, " • Norm ", ifelse(normality_ok,'OK','fail'), " (", normality_note, ") • ",
        "Lev ", ifelse(homogeneity_ok,'OK','fail'),
        if (welch_preference_applied) paste0(" • Welch: ", welch_preference_reason) else ""
      )
    }
  }
  .subtitle_minimal <- function() test_chosen

  subtitle_text <- switch(
    subtitle_style,
    full    = .subtitle_full(),
    compact = .subtitle_compact(),
    minimal = .subtitle_minimal()
  )

  p <- ggplot2::ggplot(data, ggplot2::aes(x = .x, y = .data[[.response]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6) +
    ggplot2::geom_jitter(width = jitter_width, alpha = point_alpha, size = 1.6) +
    {
      if (show_means == "point") {
        list(ggplot2::stat_summary(fun = "mean", geom = "point", size = 2.6, shape = 21, fill = "black"))
      } else if (show_means == "point+ci") {
        list(
          ggplot2::stat_summary(fun = "mean", geom = "point", size = 2.6, shape = 21, fill = "black"),
          ggplot2::stat_summary(fun.data = mean_cl_normal_local, geom = "errorbar", width = 0.2)
        )
      } else NULL
    } +
    ggplot2::geom_text(data = letters_pos,
                       ggplot2::aes(x = x, y = y, label = .group),
                       vjust = 0, size = 5) +
    ggplot2::labs(
      x = which_factor,
      y = .response,
      subtitle = subtitle_text
    ) +
    theme_base +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(face = "italic"))

  # -------------------------- Return -----------------------------------------
  list(
    test_info = list(
      test = test_chosen,
      normality_ok = normality_ok,
      normality_note = normality_note,
      homogeneity_ok = homogeneity_ok,
      levene = lev,
      n_groups_for_letters = n_groups,
      letters_skipped = skip_letters,
      welch_preference_applied = welch_preference_applied,
      welch_preference_reason  = welch_preference_reason
    ),
    model = fit,                # aov / art model; NULL for Welch/KW branches
    anova_table = anova_tbl,
    emmeans = emm,              # may be NULL when letters skipped or Welch/KW branch
    posthoc = posthoc_tbl,
    letters = letters_df,
    plot = p
  )
}
