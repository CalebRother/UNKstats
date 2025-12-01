#' Two-way factorial ANOVA (Grouped-Boxplot Version)
#'
#' Fits a 2x(2+) between-subjects factorial ANOVA (via afex).
#' Performs post-hoc pairwise comparisons via `emmeans`, produces compact letter displays,
#' and visualizes results with a grouped boxplot + jitter.
#'
#' @param data A data frame containing the variables to be analyzed.
#' @param dv Character; name of the numeric dependent variable.
#' @param factor_a Character; first factor (main effect A, used for x-axis grouping).
#' @param factor_b Character; second factor (main effect B, used for color fill).
#' @param include_interaction Logical; if TRUE includes A*B interaction term in the ANOVA table.
#'   (Note: The underlying model in afex fits the interaction to calculate Type III SS correctly).
#' @param which_factor Character; which term to perform post-hoc tests on
#'   (e.g., "Dose", "Temp", or "Dose:Temp"). Defaults to `factor_a` if unspecified.
#' @param adjust Adjustment method for pairwise comparisons
#'   ("tukey","sidak","holm","bonferroni","BH").
#' @param show_means "point" or "none" — controls whether to show mean points on top of boxes.
#' @param theme_base A ggplot2 theme (default = `ggplot2::theme_bw()`).
#'
#' @return A list containing:
#' \itemize{
#'   \item `test_info` — basic metadata on the test.
#'   \item `model` — the afex ANOVA model object.
#'   \item `anova_table` — ANOVA summary table (Type III SS).
#'   \item `posthoc` — post-hoc pairwise comparison results.
#'   \item `letters` — compact letter display results.
#'   \item `plot` — grouped boxplot (ggplot object).
#' }
#' Class `"teach_anova_result"`.
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
  
  adjust     <- match.arg(adjust)
  show_means <- match.arg(show_means)
  
  # --- Data Prep -------------------------------------------------------------
  df <- data[stats::complete.cases(data[, c(dv, factor_a, factor_b), drop = FALSE]), , drop = FALSE]
  
  if (nrow(df) < nrow(data)) {
    message("Removed ", nrow(data) - nrow(df), " row(s) with missing values.")
  }
  
  if (!is.numeric(df[[dv]])) stop("`dv` must be numeric.", call. = FALSE)
  df[[factor_a]] <- as.factor(df[[factor_a]])
  df[[factor_b]] <- as.factor(df[[factor_b]])
  
  # Dummy ID for afex (between-subjects: one row per subject)
  df$id_dummy <- factor(seq_len(nrow(df)))
  
  # --- ANOVA (Type III SS via afex) ------------------------------------------
  afex_mod <- afex::aov_ez(
    id          = "id_dummy",
    dv          = dv,
    data        = df,
    between     = c(factor_a, factor_b),
    type        = 3,
    include_aov = TRUE,
    return      = "afex_aov"
  )
  
  an_tbl <- afex::nice(afex_mod, es = "pes")
  an_tbl <- as.data.frame(an_tbl)
  
  if (!include_interaction) {
    an_tbl <- an_tbl[!grepl(":", an_tbl$Effect), ]
  }
  
  # --- Post-hoc Analysis -----------------------------------------------------
  if (is.null(which_factor)) which_factor <- factor_a
  
  emm <- try(
    emmeans::emmeans(afex_mod, specs = stats::as.formula(paste0("~", which_factor))),
    silent = TRUE
  )
  
  posthoc <- NULL
  if (!inherits(emm, "try-error")) {
    pw      <- emmeans::contrast(emm, method = "pairwise", adjust = adjust)
    posthoc <- broom::tidy(pw)
  }
  
  # --- Compact Letter Display (CLD) — OLD multcomp style ---------------------
  letters_df <- NULL
  if (!inherits(emm, "try-error")) {
    cld_res <- try(
      multcomp::cld(emm, Letters = letters, adjust = adjust),
      silent = TRUE
    )
    
    if (!inherits(cld_res, "try-error")) {
      cld_res <- as.data.frame(cld_res)
      cld_res$.group <- gsub("\\s+", "", cld_res$.group)
      
      if (grepl(":", which_factor)) {
        # interaction term, e.g. "Dose:Temp"
        facs <- strsplit(which_factor, ":")[[1]]
        cld_res$combo <- interaction(cld_res[, facs, drop = FALSE], drop = TRUE)
        letters_df <- cld_res[, c("combo", ".group")]
      } else {
        # main effect
        letters_df <- cld_res[, c(which_factor, ".group")]
        names(letters_df)[1] <- which_factor
      }
    }
  }
  
  # --- Summary Data (Means + Max for Labeling) -------------------------------
  summary_df <- df |>
    dplyr::group_by(.data[[factor_a]], .data[[factor_b]]) |>
    dplyr::summarise(
      mean    = mean(.data[[dv]], na.rm = TRUE),
      max_val = max(.data[[dv]],  na.rm = TRUE),
      se      = stats::sd(.data[[dv]], na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    )
  
  y_range <- diff(range(df[[dv]], na.rm = TRUE))
  y_pad   <- 0.05 * y_range
  
  # --- Attach Letters --------------------------------------------------------
  if (!is.null(letters_df)) {
    if ("combo" %in% names(letters_df)) {
      summary_df <- summary_df |>
        dplyr::mutate(combo = interaction(.data[[factor_a]], .data[[factor_b]], drop = TRUE))
      label_df <- dplyr::left_join(summary_df, letters_df, by = "combo")
    } else {
      join_var <- ifelse(which_factor == factor_a, factor_a, factor_b)
      label_df <- dplyr::left_join(summary_df, letters_df, by = join_var)
    }
    
    label_df <- label_df |>
      dplyr::mutate(y = max_val + y_pad)
    
  } else {
    label_df <- NULL
  }
  
  # --- Plot (Box + Jitter, Grouped) ------------------------------------------
  pd <- ggplot2::position_dodge(width = 0.7)
  
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[[factor_a]], y = .data[[dv]], fill = .data[[factor_b]])
  ) +
    ggplot2::geom_boxplot(
      position      = pd,
      outlier.shape = NA,
      width         = 0.6
    ) +
    ggplot2::geom_jitter(
      position = ggplot2::position_jitterdodge(
        jitter.width = 0.12,
        dodge.width  = 0.7
      ),
      alpha = 0.5,
      size  = 1.4
    )
  
  if (!is.null(label_df)) {
    p <- p +
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(
          x     = .data[[factor_a]],
          y     = y,
          label = .group,
          group = .data[[factor_b]]
        ),
        position = pd,
        vjust    = 0,
        size     = 5
      )
  }
  
  if (show_means == "point") {
    p <- p +
      ggplot2::geom_point(
        data = summary_df,
        ggplot2::aes(
          x     = .data[[factor_a]],
          y     = mean,
          group = .data[[factor_b]]
        ),
        position = pd,
        color    = "black",
        size     = 2
      )
  }
  
  p <- p +
    theme_base +
    ggplot2::labs(
      x = factor_a,
      y = dv,
      fill = factor_b,
      subtitle = paste0(
        "Two-way ANOVA (Type III SS)",
        if (include_interaction) " with Interaction" else " Main Effects Only"
      )
    )
  
  # --- Return ----------------------------------------------------------------
  out <- list(
    test_info = list(
      test                = "two_way_anova",
      parametric          = TRUE,
      include_interaction = include_interaction,
      which_factor        = which_factor,
      ss_type             = "Type III"
    ),
    model       = afex_mod,
    anova_table = an_tbl,
    posthoc     = posthoc,
    letters     = letters_df,
    plot        = p + theme(axis.text.x = element_text(angle=90, 
                                                       vjust=.5, 
                                                       hjust=1))
  )
  
  class(out) <- c("teach_anova_result", class(out))
  attr(out, "meta") <- list(
    design      = "two_way_between",
    dv          = dv,
    factors     = list(A = factor_a, B = factor_b),
    adjust      = adjust,
    interaction = include_interaction
  )
  
  out
}
