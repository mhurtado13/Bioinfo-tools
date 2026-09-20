#' T-test for a categorical trait association (2 groups)
#'
#' Performs Welch's t-test (unequal variances, the \code{stats::t.test()} default) between each
#' score column and a categorical clinical/experimental variable with exactly two groups.
#' Significant features are plotted as boxplots and saved as PDF files in the "Results/" directory.
#'
#' @param scores A samples x features score matrix.
#' @param coldata A data frame containing sample-level annotations including the trait to test.
#' @param trait Character. Name of the column in `coldata` used as the categorical variable;
#'   must have exactly two groups.
#' @param pval Numeric. P-value threshold for significance (default = 0.05).
#'
#' @return A matrix of only the significant features (columns) after the t-test.
#'         Returns \code{NULL} if no significant features are found.
#' @export
#'
scores.ttest.test <- function(scores, coldata, trait, pval = 0.05) {
  sig = c()
  coldata[, trait] = droplevels(factor(coldata[, trait]))
  if (nlevels(coldata[, trait]) != 2) {
    stop("'trait' must have exactly 2 groups for a t-test; found ", nlevels(coldata[, trait]),
        ". Use scores.anova.test() for 3+ groups.")
  }

  for (j in 1:ncol(scores)) {
    data = data.frame(Value = scores[, j], Trait = coldata[, trait])
    data = data[stats::complete.cases(data), ]

    res.t <- stats::t.test(Value ~ Trait, data = data)

    if (round(res.t$p.value, 5) <= pval) {
      cat("Significant p-value after t-test for", colnames(scores)[j], "\n")

      grDevices::pdf(paste0("Results/Ttest_", trait, "_", colnames(scores)[j], ".pdf"),
          width = 7, height = 6)
      print(
        ggplot2::ggplot(data, ggplot2::aes(x = Trait, y = Value, fill = Trait)) +
          ggplot2::geom_boxplot(alpha = 0.6, outlier.shape = NA) +
          ggplot2::geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "gray20") +
          ggpubr::stat_compare_means(method = "t.test", label.y.npc = "top", size = 4.5) +
          ggplot2::scale_fill_brewer(palette = "Set2") +
          ggplot2::labs(
            title    = colnames(scores)[j],
            subtitle = paste0("T-test | ", trait),
            x        = trait,
            y        = "Score"
          ) +
          ggplot2::theme_classic(base_size = 13) +
          ggplot2::theme(
            legend.position  = "none",
            axis.text        = ggplot2::element_text(size = 11, color = "black"),
            axis.title       = ggplot2::element_text(size = 13),
            plot.title       = ggplot2::element_text(size = 14, face = "bold"),
            plot.subtitle    = ggplot2::element_text(size = 10, color = "gray45"),
            panel.grid.major = ggplot2::element_line(color = "gray92", linewidth = 0.4)
          )
      )
      grDevices::dev.off()

      sig = c(sig, j)
    }
  }

  if (length(sig) == 0) {
    message("No significant features (p-value < ", pval, ") after t-test.")
    return(NULL)
  } else {
    return(scores[, sig, drop = FALSE])
  }
}

#' Wilcoxon rank-sum test for a categorical trait association (2 groups)
#'
#' Performs a Wilcoxon rank-sum test (Mann-Whitney U) between each score column and a categorical
#' clinical/experimental variable with exactly two groups -- a non-parametric alternative to
#' \code{scores.ttest.test()} that does not assume normally distributed scores. Significant
#' features are plotted as boxplots and saved to the "Results/" folder.
#'
#' @param scores A samples x features score matrix.
#' @param coldata A data frame containing sample-level annotations including the trait to test.
#' @param trait Character. Name of the column in `coldata` used as the categorical variable;
#'   must have exactly two groups.
#' @param pval Numeric. P-value threshold for significance (default = 0.05).
#'
#' @return A matrix of only the significant features (columns) after the Wilcoxon test.
#'         Returns \code{NULL} if no significant features are found.
#' @export
#'
scores.wilcox.test <- function(scores, coldata, trait, pval = 0.05) {
  sig = c()
  coldata[, trait] = droplevels(factor(coldata[, trait]))
  if (nlevels(coldata[, trait]) != 2) {
    stop("'trait' must have exactly 2 groups for a Wilcoxon test; found ", nlevels(coldata[, trait]),
        ". Use scores.kruskal.test() for 3+ groups.")
  }

  for (j in 1:ncol(scores)) {
    data = data.frame(Value = scores[, j], Trait = coldata[, trait])
    data = data[stats::complete.cases(data), ]

    res.w <- stats::wilcox.test(Value ~ Trait, data = data)

    if (round(res.w$p.value, 5) <= pval) {
      cat("Significant p-value after Wilcoxon test for", colnames(scores)[j], "\n")

      grDevices::pdf(paste0("Results/Wilcoxon_", trait, "_", colnames(scores)[j], ".pdf"),
          width = 7, height = 6)
      print(
        ggplot2::ggplot(data, ggplot2::aes(x = Trait, y = Value, fill = Trait)) +
          ggplot2::geom_boxplot(alpha = 0.6, outlier.shape = NA) +
          ggplot2::geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "gray20") +
          ggpubr::stat_compare_means(method = "wilcox.test", label.y.npc = "top", size = 4.5) +
          ggplot2::scale_fill_brewer(palette = "Set2") +
          ggplot2::labs(
            title    = colnames(scores)[j],
            subtitle = paste0("Wilcoxon test | ", trait),
            x        = trait,
            y        = "Score"
          ) +
          ggplot2::theme_classic(base_size = 13) +
          ggplot2::theme(
            legend.position  = "none",
            axis.text        = ggplot2::element_text(size = 11, color = "black"),
            axis.title       = ggplot2::element_text(size = 13),
            plot.title       = ggplot2::element_text(size = 14, face = "bold"),
            plot.subtitle    = ggplot2::element_text(size = 10, color = "gray45"),
            panel.grid.major = ggplot2::element_line(color = "gray92", linewidth = 0.4)
          )
      )
      grDevices::dev.off()

      sig = c(sig, j)
    }
  }

  if (length(sig) == 0) {
    message("No significant features (p-value < ", pval, ") after Wilcoxon test.")
    return(NULL)
  } else {
    return(scores[, sig, drop = FALSE])
  }
}

#' ANOVA for a categorical trait association (2+ groups)
#'
#' Performs a one-way ANOVA F-test (\code{stats::aov()}) between each score column and a
#' categorical clinical/experimental variable with two or more groups. Significant features are
#' plotted as boxplots and saved to the "Results/" folder.
#'
#' @param scores A samples x features score matrix.
#' @param coldata A data frame containing sample-level annotations including the trait to test.
#' @param trait Character. Name of the column in `coldata` used as the categorical variable
#'   (two or more groups).
#' @param pval Numeric. P-value threshold for significance (default = 0.05).
#'
#' @return A matrix of only the significant features (columns) after the ANOVA F-test.
#'         Returns \code{NULL} if no significant features are found.
#' @export
#'
scores.anova.test <- function(scores, coldata, trait, pval = 0.05) {
  sig = c()
  coldata[, trait] = droplevels(factor(coldata[, trait]))
  if (nlevels(coldata[, trait]) < 2) {
    stop("'trait' must have at least 2 groups for ANOVA; found ", nlevels(coldata[, trait]), ".")
  }

  for (j in 1:ncol(scores)) {
    data = data.frame(Value = scores[, j], Trait = coldata[, trait])
    data = data[stats::complete.cases(data), ]

    res.aov <- stats::aov(Value ~ Trait, data = data)
    aov.pval <- summary(res.aov)[[1]][["Pr(>F)"]][1]

    if (round(aov.pval, 5) <= pval) {
      cat("Significant p-value after ANOVA for", colnames(scores)[j], "\n")

      grDevices::pdf(paste0("Results/ANOVA_", trait, "_", colnames(scores)[j], ".pdf"),
          width = 7, height = 6)
      print(
        ggplot2::ggplot(data, ggplot2::aes(x = Trait, y = Value, fill = Trait)) +
          ggplot2::geom_boxplot(alpha = 0.6, outlier.shape = NA) +
          ggplot2::geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "gray20") +
          ggpubr::stat_compare_means(method = "anova", label.y.npc = "top", size = 4.5) +
          ggplot2::scale_fill_brewer(palette = "Set2") +
          ggplot2::labs(
            title    = colnames(scores)[j],
            subtitle = paste0("ANOVA | ", trait),
            x        = trait,
            y        = "Score"
          ) +
          ggplot2::theme_classic(base_size = 13) +
          ggplot2::theme(
            legend.position  = "none",
            axis.text        = ggplot2::element_text(size = 11, color = "black"),
            axis.title       = ggplot2::element_text(size = 13),
            plot.title       = ggplot2::element_text(size = 14, face = "bold"),
            plot.subtitle    = ggplot2::element_text(size = 10, color = "gray45"),
            panel.grid.major = ggplot2::element_line(color = "gray92", linewidth = 0.4)
          )
      )
      grDevices::dev.off()

      sig = c(sig, j)
    }
  }

  if (length(sig) == 0) {
    message("No significant features (p-value < ", pval, ") after ANOVA.")
    return(NULL)
  } else {
    return(scores[, sig, drop = FALSE])
  }
}

#' Kruskal-Wallis test for a categorical trait association (2+ groups)
#'
#' Performs a Kruskal-Wallis rank-sum test (\code{stats::kruskal.test()}) between each score
#' column and a categorical clinical/experimental variable with two or more groups -- a
#' non-parametric alternative to \code{scores.anova.test()} that does not assume normally
#' distributed scores. Significant features are plotted as boxplots and saved to the "Results/"
#' folder.
#'
#' @param scores A samples x features score matrix.
#' @param coldata A data frame containing sample-level annotations including the trait to test.
#' @param trait Character. Name of the column in `coldata` used as the categorical variable
#'   (two or more groups).
#' @param pval Numeric. P-value threshold for significance (default = 0.05).
#'
#' @return A matrix of only the significant features (columns) after the Kruskal-Wallis test.
#'         Returns \code{NULL} if no significant features are found.
#' @export
#'
scores.kruskal.test <- function(scores, coldata, trait, pval = 0.05) {
  sig = c()
  coldata[, trait] = droplevels(factor(coldata[, trait]))
  if (nlevels(coldata[, trait]) < 2) {
    stop("'trait' must have at least 2 groups for a Kruskal-Wallis test; found ", nlevels(coldata[, trait]), ".")
  }

  for (j in 1:ncol(scores)) {
    data = data.frame(Value = scores[, j], Trait = coldata[, trait])
    data = data[stats::complete.cases(data), ]

    res.kw <- stats::kruskal.test(Value ~ Trait, data = data)

    if (round(res.kw$p.value, 5) <= pval) {
      cat("Significant p-value after Kruskal-Wallis test for", colnames(scores)[j], "\n")

      grDevices::pdf(paste0("Results/Kruskal_", trait, "_", colnames(scores)[j], ".pdf"),
          width = 7, height = 6)
      print(
        ggplot2::ggplot(data, ggplot2::aes(x = Trait, y = Value, fill = Trait)) +
          ggplot2::geom_boxplot(alpha = 0.6, outlier.shape = NA) +
          ggplot2::geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "gray20") +
          ggpubr::stat_compare_means(method = "kruskal.test", label.y.npc = "top", size = 4.5) +
          ggplot2::scale_fill_brewer(palette = "Set2") +
          ggplot2::labs(
            title    = colnames(scores)[j],
            subtitle = paste0("Kruskal-Wallis test | ", trait),
            x        = trait,
            y        = "Score"
          ) +
          ggplot2::theme_classic(base_size = 13) +
          ggplot2::theme(
            legend.position  = "none",
            axis.text        = ggplot2::element_text(size = 11, color = "black"),
            axis.title       = ggplot2::element_text(size = 13),
            plot.title       = ggplot2::element_text(size = 14, face = "bold"),
            plot.subtitle    = ggplot2::element_text(size = 10, color = "gray45"),
            panel.grid.major = ggplot2::element_line(color = "gray92", linewidth = 0.4)
          )
      )
      grDevices::dev.off()

      sig = c(sig, j)
    }
  }

  if (length(sig) == 0) {
    message("No significant features (p-value < ", pval, ") after Kruskal-Wallis test.")
    return(NULL)
  } else {
    return(scores[, sig, drop = FALSE])
  }
}

#' Perform statistical analysis on scores against a categorical trait
#'
#' Unified interface for testing associations between any score matrix and a categorical
#' clinical/experimental trait, mirroring \code{scores.stat.analysis.continuous()} for the
#' categorical case.
#'
#' @param scores A samples x features score matrix.
#' @param coldata A data frame containing clinical or experimental metadata for samples.
#'   Must include the column specified in `trait`.
#' @param trait Character. The name of the column in `coldata` representing the categorical
#'   clinical or experimental trait to test against (e.g., treatment arm, response group,
#'   risk class).
#' @param method Character. Statistical test to perform. One of:
#'   \itemize{
#'     \item `"ttest"` - Welch's t-test (parametric, exactly 2 groups)
#'     \item `"wilcox"` - Wilcoxon rank-sum test (non-parametric, exactly 2 groups)
#'     \item `"anova"` - One-way ANOVA F-test (parametric, 2+ groups)
#'     \item `"kruskal"` - Kruskal-Wallis test (non-parametric, 2+ groups)
#'   }
#'   Only one method is used per call; defaults to `"ttest"` if not specified.
#' @param pval Numeric. P-value threshold for significance (default: 0.05).
#'
#' @details
#' The function automatically calls the corresponding statistical test function
#' based on the `method` argument:
#' \itemize{
#'   \item \code{scores.ttest.test()}
#'   \item \code{scores.wilcox.test()}
#'   \item \code{scores.anova.test()}
#'   \item \code{scores.kruskal.test()}
#' }
#'
#' Each test produces both a statistical result and visual outputs (PDF boxplots)
#' stored in the `"Results/"` folder.
#'
#' @return A matrix of only the significant features (columns) after the chosen test.
#'   Returns `NULL` if no significant features are found.
#'
#' @examples
#' \dontrun{
#' sig <- scores.stat.analysis.categorical(score_matrix, coldata, trait = "response_group",
#'                                         method = "kruskal", pval = 0.05)
#' }
#'
#' @export
#'
scores.stat.analysis.categorical <- function(scores, coldata, trait,
                                             method = c("ttest", "wilcox", "anova", "kruskal"),
                                             pval = 0.05) {
  method <- match.arg(method)

  message("Running ", toupper(method), " test for score comparison...\n")

  result <- switch(method,
                   ttest   = scores.ttest.test(scores, coldata, trait, pval),
                   wilcox  = scores.wilcox.test(scores, coldata, trait, pval),
                   anova   = scores.anova.test(scores, coldata, trait, pval),
                   kruskal = scores.kruskal.test(scores, coldata, trait, pval))

  if (is.null(result)) {
    message("None feature was found significant after ", method, " test (p < ", pval, ").")
  } else {
    message("Significant features found: ", ncol(result))
  }

  return(result)
}
