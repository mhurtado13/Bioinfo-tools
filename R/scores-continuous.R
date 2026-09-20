#' Pearson correlation test for continuous trait association
#'
#' Performs a Pearson correlation test between each score column and a
#' continuous clinical/experimental variable. Significant features are plotted as scatter
#' plots with a fitted regression line and saved as PDF files in the "Results/" directory.
#'
#' @param scores A samples x features score matrix.
#' @param coldata A data frame containing sample-level annotations including the trait to test.
#' @param trait Character. Name of the (numeric) column in `coldata` used as the continuous variable.
#' @param pval Numeric. P-value threshold for significance (default = 0.05).
#'
#' @return A matrix of only the significant features (columns) after the Pearson test.
#'         Returns \code{NULL} if no significant features are found.
#' @export
#'
scores.pearson.test <- function(scores, coldata, trait, pval = 0.05) {
  sig = c()
  coldata[, trait] = as.numeric(coldata[, trait])

  for (j in 1:ncol(scores)) {
    data = data.frame(Value = scores[, j], Trait = coldata[, trait])
    data = data[stats::complete.cases(data), ]

    res.cor <- stats::cor.test(data$Value, data$Trait, method = "pearson")

    if (round(res.cor$p.value, 5) <= pval) {
      cat("Significant p-value after Pearson correlation for", colnames(scores)[j], "\n")

      grDevices::pdf(paste0("Results/Pearson_", trait, "_", colnames(scores)[j], ".pdf"),
          width = 7, height = 6)
      print(
        ggplot2::ggplot(data, ggplot2::aes(x = Trait, y = Value)) +
          ggplot2::geom_point(size = 1.8, alpha = 0.6, color = "gray20") +
          ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "#3C5488",
                               fill = "#3C5488", alpha = 0.15, linewidth = 0.9) +
          ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4.5) +
          ggplot2::labs(
            title    = colnames(scores)[j],
            subtitle = paste0("Pearson correlation | ", trait),
            x        = trait,
            y        = "Score"
          ) +
          ggplot2::theme_classic(base_size = 13) +
          ggplot2::theme(
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
    message("No significant features (p-value < ", pval, ") after Pearson correlation.")
    return(NULL)
  } else {
    return(scores[, sig, drop = FALSE])
  }
}

#' Spearman correlation test for continuous trait association
#'
#' Performs a Spearman rank correlation test between each score column and a
#' continuous clinical/experimental variable, appropriate when a monotonic (non-linear) rather
#' than strictly linear relationship is expected. Significant features are plotted as scatter
#' plots and saved to the "Results/" folder.
#'
#' @param scores A samples x features score matrix.
#' @param coldata A data frame containing sample annotations and clinical traits.
#' @param trait Character. Name of the (numeric) column in `coldata` used as the continuous variable.
#' @param pval Numeric. P-value threshold for significance (default = 0.05).
#'
#' @return A matrix of only the significant features (columns), or \code{NULL} if none are found.
#' @export
#'
scores.spearman.test <- function(scores, coldata, trait, pval = 0.05) {
  sig = c()
  coldata[, trait] = as.numeric(coldata[, trait])

  for (j in 1:ncol(scores)) {
    data = data.frame(Value = scores[, j], Trait = coldata[, trait])
    data = data[stats::complete.cases(data), ]

    res.cor <- stats::cor.test(data$Value, data$Trait, method = "spearman", exact = FALSE)

    if (round(res.cor$p.value, 5) <= pval) {
      cat("Significant p-value after Spearman correlation for", colnames(scores)[j], "\n")

      grDevices::pdf(paste0("Results/Spearman_", trait, "_", colnames(scores)[j], ".pdf"),
          width = 7, height = 6)
      print(
        ggplot2::ggplot(data, ggplot2::aes(x = Trait, y = Value)) +
          ggplot2::geom_point(size = 1.8, alpha = 0.6, color = "gray20") +
          ggplot2::geom_smooth(method = "loess", formula = y ~ x, color = "#00A087",
                               fill = "#00A087", alpha = 0.15, linewidth = 0.9, se = TRUE) +
          ggpubr::stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 4.5) +
          ggplot2::labs(
            title    = colnames(scores)[j],
            subtitle = paste0("Spearman correlation | ", trait),
            x        = trait,
            y        = "Score"
          ) +
          ggplot2::theme_classic(base_size = 13) +
          ggplot2::theme(
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
    message("No significant features (p-value < ", pval, ") after Spearman correlation.")
    return(invisible(NULL))
  } else {
    return(scores[, sig, drop = FALSE])
  }
}

#' Fisher's exact test for two continuous variables (median-split)
#'
#' Both the score and the trait are binarised at their median into
#' High/Low groups, and Fisher's exact test is applied to the resulting 2x2 contingency table.
#'
#' @param scores A samples x features score matrix. Continuous scores are binarised at the
#'   median into High/Low groups.
#' @param coldata A data frame containing the clinical or experimental traits.
#' @param trait Character. Name of the (numeric) column in `coldata` to test with Fisher's
#'   exact test. It is binarised at the median into High/Low groups.
#' @param pval Numeric. P-value threshold for significance (default 0.05).
#'
#' @return A matrix of only the significant features (columns) after the Fisher test, or
#'   \code{NULL} if none are significant. Additionally, it saves corresponding barplot
#'   visualizations in the "Results/" folder.
#' @export
#'
scores.fisher.test.continuous = function(scores, coldata, trait, pval = 0.05){

  sig = c()
  coldata[, trait] = as.numeric(coldata[, trait])
  trait_level = ifelse(coldata[, trait] > stats::median(coldata[, trait], na.rm = TRUE), "High", "Low")

  for (j in 1:ncol(scores)) {
    level = scores[, j]
    score_level = ifelse(level > stats::median(level, na.rm = TRUE), "High", "Low")

    df = data.frame(score_level = score_level, trait_level = trait_level)
    df = df[stats::complete.cases(df), ]

    contingency = table(df[, "score_level"], df[, "trait_level"])
    test = stats::fisher.test(contingency)

    if (round(test$p.value, 5) <= pval) {
      cat("Significant pval after doing Fisher test for", colnames(scores)[j], "\n")
      p_label <- ifelse(test$p.value < 0.001, "< 0.001", paste0("= ", round(test$p.value, 4)))
      grDevices::pdf(paste0("Results/Fisher_", trait, "_", colnames(scores)[j], ".pdf"), width = 8, height = 6)
      print(
        ggstatsplot::ggbarstats(df, score_level, trait_level, results.subtitle = FALSE,
                                title        = colnames(scores)[j],
                                subtitle     = paste0("Fisher's exact test | p ", p_label, " | ", trait, " (median split)"),
                                xlab         = paste0(trait, " level"),
                                legend.title = "Score level") +
          ggplot2::theme(
            plot.title    = ggplot2::element_text(size = 14, face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 10, color = "gray45"),
            axis.text     = ggplot2::element_text(size = 12),
            axis.title    = ggplot2::element_text(size = 13),
            legend.title  = ggplot2::element_text(size = 12),
            legend.text   = ggplot2::element_text(size = 11)
          )
      )
      grDevices::dev.off()

      sig = c(j, sig)
    }
  }

  if (length(sig) == 0) {
    message("No significant features (pvalue < ", pval, ") after Fisher test")
    return(NULL)
  } else {
    return(scores[, sig, drop = FALSE])
  }

}

#' Perform statistical analysis on scores against a continuous trait
#'
#' Unified interface for testing associations between any score matrix and a continuous
#' (numeric) clinical/experimental trait.
#'
#' @param scores A samples x features score matrix.
#' @param coldata A data frame containing clinical or experimental metadata for samples.
#'   Must include the column specified in `trait`.
#' @param trait Character. The name of the (numeric) column in `coldata` representing the
#'   continuous clinical or experimental trait to test against (e.g., age, tumor purity,
#'   a continuous score, etc.).
#' @param method Character. Statistical test to perform. One of:
#'   \itemize{
#'     \item `"pearson"` - Pearson correlation (parametric, linear relationship)
#'     \item `"spearman"` - Spearman rank correlation (non-parametric, monotonic relationship)
#'     \item `"fisher"` - Fisher's exact test (both score and trait binarised at median)
#'   }
#'   Only one method is used per call; defaults to `"pearson"` if not specified.
#' @param pval Numeric. P-value threshold for significance (default: 0.05).
#'
#' @details
#' The function automatically calls the corresponding statistical test function
#' based on the `method` argument:
#' \itemize{
#'   \item \code{scores.pearson.test()}
#'   \item \code{scores.spearman.test()}
#'   \item \code{scores.fisher.test.continuous()}
#' }
#'
#' Each test produces both a statistical result and visual outputs (PDF plots)
#' stored in the `"Results/"` folder. These visualizations include the relevant
#' test results (p-values) annotated on the plots.
#'
#' @return A matrix of only the significant features (columns) after the chosen test.
#'   Returns `NULL` if no significant features are found.
#'
#' @examples
#' \dontrun{
#' sig <- scores.stat.analysis.continuous(score_matrix, coldata, trait = "tumor_purity",
#'                                        method = "spearman", pval = 0.05)
#' }
#'
#' @export
#'
scores.stat.analysis.continuous <- function(scores, coldata, trait,
                                            method = c("pearson", "spearman", "fisher"),
                                            pval = 0.05) {
  method <- match.arg(method)

  if (!is.numeric(coldata[, trait])) {
    stop("'trait' must be a continuous (numeric) column in 'coldata'.")
  }

  message("Running ", toupper(method), " test for score comparison...\n")

  result <- switch(method,
                   pearson  = scores.pearson.test(scores, coldata, trait, pval),
                   spearman = scores.spearman.test(scores, coldata, trait, pval),
                   fisher   = scores.fisher.test.continuous(scores, coldata, trait, pval))

  if (is.null(result)) {
    message("None feature was found significant after ", method, " test (p < ", pval, ").")
  } else {
    message("Significant features found: ", ncol(result))
  }

  return(result)
}
