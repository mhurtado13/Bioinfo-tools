#' Pearson and Spearman correlation analysis
#'
#' General-purpose correlation utility that reports both Pearson (parametric, linear) and
#' Spearman (non-parametric, monotonic) correlation together, instead of picking one method.
#' Reuses \code{stats::cor.test()} the same way \code{scores.pearson.test()}/
#' \code{scores.spearman.test()} already do for the vector and scores-vs-trait cases, and
#' \code{WGCNA::cor()}/\code{WGCNA::corPvalueStudent()} for the matrix-vs-matrix case. Handles
#' three input shapes:
#' \enumerate{
#'   \item Two numeric vectors \code{x} and \code{y} - returns a small data frame with both
#'     correlation coefficients and p-values.
#'   \item A score matrix \code{x} plus \code{coldata}/\code{trait} - calls the existing
#'     \code{scores.pearson.test()} and \code{scores.spearman.test()} on every column, unions
#'     the features significant by either method, and returns their subset plus a combined
#'     per-feature statistics table.
#'   \item Two matrices/data frames \code{x} and \code{y} (samples x features each, same
#'     samples in the same order) - correlates every column of \code{x} against every column
#'     of \code{y} with both methods and, if \code{file_name} is given, saves a heatmap per
#'     method to \code{Results/}.
#' }
#'
#' @param x A numeric vector, a score matrix (samples x features), or a
#'   samples x features matrix/data frame to correlate against \code{y}.
#' @param y Optional. A numeric vector (case 1) or a samples x features matrix/data frame
#'   (case 3) to correlate against \code{x}.
#' @param coldata Optional data frame of sample annotations (case 2), used together with \code{trait}.
#' @param trait Optional character. Column name of the continuous trait in \code{coldata} (case 2).
#' @param file_name Optional character. Base file name (without extension) used to save
#'   heatmaps to \code{Results/} in the matrix-vs-matrix case (case 3).
#' @param pval Numeric. P-value threshold for significance (default = 0.05).
#'
#' @return
#' \itemize{
#'   \item Case 1: a data frame with \code{method}, \code{estimate}, and \code{p.value} for
#'     Pearson and Spearman.
#'   \item Case 2: a list with \code{$scores} (the significant score subset, as returned by
#'     \code{scores.pearson.test()}) and \code{$stats} (a data frame of per-feature
#'     Pearson/Spearman estimates and p-values). Returns \code{NULL} (with a message) if no
#'     feature is significant by either method.
#'   \item Case 3: a list with \code{$pearson} and \code{$spearman}, each containing
#'     \code{$cor} and \code{$pval} matrices (features of \code{x} x features of \code{y}).
#' }
#' @export
#'
compute.correlation <- function(x, y = NULL, coldata = NULL, trait = NULL,
                                file_name = NULL, pval = 0.05) {

  ## ---------- Case 1: two plain numeric vectors ----------
  if (is.null(dim(x)) && is.numeric(x) && !is.null(y) && is.null(dim(y)) && is.numeric(y)) {

    res.pearson  <- stats::cor.test(x, y, method = "pearson")
    res.spearman <- stats::cor.test(x, y, method = "spearman", exact = FALSE)

    return(data.frame(
      method    = c("pearson", "spearman"),
      estimate  = c(unname(res.pearson$estimate), unname(res.spearman$estimate)),
      p.value   = c(res.pearson$p.value, res.spearman$p.value)
    ))
  }

  ## ---------- Case 2: score matrix vs a continuous trait ----------
  if (!is.null(coldata) && !is.null(trait)) {

    coldata[, trait] <- as.numeric(coldata[, trait])

    message("Running Pearson and Spearman correlation for score comparison...\n")

    sig.pearson  <- scores.pearson.test(x, coldata, trait, pval)
    sig.spearman <- scores.spearman.test(x, coldata, trait, pval)

    sig_names <- union(
      if (!is.null(sig.pearson))  colnames(sig.pearson)  else character(0),
      if (!is.null(sig.spearman)) colnames(sig.spearman) else character(0)
    )

    # Combined per-feature correlation statistics (both methods, every feature)
    stats_tab <- do.call(rbind, lapply(colnames(x), function(feat) {
      data <- data.frame(Value = x[, feat], Trait = coldata[, trait])
      data <- data[stats::complete.cases(data), ]
      rp <- stats::cor.test(data$Value, data$Trait, method = "pearson")
      rs <- stats::cor.test(data$Value, data$Trait, method = "spearman", exact = FALSE)
      data.frame(feature      = feat,
                pearson.r    = unname(rp$estimate),  pearson.p  = rp$p.value,
                spearman.rho = unname(rs$estimate),  spearman.p = rs$p.value)
    }))

    if (length(sig_names) == 0) {
      message("No significant features (p-value < ", pval, ") after Pearson or Spearman correlation.")
      return(NULL)
    }

    return(list(
      scores = x[, sig_names, drop = FALSE],
      stats  = stats_tab[stats_tab$feature %in% sig_names, ]
    ))
  }

  ## ---------- Case 3: two matrices/data frames ----------
  if (!is.null(y) && (is.matrix(y) || is.data.frame(y))) {

    matA <- data.frame(x, check.names = FALSE)
    matB <- data.frame(y, check.names = FALSE)

    if (length(rownames(matA)) == 0 || !all(rownames(matA) == rownames(matB)))
      stop("'x' and 'y' must share the same rownames (samples in the same order).")

    cor.pearson   <- WGCNA::cor(matA, matB, method = "p")
    pval.pearson  <- WGCNA::corPvalueStudent(cor.pearson, nrow(matA))
    cor.spearman  <- WGCNA::cor(matA, matB, method = "s")
    pval.spearman <- WGCNA::corPvalueStudent(cor.spearman, nrow(matA))

    if (!is.null(file_name)) {
      for (m in c("pearson", "spearman")) {
        cor_mat <- if (m == "pearson") cor.pearson  else cor.spearman
        p_mat   <- if (m == "pearson") pval.pearson else pval.spearman

        textMatrix <- paste0(signif(cor_mat, 2), "\n(", signif(p_mat, 2), ")")
        dim(textMatrix) <- dim(cor_mat)

        grDevices::pdf(paste0("Results/Correlation_", m, "_", file_name, ".pdf"), width = 8, height = 8)
        graphics::par(mar = c(15, 15, 3, 3))
        WGCNA::labeledHeatmap(Matrix = cor_mat,
                              xLabels = colnames(cor_mat),
                              yLabels = rownames(cor_mat),
                              xLabelsPosition = "top",
                              colors = WGCNA::blueWhiteRed(50),
                              textMatrix = textMatrix,
                              setStdMargins = FALSE,
                              cex.text = 0.5,
                              zlim = c(-1, 1))
        grDevices::dev.off()
      }
    }

    return(list(
      pearson  = list(cor = cor.pearson,  pval = pval.pearson),
      spearman = list(cor = cor.spearman, pval = pval.spearman)
    ))
  }

  stop("Provide either two numeric vectors ('x','y'), a score matrix 'x' with 'coldata'/'trait', or two matrices 'x' and 'y'.")
}
