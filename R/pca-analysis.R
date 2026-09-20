#' PCA analysis and visualization for a samples x features matrix
#'
#' Runs a principal component analysis on any samples x features matrix, optionally colored by
#' a grouping column from \code{coldata}, and saves a PC1 vs. PC2 scatter plot to the
#' "Results/" folder.
#'
#' @param data A samples x features numeric matrix or data frame. Zero-variance and
#'   all-NA columns are dropped automatically.
#' @param coldata Optional data frame of sample-level annotations, with row names matching
#'   the row names of \code{data}. Used to color samples in the PCA plot.
#' @param group_col Optional character. Column name in \code{coldata} used to color/group
#'   samples in the PCA plot. Ignored if \code{coldata} is \code{NULL}.
#' @param center Logical. If TRUE (default), center features before PCA (see \code{stats::prcomp}).
#' @param scale Logical. If TRUE (default), scale features to unit variance before PCA.
#' @param file.name Optional character suffix used when saving the PCA plot to \code{Results/}.
#' @param return Logical. If TRUE (default), saves the PCA plot as a PDF in \code{Results/}.
#'
#' @return A list with:
#' \describe{
#'   \item{pca}{The \code{stats::prcomp} object.}
#'   \item{coords}{A data frame with the PC1/PC2 coordinates, \code{Sample}, and \code{Group} columns.}
#'   \item{var_explained}{Percentage of variance explained by PC1 and PC2.}
#'   \item{plot}{The \code{ggplot2} PCA scatter plot object.}
#' }
#' @export
#'
compute.pca.analysis <- function(data, coldata = NULL, group_col = NULL,
                                 center = TRUE, scale = TRUE,
                                 file.name = NULL, return = TRUE) {

  m <- as.matrix(data)
  mode(m) <- "numeric"

  # Drop zero-variance / all-NA features, as needed for PCA to run
  keep <- apply(m, 2, function(x) all(!is.na(x)) && stats::var(x, na.rm = TRUE) != 0)
  m <- m[, keep, drop = FALSE]

  if (nrow(m) < 2 || ncol(m) < 2) {
    stop("Need at least 2 samples and 2 numeric features (with non-zero variance) for PCA.")
  }

  pca <- stats::prcomp(m, center = center, scale. = scale)
  pcs <- as.data.frame(pca$x[, 1:2, drop = FALSE], check.names = FALSE)
  pcs$Sample <- rownames(m)

  if (!is.null(coldata) && !is.null(group_col)) {
    cdf <- as.data.frame(coldata, check.names = FALSE)
    idx <- match(pcs$Sample, rownames(cdf))
    pcs$Group <- as.character(cdf[[group_col]][idx])
    pcs$Group[is.na(pcs$Group)] <- "Unknown"
  } else {
    pcs$Group <- "Samples"
  }

  var_exp <- summary(pca)$importance[2, 1:2] * 100
  names(var_exp) <- c("PC1", "PC2")

  p <- ggplot2::ggplot(pcs, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::geom_text(ggplot2::aes(label = Sample), vjust = -0.7, size = 3, show.legend = FALSE) +
    ggplot2::labs(
      x     = paste0("PC1 (", round(var_exp[["PC1"]], 1), "%)"),
      y     = paste0("PC2 (", round(var_exp[["PC2"]], 1), "%)"),
      color = if (!is.null(group_col)) group_col else "Group"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  if (return) {
    grDevices::pdf(paste0("Results/PCA_plot", if (!is.null(file.name)) paste0("_", file.name), ".pdf"),
        width = 7, height = 6)
    print(p)
    grDevices::dev.off()
  }

  return(list(pca = pca, coords = pcs, var_explained = var_exp, plot = p))
}
