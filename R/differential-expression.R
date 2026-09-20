#' Run differential expression analysis with edgeR/limma-voom
#'
#' Filters low-expression genes, applies TMM normalization, runs voom
#' transformation, fits a linear model, and returns the top differentially
#' expressed genes via \code{limma::topTable}. Optionally saves a volcano
#' plot (via \code{compute.volcano.plot()}) built from the full, unfiltered
#' set of tested genes.
#'
#' @param counts A raw count matrix (genes x samples).
#' @param coldata A data frame of sample metadata whose row names match
#'   the column names of \code{counts}.
#' @param group_col Character. Name of the column in \code{coldata} used as
#'   the grouping factor for differential expression.
#' @param ref_level Character or \code{NULL}. Reference level for the group
#'   factor. If \code{NULL}, the default factor ordering is used.
#' @param pval Numeric. Adjusted p-value cutoff used both to filter the
#'   returned table (via \code{topTable(p.value = pval)}) and to mark
#'   significance on the volcano plot. Default 0.05.
#' @param logFC Numeric. Absolute log2 fold-change cutoff marked on the
#'   volcano plot (does not filter the returned table). Default 1.
#' @param plot Logical; if TRUE (default), saves a volcano plot to Results/
#'   via \code{compute.volcano.plot()}.
#' @param file.name Optional character suffix used when saving the volcano
#'   plot to \code{Results/}.
#'
#' @return A data frame of differentially expressed genes (p.adj < \code{pval}) as
#'   returned by \code{limma::topTable}, with columns \code{logFC},
#'   \code{AveExpr}, \code{t}, \code{P.Value}, \code{adj.P.Val}, and \code{B}.
#'
#' @examples
#' \dontrun{
#' deg <- run_deg_analysis(counts, coldata, group_col = "Group", ref_level = "Control")
#' }
#'
#' @export
#'
run_deg_analysis <- function(counts, coldata, group_col, ref_level = NULL, pval = 0.05,
                             logFC = 1, plot = TRUE, file.name = NULL) {
  # Prepare counts
  counts_mat <- as.matrix(counts)
  mode(counts_mat) <- "numeric"
  counts_mat <- counts_mat[, rownames(coldata)]

  # Create group factor
  group <- factor(coldata[[group_col]])

  # Set reference level if provided
  if (!is.null(ref_level)) {
    group <- stats::relevel(group, ref = ref_level)
  }

  # Create DGE object and filter
  dge <- edgeR::DGEList(counts = counts_mat, group = group)
  keep <- edgeR::filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)

  # Design matrix
  design <- stats::model.matrix(~ group)

  # voom transformation
  v <- limma::voom(dge, design)

  # Fit model
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  # Extract coefficient name (second column of design)
  coef_name <- colnames(design)[2]

  # Full, unfiltered results -- used for the volcano plot so it reflects every tested gene,
  # not just the significant ones
  res_full <- limma::topTable(fit, coef = coef_name, p.value = 1, number = Inf)

  if (plot) {
    compute.volcano.plot(res_full, pval = pval, logFC = logFC,
                         title = paste0(group_col, ": ", coef_name),
                         file.name = file.name)
  }

  # Get results
  res <- res_full[res_full$adj.P.Val <= pval, ]

  return(res)
}

#' Volcano plot for differential expression results
#'
#' Plots log2 fold-change against adjusted p-value for a differential expression result (e.g. the
#' full, unfiltered output of \code{run_deg_analysis()}), highlighting genes passing both the
#' \code{pval} and \code{logFC} cutoffs, via \code{EnhancedVolcano::EnhancedVolcano()}.
#'
#' @param deg A data frame of differential expression results with a \code{logFC} column and an
#'   \code{adj.P.Val} column (e.g. \code{limma::topTable()} output), gene symbols as row names.
#'   For a meaningful volcano plot this should be the full, unfiltered set of tested genes, not
#'   only the significant ones.
#' @param pval Numeric. Adjusted p-value cutoff marked on the plot. Default 0.05.
#' @param logFC Numeric. Absolute log2 fold-change cutoff marked on the plot. Default 1.
#' @param title Character. Plot title. Default \code{NULL} (uses \code{"Volcano plot"}).
#' @param xlim Numeric vector of length 2, or \code{NULL} (default). X-axis (log2 fold-change)
#'   display range, forwarded to \code{EnhancedVolcano::EnhancedVolcano()}. If \code{NULL},
#'   \code{EnhancedVolcano}'s own default is used (data range padded by a fixed +/-1.5), which can
#'   dwarf small effect sizes and compress every point near zero. Set this explicitly (e.g.
#'   \code{c(-0.2, 0.2)}) to zoom into the actual spread of \code{logFC} in \code{deg}.
#' @param ylim Numeric vector of length 2, or \code{NULL} (default). Y-axis
#'   (\eqn{-log10} adjusted p-value) display range, forwarded to \code{EnhancedVolcano}. If
#'   \code{NULL}, \code{EnhancedVolcano}'s own default is used (\code{0} to the tallest point + 5).
#' @param top_n_labels Integer or \code{NULL}. Number of genes to label with their gene symbol,
#'   chosen as the smallest \code{adj.P.Val} (ties broken by \code{P.Value}) -- independent of
#'   \code{pval}/\code{logFC}, so gene names are still shown even when nothing clears those
#'   cutoffs (as with an underpowered comparison). Default 10. Set to \code{NULL} or \code{0} to
#'   fall back to \code{EnhancedVolcano}'s own default of labeling only genes passing both
#'   \code{pval} and \code{logFC}.
#' @param return Logical; if TRUE (default), saves the plot as a PDF in Results/.
#' @param file.name Optional character suffix used when saving the plot to \code{Results/}.
#'
#' @return The \code{ggplot} object produced by \code{EnhancedVolcano::EnhancedVolcano()}.
#'
#' @examples
#' \dontrun{
#' deg_full <- run_deg_analysis(counts, coldata, group_col = "Group", plot = FALSE)
#' compute.volcano.plot(deg_full, pval = 0.05, logFC = 1, file.name = "Group")
#'
#' # Zoom the axes to a panel with small effect sizes instead of EnhancedVolcano's
#' # default +/-1.5 log2FC padding
#' compute.volcano.plot(deg_full, pval = 0.05, logFC = 0.05,
#'                       xlim = c(-0.2, 0.2), ylim = c(0, 2), file.name = "Group_zoom")
#' }
#'
#' @export
#'
compute.volcano.plot <- function(deg, pval = 0.05, logFC = 1, title = NULL, xlim = NULL, ylim = NULL,
                                 top_n_labels = 10, return = TRUE, file.name = NULL) {

  # Always label the top genes by adjusted p-value (ties broken by raw P-value), independent of
  # pCutoff/FCcutoff: EnhancedVolcano's own default only labels genes passing BOTH cutoffs, which
  # renders with no gene names at all for an underpowered comparison where nothing is significant.
  select_lab <- NULL
  if (!is.null(top_n_labels) && top_n_labels > 0) {
    ord <- order(deg$adj.P.Val, deg$P.Value)
    select_lab <- rownames(deg)[utils::head(ord, top_n_labels)]
  }

  # legend labels are overridden because EnhancedVolcano's own defaults ("p-value", "p-value &
  # Log2 FC") are misleading here: this function always colors/thresholds on adj.P.Val (FDR),
  # never the raw P.Value column, since y = "adj.P.Val" below.
  volcano_args <- list(
    toptable              = deg,
    lab                    = rownames(deg),
    selectLab              = select_lab,
    drawConnectors         = TRUE,
    widthConnectors        = 0.3,
    max.overlaps           = Inf,
    maxoverlapsConnectors  = Inf,
    labSize                = 3.2,
    x            = "logFC",
    y            = "adj.P.Val",
    pCutoff      = pval,
    FCcutoff     = logFC,
    title        = if (is.null(title)) "Volcano plot" else title,
    subtitle     = NULL,
    ylab         = bquote(~-Log[10] ~ italic(adj.P.Val)),
    legendLabels = c("NS", "Log2 FC", "Adj. p-value", "Adj. p-value & Log2 FC")
  )
  if (!is.null(xlim)) volcano_args$xlim <- xlim
  if (!is.null(ylim)) volcano_args$ylim <- ylim

  p <- do.call(EnhancedVolcano::EnhancedVolcano, volcano_args)

  if (return) {
    grDevices::pdf(paste0("Results/Volcano_", file.name, ".pdf"), width = 10, height = 8)
    print(p)
    grDevices::dev.off()
  }

  return(p)
}
