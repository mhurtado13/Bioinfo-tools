#' Normalize raw RNA-seq counts to log-TPM and z-score scale
#'
#' Converts a raw count matrix to TPM (\code{ADImpute::NormalizeTPM()}), log-transforms it, and
#' then optionally centers/scales each gene (row) to mean 0 and unit variance across samples.
#'
#' @param raw.counts A raw gene expression count matrix, with genes as rows and samples as columns.
#' @param log Logical. If TRUE (default), log-transform the TPM values (as in \code{ADImpute::NormalizeTPM}).
#' @param zscore Logical. If TRUE, z-score scale each gene (row) across samples after
#'   TPM/log normalization. Default \code{FALSE}, since \code{compute.TFs.activity()} and
#'   \code{compute.pathway.activity()} expect log-TPM expression, not per-gene standardized values.
#' @param tr_length Optional. Transcript length information passed through to
#'   \code{ADImpute::NormalizeTPM()} (see its documentation). Default \code{NULL}.
#'
#' @return A data frame of normalized expression values (genes x samples), log-TPM and/or
#'   z-score scaled according to \code{log}/\code{zscore}.
#' @export
#'
compute.normalization <- function(raw.counts, log = TRUE, zscore = FALSE, tr_length = NULL) {

  norm.counts <- data.frame(ADImpute::NormalizeTPM(raw.counts, tr_length = tr_length, log = log),
                            check.names = FALSE)

  if (zscore) {
    # Center/scale each gene (row) across samples, preserving dimnames
    norm.counts <- as.data.frame(t(scale(t(norm.counts))), check.names = FALSE)
  }

  return(norm.counts)
}

#' Remove batch effects from raw RNA-seq counts (ComBat-seq)
#'
#' Corrects raw counts across batches/cohorts with \code{sva::ComBat_seq()} - a negative-binomial
#' model appropriate for raw RNA-seq counts - before TPM normalization:
#' \preformatted{
#' corrected_counts <- sva::ComBat_seq(as.matrix(raw_counts), batch = batch_labels)
#' counts.normalized <- ADImpute::NormalizeTPM(corrected_counts, log = TRUE)
#' }
#' This function wraps that same two-step pattern, optionally chaining into
#' \code{compute.normalization()} (log-TPM, no z-score by default) for the corrected counts.
#'
#' @param raw.counts A raw gene expression count matrix, with genes as rows and samples as columns.
#' @param batch A vector of batch/cohort labels, one per sample, in the same order as the
#'   columns of \code{raw.counts} (e.g. \code{coldata$Batch} or \code{coldata$Cohort}).
#' @param group Optional. A vector of the biological condition of interest (e.g.
#'   treatment/response), passed through to \code{sva::ComBat_seq()} so that variation of
#'   interest is preserved during correction.
#' @param covar_mod Optional. A model matrix of additional covariates to preserve, passed
#'   through to \code{sva::ComBat_seq()}.
#' @param normalize Logical. If TRUE (default), the batch-corrected counts are normalized
#'   with \code{compute.normalization()} (log-TPM, no z-score by default) before being
#'   returned, matching the ComBat_seq -> NormalizeTPM pattern shown above.
#' @param log Passed to \code{compute.normalization()} when \code{normalize = TRUE}.
#' @param zscore Passed to \code{compute.normalization()} when \code{normalize = TRUE}.
#'   Default \code{FALSE}, since \code{compute.TFs.activity()} and \code{compute.pathway.activity()}
#'   expect log-TPM expression, not per-gene standardized values.
#'
#' @return If \code{normalize = FALSE}, a data frame of batch-corrected raw counts
#'   (genes x samples). If \code{normalize = TRUE} (default), a data frame of
#'   batch-corrected, log-TPM (and optionally z-score) normalized expression values.
#'
#' @examples
#' \dontrun{
#' corrected <- compute.batch.correction(raw.counts, batch = coldata$Batch)
#' }
#'
#' @export
#'
compute.batch.correction <- function(raw.counts, batch, group = NULL, covar_mod = NULL,
                                     normalize = TRUE, log = TRUE, zscore = FALSE) {

  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("Package 'sva' is required for compute.batch.correction().")
  }

  if (length(batch) != ncol(raw.counts)) {
    stop("Length of 'batch' must match the number of samples (columns) in 'raw.counts'.")
  }

  corrected.counts <- sva::ComBat_seq(as.matrix(raw.counts), batch = batch,
                                      group = group, covar_mod = covar_mod)

  if (normalize) {
    return(compute.normalization(corrected.counts, log = log, zscore = zscore))
  }

  return(data.frame(corrected.counts, check.names = FALSE))
}
