#' BioinfoTools: General-Purpose Tools for Bulk RNA-Seq Analysis
#'
#' A collection of general-purpose R functions for bulk RNA-seq analysis, covering
#' transcription factor and pathway activity inference, differential expression,
#' gene set enrichment, survival analysis, normalization, batch correction, PCA, and
#' score-vs-trait association testing. See the package README for a full function
#' reference and a worked usage example.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

# Column/variable names referenced via non-standard evaluation (dplyr pipelines, ggplot2
# aes(), decoupleR's .data/.env pronouns) rather than as actual global variables/functions.
utils::globalVariables(c(
  ".data", ".env", "Group", "PC1", "PC2", "Regulator", "Sample", "TF", "Target", "Trait",
  "Value", "condition", "direction", "gene_symbol", "gs_name", "mor", "p_value", "score",
  "target", "tf", "weight"
))
