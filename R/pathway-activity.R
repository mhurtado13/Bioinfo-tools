#' Compute sample-level pathway activity scores
#'
#' Computes pathway activity scores from normalized gene expression data. By default, the pathway
#' resource is PROGENy (Schubert et al., 2018), but KEGG, REACTOME, or MSigDB Hallmark of Cancer gene
#' sets can be used instead (`pathway_source`). Any of these can be scored by any
#' \code{decoupleR::decouple()} statistic (`statistic`), default `"consensus"`. Optionally, it also
#' scores any user-provided gene sets with that same statistic, independently of `pathway_source`.
#'
#' @param RNA.tpm A numeric matrix of log-normalized gene expression values -- log2(TPM + pseudocount), not raw
#'   or linear-scale TPM -- with genes as rows and samples as columns (e.g. the output of
#'   \code{compute.normalization(..., log = TRUE, zscore = FALSE)}).
#' @param gene_sets A list of gene sets (e.g., hallmark signatures or user-defined sets). If provided, activity scores are
#'   computed for these sets (with `statistic`, like `pathway_source`) and saved to Results/ if `return` is TRUE, but are
#'   not part of the value returned by this function -- only the `pathway_source` matrix is returned. Default is \code{NULL}.
#' @param paths A data frame describing the pathway-gene interactions for use with PROGENy (ignored for other `pathway_source` values). If \code{NULL}, the human PROGENy resource (top 500 genes) will be used by default.
#' @param pathway_source Character. The pathway resource used for the primary pathway score. Options are `"PROGENy"` (default, signed weights),
#'   `"KEGG"`, `"REACTOME"`, or `"Hallmark"` (all three unsigned MSigDB gene set collections).
#' @param statistic Character. Which \code{decoupleR::decouple()} statistic to score `pathway_source` (and `gene_sets`, if given) with. Default `"consensus"`
#'   (an ensemble score across the top-performing methods: MLM, ULM, and normalized WSUM). Any single underlying method
#'   can be requested instead: `"aucell"`, `"udt"`, `"mdt"`, `"wmean"`, `"ulm"`, `"mlm"`, `"wsum"`, `"viper"`, `"gsva"`,
#'   `"ora"`, or `"fgsea"` (see \code{decoupleR::show_methods()}). `"mlm"` (and therefore the default `"consensus"`) fits
#'   all sources jointly and can fail with a collinearity error on large, highly redundant collections (e.g. REACTOME) --
#'   use a per-source method (e.g. `"ulm"`, `"viper"`, `"wsum"`) in that case.
#' @param min_targets_size Integer. Minimum number of target genes per source (pathway or gene set) present in `RNA.tpm`
#'   required for it to be scored; sources with fewer are dropped. Default is 5.
#' @param return Logical; if TRUE, saves matrices in Results/ folder. Default is TRUE.
#' @param file.name Optional character suffix used when writing output CSV files.
#'
#' @return A scaled matrix of pathway activity scores for `pathway_source` (samples as rows, pathways as columns).
#' If `gene_sets` was also provided, its activity scores are saved to Results/ (when `return` is TRUE) but are not
#' included in this return value.
#'
#'
#' @references
#' Schubert M, Klinger B, Kluenemann M, Sieber A, Uhlitz F, Sauer S, Garnett MJ, Bluethgen N, Saez-Rodriguez J.
#' Perturbation-response genes reveal signaling footprints in cancer gene expression. Nature Communications. 2018. \doi{10.1038/s41467-017-02391-6}
#'
#' Badia-i-Mompel P., Vélez Santiago J., Braunger J., Geiss C., Dimitrov D., Müller-Dott S., Taus P., Dugourd A., Holland C.H., Ramirez Flores R.O. and Saez-Rodriguez J. 2022. decoupleR: ensemble of computational methods to infer biological activities from omics data. Bioinformatics Advances. https://doi.org/10.1093/bioadv/vbac016
#' 
#' @examples
#' \dontrun{
#' # Compute only PROGENy activities
#' pathways <- compute.pathway.activity(counts.norm)
#'
#' # Use REACTOME instead of PROGENy
#' pathways_reactome <- compute.pathway.activity(counts.norm, pathway_source = "REACTOME")
#' }
#'
#' @export
#'
compute.pathway.activity <- function(RNA.tpm, gene_sets = NULL, paths = NULL,
                                     pathway_source = c("PROGENy", "KEGG", "REACTOME", "Hallmark"),
                                     statistic = c("consensus", "aucell", "udt", "mdt", "wmean", "ulm", "mlm", "wsum", "viper", "gsva", "ora", "fgsea"),
                                     min_targets_size = 5, return = TRUE, file.name = NULL){

  pathway_source <- match.arg(pathway_source)
  statistic <- match.arg(statistic)
  rn <- rownames(RNA.tpm)
  RNA.tpm <- apply(as.matrix(RNA.tpm), 2, as.numeric)
  rownames(RNA.tpm) <- rn
  results_list <- list()

  # Scores a source/target(/mor) network with the requested decoupleR statistic
  score_with_decoupleR <- function(net) {
    decoupleR::decouple(
        mat             = RNA.tpm,
        network         = net,
        .source         = "source",
        .target         = "target",
        statistics      = if (statistic == "consensus") NULL else statistic,
        consensus_score = (statistic == "consensus"),
        minsize         = min_targets_size
      ) %>%
      dplyr::filter(.data$statistic == .env$statistic) %>%
      tidyr::pivot_wider(id_cols = "condition", names_from = "source", values_from = "score") %>%
      tibble::column_to_rownames("condition") %>%
      as.matrix() %>%
      scale() %>%
      as.data.frame()
  }

  if (pathway_source == "PROGENy") {

    ###### PROGENy (signed weights)
    progeny_cache_file <- "Results/Pathways_collection_PROGENy.csv"
    if (is.null(paths)) {
      if (file.exists(progeny_cache_file)) {
        paths <- utils::read.csv(progeny_cache_file, row.names = 1)
        cat("Using cached PROGENy pathways collection from ", progeny_cache_file, "\n")
      } else {
        paths <- decoupleR::get_progeny(organism = "human", top = 500)
        utils::write.csv(paths, progeny_cache_file)
      }
    }
    pathway_net <- dplyr::rename(paths, mor = weight)

  } else {

    ###### KEGG / REACTOME / Hallmark (unsigned MSigDB gene set collection)
    msigdb_cache_file <- paste0("Results/Pathways_collection_", pathway_source, ".csv")
    if (file.exists(msigdb_cache_file)) {
      msigdb_sets <- utils::read.csv(msigdb_cache_file)
      cat("Using cached", pathway_source, "gene set collection from", msigdb_cache_file, "\n")
    } else {
      msigdb_cat <- if (pathway_source == "Hallmark") "H" else "C2"
      msigdb_subcat <- switch(pathway_source, KEGG = "CP:KEGG", REACTOME = "CP:REACTOME", Hallmark = "")
      msigdb_sets <- msigdbr::msigdbr(species = "Homo sapiens", category = msigdb_cat, subcategory = msigdb_subcat)
      utils::write.csv(msigdb_sets, msigdb_cache_file, row.names = FALSE)
    }

    pathway_net <- msigdb_sets %>%
      dplyr::transmute(source = gs_name, target = gene_symbol, mor = 1) %>%
      dplyr::distinct(source, target, .keep_all = TRUE)
  }

  ###### Score pathway_net with the requested decoupleR statistic
  results_list[[pathway_source]] <- score_with_decoupleR(pathway_net)

  ###### User-provided gene sets (optional), scored with the same statistic
  if (!is.null(gene_sets)) {
    gene_sets_net <- dplyr::bind_rows(lapply(names(gene_sets), function(nm) {
        data.frame(source = nm, target = gene_sets[[nm]], mor = 1)
      })) %>%
      dplyr::distinct(source, target, .keep_all = TRUE)

    results_list$gene_sets <- score_with_decoupleR(gene_sets_net)
  }

  ###### Save outputs if requested
  if (return) {
    if (!is.null(results_list[[pathway_source]])) {
      utils::write.csv(results_list[[pathway_source]], paste0("Results/Pathway_matrix_", pathway_source, "_", file.name, ".csv"))
    }
    if (!is.null(results_list$gene_sets)) {
      utils::write.csv(results_list$gene_sets, paste0("Results/Pathway_matrix_gene_sets_", file.name,".csv"))
    }
  }

  ###### Output
  sample_acts_pathway <- results_list[[pathway_source]]
  colnames(sample_acts_pathway) <- make.names(colnames(sample_acts_pathway))

  return(sample_acts_pathway)
}
