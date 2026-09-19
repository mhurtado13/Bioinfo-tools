#' Compute Transcription Factor (TF) activity
#'
#' Infers transcription factor (TF) activity from a gene expression matrix using the VIPER algorithm (Alvarez et al., 2016). The function requires a TF-target gene regulatory network, which can be provided by the user or obtained from OmnipathR resources such as CollecTRI or Dorothea. ARACNE-inferred networks are also supported.
#'
#' @param RNA.counts A gene expression matrix with genes as rows and samples as columns. The matrix should be normalized (e.g., TPM, log2CPM, etc.).
#' @param TF.collection Character. The source of the TF-target network. Options are `"CollecTRI"` (default), `"Dorothea"`, or `"ARACNE"`.
#' - `"CollecTRI"` and `"Dorothea"` use prebuilt collections from OmnipathR.
#' - `"ARACNE"` allows user input of a custom network file in a 3-column format: `regulator`, `target`, and `mutual information`.
#' @param min_targets_size Integer. Minimum number of target genes per regulon required for TF activity inference. Default is 5.
#' @param universe Optional. A user-specified data frame of TF-target interactions. If not provided, the function will fetch the relevant network based on the `TF.collection` argument.
#' @param statistic Character. Which \code{decoupleR::decouple()} statistic to extract. Default `"consensus"`
#'   (an ensemble score across the top-performing methods: MLM, ULM, and normalized WSUM). Any single
#'   underlying method can be requested instead: `"aucell"`, `"udt"`, `"mdt"`, `"wmean"`, `"ulm"`, `"mlm"`,
#'   `"wsum"`, `"viper"`, `"gsva"`, `"ora"`, or `"fgsea"` (see \code{decoupleR::show_methods()}).
#' @param cores Integer. Number of cores used by VIPER inference. Default is 4.
#' @param return Logical; if TRUE, saves matrix in Results/ folder. Default is TRUE.
#' @param file.name Optional character suffix used when writing the TF activity matrix to disk.
#'
#' @return A data frame of inferred and scaled TF activity scores, with samples as rows and TFs as columns.
#'
#' @references
#' Alvarez, M. et al. (2016). Functional characterization of somatic mutations in cancer using network-based inference of protein activity. *Nature Genetics*, 48(8), 838-847. https://doi.org/10.1038/ng.3593
#'
#' Tuerei, D., Korcsmaros, T., & Saez-Rodriguez, J. (2016). OmniPath: guidelines and gateway for literature-curated signaling pathway resources. *Nature Methods*, 13(12), 966-967. https://doi.org/10.1038/nmeth.4077
#'
#' Garcia-Alonso, L. et al. (2019). Benchmark and integration of resources for the estimation of human transcription factor activities. *Genome Research*. https://doi.org/10.1101/gr.240663.118
#'
#' Lachmann, A. et al. (2016). ARACNe-AP: gene network reverse engineering through adaptive partitioning inference of mutual information. *Bioinformatics*, 32(14), 2233-2235. https://doi.org/10.1093/bioinformatics/btw216
#'
#' Margolin, A.A. et al. (2006). ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. *BMC Bioinformatics*, 7(Suppl 1), S7. https://doi.org/10.1186/1471-2105-7-S1-S7
#'
#' @examples
#' data("counts.norm.tuto")
#' tfs_activity <- compute.TFs.activity(counts.norm.tuto, cores = 1)
#'
compute.TFs.activity <- function(RNA.counts, TF.collection = "CollecTRI", min_targets_size = 5, universe = NULL,
                                 statistic = c("consensus", "aucell", "udt", "mdt", "wmean", "ulm", "mlm", "wsum", "viper", "gsva", "ora", "fgsea"),
                                 cores = 3, return = TRUE, file.name = NULL){

  statistic <- match.arg(statistic)
  tf_cache_file <- "Results/TF_target_collection.csv"

  if(TF.collection == "ARACNE"){

    # auto-discover when only one network exists
    candidates <- list.files("input/ARACNE", pattern = "^network\\.txt$",
                             recursive = TRUE, full.names = TRUE)
    if(length(candidates) == 0)
      stop("TF.collection = 'ARACNE' requires a network.txt under input/ARACNE/")
    if(length(candidates) > 1)
      stop("Multiple ARACNe networks found under input/ARACNE/:\n",
           paste(dirname(dirname(candidates)), collapse = "\n"))
    aracne.network <- candidates[1]
    cat("Auto-detected ARACNe network:", aracne.network, "\n")

    cat("Loading ARACNe network from:", aracne.network, "\n")
    # Read network edges, filter to genes present in expression matrix
    aracne_net <- utils::read.table(aracne.network, header = TRUE, sep = "\t") %>%
      dplyr::select(source = Regulator, target = Target) %>%
      dplyr::filter(source %in% rownames(RNA.counts) & target %in% rownames(RNA.counts))

    # Compute Spearman correlation between every TF and its targets in one matrix op
    # (this is exactly what TFmode1 does internally)
    all_tfs     <- unique(aracne_net$source)
    all_targets <- unique(aracne_net$target)
    cor_mat <- suppressWarnings(
      stats::cor(t(RNA.counts[all_tfs, , drop = FALSE]),
                 t(RNA.counts[all_targets, , drop = FALSE]),
                 method = "spearman")
    )

    # mor = sign of Spearman correlation - +1 activation, -1 repression
    universe <- aracne_net %>%
      dplyr::mutate(mor = cor_mat[cbind(source, target)]) %>%
      dplyr::mutate(mor = sign(mor)) %>%
      dplyr::filter(!is.na(mor) & mor != 0)

    cat("Computing TF activities...\n")
    
    sample_acts <- decoupleR::decouple( mat     = RNA.counts,
                                        network = universe,
                                        .source = "source",
                                        .target = "target",
                                        statistics = if (statistic == "consensus") NULL else statistic,
                                        consensus_score = (statistic == "consensus"),
                                        minsize = min_targets_size
                                      ) %>%
      dplyr::filter(.data$statistic == .env$statistic) %>%
      decoupleR::pivot_wider_profile(id_cols     = source,
                                     names_from  = condition,
                                     values_from = score) %>%
      as.matrix() %>%
      t()

  } else {

    if(TF.collection == "CollecTRI"){
      if(is.null(universe)){
        if(file.exists(tf_cache_file)){
          universe = utils::read.csv(tf_cache_file, row.names = 1)
          cat("Using cached TF-target collection from", tf_cache_file, "\n")
        } else {
          universe = decoupleR::get_collectri(organism = 'human', split_complexes = F)
          utils::write.csv(universe, tf_cache_file)
        }
      }
    } else if(TF.collection == "Dorothea"){
      if(is.null(universe)){
        if(file.exists(tf_cache_file)){
          universe = utils::read.csv(tf_cache_file, row.names = 1)
        } else {
          universe = dplyr::filter(dorothea::dorothea_hs, .data$confidence %in% c("A", "B")) %>%
            dplyr::mutate(source = .data$tf) %>%
            dplyr::select(-tf)
          utils::write.csv(universe, tf_cache_file)
        }
      }
    }

    sample_acts <- decoupleR::decouple(mat     = RNA.counts,
                                       network = universe,
                                      .source = "source",
                                      .target = "target",
                                      statistics = if (statistic == "consensus") NULL else statistic,
                                      consensus_score = (statistic == "consensus"),
                                    ) %>%
                                      dplyr::filter(.data$statistic == .env$statistic) %>%
                                      decoupleR::pivot_wider_profile(id_cols     = source,
                                                                     names_from  = condition,
                                                                     values_from = score) %>%
                                      as.matrix() %>%
                                      t()

  }

  sample_acts <- sample_acts[colnames(RNA.counts), , drop = FALSE]

  if(return){
    utils::write.csv(sample_acts, paste0("Results/TF_matrix_", file.name, ".csv"))
  }

  result <- data.frame(sample_acts)
  colnames(result) <- make.names(colnames(result))
  return(result)

}

#' Computes TF-modules pathway activities scores
#'
#' This function computes pathway activity scores from normalized gene expression data.
#' By default, it uses a multivariate linear model (MLM) based on the PROGENy resource (Schubert et al., 2018).
#' Alternatively, KEGG, REACTOME, or MSigDB Hallmark of Cancer gene sets can be used instead, scored via
#' Gene Set Variation Analysis (GSVA). Optionally, it also performs GSVA using any user-provided gene sets.
#'
#' @param RNA.tpm A numeric matrix of normalized gene expression values with genes as rows and samples as columns.
#' @param gene_sets A list of gene sets (e.g., hallmark signatures or user-defined sets). If provided, GSVA scores will be computed for these sets. Default is \code{NULL}.
#' @param paths A data frame describing the pathway-gene interactions for use with PROGENy (ignored for other `pathway_source` values). If \code{NULL}, the human PROGENy resource (top 500 genes) will be used by default.
#' @param pathway_source Character. The pathway resource used for the primary pathway score. Options are `"PROGENy"` (default, MLM-based),
#'   `"KEGG"`, `"REACTOME"`, or `"Hallmark"` (all three GSVA-based, using the corresponding MSigDB collection).
#' @param return Logical; if TRUE, saves matrices in Results/ folder. Default is TRUE.
#' @param file.name Optional character suffix used when writing output CSV files.
#'
#' @return If \code{gene_sets} is \code{NULL}, a scaled matrix of pathway activity scores for `pathway_source` (samples as rows, pathways as columns).
#' If \code{gene_sets} is provided, a list with two elements:
#' \itemize{
#'   \item the `pathway_source` matrix (e.g. \code{PROGENy}): A scaled matrix of pathway activity scores.
#'   \item \code{GSVA}: A scaled matrix of GSVA scores based on the provided gene sets.
#' }
#'
#'
#' @references
#' Schubert M, Klinger B, Kluenemann M, Sieber A, Uhlitz F, Sauer S, Garnett MJ, Bluethgen N, Saez-Rodriguez J.
#' Perturbation-response genes reveal signaling footprints in cancer gene expression. Nature Communications. 2018. \doi{10.1038/s41467-017-02391-6}
#'
#' @examples
#' # Compute only PROGENy activities
#' data("counts.norm.tuto")
#' pathways <- compute.pathway.activity(counts.norm.tuto)
#'
#' # Use REACTOME instead of PROGENy
#' pathways_reactome <- compute.pathway.activity(counts.norm.tuto, pathway_source = "REACTOME")
#'
compute.pathway.activity <- function(RNA.tpm, gene_sets = NULL, paths = NULL,
                                     pathway_source = c("PROGENy", "KEGG", "REACTOME", "Hallmark"),
                                     return = TRUE, file.name = NULL) {

  pathway_source <- match.arg(pathway_source)
  rn <- rownames(RNA.tpm)
  RNA.tpm <- apply(as.matrix(RNA.tpm), 2, as.numeric)
  rownames(RNA.tpm) <- rn
  results_list <- list()

  if (pathway_source == "PROGENy") {

    ###### PROGENy (MLM)
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

    progeny <- decoupleR::run_mlm(
      mat      = RNA.tpm,
      net      = paths,
      .source  = "source",
      .target  = "target",
      .mor     = "weight",
      minsize  = 5
    )

    sample_acts_pathway <- progeny %>%
      tidyr::pivot_wider(id_cols = "condition", names_from = "source", values_from = "score") %>%
      tibble::column_to_rownames("condition") %>%
      as.matrix() %>%
      scale() %>%
      as.data.frame()

  } else {

    ###### KEGG / REACTOME / Hallmark (GSVA on the corresponding MSigDB collection)
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

    pathway_gene_sets <- split(msigdb_sets$gene_symbol, msigdb_sets$gs_name)

    pathway_gsva <- GSVA::gsva(
      RNA.tpm,
      pathway_gene_sets,
      method  = "gsva",
      kcdf    = "Gaussian",
      min.sz  = 1,
      mx.diff = TRUE,
      verbose = TRUE
    )

    sample_acts_pathway <- t(pathway_gsva) %>%
      scale() %>%
      as.data.frame()
  }

  results_list[[pathway_source]] <- sample_acts_pathway

  ###### GSVA (optional)
  if (!is.null(gene_sets)) {
    gsva_results <- GSVA::gsva(
      RNA.tpm,
      gene_sets,
      method  = "gsva",
      kcdf    = "Gaussian",
      min.sz  = 1,
      mx.diff = TRUE,
      verbose = TRUE
    )

    sample_acts_gsva <- t(gsva_results) %>%
      scale() %>%
      as.data.frame()

    results_list$GSVA <- sample_acts_gsva
  }

  ###### Save outputs if requested
  if (return) {
    if (!is.null(results_list[[pathway_source]])) {
      utils::write.csv(results_list[[pathway_source]], paste0("Results/Pathway_matrix_", pathway_source, "_", file.name, ".csv"))
    }
    if (!is.null(results_list$GSVA)) {
      utils::write.csv(results_list$GSVA, paste0("Results/Pathway_matrix_GSVA_", file.name,".csv"))
    }
  }

  ###### RETURN FORMAT LOGIC
  # If only one element -> return it directly (matrix)
  # standardise column names so they are valid R identifiers (consistent with caret internals)
  results_list <- lapply(results_list, function(df) {
    colnames(df) <- make.names(colnames(df))
    df
  })

  if (length(results_list) == 1) {
    return(results_list[[1]])
  }

  return(results_list)
}

#' Run differential expression analysis with edgeR/limma-voom
#'
#' Filters low-expression genes, applies TMM normalization, runs voom
#' transformation, fits a linear model, and returns the top differentially
#' expressed genes via \code{limma::topTable}.
#'
#' @param counts A raw count matrix (genes x samples).
#' @param coldata A data frame of sample metadata whose row names match
#'   the column names of \code{counts}.
#' @param group_col Character. Name of the column in \code{coldata} used as
#'   the grouping factor for differential expression.
#' @param ref_level Character or \code{NULL}. Reference level for the group
#'   factor. If \code{NULL}, the default factor ordering is used.
#'
#' @return A data frame of differentially expressed genes (p.adj < 0.05) as
#'   returned by \code{limma::topTable}, with columns \code{logFC},
#'   \code{AveExpr}, \code{t}, \code{P.Value}, \code{adj.P.Val}, and \code{B}.
#'
#' @keywords internal
run_deg_analysis <- function(counts, coldata, group_col, ref_level = NULL, pval = 0.05) {
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

  # Get results
  res <- limma::topTable(fit, coef = coef_name, p.value = pval, number = Inf)

  return(res)
}

#' Kaplan-Meier survival analysis on clinical groups or CellTFusion features
#'
#' Performs Kaplan-Meier survival analysis and log-rank testing, either on a predefined
#' grouping variable (e.g. a clinical risk group or a supervised cell group split), or
#' automatically for each column of a feature matrix (e.g. latent factor scores from
#' \code{compute.latent_factors()}, TF module scores from \code{compute.WTCNA()}, or
#' cell group scores from \code{construct_cell_groups()}). In the latter case, samples
#' are split into "High"/"Low" groups per feature using a quantile cutoff, and only
#' features with a significant log-rank test are returned.
#'
#' @param survival.data A data frame of clinical/survival metadata, with samples as rows.
#' @param PFS Character. Column name in \code{survival.data} with the survival/follow-up time.
#' @param PFS_event Character. Column name in \code{survival.data} with the event indicator
#'   (1 = event occurred, 0 = censored).
#' @param file_name Optional character. Suffix used when saving Kaplan-Meier plots to \code{Results/}.
#' @param features Optional. A samples x features numeric matrix or data frame (e.g.
#'   \code{latent_spaces$Z}). If provided (and \code{group_column} is \code{NULL}), each
#'   feature is tested individually. Mutually exclusive with \code{group_column}.
#' @param p.value Numeric. Log-rank test p-value threshold used to keep a feature as
#'   significant when \code{features} is used. Default is 0.05.
#' @param thres Numeric between 0 and 1. Quantile cutoff used to split each feature into
#'   High/Low groups when \code{features} is used. Default is 0.5 (median split).
#' @param group_column Optional character. Column name in \code{survival.data} defining a
#'   predefined categorical grouping (e.g. cell group membership, cluster, or trait class).
#'   Mutually exclusive with \code{features}.
#'
#' @return
#' If \code{group_column} is provided, a list with:
#' \describe{
#'   \item{km_plot}{The \code{survminer::ggsurvplot} object.}
#'   \item{median_PFS}{Median survival time per group (\code{survminer::surv_median()}).}
#'   \item{p_value}{Log-rank test p-value across groups.}
#'   \item{median_follow_up}{Median follow-up time (reverse Kaplan-Meier).}
#' }
#' If \code{features} is provided, a named list of significant \code{Surv() ~ feature} formulas
#' (one per feature with log-rank p-value below \code{p.value}); \code{NULL} (with a message)
#' if none are significant.
#'
#' A Kaplan-Meier plot (with risk table) is saved as an SVG file per significant result to
#' \code{Results/SurvPlot_<group-or-feature>_<file_name>.svg}.
#'
#' @details
#' Requires the \code{survival}, \code{survminer}, and \code{gridExtra} packages (see \code{Suggests}).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Predefined clinical/cell groups
#' compute.survival.analysis(
#'   survival.data = traitdata,
#'   PFS           = "PFS",
#'   PFS_event     = "PFS_event",
#'   group_column  = "Best.Confirmed.Overall.Response",
#'   file_name     = "Tutorial"
#' )
#'
#' # Automatic screening of latent factor scores
#' compute.survival.analysis(
#'   survival.data = traitdata,
#'   PFS           = "PFS",
#'   PFS_event     = "PFS_event",
#'   features      = res$Latent_spaces$Z,
#'   p.value       = 0.05,
#'   thres         = 0.5,
#'   file_name     = "Tutorial"
#' )
#' }
compute.survival.analysis = function(survival.data, PFS, PFS_event, file_name = NULL,
                                     features = NULL, p.value = 0.05, thres = 0.5, group_column = NULL) {

  if (!requireNamespace("survival", quietly = TRUE) ||
      !requireNamespace("survminer", quietly = TRUE) ||
      !requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Packages 'survival', 'survminer', and 'gridExtra' are required for compute.survival.analysis().")
  }

  ###### CASE 1: Use predefined groups (like patient clusters or risk groups)
  if (!is.null(group_column)) {

    data_for_model <- data.frame("time" = survival.data[, PFS],
                                 "status" = survival.data[, PFS_event])

    # Use an existing column in survival.data to define patient groups
    data_for_model$coxHL <- paste0("Group_", survival.data[[group_column]])

    # Fit Kaplan-Meier survival curves across groups
    km_fit <- survival::survfit(survival::Surv(time, status) ~ coxHL, data = data_for_model)

    # Median PFS per group
    medians <- survminer::surv_median(km_fit)

    # Median follow-up (reverse KM)
    rev_km <- survival::survfit(survival::Surv(time, 1 - status) ~ 1, data = data_for_model)
    median_follow_up <- summary(rev_km)$table["median"]

    # Log-rank test for survival difference among groups
    pval <- survminer::surv_pvalue(km_fit, data = data_for_model)$pval

    # Legend labels from strata names
    strata_names <- gsub("coxHL=", "", names(km_fit$strata))

    p <- survminer::ggsurvplot(km_fit,
                               data = data_for_model,
                               size = 1,
                               conf.int.style = "step",
                               pval = TRUE,
                               risk.table = TRUE,
                               risk.table.col = "strata",
                               legend.labs = strata_names,
                               risk.table.height = 0.3,
                               ggtheme = ggplot2::theme_grey(),
                               title = paste0("Cox PH for Surv(time, status) ~ ", group_column),
                               xlab = "Time to death/recurrence/progression")

    p$table <- p$table + ggplot2::theme(legend.position = "none")

    combined <- gridExtra::arrangeGrob(p$plot, p$table, ncol = 1, heights = c(2/3, 1/3))
    ggplot2::ggsave(
      filename = paste0("Results/SurvPlot_", group_column, "_", file_name, ".svg"),
      plot = combined,
      width = 10,
      height = 5,
      units = "in",
      device = "svg"
    )

    return(list(
      km_plot = p,
      median_PFS = medians,
      p_value = pval,
      median_follow_up = median_follow_up
    ))

  ###### CASE 2: Use features to define high/low groups automatically
  } else if (!is.null(features)) {
    significant_combinations <- list()
    contador <- 1
    n_features <- ncol(features)

    for (n in 1:n_features) {

      formula <- stats::as.formula(paste("survival::Surv(time, status) ~", colnames(features)[n]))

      data_for_model <- data.frame("time" = survival.data[, PFS],
                                   "status" = survival.data[, PFS_event])

      # Quantile cutoff used to split the feature into high vs. low groups
      quantiles <- stats::quantile(features[, n], thres)

      data_for_model$coxHL <- factor(ifelse(features[, n] >= quantiles,
                                            paste0("High_", colnames(features)[n]),
                                            paste0("Low_", colnames(features)[n])))

      km_fit <- survival::survfit(survival::Surv(time, status) ~ coxHL, data = data_for_model)

      pval <- survminer::surv_pvalue(km_fit, data = data_for_model)$pval

      if (!is.na(pval) && pval < p.value) {

        significant_combinations[[contador]] <- formula
        names(significant_combinations)[contador] <- colnames(features)[n]

        strata_names <- gsub("coxHL=", "", names(km_fit$strata))

        p <- survminer::ggsurvplot(km_fit,
                                   data = data_for_model,
                                   size = 1,
                                   palette = c("#E7B800", "#2E9FDF"),
                                   conf.int.style = "step",
                                   pval = TRUE,
                                   risk.table = TRUE,
                                   risk.table.col = "strata",
                                   legend.labs = strata_names,
                                   risk.table.height = 0.3,
                                   ggtheme = ggplot2::theme_grey(),
                                   title = paste0("Cox PH for Surv(time, status) ~", colnames(features)[n]),
                                   xlab = "Time to death/recurrence/progression")

        p$table <- p$table + ggplot2::theme(legend.position = "none")
        combined <- gridExtra::arrangeGrob(p$plot, p$table, ncol = 1, heights = c(2/3, 1/3))
        ggplot2::ggsave(
          filename = paste0("Results/SurvPlot_", names(significant_combinations)[contador], "_", file_name, ".svg"),
          plot = combined,
          width = 10,
          height = 5,
          units = "in",
          device = "svg"
        )

        contador <- contador + 1
      }
    }

    if (length(significant_combinations) == 0) {
      message("No significant combinations found.")
      return(invisible(NULL))
    } else {
      return(significant_combinations)
    }

  } else {
    stop("Please provide either 'group_column' for predefined groups or 'features' for feature-based analysis.")
  }
}

#' Pearson correlation test for continuous trait association
#'
#' Adapted from CellTFusion's \code{scores.ttest()} for a continuous (rather than binary
#' categorical) trait. Performs a Pearson correlation test between each score column and a
#' continuous clinical/experimental variable. Significant features are plotted as scatter
#' plots with a fitted regression line and saved as PDF files in the "Results/" directory.
#'
#' @param scores A list, NMF output from \code{compute.latent_factors()}, or a score matrix.
#'   When a list, the first element must be a samples x features score matrix.
#' @param coldata A data frame containing sample-level annotations including the trait to test.
#' @param trait Character. Name of the (numeric) column in `coldata` used as the continuous variable.
#' @param pval Numeric. P-value threshold for significance (default = 0.05).
#'
#' @return A list containing only significant cell groups after the Pearson test.
#'         Returns \code{NULL} if no significant groups are found.
#' @export
#'
scores.pearson.test <- function(scores, coldata, trait, pval = 0.05) {
  sig = c()
  coldata[, trait] = as.numeric(coldata[, trait])

  for (j in 1:ncol(scores[[1]])) {
    data = data.frame(Value = scores[[1]][, j], Trait = coldata[, trait])
    data = data[stats::complete.cases(data), ]

    res.cor <- stats::cor.test(data$Value, data$Trait, method = "pearson")

    if (round(res.cor$p.value, 5) <= pval) {
      cat("Significant p-value after Pearson correlation for", colnames(scores[[1]])[j], "\n")

      pdf(paste0("Results/Pearson_", trait, "_", colnames(scores[[1]])[j], ".pdf"),
          width = 7, height = 6)
      print(
        ggplot2::ggplot(data, ggplot2::aes(x = Trait, y = Value)) +
          ggplot2::geom_point(size = 1.8, alpha = 0.6, color = "gray20") +
          ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "#3C5488",
                               fill = "#3C5488", alpha = 0.15, linewidth = 0.9) +
          ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4.5) +
          ggplot2::labs(
            title    = colnames(scores[[1]])[j],
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
      dev.off()

      sig = c(sig, j)
    }
  }

  if (length(sig) == 0) {
    message("No significant features (p-value < ", pval, ") after Pearson correlation.")
    return(NULL)
  } else {
    sig.scores = list()
    sig.scores[[1]] = scores[[1]][, sig, drop = FALSE]
    if (length(scores) >= 2) sig.scores[[2]] = scores[[2]][sig]
    if (length(scores) >= 3) sig.scores[[3]] = scores[[3]][sig]
    return(sig.scores)
  }
}

#' Spearman correlation test for continuous trait association
#'
#' Adapted from CellTFusion's \code{scores.wilcox.test()} for a continuous (rather than binary
#' categorical) trait. Performs a Spearman rank correlation test between each score column and a
#' continuous clinical/experimental variable, appropriate when a monotonic (non-linear) rather
#' than strictly linear relationship is expected. Significant features are plotted as scatter
#' plots and saved to the "Results/" folder.
#'
#' @param scores A list, NMF output from \code{compute.latent_factors()}, or a score matrix.
#'   When a list, the first element must be a samples x features score matrix.
#' @param coldata A data frame containing sample annotations and clinical traits.
#' @param trait Character. Name of the (numeric) column in `coldata` used as the continuous variable.
#' @param pval Numeric. P-value threshold for significance (default = 0.05).
#'
#' @return A list containing significant features or \code{NULL} if none are found.
#' @export
#'
scores.spearman.test <- function(scores, coldata, trait, pval = 0.05) {
  sig = c()
  coldata[, trait] = as.numeric(coldata[, trait])

  for (j in 1:ncol(scores[[1]])) {
    data = data.frame(Value = scores[[1]][, j], Trait = coldata[, trait])
    data = data[stats::complete.cases(data), ]

    res.cor <- stats::cor.test(data$Value, data$Trait, method = "spearman", exact = FALSE)

    if (round(res.cor$p.value, 5) <= pval) {
      cat("Significant p-value after Spearman correlation for", colnames(scores[[1]])[j], "\n")

      pdf(paste0("Results/Spearman_", trait, "_", colnames(scores[[1]])[j], ".pdf"),
          width = 7, height = 6)
      print(
        ggplot2::ggplot(data, ggplot2::aes(x = Trait, y = Value)) +
          ggplot2::geom_point(size = 1.8, alpha = 0.6, color = "gray20") +
          ggplot2::geom_smooth(method = "loess", formula = y ~ x, color = "#00A087",
                               fill = "#00A087", alpha = 0.15, linewidth = 0.9, se = TRUE) +
          ggpubr::stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 4.5) +
          ggplot2::labs(
            title    = colnames(scores[[1]])[j],
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
      dev.off()

      sig = c(sig, j)
    }
  }

  if (length(sig) == 0) {
    message("No significant features (p-value < ", pval, ") after Spearman correlation.")
    return(invisible(NULL))
  } else {
    sig.scores = list()
    sig.scores[[1]] = scores[[1]][, sig, drop = FALSE]
    if (length(scores) >= 2) sig.scores[[2]] = scores[[2]][sig]
    if (length(scores) >= 3) sig.scores[[3]] = scores[[3]][sig]
    return(sig.scores)
  }
}

#' Kendall correlation test for continuous trait association
#'
#' Adapted from CellTFusion's \code{scores.kruskal.test()} for a continuous (rather than
#' multi-level categorical) trait. Performs a Kendall's tau rank correlation test between each
#' score column and a continuous clinical/experimental variable. Kendall's tau is a
#' non-parametric alternative to Spearman that is more robust to outliers and ties and tends
#' to be preferred for small sample sizes. Significant features are plotted as scatter plots
#' and saved to the "Results/" folder.
#'
#' @param scores A list, NMF output from \code{compute.latent_factors()}, or a score matrix.
#'   When a list, the first element must be a samples x features score matrix.
#' @param coldata A data frame containing sample annotations, including the continuous trait.
#' @param trait Character. Name of the (numeric) column in `coldata` used as the continuous variable.
#' @param pval Numeric. P-value threshold for significance (default = 0.05).
#'
#' @return A list containing only significant features after the Kendall test.
#'         Returns \code{NULL} if no significant features are found.
#' @export
#'
scores.kendall.test <- function(scores, coldata, trait, pval = 0.05) {
  sig = c()
  coldata[, trait] = as.numeric(coldata[, trait])

  for (j in 1:ncol(scores[[1]])) {
    data = data.frame(Value = scores[[1]][, j], Trait = coldata[, trait])
    data = data[stats::complete.cases(data), ]

    res.cor <- stats::cor.test(data$Value, data$Trait, method = "kendall", exact = FALSE)

    if (round(res.cor$p.value, 5) <= pval) {
      cat("Significant p-value after Kendall correlation for", colnames(scores[[1]])[j], "\n")

      pdf(paste0("Results/Kendall_", trait, "_", colnames(scores[[1]])[j], ".pdf"),
          width = 7, height = 6)
      print(
        ggplot2::ggplot(data, ggplot2::aes(x = Trait, y = Value)) +
          ggplot2::geom_point(size = 1.8, alpha = 0.6, color = "gray20") +
          ggplot2::geom_smooth(method = "loess", formula = y ~ x, color = "#DC0000",
                               fill = "#DC0000", alpha = 0.15, linewidth = 0.9, se = TRUE) +
          ggpubr::stat_cor(method = "kendall", label.x.npc = "left", label.y.npc = "top", size = 4.5) +
          ggplot2::labs(
            title    = colnames(scores[[1]])[j],
            subtitle = paste0("Kendall correlation | ", trait),
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
      dev.off()

      sig = c(sig, j)
    }
  }

  if (length(sig) == 0) {
    message("No significant features (p-value < ", pval, ") after Kendall correlation.")
    return(NULL)
  } else {
    sig.scores = list()
    sig.scores[[1]] = scores[[1]][, sig, drop = FALSE]
    if (length(scores) >= 2) sig.scores[[2]] = scores[[2]][sig]
    if (length(scores) >= 3) sig.scores[[3]] = scores[[3]][sig]
    return(sig.scores)
  }
}

#' Linear regression test for continuous trait association
#'
#' Adapted from CellTFusion's \code{scores.anova.test()} for a continuous (rather than
#' multi-level categorical) trait. Fits a simple linear model (\code{Value ~ Trait}) for each
#' score column and tests whether the trait significantly explains variation in the score
#' using the model's F-test. Significant features are plotted as scatter plots annotated with
#' the fitted regression equation and adjusted R-squared, saved to the "Results/" folder.
#'
#' @param scores A list, NMF output from \code{compute.latent_factors()}, or a score matrix.
#'   When a list, the first element must be a samples x features score matrix.
#' @param coldata A data frame containing sample annotations including the continuous trait.
#' @param trait Character. Name of the (numeric) column in `coldata` used as the continuous variable.
#' @param pval Numeric. P-value threshold for significance (default = 0.05).
#'
#' @return A list of significant features or \code{NULL} if none are significant.
#' @export
#'
scores.lm.test = function(scores, coldata, trait, pval = 0.05){
  sig = c()
  coldata[, trait] = as.numeric(coldata[, trait])

  for (j in 1:ncol(scores[[1]])) {
    data = data.frame(Value = scores[[1]][, j], Trait = coldata[, trait])
    data = data[stats::complete.cases(data), ]

    fit <- stats::lm(Value ~ Trait, data = data)
    res.lm <- stats::anova(fit)
    lm.pval <- res.lm[["Pr(>F)"]][1]

    if (round(lm.pval, 5) <= pval) {
      cat("Significant p-value after linear regression for", colnames(scores[[1]])[j], "\n")

      pdf(paste0("Results/LM_", trait, "_", colnames(scores[[1]])[j], ".pdf"),
          width = 7, height = 6)
      print(
        ggplot2::ggplot(data, ggplot2::aes(x = Trait, y = Value)) +
          ggplot2::geom_point(size = 1.8, alpha = 0.6, color = "gray20") +
          ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "#8491B4",
                               fill = "#8491B4", alpha = 0.15, linewidth = 0.9) +
          ggpubr::stat_regline_equation(
            ggplot2::aes(label = paste(ggplot2::after_stat(eq.label), ggplot2::after_stat(adj.rr.label), sep = "~~~")),
            label.x.npc = "left", label.y.npc = "top", size = 4.2
          ) +
          ggplot2::labs(
            title    = colnames(scores[[1]])[j],
            subtitle = paste0("Linear regression (F-test) | ", trait),
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
      dev.off()

      sig = c(sig, j)
    }
  }

  sig.scores = list()
  if (length(sig) == 0) {
    message("No significant features (p-value < ", pval, ") after linear regression test")
    return(NULL)
  } else {
    sig.scores[[1]] = scores[[1]][, sig, drop = FALSE]
    if (length(scores) >= 2) sig.scores[[2]] = scores[[2]][sig]
    if (length(scores) >= 3) sig.scores[[3]] = scores[[3]][sig]
    return(sig.scores)
  }
}

#' Fisher's exact test for two continuous variables (median-split)
#'
#' Adapted from CellTFusion's \code{scores.fisher.test()}, extended to a continuous (rather
#' than categorical) trait. Both the score and the trait are binarised at their median into
#' High/Low groups, and Fisher's exact test is applied to the resulting 2x2 contingency table.
#'
#' @param scores A list, NMF output from \code{compute.latent_factors()}, or a score matrix.
#'   When a list, the first element must be a samples x features score matrix.
#'   Continuous scores are binarised at the median into High/Low groups.
#' @param coldata A data frame containing the clinical or experimental traits.
#' @param trait Character. Name of the (numeric) column in `coldata` to test with Fisher's
#'   exact test. It is binarised at the median into High/Low groups.
#' @param pval Numeric. P-value threshold for significance (default 0.05).
#'
#' @return A list containing the significant features after Fisher test. Additionally,
#'         it saves corresponding barplot visualizations in the "Results/" folder.
#' @export
#'
scores.fisher.test.continuous = function(scores, coldata, trait, pval = 0.05){

  sig = c()
  coldata[, trait] = as.numeric(coldata[, trait])
  trait_level = ifelse(coldata[, trait] > stats::median(coldata[, trait], na.rm = TRUE), "High", "Low")

  for (j in 1:ncol(scores[[1]])) {
    level = scores[[1]][, j]
    score_level = ifelse(level > stats::median(level, na.rm = TRUE), "High", "Low")

    df = data.frame(score_level = score_level, trait_level = trait_level)
    df = df[stats::complete.cases(df), ]

    contingency = table(df[, "score_level"], df[, "trait_level"])
    test = stats::fisher.test(contingency)

    if (round(test$p.value, 5) <= pval) {
      cat("Significant pval after doing Fisher test for", colnames(scores[[1]])[j], "\n")
      p_label <- ifelse(test$p.value < 0.001, "< 0.001", paste0("= ", round(test$p.value, 4)))
      pdf(paste0("Results/Fisher_", trait, "_", colnames(scores[[1]])[j], ".pdf"), width = 8, height = 6)
      print(
        ggstatsplot::ggbarstats(df, score_level, trait_level, results.subtitle = FALSE,
                                title        = colnames(scores[[1]])[j],
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
      dev.off()

      sig = c(j, sig)
    }
  }

  sig.scores = list()
  if (length(sig) == 0) {
    message("No significant features (pvalue < ", pval, ") after Fisher test")
    return(NULL)
  } else {
    sig.scores[[1]] = scores[[1]][, sig, drop = FALSE]
    if (length(scores) >= 2) sig.scores[[2]] = scores[[2]][sig]
    if (length(scores) >= 3) sig.scores[[3]] = scores[[3]][sig]
    return(sig.scores)
  }

}

#' Perform statistical analysis on scores against a continuous trait
#'
#' Adapted from CellTFusion's \code{scores.stat.analysis()} to work with any continuous
#' (numeric) clinical/experimental trait instead of a categorical grouping variable. Unified
#' interface for testing associations between any score matrix and a continuous trait.
#' Accepts cell group scores, NMF latent factors (output of \code{compute.latent_factors()}),
#' or any samples x features matrix.
#'
#' @param scores A list, NMF output from \code{compute.latent_factors()}, or a
#'   score matrix (samples x features). When a list, the first element must be
#'   the score matrix; optional second and third elements are passed through to
#'   the returned result unchanged.
#' @param coldata A data frame containing clinical or experimental metadata for samples.
#'   Must include the column specified in `trait`.
#' @param trait Character. The name of the (numeric) column in `coldata` representing the
#'   continuous clinical or experimental trait to test against (e.g., age, tumor purity,
#'   a continuous score, etc.).
#' @param method Character. Statistical test to perform. One of:
#'   \itemize{
#'     \item `"pearson"` - Pearson correlation (parametric, linear relationship)
#'     \item `"spearman"` - Spearman rank correlation (non-parametric, monotonic relationship)
#'     \item `"kendall"` - Kendall's tau (non-parametric, robust to outliers/ties, small samples)
#'     \item `"lm"` - Linear regression F-test (parametric, model-based)
#'     \item `"fisher"` - Fisher's exact test (both score and trait binarised at median)
#'   }
#'   Defaults to all available options, but only one can be used per call.
#' @param pval Numeric. P-value threshold for significance (default: 0.05).
#'
#' @details
#' The function automatically calls the corresponding statistical test function
#' based on the `method` argument:
#' \itemize{
#'   \item \code{scores.pearson.test()}
#'   \item \code{scores.spearman.test()}
#'   \item \code{scores.kendall.test()}
#'   \item \code{scores.lm.test()}
#'   \item \code{scores.fisher.test.continuous()}
#' }
#'
#' Each test produces both a statistical result and visual outputs (PDF plots)
#' stored in the `"Results/"` folder. These visualizations include the relevant
#' test results (p-values) annotated on the plots.
#'
#' @return A list of significant features, where the first element contains the
#'   subset of the original score matrix for significant features. Optional
#'   second and third elements (from the input list) are subsetted accordingly.
#'   Returns `NULL` if no significant features are found.
#'
#' @examples
#' \dontrun{
#' # Cell group scores against a continuous trait
#' sig <- scores.stat.analysis.continuous(cell.groups, coldata, trait = "tumor_purity",
#'                                        method = "spearman", pval = 0.05)
#'
#' # NMF latent factors
#' nmf <- compute.latent.factors(cell.groups)
#' sig <- scores.stat.analysis.continuous(nmf, coldata, trait = "age", method = "lm")
#' }
#'
#' @export
#'
scores.stat.analysis.continuous <- function(scores, coldata, trait,
                                            method = c("pearson", "spearman", "kendall", "lm", "fisher"),
                                            pval = 0.05) {
  method <- match.arg(method)

  # Accept NMF output from compute.latent_factors() (has $Z) or a raw score matrix
  if (!is.null(scores$Z)) {
    scores <- list(scores$Z)
  } else if (is.matrix(scores) || is.data.frame(scores)) {
    scores <- list(scores)
  }

  if (!is.numeric(coldata[, trait])) {
    stop("'trait' must be a continuous (numeric) column in 'coldata'.")
  }

  message("Running ", toupper(method), " test for score comparison...\n")

  result <- switch(method,
                   pearson  = scores.pearson.test(scores, coldata, trait, pval),
                   spearman = scores.spearman.test(scores, coldata, trait, pval),
                   kendall  = scores.kendall.test(scores, coldata, trait, pval),
                   lm       = scores.lm.test(scores, coldata, trait, pval),
                   fisher   = scores.fisher.test.continuous(scores, coldata, trait, pval))

  if (is.null(result)) {
    message("None feature was found significant after ", method, " test (p < ", pval, ").")
  } else {
    message("Significant features found: ", ncol(result[[1]]))
  }

  return(result)
}

#' Normalize raw RNA-seq counts to log-TPM and z-score scale
#'
#' Adapted from the normalization step inline in CellTFusion's \code{CellTFusion()}
#' (\code{ADImpute::NormalizeTPM(raw.counts, log = TRUE)}), extended with an optional
#' per-gene z-score scaling step. Converts a raw count matrix to TPM, log-transforms it,
#' and then centers/scales each gene (row) to mean 0 and unit variance across samples.
#'
#' @param raw.counts A raw gene expression count matrix, with genes as rows and samples as columns.
#' @param log Logical. If TRUE (default), log-transform the TPM values (as in \code{ADImpute::NormalizeTPM}).
#' @param zscore Logical. If TRUE (default), z-score scale each gene (row) across samples after
#'   TPM/log normalization.
#' @param tr_length Optional. Transcript length information passed through to
#'   \code{ADImpute::NormalizeTPM()} (see its documentation). Default \code{NULL}.
#'
#' @return A data frame of normalized expression values (genes x samples), log-TPM and/or
#'   z-score scaled according to \code{log}/\code{zscore}.
#' @export
#'
compute.normalization <- function(raw.counts, log = TRUE, zscore = TRUE, tr_length = NULL) {

  norm.counts <- data.frame(ADImpute::NormalizeTPM(raw.counts, tr_length = tr_length, log = log),
                            check.names = FALSE)

  if (zscore) {
    # Center/scale each gene (row) across samples, preserving dimnames
    norm.counts <- as.data.frame(t(scale(t(norm.counts))), check.names = FALSE)
  }

  return(norm.counts)
}

#' PCA analysis and visualization for a samples x features matrix
#'
#' Adapted from the deconvolution PCA reactive in CellTFusion's Shiny app
#' (\code{inst/shiny/server.R}, \code{deconv_pca}/\code{output$deconv_pca_plot}). Runs a
#' principal component analysis on any samples x features matrix (e.g. cell-type
#' deconvolution fractions, TF activity scores, or latent factors), optionally colored by
#' a grouping column from \code{coldata}, and saves a PC1 vs. PC2 scatter plot to the
#' "Results/" folder.
#'
#' @param data A samples x features numeric matrix or data frame (e.g. deconvolution
#'   fractions, TF activity scores, or latent factor scores). Zero-variance and
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
    pdf(paste0("Results/PCA_plot", if (!is.null(file.name)) paste0("_", file.name), ".pdf"),
        width = 7, height = 6)
    print(p)
    dev.off()
  }

  return(list(pca = pca, coords = pcs, var_explained = var_exp, plot = p))
}

#' Pearson and Spearman correlation analysis
#'
#' General-purpose correlation utility that reports both Pearson (parametric, linear) and
#' Spearman (non-parametric, monotonic) correlation together, instead of picking one method.
#' Reuses \code{stats::cor.test()} the same way \code{scores.pearson.test()}/
#' \code{scores.spearman.test()} already do for the vector and scores-vs-trait cases, and
#' \code{WGCNA::cor()}/\code{WGCNA::corPvalueStudent()} the same way CellTFusion's
#' \code{compute.modules.relationship()} does for the matrix-vs-matrix case. Handles three
#' input shapes:
#' \enumerate{
#'   \item Two numeric vectors \code{x} and \code{y} - returns a small data frame with both
#'     correlation coefficients and p-values.
#'   \item A score matrix \code{x} (or NMF output with \code{$Z}) plus \code{coldata}/\code{trait}
#'     - calls the existing \code{scores.pearson.test()} and \code{scores.spearman.test()} on
#'     every column, unions the features significant by either method, and returns their
#'     subset (same shape as those functions) plus a combined per-feature statistics table.
#'   \item Two matrices/data frames \code{x} and \code{y} (samples x features each, same
#'     samples in the same order) - correlates every column of \code{x} against every column
#'     of \code{y} with both methods and, if \code{file_name} is given, saves a heatmap per
#'     method to \code{Results/}.
#' }
#'
#' @param x A numeric vector, a score matrix/NMF output (samples x features), or a
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
#'   \item Case 2: a list with the significant score subset (as returned by
#'     \code{scores.pearson.test()}) plus a \code{$stats} data frame of per-feature
#'     Pearson/Spearman estimates and p-values. Returns \code{NULL} (with a message) if no
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

    if (is.list(x) && !is.null(x$Z)) {
      scores_list <- list(x$Z)
    } else if (is.matrix(x) || is.data.frame(x)) {
      scores_list <- list(x)
    } else {
      scores_list <- x
    }

    coldata[, trait] <- as.numeric(coldata[, trait])

    message("Running Pearson and Spearman correlation for score comparison...\n")

    sig.pearson  <- scores.pearson.test(scores_list, coldata, trait, pval)
    sig.spearman <- scores.spearman.test(scores_list, coldata, trait, pval)

    sig_names <- union(
      if (!is.null(sig.pearson))  colnames(sig.pearson[[1]])  else character(0),
      if (!is.null(sig.spearman)) colnames(sig.spearman[[1]]) else character(0)
    )

    # Combined per-feature correlation statistics (both methods, every feature)
    stats_tab <- do.call(rbind, lapply(colnames(scores_list[[1]]), function(feat) {
      data <- data.frame(Value = scores_list[[1]][, feat], Trait = coldata[, trait])
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

    idx <- match(sig_names, colnames(scores_list[[1]]))
    sig.scores <- list()
    sig.scores[[1]] <- scores_list[[1]][, sig_names, drop = FALSE]
    if (length(scores_list) >= 2) sig.scores[[2]] <- scores_list[[2]][idx]
    if (length(scores_list) >= 3) sig.scores[[3]] <- scores_list[[3]][idx]
    sig.scores$stats <- stats_tab[stats_tab$feature %in% sig_names, ]

    return(sig.scores)
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

        pdf(paste0("Results/Correlation_", m, "_", file_name, ".pdf"), width = 8, height = 8)
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
        dev.off()
      }
    }

    return(list(
      pearson  = list(cor = cor.pearson,  pval = pval.pearson),
      spearman = list(cor = cor.spearman, pval = pval.spearman)
    ))
  }

  stop("Provide either two numeric vectors ('x','y'), a score matrix 'x' with 'coldata'/'trait', or two matrices 'x' and 'y'.")
}

#' Remove batch effects from raw RNA-seq counts (ComBat-seq)
#'
#' Adapted from the batch-correction step used in CellTFusion's Shiny app
#' (\code{CellTFusion_paper/man/gui.R}, the "Correct Batch Effect" button /
#' \code{output$OutputBatch}) and in the multi-cohort Vanderbilt analysis
#' (\code{LungPredict1_paper/Scripts/Vanderbilt_analysis.Rmd}). Both correct raw counts
#' across batches/cohorts with \code{sva::ComBat_seq()} - a negative-binomial model
#' appropriate for raw RNA-seq counts - before TPM normalization:
#' \preformatted{
#' GEM_corrected <- sva::ComBat_seq(as.matrix(raw_counts), batch = as.vector(traitData[, BatchCol]))
#' counts.normalized <- ADImpute::NormalizeTPM(GEM_corrected, log = TRUE)
#' }
#' This function wraps that same two-step pattern, optionally chaining into
#' \code{compute.normalization()} (log-TPM + z-score) for the corrected counts.
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
#'   with \code{compute.normalization()} (log-TPM + z-score) before being returned, matching
#'   the ComBat_seq -> NormalizeTPM pattern used in the GUI/Vanderbilt scripts.
#' @param log Passed to \code{compute.normalization()} when \code{normalize = TRUE}.
#' @param zscore Passed to \code{compute.normalization()} when \code{normalize = TRUE}.
#'
#' @return If \code{normalize = FALSE}, a data frame of batch-corrected raw counts
#'   (genes x samples). If \code{normalize = TRUE} (default), a data frame of
#'   batch-corrected, log-TPM (and optionally z-score) normalized expression values.
#' @export
#'
compute.batch.correction <- function(raw.counts, batch, group = NULL, covar_mod = NULL,
                                     normalize = TRUE, log = TRUE, zscore = TRUE) {

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
