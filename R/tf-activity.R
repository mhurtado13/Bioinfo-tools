#' Compute Transcription Factor (TF) activity
#'
#' Infers transcription factor (TF) activity from a gene expression matrix using the VIPER algorithm (Alvarez et al., 2016). The function requires a TF-target gene regulatory network, which can be provided by the user or obtained from OmnipathR resources such as CollecTRI or Dorothea. ARACNE-inferred networks are also supported.
#'
#' @param RNA.counts A gene expression matrix with genes as rows and samples as columns. The matrix should be
#'   log-normalized -- log2(TPM + pseudocount), log2CPM, etc. -- not raw or linear-scale TPM/counts (e.g. the
#'   output of \code{compute.normalization(..., log = TRUE, zscore = FALSE)}).
#' @param TF.collection Character. The source of the TF-target network. Options are `"CollecTRI"` (default), `"Dorothea"`, or `"ARACNE"`.
#' - `"CollecTRI"` and `"Dorothea"` use prebuilt collections from OmnipathR.
#' - `"ARACNE"` allows user input of a custom network file in a 3-column format: `regulator`, `target`, and `mutual information`.
#' @param min_targets_size Integer. Minimum number of target genes per regulon required for TF activity inference. Default is 5.
#' @param universe Optional. A user-specified data frame of TF-target interactions. If not provided, the function will fetch the relevant network based on the `TF.collection` argument.
#' @param confidence Character vector. Dorothea confidence levels to keep (ignored for other `TF.collection` values), from `"A"` (highest) to `"E"` (lowest). Default `c("A", "B")`.
#' @param aracne.network Character. Path to the ARACNE `network.txt` file (3 columns: `Regulator`, `Target`, mutual information), required when `TF.collection = "ARACNE"`.
#' @param statistic Character. Which \code{decoupleR::decouple()} statistic to extract. Default `"consensus"`
#'   (an ensemble score across the top-performing methods: MLM, ULM, and normalized WSUM). Any single
#'   underlying method can be requested instead: `"aucell"`, `"udt"`, `"mdt"`, `"wmean"`, `"ulm"`, `"mlm"`,
#'   `"wsum"`, `"viper"`, `"gsva"`, `"ora"`, or `"fgsea"` (see \code{decoupleR::show_methods()}).
#' @param cores Integer. Reserved for parallelizing VIPER inference; not currently used by this
#'   function. Default is 3.
#' @param return Logical; if TRUE, saves matrix in Results/ folder. Default is TRUE.
#' @param file.name Optional character suffix used when writing the TF activity matrix to disk.
#'
#' @return A data frame of inferred and scaled TF activity scores, with samples as rows and TFs as columns.
#'
#' @references
#'
#' Alvarez, M. et al. (2016). Functional characterization of somatic mutations in cancer using network-based inference of protein activity. *Nature Genetics*, 48(8), 838-847. https://doi.org/10.1038/ng.3593
#'
#' Tuerei, D., Korcsmaros, T., & Saez-Rodriguez, J. (2016). OmniPath: guidelines and gateway for literature-curated signaling pathway resources. *Nature Methods*, 13(12), 966-967. https://doi.org/10.1038/nmeth.4077
#'
#' Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. Benchmark and integration of resources for the estimation of human transcription factor activities. Genome Research. 2019. DOI: 10.1101/gr.240663.118.
#'
#' Lachmann, A. et al. (2016). ARACNe-AP: gene network reverse engineering through adaptive partitioning inference of mutual information. *Bioinformatics*, 32(14), 2233-2235. https://doi.org/10.1093/bioinformatics/btw216
#'
#' Margolin, A.A. et al. (2006). ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. *BMC Bioinformatics*, 7(Suppl 1), S7. https://doi.org/10.1186/1471-2105-7-S1-S7
#'
#' Badia-i-Mompel P., Vélez Santiago J., Braunger J., Geiss C., Dimitrov D., Müller-Dott S., Taus P., Dugourd A., Holland C.H., Ramirez Flores R.O. and Saez-Rodriguez J. 2022. decoupleR: ensemble of computational methods to infer biological activities from omics data. Bioinformatics Advances. https://doi.org/10.1093/bioadv/vbac016
#'
#' S. Hänzelmann, R. Castelo, and J. Guinney. “GSVA: gene set variation analysis for microarray and RNA-Seq data”. In: BMC Bioinformatics 14 (2013), p. 7. DOI: 10.1186/1471-2105-14-7. URL: https://doi.org/10.1186/1471-2105-14-7.
#'
#' G. Korotkevich, V. Sukhov, and A. Sergushichev. “Fast gene set enrichment analysis”. In: bioRxiv (2019). DOI: 10.1101/060012. URL: http://biorxiv.org/content/early/2016/06/20/060012.
#'
#' Müller-Dott S., Tsirvouli E., Vázquez M., Ramirez Flores R.O., Badia-i-Mompel P., Fallegger R., Lægreid A. and Saez-Rodriguez J. Expanding the coverage of regulons from high-confidence prior knowledge for accurate estimation of transcription factor activities. bioRxiv. 2023. DOI: 10.1101/2023.03.30.534849
#'
#' @examples
#' \dontrun{
#' tfs_activity <- compute.TFs.activity(counts.norm, file.name = "Tutorial")
#' }
#'
#' @export
#'
compute.TFs.activity <- function(RNA.counts, TF.collection = "CollecTRI", min_targets_size = 5, universe = NULL,
                                 confidence = c("A", "B"), aracne.network = NULL,
                                 statistic = c("consensus", "aucell", "udt", "mdt", "wmean", "ulm", "mlm", "wsum", "viper", "gsva", "ora", "fgsea"),
                                 cores = 3, return = TRUE, file.name = NULL){

  statistic <- match.arg(statistic)
  confidence <- match.arg(confidence, choices = c("A", "B", "C", "D", "E"), several.ok = TRUE)
  tf_cache_file <- "Results/TF_target_collection.csv"

  if(TF.collection == "ARACNE"){

    if(is.null(aracne.network) || !file.exists(aracne.network))
      stop("TF.collection = 'ARACNE' requires 'aracne.network' to point to a valid network.txt file.")

    cat("Loading ARACNe network from:", aracne.network, "\n")
    # Read network edges, filter to genes present in expression matrix
    aracne_net <- utils::read.table(aracne.network, header = TRUE, sep = "\t") %>%
      dplyr::select(source = Regulator, target = Target) %>%
      dplyr::filter(source %in% rownames(RNA.counts) & target %in% rownames(RNA.counts))

    # Compute Spearman correlation between every TF and its targets in one matrix 
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
          universe = dplyr::filter(dorothea::dorothea_hs, .data$confidence %in% .env$confidence) %>%
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
                                      minsize = min_targets_size,
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

#' TF activity from a DEG signature, with permutation-based significance
#'
#' Computes TF activity directly from a differential expression signature -- a single, per-gene
#' ranking statistic, rather than per-sample expression: the signature is treated as a
#' single-column "sample" and scored directly with \code{compute.TFs.activity()}. In addition to
#' that score (\code{"consensus"} by default, or whichever \code{statistic} is chosen), an
#' empirical p-value is computed by repeatedly shuffling which gene gets which signature value
#' (breaking the gene-TF regulon relationship while preserving the signature's value
#' distribution), rescoring every TF against each shuffled signature, and comparing the observed
#' score to that null distribution -- a permutation-based null model, generalized to any
#' \code{decoupleR} statistic rather than a single bootstrapped model. Also produces a top-TF bar
#' plot.
#'
#' @param signature A named numeric vector: gene symbols as names, one differential expression
#'   statistic per gene as values (unfiltered, covering every tested gene, not only the
#'   significant ones). We recommend the \code{t} statistic from \code{run_deg_analysis()} --
#'   e.g. \code{setNames(deg$t, rownames(deg))} on \code{run_deg_analysis(..., pval = 1)}'s
#'   output -- since it combines effect size and precision, but the choice of statistic (e.g.
#'   log2 fold-change) is entirely up to the caller.
#' @param TF.collection Character. The source of the TF-target network, exactly as in
#'   \code{compute.TFs.activity()}: \code{"CollecTRI"} (default), \code{"Dorothea"}, or
#'   \code{"ARACNE"}.
#' @param min_targets_size Integer. Minimum number of target genes per regulon required for TF
#'   activity inference (passed to \code{compute.TFs.activity()} and to every permutation).
#'   Default is 5.
#' @param universe Optional. A user-specified data frame of TF-target interactions (\code{source},
#'   \code{target}, \code{mor}), reused for the observed score and every permutation, exactly as
#'   in \code{compute.TFs.activity()}. If not provided, the function fetches/caches the network
#'   based on \code{TF.collection} instead. Ignored when \code{TF.collection = "ARACNE"} (that
#'   network is always rebuilt from \code{RNA.counts}, as in \code{compute.TFs.activity()}).
#' @param confidence Character vector. Dorothea confidence levels to keep (ignored for other
#'   \code{TF.collection} values). Default \code{c("A", "B")}.
#' @param aracne.network Character. Path to the ARACNE \code{network.txt} file, required when
#'   \code{TF.collection = "ARACNE"} (see \code{compute.TFs.activity()}).
#' @param RNA.counts Optional. The full genes x samples normalized expression matrix, required
#'   only when \code{TF.collection = "ARACNE"}: unlike CollecTRI/Dorothea, ARACNE's mode of
#'   regulation is not pre-signed, but derived by \code{compute.TFs.activity()} from cross-sample
#'   Spearman correlation between each TF and its targets -- which needs real multi-sample
#'   expression, not the single-column DEG signature being scored here. Ignored for other
#'   \code{TF.collection} values.
#' @param statistic Character. Which \code{decoupleR::decouple()} statistic to score the
#'   signature with. Default \code{"consensus"} (an ensemble score across the top-performing
#'   methods: MLM, ULM, and normalized WSUM), same options as \code{compute.TFs.activity()}:
#'   \code{"aucell"}, \code{"udt"}, \code{"mdt"}, \code{"wmean"}, \code{"ulm"}, \code{"mlm"},
#'   \code{"wsum"}, \code{"viper"}, \code{"gsva"}, \code{"ora"}, or \code{"fgsea"}. Note that
#'   \code{"consensus"}, \code{"udt"}, \code{"mdt"}, \code{"ora"}, \code{"fgsea"}, and
#'   \code{"gsva"} are considerably slower per permutation than \code{"viper"}/\code{"ulm"}/
#'   \code{"mlm"}/\code{"wsum"}/\code{"wmean"}/\code{"aucell"}.
#' @param n_perm Integer. Number of gene-label permutations used to build the null distribution
#'   for the empirical p-value. Permutations are internally scored in memory-bounded batches
#'   (far faster than one \code{decoupleR::decouple()} call per permutation, since the network
#'   preprocessing cost is paid once per batch, not once per permutation). Default 1000.
#' @param seed Integer. Random seed for the permutations. Default 42.
#' @param top_n Integer. Number of top up- and down-activity TFs shown in the plot. Default 10.
#' @param return Logical; if TRUE (default), saves the activity table (CSV) and plot (PDF) to
#'   \code{Results/}.
#' @param file.name Optional character suffix used when saving outputs to \code{Results/}.
#'
#' @return A list with:
#' \describe{
#'   \item{activity}{A data frame with one row per TF: \code{TF}, \code{score} (the
#'     \code{statistic} score from the observed signature), and \code{p_value} (the empirical,
#'     permutation-based two-sided p-value), sorted by \code{score} descending.}
#'   \item{plot}{A \code{ggplot2} bar plot of the top \code{top_n} up- and down-activity TFs.}
#' }
#'
#' @examples
#' \dontrun{
#' deg_full <- run_deg_analysis(counts, coldata, group_col = "Group", pval = 1, plot = FALSE)
#' signature <- setNames(deg_full$t, rownames(deg_full))
#' tf_deg <- compute.TFs.activity.DEG(signature, n_perm = 1000, file.name = "Group")
#' head(tf_deg$activity)
#' tf_deg$plot
#'
#' # ARACNE: mode of regulation derived from the full expression matrix
#' tf_deg_aracne <- compute.TFs.activity.DEG(signature, TF.collection = "ARACNE",
#'                                           aracne.network = "input/ARACNE/luad/network/network.txt",
#'                                           RNA.counts = counts.norm, n_perm = 1000)
#' }
#'
#' @export
#'
compute.TFs.activity.DEG <- function(signature, TF.collection = "CollecTRI", min_targets_size = 5, universe = NULL,
                                     confidence = c("A", "B"), aracne.network = NULL, RNA.counts = NULL,
                                     statistic = c("consensus", "aucell", "udt", "mdt", "wmean", "ulm", "mlm",
                                                  "wsum", "viper", "gsva", "ora", "fgsea"),
                                     n_perm = 1000, seed = 42,
                                     top_n = 10, return = TRUE, file.name = NULL) {

  statistic <- match.arg(statistic)
  confidence <- match.arg(confidence, choices = c("A", "B", "C", "D", "E"), several.ok = TRUE)

  # Permutations are scored in memory-bounded batches rather than all at once or one at a time
  # (not user-configurable: 100 is a safe, validated batch size)
  chunk_size <- 100

  if (TF.collection == "ARACNE" && is.null(RNA.counts)) {
    stop("TF.collection = 'ARACNE' requires 'RNA.counts' (the full genes x samples expression matrix) ",
        "to derive each edge's mode of regulation via cross-sample correlation.")
  }

  if (!is.numeric(signature) || is.null(names(signature))) {
    stop("'signature' must be a named numeric vector (names = gene symbols), e.g. ",
        "setNames(deg$t, rownames(deg)).")
  }

  sig_vec <- signature[!is.na(signature)]
  sig_mat <- matrix(sig_vec, ncol = 1, dimnames = list(names(sig_vec), "signature"))

  ###### Build the TF-target network once, reused for both the observed score and every permutation.
  # Built here (rather than left to compute.TFs.activity()) because its ARACNE branch derives mor from
  # cross-sample correlation on whatever matrix it is asked to score -- which would be this signature's
  # single column, not real multi-sample expression, if left to compute internally.
  tf_cache_file <- "Results/TF_target_collection.csv"
  if (TF.collection == "ARACNE") {
    aracne_net <- utils::read.table(aracne.network, header = TRUE, sep = "\t") %>%
      dplyr::select(source = Regulator, target = Target) %>%
      dplyr::filter(source %in% rownames(RNA.counts) & target %in% rownames(RNA.counts))

    all_tfs <- unique(aracne_net$source)
    all_targets <- unique(aracne_net$target)
    cor_mat <- suppressWarnings(
      stats::cor(t(RNA.counts[all_tfs, , drop = FALSE]),
                t(RNA.counts[all_targets, , drop = FALSE]),
                method = "spearman")
    )

    perm_universe <- aracne_net %>%
      dplyr::mutate(mor = cor_mat[cbind(source, target)]) %>%
      dplyr::mutate(mor = sign(mor)) %>%
      dplyr::filter(!is.na(mor) & mor != 0)

  } else if (!is.null(universe)) {
    perm_universe <- universe
  } else if (TF.collection == "CollecTRI") {
    if (file.exists(tf_cache_file)) {
      perm_universe <- utils::read.csv(tf_cache_file, row.names = 1)
    } else {
      perm_universe <- decoupleR::get_collectri(organism = 'human', split_complexes = FALSE)
      utils::write.csv(perm_universe, tf_cache_file)
    }
  } else if (TF.collection == "Dorothea") {
    perm_universe <- dplyr::filter(dorothea::dorothea_hs, .data$confidence %in% .env$confidence) %>%
      dplyr::mutate(source = .data$tf) %>%
      dplyr::select(-tf)
  } else {
    stop("Unknown TF.collection: '", TF.collection, "'. Use 'CollecTRI', 'Dorothea', or 'ARACNE'.")
  }

  ###### Observed TF activity from the real DEG signature, scored against the network built above
  observed <- compute.TFs.activity(sig_mat, TF.collection = "CollecTRI", min_targets_size = min_targets_size,
                                   universe = perm_universe, statistic = statistic, return = FALSE)
  obs_score <- stats::setNames(as.numeric(observed[1, ]), colnames(observed))

  ###### Null distribution: shuffle which gene gets which signature value, in memory-safe batches
  set.seed(seed)
  null_scores <- matrix(NA_real_, nrow = 0, ncol = length(obs_score), dimnames = list(NULL, names(obs_score)))

  n_done <- 0
  while (n_done < n_perm) {
    this_chunk <- min(chunk_size, n_perm - n_done)
    perm_mat <- replicate(this_chunk, sample(sig_vec))
    rownames(perm_mat) <- names(sig_vec)

    chunk_res <- decoupleR::decouple(mat = perm_mat, network = perm_universe, .source = "source", .target = "target",
                                     statistics = if (statistic == "consensus") NULL else statistic,
                                     consensus_score = (statistic == "consensus"), minsize = min_targets_size) %>%
      dplyr::filter(.data$statistic == .env$statistic) %>%
      tidyr::pivot_wider(id_cols = "condition", names_from = "source", values_from = "score") %>%
      tibble::column_to_rownames("condition") %>%
      as.matrix()
    colnames(chunk_res) <- make.names(colnames(chunk_res))

    pad <- matrix(NA_real_, nrow = nrow(chunk_res), ncol = length(obs_score), dimnames = list(NULL, names(obs_score)))
    common <- intersect(colnames(chunk_res), names(obs_score))
    pad[, common] <- chunk_res[, common]
    null_scores <- rbind(null_scores, pad)

    n_done <- n_done + this_chunk
  }

  ###### Empirical p-value: two-sided, fraction of |null score| >= |observed score|
  emp_pval <- sapply(names(obs_score), function(tf) {
    null_vals <- stats::na.omit(null_scores[, tf])
    if (length(null_vals) == 0) return(NA_real_)
    (sum(abs(null_vals) >= abs(obs_score[[tf]])) + 1) / (length(null_vals) + 1)
  })

  result <- data.frame(TF = names(obs_score), score = unname(obs_score), p_value = unname(emp_pval)) %>%
    dplyr::arrange(dplyr::desc(score))

  ###### Plot: top/bottom TFs by score, labeled with the empirical p-value
  # tail() excludes rows already in head() so a small network (fewer TFs than 2*top_n) never
  # lists the same TF as both "Up" and "Down"
  top_up <- utils::head(result, top_n)
  remaining <- result[!(result$TF %in% top_up$TF), ]
  top_down <- utils::tail(remaining, min(top_n, nrow(remaining)))

  plot_data <- dplyr::bind_rows(
    top_up %>% dplyr::mutate(direction = "Up"),
    top_down %>% dplyr::mutate(direction = "Down")
  )
  plot_data$TF <- factor(plot_data$TF, levels = rev(plot_data$TF))

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = score, y = TF, fill = direction)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(ggplot2::aes(label = signif(p_value, 2), hjust = ifelse(score > 0, -0.1, 1.1)), size = 3) +
    ggplot2::scale_fill_manual(values = c(Up = "firebrick", Down = "steelblue")) +
    ggplot2::labs(x = paste0(statistic, " score"), y = NULL, fill = NULL,
                 title = "Top TFs from DEG signature",
                 subtitle = paste0("Labels: empirical p-value from ", n_perm, " gene-label permutations")) +
    ggplot2::theme_bw()

  if (return) {
    utils::write.csv(result, paste0("Results/TF_DEG_activity_", file.name, ".csv"), row.names = FALSE)
    grDevices::pdf(paste0("Results/TF_DEG_activity_", file.name, ".pdf"), width = 8, height = 6)
    print(p)
    grDevices::dev.off()
  }

  return(list(activity = result, plot = p))
}
