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
#' @param return Logical; if TRUE (default), saves the plot as a PDF in Results/.
#' @param file.name Optional character suffix used when saving the plot to \code{Results/}.
#'
#' @return The \code{ggplot} object produced by \code{EnhancedVolcano::EnhancedVolcano()}.
#'
#' @examples
#' \dontrun{
#' deg_full <- run_deg_analysis(counts, coldata, group_col = "Group", plot = FALSE)
#' compute.volcano.plot(deg_full, pval = 0.05, logFC = 1, file.name = "Group")
#' }
#'
#' @export
#'
compute.volcano.plot <- function(deg, pval = 0.05, logFC = 1, title = NULL, return = TRUE, file.name = NULL) {

  p <- EnhancedVolcano::EnhancedVolcano(deg,
                                        lab       = rownames(deg),
                                        x         = "logFC",
                                        y         = "adj.P.Val",
                                        pCutoff   = pval,
                                        FCcutoff  = logFC,
                                        title     = if (is.null(title)) "Volcano plot" else title,
                                        subtitle  = NULL,
                                        ylab      = bquote(~-Log[10] ~ italic(adj.P.Val)))

  if (return) {
    grDevices::pdf(paste0("Results/Volcano_", file.name, ".pdf"), width = 10, height = 8)
    print(p)
    grDevices::dev.off()
  }

  return(p)
}

#' Gene set enrichment analysis (ORA or GSEA) across multiple databases
#'
#' Tests a gene list for over-representation (ORA) or a ranked gene list for gene set enrichment
#' (GSEA) against Gene Ontology (Biological Process, Cellular Component, or Molecular Function),
#' KEGG, Reactome, PROGENy, or MSigDB Hallmark. Wraps \code{clusterProfiler}/\code{ReactomePA},
#' following the same ORA/GSEA workflow already used elsewhere (\code{enrichKEGG()}/\code{enrichPathway()}
#' on DEG targets, \code{GSEA()} with a Hallmark \code{TERM2GENE} on a ranked gene list).
#'
#' @param genes For `method = "ORA"`, a character vector of gene SYMBOLs to test (e.g.
#'   \code{rownames(run_deg_analysis(...))}). For `method = "GSEA"`, a named numeric vector of a ranking
#'   statistic (e.g. \code{logFC} or \code{t}) with gene SYMBOLs as names, for every tested gene --
#'   not pre-filtered to significant ones (e.g. from \code{run_deg_analysis(..., pval = 1)}).
#' @param database Character. Resource to test against: `"GO_BP"` (default), `"GO_CC"`, `"GO_MF"`,
#'   `"KEGG"`, `"REACTOME"`, `"PROGENy"`, or `"Hallmark"`.
#' @param method Character. `"ORA"` (over-representation, default) or `"GSEA"` (gene set enrichment).
#' @param universe Optional character vector of background gene SYMBOLs, used only for `method = "ORA"`
#'   (e.g. all genes tested in the DEG analysis, \code{rownames(counts)}). If \code{NULL}, each tool's
#'   own default universe is used.
#' @param pval Numeric. Adjusted p-value cutoff for both ORA and GSEA. Default 0.05.
#' @param top_n Integer. Number of top terms shown in the dot plot. Default 20.
#' @param return Logical; if TRUE (default), saves a dot plot as a PDF in Results/.
#' @param file.name Optional character suffix used when saving the plot.
#'
#' @return A \code{clusterProfiler}/\code{ReactomePA} enrichment result object (\code{enrichResult} for
#'   ORA, \code{gseaResult} for GSEA), or \code{NULL} (with a message) if no term passed `pval`.
#'
#' @examples
#' \dontrun{
#' deg <- run_deg_analysis(counts, coldata, group_col = "Group")
#' ora_kegg <- compute.enrichment(rownames(deg), database = "KEGG", universe = rownames(counts))
#'
#' deg_full <- run_deg_analysis(counts, coldata, group_col = "Group", pval = 1)
#' ranks <- setNames(deg_full$t, rownames(deg_full))
#' gsea_hallmark <- compute.enrichment(ranks, database = "Hallmark", method = "GSEA")
#' }
#'
#' @references
#' Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012 May;16(5):284-7. doi: 10.1089/omi.2011.0118. Epub 2012 Mar 28. PMID: 22455463; PMCID: PMC3339379.
#'
#' Yu G, He QY. ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. Molecular BioSystems. 2016, 12(2):477-479. doi: 10.1039/C5MB00663E.
#'
#' @export
#'
compute.enrichment <- function(genes, database = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "REACTOME", "PROGENy", "Hallmark"),
                               method = c("ORA", "GSEA"), universe = NULL, pval = 0.05, top_n = 20,
                               return = TRUE, file.name = NULL) {

  database <- match.arg(database)
  method <- match.arg(method)

  # SYMBOL -> ENTREZID, for the databases keyed on Entrez (KEGG, REACTOME)
  to_entrez <- function(symbols) {
    map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = unique(symbols),
                                                  columns = "ENTREZID", keytype = "SYMBOL"))
    unique(stats::na.omit(map$ENTREZID))
  }

  # Re-key a SYMBOL-named ranked vector to ENTREZID, for gseKEGG/gsePathway
  to_entrez_ranked <- function(ranked) {
    map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = names(ranked),
                                                  columns = "ENTREZID", keytype = "SYMBOL"))
    map <- map[!is.na(map$ENTREZID) & !duplicated(map$SYMBOL), ]
    sort(stats::setNames(ranked[map$SYMBOL], map$ENTREZID), decreasing = TRUE)
  }

  # TERM2GENE for the resources without a dedicated clusterProfiler/ReactomePA wrapper,
  # reusing the same Results/ cache as compute.pathway.activity()
  progeny_term2gene <- function() {
    progeny_cache_file <- "Results/Pathways_collection_PROGENy.csv"
    if (file.exists(progeny_cache_file)) {
      net <- utils::read.csv(progeny_cache_file, row.names = 1)
    } else {
      net <- decoupleR::get_progeny(organism = "human", top = 500)
      utils::write.csv(net, progeny_cache_file)
    }
    dplyr::transmute(net, term = source, gene = target)
  }

  hallmark_term2gene <- function() {
    hallmark_cache_file <- "Results/Pathways_collection_Hallmark.csv"
    if (file.exists(hallmark_cache_file)) {
      sets <- utils::read.csv(hallmark_cache_file)
    } else {
      sets <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
      utils::write.csv(sets, hallmark_cache_file, row.names = FALSE)
    }
    dplyr::distinct(dplyr::transmute(sets, term = gs_name, gene = gene_symbol))
  }

  if (method == "ORA") {

    ###### Over-representation analysis
    gene_list <- if (is.numeric(genes)) names(genes) else genes

    result <- switch(database,
      GO_BP    = clusterProfiler::enrichGO(gene = gene_list, universe = universe, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                           keyType = "SYMBOL", ont = "BP", pvalueCutoff = pval),
      GO_CC    = clusterProfiler::enrichGO(gene = gene_list, universe = universe, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                           keyType = "SYMBOL", ont = "CC", pvalueCutoff = pval),
      GO_MF    = clusterProfiler::enrichGO(gene = gene_list, universe = universe, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                           keyType = "SYMBOL", ont = "MF", pvalueCutoff = pval),
      KEGG     = clusterProfiler::setReadable(
                   clusterProfiler::enrichKEGG(gene = to_entrez(gene_list), organism = "hsa",
                                               universe = if (!is.null(universe)) to_entrez(universe) else NULL,
                                               pvalueCutoff = pval),
                   OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID"),
      REACTOME = ReactomePA::enrichPathway(gene = to_entrez(gene_list), organism = "human",
                                           universe = if (!is.null(universe)) to_entrez(universe) else NULL,
                                           pvalueCutoff = pval, readable = TRUE),
      PROGENy  = clusterProfiler::enricher(gene = gene_list, universe = universe, pvalueCutoff = pval,
                                           TERM2GENE = progeny_term2gene()),
      Hallmark = clusterProfiler::enricher(gene = gene_list, universe = universe, pvalueCutoff = pval,
                                           TERM2GENE = hallmark_term2gene())
    )

  } else {

    ###### Gene set enrichment analysis
    if (!is.numeric(genes) || is.null(names(genes)))
      stop("method = 'GSEA' requires 'genes' to be a named numeric ranking vector (names = gene symbols).")

    gene_list <- sort(genes, decreasing = TRUE)

    result <- switch(database,
      GO_BP    = clusterProfiler::gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                        keyType = "SYMBOL", ont = "BP", pvalueCutoff = pval),
      GO_CC    = clusterProfiler::gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                        keyType = "SYMBOL", ont = "CC", pvalueCutoff = pval),
      GO_MF    = clusterProfiler::gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                        keyType = "SYMBOL", ont = "MF", pvalueCutoff = pval),
      KEGG     = clusterProfiler::setReadable(
                   clusterProfiler::gseKEGG(geneList = to_entrez_ranked(gene_list), organism = "hsa", pvalueCutoff = pval),
                   OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID"),
      REACTOME = ReactomePA::gsePathway(geneList = to_entrez_ranked(gene_list), organism = "human", pvalueCutoff = pval),
      PROGENy  = clusterProfiler::GSEA(geneList = gene_list, TERM2GENE = progeny_term2gene(), pvalueCutoff = pval),
      Hallmark = clusterProfiler::GSEA(geneList = gene_list, TERM2GENE = hallmark_term2gene(), pvalueCutoff = pval)
    )
  }

  if (is.null(result) || nrow(result@result) == 0) {
    message("No significant terms found for database = '", database, "', method = '", method, "' (p < ", pval, ").")
    return(invisible(NULL))
  }

  if (return) {
    grDevices::pdf(paste0("Results/Enrichment_", method, "_", database, "_", file.name, ".pdf"), width = 8, height = 10)
    print(enrichplot::dotplot(result, showCategory = top_n, title = paste(database, method)))
    grDevices::dev.off()
  }

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

#' Kaplan-Meier survival analysis on clinical groups
#'
#' Performs Kaplan-Meier survival analysis and log-rank testing, either on a predefined
#' grouping variable (e.g. a clinical risk group or a supervised cell group split), or
#' automatically for each column of a feature matrix. In the latter case, samples
#' are split into "High"/"Low" groups per feature using a quantile cutoff, and only
#' features with a significant log-rank test are returned.
#'
#' @param survival.data A data frame of clinical/survival metadata, with samples as rows.
#' @param PFS Character. Column name in \code{survival.data} with the survival/follow-up time.
#' @param PFS_event Character. Column name in \code{survival.data} with the event indicator
#'   (1 = event occurred, 0 = censored).
#' @param file_name Optional character. Suffix used when saving Kaplan-Meier plots to \code{Results/}.
#' @param features Optional. A samples x features numeric matrix or data frame. If provided 
#' (and \code{group_column} is \code{NULL}), each feature is tested individually. 
#'  Mutually exclusive with \code{group_column}.
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
#' Requires the \code{survival}, \code{survminer}, and \code{gridExtra} packages to be installed.
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
#' # Automatic screening of a feature matrix (e.g. pathway activity scores)
#' compute.survival.analysis(
#'   survival.data = traitdata,
#'   PFS           = "PFS",
#'   PFS_event     = "PFS_event",
#'   features      = pathway_scores,
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

      pdf(paste0("Results/Pearson_", trait, "_", colnames(scores)[j], ".pdf"),
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
      dev.off()

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

      pdf(paste0("Results/Spearman_", trait, "_", colnames(scores)[j], ".pdf"),
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
      dev.off()

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
      pdf(paste0("Results/Fisher_", trait, "_", colnames(scores)[j], ".pdf"), width = 8, height = 6)
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
      dev.off()

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

      pdf(paste0("Results/Ttest_", trait, "_", colnames(scores)[j], ".pdf"),
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
      dev.off()

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

      pdf(paste0("Results/Wilcoxon_", trait, "_", colnames(scores)[j], ".pdf"),
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
      dev.off()

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

      pdf(paste0("Results/ANOVA_", trait, "_", colnames(scores)[j], ".pdf"),
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
      dev.off()

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

      pdf(paste0("Results/Kruskal_", trait, "_", colnames(scores)[j], ".pdf"),
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
      dev.off()

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
