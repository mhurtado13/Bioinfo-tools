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
