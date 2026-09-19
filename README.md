# Bioinfo-tools

Collection of general-purpose R functions for bulk RNA-seq analysis, gathered in a single
[`functions.r`](functions.r) file and shared across projects. Covers transcription factor
and pathway activity inference, differential expression, survival analysis, and
score-vs-trait association testing.

## Functions

| Function | Description |
| --- | --- |
| `compute.TFs.activity()` | Infers transcription factor activity from a normalized expression matrix, with a TF-target network from CollecTRI, Dorothea, or a user-supplied ARACNE network, scored by a `decoupleR` statistic (`"consensus"` by default, or a single method such as `"mlm"`/`"viper"`). |
| `compute.pathway.activity()` | Computes pathway activity scores, from PROGENy (MLM, default) or KEGG/REACTOME/Hallmark (GSVA, via `pathway_source`), optionally combined with GSVA scores on custom gene sets. |
| `run_deg_analysis()` | Differential expression analysis between two groups with edgeR/limma-voom (filtering, TMM normalization, voom, `topTable`). |
| `compute.survival.analysis()` | Kaplan-Meier analysis and log-rank testing, either on a predefined clinical grouping column or automatically screened across every column of a feature matrix (High/Low quantile split). |
| `compute.normalization()` | Normalizes raw counts to log-TPM (`ADImpute::NormalizeTPM`), with optional per-gene z-score scaling. |
| `compute.pca.analysis()` | PCA on a samples x features matrix, optionally colored by a clinical grouping column, with the PC1/PC2 plot saved to `Results/`. |
| `compute.correlation()` | Pearson + Spearman correlation, dispatching on input shape: two vectors, a score matrix vs. a continuous trait, or two matrices (with correlation heatmaps saved to `Results/`). |
| `scores.pearson.test()` | Pearson correlation test between each column of a score matrix and a continuous trait. |
| `scores.spearman.test()` | Spearman rank correlation test between each column of a score matrix and a continuous trait. |
| `scores.kendall.test()` | Kendall's tau correlation test between each column of a score matrix and a continuous trait. |
| `scores.lm.test()` | Linear regression F-test between each column of a score matrix and a continuous trait. |
| `scores.fisher.test.continuous()` | Fisher's exact test on a score matrix and a continuous trait, both binarized at the median. |
| `scores.stat.analysis.continuous()` | Unified dispatcher for the `scores.*.test()` family above, selected via a `method` argument. |

Each function is documented with a roxygen block above its definition in `functions.r`
(parameters, return value, and usage example), so refer there for full argument details.

A ready-to-run example of every function and its documented variations, tested on the
Vanderbilt lung cancer cohort from the LungPredict1 paper, is in
[`test_functions.qmd`](test_functions.qmd).

## Usage example

```r
library(dplyr)
source("functions.r")
dir.create("Results", showWarnings = FALSE)

# raw.counts: genes x samples raw count matrix
# coldata:    sample metadata, row names matching colnames(raw.counts)

counts.norm <- edgeR::cpm(raw.counts, log = TRUE)

tfs      <- compute.TFs.activity(counts.norm, TF.collection = "CollecTRI")
pathways <- compute.pathway.activity(counts.norm)
deg      <- run_deg_analysis(raw.counts, coldata, group_col = "Group", ref_level = "Control")

km <- compute.survival.analysis(coldata, PFS = "PFS", PFS_event = "PFS_event",
                                 features = pathways, p.value = 0.05)
```

## Contributing

New functions are welcome — just add them to `functions.r` documented with the same
[roxygen2](https://roxygen2.r-lib.org/) style already used throughout the file:

```r
#' One-line title
#'
#' A short paragraph describing what the function does and why.
#'
#' @param x Description of each parameter, including type and default behavior.
#' @param ... (one @param per argument)
#'
#' @return What the function returns.
#'
#' @examples
#' my_new_function(x)
#'
my_new_function <- function(x, ...) {
  ...
}
```

Save outputs (plots, matrices) to a `Results/` folder, following the existing functions'
convention, and add a usage example to `test_functions.qmd` if the function can be
demonstrated on the Vanderbilt cohort already loaded there.

## Acknowledgements

This repository is maintained by [Marcelo Hurtado](https://github.com/mhurtado13) in the
[Network Biology for Immuno-oncology (NetB(IO)²)](https://www.crct-inserm.fr/en/netbio2_en/)
group at the Cancer Research Center of Toulouse, in supervision of
[Vera Pancaldi](https://github.com/VeraPancaldi).
