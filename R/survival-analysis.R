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
