#' Visualization of Cell Death ssGSEA Scores (Heatmap)
#'
#' Generate a heatmap of ssGSEA scores (pathways x samples) with optional
#' group annotation. This is the most common visualization for ssGSEA
#' results in publications.
#' Automatically adapt for two modes:
#' 1. Global mode: multiple pathways heatmap with row scaling
#' 2. Targeted mode: single-pathway heatmap without row scaling
#'
#' @param ssgsea_result A \code{celldeath_ssgsea} object returned by \code{celldeath_ssgsea()}.
#' @param show_category Number of top pathways to display (ordered by mean
#' absolute score). Default: 30.   # heatmap；bar 里 Default: 20
#' @param scale_rows Logical. Whether to z-score scale each pathway across samples
#' (recommended for multi-pathway mode). Default: TRUE.
#' @param cluster_cols Logical. Whether to cluster samples. Default: TRUE.
#' Ignored automatically if only one pathway is present.
#' @param palette Color palette for score values (low to high).
#' Default: c("#2E8B57", "white", "#E74C3C").
#' @param group_palette Named or unnamed colors for group annotation.
#' Default: c("#2E8B57", "#F39C12", "#E74C3C", "#8E44AD").
#' @param title Plot title. Default: "Cell Death Pathway ssGSEA Scores".
#' @param filename If provided, save plot to file (e.g., "heatmap.png"). Default: NULL.
#' @param width Plot width in inches when saving. Default: 10.
#' @param height Plot height in inches when saving. Default: 8.
#' @return Invisibly returns the \code{pheatmap} object.
#' @export
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
plot_death_ssgsea_heatmap <- function(ssgsea_result,
                                      show_category = 30,
                                      scale_rows = TRUE,
                                      cluster_cols = TRUE,
                                      palette = c("#2E8B57", "white", "#E74C3C"),
                                      group_palette = c("#2E8B57", "#F39C12", "#E74C3C", "#8E44AD"),
                                      title = "Cell Death Pathway ssGSEA Scores",
                                      filename = NULL, width = 10, height = 8) {
  # ========== Input validity check ==========
  if (!inherits(ssgsea_result, "celldeath_ssgsea")) {
    stop("Error: Parameter ssgsea_result must be the output of celldeath_ssgsea()!")
  }
  score_mat <- ssgsea_result$score
  group <- ssgsea_result$group

  if (nrow(score_mat) == 0 || ncol(score_mat) == 0) {
    stop("Error: The ssGSEA score matrix is empty!")
  }

  # Keep top pathways by mean absolute score to avoid label overlap
  if (nrow(score_mat) > show_category) {
    keep <- names(sort(rowMeans(abs(score_mat)), decreasing = TRUE))[1:show_category]
    score_mat <- score_mat[keep, , drop = FALSE]
  }

  # ========== Prepare plot matrix (row z-score for multi-pathway mode) ==========
  single_pathway <- nrow(score_mat) == 1
  plot_mat <- score_mat
  if (scale_rows && !single_pathway) {
    plot_mat <- t(scale(t(score_mat)))
    plot_mat[is.na(plot_mat)] <- 0   # constant-score pathway -> NA after scaling
  }

  # ========== Single-pathway adjustments ==========
  # Put the pathway name into the title (row label is hidden later)
  if (single_pathway) {
    title <- paste0(title, "\n", rownames(score_mat),
                    if (scale_rows) " (z-scored)" else "")
    # No clustering on a single row: order samples by group instead
    if (!is.null(group)) cluster_cols <- FALSE
  }

  # ========== Build group annotation ==========
  annotation_col <- NULL
  annotation_colors <- NULL
  if (!is.null(group)) {
    annotation_col <- data.frame(Group = group, row.names = colnames(plot_mat))
    group_levels <- levels(group)
    anno_cols <- group_palette[seq_len(length(group_levels))]
    names(anno_cols) <- group_levels
    annotation_colors <- list(Group = anno_cols)
    # Order samples by group for cleaner layout (disable clustering in this case)
    if (!cluster_cols) {
      plot_mat <- plot_mat[, order(group), drop = FALSE]
      annotation_col <- annotation_col[colnames(plot_mat), , drop = FALSE]
    }
  }

  # ========== Draw heatmap ==========
  col_fun <- grDevices::colorRampPalette(palette)(100)
  p <- pheatmap::pheatmap(
    plot_mat,
    color = col_fun,
    scale = "none",                      # scaling handled manually above
    cluster_rows = !single_pathway,
    cluster_cols = cluster_cols && ncol(plot_mat) > 2,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_colnames = ncol(plot_mat) <= 50,
    show_rownames = !single_pathway,
    fontsize_row = 9,
    main = title,
    silent = TRUE
  )

  # Save plot if filename is provided
  if (!is.null(filename)) {
    grDevices::png(filename, width = width, height = height, units = "in", res = 300)
    grid::grid.draw(p$gtable)
    grDevices::dev.off()
  } else {
    grid::grid.draw(p$gtable)
  }

  invisible(p)
}


#' Visualization of Cell Death ssGSEA Scores (Grouped Boxplot)
#'
#' Generate boxplots comparing ssGSEA score distributions across groups.
#' Best companion for \code{celldeath_ssgsea_diff()} results.
#' Automatically adapt for two modes:
#' 1. Global mode: pathways on x-axis (faceted or dodged by group)
#' 2. Targeted mode: single-pathway group comparison layout
#'
#' @param ssgsea_result A \code{celldeath_ssgsea} object returned by
#' \code{celldeath_ssgsea()} (must have non-NULL \code{group}).
#' @param show_category Number of top pathways to display (ordered by mean
#' absolute score). Default: 10.
#' @param palette Colors for groups. Default: c("#2E8B57", "#F39C12", "#E74C3C", "#8E44AD").
#' @param add_points Logical. Whether to overlay jittered sample points. Default: TRUE.
#' @param title Plot title. Default: "Cell Death Pathway ssGSEA Score Comparison".
#' @param filename If provided, save plot to file (e.g., "boxplot.png"). Default: NULL.
#' @param return_data If TRUE, return the data frame used for plotting. Default: FALSE.
#' @return A \code{ggplot2} object (or data frame if return_data = TRUE).
#' @export
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter facet_wrap scale_fill_manual
#' @importFrom ggplot2 labs theme_minimal theme element_text ggsave
plot_death_ssgsea_boxplot <- function(ssgsea_result,
                                      show_category = 10,
                                      palette = c("#2E8B57", "#F39C12", "#E74C3C", "#8E44AD"),
                                      add_points = TRUE,
                                      title = "Cell Death Pathway ssGSEA Score Comparison",
                                      filename = NULL, return_data = FALSE) {
  # ========== Input validity check ==========
  if (!inherits(ssgsea_result, "celldeath_ssgsea")) {
    stop("Error: Parameter ssgsea_result must be the output of celldeath_ssgsea()!")
  }
  if (is.null(ssgsea_result$group)) {
    stop("Error: No group information found! Please provide the group parameter when running celldeath_ssgsea().")
  }
  score_mat <- ssgsea_result$score
  group <- ssgsea_result$group

  # ========== Convert to long format ==========
  df <- data.frame(
    Sample = rep(colnames(score_mat), each = nrow(score_mat)),
    Pathway = rep(rownames(score_mat), times = ncol(score_mat)),
    Score = as.vector(score_mat),
    Group = rep(group, each = nrow(score_mat)),
    stringsAsFactors = FALSE
  )

  # Select top pathways by mean absolute score
  if (nrow(score_mat) > show_category) {
    keep <- names(sort(rowMeans(abs(score_mat)), decreasing = TRUE))[1:show_category]
    df <- df[df$Pathway %in% keep, ]
  }
  df$Pathway <- factor(df$Pathway,
                       levels = rownames(score_mat),
                       labels = .wrap_pathway_labels(rownames(score_mat), width = 22)
                       )

  # ========== Draw boxplot ==========
  single_pathway <- length(unique(df$Pathway)) == 1

  # Friendly message for targeted single-pathway mode
  if (single_pathway) {
    message("INFO: Current result contains a single pathway (targeted mode); showing group comparison for this pathway only.")
  }

  if (!single_pathway) {
    # Global mode: one facet per pathway, groups on x-axis
    p <- ggplot2::ggplot(df, ggplot2::aes(x = Group, y = Score, fill = Group)) +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      ggplot2::facet_wrap(~ Pathway, scales = "free_y") +
      ggplot2::scale_fill_manual(values = palette, name = "Group") +
      ggplot2::labs(x = NULL, y = "ssGSEA Score", title = title) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        strip.text = ggplot2::element_text(size = 8),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
  } else {
    # Targeted single-pathway mode: simple grouped boxplot with legend
    p <- ggplot2::ggplot(df, ggplot2::aes(x = Group, y = Score, fill = Group)) +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
      ggplot2::scale_fill_manual(values = palette, name = "Group") +
      ggplot2::labs(x = "Group", y = "ssGSEA Score",
                    title = paste0(title, " (", unique(df$Pathway), ")")) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.x = ggplot2::element_text(size = 10)
      )
  }

  # Overlay jittered sample points
  if (add_points) {
    p <- p + ggplot2::geom_jitter(width = 0.15, size = 1.2, alpha = 0.6)
  }

  # Save plot if filename is provided
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, p, width = 10, height = 8, dpi = 300)
  }

  # Return data or plot
  if (return_data) return(df) else return(p)
}


#' Visualization of Cell Death ssGSEA Scores (Mean Barplot)
#'
#' Generate a horizontal barplot of mean ssGSEA scores per group,
#' with SEM (standard error of mean) error bars.
#' Suitable for comparing overall pathway activity across groups
#' when many pathways are involved.
#'
#' @param ssgsea_res A \code{celldeath_ssgsea} object returned by \code{celldeath_ssgsea()}
#' (must have non-NULL \code{group}).
#' @param show_category Number of top pathways to display (ordered by mean
#' absolute score). Default: 30.   # heatmap；bar 里 Default: 20
#' @param palette Group colors. Default: c("#2E8B57", "#F39C12", "#E74C3C", "#8E44AD").
#' @param title Plot title. Default: "Mean Cell Death Pathway ssGSEA Score".
#' @param filename If provided, save plot to file (e.g., "barplot.png"). Default: NULL.
#' @param return_data If TRUE, return the data frame used for plotting. Default: FALSE.
#' @return A \code{ggplot2} object (or data frame if return_data = TRUE).
#' @export
#' @importFrom ggplot2 ggplot aes geom_bar geom_errorbar scale_fill_manual
#' @importFrom ggplot2 labs theme_minimal theme element_text ggsave coord_flip position_dodge
#' @importFrom stringr str_wrap
plot_death_ssgsea_bar <- function(ssgsea_res,
                                  show_category = 20,
                                  palette = c("#2E8B57", "#F39C12", "#E74C3C", "#8E44AD"),
                                  title = "Mean Cell Death Pathway ssGSEA Score",
                                  filename = NULL, return_data = FALSE) {
  # ========== Input validity check ==========
  if (!inherits(ssgsea_res, "celldeath_ssgsea")) {
    stop("Error: Parameter ssgsea_res must be the output of celldeath_ssgsea()!")
  }
  if (is.null(ssgsea_res$group)) {
    stop("Error: Mean barplot requires group information! Please provide the group parameter when running celldeath_ssgsea().")
  }

  score_mat <- ssgsea_res$score
  group <- ssgsea_res$group

  # Friendly message for targeted single-pathway mode
  if (nrow(score_mat) == 1) {
    message("INFO: Current result contains a single pathway (targeted mode); the barplot shows only this pathway. Consider plot_death_ssgsea_boxplot() for a clearer group comparison.")
  }

  # ========== Compute mean and SEM per pathway per group ==========
  df <- do.call(rbind, lapply(rownames(score_mat), function(pw) {
    do.call(rbind, lapply(levels(group), function(g) {
      vals <- score_mat[pw, group == g]
      data.frame(
        Pathway = pw,
        Group = g,
        Mean = mean(vals),
        SEM = stats::sd(vals) / sqrt(length(vals)),
        stringsAsFactors = FALSE
      )
    }))
  }))
  df$Pathway <- stringr::str_wrap(df$Pathway, width = 30)
  # Order pathways by overall mean score
  pathway_order <- names(sort(tapply(df$Mean, df$Pathway, mean)))
  df$Pathway <- factor(df$Pathway, levels = pathway_order)

  # Keep top pathways by mean absolute score
  unique_pathways <- unique(df$Pathway)
  if (length(unique_pathways) > show_category) {
    mean_abs <- tapply(df$Mean, df$Pathway, function(x) mean(abs(x)))
    keep <- names(sort(mean_abs, decreasing = TRUE))[1:show_category]
    df <- df[df$Pathway %in% keep, ]
  }

  # ========== Draw horizontal grouped barplot ==========
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Pathway, y = Mean, fill = Group)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.8),
                      alpha = 0.8, width = 0.7) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean - SEM, ymax = Mean + SEM),
                           position = ggplot2::position_dodge(width = 0.8), width = 0.25) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(x = "Cell Death Pathway", y = "Mean ssGSEA Score", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 9),
      legend.position = "top"
    )

  # Save plot if filename is provided
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, p, width = 9, height = 7, dpi = 300)
  }

  # Return data or plot
  if (return_data) return(df) else return(p)
}

