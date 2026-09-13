#' Visualization of Multi-Group Cell Death Enrichment Comparison
#'
#' Generate a dotplot for multi-group cell death pathway enrichment comparison.
#' Support two analysis modes from celldeath_compare_enrich():
#' 1. Global mode: all groups compare against all cell death pathways
#' 2. Targeted mode: all groups compare against one single specified pathway
#'
#' @param compare_res A \code{compareClusterResult} object returned by \code{celldeath_compare_enrich()}.
#' @param show_category Number of top pathways to display per group (sorted by p-value).
#' @param palette Color palette for p.adjust (low to high). Default: c("#2E8B57", "#F39C12", "#E74C3C").
#' @param title Plot title. Default: "Cell Death Pathway Enrichment (Multi-Group)".
#' @param filename If provided, save plot to file (e.g., "compare_bubble.png"). Default: NULL.
#' @param return_data If TRUE, return the data frame used for plotting. Default: FALSE.
#' @return A \code{ggplot2} object (or data frame if return_data = TRUE).
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn scale_size_continuous
#' @importFrom ggplot2 facet_grid vars labs theme_minimal theme element_text unit ggsave
#' @importFrom stringr str_wrap
#' @examples
#' data("demo_deg_list", package = "CDEnrich")
#' data("demo_all_genes", package = "CDEnrich")
#' data("representative_genes", package = "CDEnrich")
#'
#' compare_result <- celldeath_compare_enrich(
#'   deg_list_list = demo_deg_list,
#'   universe = demo_all_genes
#' )
#'
#' plot_death_compare_bubble(
#'   compare_result,
#'   show_category = 10
#' )
plot_death_compare_bubble <- function(compare_res, show_category = 10,
                                      palette = c("#2E8B57", "#F39C12", "#E74C3C"),
                                      title = "Cell Death Pathway Enrichment (Multi-Group)",
                                      filename = NULL, return_data = FALSE) {
  # 1. Input validity check
  if (is.null(compare_res)) {
    stop("Please run celldeath_compare_enrich() first to obtain valid multi-group enrichment results!")
  }
  df <- as.data.frame(compare_res)
  if (nrow(df) == 0) {
    stop("Multi-group enrichment result is empty, cannot generate visualization!")
  }

  # 2. Subset top pathways for each group independently
  df_list <- split(df, df$Cluster)
  top_df <- do.call(rbind, lapply(df_list, function(x) {
    x <- x[order(x$p.adjust), ][seq_len(min(show_category, nrow(x))), ]
    x$Description <- stringr::str_wrap(x$Description, width = 40)
    return(x)
  }))

  # Friendly message for single targeted pathway mode
  unique_path_num <- length(unique(top_df$Description))
  if (unique_path_num == 1) {
    message("INFO: Current multi-group comparison uses targeted single-pathway mode; plot only shows enrichment results of this specified pathway across groups.")
  }

  # 3. Preserve original group order as factor levels
  top_df$Cluster <- factor(top_df$Cluster, levels = unique(df$Cluster))

  # 4. Draw faceted bubble plot
  p <- ggplot2::ggplot(top_df, ggplot2::aes(x = GeneRatio, y = stats::reorder(Description, -p.adjust))) +
    ggplot2::geom_point(ggplot2::aes(size = Count, color = p.adjust), alpha = 0.8) +
    ggplot2::scale_color_gradientn(colors = palette, trans = "log10",
                                   name = "Adjusted P-value") +
    ggplot2::scale_size_continuous(range = c(2, 8), name = "Gene Count") +
    ggplot2::facet_grid(rows = ggplot2::vars(Cluster), scales = "free_y", space = "free_y") +
    ggplot2::labs(x = "Gene Ratio", y = "Cell Death Pathway", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 9),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      panel.spacing = ggplot2::unit(0.5, "cm")
    )

  # 5. Save plot
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, p, width = 12, height = 10, dpi = 300)
  }

  # 6. Return data or plot
  if (return_data) return(top_df) else return(p)
}





#' Visualization of Multi-Group Cell Death Enrichment Comparison (Heatmap)
#'
#' Generate a heatmap for multi-group cell death pathway enrichment comparison,
#' showing -log10(p.adjust) for each pathway-group pair.
#' Compatible with two analysis modes from celldeath_compare_enrich():
#' 1. Global mode: all groups compare against all cell death pathways
#' 2. Targeted mode: all groups compare against one single specified pathway
#'
#' @param compare_res A \code{compareClusterResult} object returned by \code{celldeath_compare_enrich()}.
#' @param show_category Number of top pathways to display (selected by lowest p.adjust across groups).
#' @param palette Color palette for fill (low to high). Default: c("#2E8B57", "#F39C12", "#E74C3C").
#' @param title Plot title. Default: "Cell Death Pathway Enrichment (Multi-Group Heatmap)".
#' @param filename If provided, save plot to file (e.g., "compare_heatmap.png"). Default: NULL.
#' @param return_data If TRUE, return the data frame used for plotting. Default: FALSE.
#' @return A \code{ggplot2} object (or data frame if return_data = TRUE).
#' @export
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn labs theme_minimal
#' @importFrom ggplot2 theme element_text element_rect unit ggsave
#' @importFrom stringr str_wrap
#' @examples
#' data("demo_deg_list", package = "CDEnrich")
#' data("demo_all_genes", package = "CDEnrich")
#' data("representative_genes", package = "CDEnrich")
#'
#' compare_result <- celldeath_compare_enrich(
#'   deg_list_list = demo_deg_list,
#'   universe = demo_all_genes
#' )
#'
#' plot_death_compare_heatmap(
#'   compare_result,
#'   show_category = 10
#' )
plot_death_compare_heatmap <- function(compare_res, show_category = 10,
                                       palette = c("#2E8B57", "#F39C12", "#E74C3C"),
                                       title = "Cell Death Pathway Enrichment (Multi-Group Heatmap)",
                                       filename = NULL, return_data = FALSE) {
  # 1. Input validity check
  if (is.null(compare_res)) {
    stop("Please run celldeath_compare_enrich() first to obtain valid multi-group enrichment results!")
  }
  df <- as.data.frame(compare_res)
  if (nrow(df) == 0) {
    stop("Multi-group enrichment result is empty, cannot generate visualization!")
  }

  # 2. Calculate -log10(p.adjust) as fill value
  df$logp <- -log10(df$p.adjust)

  # 3. Select top pathways
  min_p_by_pathway <- tapply(df$p.adjust, df$Description, min, na.rm = TRUE)
  sorted_pathways <- names(sort(min_p_by_pathway))
  pathway_order <- sorted_pathways[seq_len(min(show_category, length(sorted_pathways)))]

  # 4. Filter data and wrap pathway labels
  heat_df <- df[df$Description %in% pathway_order, ]
  heat_df$Description <- stringr::str_wrap(heat_df$Description, width = 40)

  # Friendly message for single targeted pathway mode
  unique_path_count <- length(unique(heat_df$Description))
  if (unique_path_count == 1) {
    message("INFO: Current multi-group comparison heatmap uses targeted single-pathway mode; chart only shows significance tile of this specified pathway across groups.")
  }

  # 5. Set factor levels
  heat_df$Description <- factor(heat_df$Description, levels = stringr::str_wrap(pathway_order, width = 40))
  heat_df$Cluster <- factor(heat_df$Cluster, levels = unique(df$Cluster))

  # 6. Draw heatmap
  p <- ggplot2::ggplot(heat_df, ggplot2::aes(x = Cluster, y = Description, fill = logp)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradientn(
      colors = palette,
      name = expression(-log[10]("p.adjust")),
      na.value = "grey90"
    ) +
    ggplot2::labs(x = "Group", y = "Cell Death Pathway", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(size = 9),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.grid = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 9)
    )

  # 7. Save plot
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, p, width = 10, height = max(6, nrow(heat_df) * 0.3), dpi = 300)
  }

  # 8. Return data or plot
  if (return_data) return(heat_df) else return(p)
}





#' Plot GSEA Curves for Multiple Group Cell Death Pathways
#'
#' @param gsea_res_list Named list of GSEA results from celldeath_gsea_multiple().
#' @param pathway Character. Either an **exact pathway ID** present in the
#' results, or a **cell death type name** (e.g. "Ferroptosis"), which is
#' resolved to the recommended representative gene set when available
#' (see \code{\link{get_representative_pathway}}).
#' @param use_recommended Logical. Only takes effect when \code{pathway} is
#' a cell death type with a curated recommendation: TRUE (default) plots
#' the recommended representative gene set; FALSE plots ALL gene sets of
#' this type present in each group's result.
#' @param filename_prefix If provided, save each plot as "prefix_groupname.png". Default: NULL.
#' @param ... Additional arguments passed to \code{plot_death_gsea_curve}.
#'
#' @return A named list of ggplot2 objects, each containing one GSEA curve.
#' @export
#' @importFrom ggplot2 ggsave
#' @examples
#' data("demo_rank_apop_list", package = "CDEnrich")
#' data("representative_pathways", package = "CDEnrich")
#'
#' apoptosis_id <- get_representative_pathway("Apoptosis")$Pathway
#'
#' gsea_groups <- celldeath_gsea_multiple(
#'   geneList_list = demo_rank_apop_list,
#'   pvalueCutoff = 1
#' )
#'
#' curves <- plot_death_gsea_curve_multiple(
#'   gsea_groups,
#'   pathway = apoptosis_id
#' )
#'
#' curves$Control
#' curves$Treatment
plot_death_gsea_curve_multiple <- function(gsea_res_list, pathway,
                                           use_recommended = TRUE,
                                           filename_prefix = NULL, ...) {
  # 1. Check input is a named list
  if (is.null(names(gsea_res_list)) || any(names(gsea_res_list) == "")) {
    stop("gsea_res_list must be a NAMED list!")
  }

  # 2. Collect the union of pathway IDs across all non-empty group results
  ids_union <- unique(unlist(lapply(gsea_res_list, function(res) {
    if (is.null(res)) return(NULL)
    df_res <- if (is.data.frame(res)) res else as.data.frame(res)
    if (nrow(df_res) == 0) return(NULL)
    unique(df_res$ID)
  }), use.names = FALSE))

  if (length(ids_union) == 0) {
    stop("All GSEA results are empty, nothing to plot.")
  }

  # 3. Resolve pathway ONCE against the union of result IDs
  #    (exact ID -> used directly; death type -> recommended set, with notice)
  resolved_pathways <- .resolve_plot_pathway(pathway, ids_union,
                                             use_recommended = use_recommended)

  plot_list <- list()

  for (group_name in names(gsea_res_list)) {
    message("Plotting GSEA curve for group: ", group_name, " ...")

    current_res <- gsea_res_list[[group_name]]

    # Support data.frame and gseaResult object
    if (is.data.frame(current_res)) {
      df_res <- current_res
    } else {
      df_res <- as.data.frame(current_res)
    }

    # Check empty result
    if (is.null(current_res) || nrow(df_res) == 0) {
      warning("Group [", group_name, "] GSEA result is empty, skip plotting.")
      next
    }

    # Check which resolved pathways are present in this group
    present <- resolved_pathways[resolved_pathways %in% df_res$ID]
    if (length(present) == 0) {
      warning("Group [", group_name, "] does not contain pathway [",
                     paste(resolved_pathways, collapse = ", "), "], skip plotting.")
      next
    }

    # Pass the already-resolved exact IDs: skips per-group re-resolution,
    # respects use_recommended, and avoids duplicate messages
    p <- plot_death_gsea_curve(current_res, pathway = present, ...)

    plot_list[[group_name]] <- p

    # Batch save plots if prefix provided
    if (!is.null(filename_prefix)) {
      filename <- paste0(filename_prefix, "_", group_name, ".png")
      ggplot2::ggsave(filename, p, width = 8, height = 6, dpi = 300)
      message("Plot saved: ", filename)
    }
  }

  message("All group GSEA plots completed!")
  return(plot_list)
}




#' Heatmap for Multi-Group GSEA Results
#'
#' Generate a heatmap to compare pathway statistics (e.g., NES, -log10 p.adjust)
#' across multiple GSEA groups.
#' Compatible with two modes from celldeath_gsea():
#' 1. Global mode: ranked list vs all cell death pathways
#' 2. Targeted mode: ranked list vs one specific single pathway
#' Missing pathway in certain group will be shown as grey tile.
#'
#' @param gsea_res_list A named list of GSEA result objects (from \code{celldeath_gsea_multiple}).
#' @param show_category Number of top pathways to display. Default: 10.
#' @param statistic Which statistic to plot. Options: "NES" (default) or "p.adjust".
#' @param palette Color palette for fill (low to high). Default: c("#2E8B57", "#F39C12", "#E74C3C").
#' @param title Plot title. Default: "GSEA Multi-Group Enrichment Heatmap".
#' @param filename If provided, save plot to file. Default: NULL.
#' @param return_data If TRUE, return the data frame used for plotting. Default: FALSE.
#' @return A \code{ggplot2} object (or data frame if return_data = TRUE).
#' @export
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn labs theme_minimal
#' @importFrom ggplot2 theme element_text element_rect unit ggsave
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_wrap
#' @examples
#' data("demo_rank_apop_list", package = "CDEnrich")
#' data("representative_genes", package = "CDEnrich")
#'
#' gsea_groups <- celldeath_gsea_multiple(
#'   geneList_list = demo_rank_apop_list,
#'   pvalueCutoff = 1
#' )
#'
#' plot_death_gsea_heatmap_multiple(
#'   gsea_groups,
#'   statistic = "NES",
#'   show_category = 10
#' )
plot_death_gsea_heatmap_multiple <- function(gsea_res_list,
                                            show_category = 10,
                                            statistic = c("NES", "p.adjust"),
                                            palette = c("#2E8B57", "#F39C12", "#E74C3C"),
                                            title = "GSEA Multi-Group Enrichment Heatmap",
                                            filename = NULL,
                                            return_data = FALSE) {

  statistic <- match.arg(statistic)

  # 1. Input validity check
  if (is.null(names(gsea_res_list)) || any(names(gsea_res_list) == "")) {
    stop("gsea_res_list must be a NAMED list (each group must have a name)!")
  }

  # 2. Extract statistics from each group
  group_names <- names(gsea_res_list)
  all_data <- list()
  for (grp in group_names) {
    res <- gsea_res_list[[grp]]

    if (is.data.frame(res)) {
      df_res <- res
    } else {
      df_res <- as.data.frame(res)
    }
    if (is.null(res) || nrow(df_res) == 0) {
      warning("Group ", grp, " has no result, skipped.")
      next
    }
    df_grp <- df_res[, c("ID", "Description", statistic)]
    colnames(df_grp)[3] <- grp
    all_data[[grp]] <- df_grp
  }
  if (length(all_data) == 0) stop("No valid GSEA results.")
  group_names <- names(all_data)   # keep only groups with valid results

  # Merge all groups
  merged <- Reduce(function(x, y) merge(x, y, by = c("ID", "Description"), all = TRUE), all_data)
  rownames(merged) <- merged$ID

  # 3. Filter top pathways
  if (statistic == "p.adjust") {
    p_cols <- group_names[group_names %in% colnames(merged)]
    merged$min_p <- apply(merged[, p_cols, drop = FALSE], 1, min, na.rm = TRUE)
    merged <- merged[order(merged$min_p), ]
    for (grp in group_names) {
      merged[[grp]] <- -log10(merged[[grp]])
    }
    legend_name <- expression(-log[10]("p.adjust"))
  } else { # NES
    nes_cols <- group_names[group_names %in% colnames(merged)]
    merged$max_abs <- apply(abs(merged[, nes_cols, drop = FALSE]), 1, max, na.rm = TRUE)
    merged <- merged[order(-merged$max_abs), ]
    legend_name <- "NES"
  }

  pathway_order <- utils::head(rownames(merged), show_category)
  heat_df <- merged[pathway_order, group_names, drop = FALSE]

  # 4. Reshape to long format
  heat_df_with_id <- data.frame(ID = rownames(heat_df), heat_df, check.names = FALSE)
  heat_long <- tidyr::pivot_longer(heat_df_with_id,
                                   cols = -ID,
                                   names_to = "Group",
                                   values_to = "Value")

  # Attach wrapped pathway description
  desc_map <- merged[pathway_order, c("ID", "Description")]
  heat_long$Description <- desc_map$Description[match(heat_long$ID, desc_map$ID)]
  heat_long$Description <- stringr::str_wrap(heat_long$Description, width = 40)

  # Friendly message for single-pathway targeted comparison
  unique_path_count <- length(unique(heat_long$Description))
  if (unique_path_count == 1) {
    message("INFO: This multi-group GSEA comparison contains only one pathway, corresponding to targeted single-pathway cross-group analysis.")
  }

  # 5. Factor ordering
  heat_long$Description <- factor(heat_long$Description,
                                  levels = stringr::str_wrap(pathway_order, width = 40))
  heat_long$Group <- factor(heat_long$Group, levels = group_names)

  # 6. Draw heatmap
  p <- ggplot2::ggplot(heat_long, ggplot2::aes(x = Group, y = Description, fill = Value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradientn(
      colors = palette,
      name = legend_name,
      na.value = "grey90"
    ) +
    ggplot2::labs(x = "Group", y = "Cell Death Pathway", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(size = 9),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.grid = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 9)
    )

  # 7. Save plot
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, p, width = 10, height = max(6, nrow(heat_long) * 0.3), dpi = 300)
  }

  if (return_data) return(heat_long) else return(p)
}
