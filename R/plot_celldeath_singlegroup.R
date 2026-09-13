#' Visualization of Cell Death Enrichment Results (Bubble Plot)
#'
#' Generate a dotplot for single-group cell death pathway enrichment analysis.
#' Automatically adapt for two modes:
#' 1. Global mode: multiple pathways bubble plot
#' 2. Targeted mode: single-pathway ORA plot layout adjustment
#'
#' @param enrich_res An \code{enrichResult} object returned by \code{celldeath_enrich()}.
#' @param show_category Number of top pathways to display (sorted by p-value).
#' @param palette Color palette for p.adjust (low to high). Default: c("#2E8B57", "#F39C12", "#E74C3C").
#' @param size_range Range of point size (min, max). Default: c(2, 8).
#' @param title Plot title. Default: "Cell Death Pathway Enrichment".
#' @param filename If provided, save plot to file (e.g., "bubble.png"). Default: NULL.
#' @param return_data If TRUE, return the data frame used for plotting. Default: FALSE.
#' @return A \code{ggplot2} object (or data frame if return_data = TRUE).
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn scale_size_continuous
#' @importFrom ggplot2 labs theme_minimal theme element_text ggsave
#' @importFrom stringr str_wrap
#' @examples
#' data("demo_deg", package = "CDEnrich")
#' data("demo_all_genes", package = "CDEnrich")
#' data("representative_genes", package = "CDEnrich")
#'
#' ora_result <- celldeath_enrich(
#'   deg = demo_deg,
#'   universe = demo_all_genes
#' )
#'
#' plot_death_enrich_bubble(
#'   ora_result,
#'   show_category = 10
#' )
plot_death_enrich_bubble <- function(enrich_res, show_category = 10,
                                     palette = c("#2E8B57", "#F39C12", "#E74C3C"),
                                     size_range = c(2, 8),
                                     title = "Cell Death Pathway Enrichment",
                                     filename = NULL, return_data = FALSE) {
  # ========== Input validity check ==========
  if (is.null(enrich_res)) {
    stop("Please run celldeath_enrich() first to obtain valid enrichment results!")
  }

  # Prepare data
  df_all <- as.data.frame(enrich_res)

  if (nrow(df_all) == 0) {
    stop("Enrichment result is empty!")
  }

  is_targeted <- nrow(df_all) == 1

  if (is_targeted) {
    message(
      "Targeted single-pathway result detected; ",
      "using plot_death_enrich_targeted() instead of a multi-pathway bubble plot."
    )

    return(plot_death_enrich_targeted(
      enrich_res = enrich_res,
      color = palette[2],
      title = title,
      filename = filename,
      return_data = return_data
    ))
  }

  df <- df_all
  df <- df[order(df$p.adjust), ][seq_len(min(show_category, nrow(df))), ]
  df$Description <- stringr::str_wrap(df$Description, width = 40)

  # ========== Draw bubble plot ==========
  # Global ORA with multiple pathways: standard ordered bubble plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = GeneRatio, y = stats::reorder(Description, -p.adjust))) +
    ggplot2::geom_point(ggplot2::aes(size = Count, color = p.adjust), alpha = 0.8) +
    ggplot2::scale_color_gradientn(colors = palette, trans = "log10",
                                   name = "Adjusted P-value") +
    ggplot2::scale_size_continuous(range = size_range, name = "Gene Count") +
    ggplot2::labs(x = "Gene Ratio", y = "Cell Death Pathway", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 9)
    )

  # Save plot if filename is provided
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, p, width = 10, height = 8, dpi = 300)
  }

  # Return data or plot
  if (return_data) return(df) else return(p)
}




#' Targeted Single-Pathway ORA Summary Plot
#'
#' Generate a compact summary plot for a targeted ORA result containing
#' exactly one pathway. The bar shows the number of overlapping genes and
#' the annotation reports GeneRatio, BgRatio, and adjusted P-value.
#'
#' @param enrich_res An \code{enrichResult} object or data frame returned by
#'   \code{celldeath_enrich()}.
#' @param color Bar color. Default: "#F39C12".
#' @param title Plot title. Default: "Targeted Cell Death Pathway Enrichment".
#' @param filename If provided, save plot to file. Default: NULL.
#' @param return_data If TRUE, return the data frame used for plotting.
#'   Default: FALSE.
#'
#' @return A \code{ggplot2} object (or data frame if \code{return_data = TRUE}).
#' @export
#' @importFrom stringr str_wrap
#' @importFrom rlang .data
#' @examples
#' data("demo_targeted_deg", package = "CDEnrich")
#' data("demo_all_genes", package = "CDEnrich")
#' data("representative_pathways", package = "CDEnrich")
#'
#' ora_ferroptosis <- celldeath_enrich(
#'   deg = demo_targeted_deg,
#'   universe = demo_all_genes,
#'   term_select = "Ferroptosis"
#' )
#'
#' plot_death_enrich_targeted(ora_ferroptosis)
plot_death_enrich_targeted <- function(enrich_res,
                                       color = "#F39C12",
                                       title = "Targeted Cell Death Pathway Enrichment",
                                       filename = NULL,
                                       return_data = FALSE) {
  if (is.null(enrich_res)) {
    stop("Please run celldeath_enrich() first to obtain valid enrichment results!")
  }

  df <- as.data.frame(enrich_res)

  if (nrow(df) == 0) {
    stop("Enrichment result is empty!")
  }

  if (nrow(df) != 1) {
    stop(
      "plot_death_enrich_targeted() requires exactly one enrichment result. ",
      "For multiple pathways, use plot_death_enrich_bubble(), ",
      "plot_death_enrich_network(), or plot_death_volcano()."
    )
  }

  required_cols <- c(
    "Description", "GeneRatio", "BgRatio",
    "p.adjust", "Count"
  )
  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  df$Description <- stringr::str_wrap(df$Description, width = 40)
  df$Label <- sprintf(
    "GeneRatio: %s\nBgRatio: %s\nAdjusted P-value: %.3e",
    df$GeneRatio,
    df$BgRatio,
    df$p.adjust
  )

  if (return_data) {
    return(df)
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = Count, y = Description)
  ) +
    ggplot2::geom_col(fill = color, width = 0.45) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$Label),
      hjust = 0,
      nudge_x = max(1, max(df$Count) * 0.03),
      size = 4
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.65))
    ) +
    ggplot2::coord_cartesian(
      clip = "off"
    ) +
    ggplot2::labs(
      x = "Overlapping Gene Count",
      y = "Targeted Cell Death Pathway",
      title = title
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.ticks.y = ggplot2::element_blank()
    )

  if (!is.null(filename)) {
    ggplot2::ggsave(
      filename = filename,
      plot = p,
      width = 11,
      height = 4,
      dpi = 300
    )
  }

  return(p)
}





#' Visualization of Cell Death Enrichment Results (Enrichment Network)
#'
#' Generate an enrichment network plot for single-group cell death pathway enrichment analysis.
#' Automatically adapt two analysis modes from celldeath_enrich():
#' 1. Global mode: multiple pathways, full gene overlap network with edges
#' 2. Targeted mode: single specified pathway, simple single-node plot
#' Nodes represent pathways, edges represent gene overlaps (similarity)(only for multi-pathway).
#'
#' @param enrich_res An \code{enrichResult} object returned by \code{celldeath_enrich()}.
#' @param show_category Number of top pathways to display (sorted by p-value). Default: 30.
#' @param color Which variable to use for node color. Options: "p.adjust" (default), "pvalue", "qvalue", "Count", etc.
#' @param filename If provided, save plot to file (e.g., "enrich_network.png"). Default: NULL.
#' @param ... Additional arguments passed to \code{enrichplot::emapplot} (e.g., node_label, cex_label_group).
#'
#' @return A \code{ggplot2} object.
#' @export
#' @importFrom enrichplot pairwise_termsim emapplot
#' @importFrom stringr str_wrap
#' @importFrom ggplot2 ggsave ggplot aes geom_point scale_color_gradientn
#' @importFrom ggplot2 scale_size_continuous labs theme_minimal element_text
#' @examples
#' data("demo_deg", package = "CDEnrich")
#' data("demo_all_genes", package = "CDEnrich")
#' data("representative_genes", package = "CDEnrich")
#'
#' ora_result <- celldeath_enrich(
#'   deg = demo_deg,
#'   universe = demo_all_genes
#' )
#'
#' plot_death_enrich_network(
#'   ora_result,
#'   show_category = 10
#' )
plot_death_enrich_network <- function(enrich_res, show_category = 30,
                                      color = "p.adjust", filename = NULL, ...) {
  # ========== Input validity check ==========
  if (is.null(enrich_res)) {
    stop("Please run celldeath_enrich() first to obtain valid enrichment results!")
  }
  df <- as.data.frame(enrich_res)
  if (nrow(df) == 0) {
    stop("Enrichment result is empty, cannot draw enrichment network!")
  }

  # Targeted single-pathway ORA has no inter-pathway edges.
  # Use the compact targeted summary plot instead.
  if (nrow(df) == 1) {
    message(
      "Targeted single-pathway result detected; network plot requires at least two pathways. ",
      "Using plot_death_enrich_targeted() instead."
    )

    return(plot_death_enrich_targeted(
      enrich_res = enrich_res,
      title = "Single Targeted Pathway Enrichment",
      filename = filename
    ))
  }

  # Subset top significant pathways
  df <- df[order(df$p.adjust), ][seq_len(min(show_category, nrow(df))), ]

  # Validate color argument
  if (!color %in% colnames(df)) {
    warning("Parameter 'color' '", color, "' is not present in result columns. 'p.adjust' will be used for color mapping.")
    color <- "p.adjust"
  }

  # ========== Draw network plot ==========
  # Global multi-pathway mode: calculate similarity and generate connected network
  termsim <- enrichplot::pairwise_termsim(enrich_res)
  p <- enrichplot::emapplot(termsim,
                            showCategory = show_category,
                            color = color,
                            ...)

  # ========== Save plot if filename is provided ==========
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, p, width = 10, height = 8, dpi = 300)
  }

  return(p)
}




#' Volcano Plot for Cell Death Enrichment Results
#'
#' Generate a volcano-style plot for single-group enrichment analysis (ORA or GSEA),
#' showing effect size (GeneRatio for ORA, NES for GSEA) versus significance (-log10 p.adjust).
#' Support both global multi-pathway mode and targeted single-pathway mode.
#'
#' @param enrich_res An \code{enrichResult} (ORA) or \code{gseaResult} (GSEA) object.
#' @param label_top Number of top significant pathways to label (sorted by p.adjust). Default: 5.
#' @param color_point Color of non-significant points. Default: "#D3D3D3".
#' @param color_signif Color of significant points (p.adjust < 0.05). Default: "#D62728".
#' @param point_size Size of points. Default: 2.
#' @param alpha Transparency of points. Default: 0.8.
#' @param label_size Size of pathway labels. Default: 3.5.
#' @param x_lab X-axis label. If NULL, auto-detected ("Gene Ratio" for ORA, "NES" for GSEA). Default: NULL.
#' @param y_lab Y-axis label. Default: "-log10(Adjusted P-value)".
#' @param title Plot title. If NULL, auto-generated. Default: NULL.
#' @param filename If provided, save plot to file (e.g., "volcano.png"). Default: NULL.
#' @param ... Additional arguments passed to \code{ggrepel::geom_text_repel} (e.g., max.overlaps, force).
#'
#' @return A \code{ggplot2} object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs theme_bw element_text
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggsave
#' @examples
#' data("demo_deg", package = "CDEnrich")
#' data("demo_all_genes", package = "CDEnrich")
#' data("representative_genes", package = "CDEnrich")
#'
#' ora_result <- celldeath_enrich(
#'   deg = demo_deg,
#'   universe = demo_all_genes
#' )
#'
#' plot_death_volcano(
#'   ora_result,
#'   label_top = 5
#' )
plot_death_volcano <- function(enrich_res, label_top = 5,
                               color_point = "lightblue",
                               color_signif = "darkred",
                               point_size = 2,
                               alpha = 0.8,
                               label_size = 3.5,
                               x_lab = NULL,
                               y_lab = expression(-log[10]("Adjusted P-value")),
                               title = NULL,
                               filename = NULL,
                               ...) {
  # 1. Input validity check
  if (is.null(enrich_res)) {
    stop("enrich_res cannot be NULL, please run celldeath_enrich() or celldeath_gsea() first!")
  }

  # 2. Extract data frame (supports enrichResult/gseaResult and plain data frame)
  if (is.data.frame(enrich_res)) {
    df <- enrich_res
  } else {
    df <- as.data.frame(enrich_res)
  }
  if (nrow(df) == 0) {
    stop("Enrichment result is empty, cannot draw volcano plot!")
  }

  # 3. Detect analysis type and define effect size column
  if ("GeneRatio" %in% colnames(df)) {
    effect_col <- "GeneRatio"
    default_x_lab <- "Gene Ratio"
    type_name <- "ORA"
  } else if ("NES" %in% colnames(df)) {
    effect_col <- "NES"
    default_x_lab <- "Normalized Enrichment Score (NES)"
    type_name <- "GSEA"
  } else if (inherits(enrich_res, "enrichResult")) {
    effect_col <- "GeneRatio"
    default_x_lab <- "Gene Ratio"
    type_name <- "ORA"
  } else if (inherits(enrich_res, "gseaResult")) {
    effect_col <- "NES"
    default_x_lab <- "Normalized Enrichment Score (NES)"
    type_name <- "GSEA"
  } else {
    stop("Cannot identify enrichment result type: column 'GeneRatio' (ORA) or 'NES' (GSEA) missing, and input is not enrichResult/gseaResult object.")
  }

  # Check effect column exists
  if (!effect_col %in% colnames(df)) {
    stop("Cannot find effect size column:", effect_col)
  }

  # Set x-axis label
  if (is.null(x_lab)) x_lab <- default_x_lab

  # Single-pathway ORA results are better shown as a compact summary plot
  if (nrow(df) == 1 && effect_col == "GeneRatio") {
    message(
      "Targeted single-pathway ORA result detected; ",
      "using plot_death_enrich_targeted() instead of a one-point volcano plot."
    )

    targeted_color <- if (isTRUE(df$p.adjust[1] < 0.05)) {
      color_signif
    } else {
      color_point
    }

    return(plot_death_enrich_targeted(
      enrich_res = enrich_res,
      color = targeted_color,
      title = if (is.null(title)) {
        "Targeted Cell Death Pathway Enrichment"
      } else {
        title
      },
      filename = filename
    ))
  }

  # Auto-generate title
  if (is.null(title)) {
    title <- paste0(type_name, " Enrichment Volcano Plot")
  }

  # 4. Compute -log10(p.adjust)
  df$logp <- -log10(df$p.adjust)
  df$signif <- df$p.adjust < 0.05

  # 5. Prepare label data (top label_top entries sorted by p.adjust)
  label_df <- df[order(df$p.adjust), ][seq_len(min(label_top, nrow(df))), ]

  # Friendly message for single targeted GSEA pathway mode
  if (nrow(df) == 1) {
    message(
      "INFO: Targeted single-pathway GSEA result detected; ",
      "volcano plot displays only one data point for this specified pathway."
    )
  }

  # 6. Draw volcano plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[effect_col]], y = logp)) +
    ggplot2::geom_point(
      ggplot2::aes(color = signif),
      size = point_size,
      alpha = alpha,
      stroke = 0.2
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = color_point, "TRUE" = color_signif),
      name = "Significance",
      labels = c("p.adjust >= 0.05", "p.adjust < 0.05")
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = Description),
      size = label_size,
      box.padding = 0.35,
      point.padding = 0.5,
      segment.color = "grey50",
      show.legend = FALSE,
      ...
    ) +
    ggplot2::labs(x = x_lab, y = y_lab, title = title) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(colour = "black"),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 0.8)
    )

  # 7. Save plot
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, p, width = 8, height = 6, dpi = 300)
  }

  return(p)
}





#' Plot GSEA Curve for Cell Death Pathways
#'
#' Generate standard three-panel GSEA enrichment plot for cell death pathway.
#' Support both global multi-pathway mode and targeted single-pathway mode.
#'
#' @param gsea_res Result from \code{celldeath_gsea()}.
#' @param pathway Character. Either an **exact pathway ID** present in the
#' result, or a **cell death type name** (e.g. "Ferroptosis"), which is
#' resolved to the recommended representative gene set when available
#' (see \code{\link{get_representative_pathway}}).
#' @param use_recommended Logical. Only takes effect when \code{pathway} is
#' a cell death type with a curated recommendation: TRUE (default) plots
#' the recommended representative gene set; FALSE plots ALL gene sets of
#' this type present in the result.
#'
#' @return A ggplot2 object.
#' @export
#' @importFrom enrichplot gseaplot2
#' @importFrom ggplot2 theme element_text
#' @examples
#' data("demo_rank_apop", package = "CDEnrich")
#' data("representative_pathways", package = "CDEnrich")
#'
#' apoptosis_id <- get_representative_pathway("Apoptosis")$Pathway
#'
#' gsea_result <- celldeath_gsea(
#'   geneList = demo_rank_apop
#' )
#'
#' plot_death_gsea_curve(
#'   gsea_result,
#'   pathway = apoptosis_id
#' )
plot_death_gsea_curve <- function(gsea_res, pathway, use_recommended = TRUE) {
  # Convert to data frame for validation
  if (is.data.frame(gsea_res)) {
    df <- gsea_res
  } else {
    df <- as.data.frame(gsea_res)
  }
  # 1. Check empty result
  if (is.null(gsea_res) || nrow(df) == 0) {
    stop("GSEA result is empty, please check gene list or threshold.")
  }
  # 2. Resolve target pathway(s) element-wise: exact IDs pass through
  #    directly; a death type resolves to its recommended set(s)
  avail_ids <- unique(df$ID)
  pathway <- unlist(lapply(pathway, function(p) {
    .resolve_plot_pathway(p, avail_ids, use_recommended = use_recommended)
  }), use.names = FALSE)

  # 3. Draw standard three-panel GSEA plot
  p <- enrichplot::gseaplot2(
    gsea_res,
    geneSetID = pathway,
    title = paste(pathway, collapse = " vs "),
    color = "#2E8B57",
    base_size = 10
  )

  return(p)
}




#' Plot GSEA Ranking Score Curve for Cell Death Pathways
#'
#' Generate a running score plot (ranked list metric) for a single GSEA result.
#' Compatible with two modes from celldeath_gsea():
#' 1. Global mode: ranked list vs all cell death pathways
#' 2. Targeted mode: ranked list vs one specific single pathway
#' This plot shows the enrichment score as a function of the ranked gene list.
#'
#' @param gsea_res Result from \code{celldeath_gsea()}.
#' @param pathway Character. Either an **exact pathway ID** present in the
#' result, or a **cell death type name** (e.g. "Ferroptosis"), which is
#' resolved to the recommended representative gene set when available
#' (see \code{\link{get_representative_pathway}}).
#' @param use_recommended Logical. Only takes effect when \code{pathway} is
#' a cell death type with a curated recommendation: TRUE (default) plots
#' the recommended representative gene set; FALSE plots ALL gene sets of
#' this type present in the result.
#' @param color Color of the running score line. Default: "#2E8B57".
#' @param base_size Base font size. Default: 10.
#' @param title Plot title. If NULL, uses the pathway name. Default: NULL.
#' @param filename If provided, save plot to file (e.g., "gsea_ranked.png"). Default: NULL.
#' @param ... Additional arguments passed to \code{enrichplot::gseaplot}.
#'
#' @return A \code{ggplot2} object.
#' @export
#' @importFrom enrichplot gseaplot
#' @importFrom ggplot2 ggsave
#' @examples
#' data("demo_rank_apop", package = "CDEnrich")
#' data("representative_pathways", package = "CDEnrich")
#'
#' apoptosis_id <- get_representative_pathway("Apoptosis")$Pathway
#'
#' gsea_result <- celldeath_gsea(
#'   geneList = demo_rank_apop
#' )
#'
#' plot_death_gsea_ranked(
#'   gsea_result,
#'   pathway = apoptosis_id
#' )
#'
plot_death_gsea_ranked <- function(gsea_res, pathway,
                                   use_recommended = TRUE,
                                   color = "#2E8B57",
                                   base_size = 10,
                                   title = NULL,
                                   filename = NULL,
                                   ...) {

  # Convert to data frame for validation
  if (is.data.frame(gsea_res)) {
    df <- gsea_res
  } else {
    df <- as.data.frame(gsea_res)
  }

  # 1. Input validity check
  if (is.null(gsea_res) || nrow(df) == 0) {
    stop("GSEA result is empty, please check gene list or threshold.")
  }

  # 2. Resolve target pathway (exact ID, or death type -> recommended set)
  avail_path <- unique(df$ID)
  pathway <- unlist(lapply(pathway, function(p) {
    .resolve_plot_pathway(p, avail_path, use_recommended = use_recommended)
  }), use.names = FALSE)

  # 3. Set plot title
  plot_title <- if (is.null(title)) paste(pathway, collapse = " vs ") else title

  # 4. Draw running score curve
  p <- enrichplot::gseaplot(gsea_res,
                            geneSetID = pathway,
                            by = "runningScore",
                            color.line = color,
                            base_size = base_size,
                            title = plot_title,
                            ...)

  # 5. Save plot
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, p, width = 7, height = 5, dpi = 300)
  }

  return(p)
}
