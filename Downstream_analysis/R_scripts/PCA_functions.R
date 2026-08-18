# PCA and PERMANOVA analysis of overlapping taxa
# Input: taxon count matrices from Kraken2 or DIAMOND-MEGAN6
# Analysis: Hellinger transformation, PCA, PERMANOVA, and visualization

library(vegan)     # Hellinger transformation, PCA, PERMANOVA
library(ggplot2)   # plotting
library(ggrepel)   # non-overlapping labels
library(ggalt)     # geom_encircle

run_PCA_analysis <- function(
    count_file,
    overlap_file,
    metadata,
    permanova_rhs,
    permanova_by = NULL
) {
  # --------------------------------------------------------
  # 1. Read count matrix
  # --------------------------------------------------------
  
  micro_data <- read.delim(
    count_file,
    sep = "\t",
    header = FALSE,
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = ""
  )
  
  # Standardize taxon names
  rownames(micro_data) <- gsub("_", " ", rownames(micro_data))
  
  
  # --------------------------------------------------------
  # 2. Assign sample names
  # --------------------------------------------------------
  
  if (!"sample_id" %in% colnames(metadata)) {
    stop("metadata must contain a 'sample_id' column.")
  }
  
  if (ncol(micro_data) != nrow(metadata)) {
    stop(
      "The number of count-matrix columns does not match ",
      "the number of metadata rows."
    )
  }
  
  colnames(micro_data) <- metadata$sample_id
  
  
  # --------------------------------------------------------
  # 3. Retain overlapping taxa
  # --------------------------------------------------------
  
  overlap <- read.delim(
    overlap_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )[[1]]
  
  missing_taxa <- setdiff(overlap, rownames(micro_data))
  
  if (length(missing_taxa) > 0) {
    stop(
      "Some overlapping taxa were not found in the count matrix: ",
      paste(missing_taxa, collapse = ", ")
    )
  }
  
  micro_data <- micro_data[overlap, , drop = FALSE]
  
  
  # --------------------------------------------------------
  # 4. Transpose to sample × taxon
  # --------------------------------------------------------
  
  micro_data <- t(micro_data)
  
  
  # --------------------------------------------------------
  # 5. Match metadata to sample order
  # --------------------------------------------------------
  
  metadata_matched <- metadata[
    match(rownames(micro_data), metadata$sample_id),
    ,
    drop = FALSE
  ]
  
  if (anyNA(metadata_matched$sample_id)) {
    stop("Some samples could not be matched to metadata.")
  }
  
  
  # --------------------------------------------------------
  # 6. Hellinger transformation and PCA
  # --------------------------------------------------------
  
  micro_hel <- decostand(
    micro_data,
    method = "hellinger"
  )
  
  PCA <- rda(
    micro_hel,
    scale = FALSE
  )
  
  PCA_exp <- PCA$CA$eig / sum(PCA$CA$eig)
  
  site_scores <- as.data.frame(
    scores(
      PCA,
      display = "sites",
      scaling = 1,
      choices = 1:2
    )
  )
  
  taxa_scores <- as.data.frame(
    scores(
      PCA,
      display = "species",
      scaling = 2,
      choices = 1:2
    )
  )
  
  
  # --------------------------------------------------------
  # 7. Add metadata to sample scores
  # --------------------------------------------------------
  
  metadata_columns <- setdiff(
    colnames(metadata_matched),
    "sample_id"
  )
  
  site_scores <- cbind(
    site_scores,
    metadata_matched[, metadata_columns, drop = FALSE]
  )
  
  
  # --------------------------------------------------------
  # 8. PERMANOVA
  # --------------------------------------------------------
  
  if (
    !is.character(permanova_rhs) ||
    length(permanova_rhs) != 1 ||
    !nzchar(trimws(permanova_rhs))
  ) {
    stop("permanova_rhs must be a single non-empty character string.")
  }
  
  if (
    !is.null(permanova_by) &&
    !permanova_by %in% c("terms", "margin", "onedf")
  ) {
    stop(
      "permanova_by must be NULL, 'terms', 'margin', or 'onedf'."
    )
  }
  
  permanova_variables <- all.vars(
    as.formula(
      paste("~", permanova_rhs)
    )
  )
  
  missing_variables <- setdiff(
    permanova_variables,
    colnames(metadata_matched)
  )
  
  if (length(missing_variables) > 0) {
    stop(
      "PERMANOVA variable(s) not found in metadata: ",
      paste(missing_variables, collapse = ", ")
    )
  }
  
  permanova_formula <- as.formula(
    paste("micro_data ~", permanova_rhs)
  )
  
  permanova <- adonis2(
    permanova_formula,
    data = metadata_matched,
    permutations = 9999,
    method = "bray",
    by = permanova_by
  )
  
  
  # --------------------------------------------------------
  # 9. Return analysis results
  # --------------------------------------------------------
  
  return(list(
    count_matrix = micro_data,
    hellinger_matrix = micro_hel,
    PCA = PCA,
    explained_variance = PCA_exp,
    site_scores = site_scores,
    taxa_scores = taxa_scores,
    permanova = permanova,
    permanova_rhs = permanova_rhs,
    permanova_by = permanova_by
  ))
}


# ============================================================
# Create a PCA plot
# ============================================================

plot_PCA_result <- function(
    analysis_result,
    color_col,
    pipeline_name,
    organism,
    taxonomic_level,
    group_colors = NULL,
    shape_col = NULL,
    enclosure_method = "encircle",
    polygon_col = NULL,
    shape_values = NULL,
    show_sample_labels = TRUE
) {
  plot_data <- analysis_result$site_scores
  
  if (
    !is.logical(show_sample_labels) ||
    length(show_sample_labels) != 1 ||
    is.na(show_sample_labels)
  ) {
    stop("show_sample_labels must be TRUE or FALSE.")
  }

  if (!all(c("PC1", "PC2") %in% colnames(plot_data))) {
    stop("PC1 and PC2 were not found in the PCA site scores.")
  }
  
  if (!color_col %in% colnames(plot_data)) {
    stop(
      "Color column not found in PCA site scores: ",
      color_col
    )
  }
  
  if (
    !is.null(shape_col) &&
    !shape_col %in% colnames(plot_data)
  ) {
    stop(
      "Shape column not found in PCA site scores: ",
      shape_col
    )
  }
  
  if (
    !is.null(polygon_col) &&
    !polygon_col %in% colnames(plot_data)
  ) {
    stop(
      "Polygon column not found in PCA site scores: ",
      polygon_col
    )
  }
  
  plot_data$sample_id <- rownames(plot_data)
  plot_data$plot_color <- plot_data[[color_col]]
  
  if (!is.null(shape_col)) {
    plot_data$plot_shape <- plot_data[[shape_col]]
  }
  
  if (!is.null(polygon_col)) {
    plot_data$plot_polygon <- plot_data[[polygon_col]]
  }
  
  pc1_label <- paste0(
    "PC1: ",
    round(analysis_result$explained_variance[1] * 100, 2),
    "%"
  )
  
  pc2_label <- paste0(
    "PC2: ",
    round(analysis_result$explained_variance[2] * 100, 2),
    "%"
  )
  
  # --------------------------------------------------------
  # Build PERMANOVA subtitle text
  # --------------------------------------------------------
  
  permanova_table <- analysis_result$permanova
  
  if (is.null(analysis_result$permanova_by)) {
    
    model_row <- if (
      "Model" %in% rownames(permanova_table)
    ) {
      "Model"
    } else {
      rownames(permanova_table)[1]
    }
    
    permanova_text <- paste0(
      "PERMANOVA (",
      analysis_result$permanova_rhs,
      "): R2 = ",
      round(
        permanova_table[model_row, "R2"] * 100,
        2
      ),
      "%; p = ",
      signif(
        permanova_table[model_row, "Pr(>F)"],
        3
      )
    )
    
  } else {
    
    term_rows <- setdiff(
      rownames(permanova_table),
      c("Residual", "Total")
    )
    
    term_text <- vapply(
      term_rows,
      function(term) {
        paste0(
          term,
          ": R2 = ",
          round(
            permanova_table[term, "R2"] * 100,
            2
          ),
          "%; p = ",
          signif(
            permanova_table[term, "Pr(>F)"],
            3
          )
        )
      },
      character(1)
    )
    
    permanova_text <- paste(
      c(
        paste0(
          "PERMANOVA (",
          analysis_result$permanova_rhs,
          "):"
        ),
        term_text
      ),
      collapse = "\n"
    )
  }
  
  
  # --------------------------------------------------------
  # Check plotting colors
  # --------------------------------------------------------
  
  if (!is.null(group_colors)) {
    missing_colors <- setdiff(
      unique(as.character(plot_data$plot_color)),
      names(group_colors)
    )
    
    if (length(missing_colors) > 0) {
      stop(
        "No plotting colors were defined for: ",
        paste(missing_colors, collapse = ", ")
      )
    }
  }
  
  
  # --------------------------------------------------------
  # Create PCA plot
  # --------------------------------------------------------
  
  if (is.null(shape_col)) {
    point_layer <- geom_point(
      aes(color = plot_color),
      size = 1.5,
      alpha = 0.8
    )
  } else {
    point_layer <- geom_point(
      aes(
        color = plot_color,
        shape = plot_shape
      ),
      size = 1.5,
      alpha = 0.8
    )
  }
  
  p <- ggplot(
    plot_data,
    aes(PC1, PC2)
  ) +
    theme(
      panel.grid = element_line(
        color = "gray",
        linetype = 2,
        size = 0.1
      ),
      panel.background = element_rect(
        color = "black",
        fill = "transparent"
      ),
      legend.key = element_rect(
        fill = "transparent"
      )
    ) +
    labs(
      x = pc1_label,
      y = pc2_label,
      color = color_col,
      fill = color_col,
      shape = shape_col,
      title = paste0(
        organism,
        " (",
        pipeline_name,
        ")"
      ),
      subtitle = paste0(
        taxonomic_level,
        " level\n",
        permanova_text
      )
    ) +
    geom_vline(
      xintercept = 0,
      color = "gray",
      size = 0.4
    ) +
    geom_hline(
      yintercept = 0,
      color = "gray",
      size = 0.4
    )

  # Optional polygon layer (used for Place x genotype groups)
  if (!is.null(polygon_col)) {

    p <- p +
      geom_polygon(
        aes(
          group = plot_polygon,
          fill = plot_color
        ),
        alpha = 0.2,
        color = NA,
        show.legend = FALSE
      )
  }

  # Add points
  p <- p + point_layer

  # Optional enclosure
  if (enclosure_method == "encircle") {

    p <- p +
      ggalt::geom_encircle(
        aes(
          group = plot_color,
          fill = plot_color
        ),
        alpha = 0.2,
        show.legend = FALSE
      )

  } else if (enclosure_method == "ellipse") {

    p <- p +
      ggforce::geom_mark_ellipse(
        aes(
          group = plot_color,
          color = plot_color
        ),
        show.legend = FALSE
      )

  } else if (enclosure_method != "none") {

    stop(
      "enclosure_method must be 'encircle', 'ellipse', or 'none'."
    )
  }

  if (!is.null(group_colors)) {
    p <- p +
      scale_color_manual(values = group_colors) +
      scale_fill_manual(values = group_colors)
  }

  if (!is.null(shape_values)) {
    p <- p +
      scale_shape_manual(values = shape_values)
  }

  # Add sample labels last, if requested
  if (show_sample_labels) {
    p <- p +
      geom_text_repel(
        aes(label = sample_id),
        size = 1.5
      )
  }

  return(p)
}
# ============================================================
# Add top-abundance taxon arrows to an existing PCA plot
# ============================================================
#
# Taxa are ranked by total abundance across the same overlapping
# count matrix and the same samples used in the PCA.
#
# This reproduces the original all-potato arrow-selection logic
# without pipeline-specific deletion of sample columns.
#
add_taxa_arrows <- function(
    p,
    analysis_result,
    top_n = 10,
    arrow_scale = 0.9,
    label_x_scale = 0.75,
    label_y_scale = 1,
    arrow_color = "red",
    arrow_linewidth = 0.3,
    label_size = 2
) {

  # --------------------------------------------------------
  # 1. Validate inputs
  # --------------------------------------------------------

  if (!inherits(p, "ggplot")) {
    stop("p must be a ggplot object.")
  }

  if (
    !is.numeric(top_n) ||
    length(top_n) != 1 ||
    is.na(top_n) ||
    top_n < 1
  ) {
    stop("top_n must be a positive number.")
  }

  top_n <- as.integer(top_n)

  if (is.null(analysis_result$count_matrix)) {
    stop("analysis_result does not contain count_matrix.")
  }

  if (is.null(analysis_result$taxa_scores)) {
    stop("analysis_result does not contain taxa_scores.")
  }

  count_matrix <- analysis_result$count_matrix
  taxa_scores <- analysis_result$taxa_scores

  if (is.null(colnames(count_matrix))) {
    stop("count_matrix must contain taxon column names.")
  }

  if (!all(c("PC1", "PC2") %in% colnames(taxa_scores))) {
    stop("PC1 and PC2 were not found in taxa_scores.")
  }


  # --------------------------------------------------------
  # 2. Select the most abundant overlapping taxa
  # --------------------------------------------------------

  taxon_abundance <- colSums(
    count_matrix,
    na.rm = TRUE
  )

  ranked_taxa <- names(
    sort(
      taxon_abundance,
      decreasing = TRUE
    )
  )

  top_n <- min(
    top_n,
    length(ranked_taxa)
  )

  top_taxa <- ranked_taxa[
    seq_len(top_n)
  ]

  missing_scores <- setdiff(
    top_taxa,
    rownames(taxa_scores)
  )

  if (length(missing_scores) > 0) {
    stop(
      "PCA taxa scores were not found for: ",
      paste(missing_scores, collapse = ", ")
    )
  }


  # --------------------------------------------------------
  # 3. Prepare arrow coordinates
  # --------------------------------------------------------

  arrow_data <- taxa_scores[
    top_taxa,
    c("PC1", "PC2"),
    drop = FALSE
  ]

  arrow_data$taxon <- rownames(
    arrow_data
  )


  # --------------------------------------------------------
  # 4. Add arrows and taxon labels
  # --------------------------------------------------------

  p <- p +
    geom_segment(
      data = arrow_data,
      aes(
        x = 0,
        y = 0,
        xend = PC1 * arrow_scale,
        yend = PC2 * arrow_scale
      ),
      inherit.aes = FALSE,
      arrow = grid::arrow(
        length = grid::unit(
          0.1,
          "cm"
        )
      ),
      linewidth = arrow_linewidth,
      color = arrow_color
    ) +
    ggrepel::geom_text_repel(
      data = arrow_data,
      aes(
        x = PC1 * label_x_scale,
        y = PC2 * label_y_scale,
        label = taxon
      ),
      inherit.aes = FALSE,
      color = arrow_color,
      size = label_size,
      box.padding = 0.3,
      point.padding = 0.2,
      max.overlaps = Inf
    )

  return(p)
}
