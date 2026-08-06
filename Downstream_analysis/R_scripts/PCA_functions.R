# PCA and PERMANOVA analysis of overlapping taxa
# Input: taxon count matrices from Kraken2 or DIAMOND-MEGAN6
# Analysis: Hellinger transformation, PCA, PERMANOVA, and visualization

library(vegan)     # Hellinger transformation, PCA, PERMANOVA
library(ggplot2)   # plotting
library(ggrepel)   # non-overlapping labels
library(patchwork) # plot combination
library(ggalt)     # geom_encircle
library(ggpubr)    # ggarrange


run_PCA_analysis <- function(
    count_file,
    overlap_file,
    metadata,
    group_col
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
  
  if (!group_col %in% colnames(metadata_matched)) {
    stop(
      "Grouping column not found in metadata: ",
      group_col
    )
  }
  
  permanova_data <- data.frame(
    group = metadata_matched[[group_col]]
  )
  
  permanova <- adonis2(
    micro_data ~ group,
    data = permanova_data,
    permutations = 9999,
    method = "bray"
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
    permanova = permanova
  ))
}

# ============================================================
# Create a PCA plot
# ============================================================

plot_PCA_result <- function(
    analysis_result,
    group_col,
    pipeline_name,
    organism,
    taxonomic_level,
    group_colors = NULL
) {
  plot_data <- analysis_result$site_scores

  if (!all(c("PC1", "PC2") %in% colnames(plot_data))) {
    stop("PC1 and PC2 were not found in the PCA site scores.")
  }

  if (!group_col %in% colnames(plot_data)) {
    stop(
      "Grouping column not found in PCA site scores: ",
      group_col
    )
  }

  plot_data$sample_id <- rownames(plot_data)
  plot_data$plot_group <- plot_data[[group_col]]

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

  permanova_r2 <- round(
    analysis_result$permanova$R2[1] * 100,
    2
  )

  permanova_p <- signif(
    analysis_result$permanova$`Pr(>F)`[1],
    3
  )

  if (!is.null(group_colors)) {
    missing_colors <- setdiff(
      unique(as.character(plot_data$plot_group)),
      names(group_colors)
    )

    if (length(missing_colors) > 0) {
      stop(
        "No plotting colors were defined for: ",
        paste(missing_colors, collapse = ", ")
      )
    }
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
      color = group_col,
      fill = group_col,
      title = paste0(
        organism,
        " (",
        pipeline_name,
        ")"
      ),
      subtitle = paste0(
        taxonomic_level,
        " level\n",
        "PERMANOVA: R2 = ",
        permanova_r2,
        "%; p = ",
        permanova_p
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
    ) +
    geom_point(
      aes(color = plot_group),
      size = 1.5,
      alpha = 0.8
    ) +
    geom_text_repel(
      aes(label = sample_id),
      size = 1.5
    ) +
    geom_encircle(
      aes(
        group = plot_group,
        fill = plot_group
      ),
      alpha = 0.2,
      show.legend = FALSE
    )

  if (!is.null(group_colors)) {
    p <- p +
      scale_color_manual(values = group_colors) +
      scale_fill_manual(values = group_colors)
  }

  return(p)
}

