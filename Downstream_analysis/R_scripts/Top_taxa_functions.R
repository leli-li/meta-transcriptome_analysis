# Rank and plot taxa by group using one count matrix at a time.
# Sum counts across samples in each group.
# The overlap list marks shared taxa but does not filter the matrix.

summarize_top_taxa <- function(
    count_file,
    metadata,
    group_col,
    overlap_file = NULL,
    sample_id_transform = identity,
    top_n = 10
) {
  # 1. Read taxon names from column 1 and sample IDs from the header.
  counts <- read.delim(
    count_file, header = TRUE, row.names = 1,
    sep = "\t", check.names = FALSE, quote = "",
    stringsAsFactors = FALSE
  )
  rownames(counts) <- gsub("_", " ", rownames(counts))
  colnames(counts) <- as.character(sample_id_transform(colnames(counts)))

  if (nrow(counts) == 0 || ncol(counts) == 0 ||
      anyDuplicated(rownames(counts)) || any(trimws(rownames(counts)) == "")) {
    stop("Count matrix must contain unique, non-empty taxon names.")
  }
  if (anyNA(colnames(counts)) || any(trimws(colnames(counts)) == "") ||
      anyDuplicated(colnames(counts))) {
    stop("Count-matrix sample IDs must be unique and non-missing.")
  }
  if (!all(vapply(counts, is.numeric, logical(1)))) {
    stop("All sample columns must contain numeric counts.")
  }
  counts <- as.matrix(counts)
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("Counts must be finite and non-negative.")
  }
  if (length(top_n) != 1 || !is.numeric(top_n) || !is.finite(top_n) ||
      top_n < 1 || top_n != floor(top_n)) {
    stop("top_n must be a positive integer.")
  }

  # 2. Match metadata by sample ID.
  if (!all(c("sample_id", group_col) %in% colnames(metadata))) {
    stop("Metadata must contain sample_id and the requested group column.")
  }
  for (column in c("sample_id", group_col)) {
    if (anyNA(metadata[[column]]) || any(trimws(metadata[[column]]) == "")) {
      stop("Missing or empty metadata values in: ", column)
    }
  }
  if (anyDuplicated(metadata$sample_id) ||
      !setequal(colnames(counts), metadata$sample_id)) {
    stop("Count-matrix sample IDs must match unique metadata sample_id values.")
  }
  groups <- unique(as.character(metadata[[group_col]]))
  metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]

  # 3. Mark shared taxa, or use NA if there is no list.
  shared <- rep(NA, nrow(counts))
  if (!is.null(overlap_file)) {
    overlap <- trimws(readLines(overlap_file, warn = FALSE))
    overlap <- gsub("_", " ", overlap[nzchar(overlap)])
    if (!all(overlap %in% rownames(counts))) {
      stop("Some overlap taxa are absent from this count matrix.")
    }
    shared <- rownames(counts) %in% overlap
  }

  # 4. Sum and rank counts within each group.
  ranked_list <- list()
  top_list <- list()
  plot_list <- list()

  for (group in groups) {
    sample_ids <- metadata$sample_id[metadata[[group_col]] == group]
    total_count <- rowSums(counts[, sample_ids, drop = FALSE])
    group_total <- sum(total_count)
    if (group_total == 0) {
      stop("All counts are zero in group: ", group)
    }

    group_table <- data.frame(
      group = group,
      taxon = rownames(counts),
      total_count = unname(total_count),
      percentage = unname(100 * total_count / group_total),
      shared_between_pipelines = shared,
      stringsAsFactors = FALSE
    )
    # Break ties by taxon name.
    group_table <- group_table[
      order(-group_table$total_count, group_table$taxon), , drop = FALSE
    ]
    group_table$rank <- seq_len(nrow(group_table))
    ranked_list[[group]] <- group_table

    # Exclude zero-count taxa from the Top N.
    top_table <- head(group_table[group_table$total_count > 0, ], top_n)
    top_list[[group]] <- top_table

    # Use all taxa in the group to calculate percentages.
    plot_table <- top_table[, c("group", "taxon", "total_count", "percentage")]
    other_count <- group_total - sum(top_table$total_count)
    if (other_count > 0) {
      plot_table <- rbind(plot_table, data.frame(
        group = group, taxon = "Others", total_count = other_count,
        percentage = 100 * other_count / group_total
      ))
    }
    plot_list[[group]] <- plot_table
  }

  ranked_counts <- do.call(rbind, ranked_list)
  top_taxa <- do.call(rbind, top_list)
  plot_data <- do.call(rbind, plot_list)
  rownames(ranked_counts) <- rownames(top_taxa) <- rownames(plot_data) <- NULL

  return(list(
    ranked_counts = ranked_counts,
    top_taxa = top_taxa,
    plot_data = plot_data,
    top_n = top_n
  ))
}


plot_top_taxa <- function(result, title = "", count_unit = "counts", taxon_colors = NULL) {
  # "Others" can contain different taxa in different groups.
  # Label counts as read pairs for Kraken2 or contigs for DM.
  plot_data <- result$plot_data
  taxa <- unique(result$top_taxa$taxon)
  plot_data$group <- factor(plot_data$group, levels = unique(plot_data$group))
  plot_data$taxon <- factor(plot_data$taxon, levels = c(taxa, "Others"))

  if (is.null(taxon_colors)) {
    taxon_colors <- setNames(grDevices::hcl.colors(length(taxa), "Dark 3"), taxa)
  }
  if (!all(taxa %in% names(taxon_colors))) {
    stop("taxon_colors must supply a named color for every displayed taxon.")
  }
  taxon_colors["Others"] <- "#D9D9D9"

  ggplot2::ggplot(plot_data, ggplot2::aes(x = group, y = percentage, fill = taxon)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = taxon_colors, drop = TRUE) +
    ggplot2::scale_y_continuous(breaks = seq(0, 100, 25)) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Top ", result$top_n, " taxa per group + Others"),
      x = NULL, y = paste0("Percentage of ", count_unit), fill = "Taxon"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2, byrow = TRUE))
}
