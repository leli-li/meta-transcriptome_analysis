# Calculate Shannon diversity for shared bacterial and fungal taxa.
# Combine their counts within each taxonomic rank.
# Use counts directly, without a Hellinger transformation.

library(vegan)
library(ggplot2)

read_Shannon_count_matrix <- function(
    count_file,
    sample_ids,
    count_header = FALSE,
    sample_id_transform = identity
) {
  count_data <- read.delim(
    file = count_file,
    sep = "\t",
    header = count_header,
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = ""
  )

  rownames(count_data) <- gsub("_", " ", rownames(count_data))

  if (ncol(count_data) != length(sample_ids)) {
    stop(
      "The number of count-matrix columns does not match the number ",
      "of supplied sample IDs for: ",
      count_file
    )
  }

  if (count_header) {
    # Get sample IDs from the header.
    ids <- sample_id_transform(colnames(count_data))

    if (
      length(ids) != ncol(count_data) ||
      anyNA(ids) ||
      any(ids == "") ||
      anyDuplicated(ids) ||
      !setequal(ids, sample_ids)
    ) {
      stop(
        "Count-matrix sample IDs do not match sample_ids: ",
        count_file
      )
    }

    # Put columns in sample ID order.
    colnames(count_data) <- ids
    count_data <- count_data[, sample_ids, drop = FALSE]

  } else {
    # Without headers, sample_ids must match the column order.
    colnames(count_data) <- sample_ids
  }

  if (anyNA(count_data)) {
    stop("Missing values were found in the count matrix: ", count_file)
  }

  if (any(as.matrix(count_data) < 0)) {
    stop("Negative values were found in the count matrix: ", count_file)
  }

  return(count_data)
}


retain_Shannon_overlap <- function(
    count_data,
    overlap_file
) {
  overlap <- read.delim(
    overlap_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )[[1]]

  overlap <- gsub("_", " ", overlap)

  missing_taxa <- setdiff(
    overlap,
    rownames(count_data)
  )

  if (length(missing_taxa) > 0) {
    stop(
      "Some overlapping taxa were not found in the count matrix: ",
      paste(missing_taxa, collapse = ", ")
    )
  }

  count_data <- count_data[
    overlap,
    ,
    drop = FALSE
  ]

  return(count_data)
}


run_Shannon_analysis <- function(
    bac_count_file,
    fungi_count_file,
    bac_overlap_file,
    fungi_overlap_file,
    sample_ids,
    metadata,
    count_header = FALSE,
    sample_id_transform = identity
) {

  if (!"sample_id" %in% colnames(metadata)) {
    stop("metadata must contain a 'sample_id' column.")
  }

  if (anyDuplicated(sample_ids)) {
    stop("sample_ids contains duplicated sample IDs.")
  }

  if (length(sample_ids) != nrow(metadata)) {
    stop(
      "The number of supplied sample IDs does not match ",
      "the number of metadata rows."
    )
  }

  missing_metadata <- setdiff(
    sample_ids,
    metadata$sample_id
  )

  extra_metadata <- setdiff(
    metadata$sample_id,
    sample_ids
  )

  if (length(missing_metadata) > 0) {
    stop(
      "Some supplied sample IDs were not found in metadata: ",
      paste(missing_metadata, collapse = ", ")
    )
  }

  if (length(extra_metadata) > 0) {
    stop(
      "Some metadata sample IDs were not found in supplied sample_ids: ",
      paste(extra_metadata, collapse = ", ")
    )
  }

  bac_counts <- read_Shannon_count_matrix(
    count_file = bac_count_file,
    sample_ids = sample_ids,
    count_header = count_header,
    sample_id_transform = sample_id_transform
  )

  fungi_counts <- read_Shannon_count_matrix(
    count_file = fungi_count_file,
    sample_ids = sample_ids,
    count_header = count_header,
    sample_id_transform = sample_id_transform
  )

  bac_overlap_counts <- retain_Shannon_overlap(
    count_data = bac_counts,
    overlap_file = bac_overlap_file
  )

  fungi_overlap_counts <- retain_Shannon_overlap(
    count_data = fungi_counts,
    overlap_file = fungi_overlap_file
  )

  combined_counts <- rbind(
    bac_overlap_counts,
    fungi_overlap_counts
  )

  combined_counts <- t(combined_counts)

  if (any(rowSums(combined_counts) == 0)) {
    zero_samples <- rownames(combined_counts)[
      rowSums(combined_counts) == 0
    ]

    stop(
      "No retained reads were found for sample(s): ",
      paste(zero_samples, collapse = ", ")
    )
  }

  # diversity() converts counts to proportions within each sample.
  shannon <- diversity(
    combined_counts,
    index = "shannon",
    MARGIN = 1,
    base = exp(1)
  )

  metadata_matched <- metadata[
    match(
      rownames(combined_counts),
      metadata$sample_id
    ),
    ,
    drop = FALSE
  ]

  if (anyNA(metadata_matched$sample_id)) {
    stop("Some samples could not be matched to metadata.")
  }

  if (!identical(
    as.character(metadata_matched$sample_id),
    rownames(combined_counts)
  )) {
    stop(
      "Metadata order could not be matched exactly to the count matrix."
    )
  }

  shannon_data <- cbind(
    data.frame(
      sample_id = rownames(combined_counts),
      Shannon = as.numeric(shannon),
      stringsAsFactors = FALSE
    ),
    metadata_matched[
      ,
      setdiff(colnames(metadata_matched), "sample_id"),
      drop = FALSE
    ]
  )

  return(
    list(
      count_matrix = combined_counts,
      shannon = shannon,
      shannon_data = shannon_data,
      bacterial_taxa = rownames(bac_overlap_counts),
      fungal_taxa = rownames(fungi_overlap_counts)
    )
  )
}


plot_Shannon_result <- function(
    analysis_result,
    group_col,
    pipeline_name,
    taxonomic_level,
    group_levels = NULL,
    group_labels = NULL,
    y_limits = NULL
) {

  plot_data <- analysis_result$shannon_data

  if (!group_col %in% colnames(plot_data)) {
    stop(
      "Grouping column not found in Shannon data: ",
      group_col
    )
  }

  plot_data$plot_group <- plot_data[[group_col]]

  if (!is.null(group_levels)) {
    missing_groups <- setdiff(
      unique(as.character(plot_data$plot_group)),
      group_levels
    )

    if (length(missing_groups) > 0) {
      stop(
        "group_levels is missing group(s): ",
        paste(missing_groups, collapse = ", ")
      )
    }

    plot_data$plot_group <- factor(
      plot_data$plot_group,
      levels = group_levels
    )
  }

  p <- ggplot(
    plot_data,
    aes(
      x = plot_group,
      y = Shannon
    )
  ) +
    geom_boxplot() +
    labs(
      x = NULL,
      y = "Shannon index",
      title = pipeline_name,
      subtitle = paste0(
        taxonomic_level,
        " level"
      )
    ) +
    theme(
      panel.grid = element_line(
        color = "gray",
        linetype = 2,
        linewidth = 0.1
      ),
      panel.background = element_rect(
        color = "black",
        fill = "transparent"
      ),
      legend.key = element_rect(
        fill = "transparent"
      )
    )

  if (!is.null(group_labels)) {
    p <- p +
      scale_x_discrete(
        labels = group_labels
      )
  }

  if (!is.null(y_limits)) {
    if (
      !is.numeric(y_limits) ||
      length(y_limits) != 2 ||
      anyNA(y_limits)
    ) {
      stop("y_limits must be a numeric vector of length 2.")
    }

    p <- p +
      coord_cartesian(
        ylim = y_limits
      )
  }

  return(p)
}
