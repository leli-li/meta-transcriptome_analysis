# Rank and plot top taxa for controlled-climate potato.
# Run from Downstream_analysis/ with Rscript.
# Put metadata.tsv, kk2/, dm/, and overlap/ in dataset_dir.
# Packages: ggplot2, patchwork.
# Reruns replace existing output tables and PDFs.

# ============================================================
# 1. Load shared functions and set paths
# ============================================================

source("R_scripts/Top_taxa_functions.R")
library(ggplot2)

dataset_dir <- file.path("Datasets", "controlled-climate-potato")
output_dir <- file.path(dataset_dir, "Top_taxa_results")
table_dir <- file.path(output_dir, "tables")
plot_dir <- file.path(output_dir, "plots")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

top_n <- 10

# ============================================================
# 2. Read metadata and define the groups to summarize
# ============================================================

metadata <- read.delim(
  file.path(dataset_dir, "metadata.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE,
  check.names = FALSE
)

# Select the Top 10 independently for each genotype (three samples per group).
group_levels <- c("Bintje", "Desiree", "KingEdward")
group_labels <- c(Bintje = "Bintje", Desiree = "Desiree", KingEdward = "King Edward")

if (!"genotype" %in% names(metadata) ||
    anyNA(metadata$genotype) ||
    any(!metadata$genotype %in% group_levels)) {
  stop("Metadata must contain valid genotype groups.")
}

# Set the plot order: Bintje, Desiree, then King Edward.
metadata <- metadata[
  order(match(metadata$genotype, group_levels)), , drop = FALSE
]

# Remove workflow suffixes from the column names to get sample IDs.
sample_id_transform <- function(x) {
  sub("(_kk2_nr[.]report|_Unmapped[.]out[.]trinity)$", "", x)
}

# ============================================================
# 3. List the input tables and workflow labels
# ============================================================

# The six tables cover bacteria and fungi at the family, genus, and species levels.
# The kk2, dm, and overlap directories use the same filenames.
file_ids <- c("bac_f", "bac_g", "bac_sp", "fungi_f", "fungi_g", "fungi_sp")
table_titles <- c(
  bac_f = "Bacteria - Family level",
  bac_g = "Bacteria - Genus level",
  bac_sp = "Bacteria - Species level",
  fungi_f = "Fungi - Family level",
  fungi_g = "Fungi - Genus level",
  fungi_sp = "Fungi - Species level"
)

pipelines <- c("kk2", "dm")
pipeline_titles <- c(kk2 = "Kraken2", dm = "DIAMOND-MEGAN6")
count_units <- c(kk2 = "read-pair counts", dm = "contig counts")

# ============================================================
# 4. Rank taxa and save tables
# ============================================================

# Store 12 results: 2 workflows x 2 kingdoms x 3 taxonomic ranks.
# Each result contains all three genotypes, identified by the group column,
# with complete rankings, Top 10 tables, and Top 10 + Others data for plotting.
results <- list()

for (pipeline in pipelines) {
  for (file_id in file_ids) {
    result_name <- paste(pipeline, file_id, sep = "_")
    message("Summarizing: ", result_name)

    # Rank by total count in the group, not mean relative abundance.
    # The overlap list marks shared taxa without filtering them.
    result <- summarize_top_taxa(
      count_file = file.path(dataset_dir, pipeline, paste0(file_id, ".txt")),
      metadata = metadata,
      group_col = "genotype",
      overlap_file = file.path(dataset_dir, "overlap", paste0(file_id, ".txt")),
      sample_id_transform = sample_id_transform,
      top_n = top_n
    )
    results[[result_name]] <- result

    # ranked_counts: counts, percentages, shared-taxon flags, and ranks for all taxa.
    # top_taxa: up to 10 taxa with positive counts per genotype.
    # plot_data: Top 10 + Others. Percentages use the total counts of all taxa
    # in this group, workflow, kingdom, and rank, rather than only the Top 10.
    write.table(
      result$ranked_counts,
      file.path(table_dir, paste0(result_name, "_ranked_counts.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    write.table(
      result$top_taxa,
      file.path(table_dir, paste0(result_name, "_top", top_n, ".tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    write.table(
      result$plot_data,
      file.path(table_dir, paste0(result_name, "_plot_data.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
}

# ============================================================
# 5. Plot Top 10 + Others and save individual PDFs
# ============================================================

plots <- list()

for (file_id in file_ids) {
  # Use the same taxon colors across workflows for a given kingdom and rank.
  # Combine taxon names only to assign colors; each workflow retains its own
  # counts and independently selected Top 10.
  displayed_taxa <- union(
    results[[paste0("kk2_", file_id)]]$top_taxa$taxon,
    results[[paste0("dm_", file_id)]]$top_taxa$taxon
  )
  taxon_colors <- setNames(
    grDevices::hcl.colors(length(displayed_taxa), palette = "Dark 3"),
    displayed_taxa
  )

  for (pipeline in pipelines) {
    result_name <- paste(pipeline, file_id, sep = "_")
    plots[[result_name]] <- plot_top_taxa(
      result = results[[result_name]],
      title = paste(pipeline_titles[[pipeline]], table_titles[[file_id]], sep = " | "),
      count_unit = count_units[[pipeline]],
      taxon_colors = taxon_colors
    ) +
      scale_x_discrete(limits = group_levels, labels = group_labels) +
      labs(x = "Genotype")

    # Each bar has its own Top 10; Others can contain different taxa across groups.
    # Allow extra legend space when genotypes have different Top 10 taxa.
    n_legend_items <- length(unique(results[[result_name]]$plot_data$taxon))
    plot_height <- max(8, 6 + 0.25 * ceiling(n_legend_items / 2))

    # Save each plot separately and specify plot explicitly to save the intended figure.
    # The workflows use different counting units; expressing counts as percentages
    # does not remove this methodological difference.
    ggsave(
      filename = file.path(plot_dir, paste0(result_name, "_top", top_n, ".pdf")),
      plot = plots[[result_name]],
      width = 10, height = plot_height, units = "in"
    )
  }
}

# ============================================================
# 6. Inspect a result interactively (optional)
# ============================================================

# After running the script, uncomment any line below to inspect bacterial
# species-level results or display the corresponding plot.
# results$kk2_bac_sp$top_taxa
# results$dm_bac_sp$top_taxa
# print(plots$kk2_bac_sp)
# print(plots$dm_bac_sp)

# ============================================================
# 7. Combine ranks and save four panel PDFs
# ============================================================

# If plots and output_dir are still available from this dataset's run,
# this section can be run on its own without repeating the analysis.
# patchwork combines the existing ggplot objects; no PDF inputs are needed.
panel_dir <- file.path(output_dir, "panels")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

# Each panel has three columns: Family, Genus, and Species.
# Keep separate legends because the ranks contain different taxa.
panel_sets <- list(
  kk2_bacteria = c("kk2_bac_f", "kk2_bac_g", "kk2_bac_sp"),
  kk2_fungi = c("kk2_fungi_f", "kk2_fungi_g", "kk2_fungi_sp"),
  dm_bacteria = c("dm_bac_f", "dm_bac_g", "dm_bac_sp"),
  dm_fungi = c("dm_fungi_f", "dm_fungi_g", "dm_fungi_sp")
)

panels <- list()
panel_width <- 30  # Inches; allow space for all three plots and their legends.

for (panel_name in names(panel_sets)) {
  panel_plots <- plots[panel_sets[[panel_name]]]
  panels[[panel_name]] <- patchwork::wrap_plots(
    panel_plots, ncol = 3, guides = "keep"
  )

  # Use the longest legend to set the height of the combined figure.
  # Each plot uses two legend columns; increase the height for more taxa.
  legend_sizes <- vapply(
    panel_plots,
    function(p) length(unique(p$data$taxon)),
    integer(1)
  )
  panel_height <- max(8, 6 + 0.25 * ceiling(max(legend_sizes) / 2))

  ggplot2::ggsave(
    filename = file.path(panel_dir, paste0(panel_name, ".pdf")),
    plot = panels[[panel_name]],
    width = panel_width, height = panel_height, units = "in"
  )
}

# Optional: display a combined figure in the RStudio plot pane.
# print(panels$kk2_bacteria)

message("Done. Tables, individual plots, and panels saved in: ", output_dir)
