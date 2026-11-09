# Run PCA and PERMANOVA for controlled-climate potato.
# Run from Downstream_analysis/ with Rscript.
# Put metadata.tsv, kk2/, dm/, and overlap/ in dataset_dir.
# Packages: patchwork, vegan, ggplot2, ggrepel, ggalt, ggforce.
library(patchwork)

# ============================================================
# 1. Load shared functions
# ============================================================
source("R_scripts/PCA_functions.R")

# ============================================================
# 2. Dataset-specific settings
# ============================================================
dataset_name <- "controlled-climate-potato"
dataset_dir <- file.path("Datasets", dataset_name)

metadata_file <- file.path(
  dataset_dir,
  "metadata.tsv"
)

# PERMANOVA model
permanova_rhs <- "genotype"
permanova_by <- NULL

# Extract sample IDs from Kraken2 and DIAMOND-MEGAN6 headers.
sample_id_transform <- function(x) {
  sub(
    "(_kk2_nr[.]report|_Unmapped[.]out[.]trinity)$",
    "",
    x
  )
}

# PCA plot aesthetics
color_col <- "genotype"
shape_col <- NULL
enclosure_method <- "ellipse"

# Use default ggplot2 colors
group_colors <- NULL

# ============================================================
# 3. Load sample metadata
# ============================================================

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

metadata <- read.delim(
  metadata_file,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 4. Define analysis sets and construct input paths
# ============================================================

analysis_plan <- data.frame(
  file_id = c(
    "bac_f",
    "bac_g",
    "bac_sp",
    "fungi_f",
    "fungi_g",
    "fungi_sp"
  ),
  organism = c(
    "Bacteria",
    "Bacteria",
    "Bacteria",
    "Fungi",
    "Fungi",
    "Fungi"
  ),
  taxonomic_level = c(
    "Family",
    "Genus",
    "Species",
    "Family",
    "Genus",
    "Species"
  ),
  stringsAsFactors = FALSE
)

pipelines <- c(
  kk2 = "Kraken2",
  dm = "MEGAN6"
)

# Expand the six taxonomic analysis sets across both pipelines
n_analysis_sets <- nrow(analysis_plan)

analysis_plan <- analysis_plan[
  rep(seq_len(n_analysis_sets), times = length(pipelines)),
  ,
  drop = FALSE
]

analysis_plan$pipeline_dir <- rep(
  names(pipelines),
  each = n_analysis_sets
)

analysis_plan$pipeline_name <- unname(
  pipelines[analysis_plan$pipeline_dir]
)

analysis_plan$count_file <- file.path(
  dataset_dir,
  analysis_plan$pipeline_dir,
  paste0(analysis_plan$file_id, ".txt")
)

analysis_plan$overlap_file <- file.path(
  dataset_dir,
  "overlap",
  paste0(analysis_plan$file_id, ".txt")
)

analysis_plan$result_name <- paste(
  analysis_plan$pipeline_dir,
  analysis_plan$file_id,
  sep = "_"
)

rownames(analysis_plan) <- NULL

analysis_plan$count_file_exists <- file.exists(
  analysis_plan$count_file
)

analysis_plan$overlap_file_exists <- file.exists(
  analysis_plan$overlap_file
)

cat("\nConstructed analysis plan:\n\n")

print(
  analysis_plan[
    ,
    c(
      "result_name",
      "pipeline_name",
      "organism",
      "taxonomic_level",
      "count_file",
      "count_file_exists",
      "overlap_file",
      "overlap_file_exists"
    )
  ],
  row.names = FALSE
)

# Check all input files before starting the analyses
required_files <- unique(
  c(
    analysis_plan$count_file,
    analysis_plan$overlap_file
  )
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Required input files were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}


# ============================================================
# 5. Run PCA and PERMANOVA analyses and create plots
# ============================================================

results <- list()
plots <- list()

for (i in seq_len(nrow(analysis_plan))) {
  
  current_analysis <- analysis_plan[i, ]
  result_name <- current_analysis$result_name
  
  results[[result_name]] <- run_PCA_analysis(
    count_file = current_analysis$count_file,
    overlap_file = current_analysis$overlap_file,
    metadata = metadata,
    permanova_rhs = permanova_rhs,
    permanova_by = permanova_by,
    count_header = TRUE,
    sample_id_transform = sample_id_transform
  )
  
  plots[[result_name]] <- plot_PCA_result(
    analysis_result = results[[result_name]],
    color_col = color_col,
    pipeline_name = current_analysis$pipeline_name,
    organism = current_analysis$organism,
    taxonomic_level = current_analysis$taxonomic_level,
    group_colors = group_colors,
    shape_col = shape_col,
    enclosure_method = enclosure_method
  )
  
  # Add arrows for the top 10 most abundant overlapping taxa
  # to species-level PCA plots
  if (current_analysis$taxonomic_level == "Species") {
    
    plots[[result_name]] <- add_taxa_arrows(
      p = plots[[result_name]],
      analysis_result = results[[result_name]],
      top_n = 10,
      arrow_color = "purple"
    )
  }
}

# ============================================================
# 6. Save analysis results and PCA plots
# ============================================================

output_dir <- file.path(dataset_dir, "PCA_results")

plot_dir <- file.path(
  output_dir,
  "PCA_plots"
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

results_file <- file.path(
  output_dir,
  "PCA_permanova_results.rds"
)

saveRDS(
  results,
  file = results_file
)

for (plot_name in names(plots)) {
  
  plot_file <- file.path(
    plot_dir,
    paste0(plot_name, ".pdf")
  )
  
  ggsave(
    filename = plot_file,
    plot = plots[[plot_name]],
    width = 6,
    height = 5
  )
}

message("Analysis results saved to: ", results_file)
message(length(plots), " PCA plots saved to: ", plot_dir)


# Save combined PCA panels (3 plots per PDF)
combined_plot_dir <- file.path(
  output_dir,
  "PCA_panels"
)

dir.create(
  combined_plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

combined_plots <- list(
  kk2_bacteria = (plots$kk2_bac_f + plots$kk2_bac_g + plots$kk2_bac_sp) +
    plot_layout(guides = "collect"),
  
  kk2_fungi = (plots$kk2_fungi_f + plots$kk2_fungi_g + plots$kk2_fungi_sp) +
    plot_layout(guides = "collect"),
  
  dm_bacteria = (plots$dm_bac_f + plots$dm_bac_g + plots$dm_bac_sp) +
    plot_layout(guides = "collect"),
  
  dm_fungi = (plots$dm_fungi_f + plots$dm_fungi_g + plots$dm_fungi_sp) +
    plot_layout(guides = "collect")
)

for (panel_name in names(combined_plots)) {
  
  panel_file <- file.path(
    combined_plot_dir,
    paste0(panel_name, ".pdf")
  )
  
  ggsave(
    filename = panel_file,
    plot = combined_plots[[panel_name]],
    width = 14.4,
    height = 4.71
  )
}

message(length(combined_plots), " combined PCA panels saved to: ", combined_plot_dir)


# ============================================================
# 7. Inspect results and plots
# ============================================================

results$kk2_bac_f$explained_variance[1:2]
results$kk2_bac_f$permanova
plots$kk2_bac_f

results$dm_bac_f$explained_variance[1:2]
results$dm_bac_f$permanova
plots$dm_bac_f
