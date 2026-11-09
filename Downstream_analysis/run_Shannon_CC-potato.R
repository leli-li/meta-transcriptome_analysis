# Run Shannon diversity for controlled-climate potato.
# Run from Downstream_analysis/ with Rscript.
# Put metadata.tsv, kk2/, dm/, and overlap/ in dataset_dir.
# Packages: patchwork, vegan, ggplot2.
# Combine shared bacterial and fungal counts within each taxonomic rank.

# ============================================================
# 1. Load shared analysis functions
# ============================================================

source("R_scripts/Shannon_functions.R")
library(patchwork)

# ============================================================
# 2. Dataset settings
# ============================================================

dataset_name <- "controlled-climate-potato"

dataset_dir <- file.path("Datasets", dataset_name)

metadata_file <- file.path(
  dataset_dir,
  "metadata.tsv"
)

group_col <- "genotype"

group_levels <- c(
  "Bintje",
  "Desiree",
  "KingEdward"
)

group_labels <- group_levels

# ============================================================
# 3. Load metadata
# ============================================================

if (!file.exists(metadata_file)) {
  stop(
    "Metadata file not found: ",
    metadata_file
  )
}

metadata <- read.delim(
  metadata_file,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"sample_id" %in% colnames(metadata)) {
  stop("metadata must contain a 'sample_id' column.")
}

if (!group_col %in% colnames(metadata)) {
  stop(
    "Grouping column not found in metadata: ",
    group_col
  )
}

# ============================================================
# 4. Match samples by ID
# ============================================================

# Match matrix columns to the sample IDs in metadata.
sample_ids <- as.character(metadata$sample_id)

# Remove pipeline-specific suffixes from count-matrix column names.
sample_id_transform <- function(x) {
  sub(
    "(_kk2_nr[.]report|_Unmapped[.]out[.]trinity)$",
    "",
    x
  )
}


# ============================================================
# 5. Define analysis sets and construct input paths
# ============================================================

analysis_plan <- data.frame(
  taxon_id = c(
    "f",
    "g",
    "sp"
  ),
  taxonomic_level = c(
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

n_analysis_sets <- nrow(analysis_plan)

analysis_plan <- analysis_plan[
  rep(
    seq_len(n_analysis_sets),
    times = length(pipelines)
  ),
  ,
  drop = FALSE
]

analysis_plan$pipeline_dir <- rep(
  names(pipelines),
  each = n_analysis_sets
)

analysis_plan$pipeline_name <- unname(
  pipelines[
    analysis_plan$pipeline_dir
  ]
)

analysis_plan$bac_count_file <- file.path(
  dataset_dir,
  analysis_plan$pipeline_dir,
  paste0(
    "bac_",
    analysis_plan$taxon_id,
    ".txt"
  )
)

analysis_plan$fungi_count_file <- file.path(
  dataset_dir,
  analysis_plan$pipeline_dir,
  paste0(
    "fungi_",
    analysis_plan$taxon_id,
    ".txt"
  )
)

analysis_plan$bac_overlap_file <- file.path(
  dataset_dir,
  "overlap",
  paste0(
    "bac_",
    analysis_plan$taxon_id,
    ".txt"
  )
)

analysis_plan$fungi_overlap_file <- file.path(
  dataset_dir,
  "overlap",
  paste0(
    "fungi_",
    analysis_plan$taxon_id,
    ".txt"
  )
)

analysis_plan$result_name <- paste(
  analysis_plan$pipeline_dir,
  analysis_plan$taxon_id,
  sep = "_"
)

rownames(analysis_plan) <- NULL

cat(
  "\nConstructed all-potato Shannon analysis plan:\n\n"
)

print(
  analysis_plan[
    ,
    c(
      "result_name",
      "pipeline_name",
      "taxonomic_level",
      "bac_count_file",
      "fungi_count_file",
      "bac_overlap_file",
      "fungi_overlap_file"
    )
  ],
  row.names = FALSE
)

required_files <- unique(
  c(
    analysis_plan$bac_count_file,
    analysis_plan$fungi_count_file,
    analysis_plan$bac_overlap_file,
    analysis_plan$fungi_overlap_file
  )
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Required input files were not found:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

# ============================================================
# 6. Run Shannon analyses
# ============================================================

results <- list()

for (i in seq_len(nrow(analysis_plan))) {
  
  current_analysis <- analysis_plan[
    i,
    ,
    drop = FALSE
  ]
  
  result_name <- current_analysis$result_name
  
  results[[result_name]] <- run_Shannon_analysis(
    bac_count_file = current_analysis$bac_count_file,
    fungi_count_file = current_analysis$fungi_count_file,
    bac_overlap_file = current_analysis$bac_overlap_file,
    fungi_overlap_file = current_analysis$fungi_overlap_file,
    sample_ids = sample_ids,
    metadata = metadata,
    count_header = TRUE,
    sample_id_transform = sample_id_transform
  )
}

# Summarize per-sample Shannon indices in the same group order as the plots.
# n is the number of samples in each group, not the number of reads or taxa.
shannon_summary_rows <- list()

for (i in seq_len(nrow(analysis_plan))) {

  current_analysis <- analysis_plan[i, , drop = FALSE]
  shannon_data <- results[[current_analysis$result_name]]$shannon_data

  for (j in seq_along(group_levels)) {

    current_group <- group_levels[j]
    values <- shannon_data$Shannon[
      shannon_data[[group_col]] == current_group
    ]

    shannon_summary_rows[[length(shannon_summary_rows) + 1]] <- data.frame(
      pipeline = current_analysis$pipeline_name,
      taxonomic_level = current_analysis$taxonomic_level,
      group = gsub("\n", " - ", group_labels[j]),
      n = length(values),
      median_Shannon = median(values),
      stringsAsFactors = FALSE
    )
  }
}

shannon_summary <- do.call(rbind, shannon_summary_rows)
rownames(shannon_summary) <- NULL

# Round only the console display; keep full precision in the exported table.
cat("\nShannon group medians:\n\n")
print(
  transform(shannon_summary, median_Shannon = round(median_Shannon, 3)),
  row.names = FALSE
)

# ============================================================
# 7. Determine matched y-axis ranges for Kraken2 and MEGAN6
# ============================================================

get_Shannon_limits <- function(
    results,
    result_names
) {

  values <- unlist(
    lapply(
      result_names,
      function(result_name) {
        results[[result_name]]$shannon_data$Shannon
      }
    ),
    use.names = FALSE
  )

  range(
    values,
    na.rm = TRUE
  )
}

shannon_limits <- list(
  Family = get_Shannon_limits(
    results,
    c("kk2_f", "dm_f")
  ),
  Genus = get_Shannon_limits(
    results,
    c("kk2_g", "dm_g")
  ),
  Species = get_Shannon_limits(
    results,
    c("kk2_sp", "dm_sp")
  )
)


# ============================================================
# 8. Create Shannon plots
# ============================================================

plots <- list()

for (i in seq_len(nrow(analysis_plan))) {

  current_analysis <- analysis_plan[
    i,
    ,
    drop = FALSE
  ]

  result_name <- current_analysis$result_name
  taxonomic_level <- current_analysis$taxonomic_level

  plots[[result_name]] <- plot_Shannon_result(
    analysis_result = results[[result_name]],
    group_col = group_col,
    pipeline_name = current_analysis$pipeline_name,
    taxonomic_level = taxonomic_level,
    group_levels = group_levels,
    group_labels = group_labels,
    y_limits = shannon_limits[[taxonomic_level]]
  )
}


# ============================================================
# 9. Save analysis results
# ============================================================

output_dir <- file.path(
  dataset_dir,
  "Shannon_results"
)

plot_dir <- file.path(
  output_dir,
  "Shannon_plots"
)

panel_dir <- file.path(
  output_dir,
  "Shannon_panels"
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  panel_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

results_file <- file.path(
  output_dir,
  "Shannon_results.rds"
)

saveRDS(
  results,
  file = results_file
)

message(
  "Shannon results saved to: ",
  results_file
)

# Export the same group summary shown in the console.
summary_file <- file.path(output_dir, "Shannon_group_summary.tsv")

write.table(
  shannon_summary,
  file = summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Shannon group summary saved to: ", summary_file)


# ============================================================
# 10. Save individual Shannon plots
# ============================================================

for (plot_name in names(plots)) {

  plot_file <- file.path(
    plot_dir,
    paste0(
      plot_name,
      "_Shannon.pdf"
    )
  )

  ggsave(
    filename = plot_file,
    plot = plots[[plot_name]],
    width = 6,
    height = 5
  )
}

message(
  length(plots),
  " Shannon plots saved to: ",
  plot_dir
)


# ============================================================
# 11. Save combined Shannon panel
# ============================================================

# Rows = Family / Genus / Species
# Columns = Kraken2 / MEGAN6

CC_potato_panel <- (
  plots$kk2_f +
    plots$dm_f
) / (
  plots$kk2_g +
    plots$dm_g
) / (
  plots$kk2_sp +
    plots$dm_sp
)

panel_file <- file.path(
  panel_dir,
  paste0(dataset_name, "_Shannon.pdf")
)

ggsave(
  filename = panel_file,
  plot = CC_potato_panel,
  width = 12.4,
  height = 16
)

message(
  "Combined all-potato Shannon panel saved to: ",
  panel_file
)


# ============================================================
# 12. Inspect results and plots
# ============================================================

results$kk2_f$shannon_data
results$dm_f$shannon_data

plots$kk2_f
plots$dm_f

CC_potato_panel
