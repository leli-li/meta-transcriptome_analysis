# Check required packages
needed_pkgs <- c("polyester", "Biostrings")
missing <- needed_pkgs[!sapply(needed_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing) > 0) {
  stop(
    "Missing R package(s): ", paste(missing, collapse = ", "), "\n",
    "Please install the missing package(s) by running:\n",
    "  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
    "  BiocManager::install(c('", paste(missing, collapse = "', '"), "'))\n",
    call. = FALSE
  )
}

# Load packages 
library(polyester)
library(Biostrings)



# === 0. Setup ===
# Set working directory
setwd("path/to/working directory")

# Set input paths
potato_fasta   <- "Data/random_potato_100.fna"
microbe_fasta  <- "Data/microbes_all.fna"
# Create output directory
outdir <- "polyester_out"
if (dir.exists(outdir) && length(list.files(outdir)) > 0L) {
  stop("Output directory is not empty: ", outdir, call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Set the random seed for reproducibility.
seed <- 20261326



# === 1. Read transcript sequences ===
potato_seqs   <- readDNAStringSet(potato_fasta)
microbe_seqs  <- readDNAStringSet(microbe_fasta)
all_seqs      <- c(potato_seqs, microbe_seqs)
tx_names      <- names(all_seqs)
N             <- length(tx_names)



# === 2. Identify host (potato) vs microbial transcripts ===
is_potato <- seq_along(all_seqs) <= length(potato_seqs)
stopifnot(!anyDuplicated(tx_names))
stopifnot(all(width(potato_seqs) >= 150))



# === 3. Define the simulation design ===
sample_names  <- c("profile1_contam0", "profile1_contam0.5", "profile1_contam5",
                   "profile2_contam0", "profile2_contam0.5", "profile2_contam5") 
profile_ids   <- rep(c("profile1", "profile2"), each = 3)
contam_pct    <- rep(c(0, 0.005, 0.05), times = 2)
samples       <- length(sample_names)
total_reads   <- 5e5 # each sample will contain 5e5 read pairs

metadata <- data.frame(
  sample          = sample_names,
  polyester_id    = sprintf("sample_%02d", seq_len(samples)),
  profile         = profile_ids,
  host_pct        = contam_pct * 100,
  expected_reads  = total_reads,
  read1_file      = paste0(sample_names, "_1.fasta"),
  read2_file      = paste0(sample_names, "_2.fasta"),
  potato_reads    = NA_real_,
  microbe_reads   = NA_real_,
  actual_reads    = NA_real_,
  actual_host_pct = NA_real_,
  stringsAsFactors = FALSE
) 

# Initialize the transcript-by-sample count matrix.
expr <- matrix(0, nrow = N, ncol = samples,
               dimnames = list(tx_names, sample_names))

set.seed(seed)

# Keep the relative weights of host transcripts constant across all mock samples
potato_profile <- sample(1:10, length(potato_seqs), replace = TRUE)
# Host-transcript relative weights are sampled from 1 to 10.

# Generate two microbial expression profiles.
microbe_profiles <- list(
  profile1 = sample(10:100, length(microbe_seqs), replace = TRUE),
  profile2 = sample(10:100, length(microbe_seqs), replace = TRUE)
)
# 10:100 means the range of relative weights of microbial transcript. It can be adjusted. 

# Calculate the proportion of transcripts for each simulated sample
for (i in seq_len(samples)) {
  pct_p <- contam_pct[i]
  reads_potato  <- round(total_reads * pct_p)
  reads_microbe <- total_reads - reads_potato
  
  potato_expr  <- potato_profile
  microbe_expr <- microbe_profiles[[profile_ids[i]]]
  
  potato_scaled  <- round(potato_expr / sum(potato_expr) * reads_potato)
  microbe_scaled <- round(microbe_expr / sum(microbe_expr) * reads_microbe)
  
  expr[is_potato, i]  <- potato_scaled
  expr[!is_potato, i] <- microbe_scaled
  
  metadata$potato_reads[i]     <- sum(potato_scaled)
  metadata$microbe_reads[i]    <- sum(microbe_scaled)
  metadata$actual_reads[i]     <- sum(expr[, i])
  metadata$actual_host_pct[i]  <- metadata$potato_reads[i] / metadata$actual_reads[i] * 100
}



# === 4. Pre-simulation validation ===
stopifnot(nrow(expr) == N)
stopifnot(ncol(expr) == samples)
stopifnot(all(dimnames(expr)[[2]] == sample_names))
stopifnot(all(expr >= 0))
stopifnot(all(profile_ids %in% c("profile1", "profile2")))
stopifnot(all(abs(metadata$actual_reads - total_reads) <= 0.01 * total_reads))
stopifnot(all(expr[is_potato, contam_pct == 0, drop = FALSE] == 0))
stopifnot(all(
  abs(metadata$actual_host_pct - metadata$host_pct) <= 0.01
))

message("Pre-simulation checks passed!")



# === 5. Simulation ===
simulate_experiment_countmat(
  fasta       = c(potato_fasta, microbe_fasta),
  readmat     = expr,
  outdir      = outdir,
  paired      = TRUE,
  seed        = seed,
  readlen     = 150,
  fraglen     = 300,
  fragsd      = 50,
  error_model = "uniform",
  error_rate  = 0.005
)


# === 6. Rename Polyester output files and save simulation truth ===
# Column order defines the sample order in Polyester output.
expected_order <- colnames(expr)

# Locate the generated mate files.
files_1 <- list.files(
  outdir,
  pattern = "^sample_\\d+_1\\.fasta$",
  full.names = TRUE
)
files_2 <- list.files(
  outdir,
  pattern = "^sample_\\d+_2\\.fasta$",
  full.names = TRUE
)

# Extract sample numbers and sort files numerically.
get_num <- function(x) {
  as.integer(sub(".*sample_(\\d+)_.*", "\\1", basename(x)))
}

files_1 <- files_1[order(get_num(files_1))]
files_2 <- files_2[order(get_num(files_2))]

# Confirm that Polyester generated the expected number of mate files.
stopifnot(length(files_1) == length(expected_order))
stopifnot(length(files_2) == length(expected_order))

# Rename files and confirm that every operation succeeded.
renamed_1 <- file.rename(
  files_1,
  file.path(outdir, paste0(expected_order, "_1.fasta"))
)
renamed_2 <- file.rename(
  files_2,
  file.path(outdir, paste0(expected_order, "_2.fasta"))
)
stopifnot(all(renamed_1), all(renamed_2))


# Save
write.csv(metadata, file.path(outdir, "metadata.csv"), row.names = FALSE)
write.csv(expr, file.path(outdir, "truth_counts.csv"), row.names = TRUE)
saveRDS(expr, file.path(outdir, "truth_counts.rds"))

message("Simulation completed; output files were renamed to match sample IDs.\n Next step: follow the validation commands in ",
        "03_Mock_RNA_check.md.")

# NEXT-STEP: Check the simulated samples
## Find the 03_Mock_RNA_check.md and follow the commands in it.