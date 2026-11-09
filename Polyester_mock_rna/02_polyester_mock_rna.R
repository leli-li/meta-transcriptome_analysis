# Simulate paired-end RNA-seq reads from potato, 22 bacterial and 8 fungal species.
# Profiles vary transcript weights within species, not species-level proportions.
# Run from the project directory containing data/ and fungi/.

# === 0. Packages, paths and simulation settings ===
library(polyester)
library(Biostrings)

potato_fasta <- "data/random_potato_100.fna"
bacteria_dir <- "data"
fungi_dir <- "fungi/fungal_cds/ncbi_dataset/data"
outdir <- "polyester_out_bacteria_fungi"

seed <- 20261326
total_reads <- 5e5                  # Read PAIRS per sample, not individual mates.
readlen <- 150
fraglen <- 300
fragsd <- 50
error_rate <- 0.005
min_template_length <- readlen      # Exclude CDS shorter than a full read.

# Refuse to overwrite an existing run.
if (file.exists(outdir) && !dir.exists(outdir)) {
  stop("Output path is an existing file: ", outdir)
}
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE,
                                          no.. = TRUE)) > 0L) {
  stop("Output directory is not empty: ", outdir,
       ". Choose a different outdir for a new run.", call. = FALSE)
}


# === 1. Map source assemblies to species ===
# Assemblies with the same species TaxID share one read-pair budget.
sources <- read.table(text = "
kingdom assembly species_taxid species
bacteria GCA_000006785.2 1314 Streptococcus_pyogenes
bacteria GCA_000013285.1 1502 Clostridium_perfringens
bacteria GCA_000013785.1 316 Stutzerimonas_stutzeri
bacteria GCA_000015865.1 1515 Acetivibrio_thermocellus
bacteria GCA_000022345.1 562 Escherichia_coli
bacteria GCA_000023785.1 2724 Hirschia_baltica
bacteria GCA_000025905.1 395922 Coraliomargarita_akajimensis
bacteria GCA_000092125.1 52022 Allomeiothermus_silvanus
bacteria GCA_000092825.1 286802 Segniliparus_rotundus
bacteria GCA_000092985.1 2014 Nocardiopsis_dassonvillei
bacteria GCA_000143845.1 133926 Olsenella_uli
bacteria GCA_000143985.1 55206 Sediminispirochaeta_smaragdinae
bacteria GCA_000227705.3 377615 Thermobacillus_composti
bacteria GCA_000231385.3 79209 Desulfosporosinus_meridiei
bacteria GCA_000233715.3 102134 Desulfoscipio_gibsoniae
bacteria GCA_000235405.3 93466 Fervidobacterium_pennivorans
bacteria GCA_000242255.3 81475 Frateuria_aurantia
bacteria GCA_000255115.3 885581 Desulfosporosinus_acidiphilus
bacteria GCA_000265425.1 392734 Terriglobus_roseus
bacteria GCA_000325705.1 390884 Echinicola_vietnamensis
bacteria GCA_000439255.1 54736 Salmonella_bongori
bacteria GCA_000445015.1 1718 Corynebacterium_glutamicum
bacteria GCA_001692755.1 1515 Acetivibrio_thermocellus
bacteria GCA_001877055.1 2014 Nocardiopsis_dassonvillei
fungi GCF_023522795.1 169388 Fusarium_solani
fungi GCF_000149615.1 33178 Aspergillus_terreus
fungi GCF_033870435.1 5147 Sordaria_macrospora
fungi GCF_000150705.2 1048829 Paracoccidioides_lutzii
fungi GCF_003942205.1 105984 Apiotrichum_porosum
fungi GCF_000271585.1 5325 Trametes_versicolor
fungi GCF_001720195.1 5208 Cryptococcus_depauperatus
fungi GCF_000320585.1 984962 Heterobasidion_irregulare
", header = TRUE, stringsAsFactors = FALSE)
sources$species <- gsub("_", " ", sources$species, fixed = TRUE)

# Read bacterial CDS files separately and keep the NCBI directory layout for fungi.
sources$fasta <- ifelse(
  sources$kingdom == "bacteria",
  file.path(bacteria_dir, paste0(sources$assembly, ".cds.fna")),
  file.path(fungi_dir, sources$assembly, "cds_from_genomic.fna")
)
for (j in seq_len(nrow(sources))) {
  if (!file.exists(sources$fasta[j]) &&
      file.exists(paste0(sources$fasta[j], ".gz"))) {
    sources$fasta[j] <- paste0(sources$fasta[j], ".gz")
  }
}

host_source <- data.frame(
  kingdom = "host", assembly = "GCF_000226075.1", species_taxid = 4113L,
  species = "Solanum tuberosum", fasta = potato_fasta,
  stringsAsFactors = FALSE
)
sources <- rbind(host_source, sources)
stopifnot(!anyDuplicated(sources$assembly))
missing_files <- sources$fasta[!file.exists(sources$fasta)]
if (length(missing_files) > 0L) {
  stop("Input FASTA file(s) not found:\n", paste(missing_files, collapse = "\n"),
       "\nCheck the paths at the top of this script; do not merge inputs first.",
       call. = FALSE)
}


# === 2. Read templates and retain their source information ===
sequence_parts <- vector("list", nrow(sources))
mapping_parts <- vector("list", nrow(sources))
sources$input_sequences <- sources$retained_sequences <- sources$short_removed <- 0L
sources$input_md5 <- unname(tools::md5sum(sources$fasta))

for (j in seq_len(nrow(sources))) {
  seqs <- readDNAStringSet(sources$fasta[j])
  sources$input_sequences[j] <- length(seqs)
  keep <- width(seqs) >= min_template_length
  if (sources$kingdom[j] == "host" &&
      (length(seqs) != 100L || !all(keep))) {
    stop("Expected the original 100 potato transcripts, all >= readlen nt.")
  }
  sources$short_removed[j] <- sum(!keep)
  original_index <- which(keep)
  seqs <- seqs[keep]
  sources$retained_sequences[j] <- length(seqs)
  if (length(seqs) == 0L) {
    stop("No usable templates remain in: ", sources$fasta[j])
  }

  # Assign unique template IDs and retain the original headers for traceability.
  original_headers <- names(seqs)
  new_ids <- paste0(sources$kingdom[j], "_", sources$assembly[j],
                    "_tx", sprintf("%07d", original_index))
  names(seqs) <- new_ids
  sequence_parts[[j]] <- seqs
  mapping_parts[[j]] <- data.frame(
    transcript_id = new_ids,
    kingdom = sources$kingdom[j],
    assembly = sources$assembly[j],
    species_taxid = sources$species_taxid[j],
    species = sources$species[j],
    source_fasta = sources$fasta[j],
    original_header = original_headers,
    length_nt = width(seqs),
    stringsAsFactors = FALSE
  )
}

all_seqs <- do.call(c, sequence_parts)
transcript_sources <- do.call(rbind, mapping_parts)
rownames(transcript_sources) <- NULL
tx_names <- names(all_seqs)
N <- length(all_seqs)
stopifnot(!anyDuplicated(tx_names),
          identical(tx_names, transcript_sources$transcript_id))

is_potato <- transcript_sources$kingdom == "host"
is_bacteria <- transcript_sources$kingdom == "bacteria"
is_fungi <- transcript_sources$kingdom == "fungi"
microbial_species <- unique(transcript_sources$species_taxid[!is_potato])
stopifnot(length(unique(transcript_sources$species_taxid[is_bacteria])) == 22L,
          length(unique(transcript_sources$species_taxid[is_fungi])) == 8L,
          length(microbial_species) == 30L)

species_rows <- lapply(microbial_species, function(id) {
  which(!is_potato & transcript_sources$species_taxid == id)
})
names(species_rows) <- as.character(microbial_species)


# === 3. Allocate read pairs for the six samples ===
# Largest-remainder allocation preserves the exact read-pair total.
allocate_counts <- function(weights, total) {
  stopifnot(length(weights) > 0L, all(is.finite(weights)), all(weights >= 0),
            sum(weights) > 0, length(total) == 1L, is.finite(total),
            total >= 0, total == round(total))
  expected <- weights / sum(weights) * total
  counts <- floor(expected)
  remainder <- as.integer(total - sum(counts))
  if (remainder > 0L) {
    add_to <- order(expected - counts, decreasing = TRUE)[seq_len(remainder)]
    counts[add_to] <- counts[add_to] + 1
  }
  stopifnot(sum(counts) == total)
  counts
}

sample_names <- c("profile1_contam0", "profile1_contam0.5", "profile1_contam5",
                  "profile2_contam0", "profile2_contam0.5", "profile2_contam5")
profile_ids <- rep(c("profile1", "profile2"), each = 3)
contam_pct <- rep(c(0, 0.005, 0.05), times = 2)
samples <- length(sample_names)
metadata <- data.frame(
  sample = sample_names,
  polyester_id = sprintf("sample_%02d", seq_len(samples)),
  profile = profile_ids,
  host_pct = contam_pct * 100,        # INPUT host proportion, before filtering.
  expected_reads = total_reads,       # All *_reads columns count read PAIRS.
  read1_file = paste0(sample_names, "_1.fasta"),
  read2_file = paste0(sample_names, "_2.fasta"),
  potato_reads = NA_real_, bacteria_reads = NA_real_, fungi_reads = NA_real_,
  microbe_reads = NA_real_, assigned_reads = NA_real_,
  assigned_host_pct = NA_real_, fungi_pct_of_microbes = NA_real_,
  stringsAsFactors = FALSE
)
expr <- matrix(0, nrow = N, ncol = samples,
               dimnames = list(tx_names, sample_names))

set.seed(seed)
potato_profile <- sample(1:10, sum(is_potato), replace = TRUE)
microbe_profiles <- list(
  profile1 = sample(10:100, sum(!is_potato), replace = TRUE),
  profile2 = sample(10:100, sum(!is_potato), replace = TRUE)
)
profile_weights <- lapply(microbe_profiles, function(x) {
  full <- numeric(N)
  full[!is_potato] <- x
  full
})

for (i in seq_len(samples)) {
  reads_potato <- round(total_reads * contam_pct[i])
  reads_microbe <- total_reads - reads_potato
  expr[is_potato, i] <- allocate_counts(potato_profile, reads_potato)
  species_budgets <- allocate_counts(rep(1, length(species_rows)), reads_microbe)
  weights <- profile_weights[[profile_ids[i]]]
  for (s in seq_along(species_rows)) {
    rows <- species_rows[[s]]
    expr[rows, i] <- allocate_counts(weights[rows], species_budgets[s])
  }
  metadata$potato_reads[i] <- sum(expr[is_potato, i])
  metadata$bacteria_reads[i] <- sum(expr[is_bacteria, i])
  metadata$fungi_reads[i] <- sum(expr[is_fungi, i])
  metadata$microbe_reads[i] <- sum(expr[!is_potato, i])
  metadata$assigned_reads[i] <- sum(expr[, i])
  metadata$assigned_host_pct[i] <- sum(expr[is_potato, i]) / total_reads * 100
  metadata$fungi_pct_of_microbes[i] <- metadata$fungi_reads[i] / reads_microbe * 100
}


# === 4. Validate and save the simulation truth ===
stopifnot(nrow(expr) == N, ncol(expr) == samples,
          all(is.finite(expr)), all(expr >= 0), all(expr == floor(expr)),
          all(colSums(expr) == total_reads),
          all(expr[is_potato, contam_pct == 0, drop = FALSE] == 0),
          all(abs(metadata$assigned_host_pct - metadata$host_pct) < 1e-8))

# The truth table includes potato; microbial scoring excludes it.
species_info <- unique(transcript_sources[c("kingdom", "species_taxid", "species")])
species_counts <- t(vapply(species_info$species_taxid, function(id) {
  colSums(expr[transcript_sources$species_taxid == id, , drop = FALSE])
}, numeric(samples)))
colnames(species_counts) <- sample_names
for (i in seq_len(samples)) {
  microbial_counts <- species_counts[species_info$kingdom != "host", i]
  stopifnot(all(microbial_counts > 0),
            max(microbial_counts) - min(microbial_counts) <= 1)
}
species_truth <- cbind(species_info, as.data.frame(species_counts))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
combined_fasta <- file.path(outdir, "combined_transcripts.fa")
writeXStringSet(all_seqs, combined_fasta)
write.csv(sources, file.path(outdir, "source_assemblies.csv"), row.names = FALSE)
write.csv(transcript_sources, file.path(outdir, "transcript_sources.csv"), row.names = FALSE)
write.csv(species_truth, file.path(outdir, "species_truth_counts.csv"), row.names = FALSE)
write.csv(metadata, file.path(outdir, "metadata.csv"), row.names = FALSE)
write.csv(expr, file.path(outdir, "truth_counts.csv"), row.names = TRUE)
saveRDS(expr, file.path(outdir, "truth_counts.rds"))
saveRDS(list(seed = seed, total_read_pairs = total_reads, readlen = readlen,
             fraglen = fraglen, fragsd = fragsd, error_rate = error_rate,
             min_template_length = min_template_length,
             allocation = "Equal microbial species budgets; within-species transcript weights 10:100",
             host_weights = potato_profile, microbial_weights = microbe_profiles,
             polyester_version = as.character(packageVersion("polyester")),
             biostrings_version = as.character(packageVersion("Biostrings"))),
        file.path(outdir, "simulation_settings.rds"))
capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))
message("Pre-simulation checks passed: 22 bacterial + 8 fungal species.")
message("Fungi receive approximately 8/30 of non-host read pairs, by design.")


# === 5. Simulate paired-end FASTA reads ===
# The combined FASTA keeps template order aligned with readmat rows.
simulate_experiment_countmat(
  fasta = combined_fasta,
  readmat = expr,
  outdir = outdir,
  paired = TRUE,
  seed = seed,
  readlen = readlen,
  fraglen = fraglen,
  fragsd = fragsd,
  error_model = "uniform",
  error_rate = error_rate
)


# === 6. Rename outputs and verify read-pair counts ===
get_num <- function(x) {
  as.integer(sub("^sample_([0-9]+)_[12][.]fasta$", "\\1", basename(x)))
}
files_1 <- list.files(outdir, "^sample_[0-9]+_1[.]fasta$", full.names = TRUE)
files_2 <- list.files(outdir, "^sample_[0-9]+_2[.]fasta$", full.names = TRUE)
files_1 <- files_1[order(get_num(files_1))]
files_2 <- files_2[order(get_num(files_2))]
stopifnot(identical(get_num(files_1), seq_len(samples)),
          identical(get_num(files_2), seq_len(samples)))
final_1 <- file.path(outdir, metadata$read1_file)
final_2 <- file.path(outdir, metadata$read2_file)
stopifnot(!any(file.exists(c(final_1, final_2))))
stopifnot(all(file.rename(files_1, final_1)),
          all(file.rename(files_2, final_2)))

# Count FASTA headers without loading whole files.
count_fasta_reads <- function(path) {
  con <- file(path, open = "r")
  on.exit(close(con))
  n <- 0
  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (length(lines) == 0L) break
    n <- n + sum(startsWith(lines, ">"))
  }
  n
}
metadata$observed_mate1_reads <- vapply(final_1, count_fasta_reads, numeric(1))
metadata$observed_mate2_reads <- vapply(final_2, count_fasta_reads, numeric(1))
write.csv(metadata, file.path(outdir, "metadata.csv"), row.names = FALSE)
stopifnot(all(metadata$observed_mate1_reads == metadata$assigned_reads),
          all(metadata$observed_mate2_reads == metadata$assigned_reads))

message("Simulation and read-count checks completed. Output: ",
        normalizePath(outdir))
# Continue with the FASTA-compatible steps in upstream-analysis.md.
