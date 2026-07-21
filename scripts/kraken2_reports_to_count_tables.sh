#!/bin/bash

# Generate a combined read-count matrix and rank-specific count tables from multiple Kraken2 report files.
# Run this script in the directory containing the Kraken2 *.report files.

# All report files should have been generated using:
# 1. The same Kraken2 database build.
# 2. The same Kraken2 classification settings.

# All MPA files generated below use the same kreport2mpa.py options.



# Read the KrakenTools path from the command line
if [[ $# -ne 1 ]]; then
  echo "Usage: bash $0 /path/to/KrakenTools" >&2
  exit 1
fi

KRAKENTOOLS_DIR="$1"
# Check the required KrakenTools scripts
if [[ ! -f "$KRAKENTOOLS_DIR/kreport2mpa.py" ]]; then
  echo "Cannot find: $KRAKENTOOLS_DIR/kreport2mpa.py" >&2
  exit 1
fi

if [[ ! -f "$KRAKENTOOLS_DIR/combine_mpa.py" ]]; then
  echo "Cannot find: $KRAKENTOOLS_DIR/combine_mpa.py" >&2
  exit 1
fi

### Check for Kraken2 report files
if ! compgen -G "*.report" > /dev/null; then
  echo "No Kraken2 *.report files were found in the current directory." >&2
  exit 1
fi

### Create an output directory for MPA-style reports
mkdir -p "mpa_report"

# Use an empty mpa_report directory when rerunning this script.
# Existing *.mpa files from previous runs could otherwise be included
# when combine_mpa.py processes all files matching *.mpa.

# Step 1: Convert each Kraken2 report to MPA format
echo "[Step 1] Converting .report files to .mpa format..."

for file in *.report; do
  # Remove the .report suffix while preserving the rest of the filename.
  sample="${file%.report}"

  python "$KRAKENTOOLS_DIR/kreport2mpa.py" \
    -r "$file" \
    -o "mpa_report/${sample}.mpa" \
    --display-header

# Step 2: Combine all MPA files into one read-count matrix
echo "[Step 2] Combining all .mpa files into combine.mpa..."

cd "mpa_report" || exit 1

python "$KRAKENTOOLS_DIR/combine_mpa.py" \
  -i *.mpa \
  -o "combine.mpa"

# Step 3: Extract rank-specific bacterial and fungal count tables
echo "[Step 3] Extracting taxonomic levels from combine.mpa..."
# These rank-specific tables are intentionally written without header lines.
# Sample names are added later during downstream processing in R.

awk -F "|" '/k__Bacteria/ {print $NF}' combine.mpa \
  | grep "^f__" \
  | sed 's/^f__//' \
  | sort -u > "bac_f.txt"

awk -F "|" '/k__Bacteria/ {print $NF}' combine.mpa \
  | grep "^g__" \
  | sed 's/^g__//' \
  | sort -u > "bac_g.txt"

awk -F "|" '/k__Bacteria/ {print $NF}' combine.mpa \
  | grep "^s__" \
  | sed 's/^s__//' \
  | sort -u > "bac_sp.txt"

awk -F "|" '/k__Fungi/ {print $NF}' combine.mpa \
  | grep "^f__" \
  | sed 's/^f__//' \
  | sort -u > "fungi_f.txt"

awk -F "|" '/k__Fungi/ {print $NF}' combine.mpa \
  | grep "^g__" \
  | sed 's/^g__//' \
  | sort -u > "fungi_g.txt"

awk -F "|" '/k__Fungi/ {print $NF}' combine.mpa \
  | grep "^s__" \
  | sed 's/^s__//' \
  | sort -u > "fungi_sp.txt"

echo "All done. Output files are saved in: $(pwd)"
