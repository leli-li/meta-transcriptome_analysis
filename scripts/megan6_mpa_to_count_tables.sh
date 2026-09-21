#!/bin/bash

# Generate a combined contig-count matrix and rank-specific count tables
# from multiple MEGAN6 MPA files.
# Run this script in the directory containing the *.megan.mpa files.


# Read the KrakenTools path from the command line
if [[ $# -ne 1 ]]; then
  echo "Usage: bash $0 /path/to/KrakenTools" >&2
  exit 1
fi

KRAKENTOOLS_DIR="$1"

# Check the required KrakenTools script
if [[ ! -f "$KRAKENTOOLS_DIR/combine_mpa.py" ]]; then
  echo "Cannot find: $KRAKENTOOLS_DIR/combine_mpa.py" >&2
  exit 1
fi

### Check for MEGAN6 MPA files
if ! compgen -G "*.megan.mpa" > /dev/null; then
  echo "No MEGAN6 *.megan.mpa files were found in the current directory." >&2
  exit 1
fi

# Step 1: Add the filename as the sample-name header
echo "[Step 1] Adding sample-name headers to MEGAN6 MPA files..."

for file in *.megan.mpa; do
  sample="${file%.megan.mpa}"

  if ! head -n 1 "$file" | grep -q '^#Classification'; then
    {
      printf '#Classification\t%s\n' "$sample"
      cat "$file"
    } > "${file}.tmp"

    mv "${file}.tmp" "$file"
  fi
done

# Step 2: Combine all MPA files into one contig-count matrix
echo "[Step 2] Combining all .megan.mpa files into combine.mpa..."

python3 "$KRAKENTOOLS_DIR/combine_mpa.py" \
  -i *.megan.mpa \
  -o "combine.mpa"

# Step 3: Extract rank-specific bacterial and fungal count tables
echo "[Step 3] Extracting taxonomic levels from combine.mpa..."

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

# Add the header from combine.mpa to each rank-specific table
for table in bac_f.txt bac_g.txt bac_sp.txt fungi_f.txt fungi_g.txt fungi_sp.txt; do
  {
    head -n 1 combine.mpa
    cat "$table"
  } > "${table}.tmp"

  mv "${table}.tmp" "$table"
done

echo "All done. Output files are saved in: $(pwd)"
