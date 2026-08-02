#!/usr/bin/env bash
set -euo pipefail

diamond_dir="02.1_secondary_mapping"
ids_dir="${diamond_dir}/read_ids"
temporary_dir="${diamond_dir}/temporary_individually_filtered_reads"
filtered_dir="${diamond_dir}/filtered_reads"
summary_file="${diamond_dir}/secondary_filtering_summary.tsv"

mkdir -p "$ids_dir" "$temporary_dir" "$filtered_dir"

# Step 1: Summarize direct Viridiplantae BLASTX hits
echo "[1/4] Summarizing Viridiplantae BLASTX hits ..."

echo -e \
  "read_file\ttotal_reads\tunique_viridiplantae_hits\thit_rate_percent" \
  > "$summary_file"

for f in *mate1* *mate2*; do
  hit_file="${diamond_dir}/${f}.viridiplantae.blastx"

  if [[ ! -f "$hit_file" ]]; then
    echo "ERROR: DIAMOND output was not found: $hit_file" >&2
    exit 1
  fi

  total_reads=$(( $(wc -l < "$f") / 4 ))
  hit_reads=$(cut -f1 "$hit_file" | sort -u | wc -l)
  hit_rate=$(awk \
    -v hits="$hit_reads" \
    -v total="$total_reads" \
    'BEGIN {
      if (total > 0) {
        printf "%.2f", hits / total * 100
      } else {
        printf "NA"
      }
    }')

  echo -e \
    "${f}\t${total_reads}\t${hit_reads}\t${hit_rate}" \
    >> "$summary_file"
done

# Step 2: Extract IDs of Viridiplantae-hit reads
echo "[2/4] Extracting Viridiplantae-hit read IDs ..."

for hit_file in "${diamond_dir}"/*.viridiplantae.blastx; do
  read_file="$(basename "${hit_file%.viridiplantae.blastx}")"
  ids_file="${ids_dir}/${read_file}.viridiplantae_hit_ids.txt"

  cut -f1 "$hit_file" \
    | sort -u \
    > "$ids_file"
done

# Step 3: Remove Viridiplantae-hit reads from STAR-unmapped reads
echo "[3/4] Removing Viridiplantae-hit reads ..."

for f in *mate1* *mate2*; do
  ids_file="${ids_dir}/${f}.viridiplantae_hit_ids.txt"
  temporary_file="${temporary_dir}/${f}_nonplant.fq"

  if [[ ! -f "$ids_file" ]]; then
    echo "ERROR: Read-ID file was not found: $ids_file" >&2
    exit 1
  fi

  if [[ -s "$ids_file" ]]; then
    seqkit grep \
      -v \
      -f "$ids_file" \
      "$f" \
      -o "$temporary_file"
  else
    cp "$f" "$temporary_file"
  fi
done

# Step 4: Retain complete non-plant read pairs
echo "[4/4] Retaining complete non-plant read pairs ..."

for mate1 in "${temporary_dir}"/*mate1*_nonplant.fq; do
  mate2="${mate1/mate1/mate2}"

  if [[ ! -f "$mate2" ]]; then
    echo "ERROR: Corresponding mate2 file was not found: $mate2" >&2
    exit 1
  fi

  seqkit pair \
    -1 "$mate1" \
    -2 "$mate2" \
    -O "$filtered_dir" \
    --save-unpaired
done

# Remove unpaired reads and temporary independently filtered files
echo "Removing intermediate and unpaired reads ..."
rm -f "${filtered_dir}"/*.unpaired.fq
rm -rf "$temporary_dir"

echo "Done."
echo "Summary: $summary_file"
echo "Final paired non-plant reads: $filtered_dir/"