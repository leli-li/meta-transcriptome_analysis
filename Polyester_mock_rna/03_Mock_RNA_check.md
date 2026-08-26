# Validation of Polyester-generated mock RNA-seq reads
## Run the following checks from the Polyester output directory.
cd path/to/output_directory



# 1. Inspect file count, format, and read lengths
find . -maxdepth 1 -name "*.fasta" | sort
find . -maxdepth 1 -name "*.fasta" | wc -l

seqkit stats -a -T ./*.fasta \
  > "fasta_stats.tsv"

column -t -s $'\t' "fasta_stats.tsv"



# 2. Verify mate counts and read-ID pairing
for r1 in ./*_1.fasta; do
    r2="${r1%_1.fasta}_2.fasta"

    n1=$(grep -c '^>' "$r1")
    n2=$(grep -c '^>' "$r2")

    printf "%s\tR1=%s\tR2=%s\n" "$(basename "${r1%_1.fasta}")" "$n1" "$n2"

    test "$n1" -eq "$n2" || {
        echo "ERROR: mate counts differ"
        exit 1
    }

    cmp -s \
      <(seqkit seq -n -i "$r1") \
      <(seqkit seq -n -i "$r2") || {
        echo "ERROR: mate IDs differ for $r1"
        exit 1
    }
done


# 3. Verify FASTA record counts against metadata
tail -n +2 metadata.csv | sed 's/"//g' | cut -d',' -f1,10 | tr ',' '\t' | sort > meta_pairs.txt
tail -n +2 fasta_stats.tsv | cut -f1,4 | sed -E 's#^\./##; s/_[12]\.fasta\t/\t/' | sort -u > fasta_pairs.txt
if comm -3 meta_pairs.txt fasta_pairs.txt | grep -q .; then
    echo "ERROR: FASTA record counts do not match metadata." >&2
    exit 1
else
    echo "PASS: all FASTA record counts match metadata."
fi



# 4. Verify host-read proportions from the final FASTA headers
seqkit seq -n -i ../Data/random_potato_100.fna \
  | sort -u \
  > potato_transcript_ids.txt

for r1 in ./*_1.fasta; do
    total=$(grep -c '^>' "$r1")

    host=$(
      seqkit seq -n -i "$r1" \
        | sed -E 's#^[^/]*/##; s#;.*$##' \
        | awk '
            NR == FNR { host[$1] = 1; next }
            ($1 in host) { n++ }
            END { print n + 0 }
          ' potato_transcript_ids.txt -
    )

    pct=$(awk -v h="$host" -v n="$total" \
      'BEGIN { printf "%.6f", 100 * h / n }')

    printf "%s\ttotal=%s\thost=%s\thost_pct=%s%%\n" \
      "$(basename "${r1%_1.fasta}")" "$total" "$host" "$pct"
done