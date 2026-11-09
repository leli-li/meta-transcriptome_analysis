Replace the `/path/to/...` placeholders with your local software and database paths. Run the command blocks with Bash. Wait for each background job to finish successfully before starting a step that uses its outputs.

# 01_trimming
Simulated reads are in FASTA format, so Trimmomatic trimming is skipped.



# 02_mapping
Run from the directory containing the simulation output folder.

```bash
cd polyester_out_bacteria_fungi
mkdir -p 02_mapping

nohup bash -c 'set -e; for file in *_1.fasta; do sample="${file%_1.fasta}"; STAR --runThreadN 16 --genomeDir "/path/to/potato_star_index" --readFilesIn "$file" "${sample}_2.fasta" --outFileNamePrefix "02_mapping/${sample}_" --outSAMtype None --outReadsUnmapped Fastx; done' > STAR.log 2>&1 &
```

# 02.1_secondary_filtering
After STAR has completed for all six samples, enter `02_mapping/`. Run both secondary-filtering blocks from that directory.

```bash
cd 02_mapping
```

## DIAMOND blastx
```bash
#!/bin/bash
set -euo pipefail

OUTDIR="02.1_secondary_mapping"
mkdir -p "$OUTDIR"

for mate1 in *mate1*; do
    mate2="${mate1/mate1/mate2}"

    if [[ ! -f "$mate2" ]]; then
        echo "$(date '+%F %T') [ERROR] Paired mate2 file was not found for: $mate1" >&2
        exit 1
    fi

    echo "$(date '+%F %T') [INFO] Processing: $mate1 and $mate2"

    "/path/to/diamond" blastx \
        --db "/path/to/viridiplantae_swissprot.dmnd" \
        --query "$mate1" \
        --out "${OUTDIR}/${mate1}.viridiplantae.blastx" \
        --evalue 1e-5 \
        --threads 8 \
        --outfmt 6 qseqid sseqid pident evalue bitscore

    "/path/to/diamond" blastx \
        --db "/path/to/viridiplantae_swissprot.dmnd" \
        --query "$mate2" \
        --out "${OUTDIR}/${mate2}.viridiplantae.blastx" \
        --evalue 1e-5 \
        --threads 8 \
        --outfmt 6 qseqid sseqid pident evalue bitscore

    echo "$(date '+%F %T') [INFO] Done: $mate1 / $mate2"
done

echo "$(date '+%F %T') [INFO] All pairs finished."
```

## Filtering
```bash
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

  total_reads=$(seqkit seq -n "$f" | wc -l)
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
  temporary_file="${temporary_dir}/${f}_nonplant.fasta"

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

for mate1 in "${temporary_dir}"/*mate1*_nonplant.fasta; do
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
rm -f "${filtered_dir}"/*.unpaired.fasta
rm -rf "$temporary_dir"

echo "Done."
echo "Summary: $summary_file"
echo "Final paired non-plant reads: $filtered_dir/"
```

## Check paired reads
Proceed only if the paired files have matching IDs in the same order.

```bash
cd 02.1_secondary_mapping/filtered_reads
for mate1 in *mate1_nonplant.fasta; do
    mate2="${mate1%mate1_nonplant.fasta}mate2_nonplant.fasta"
    if [[ ! -f "$mate2" ]]; then
        echo "Missing mate2: $mate2" >&2
        exit 1
    fi
    diff -q <(grep '^>' "$mate1" | cut -d' ' -f1) \
            <(grep '^>' "$mate2" | cut -d' ' -f1) || exit 1
    echo "PASS: $mate1"
done
```

# Kraken2
Run Kraken2 and Trinity from the paired non-plant FASTA directory.

```bash
mkdir -p 03_kraken2_output
nohup bash -c 'set -e; for mate1 in *_Unmapped.out.mate1_nonplant.fasta; do sample="${mate1%_Unmapped.out.mate1_nonplant.fasta}"; kraken2 --db "/path/to/nr_bac_and_fungi" --threads 16 --paired --fasta-input "$mate1" "${sample}_Unmapped.out.mate2_nonplant.fasta" --report "03_kraken2_output/${sample}_kk2_nr.report" --output "03_kraken2_output/${sample}_kk2_nr.out" --confidence 0.05; done' > kraken2_batch_new.log 2>&1 &
```

# Trinity
```bash
mkdir -p 04_assembled_contigs
nohup bash -c 'for file in *mate1*_nonplant.fasta; do prefix="${file%mate1*}"; suffix="${file#*mate1}"; mate2="${prefix}mate2${suffix}"; if [[ ! -f "$mate2" ]]; then echo "Paired mate2 file was not found for: $file" >&2; exit 1; fi; Trinity --seqType fa --left "$file" --right "$mate2" --CPU 16 --max_memory 32G --NO_SEQTK --no_normalize_reads --output "04_assembled_contigs/${prefix}trinity" --full_cleanup; done' > 'trinity.log' 2>&1 &
```

# DIAMOND blastx
After Trinity has completed, enter the directory containing the assembled contigs.

```bash
cd 04_assembled_contigs
mkdir -p DIAMOND_output
nohup bash -c 'set -e; for file in *.Trinity.fasta; do contig_stem="${file%.Trinity.fasta}"; echo "Processing: $file"; "/path/to/diamond" blastx --db "/path/to/nr_bac_and_fungi.dmnd" --query "$file" --out "DIAMOND_output/${contig_stem}.daa" --outfmt 100 --long-reads --threads 8; done' > "diamond_blastx.log" 2>&1 &
```

# Meganize
After DIAMOND has completed, run the following in its output directory. The manual commands meganize DAA files in place.

```bash
cd DIAMOND_output
for file in *.daa; do
    "/path/to/megan/tools/daa-meganizer" -i "$file" \
        -mdb "/path/to/megan-map-Feb2022.db" -t 8 -lg -ram readCount -v || exit 1
done
```

## Export NCBI taxonomic counts
Run from the same directory after meganization has completed.

```bash
mkdir -p NCBI_c2c
for file in *.daa; do
    "/path/to/megan/tools/daa2info" -i "$file" -c2c Taxonomy -n -r \
        -o "NCBI_c2c/${file%.daa}_NCBI_c2c.txt" || exit 1
done
```

## Optional: Generate a MEGAN comparison file
This file is not needed for the c2c count-table route.

```bash
"/path/to/megan/tools/compute-comparison" -i *.daa -o mock-rna_absolute.megan -n false
```

## Next steps
Use the Kraken2 `.report` files and MEGAN `NCBI_c2c/` exports to generate the respective count matrices (Workflow 03 and 04). Overlap extraction and benchmark scoring are separate downstream steps.
