# Validation of Polyester-generated mock RNA-seq reads

Run these checks with Bash from the mock working directory created in [step 01](01_RNA_simulation_preparations.md). The simulation must have finished writing its paired FASTA files, `metadata.csv`, and `transcript_sources.csv`.

## Inspect the simulated FASTA files

There should be six samples and 12 mate files. The combined template FASTA is not part of this check.

```bash
cd polyester_out_bacteria_fungi
seqkit stats -a -T ./*_[12].fasta > fasta_stats.tsv
column -t -s $'\t' fasta_stats.tsv
```

## Verify paired IDs and source counts

Run this block from the simulation output directory. It reads filenames and counts by column name from `metadata.csv`, checks the IDs of both mates in order, and uses the template ID in each read header to look up its source in `transcript_sources.csv`. Count each pair once, using mate 1 after confirming the IDs match.

```bash
python3 - <<'PY'
import csv
from collections import Counter
from itertools import zip_longest

with open("metadata.csv", newline="", encoding="utf-8") as handle:
    samples = list(csv.DictReader(handle))
with open("transcript_sources.csv", newline="", encoding="utf-8") as handle:
    sources = {
        row["transcript_id"]: row["kingdom"] for row in csv.DictReader(handle)
    }

if len(samples) != 6 or len({row["sample"] for row in samples}) != 6:
    raise SystemExit("Expected six distinct mock samples in metadata.csv.")


def fasta_ids(path):
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                yield line[1:].split()[0]


source_columns = {
    "host": "potato_reads",
    "bacteria": "bacteria_reads",
    "fungi": "fungi_reads",
}
for row in samples:
    sample = row["sample"]
    counts = Counter()
    for r1, r2 in zip_longest(fasta_ids(row["read1_file"]), fasta_ids(row["read2_file"])):
        if r1 != r2:
            raise SystemExit(f"{sample}: mate IDs or record counts differ.")
        # Polyester uses read_id/template_id;... in each simulated FASTA header.
        template = r1.partition("/")[2].split(";", 1)[0]
        kingdom = sources.get(template)
        if kingdom not in source_columns:
            raise SystemExit(f"{sample}: unknown source template in {r1}")
        counts[kingdom] += 1

    pairs = sum(counts.values())
    if pairs == 0 or pairs != float(row["assigned_reads"]):
        raise SystemExit(f"{sample}: read-pair count differs from metadata.")
    for kingdom, column in source_columns.items():
        if counts[kingdom] != float(row[column]):
            raise SystemExit(f"{sample}: {kingdom} read-pair count differs from metadata.")

    print(f"PASS: {sample}: {pairs:,} pairs; "
          f"potato={counts['host']:,} ({100 * counts['host'] / pairs:.2f}%); "
          f"bacteria={counts['bacteria']:,}; fungi={counts['fungi']:,}")
PY
```

These checks validate the simulated input files, not classifier accuracy. Once all six samples pass, continue with [04_Upstream-analysis.md](04_Upstream-analysis.md).
