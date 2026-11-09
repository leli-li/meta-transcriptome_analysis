# Preparation of reference sequences for mock RNA-seq simulation
Prepare the bacterial and fungal CDS templates and potato transcripts used by the Polyester simulation.

## Install the required tools

```bash
conda create -n mock_rna \
    -c conda-forge \
    -c bioconda \
    ncbi-datasets-cli=16.6.0 \
    seqkit=2.13.0
conda activate mock_rna
```

## Prepare the working directory

Use a new working directory outside the Git checkout. Set the two paths below and run the command blocks in the same Bash session.

```bash
repo_dir="/path/to/meta-transcriptome_analysis"
mock_dir="/path/to/mock_rna_project"
mkdir -p "$mock_dir/data" "$mock_dir/fungi"
cd "$mock_dir"
```

## Download bacterial CDSs

The original bacterial CDS files came from a local CAMISIM genome collection. The commands below retrieve CDSs for the same versioned assembly accessions listed in [microbial_assembly_accessions.txt](Data/microbial_assembly_accessions.txt). Keep each assembly in a separate file, as required by the simulation script.

```bash
bacterial_accessions="$repo_dir/Polyester_mock_rna/Data/microbial_assembly_accessions.txt"
datasets download genome accession \
    --inputfile "$bacterial_accessions" --include cds --filename data/bacterial_cds.zip
unzip data/bacterial_cds.zip -d data/ncbi_cds

while IFS= read -r accession || [[ -n "$accession" ]]; do
    [[ -n "$accession" ]] || continue
    cds="data/ncbi_cds/ncbi_dataset/data/$accession/cds_from_genomic.fna"
    if [[ ! -s "$cds" ]]; then
        echo "Missing or empty bacterial CDS file: $cds" >&2
        exit 1
    fi
    cp "$cds" "data/${accession}.cds.fna" || exit 1
done < "$bacterial_accessions"
```

## Download fungal CDSs

Save the following eight versioned assembly accessions, one per line, as `fungi/fungal_assembly_accessions.txt`:
```text
GCF_023522795.1
GCF_000149615.1
GCF_033870435.1
GCF_000150705.2
GCF_003942205.1
GCF_000271585.1
GCF_001720195.1
GCF_000320585.1
```
Keep the extracted NCBI directory structure under `fungi/fungal_cds/`.

```bash
datasets download genome accession \
    --inputfile fungi/fungal_assembly_accessions.txt --include cds \
    --filename fungi/fungal_cds.zip
unzip fungi/fungal_cds.zip -d fungi/fungal_cds
```

## Prepare potato transcripts

```bash
cd data
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Solanum_tuberosum/latest_assembly_versions/GCF_000226075.1_SolTub_3.0/GCF_000226075.1_SolTub_3.0_rna.fna.gz
## Decompress
gunzip GCF_000226075.1_SolTub_3.0_rna.fna.gz

# Inspect the host reference transcriptome
seqkit stats GCF_000226075.1_SolTub_3.0_rna.fna

# Retain transcripts at least 150 nt long
seqkit seq -m 150 \
    GCF_000226075.1_SolTub_3.0_rna.fna \
    -o potato_rna_min150.fna

# Randomly sample 100 host transcripts
seqkit sample2 -n 100 -s 56 -2 \
    potato_rna_min150.fna \
    -o random_potato_100.fna

# Verify the sampled transcript count and length distribution
seqkit fx2tab -n -l random_potato_100.fna \
    | sort -t $'\t' -k2,2n \
    | head
cd ..
```

## Run the simulation

The simulation expects `data/random_potato_100.fna`, one `data/<accession>.cds.fna` per bacterial assembly, and `fungi/fungal_cds/ncbi_dataset/data/<accession>/cds_from_genomic.fna` per fungal assembly. Do not merge the bacterial CDS files.

With `polyester` and `Biostrings` installed, run from the working directory containing `data/` and `fungi/`:

```bash
Rscript "$repo_dir/Polyester_mock_rna/02_polyester_mock_rna.R"
```

Check the resulting files using [03_Mock_RNA_check.md](03_Mock_RNA_check.md).
