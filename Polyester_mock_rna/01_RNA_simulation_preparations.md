# Preparation of reference sequences for mock RNA-seq simulation
## This document describes the preparation of microbial CDS templates and host transcript sequences used for Polyester simulation. 



# Install necessary tools
conda create -n mock_rna \
    -c conda-forge \
    -c bioconda \
    ncbi-datasets-cli=16.6.0 \
    seqkit=2.13.0


# Create directory for storing source sequences
mkdir -p Data
cd Data

# Download microbial sequences
## The microbial CDS FASTA files used in the original simulation were obtained from a local CAMISIM genome collection. The NCBI Datasets commands below provide a reproducible route for retrieving CDSs corresponding to the same versioned assembly accessions.
## Save the following 24 versioned assembly accessions, one per line, as microbial_assembly_accessions.txt
<!--
GCA_000006785.2
GCA_000013285.1
GCA_000013785.1
GCA_000015865.1
GCA_000022345.1
GCA_000023785.1
GCA_000025905.1
GCA_000092125.1
GCA_000092825.1
GCA_000092985.1
GCA_000143845.1
GCA_000143985.1
GCA_000227705.3
GCA_000231385.3
GCA_000233715.3
GCA_000235405.3
GCA_000242255.3
GCA_000255115.3
GCA_000265425.1
GCA_000325705.1
GCA_000439255.1
GCA_000445015.1
GCA_001692755.1
GCA_001877055.1
-->

## Fetch using NCBI tool
datasets download genome accession \
    --inputfile microbial_assembly_accessions.txt \
    --include cds  
## Decompress
unzip ncbi_dataset.zip -d ncbi_cds
## Combine microbial sequence
cat ncbi_cds/ncbi_dataset/data/GCA_*/*cds_from_genomic.fna  > microbes_all.fna


# Download the host reference transcriptome
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


# NEXT-STEP: Run the simulation
## Open `02_polyester_mock_rna.R` in R or RStudio. Set the working directory to the project root directory containing the `Data` folder, and run the script section by section.