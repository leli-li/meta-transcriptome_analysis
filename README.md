# Meta-transcriptome analysis workflow

This repository documents the analysis workflow used for host-filtered metatranscriptomic profiling of crop-associated microbial communities.

The workflow includes read trimming, host read removal, secondary plant-read filtering, taxonomic classification, assembly-based classification, and downstream functional annotation.

## Workflow overview

The main analysis steps are:

1. Quality trimming of paired-end RNA-seq reads using Trimmomatic.
2. Host read removal using STAR.
3. Secondary plant-read filtering using DIAMOND blastx against reviewed Viridiplantae Swiss-Prot proteins.
4. Read-level taxonomic classification using Kraken2.
5. De novo assembly of non-plant reads using Trinity.
6. Contig-level taxonomic classification using DIAMOND and MEGAN6.
7. Functional annotation using FAPROTAX and FUNGuild.

## Repository scope

This repository provides workflow documentation and curated analysis command records for the associated manuscript.

It is not a fully automated software package. Some steps, especially database preparation and MEGAN6-based taxonomic export, may require adaptation to the local computing environment.

## Notes

Large sequencing files, intermediate outputs, and reference databases are not included in this repository.

