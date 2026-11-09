# Meta-transcriptome analysis workflow
This repository documents the analysis workflow used for host-filtered metatranscriptomic profiling of crop-associated microbial communities.

The workflow includes read trimming, host read removal, secondary plant-read filtering, taxonomic classification, assembly-based classification, and downstream community analyses.

## Workflow overview
The main analysis steps are:

1. Quality trimming of paired-end RNA-seq reads using Trimmomatic. [Trimming and STAR commands](workflow/01_Trimming-02_Host_read_removal_with_STAR.md)
2. Host read removal using STAR.
3. Secondary plant-read filtering using DIAMOND blastx against reviewed Viridiplantae Swiss-Prot proteins. [Secondary filtering commands](workflow/02.1_Secondary_plant-read_filtering_with_DIAMOND.md)
4. Read-level taxonomic classification using Kraken2 and generation of read-pair count matrices. [Kraken2 workflow](Workflow/03_Kraken2_classification.md)
5. De novo assembly of non-plant reads using Trinity.
6. Contig-level taxonomic classification using DIAMOND and MEGAN6 and generation of contig-count matrices. [Assembly and contig classification workflow](workflow/04_Assembly_and_contig-level_taxonomic_classification.md)
7. Identification of taxa shared between the two workflows, followed by downstream analyses including Shannon diversity and PCA/PERMANOVA.

## Mock RNA-seq
The mock dataset contains simulated paired-end reads from potato, 22 bacterial species, and eight fungal species. The mock workflow uses paired FASTA files and skips Trimmomatic.

- [Prepare the reference sequences](Polyester_mock_rna/01_RNA_simulation_preparations.md)
- [Simulate reads with Polyester](Polyester_mock_rna/02_polyester_mock_rna.R)
- [Check the simulated reads](Polyester_mock_rna/03_Mock_RNA_check.md)
- [Process and classify the simulated reads](Polyester_mock_rna/04_Upstream-analysis.md)

## Downstream analysis
The [overlap script](scripts/find_overlap.py) identifies shared taxon names from the Kraken2 and DIAMOND–MEGAN6 count tables. Shared [R functions](Downstream_analysis/R_scripts/) are provided for Shannon diversity and PCA/PERMANOVA. To apply them to another dataset, prepare the count matrices and sample metadata, and adapt the sample-ID mapping and model formula to the study design.

## Repository scope
This repository provides workflow documentation and curated analysis command records for the associated manuscript.

It is not a fully automated software package. Some steps, especially database preparation and MEGAN6-based taxonomic export, may require adaptation to the local computing environment.

## Notes
Large sequencing files, intermediate outputs, and reference databases are not included in this repository.