# 04_Assembly_and_contig-level_taxonomic_classification
# This workflow uses the paired non-plant reads generated in Step 02.1. The reads are assembled with Trinity, the assembled contigs are aligned against bacterial and fungal protein sequences with DIAMOND, and the resulting DAA files are taxonomically assigned and exported with MEGAN6 for overlap extraction and downstream analysis.

## 04.1_De_novo_assembly_with_Trinity
### Trinity installation
# Trinity v2.14.0 was used in this workflow.
#### Download Trinity package
wget "https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.14.0/trinityrnaseq-v2.14.0.FULL_with_extendedTestData.tar.gz"
#### Decompress
tar -xzf "trinityrnaseq-v2.14.0.FULL_with_extendedTestData.tar.gz"
#### Build Trinity (with plugins) from source 
cd "/path/to/trinityrnaseq-v2.14.0"
make
make plugins

### Run Trinity to assemble paired non-plant reads
# This step uses the paired-end non-plant reads generated after host-read removal and secondary plant-read filtering.
# The CPU and memory settings shown below were used in this workflow and should be adjusted according to dataset size and available computing resources.
# --NO_SEQTK was used in the original workflow as a workaround for a seqtk-related error. It disables seqtk-based FASTQ-to-FASTA conversion  and uses Trinity's slower Perl implementation instead. This option may be omitted if the default seqtk-based conversion works.
# --full_cleanup removes intermediate assembly files after a successful run and retains the final assembled FASTA file.

#### Enter the directory containing the paired non-plant reads
cd "/path/to/STAR_output/02.1_secondary_mapping/filtered_reads"
#### Create the assembly output directory
mkdir -p "/path/to/04_assembled_contigs"

#### Run Trinity
for file in *mate1*_nonplant.fq; do
  prefix="${file%mate1*}"
  suffix="${file#*mate1}"
  mate2="${prefix}mate2${suffix}"

  if [[ ! -f "$mate2" ]]; then
    echo "Paired mate2 file was not found for: $file" >&2
    exit 1
  fi

  "/path/to/software/trinityrnaseq-v2.14.0/Trinity" \
    --seqType fq \
    --left "$file" \
    --right "$mate2" \
    --CPU 16 \
    --max_memory 32G \
    --NO_SEQTK \
    --output "/path/to/04_assembled_contigs/${prefix}trinity" \
    --full_cleanup
done




## 04.2_Protein_alignment_of_assembled_contigs_with_DIAMOND
### Download and extract DIAMOND
# DIAMOND v2.0.15 was used in this workflow.
# The command below downloads the precompiled Linux binary.
mkdir -p "/path/to/software/diamond-v2.0.15"
cd "/path/to/software/diamond-v2.0.15"
wget "https://github.com/bbuchfink/diamond/releases/download/v2.0.15/diamond-linux64.tar.gz"
tar -xzf "diamond-linux64.tar.gz"
# Check the installed version.
./diamond version

### Build a DIAMOND-formatted protein database
# The bacterial and fungal protein FASTA files were generated during Kraken2 database preparation.
# The two FASTA files are concatenated with the cat command and streamed directly to the DIAMOND makedb command through standard input.
mkdir -p "/path/to/DIAMOND_database"
cat \
  "/path/to/kraken2_database_preparation/nr_sequences/nr_bacteria.fa" \
  "/path/to/kraken2_database_preparation/nr_sequences/nr_fungi.fa" \
  | "/path/to/software/diamond-v2.0.15/diamond" makedb \
    --in - \
    --db "/path/to/DIAMOND_database/nr_bacteria_and_fungi" \
    --threads 8
# Output database: /path/to/DIAMOND_database/nr_bacteria_and_fungi.dmnd

### Align assembled contigs against the DIAMOND-formatted protein database
# The --long-reads option was used because assembled contigs are long query sequences that may contain multiple coding regions.
# The --outfmt 100 option writes the alignments in DAA format for subsequent processing with MEGAN6.
# Standard output and error messages from the complete batch run are saved in a log file for checking run completion and troubleshooting.
cd "/path/to/04_assembled_contigs"
mkdir -p "/path/to/DIAMOND_output"

for file in *.Trinity.fasta; do
  contig_stem="${file%.Trinity.fasta}"

  echo "Processing: $file"

  "/path/to/software/diamond-v2.0.15/diamond" blastx \
    --db "/path/to/DIAMOND_database/nr_bacteria_and_fungi.dmnd" \
    --query "$file" \
    --out "/path/to/DIAMOND_output/${contig_stem}.daa" \
    --outfmt 100 \
    --long-reads \
    --threads 8
done > "/path/to/DIAMOND_output/diamond_blastx.log" 2>&1



## 04.3_Meganize_DIAMOND_DAA_files_for_taxonomic_assignment
### Download and install MEGAN Community Edition
# MEGAN Community Edition v6.23.4 was used in this workflow.
# The installer below is intended for Linux and Unix systems.
# Users of other operating systems should download the corresponding installer from the MEGAN6 download page.
cd "/path/to/software"
wget "https://software-ab.informatik.uni-tuebingen.de/download/megan6/MEGAN_Community_unix_6_23_4.sh"
# Run the interactive Linux/Unix installer.
sh "MEGAN_Community_unix_6_23_4.sh"

### Download and extract the MEGAN accession mapping database
# The February 2022 MEGAN mapping database was used in this workflow.
# It maps NCBI-nr accessions to taxonomic and functional classifications,
# including NCBI taxonomy, GTDB, EC, eggNOG, InterPro2GO, and SEED.
mkdir -p "/path/to/MEGAN_mapping_database"
cd "/path/to/MEGAN_mapping_database"
wget "https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-map-Feb2022.db.zip"
unzip "megan-map-Feb2022.db.zip"
# Extracted mapping database:
# /path/to/MEGAN_mapping_database/megan-map-Feb2022.db

### Meganize DIAMOND DAA files
# The daa-meganizer program is located in the tools directory of the MEGAN6 package.
# Run this step in the directory containing the DIAMOND DAA files.
# The -lg option enables long-read analysis mode, which was used because the input sequences are assembled contigs.
# The -t option sets the number of threads used by daa-meganizer.
# The -v option enables verbose output.
# daa-meganizer appends the classification results and the required indices directly to each DAA file. It does not generate a separate output file.
cd "/path/to/DIAMOND_output"
"/path/to/megan/tools/daa-meganizer" \
  -i *.daa \
  -mdb "/path/to/MEGAN_mapping_database/megan-map-Feb2022.db" \
  -t 8 \
  -lg \
  -v


### Generating read count matrix using MEGAN6
### [In GUI]Generate a rank-level taxonomic count matrix using the MEGAN6 GUI
# Step 1: Create a comparison document
Select **File → Compare...** and add all meganized DAA files to be
included in the comparison.
Select **Use Absolute Counts** to retain the original taxonomic counts
for each sample, and then create the comparison document.
⬇️
# Step 2: Project taxonomic assignments to the required rank
In the comparison document, select:
**Options → Project Assignments To Rank...**
Choose the taxonomic rank used in the downstream analysis.
Repeat the projection and export steps separately for the Family, Genus, and Species ranks.
# Taxonomic assignments in MEGAN may occur at different levels of the taxonomy. Projection generates a taxonomic profile represented at a consistent rank for subsequent count-matrix export and sample comparison.
⬇️
# Step 3: Select bacterial or fungal taxa
In the projected taxonomy document, select either the **Bacteria** node or the **Fungi** node. Bacterial and fungal count matrices are exported separately.
⬇️
# Step 4: Select the taxa below the selected nodes
Select:
**Select → Leaves Below**
This selects the leaf taxa displayed below the selected Bacteria and
Fungi nodes.
⬇️
# Step 5: Export the taxonomic count matrix
Select:
**File → Export → Text (CSV) Format...**
Export the selected taxonomic counts using a tab-separated text format.
Repeat Steps 2–5 for each taxonomic rank and organism group to generate six tab-separated count matrices: family-, genus-, and species-level matrices for bacteria and fungi.