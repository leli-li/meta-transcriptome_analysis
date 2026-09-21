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

  Trinity \
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

### Generating contig-count matrices
This workflow exports taxonomic counts from meganized DAA files, converts them to MPA format, and combines the per-sample reports into rank-specific count matrices. Because the DIAMOND queries are assembled contigs, the resulting matrices contain contig counts rather than read-abundance estimates.
#### Create Conda environment
conda create -n megan2mpa -c conda-forge python=3.10 ete3 pandas -y
conda activate megan2mpa
#### Export NCBI taxonomic class-to-count files
# Run the following commands in the directory containing the meganized DAA files:
mkdir -p NCBI_c2c
for file in *.daa; do
    /path/to/daa2info \
        -i "$file" \
        -c2c Taxonomy \
        -n \
        -r \
        -o "NCBI_c2c/${file%.daa}_NCBI_c2c.txt"
done

### Convert MEGAN c2c files to MPA format
# MEGAN taxonomic class-to-count (c2c) files were converted to MPA-style count tables using a locally modified version of `Convert_MEGAN-NCBI-c2c_to_kreport-mpa.py` from [PacBio’s pb-metagenomics-tools](https://github.com/PacificBiosciences/pb-metagenomics-tools). 
# Find the The modified script, `Convert_MEGAN-NCBI-c2c_to_kreport-mpa_revised.py` in the `Scripts/` directory. 

#### Run converter
# Use an empty mpa_report/ directory before running the converter, as the conversion script appends to existing output files.
cd NCBI_c2c
mkdir -p mpa_report

for file in *_NCBI_c2c.txt; do
    sample="$(basename "${file%_NCBI_c2c.txt}")"

    python3 /path/to/Scripts/Convert_MEGAN-NCBI-c2c_to_kreport-mpa_revised.py \
        -i "$file" \
        -c three \
        -l "mpa_report/${sample}" \
        -r 1
    
    # Remove intermediate files and the unused kreport output.
    rm -f \
        "mpa_report/${sample}.intermediate.names.txt" \
        "mpa_report/${sample}.intermediate.taxid.txt" \
        "mpa_report/${sample}.megan-c2c.kreport.txt"
    
    # Rename the MPA output for downstream processing.
    mv \
        "mpa_report/${sample}.megan-c2c.mpa.txt" \
        "mpa_report/${sample}.megan.mpa"
done
# The -r argument supplies the denominator for percentages in the kreport output. Here, -r 1 is used as a placeholder because the kreport output is discarded; it does not affect the MPA count values.

#### Generate bacterial and fungal count matrices at the family, genus, and species levels
# Combine the per-sample MPA reports and extract six contig-count matrices: bacterial and fungal matrices at the family, genus, and species levels. Each matrix contains taxa as rows and samples as columns.
cd mpa_report
# Run this command in the directory containing all .megan.mpa files:
bash /path/to/megan6_mpa_to_count_tables.sh /path/to/KrakenTools