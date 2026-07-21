# 03_Kraken2_classification
# This step uses the non-plant reads generated after host-read removal and secondary plant-read filtering.
# Kraken2 is used for read-level taxonomic classification.
# In this workflow, a custom Kraken2 protein database was built from bacterial and fungal sequences extracted from the NCBI nr database.




## 03.1_Prepare_custom_bacterial_and_fungal_protein_database
### Download the NCBI nr database in parallel manually
#### build directory for NCBI-nr
mkdir -p "/path/to/ncbi_nr_blastdb"
cd "/path/to/ncbi_nr_blastdb"
#### download the multi-part compressed files
# Replace 116 with the final nr volume number available for the database snapshot being downloaded.
# The -P 4 option allows up to four files to be downloaded concurrently.
seq -w 0 116 | xargs -P 4 -I {} \
  wget -c "https://ftp.ncbi.nlm.nih.gov/blast/db/nr.{}.tar.gz" 
#### download the corresponding md5 files, which are used to check the integrity of the files
seq -w 0 116 | xargs -P 4 -I {} \
  wget -c "https://ftp.ncbi.nlm.nih.gov/blast/db/nr.{}.tar.gz.md5" 
#### Verify downloaded archives using MD5 checksums
md5sum --check nr.*.tar.gz.md5
#### Decompress the downloaded nr volumes
# Each archive is deleted only after successful extraction. 
for file in nr.*.tar.gz; do
  echo "Extracting $file..."

  if tar -xzf "$file"; then
    rm "$file"
    echo "$file extracted and deleted."
  else
    echo "Failed to extract $file." >&2
    exit 1
  fi
done
### Alternative: automatic download with update_blastdb.pl
# This official BLAST+ utility automatically identifies, downloads,verifies, and decompresses the required database volumes.
update_blastdb.pl --decompress nr

### Prepare taxid lists for bacteria and fungi
# This step uses get_species_taxids.sh to retrieve species-level taxids under each taxonomic group.
# Taxid 2 corresponds to Bacteria.
# Taxid 4751 corresponds to Fungi.
mkdir -p "/path/to/Kraken2_database_preparation/taxid_lists"
get_species_taxids.sh -t 2 > "/path/to/Kraken2_database_preparation/taxid_lists/bacteria_taxid_list.txt"
get_species_taxids.sh -t 4751 > "/path/to/Kraken2_database_preparation/taxid_lists/fungi_taxid_list.txt"

### Extract bacterial and fungal protein sequences from NCBI-nr using blastdbcmd command.
# The -taxidlist  restricts extraction to the taxids listed in the provided file.
# The -outfmt %f  outputs sequences in FASTA format.
# The -target_only option is used for the non-redundant nr database to limit output to target entries matching the requested taxids.
mkdir -p /path/to/Kraken2_database_preparation/nr_sequences
blastdbcmd \
  -db nr \
  -dbtype prot \
  -taxidlist "/path/to/Kraken2_database_preparation/taxid_lists/bacteria_taxid_list.txt" \
  -outfmt %f \
  -target_only \
  > "/path/to/Kraken2_database_preparation/nr_sequences/nr_bacteria.fa"
blastdbcmd \
  -db nr \
  -dbtype prot \
  -taxidlist "/path/to/Kraken2_database_preparation/taxid_lists/fungi_taxid_list.txt" \
  -outfmt %f \
  -target_only \
  > "/path/to/Kraken2_database_preparation/nr_sequences/nr_fungi.fa"




## 03.2_Build_custom_Kraken2_protein_database
### Kraken2 (v2.1.2) installation
conda install -c bioconda kraken2=2.1.2

### Build Kraken2 database
#### Create Kraken2 database directory
mkdir -p "/path/to/kraken2_database/nr_bacteria_and_fungi"
#### Download NCBI taxonomy using Kraken2-build function
# The --protein option is used because this custom database is built from protein sequences.
# Kraken2 will download (为什么我觉得这里应该用将来时？？？？) taxonomy files and accession-to-taxid mapping files into the database directory.
kraken2-build \
  --download-taxonomy \
  --db "/path/to/kraken2_database/nr_bacteria_and_fungi" \
  --protein
#### Add bacterial and fungal protein FASTA files to the Kraken2 database library
# The input FASTA files were extracted from the NCBI nr protein database.
# The --protein option indicates that the added library files contain amino acid sequences.
# The --no-masking option was used in this workflow to skip low-complexity masking during library addition.
# The two add-to-library commands are run sequentially to avoid concurrent writing to the same Kraken2 database directory.(感觉写两遍就是分开来运行的，也许可以删掉这句？)
kraken2-build \
  --add-to-library "/path/to/database_preparation/nr_bacteria.fa" \
  --db "/path/to/kraken2_database/nr_bacteria_and_fungi" \
  --protein \
  --no-masking \
  --threads 8 
kraken2-build 
  --add-to-library "/path/to/database_preparation/nr_fungi.fa" \
  --db "/path/to/kraken2_database/nr_bacteria_and_fungi" \
  --protein \
  --no-masking \
  --threads 8
#### Build the custom Kraken2 protein database
# The --protein option is used because the database is built from protein sequences.
# The --fast-build option may reduce database build time but should be recorded because it affects database construction.
kraken2-build 
--build \
--db "/path/to/kraken2_database/nr_bacteria_and_fungi" \
--protein \
--fast-build \
--threads 8 

### Clean intermediate Kraken2 database files
# This removes intermediate files generated during database construction while keeping the final Kraken2 database files.
kraken2-build \
--clean \
--db "/path/to/kraken2_database/nr_bacteria_and_fungi"




## 03.3_Run_Kraken2_classification
### Create output directory for Kraken2 classification results
mkdir -p "/path/to/kraken2_output"

### Run Kraken2 classification on paired non-plant reads
# [事后修改]⭐️The exact input filenames should match the final outputs documented in section 02.1. 

# The previous step generates paired-end non-plant reads in separate mate1 and mate2 files, which are used as input to Kraken2.
# The --paired option instructs Kraken2 to interpret mate1 and mate2 as reads originating from the same sequenced fragment.
# Kraken2 uses classification evidence from both mates to assign one classification to each read pair.
# The --report file summarizes taxonomic counts and percentages across all read pairs in the sample.
# The --output file contains one detailed classification record for each read pair.
# A confidence threshold of 0.05 was used in this workflow.
kraken2 \
  --db "/path/to/kraken2_database/nr_bacteria_and_fungi" \
  --threads 16 \
  --confidence 0.05 \
  --report "/path/to/kraken2_output/sample_kk2_nr.report" \
  --output "/path/to/kraken2_output/sample_kk2_nr.out" \
  --paired \
  "/path/to/nonplant_reads/sample_Unmapped.out.mate1_nonplant.fq" \
  "/path/to/nonplant_reads/sample_Unmapped.out.mate2_nonplant.fq"

### Process Kraken2 reports to generate a read-count matrix
# In this workflow, individual Kraken2 report files are first converted to MPA-style count tables and then combined into a cross-sample matrix.
# KrakenTools provides the kreport2mpa.py and combine_mpa.py scripts used for these two steps.
#### Download KrakenTools
# No separate installation is required; the KrakenTools scripts can be run directly with Python.
# The following git clone command creates the directory /path/to/software/KrakenTools and downloads the repository into it.
mkdir -p "/path/to/software"
git clone \
  "https://github.com/jenniferlu717/KrakenTools.git" \
  "/path/to/software/KrakenTools"

#### Run the report-processing script
# Enter the directory containing the Kraken2 *.report files.
cd "/path/to/kraken2_output"
# [事后修改]⭐️ Download the `kraken2_reports_to_count_tables.sh` script from the `scripts` directory of this repository, place it in the directory containing the Kraken2 `*.report` files, and run it from that directory.
# Specify the local KrakenTools directory as the first command-line argument.
# Run the following commands in a terminal.
bash "kraken2_reports_to_count_tables.sh" "/path/to/software/KrakenTools"
