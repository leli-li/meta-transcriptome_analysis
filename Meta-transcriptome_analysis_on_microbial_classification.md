# Meta-transcriptome analysis on microbial classification
# >>> STEP 1: Quality Control <<<
for r1 in *R1.fastq.gz; do
  base=${r1%%_R1*}
  trimmomatic PE -threads $THREADS \
    $r1 ${base}_R2.fastq.gz \
    ${base}_1P.fq.gz ${base}_R1.fq.gz \
    ${base}_2P.fq.gz ${base}_R2.fq.gz \
    ILLUMINACLIP:${REF_DIR}/adapters.fa:2:30:10 \
    SLIDINGWINDOW:5:20 MINLEN:70
done

# >>> STEP 2: Host Sequence Removal <<<
STAR --genomeDir ${REF_DIR}/host_index \
  --readFilesIn *R1.fq.gz \
  --readFilesCommand zcat \
  --outReadsUnmapped Fastx \
  --runThreadN $THREADS \
  --outSAMtype None 

# >>> STEP 3: Microbial Classification Using Kraken2 <<<
# Custom database based on NCBI-nr
kraken2 --db ${REF_DIR}/microbial_db \
  --paired *_unmapped_1.fq *_unmapped_2.fq \
  --threads $THREADS \
  --confidence $SENSITIVITY \
  --report kraken_report.tsv

# >>> STEP 4: Functional Annotation <<<
# DIAMOND+MEGAN6 pipeline for Microbial Classification
diamond blastx -d ${REF_DIR}/nr_fungal.dmnd \
  -q contigs.fa -o annotations.daa
daa-meganizer -i annotations.daa \
  -mdb ${REF_DIR}/megan-mapping.db