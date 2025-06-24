#!/bin/bash

# 1. Quality Control with FastQC
# This step assesses the quality of raw sequencing reads.
mkdir -p fastqc1
fastqc rawreads/*fastq -o fastqc1/

# 2. MultiQC Report
# Aggregates all FastQC reports into one summary HTML file.
multiqc fastqc1/

# 3. Adapter Trimming with Cutadapt
# Removes adapter sequences from small RNA reads.
mkdir -p trimmed
for file in fastq/*.fastq; do
    base=$(basename "$file" .fastq)
    cutadapt -a AGATCGGAAGAGC -o trimmed/${base}_trimmed.fastq "$file"
done

# 5. Build Bowtie Index from Arabidopsis Genome
# Required for mapping short reads in miRDeep2.
bowtie-build references/fna/GCF_000001735.4_TAIR10.1_genomic.fna references/fna/fna.index

# 6. Collapse and Map Reads using mapper.pl (miRDeep2)
# -e: convert to fasta
# -h: collapse reads
# -j: remove reads with non-canonical letters
# -m: discard reads shorter than min length
# -k: 3' adapter sequence
# -s: output collapsed reads
# -t: output ARF mapping format
# -p: bowtie index
# -v -n -o -r -l: various mapping options
mkdir -p collapsed
for file in trimmed/*trimmed.fastq; do
    sample=$(basename "$file" _trimmed.fastq)

    mapper.pl "$file" \
        -e -h -j -m \
        -k TGGAATTCTCGGGTGCCAAGG \
        -s collapsed/${sample}_collapsed.fa \
        -t collapsed/${sample}.arf \
        -p references/fna/fna.index \
        -v -n -o 4 -r 4 -l 18
done

# 7. Quantify known miRNAs using quantifier.pl
# This matches the collapsed reads to known mature and hairpin miRNAs.
# All output files (CSV, HTML, etc.) are written to quantifier_out/
mkdir -p quantifier_out
for file in collapsed/*_collapsed.fa; do
    sample=$(basename "$file" _collapsed.fa)

    (
        cd quantifier_out

        quantifier.pl \
            -p ../references/hairpin_ath.fa \
            -m ../references/mature_ath.fa \
            -r ../"$file" \
            -y "$sample"
    )
done

# Done! 🎉 Your quantification output is now in quantifier_out/
