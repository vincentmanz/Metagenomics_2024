#!/bin/bash

# Define list of samples
SAMPLES='ERR315858'

# Activate the conda environment
# mamba activate QC_env

# Create output folder
mkdir -p TRIMMEDDATA

# Run Cutadapt
for SAMPLE in $SAMPLES; do 
  cutadapt RAWDATA/"$SAMPLE"_1.fastq.gz \
           RAWDATA/"$SAMPLE"_2.fastq.gz \
           -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
           -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
           -o TRIMMED/$SAMPLE.R1.fastq.gz \
           -p TRIMMED/$SAMPLE.R2.fastq.gz \
           --minimum-length 50 \
           --quality-cutoff 20 \
           --cores 4 > TRIMMED/$SAMPLE.cutadapt.log.txt
done