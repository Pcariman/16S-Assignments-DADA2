#!/bin/bash

# ------------------------------------------------------------------------------
# Script: 01_cutadapt.sh
# Description: Trims primers from paired-end reads using Cutadapt.
# Project: 16S-Assignments-DADA2
# ------------------------------------------------------------------------------


# Define primers
F_PRIMER="CCTACGGGAGGCAGCAG"
R_PRIMER="GACTACAAGGGTATCTAATCC"
LENGTH="-m 1"
DISCARD="--discard-untrimmed"  # Mantener lecturas no recortadas
ERROR_RATE="-e 0.12"
MIN_OVERLAP="--overlap 5"

# Directory containing the FASTQ files
DIRECTORY="/raw_data/"

# Create a subdirectory for output files
OUTPUT_DIR="$DIRECTORY/2_trimmed_cutadapt"
mkdir -p "$OUTPUT_DIR"

# Loop through each pair of R1 and R2 files in the directory
for R1_INPUT in "$DIRECTORY"/*_R1_001.fastq; do
    R2_INPUT="${R1_INPUT/_R1_/_R2_}"
    BASENAME=$(basename "$R1_INPUT" _R1_001.fastq)
    R1_OUTPUT="$OUTPUT_DIR/${BASENAME}_R1_001.trimmed.fastq"
    R2_OUTPUT="$OUTPUT_DIR/${BASENAME}_R2_001.trimmed.fastq"

    cutadapt -g $F_PRIMER -G $R_PRIMER -o $R1_OUTPUT -p $R2_OUTPUT $R1_INPUT $R2_INPUT $LENGTH $DISCARD $ERROR_RATE $MIN_OVERLAP
done
