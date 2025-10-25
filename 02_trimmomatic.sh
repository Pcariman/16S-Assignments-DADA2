#!/bin/bash

# ------------------------------------------------------------------------------
# Script: 02_trimmomatic.sh
# Description: Performs quality trimming on paired-end reads using Trimmomatic.
# Project: 16S-Assignments-DADA2
# ------------------------------------------------------------------------------

# Define input and output directories (relative to the project root)
input_dir="/raw_data/2_trimmed_cutadapt/"
output_dir="/raw_data/2_trimmed_cutadapt/trimmomatic"

# Create output directory if it does not exist
mkdir -p $output_dir

# Loop over the files in the directory
for f1 in $input_dir/*_R1_*.fastq; do
# Find the corresponding paired file
  f2=${f1/_R1_/_R2_}

# Create output file names
  base_name=$(basename $f1 _R1_001.trimmed.fastq)
  newf1_paired="${output_dir}/${base_name}_R1_trimmomatic_paired.fastq"
  newf1_unpaired="${output_dir}/${base_name}_R1_trimmomatic_unpaired.fastq"
  newf2_paired="${output_dir}/${base_name}_R2_trimmomatic_paired.fastq"
  newf2_unpaired="${output_dir}/${base_name}_R2_trimmomatic_unpaired.fastq"

# Run Trimmomatic
  trimmomatic PE -threads 5 -phred33 -trimlog ${output_dir}/${base_name}_trimLogFile -summary ${output_dir}/${base_name}_statsSummaryFile \
    $f1 $f2 \
    $newf1_paired $newf1_unpaired \
    $newf2_paired $newf2_unpaired \
    LEADING:3 TRAILING:3 MINLEN:100
done
