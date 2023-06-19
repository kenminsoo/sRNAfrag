#!/bin/bash

# This script processes the following annotation:
# https://www.mirbase.org/ftp/CURRENT/genomes/hsa.gff3 #
# To be ready for the sRNAfrag pipeline
# For peak detection benchmarking

## ENSURE THAT THE CONDA ENVIRONMENT IS ACTIVATED
## FROM THE GITHUB REPO

## CONFIG OPTIONS ##
outdir=/Users/kenminsoo/Desktop/unprocessed-annotations
gtf_location=/Users/kenminsoo/Desktop/unprocessed-annotations/mirbase.org_ftp_CURRENT_genomes_hsa.gff3
gtf_basename=hsa_mir
ref_genome=/Users/kenminsoo/Desktop/unprocessed-annotations/hg38_std.fa

# First convert to gtf

python base/conversion_tools.py gff3_to_gtf $gtf_location $outdir"/"$gtf_basename.gtf

# Now select pre-miRNAs

python base/gtf_modifiers.py select_column $outdir"/"$gtf_basename.gtf $outdir"/"$gtf_basename"_1.gtf" 2 miRNA_primary_transcript

# Change NAME to transcript_id

python base/gtf_modifiers.py standardize_attributes $outdir"/"$gtf_basename"_1.gtf" $outdir"/"$gtf_basename"_2.gtf" {Name:transcript_id}

# Add sequences to annotation

python base/gtf_modifiers.py add_sequence_gtf $outdir"/"$gtf_basename"_2.gtf" $ref_genome $outdir"/"$gtf_basename"_3.gtf" sequence