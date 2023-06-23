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

#python ../conversion_tools.py gff3_to_gtf $gtf_location $outdir"/"$gtf_basename.gtf

# Now select pre-miRNAs

#python ../gtf_modifiers.py select_column $outdir"/"$gtf_basename.gtf $outdir"/"$gtf_basename"_1.gtf" 2 miRNA_primary_transcript

# Change NAME to transcript_id

#python ../gtf_modifiers.py standardize_attributes $outdir"/"$gtf_basename"_1.gtf" $outdir"/"$gtf_basename"_2.gtf" {Name:transcript_id}

# Add sequences to annotation

#python ../gtf_modifiers.py add_sequence_gtf $outdir"/"$gtf_basename"_2.gtf" $ref_genome $outdir"/"$gtf_basename"_3.gtf" sequence

# Now filter for mature miRNAs

#python ../gtf_modifiers.py select_column $outdir"/"$gtf_basename.gtf $outdir"/"$gtf_basename"_mature.gtf" 2 miRNA

# Make csv files for mature miRNAs and pre miRNAs to join in R or Python to get true start and ends

#python ../conversion_tools.py gtf_to_tsv $outdir"/"$gtf_basename"_mature.gtf" $outdir/mature.csv [0,3,4,6] [5,7] [chr,start,end,strand,name,derived_id]

#python ../conversion_tools.py gtf_to_tsv $outdir"/"$gtf_basename"_1.gtf" $outdir/pre.csv [0,3,4,6] [1,5] [chr,start,end,strand,name,derived_id]

# first combine the elements that we need to match up in mature dataset
# python ../gtf_modifiers.py merge_two_attributes $outdir"/"$gtf_basename"_mature.gtf" Derives_from Name deriv_id $outdir"/mirna_aasra.gtf"

# python ../conversion_tools.py gtf_to_fasta $outdir"/mirna_aasra.gtf" $outdir"/AASRA_mirna.fa" $ref_genome deriv_id

# From this point on
# Copy the miR_bench.R file wherever gtfs were saved & run
# Copy the peak files to the same location in the output after running pipeline