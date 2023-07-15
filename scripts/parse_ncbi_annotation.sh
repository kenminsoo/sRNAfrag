#!/bin/bash

# This is for the species
# M. musculus

# First download annotation and reference genome from
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/ #

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) gtf_location="$2"; shift ;;
        -f|--fasta) ref_genome="$2"; shift ;;
        -b|--basename) gtf_basename="$2"; shift ;;
        -o|--output) outdir="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# First remove comments
python ../basics.py remove_comments $gtf_location $outdir"/"$gtf_basename"_1.gtf"

# First select transcript from the annotaion file
# Uses 0 as first index
python ../gtf_modifiers.py select_column $outdir"/"$gtf_basename"_1.gtf" $outdir"/"$gtf_basename"_2.gtf" 2 transcript

# Now standardize the attributes to make it easier to work with
python ../gtf_modifiers.py standardize_attributes $outdir"/"$gtf_basename"_2.gtf" $outdir"/"$gtf_basename"_3.gtf" '{"transcript_id":"transcript_id", "transcript_biotype":"biotype"}'

# Pull rRNA
python ../gtf_modifiers.py select_gtf $outdir"/"$gtf_basename"_3.gtf" $outdir"/"$gtf_basename"_rRNA.gtf" biotype rRNA
 
# Pull snoRNA
python ../gtf_modifiers.py select_gtf $outdir"/"$gtf_basename"_3.gtf" $outdir"/"$gtf_basename"_snoRNA.gtf" biotype snoRNA

# Pull snRNA
python ../gtf_modifiers.py select_gtf $outdir"/"$gtf_basename"_3.gtf" $outdir"/"$gtf_basename"_snRNA.gtf" biotype snRNA

# add sequences
python ../gtf_modifiers.py add_sequence_gtf $outdir"/"$gtf_basename"_rRNA.gtf" $ref_genome $outdir"/"$gtf_basename"_rRNA_final.gtf" sequence
python ../gtf_modifiers.py add_sequence_gtf $outdir"/"$gtf_basename"_snoRNA.gtf" $ref_genome $outdir"/"$gtf_basename"_snoRNA_final.gtf" sequence
python ../gtf_modifiers.py add_sequence_gtf $outdir"/"$gtf_basename"_snRNA.gtf" $ref_genome $outdir"/"$gtf_basename"_snRNA_final.gtf" sequence