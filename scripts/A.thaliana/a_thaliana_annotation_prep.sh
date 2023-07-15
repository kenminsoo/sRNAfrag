#!/bin/bash

# This is for the species
# A. thaliana

# download reference genome from NCBI, manually please
# https://www.ncbi.nlm.nih.gov/datasets/taxonomy/3702/
# Genbank - GTF and Fasta


#== Config == #
outdir=/Volumes/Extreme_SSD/a_thaliana
gtf_location=/Volumes/Extreme_SSD/a_thaliana/genomic.gtf
gtf_basename=tair
ref_genome=/Volumes/Extreme_SSD/a_thaliana/GCA_000001735.2_TAIR10.1_genomic.fna
#== Config End == #

# First remove comments
# python ../basics.py remove_comments $gtf_location $outdir"/"$gtf_basename"_1.gtf"

# First select transcript from the annotaion file
# Uses 0 as first index
# python ../gtf_modifiers.py select_column $outdir"/"$gtf_basename"_1.gtf" $outdir"/"$gtf_basename"_2.gtf" 2 transcript

# Now standardize the attributes to make it easier to work with
# python ../gtf_modifiers.py standardize_attributes $outdir"/"$gtf_basename"_2.gtf" $outdir"/"$gtf_basename"_3.gtf" '{"transcript_id":"transcript_id", "transcript_biotype":"biotype"}'

# Pull rRNA
# python ../gtf_modifiers.py select_gtf $outdir"/"$gtf_basename"_3.gtf" $outdir"/"$gtf_basename"_rRNA.gtf" biotype rRNA
 
# Pull snoRNA
# python ../gtf_modifiers.py select_gtf $outdir"/"$gtf_basename"_3.gtf" $outdir"/"$gtf_basename"_snoRNA.gtf" biotype snoRNA

# Pull snRNA
# python ../gtf_modifiers.py select_gtf $outdir"/"$gtf_basename"_3.gtf" $outdir"/"$gtf_basename"_snRNA.gtf" biotype snRNA

# add sequences
# python ../gtf_modifiers.py add_sequence_gtf $outdir"/"$gtf_basename"_rRNA.gtf" $ref_genome $outdir"/"$gtf_basename"_rRNA_final.gtf" sequence
# python ../gtf_modifiers.py add_sequence_gtf $outdir"/"$gtf_basename"_snoRNA.gtf" $ref_genome $outdir"/"$gtf_basename"_snoRNA_final.gtf" sequence
# python ../gtf_modifiers.py add_sequence_gtf $outdir"/"$gtf_basename"_snRNA.gtf" $ref_genome $outdir"/"$gtf_basename"_snRNA_final.gtf" sequence