#!/bin/bash

# This is for the species
# C. elegans

# First download annotation and ref genome from
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002985.6/ #

#== Config == #
outdir=/Volumes/Extreme_SSD/c_elegans
gtf_location=/Volumes/Extreme_SSD/c_elegans/genomic.gtf
gtf_basename=cel235
ref_genome=/Volumes/Extreme_SSD/c_elegans/GCA_000002985.3_WBcel235_genomic.fna
#== Config End == #

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